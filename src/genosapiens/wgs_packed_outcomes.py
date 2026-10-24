"""Disposable bounded chunks while a complete outcome catalog is discovered."""
from array import array
import hashlib
import sys
import zlib

from genosapiens.wgs_evidence_store import EvidenceStoreError
from genosapiens.wgs_pair_outcomes import CHUNK_PAIRS


class PackedOutcomes:
    def __init__(self, connection, name):
        self.connection = connection
        self.table = "packed_" + hashlib.sha256(name.encode()).hexdigest()
        connection.execute(f"CREATE TABLE {self.table} (ordinal INTEGER PRIMARY KEY,count INTEGER,content BLOB)")
        self.pending, self.chunks, self.count, self.closed = array("I"), 0, 0, False

    def append(self, value):
        if self.closed or type(value) is not int or not 0 <= value < 2**32:
            raise EvidenceStoreError("invalid private outcome code")
        self.pending.append(value)
        self.count += 1
        if len(self.pending) == CHUNK_PAIRS:
            self.flush()

    def flush(self):
        if self.closed:
            raise EvidenceStoreError("private outcome buffer is closed")
        if not self.pending:
            return
        if sys.byteorder != "little":
            self.pending.byteswap()
        self.chunks += 1
        self.connection.execute(f"INSERT INTO {self.table} VALUES (?,?,?)",
            (self.chunks,len(self.pending),zlib.compress(self.pending.tobytes(),6)))
        self.pending = array("I")

    def values(self):
        self.flush()
        emitted = 0
        for count, content in self.connection.execute(f"SELECT count,content FROM {self.table} ORDER BY ordinal"):
            if type(count) is not int or not 1 <= count <= CHUNK_PAIRS:
                raise EvidenceStoreError("private outcome chunk capacity changed")
            decoder = zlib.decompressobj()
            raw = decoder.decompress(content,count*4+1)
            if len(raw) != count*4 or not decoder.eof or decoder.unused_data or decoder.unconsumed_tail:
                raise EvidenceStoreError("private outcome chunk changed")
            values = array("I")
            values.frombytes(raw)
            if sys.byteorder != "little":
                values.byteswap()
            emitted += len(values)
            yield from values
        if emitted != self.count:
            raise EvidenceStoreError("private outcome coverage changed")

    def value_at(self, index):
        self.flush()
        if type(index) is not int or not 0 <= index < self.count:
            raise EvidenceStoreError("private outcome position is unavailable")
        count,content = self.connection.execute(f"SELECT count,content FROM {self.table} WHERE ordinal=?",
                                               (index//CHUNK_PAIRS+1,)).fetchone()
        return outcome_chunk_values(count,content)[index%CHUNK_PAIRS]

    def close(self):
        if not self.closed:
            self.connection.execute("DROP TABLE " + self.table)
        self.pending = array("I")
        self.closed = True


def outcome_chunk_values(count, content):
    if type(count) is not int or not 1 <= count <= CHUNK_PAIRS:
        raise EvidenceStoreError("private outcome chunk capacity changed")
    decoder = zlib.decompressobj()
    raw = decoder.decompress(content,count*4+1)
    if len(raw) != count*4 or not decoder.eof or decoder.unused_data or decoder.unconsumed_tail:
        raise EvidenceStoreError("private outcome chunk changed")
    values = array("I")
    values.frombytes(raw)
    if sys.byteorder != "little":
        values.byteswap()
    return values
