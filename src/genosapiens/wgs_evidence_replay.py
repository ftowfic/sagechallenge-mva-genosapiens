"""Compare computed component roots with one authenticated saved owner.

This sink writes no evidence. Its add methods close every computed logical
value against the saved root and return that root's actual digest. It permits
component replay to reuse the original construction code, with bounded rows.
"""
import hashlib
from itertools import zip_longest

from genosapiens.evidence_serialization import _canonical
from genosapiens.wgs_evidence_store import EvidenceReader, EvidenceStoreError, MAX_RECONSTRUCT_BYTES


class EvidenceReplay:
    def __init__(self, reader):
        if type(reader) is not EvidenceReader or not reader.verified:
            raise EvidenceStoreError("authenticated evidence reader required for replay")
        self.reader, self.connection, self.temporary = reader, reader.connection, reader.temporary
        self.failed = False

    def add(self, name, kind, value):
        reference = self.reader.reference(name, kind=kind)
        if self.reader._skip_verification:
            return reference.sha256
        content = _canonical(value)
        if len(content) > MAX_RECONSTRUCT_BYTES:
            raise EvidenceStoreError("replayed evidence window exceeds capacity", capacity=True)
        if content != _canonical(self.reader.reconstruct(reference)):
            raise EvidenceStoreError("saved evidence root does not replay")
        return reference.sha256

    def add_sequence(self, name, kind, values):
        reference = self.reader.reference(name, kind=kind)
        if self.reader._skip_verification:
            # Constructors can populate private projection tables while their
            # sequence runs. Keep those side effects, but omit serialization,
            # saved-row reconstruction, and the comparison itself.
            iterator = iter(values)
            try:
                for _ in iterator:
                    pass
            finally:
                close = getattr(iterator, 'close', None)
                if close is not None:
                    close()
            return reference.sha256
        absent = object()
        expected_rows,actual_rows = iter(values),self.reader.iter_array(reference)
        try:
            for expected, actual in zip_longest(expected_rows, actual_rows, fillvalue=absent):
                if expected is absent or actual is absent:
                    raise EvidenceStoreError("saved evidence sequence coverage does not replay")
                content = _canonical(expected)
                if len(content) > MAX_RECONSTRUCT_BYTES:
                    raise EvidenceStoreError("replayed evidence row exceeds capacity", capacity=True)
                if content != _canonical(actual):
                    raise EvidenceStoreError("saved evidence sequence does not replay")
        finally:
            # A mismatch may suspend a source cursor. Close both iterators
            # before context owners drop their private SQLite tables, even
            # when the exception traceback still retains this frame.
            for iterator in (expected_rows,actual_rows):
                close = getattr(iterator,"close",None)
                if close is not None:
                    close()
        return reference.sha256


class ScopedSequenceReplay(EvidenceReplay):
    """Reuse bounded row commitments within an enclosing authenticated read.

    The enclosing scope must check all source bytes before accepting results.
    Only sequences observed during their first necessary read are reused.
    """
    def __init__(self, reader, *, max_sequences, max_item_bytes):
        super().__init__(reader)
        self.depth = reader._read_scope_depth
        self.derivative = reader._active_derivative
        self.rows = {}
        self.max_sequences, self.max_item_bytes = max_sequences, max_item_bytes
        self._require_scope()

    def _require_scope(self):
        if (not self.reader.verified or self.reader.connection is not self.connection
                or self.reader._active_derivative is not self.derivative
                or self.reader._read_scope_depth < self.depth
                or not (self.depth or self.derivative is not None)):
            raise EvidenceStoreError('sequence replay input scope changed')

    def observe(self, reference, values):
        self._require_scope()
        if self.reader._skip_verification:
            yield from values
            return
        if reference.name in self.rows or len(self.rows) >= self.max_sequences:
            raise EvidenceStoreError('sequence replay commitments exceed capacity', capacity=True)
        digest, count = hashlib.sha256(), 0
        for value in values:
            digest.update((_canonical(value) + '\n').encode('ascii'))
            count += 1
            yield value
        self.rows[reference.name] = (reference, count, digest.digest())

    def add_sequence(self, name, kind, values):
        saved = self.rows.get(name)
        if saved is None:
            return super().add_sequence(name, kind, values)
        self._require_scope()
        reference, count, expected = saved
        if self.reader.reference(name, kind=kind) != reference:
            raise EvidenceStoreError('replay sequence changed')
        digest, actual_count = hashlib.sha256(), 0
        iterator = iter(values)
        try:
            for value in iterator:
                content = _canonical(value)
                if len(content) > self.max_item_bytes:
                    raise EvidenceStoreError('replay item exceeds capacity', capacity=True)
                digest.update((content + '\n').encode('ascii'))
                actual_count += 1
            if actual_count != count or digest.digest() != expected:
                raise EvidenceStoreError('saved sequence does not replay')
        finally:
            close = getattr(iterator, 'close', None)
            if close is not None:
                close()
        return reference.sha256
