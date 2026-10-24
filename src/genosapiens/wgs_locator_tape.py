"""Complete fixed-width locator sequences with bounded decoding and seeking.

Locators address evaluated evidence; they are not an unevaluated enumeration
recipe. The owning ranking/selection component verifies their exact membership
and order. This codec verifies every stored occurrence and byte boundary.
"""
from array import array
import base64
import hashlib
import sys
import zlib

from genosapiens.wgs_evidence_store import EvidenceStoreError

CHUNK_ROWS = 4096
MAX_LOCATORS = 2**63 - 1


def _decode(content, count, width):
    if type(count) is not int or not 1 <= count <= CHUNK_ROWS or type(content) is not bytes:
        raise EvidenceStoreError("invalid locator chunk")
    decoder = zlib.decompressobj()
    raw = decoder.decompress(content,count*width*4+1)
    if len(raw) != count*width*4 or not decoder.eof or decoder.unused_data or decoder.unconsumed_tail:
        raise EvidenceStoreError("locator chunk expansion changed")
    values = array("I")
    values.frombytes(raw)
    if sys.byteorder != "little":
        values.byteswap()
    return raw,values


class LocatorTape:
    def __init__(self, writer, rows, *, width, name):
        self._initialize(writer,width,name)
        digest = hashlib.sha256()
        try:
            def chunks():
                pending = array("I")
                for row in rows:
                    if (type(row) is not tuple or len(row) != width
                            or any(type(v) is not int or not 0 <= v < 2**32 for v in row)):
                        raise EvidenceStoreError("invalid evaluated evidence locator")
                    self.count += 1
                    if self.count > MAX_LOCATORS:
                        raise EvidenceStoreError("locator sequence exceeds capacity",capacity=True)
                    pending.extend(row)
                    if len(pending) == CHUNK_ROWS*width:
                        yield self._store(pending,digest)
                        pending = array("I")
                if pending:
                    yield self._store(pending,digest)
            chunks_sha = writer.add_sequence(name + ".chunks","locator-chunks-v1",chunks())
            self.stream_sha256 = digest.hexdigest()
            self.sha256 = writer.add(name,"locator-tape-v1",{
                "schema":"genosapiens.wgs-locator-tape","version":1,"width":width,
                "row_count":self.count,"chunk_count":self.chunk_count,"chunk_rows":CHUNK_ROWS,
                "stream_sha256":self.stream_sha256,
                "chunks":{"name":name + ".chunks","sha256":chunks_sha}})
        except BaseException:
            writer.failed = True
            self.close()
            raise

    def _initialize(self, owner, width, name):
        if type(width) is not int or width not in (2,3):
            raise EvidenceStoreError("invalid locator width")
        self.connection,self.width,self.name = owner.connection,width,name
        self.table = "locators_" + hashlib.sha256(name.encode()).hexdigest()
        self.connection.execute(f"CREATE TABLE {self.table} (ordinal INTEGER PRIMARY KEY,count INTEGER,content BLOB)")
        self.closed,self.count,self.chunk_count = False,0,0

    def _store(self, values, digest):
        if sys.byteorder != "little":
            values.byteswap()
        raw = values.tobytes()
        digest.update(raw)
        self.chunk_count += 1
        count = len(values)//self.width
        content = zlib.compress(raw,6)
        self.connection.execute(f"INSERT INTO {self.table} VALUES (?,?,?)",(self.chunk_count,count,content))
        return {"ordinal":self.chunk_count,"count":count,"codec":"uint32le-zlib-base64",
                "sha256":hashlib.sha256(raw).hexdigest(),"data":base64.b64encode(content).decode("ascii")}

    def rows(self, *, start=0, count=None):
        if self.closed:
            raise EvidenceStoreError("locator tape is closed")
        count = self.count-start if count is None else count
        if (type(start) is not int or type(count) is not int or min(start,count) < 0
                or start+count > self.count):
            raise EvidenceStoreError("invalid locator range")
        emitted = 0
        if count:
            first,last = start//CHUNK_ROWS+1,(start+count-1)//CHUNK_ROWS+1
            # Release each cursor before yielding. A consumer may reconstruct
            # a context in the same connection and release its scratch tables.
            for ordinal in range(first,last+1):
                row = self.connection.execute(f"SELECT count,content FROM {self.table} WHERE ordinal=?",(ordinal,)).fetchone()
                if row is None:
                    raise EvidenceStoreError("locator chunk is missing")
                size,content = row
                _,values = _decode(content,size,self.width)
                lower = max(0,start-(ordinal-1)*CHUNK_ROWS)
                upper = min(size,start+count-(ordinal-1)*CHUNK_ROWS)
                for i in range(lower,upper):
                    emitted += 1
                    yield tuple(values[i*self.width:(i+1)*self.width])
        if emitted != count:
            raise EvidenceStoreError("locator range coverage changed")

    def close(self):
        if not self.closed:
            self.connection.execute("DROP TABLE " + self.table)
            self.closed = True

    def __enter__(self):
        return self

    def __exit__(self,*_):
        self.close()


class VerifiedLocatorTape(LocatorTape):
    def __init__(self, reader, reference, *, width):
        if reference.kind != "locator-tape-v1":
            raise EvidenceStoreError("foreign locator tape root")
        header = reader.reconstruct(reference)
        if (set(header) != {"schema","version","width","row_count","chunk_count","chunk_rows","stream_sha256","chunks"}
                or header["schema"] != "genosapiens.wgs-locator-tape"
                or type(header["version"]) is not int or header["version"] != 1
                or type(header["width"]) is not int or header["width"] != width
                or type(header["row_count"]) is not int or not 0 <= header["row_count"] <= MAX_LOCATORS
                or type(header["chunk_count"]) is not int
                or header["chunk_count"] != (header["row_count"]+CHUNK_ROWS-1)//CHUNK_ROWS
                or type(header["chunk_rows"]) is not int or header["chunk_rows"] != CHUNK_ROWS
                or set(header["chunks"]) != {"name","sha256"}
                or header["chunks"]["name"] != reference.name + ".chunks"):
            raise EvidenceStoreError("invalid locator tape header")
        child = reader.reference(header["chunks"]["name"],kind="locator-chunks-v1")
        if child.sha256 != header["chunks"]["sha256"]:
            raise EvidenceStoreError("locator chunk binding changed")
        self._initialize(reader,width,reference.name)
        digest = hashlib.sha256()
        try:
            for chunk in reader.iter_array(child,max_item_bytes=128*1024):
                if (set(chunk) != {"ordinal","count","codec","sha256","data"}
                        or type(chunk["ordinal"]) is not int or chunk["ordinal"] != self.chunk_count+1
                        or type(chunk["count"]) is not int
                        or chunk["count"] != min(CHUNK_ROWS,header["row_count"]-self.count)
                        or chunk["codec"] != "uint32le-zlib-base64" or type(chunk["data"]) is not str):
                    raise EvidenceStoreError("invalid locator chunk")
                try:
                    content = base64.b64decode(chunk["data"],validate=True)
                    raw,_ = _decode(content,chunk["count"],width)
                except (ValueError,zlib.error):
                    raise EvidenceStoreError("invalid locator chunk encoding") from None
                if (base64.b64encode(content).decode("ascii") != chunk["data"]
                        or hashlib.sha256(raw).hexdigest() != chunk["sha256"]):
                    raise EvidenceStoreError("locator chunk digest changed")
                digest.update(raw)
                self.count += chunk["count"]
                self.chunk_count += 1
                self.connection.execute(f"INSERT INTO {self.table} VALUES (?,?,?)",(self.chunk_count,chunk["count"],content))
            if (self.count != header["row_count"] or self.chunk_count != header["chunk_count"]
                    or digest.hexdigest() != header["stream_sha256"]):
                raise EvidenceStoreError("locator tape coverage changed")
            self.sha256,self.stream_sha256 = reference.sha256,digest.hexdigest()
        except BaseException:
            self.close()
            raise
