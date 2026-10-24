"""Bounded external merge of complete evaluated context rank families.

Only a fixed fan-in of private binary runs is decoded at once. Consumed input
chunks are deleted while merging so SQLite can reuse their pages. Portable
outputs retain all ranked occurrences as twelve-byte locators, plus their
computed total/tie families and every context-order dependency.
"""
import hashlib
from genosapiens.wgs_context_domain import next_context_ordinal, require_complete_domain
import heapq
import json
import struct
import zlib
from fractions import Fraction

from genosapiens.evidence_fusion import FusionPolicy, policy_as_json
from genosapiens.evidence_model import _typed
from genosapiens.evidence_serialization import _canonical
from genosapiens.wgs_context_order import MAX_ORDER_KEY_BYTES, _identity_bytes, order_sql_clauses
from genosapiens.wgs_evidence_replay import EvidenceReplay
from genosapiens.wgs_evidence_store import EvidenceStoreError
from genosapiens.wgs_locator_tape import LocatorTape

MERGE_FAN_IN = 64
RUN_CHUNK_ROWS = 4096
_RECORD = struct.Struct("<32sIII")


def _context_order_rows(ordinal, order, family):
    for locator in order.tape.rows(start=family["start"],count=family["count"]):
        yield _identity_bytes(order.identity(locator)),ordinal,*locator


def _context_order_header(ordinal, order):
    return {"ordinal":ordinal,"context_id":order.scores.universe.members.context_id,"hypothesis_count":order.count,
        "order":{"name":order.name,"sha256":order.sha256},
        "scores":{"name":order.scores.name,"sha256":order.scores.sha256}}


def _context_order_families(ordinal, order):
    for family in order.families():
        yield family,_context_order_rows(ordinal,order,family)


class CompleteCaseOrder:
    def __init__(self, writer, policy, *, policy_sha256, name, context_domain=None):
        self.context_domain = context_domain
        self.writer,self.connection,self.name = writer,writer.connection,name
        self.policy = _typed(policy,FusionPolicy)
        if hashlib.sha256(policy_as_json(self.policy)).hexdigest() != policy_sha256:
            raise EvidenceStoreError("case order policy changed")
        self.policy_sha256,self.rules = policy_sha256,self.policy.ties[:-1]
        self.prefix = "case_order_" + hashlib.sha256(name.encode()).hexdigest()
        self.closed,self.published,self.tables,self.tape = False,False,[],None
        self.context_count,self.hypothesis_count,self.family_count,self.run_count,self.merge_count = 0,0,0,0,0
        self.connection.create_collation("EXACT_FRACTION",lambda a,b:(Fraction(a)>Fraction(b))-(Fraction(a)<Fraction(b)))
        try:
            columns = ",".join(f"tie_{i} TEXT" for i in range(len(self.rules)))
            for suffix,definition in (
                ("contexts","ordinal INTEGER PRIMARY KEY,id TEXT UNIQUE,content TEXT"),
                ("families","id INTEGER PRIMARY KEY,content TEXT UNIQUE,total TEXT,count INTEGER,start INTEGER,rank INTEGER," + columns),
                ("runs","id INTEGER PRIMARY KEY,family INTEGER,count INTEGER"),
                ("chunks","run INTEGER,ordinal INTEGER,count INTEGER,digest BLOB,content BLOB,PRIMARY KEY(run,ordinal)"),
            ):
                table = self.prefix + "_" + suffix
                self.connection.execute(f"CREATE TABLE {table} ({definition})")
                self.tables.append(table)
            self.connection.execute(f"CREATE INDEX {self.prefix}_run_families ON {self.prefix}_runs(family,id)")
            self.connection.execute(f"CREATE UNIQUE INDEX {self.prefix}_family_rank ON {self.prefix}_families(rank)")
            self.connection.execute(f"CREATE UNIQUE INDEX {self.prefix}_family_start ON {self.prefix}_families(start)")
            self.connection.execute(f"CREATE INDEX {self.prefix}_family_values ON {self.prefix}_families(" + order_sql_clauses(self.rules) + ")")
        except BaseException:
            writer.failed = True
            self.close()
            raise

    def add_context(self, ordinal, order):
        if (self.closed or self.published or type(ordinal) is not int or ordinal != next_context_ordinal(self)
                or not ordinal < 2**32 or order.closed or order.connection is not self.connection
                or order.scores.policy_sha256 != self.policy_sha256):
            raise EvidenceStoreError("foreign context order or incomplete case sequence")
        try:
            self._append_context(_context_order_header(ordinal,order),_context_order_families(ordinal,order))
        except BaseException:
            self.writer.failed = True
            raise

    def _append_context(self, header, families):
        ordinal = header["ordinal"]
        if self.closed or self.published or type(ordinal) is not int or ordinal != next_context_ordinal(self) or ordinal >= 2**32:
            raise EvidenceStoreError("foreign context order or incomplete case sequence")
        try:
            self.connection.execute(f"INSERT INTO {self.prefix}_contexts VALUES (?,?,?)",
                (ordinal,header["context_id"],_canonical(header)))
            count = 0
            for family,rows in families:
                family_id = self._family(family["total"],tuple(family["ties"]))
                _,size = self._write_run(family_id,rows)
                if size != family["count"]:
                    raise EvidenceStoreError("context rank family count changed")
                self.connection.execute(f"UPDATE {self.prefix}_families SET count=count+? WHERE id=?",(size,family_id))
                count += size
            if count != header["hypothesis_count"]:
                raise EvidenceStoreError("case order context coverage changed")
            self.context_count += 1
            self.hypothesis_count += count
        except BaseException:
            self.writer.failed = True
            raise

    def _family(self, total, ties):
        content = _canonical((total,ties))
        if len(content) > MAX_ORDER_KEY_BYTES or len(ties) != len(self.rules):
            raise EvidenceStoreError("case order key exceeds capacity",capacity=True)
        row = self.connection.execute(f"SELECT id FROM {self.prefix}_families WHERE content=?",(content,)).fetchone()
        if row is not None:
            return row[0]
        self.family_count += 1
        self.connection.execute(f"INSERT INTO {self.prefix}_families VALUES (" + ",".join("?" for _ in range(6+len(ties))) + ")",
                                (self.family_count,content,total,0,None,None,*ties))
        return self.family_count

    def _write_run(self, family, rows):
        self.run_count += 1
        run,count,chunk,previous,pending = self.run_count,0,0,None,bytearray()
        def flush():
            nonlocal pending,chunk
            chunk += 1
            raw = bytes(pending)
            self.connection.execute(f"INSERT INTO {self.prefix}_chunks VALUES (?,?,?,?,?)",
                (run,chunk,len(raw)//_RECORD.size,hashlib.sha256(raw).digest(),zlib.compress(raw,6)))
            pending = bytearray()
        for identity,context,left,right in rows:
            if (type(identity) is not bytes or len(identity) != 32
                    or previous is not None and identity <= previous
                    or any(type(v) is not int or not 0 <= v < 2**32 for v in (context,left,right))
                    or context == 0 or left == 0):
                raise EvidenceStoreError("case rank identity order or locator changed")
            previous = identity
            pending.extend(_RECORD.pack(identity,context,left,right))
            count += 1
            if len(pending) == RUN_CHUNK_ROWS*_RECORD.size:
                flush()
        if pending:
            flush()
        if count == 0:
            raise EvidenceStoreError("empty rank family run")
        self.connection.execute(f"INSERT INTO {self.prefix}_runs VALUES (?,?,?)",(run,family,count))
        return run,count

    def _consume_run(self, run):
        expected = self.connection.execute(f"SELECT count FROM {self.prefix}_runs WHERE id=?",(run,)).fetchone()
        if expected is None:
            raise EvidenceStoreError("case rank run is unavailable")
        count,ordinal,previous = 0,0,None
        while True:
            row = self.connection.execute(f"SELECT ordinal,count,digest,content FROM {self.prefix}_chunks WHERE run=? AND ordinal=?",
                                          (run,ordinal+1)).fetchone()
            if row is None:
                break
            ordinal,size,digest,content = row
            if type(size) is not int or size != min(RUN_CHUNK_ROWS,expected[0]-count):
                raise EvidenceStoreError("case rank chunk count changed")
            decoder = zlib.decompressobj()
            raw = decoder.decompress(content,size*_RECORD.size+1)
            if (len(raw) != size*_RECORD.size or not decoder.eof or decoder.unused_data or decoder.unconsumed_tail
                    or hashlib.sha256(raw).digest() != digest):
                raise EvidenceStoreError("case rank chunk changed")
            # The bounded decoded window owns the consumed bytes from here.
            self.connection.execute(f"DELETE FROM {self.prefix}_chunks WHERE run=? AND ordinal=?",(run,ordinal))
            for value in _RECORD.iter_unpack(raw):
                if previous is not None and value[0] <= previous:
                    raise EvidenceStoreError("case rank run ordering changed")
                previous = value[0]
                count += 1
                yield value
        if count != expected[0]:
            raise EvidenceStoreError("case rank run coverage changed")
        self.connection.execute(f"DELETE FROM {self.prefix}_runs WHERE id=?",(run,))

    def _merged_family(self, family):
        while True:
            runs = tuple(row[0] for row in self.connection.execute(
                f"SELECT id FROM {self.prefix}_runs WHERE family=? ORDER BY id LIMIT ?",(family,MERGE_FAN_IN)))
            if len(runs) == 1:
                yield from self._consume_run(runs[0])
                return
            if not runs:
                raise EvidenceStoreError("case rank family has no run")
            self._write_run(family,heapq.merge(*(self._consume_run(run) for run in runs),key=lambda v:v[0]))
            self.merge_count += 1

    def publish(self):
        require_complete_domain(self)
        if self.closed or self.published:
            raise EvidenceStoreError("case order is unavailable")
        try:
            start = 0
            for rank,(family,count) in enumerate(self.connection.execute(
                    f"SELECT id,count FROM {self.prefix}_families ORDER BY " + order_sql_clauses(self.rules)),1):
                self.connection.execute(f"UPDATE {self.prefix}_families SET rank=?,start=? WHERE id=?",(rank,start,family))
                start += count
            if start != self.hypothesis_count:
                raise EvidenceStoreError("case order family totals changed")
            def rows():
                for family,expected in self.connection.execute(f"SELECT id,count FROM {self.prefix}_families ORDER BY rank"):
                    count = 0
                    for _,context,left,right in self._merged_family(family):
                        count += 1
                        yield context,left,right
                    if count != expected:
                        raise EvidenceStoreError("complete case rank family coverage changed")
            self.tape = LocatorTape(self.writer,rows(),width=3,name=self.name + ".locators")
            if self.tape.count != self.hypothesis_count:
                raise EvidenceStoreError("complete case order coverage changed")
            contexts_sha = self.writer.add_sequence(self.name + ".contexts","case-order-contexts-v1",
                (json.loads(row[0]) for row in self.connection.execute(f"SELECT content FROM {self.prefix}_contexts ORDER BY ordinal")))
            families_sha = self.writer.add_sequence(self.name + ".families","case-order-families-v1",self.families())
            self.sha256 = self.writer.add(self.name,"complete-case-order-v1",{
                "schema":"genosapiens.wgs-complete-case-order","version":1,
                "policy_sha256":self.policy_sha256,"tie_rules":self.policy.ties,
                "hypothesis_count":self.hypothesis_count,"context_count":self.context_count,"family_count":self.family_count,
                "contexts":{"name":self.name + ".contexts","sha256":contexts_sha},
                "families":{"name":self.name + ".families","sha256":families_sha},
                "locators":{"name":self.tape.name,"sha256":self.tape.sha256}})
            self.published = True
            return self.sha256
        except BaseException:
            self.writer.failed = True
            raise

    def families(self):
        if self.closed:
            raise EvidenceStoreError("case order is closed")
        for rank in range(1,self.family_count+1):
            row = self.connection.execute(f"SELECT content,start,count FROM {self.prefix}_families WHERE rank=?",(rank,)).fetchone()
            if row is None:
                raise EvidenceStoreError("case rank family is missing")
            content,start,count = row
            total,ties = json.loads(content)
            yield {"ordinal":rank,"total":total,"ties":ties,"start":start,"count":count}

    def at_rank(self, rank):
        if (self.closed or not self.published or type(rank) is not int
                or not 1 <= rank <= self.hypothesis_count):
            raise EvidenceStoreError("case rank is unavailable")
        row = self.connection.execute(f"SELECT content,start,count FROM {self.prefix}_families "
            "WHERE start<=? ORDER BY start DESC LIMIT 1",(rank-1,)).fetchone()
        if row is None or not row[1] <= rank-1 < row[1]+row[2]:
            raise EvidenceStoreError("case rank family coverage changed")
        total,ties = json.loads(row[0])
        locator, = self.tape.rows(start=rank-1,count=1)
        return locator,total,tuple(ties)

    @property
    def reference(self):
        if self.closed or not self.published:
            raise EvidenceStoreError("case order is unavailable")
        return {"name":self.name,"kind":"complete-case-order-v1","sha256":self.sha256}

    def context(self, ordinal):
        if self.closed or not self.published or type(ordinal) is not int:
            raise EvidenceStoreError("case order context is unavailable")
        row = self.connection.execute(f"SELECT content FROM {self.prefix}_contexts WHERE ordinal=?",(ordinal,)).fetchone()
        if row is None:
            raise EvidenceStoreError("case order context is unavailable")
        return json.loads(row[0])

    def close(self):
        if self.tape is not None:
            self.tape.close()
        for table in reversed(self.tables):
            self.connection.execute("DROP TABLE " + table)
        self.tables.clear()
        self.closed = True

    def __enter__(self):
        return self

    def __exit__(self,*_):
        self.close()


class CompletedCaseOrder(CompleteCaseOrder):
    """Load saved rank families and locators; never sort or score contexts."""
    def __init__(self, reader, reference, *, policy, policy_sha256, progress=None, context_domain=None):
        from genosapiens.wgs_completed_science import child, fields, header, integer, require, observed_rows
        from genosapiens.wgs_locator_tape import VerifiedLocatorTape
        from genosapiens.evidence_model import _exact
        value = header(reader, reference, kind="complete-case-order-v1",
            schema="genosapiens.wgs-complete-case-order", keys=("policy_sha256", "tie_rules",
                "hypothesis_count", "context_count", "family_count", "contexts", "families", "locators"))
        require(value["policy_sha256"] == policy_sha256 and _canonical(value["tie_rules"]) == _canonical(policy.ties))
        for key in ("hypothesis_count", "context_count", "family_count"):
            integer(value[key])
        super().__init__(EvidenceReplay(reader), policy, policy_sha256=policy_sha256, name=reference.name, context_domain=context_domain)
        try:
            for row in observed_rows(reader.iter_array(child(reader, value["contexts"], "case-order-contexts-v1")),
                    progress, 'saved_order', value['context_count']):
                fields(row, ("ordinal", "context_id", "order", "scores", "hypothesis_count"))
                require(integer(row["ordinal"], minimum=1) == next_context_ordinal(self))
                require(type(row["context_id"]) is str and bool(row["context_id"]))
                child(reader, row["order"], "complete-context-order-v1")
                child(reader, row["scores"], "complete-context-scores-v1")
                self.context_count += 1
                self.hypothesis_count += integer(row["hypothesis_count"])
                self.connection.execute(f"INSERT INTO {self.prefix}_contexts VALUES (?,?,?)",
                    (row["ordinal"], row["context_id"], _canonical(row)))
            total = 0
            for row in reader.iter_array(child(reader, value["families"], "case-order-families-v1")):
                fields(row, ("ordinal", "total", "ties", "start", "count"))
                ordinal = integer(row["ordinal"], minimum=1)
                require(ordinal == self.family_count + 1 and integer(row["start"]) == total)
                count = integer(row["count"], minimum=1)
                _exact(row["total"])
                require(type(row["ties"]) is list and len(row["ties"]) == len(self.rules))
                for tie, rule in zip(row["ties"], self.rules, strict=True):
                    require(tie is None or type(tie) is str)
                    if tie is not None and rule.kind == "exact":
                        _exact(tie)
                columns = ",".join(f"tie_{i}" for i in range(len(self.rules)))
                values = (ordinal, _canonical((row["total"], row["ties"])), row["total"], count,
                          total, ordinal, *row["ties"])
                self.connection.execute(f"INSERT INTO {self.prefix}_families "
                    f"(id,content,total,count,start,rank,{columns}) VALUES ({','.join('?' for _ in values)})", values)
                self.family_count += 1
                total += count
            self.tape = VerifiedLocatorTape(reader, child(reader, value["locators"], "locator-tape-v1"), width=3)
            require(all(getattr(self, key) == value[key] for key in ("context_count", "hypothesis_count", "family_count"))
                and total == self.hypothesis_count == self.tape.count, "completed case order coverage changed")
            self.sha256, self.published = reference.sha256, True
        except BaseException:
            self.close()
            raise


class VerifiedCaseOrder(CompleteCaseOrder):
    def __init__(self, reader, reference, context_orders, *, policy, policy_sha256):
        if reference.kind != "complete-case-order-v1":
            raise EvidenceStoreError("foreign complete case order root")
        reader._require_reference(reference)
        super().__init__(EvidenceReplay(reader),policy,policy_sha256=policy_sha256,name=reference.name)
        try:
            for ordinal,order in context_orders:
                self.add_context(ordinal,order)
            self.publish()
        except BaseException:
            self.close()
            raise
