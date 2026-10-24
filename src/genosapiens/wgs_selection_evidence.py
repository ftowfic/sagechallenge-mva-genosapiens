"""Complete compact selection audit with the existing ordered group rules.

The source adapter supplies each authenticated occurrence in the complete case
order. Identity keys are case-wide interned rendered alleles; context member
keys remain separate for the original complete-pair/partial-member rule.
Only final selected ranks live in memory. Every audit occurrence is retained
as a decision code and exact primary rank, without repeating rich evidence.
"""
from collections import OrderedDict
from dataclasses import dataclass
import hashlib
import json
import struct

from genosapiens.evidence_serialization import _canonical
from genosapiens.wgs_evidence_store import EvidenceStoreError
from genosapiens.wgs_evidence_replay import EvidenceReplay
from genosapiens.wgs_locator_tape import LocatorTape, MAX_LOCATORS
from genosapiens.wgs_packed_outcomes import PackedOutcomes

MAX_SELECTION_CLASSES = 1_000_000
MAX_SELECTION_REASON_BYTES = 64 * 1024
_COMMON = (
    ("submitted",()),
    ("secondary",("outside-primary-row-limit",)),
    ("secondary",("alternate-interpretation-of-submitted-identity",)),
    ("secondary",("alternate-interpretation-of-secondary-identity",)),
    ("excluded",("complete-pair-not-split",)),
)


@dataclass(frozen=True, slots=True)
class SelectionOccurrence:
    locator: tuple[int,int,int]
    member_keys: tuple[int,...]
    rendered_keys: tuple[int,...]
    complete: bool
    partial: bool
    reasons: tuple[str,...]


class CompleteSelectionEvidence:
    def __init__(self, writer, order, source_reference, occurrences, *, name):
        self.writer,self.connection,self.name = writer,writer.connection,name
        self.prefix = "selection_" + hashlib.sha256(name.encode()).hexdigest()
        self.tables,self.closed,self.buffer,self.tape,self.hot = [],False,None,None,OrderedDict()
        self.catalog_count,self.count,self.tail_count,self.excluded_count,self.alternative_count = 0,0,0,0,0
        self.selected,self.removed = [],set()
        try:
            if order.closed or not order.published or order.connection is not self.connection:
                raise EvidenceStoreError("foreign selection order owner")
            if (type(source_reference) is not dict or set(source_reference) != {"name","kind","sha256"}
                    or source_reference["kind"] != "case-submission-inputs-v1"
                    or self.connection.execute("SELECT kind,digest FROM roots WHERE name=?",(source_reference["name"],)).fetchone()
                    != (source_reference["kind"],source_reference["sha256"])):
                raise EvidenceStoreError("foreign selection source root")
            for suffix,definition in (("groups","identity BLOB PRIMARY KEY,rank INTEGER,selected INTEGER"),
                                      ("catalog","id INTEGER PRIMARY KEY,content TEXT UNIQUE")):
                table = self.prefix + "_" + suffix
                self.connection.execute(f"CREATE TABLE {table} ({definition})")
                self.tables.append(table)
            for expected,(category,reasons) in enumerate(_COMMON):
                if self._decision(category,reasons) != expected:
                    raise EvidenceStoreError("selection decision catalog changed")
            self.buffer = PackedOutcomes(self.connection,name)
            selected_pairs,selected_singletons = set(),{}
            for row,locator in zip(occurrences,order.tape.rows(),strict=True):
                self.count += 1
                if self.count > MAX_LOCATORS:
                    raise EvidenceStoreError("selection count exceeds capacity",capacity=True)
                self._require_occurrence(row,locator)
                rank,primary = self.count,0
                context = locator[0]
                if row.reasons:
                    code = self._decision("excluded",row.reasons)
                elif len(row.member_keys) == 1 and row.partial and (context,row.member_keys[0]) in selected_pairs:
                    code = 4
                else:
                    identity = struct.pack("<"+"I"*len(row.rendered_keys),*row.rendered_keys)
                    existing = self.connection.execute(f"SELECT rank,selected FROM {self.prefix}_groups WHERE identity=?",(identity,)).fetchone()
                    if existing is not None:
                        primary,state = existing
                        if state == -1:
                            code,primary = 4,0
                        else:
                            code = 2 if state else 3
                    else:
                        selected = len(self.selected) < 10
                        self.connection.execute(f"INSERT INTO {self.prefix}_groups VALUES (?,?,?)",(identity,rank,int(selected)))
                        if selected:
                            if len(row.member_keys) == 2 and row.complete:
                                for member in row.member_keys:
                                    key = context,member
                                    previous = selected_singletons.pop(key,None)
                                    if previous is not None:
                                        earlier,earlier_identity = previous
                                        self.selected.remove(earlier)
                                        self.removed.add(earlier)
                                        self.connection.execute(f"UPDATE {self.prefix}_groups SET selected=-1 WHERE identity=?",(earlier_identity,))
                                    selected_pairs.add(key)
                            self.selected.append(rank)
                            if len(row.member_keys) == 1 and row.partial:
                                selected_singletons[context,row.member_keys[0]] = rank,identity
                            code = 0
                        else:
                            code = 1
                self.buffer.append(code)
                self.buffer.append(primary & 0xffffffff)
                self.buffer.append(primary >> 32)
            if self.count != order.hypothesis_count or len(self.removed) > 20 or len(selected_pairs) > 20:
                raise EvidenceStoreError("complete selection coverage changed")
            def final_rows():
                values = self.buffer.values()
                for rank,(code,low,high) in enumerate(zip(values,values,values,strict=True),1):
                    primary = low+(high<<32)
                    if rank in self.removed or primary in self.removed:
                        code,low,high,primary = 4,0,0,0
                    self.excluded_count += code >= 4
                    self.tail_count += code == 1
                    self.alternative_count += code in (2,3)
                    yield code,low,high
            self.tape = LocatorTape(writer,final_rows(),width=3,name=name + ".audit")
            if (self.tape.count != self.count or len(self.selected)+self.tail_count+self.alternative_count+self.excluded_count != self.count):
                raise EvidenceStoreError("selection final accounting changed")
            self.buffer.close()
            self.buffer = None
            catalog_sha = writer.add_sequence(name + ".catalog","selection-decisions-v1",self.catalog())
            self.sha256 = writer.add(name,"complete-submission-selection-v1",{
                "schema":"genosapiens.wgs-submission-selection","version":1,"order":order.reference,
                "inputs":source_reference,"selected_ranks":tuple(self.selected),"hypothesis_count":self.count,
                "selected_count":len(self.selected),"secondary_group_count":self.tail_count,
                "excluded_count":self.excluded_count,"alternate_interpretation_count":self.alternative_count,
                "catalog":{"name":name + ".catalog","sha256":catalog_sha},
                "audit":{"name":self.tape.name,"sha256":self.tape.sha256}})
        except BaseException:
            writer.failed = True
            self.close()
            raise

    @staticmethod
    def _require_occurrence(row, locator):
        if (type(row) is not SelectionOccurrence or row.locator != locator or type(row.complete) is not bool
                or type(row.partial) is not bool or row.partial and row.complete
                or type(row.member_keys) is not tuple or type(row.rendered_keys) is not tuple
                or any(type(v) is not int or not 0 < v < 2**32 for v in (*row.member_keys,*row.rendered_keys))
                or len(set(row.member_keys)) != len(row.member_keys)
                or type(row.reasons) is not tuple
                or any(type(v) is not str or not v for v in row.reasons)
                or len(set(row.reasons)) != len(row.reasons)
                or not row.reasons and (len(row.member_keys) not in (1,2) or len(row.rendered_keys) != len(row.member_keys))):
            raise EvidenceStoreError("selection occurrence dependencies changed")

    def _decision(self, category, reasons):
        content = _canonical((category,reasons))
        if len(content) > MAX_SELECTION_REASON_BYTES:
            raise EvidenceStoreError("selection reasons exceed capacity",capacity=True)
        cached = self.hot.get(content)
        if cached is not None:
            return cached
        row = self.connection.execute(f"SELECT id FROM {self.prefix}_catalog WHERE content=?",(content,)).fetchone()
        if row is None:
            if self.catalog_count == MAX_SELECTION_CLASSES:
                raise EvidenceStoreError("selection decisions exceed capacity",capacity=True)
            code = self.catalog_count
            self.connection.execute(f"INSERT INTO {self.prefix}_catalog VALUES (?,?)",(code,content))
            self.catalog_count += 1
        else:
            code, = row
        self.hot[content] = code
        if len(self.hot) > 256:
            self.hot.popitem(last=False)
        return code

    def catalog(self):
        for code,content in self.connection.execute(f"SELECT * FROM {self.prefix}_catalog ORDER BY id"):
            category,reasons = json.loads(content)
            yield {"ordinal":code,"category":category,"reasons":reasons}

    def audit(self):
        if self.closed:
            raise EvidenceStoreError("selection audit is closed")
        for rank,(code,low,high) in enumerate(self.tape.rows(),1):
            content, = self.connection.execute(f"SELECT content FROM {self.prefix}_catalog WHERE id=?",(code,)).fetchone()
            category,reasons = json.loads(content)
            yield rank,category,tuple(reasons),low+(high<<32)

    def close(self):
        if self.buffer is not None:
            self.buffer.close()
        if self.tape is not None:
            self.tape.close()
        for table in reversed(self.tables):
            self.connection.execute("DROP TABLE " + table)
        self.tables.clear()
        self.hot.clear()
        self.closed = True

    def __enter__(self):
        return self

    def __exit__(self,*_):
        self.close()


class CompletedSelectionEvidence(CompleteSelectionEvidence):
    """Consume the saved decision audit without running selection again."""
    def __init__(self, reader, reference, order, source_reference):
        from genosapiens.wgs_completed_science import child, fields, header, integer, require
        from genosapiens.wgs_locator_tape import VerifiedLocatorTape
        value = header(reader, reference, kind="complete-submission-selection-v1",
            schema="genosapiens.wgs-submission-selection", keys=("order", "inputs", "selected_ranks",
                "hypothesis_count", "selected_count", "secondary_group_count", "excluded_count",
                "alternate_interpretation_count", "catalog", "audit"))
        require(value["order"] == order.reference and value["inputs"] == source_reference)
        require(type(value["selected_ranks"]) is list and len(value["selected_ranks"]) <= 10)
        self.selected = value["selected_ranks"]
        require(all(integer(rank, minimum=1) <= order.hypothesis_count for rank in self.selected)
            and self.selected == sorted(set(self.selected)))
        self.writer, self.connection, self.name = reader, reader.connection, reference.name
        self.prefix = "selection_" + hashlib.sha256(self.name.encode()).hexdigest()
        self.tables, self.closed, self.buffer, self.tape, self.hot = [], False, None, None, OrderedDict()
        self.catalog_count = 0
        self.count = integer(value["hypothesis_count"])
        self.tail_count = integer(value["secondary_group_count"])
        self.excluded_count = integer(value["excluded_count"])
        self.alternative_count = integer(value["alternate_interpretation_count"])
        require(self.count == order.hypothesis_count and integer(value["selected_count"]) == len(self.selected)
            and len(self.selected) + self.tail_count + self.excluded_count + self.alternative_count == self.count)
        try:
            table = self.prefix + "_catalog"
            self.connection.execute(f"CREATE TABLE {table} (id INTEGER PRIMARY KEY,content TEXT UNIQUE)")
            self.tables.append(table)
            for row in reader.iter_array(child(reader, value["catalog"], "selection-decisions-v1"),
                    max_item_bytes=MAX_SELECTION_REASON_BYTES):
                fields(row, ("ordinal", "category", "reasons"))
                require(integer(row["ordinal"]) == self.catalog_count < MAX_SELECTION_CLASSES)
                require(row["category"] in ("submitted", "secondary", "excluded")
                    and type(row["reasons"]) is list and all(type(v) is str and v for v in row["reasons"])
                    and len(row["reasons"]) == len(set(row["reasons"])))
                if self.catalog_count < len(_COMMON):
                    require((row["category"], tuple(row["reasons"])) == _COMMON[self.catalog_count])
                else:
                    require(row["category"] == "excluded" and row["reasons"])
                self.connection.execute(f"INSERT INTO {table} VALUES (?,?)",
                    (row["ordinal"], _canonical((row["category"], row["reasons"]))))
                self.catalog_count += 1
            require(self.catalog_count >= len(_COMMON))
            self.tape = VerifiedLocatorTape(reader, child(reader, value["audit"], "locator-tape-v1"), width=3)
            selected, tail, excluded, alternatives = [], 0, 0, 0
            for rank, (code, low, high) in enumerate(self.tape.rows(), 1):
                primary_rank = low + (high << 32)
                require(code < self.catalog_count and (0 < primary_rank < rank if code in (2, 3) else primary_rank == 0))
                if code == 0:
                    selected.append(rank)
                    require(len(selected) <= 10)
                tail += code == 1
                alternatives += code in (2, 3)
                excluded += code >= 4
            require(self.tape.count == self.count and selected == self.selected
                and (tail, excluded, alternatives) == (self.tail_count, self.excluded_count, self.alternative_count))
            self.sha256 = reference.sha256
        except BaseException:
            self.close()
            raise


class VerifiedSelectionEvidence(CompleteSelectionEvidence):
    def __init__(self, reader, reference, order, source_reference, occurrences):
        if reference.kind != "complete-submission-selection-v1":
            raise EvidenceStoreError("foreign complete selection root")
        reader._require_reference(reference)
        super().__init__(EvidenceReplay(reader),order,source_reference,occurrences,name=reference.name)
