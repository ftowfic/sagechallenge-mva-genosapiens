"""Concrete complete context orders from evaluated scalar score outcomes.

Each rank family has identical total and all ties except the final hypothesis
ID. Every hypothesis ID is computed, collision checked and actually sorted.
The portable result retains the full ordered eight-byte locator stream, not a
recipe to enumerate or sort pairs later. No rich pair hypothesis is constructed.
"""
from array import array
from dataclasses import dataclass
from fractions import Fraction
import hashlib
import json
import sys

from genosapiens.evidence_fusion import scientific_tie_key
from genosapiens.evidence_serialization import _canonical
from genosapiens.genotype_hypotheses import hypothesis_id_for_members
from genosapiens.wgs_evidence_replay import EvidenceReplay
from genosapiens.wgs_evidence_store import EvidenceStoreError
from genosapiens.wgs_locator_tape import LocatorTape

MAX_ORDER_FAMILIES = 1_000_000
MAX_ORDER_METADATA_BYTES = 64 * 1024**2
MAX_ORDER_KEY_BYTES = 64 * 1024


@dataclass(frozen=True, slots=True)
class _TieComponent:
    normalized: str | None
    neutral_value: str


def _identity_bytes(value):
    if type(value) is not str or len(value) != 75 or not value.startswith("hypothesis-"):
        raise EvidenceStoreError("invalid complete hypothesis identity")
    try:
        raw = bytes.fromhex(value[11:])
    except ValueError:
        raise EvidenceStoreError("invalid complete hypothesis identity") from None
    if raw.hex() != value[11:]:
        raise EvidenceStoreError("invalid complete hypothesis identity")
    return raw


def scalar_rank_values(score, inputs, context, policy):
    """Original tie kernel over computed values, without rich source ledgers."""
    neutral = {c.name:c.neutral_value for c in policy.components}
    neutral["phenotype"] = policy.phenotype.neutral
    components = {c.name:_TieComponent(c.normalized,neutral[c.name]) for c in score.components}
    return (score.total,tuple(scientific_tie_key(rule,context=context,
        completeness=inputs.completeness,phase=inputs.phase,hypothesis_id=None,components=components).value
        for rule in policy.ties[:-1]))


def order_sql_clauses(rules):
    clauses = ["total COLLATE EXACT_FRACTION DESC"]
    for i,rule in enumerate(rules):
        clauses.append(f"CASE WHEN tie_{i} IS NULL THEN 1 ELSE 0 END ASC")
        collation = " COLLATE EXACT_FRACTION" if rule.kind == "exact" else ""
        direction = "DESC" if rule.direction == "descending" else "ASC"
        clauses.append(f"tie_{i}{collation} {direction}")
    return ",".join(clauses)


class CompleteContextOrder:
    def __init__(self, writer, scores, *, name, retain_identity_index=False):
        self.connection,self.scores,self.name = writer.connection,scores,name
        self.prefix = "context_order_" + hashlib.sha256(name.encode()).hexdigest()
        self.closed,self.tables,self.tape = False,[],None
        self.count,self.family_count,self.metadata_bytes = 0,0,0
        self.members,self.member_ordinals,self.family_counts = [None],{},array("Q")
        if scores.closed or scores.connection is not self.connection:
            raise EvidenceStoreError("foreign context score owner for ordering")
        policy = scores.policy
        self.rules = policy.ties[:-1]
        self.connection.create_collation("EXACT_FRACTION",lambda a,b:(Fraction(a)>Fraction(b))-(Fraction(a)<Fraction(b)))
        try:
            columns = ",".join(f"tie_{i} TEXT" for i in range(len(self.rules)))
            for suffix,definition in (
                ("families","id INTEGER PRIMARY KEY,content TEXT UNIQUE,total TEXT,rank INTEGER,start INTEGER,count INTEGER," + columns),
                ("entries","hid BLOB PRIMARY KEY,family INTEGER,left_id INTEGER,right_id INTEGER"),
                ("singletons","ordinal INTEGER PRIMARY KEY,id TEXT UNIQUE"),
            ):
                table = self.prefix + "_" + suffix
                self.connection.execute(f"CREATE TABLE {table} ({definition})")
                self.tables.append(table)
            self.connection.execute(f"CREATE INDEX {self.prefix}_entry_order ON {self.prefix}_entries(family,hid)")
            self.connection.execute(f"CREATE UNIQUE INDEX {self.prefix}_family_rank ON {self.prefix}_families(rank)")
            self.connection.execute(f"CREATE INDEX {self.prefix}_family_values ON {self.prefix}_families(" + order_sql_clauses(self.rules) + ")")
            for ordinal,ledger in enumerate(scores.singleton_ledgers(),1):
                if (tuple((k.name,k.direction) for k in ledger.tie_keys) != tuple((r.field,r.direction) for r in policy.ties)
                        or ledger.tie_keys[-1].value != ledger.hypothesis_id):
                    raise EvidenceStoreError("singleton order tie contract changed")
                values = self._ordering_total(ledger.total),tuple(k.value for k in ledger.tie_keys[:-1])
                family = self._family(values)
                self.connection.execute(f"INSERT INTO {self.prefix}_singletons VALUES (?,?)",(ordinal,ledger.hypothesis_id))
                self._append(ledger.hypothesis_id,family,(ordinal,0))
            if scores.pairs is not None:
                pairs = scores.universe.pairs
                for ordinal,key,_,_,_ in pairs.domain.members():
                    if ordinal != len(self.members):
                        raise EvidenceStoreError("pair order member coverage changed")
                    self.metadata_bytes += sys.getsizeof(key)+128
                    if self.metadata_bytes > MAX_ORDER_METADATA_BYTES:
                        raise EvidenceStoreError("order member metadata exceeds capacity",capacity=True)
                    self.members.append(key)
                    self.member_ordinals[key] = ordinal
                families = array("I")
                for ordinal,inputs,decision in scores.pairs.catalog.records():
                    if ordinal != len(families):
                        raise EvidenceStoreError("pair score catalog ordering changed")
                    total,ties = scalar_rank_values(decision,inputs,pairs.context,policy)
                    families.append(self._family((self._ordering_total(total),ties)))
                for pair,code in zip(pairs.domain.pairs(),scores.pairs.buffer.values(),strict=True):
                    if scores.pairs.offset and code == 0:
                        continue
                    index = code-scores.pairs.offset
                    if not 0 <= index < len(families):
                        raise EvidenceStoreError("pair score assignment is unavailable")
                    identity = hypothesis_id_for_members(pairs.context.context_id,pair,pairs.facts.header.sample_name)
                    self._append(identity,families[index],tuple(self.member_ordinals[k] for k in pair))
            if self.count != scores.hypothesis_count or any(count == 0 for count in self.family_counts):
                raise EvidenceStoreError("complete order coverage changed")
            start = 0
            for rank,(family,) in enumerate(self.connection.execute(
                    f"SELECT id FROM {self.prefix}_families ORDER BY " + order_sql_clauses(self.rules)),1):
                count = self.family_counts[family]
                self.connection.execute(f"UPDATE {self.prefix}_families SET rank=?,start=?,count=? WHERE id=?",
                                        (rank,start,count,family))
                start += count
            self.tape = LocatorTape(writer,self._ordered_locators(),width=2,name=name + ".locators")
            if self.tape.count != self.count:
                raise EvidenceStoreError("complete order locator coverage changed")
            families_sha = writer.add_sequence(name + ".families","context-order-families-v1",self.families())
            self.sha256 = writer.add(name,"complete-context-order-v1",{
                "schema":"genosapiens.wgs-complete-context-order","version":1,
                "scores":{"name":scores.name,"sha256":scores.sha256},
                "policy_sha256":scores.policy_sha256,"tie_rules":policy.ties,
                "family_count":self.family_count,"hypothesis_count":self.count,
                "families":{"name":name + ".families","sha256":families_sha},
                "locators":{"name":self.tape.name,"sha256":self.tape.sha256}})
            # The only O(H) scratch relation is scoped to this context.
            if not retain_identity_index:
                self.connection.execute(f"DROP TABLE {self.prefix}_entries")
                self.tables.remove(self.prefix + "_entries")
        except BaseException:
            writer.failed = True
            self.close()
            raise

    def _family(self, values):
        content = _canonical(values)
        if len(content) > MAX_ORDER_KEY_BYTES:
            raise EvidenceStoreError("scientific order key exceeds capacity",capacity=True)
        row = self.connection.execute(f"SELECT id FROM {self.prefix}_families WHERE content=?",(content,)).fetchone()
        if row is not None:
            return row[0]
        if self.family_count == MAX_ORDER_FAMILIES:
            raise EvidenceStoreError("context order families exceed capacity",capacity=True)
        ordinal = self.family_count
        self.connection.execute(f"INSERT INTO {self.prefix}_families VALUES (" + ",".join("?" for _ in range(6+len(values[1]))) + ")",
            (ordinal,content,values[0],None,None,None,*values[1]))
        self.family_count += 1
        self.family_counts.append(0)
        return ordinal

    def _ordering_total(self, total):
        return total if not hasattr(self.scores,"ordering_total") else self.scores.ordering_total(total)

    def _append(self, identity, family, locator):
        self.connection.execute(f"INSERT INTO {self.prefix}_entries VALUES (?,?,?,?)",(_identity_bytes(identity),family,*locator))
        self.family_counts[family] += 1
        self.count += 1

    def _ordered_locators(self):
        previous_rank,emitted = 0,0
        # Both traversals use incrementally maintained disk indexes. A joined
        # ORDER BY across all H occurrences would require an unbounded SQLite
        # temporary sort, potentially on the runtime's small tmpfs.
        for family,rank,count in self.connection.execute(
                f"SELECT id,rank,count FROM {self.prefix}_families ORDER BY rank"):
            if rank < previous_rank:
                raise EvidenceStoreError("context order family order changed")
            previous_rank = rank
            observed = 0
            for left,right in self.connection.execute(
                    f"SELECT left_id,right_id FROM {self.prefix}_entries WHERE family=? ORDER BY hid",(family,)):
                observed += 1
                emitted += 1
                yield left,right
            if observed != count:
                raise EvidenceStoreError("context order family coverage changed")
        if emitted != self.count:
            raise EvidenceStoreError("context order enumeration changed")

    def families(self):
        if self.closed:
            raise EvidenceStoreError("context order is closed")
        for rank,content,start,count in self.connection.execute(
                f"SELECT rank,content,start,count FROM {self.prefix}_families ORDER BY rank"):
            total,ties = json.loads(content)
            yield {"ordinal":rank,"total":total,"ties":ties,"start":start,"count":count}

    def identity(self, locator):
        if self.closed or type(locator) is not tuple or len(locator) != 2 or any(type(v) is not int for v in locator):
            raise EvidenceStoreError("invalid context order locator")
        left,right = locator
        if right == 0:
            row = self.connection.execute(f"SELECT id FROM {self.prefix}_singletons WHERE ordinal=?",(left,)).fetchone()
            if row is None:
                raise EvidenceStoreError("unknown singleton order locator")
            return row[0]
        if not 1 <= left < right < len(self.members):
            raise EvidenceStoreError("unknown pair order locator")
        pairs = self.scores.universe.pairs
        return hypothesis_id_for_members(pairs.context.context_id,(self.members[left],self.members[right]),pairs.facts.header.sample_name)

    def identity_locators(self):
        """Original canonical hypothesis-ID order for explicit rich replay."""
        if self.closed or self.prefix + "_entries" not in self.tables:
            raise EvidenceStoreError("context identity traversal is unavailable")
        cursor = self.connection.execute(f"SELECT hid,left_id,right_id FROM {self.prefix}_entries ORDER BY hid")
        count = 0
        try:
            for digest,left,right in cursor:
                count += 1
                yield "hypothesis-" + digest.hex(),(left,right)
            if count != self.count:
                raise EvidenceStoreError("context identity traversal coverage changed")
        finally:
            cursor.close()

    def close(self):
        if self.tape is not None:
            self.tape.close()
        for table in reversed(self.tables):
            self.connection.execute("DROP TABLE " + table)
        self.tables.clear()
        self.members.clear()
        self.member_ordinals.clear()
        self.family_counts = array("Q")
        self.closed = True

    def __enter__(self):
        return self

    def __exit__(self,*_):
        self.close()


class VerifiedContextOrder(CompleteContextOrder):
    def __init__(self, reader, reference, scores, *, retain_identity_index=False):
        if reference.kind != "complete-context-order-v1":
            raise EvidenceStoreError("foreign context order root")
        reader._require_reference(reference)
        super().__init__(EvidenceReplay(reader),scores,name=reference.name,retain_identity_index=retain_identity_index)
