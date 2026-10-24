"""Disk-backed source/member genotype facts for a globally closed context.

Calls are grouped on disk before interpretation. The shared scientific kernels
own source GT, exact-event copy aggregation and phase; this owner only supplies
their complete inputs. Rich admission/annotation evidence is never accumulated
for the whole context. Pair evaluation is a subsequent complete-domain stage.
"""
from __future__ import annotations

from collections import OrderedDict
from dataclasses import dataclass
import hashlib
import json
from types import SimpleNamespace

from genosapiens.allele_identity import MosaicAllele, RepeatAllele
from genosapiens.evidence_model import SourcedRecord, _sha, _typed
from genosapiens.evidence_serialization import _canonical
from genosapiens.genotype_properties import (
    POLICY_SHA256, AlleleGenotypeProperties, SourceGenotypeProperties, SampleCallQc,
    _LedgerSourceView, _allele_event, _identity_contig, ledger_call_binding,
    evaluate_source_genotype, evaluate_allele_genotype, assess_pair_phase,
)
from genosapiens import longtr_admission as longtr
from genosapiens.wgs_contexts import ClosedContext
from genosapiens.wgs_evidence_store import EvidenceStoreError

MAX_FACT_BYTES = 64 * 1024**2
MAX_SOURCE_WINDOW_BYTES = 64 * 1024**2
FACT_CACHE_BYTES = 16 * 1024**2
MAX_PHASE_BASES = 100_000
INVENTORY_SCHEMA = "genosapiens.wgs-context-genotype-inventory"


@dataclass(frozen=True, slots=True, repr=False)
class GenotypeFactHeader:
    inventory_sha256: str
    input_sha256: str
    reference_manifest_sha256: str
    sample_name: str
    policy_sha256: str
    ledger_source_samples: tuple[tuple[str,str], ...]
    source_set: None = None
    inventory_schema: str = INVENTORY_SCHEMA
    inventory_version: int = 1

    def __post_init__(self):
        for value in (self.inventory_sha256,self.input_sha256,self.reference_manifest_sha256,self.policy_sha256):
            _sha(value)
        for value,sample in self.ledger_source_samples:
            _sha(value)
            if type(sample) is not str or not sample:
                raise EvidenceStoreError("invalid genotype source sample")
        if (self.inventory_schema != INVENTORY_SCHEMA or type(self.inventory_version) is not int or self.inventory_version != 1
                or self.policy_sha256 != POLICY_SHA256 or not self.sample_name or self.source_set is not None
                or self.ledger_source_samples != tuple(sorted(self.ledger_source_samples))
                or len({k for k,_ in self.ledger_source_samples}) != len(self.ledger_source_samples)):
            raise EvidenceStoreError("invalid genotype fact header")


class ContextGenotypeFacts:
    """One closed context's bounded, disposable query projection."""
    def __init__(self, writer, closed, *, sample_name, source_manifest_sha256, reference_manifest_sha256):
        from genosapiens.wgs_context_evidence import EnrichedContextEvidence, VerifiedContextEvidence
        if type(closed) not in (ClosedContext,EnrichedContextEvidence,VerifiedContextEvidence) or not sample_name:
            raise EvidenceStoreError("complete genotype context is required")
        # A globally closed empty context retains its zero-member inventory
        # and source parents. It evaluates no genotypes and asserts no absence
        # of partners outside the admitted context.
        if type(closed) in (EnrichedContextEvidence,VerifiedContextEvidence) and (closed.closed or closed.writer is not writer):
            raise EvidenceStoreError("foreign enriched genotype context")
        _sha(source_manifest_sha256)
        _sha(reference_manifest_sha256)
        self.writer, self.connection = writer, writer.connection
        self.prefix = "genotypes_" + hashlib.sha256(closed.context_id.encode()).hexdigest()
        self.tables, self.cache, self.cache_bytes = [], OrderedDict(), 0
        self.closed = False
        self.header = None
        self.published = None
        self.context_id = closed.context_id
        self.input_evidence = None if type(closed) is ClosedContext else {"name":closed.name,"sha256":closed.sha256}
        self.sample_name, self.source_manifest_sha256 = sample_name, source_manifest_sha256
        self.reference_manifest_sha256 = reference_manifest_sha256
        try:
            for name, definition in (
                ("members", "id TEXT PRIMARY KEY,native TEXT,event TEXT,sequence INTEGER,locus INTEGER,screened INTEGER,parent INTEGER"),
                ("calls", "source_id TEXT,member_id TEXT,read_required INTEGER,PRIMARY KEY(source_id,member_id)"),
                ("sources", "id TEXT PRIMARY KEY,record_id TEXT,allele_number INTEGER,input_sha TEXT,sample_name TEXT,sample TEXT"),
                ("records", "id TEXT PRIMARY KEY,contig TEXT,input_sha TEXT,release TEXT,sample TEXT"),
                ("samples", "input_sha TEXT PRIMARY KEY,sample_name TEXT"),
                ("expected", "id TEXT PRIMARY KEY"),
                ("source_facts", "id TEXT PRIMARY KEY,content TEXT"),
                ("source_proofs", "id TEXT PRIMARY KEY,content TEXT"),
                ("allele_facts", "id TEXT PRIMARY KEY,content TEXT"),
            ):
                table = self.prefix + "_" + name
                self.connection.execute(f"CREATE TABLE {table} ({definition})")
                self.tables.append(table)
            self.connection.execute(f"CREATE INDEX {self.prefix}_events ON {self.prefix}_members(event)")
            self.connection.execute(f"CREATE INDEX {self.prefix}_member_calls ON {self.prefix}_calls(member_id,source_id)")
            self.connection.execute(f"CREATE INDEX {self.prefix}_record_inputs ON {self.prefix}_records(input_sha,id)")
            digest = hashlib.sha256(b'{"members":[')
            count = 0
            for member in closed.iter_members():
                if member.parent_source_manifest_sha256 != source_manifest_sha256 or member.allele is None:
                    raise EvidenceStoreError("genotype member source binding changed")
                allele = member.allele
                native = next(i for i in allele.identities if i.identity_id == allele.native_identity_id)
                if native.reference_manifest_sha256 != reference_manifest_sha256:
                    raise EvidenceStoreError("genotype member reference binding changed")
                # D28 hashes the explicit logical member object. It does not
                # claim the old expanded ledger-string inventory identity.
                content = _canonical(member).encode("ascii")
                if len(content) > MAX_FACT_BYTES:
                    raise EvidenceStoreError("genotype member exceeds record capacity",capacity=True)
                digest.update((b"," if count else b"") + content)
                count += 1
                event = _allele_event(allele)
                projected = longtr.is_projection(allele)
                self.connection.execute(f"INSERT INTO {self.prefix}_members VALUES (?,?,?,?,?,?,?)",
                    (allele.allele_id,allele.native_identity_id,_canonical(event),int(event[0]=="sequence"),
                     int(isinstance(native.identity,RepeatAllele) and not projected),
                     int(isinstance(native.identity,MosaicAllele)),int(longtr.is_parent(allele))))
                self.connection.executemany(f"INSERT OR IGNORE INTO {self.prefix}_expected VALUES (?)",
                    ((key,) for key in member.source_allele_ids))
                for call in allele.calls or ():
                    source_id,record_id,number,input_sha = ledger_call_binding(allele,call,
                        sample_name=sample_name,source_manifest_sha256=source_manifest_sha256)
                    sample = _canonical(call.sample)
                    if len(sample) > MAX_FACT_BYTES:
                        raise EvidenceStoreError("genotype source sample exceeds capacity",capacity=True)
                    prior = self.connection.execute(f"SELECT sample FROM {self.prefix}_records WHERE id=?",(record_id,)).fetchone()
                    if prior is not None and prior[0] != sample:
                        raise EvidenceStoreError("genotype record sample conflicts")
                    self.connection.execute(f"INSERT OR IGNORE INTO {self.prefix}_records VALUES (?,?,?,?,?)",
                        (record_id,_identity_contig(allele,source_id),input_sha,call.source.source_release,sample))
                    prior = self.connection.execute(f"SELECT sample,record_id FROM {self.prefix}_sources WHERE id=?",(source_id,)).fetchone()
                    if prior is not None and prior[0] != sample:
                        raise EvidenceStoreError("genotype source sample conflicts")
                    if prior is not None and prior[1] != record_id:
                        # The former complete-context relevance owner refused
                        # one source allele assigned to different records.
                        raise EvidenceStoreError("genotype source record conflicts")
                    self.connection.execute(f"INSERT OR IGNORE INTO {self.prefix}_sources VALUES (?,?,?,?,?,?)",
                        (source_id,record_id,number,input_sha,call.sample.sample_name,sample))
                    self.connection.execute(f"INSERT OR REPLACE INTO {self.prefix}_calls VALUES (?,?,?)",
                        (source_id,allele.allele_id,int(bool(call.source.source_fields["read_validation_required"]))))
                    proof = _canonical(call.source)
                    if len(proof) > MAX_FACT_BYTES:
                        raise EvidenceStoreError("genotype source proof exceeds capacity",capacity=True)
                    self.connection.execute(f"INSERT OR REPLACE INTO {self.prefix}_source_proofs VALUES (?,?)",(source_id,proof))
                    for key in (call.source.source_sha256,input_sha):
                        previous = self.connection.execute(f"SELECT sample_name FROM {self.prefix}_samples WHERE input_sha=?",(key,)).fetchone()
                        if previous is not None and previous[0] != call.sample.sample_name:
                            raise EvidenceStoreError("genotype input sample conflicts")
                        self.connection.execute(f"INSERT OR IGNORE INTO {self.prefix}_samples VALUES (?,?)",(key,call.sample.sample_name))
            if count != closed.member_count or self.connection.execute(
                    f"SELECT id FROM {self.prefix}_sources EXCEPT SELECT id FROM {self.prefix}_expected").fetchone() or self.connection.execute(
                    f"SELECT id FROM {self.prefix}_expected EXCEPT SELECT id FROM {self.prefix}_sources").fetchone():
                raise EvidenceStoreError("genotype source/member closure is incomplete")
            digest.update(b'],"reference_manifest_sha256":'+_canonical(reference_manifest_sha256).encode("ascii")
                +b',"schema":'+_canonical(INVENTORY_SCHEMA).encode("ascii")+b',"source_manifest_sha256":'
                +_canonical(source_manifest_sha256).encode("ascii")+b',"version":1}')
            samples, sample_bytes = [], 2
            for row in self.connection.execute(f"SELECT input_sha,sample_name FROM {self.prefix}_samples ORDER BY input_sha"):
                sample_bytes += len(_canonical(row)) + bool(samples)
                if sample_bytes > MAX_FACT_BYTES:
                    raise EvidenceStoreError("genotype sample metadata exceeds capacity",capacity=True)
                samples.append(row)
            samples = tuple(samples)
            self.header = GenotypeFactHeader(digest.hexdigest(),source_manifest_sha256,reference_manifest_sha256,
                sample_name,POLICY_SHA256,samples)
            self._source_facts()
            self._allele_facts()
        except BaseException:
            writer.failed = True
            self.close()
            raise

    def _write(self, name, key, value):
        content = _canonical(value)
        if len(content) > MAX_FACT_BYTES:
            raise EvidenceStoreError("genotype fact exceeds record capacity",capacity=True)
        self.connection.execute(f"INSERT INTO {self.prefix}_{name} VALUES (?,?)",(key,content))

    def _source_facts(self):
        for source_id,record_id,number,input_sha,sample_name in self.connection.execute(
                f"SELECT id,record_id,allele_number,input_sha,sample_name FROM {self.prefix}_sources ORDER BY id"):
            contig,record_input,_,content = self.connection.execute(
                f"SELECT contig,input_sha,release,sample FROM {self.prefix}_records WHERE id=?",(record_id,)).fetchone()
            # Match the original inventory's input->release dictionary exactly.
            release, = self.connection.execute(f"SELECT release FROM {self.prefix}_records WHERE input_sha=? ORDER BY id DESC LIMIT 1",
                (input_sha,)).fetchone()
            derived = tuple(row[0] for row in self.connection.execute(
                f"SELECT member_id FROM {self.prefix}_calls WHERE source_id=? ORDER BY member_id",(source_id,)))
            read_required, = self.connection.execute(f"SELECT max(read_required) FROM {self.prefix}_calls WHERE source_id=?",(source_id,)).fetchone()
            sample = _typed(json.loads(content),SampleCallQc,decode=True)
            if sample.sample_name != sample_name:
                raise EvidenceStoreError("genotype source/record sample conflicts")
            value = evaluate_source_genotype(original=_LedgerSourceView(record_id,source_id,number,input_sha,sample_name),
                record=SimpleNamespace(contig=contig),sample=sample,derived=derived,source_status="derived",
                release=release,read_required=bool(read_required))
            self._write("source_facts",source_id,value)

    def _load(self, name, key, cls):
        if self.closed:
            raise EvidenceStoreError("genotype fact owner is closed")
        cache_key = (name,key)
        if cache_key in self.cache:
            self.cache.move_to_end(cache_key)
            return self.cache[cache_key][0]
        row = self.connection.execute(f"SELECT content FROM {self.prefix}_{name} WHERE id=?",(key,)).fetchone()
        if row is None:
            raise EvidenceStoreError("genotype fact is missing")
        value = _typed(json.loads(row[0]),cls,decode=True)
        size = len(row[0])
        if size <= FACT_CACHE_BYTES:
            while self.cache and self.cache_bytes + size > FACT_CACHE_BYTES:
                _,(_,old_size) = self.cache.popitem(last=False)
                self.cache_bytes -= old_size
            self.cache[cache_key] = (value,size)
            self.cache_bytes += size
        return value

    def source(self, key):
        return self._load("source_facts",key,SourceGenotypeProperties)

    def source_proof(self, key):
        return self._load("source_proofs",key,SourcedRecord)

    def allele(self, key):
        return self._load("allele_facts",key,AlleleGenotypeProperties)

    def _sources(self, cursor):
        values,size = [],0
        for (key,) in cursor:
            length, = self.connection.execute(f"SELECT length(content) FROM {self.prefix}_source_facts WHERE id=?",(key,)).fetchone()
            size += length
            if size > MAX_SOURCE_WINDOW_BYTES:
                raise EvidenceStoreError("genotype contributor window exceeds capacity",capacity=True)
            values.append(self.source(key))
        return tuple(values)

    def _require_open(self):
        if self.closed:
            raise EvidenceStoreError("genotype fact owner is closed")

    def own_sources(self, key):
        self._require_open()
        return self._sources(self.connection.execute(f"SELECT source_id FROM {self.prefix}_calls WHERE member_id=? ORDER BY source_id",(key,)))

    def event_sources(self, key):
        self._require_open()
        row = self.connection.execute(f"SELECT event,sequence FROM {self.prefix}_members WHERE id=?",(key,)).fetchone()
        if row is None:
            raise EvidenceStoreError("genotype member is missing")
        event,sequence = row
        if not sequence:
            return self.own_sources(key)
        return self._sources(self.connection.execute(f"SELECT DISTINCT c.source_id FROM {self.prefix}_members m "
            f"JOIN {self.prefix}_calls c ON c.member_id=m.id WHERE m.event=? ORDER BY c.source_id",(event,)))

    def _allele_facts(self):
        for key,native,locus,screened in self.connection.execute(f"SELECT id,native,locus,screened FROM {self.prefix}_members ORDER BY id"):
            own = self.own_sources(key)
            if not own:
                continue
            value = evaluate_allele_genotype(allele_id=key,native_identity_id=native,own_calls=own,
                calls=self.event_sources(key),locus_call=bool(locus),screened=bool(screened),maximum_phase_bases=MAX_PHASE_BASES)
            self._write("allele_facts",key,value)

    def phase(self, left, right):
        self._require_open()
        if left >= right:
            raise EvidenceStoreError("genotype pair order is invalid")
        members = [self.connection.execute(f"SELECT event,locus,parent FROM {self.prefix}_members WHERE id=?",(key,)).fetchone()
                   for key in (left,right)]
        if any(row is None or row[2] for row in members) or members[0][0] == members[1][0]:
            raise EvidenceStoreError("genotype pair domain is invalid")
        return assess_pair_phase(left,right,left_properties=self.allele(left),right_properties=self.allele(right),
            left_sources=self.event_sources(left),right_sources=self.event_sources(right),
            unprojected_repeats={key for key,row in zip((left,right),members,strict=True) if row[1]},maximum_phase_bases=MAX_PHASE_BASES)

    def source_rows(self):
        self._require_open()
        for (key,) in self.connection.execute(f"SELECT id FROM {self.prefix}_source_facts ORDER BY id"):
            yield self.source(key)

    def allele_rows(self):
        self._require_open()
        for (key,) in self.connection.execute(f"SELECT id FROM {self.prefix}_allele_facts ORDER BY id"):
            yield self.allele(key)

    def publish(self, name):
        self._require_open()
        if self.published is not None:
            raise EvidenceStoreError("genotype facts are already published")
        source_digest = self.writer.add_sequence(name+".sources","genotype-source-facts-v1",self.source_rows())
        allele_digest = self.writer.add_sequence(name+".alleles","genotype-member-facts-v1",self.allele_rows())
        digest = self.writer.add(name,"genotype-facts-v1",{
            "schema":"genosapiens.wgs-context-genotype-facts","version":1,"context_id":self.context_id,"header":self.header,
            "input_evidence":self.input_evidence,
            "source_count":self.connection.execute(f"SELECT count(*) FROM {self.prefix}_source_facts").fetchone()[0],
            "member_count":self.connection.execute(f"SELECT count(*) FROM {self.prefix}_allele_facts").fetchone()[0],
            "sources_name":name+".sources","sources_sha256":source_digest,
            "alleles_name":name+".alleles","alleles_sha256":allele_digest,
            "scope":"complete-source-and-member-facts;pair-evidence-separately-required"})
        self.published = {"name":name,"kind":"genotype-facts-v1","sha256":digest}
        return digest

    def close(self):
        self.cache.clear()
        self.cache_bytes = 0
        for table in reversed(self.tables):
            self.connection.execute("DROP TABLE "+table)
        self.tables.clear()
        self.closed = True

    def __enter__(self):
        return self

    def __exit__(self,*_):
        self.close()


def verify_context_genotype_facts(reader, context, reference):
    """Replay original closed source/member facts and compare every saved row."""
    from itertools import zip_longest
    from genosapiens.wgs_context_evidence import VerifiedContextEvidence
    from genosapiens.wgs_evidence_store import EvidenceReader
    if (type(reader) is not EvidenceReader or type(context) is not VerifiedContextEvidence
            or context.writer is not reader or reference.kind != "genotype-facts-v1"):
        raise EvidenceStoreError("foreign genotype fact owner")
    value = reader.reconstruct(reference,max_bytes=MAX_FACT_BYTES)
    keys = {"schema","version","context_id","header","input_evidence","source_count","member_count",
            "sources_name","sources_sha256","alleles_name","alleles_sha256","scope"}
    if (set(value) != keys or value["schema"] != "genosapiens.wgs-context-genotype-facts"
            or type(value["version"]) is not int or value["version"] != 1
            or value["context_id"] != context.context_id
            or value["input_evidence"] != {"name":context.name,"sha256":context.sha256}
            or value["scope"] != "complete-source-and-member-facts;pair-evidence-separately-required"
            or any(type(value[k]) is not int or value[k] < 0 for k in ("source_count","member_count"))):
        raise EvidenceStoreError("invalid genotype fact closure")
    header = _typed(value["header"],GenotypeFactHeader,decode=True)
    facts = ContextGenotypeFacts(reader,context,sample_name=header.sample_name,
        source_manifest_sha256=header.input_sha256,reference_manifest_sha256=header.reference_manifest_sha256)
    try:
        if facts.header != header:
            raise EvidenceStoreError("genotype fact source header does not replay")
        for suffix,kind,cls,rows,count_key in (
            ("sources","genotype-source-facts-v1",SourceGenotypeProperties,facts.source_rows,"source_count"),
            ("alleles","genotype-member-facts-v1",AlleleGenotypeProperties,facts.allele_rows,"member_count"),
        ):
            if value[suffix+"_name"] != reference.name+"."+suffix:
                raise EvidenceStoreError("genotype fact root name is foreign")
            ref = reader.reference(value[suffix+"_name"],kind=kind)
            if ref.sha256 != value[suffix+"_sha256"]:
                raise EvidenceStoreError("genotype fact root digest changed")
            count = 0
            for expected,actual in zip_longest(rows(),reader.iter_array(ref,max_item_bytes=MAX_FACT_BYTES)):
                if expected is None or actual is None or _typed(actual,cls,decode=True) != expected:
                    raise EvidenceStoreError("genotype facts do not replay")
                count += 1
            if count != value[count_key]:
                raise EvidenceStoreError("genotype fact count changed")
        facts.published = {"name":reference.name,"kind":reference.kind,"sha256":reference.sha256}
        return facts
    except BaseException:
        facts.close()
        raise
