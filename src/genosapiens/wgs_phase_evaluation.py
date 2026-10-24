"""Exact dependency caching of the existing phase evaluator.

For members with one source and no member conflict bases, source identity,
record, input and contig participate only through equality or output labels.
The key preserves those four cross-member relations and every remaining field
read by assess_pair_phase. A miss calls that kernel; a hit reuses its evaluated
slot pattern. Complex source sets keep the complete ordinary evaluation path.
"""
from collections import OrderedDict
from dataclasses import dataclass

from genosapiens.evidence_serialization import _canonical
from genosapiens.wgs_evidence_store import EvidenceStoreError
from genosapiens.wgs_phase_patterns import phase_pattern

MAX_METADATA_BYTES = 64 * 1024**2
MAX_LOCAL_CLASSES = 100_000
HOT_PHASES = 4096
MAX_CACHED_PATTERN_BYTES = 16 * 1024**2


@dataclass(frozen=True,slots=True)
class _SimpleMember:
    local_class: int
    source_id: str
    record_id: str
    input_sha256: str | None
    contig: str
    event: str


class CachedPairPhases:
    def __init__(self, facts):
        facts._require_open()
        self.facts,self.members,self.hot = facts,{},OrderedDict()
        self.evaluations,self.hits,self.metadata_bytes = 0,0,0
        self.cache_bytes = 0
        classes = {}
        for prop in facts.allele_rows():
            if len(prop.event_source_allele_ids) != 1 or prop.phase_conflict_bases:
                continue
            source = facts.source(prop.event_source_allele_ids[0])
            locus,parent,event = facts.connection.execute(f"SELECT locus,parent,event FROM {facts.prefix}_members WHERE id=?",(prop.allele_id,)).fetchone()
            if parent:
                continue
            local = (source.reported_ploidy,source.minimum_alt_copies,source.maximum_alt_copies,
                source.phase_conflict,source.allele_number,source.orientations,
                prop.conflicting_source_genotypes,bool(locus))
            if local not in classes:
                if len(classes) >= MAX_LOCAL_CLASSES:
                    raise EvidenceStoreError("phase local input class capacity exceeded",capacity=True)
                classes[local] = len(classes)
                self.metadata_bytes += len(_canonical(local))
            member = _SimpleMember(classes[local],source.source_allele_id,source.source_record_id,source.input_sha256,source.contig,event)
            self.metadata_bytes += len(_canonical((prop.allele_id,member)))
            if self.metadata_bytes > MAX_METADATA_BYTES:
                raise EvidenceStoreError("phase relation metadata exceeds capacity",capacity=True)
            self.members[prop.allele_id] = member

    def pattern(self, left, right):
        self.facts._require_open()
        if left >= right:
            raise EvidenceStoreError("phase pair order is invalid")
        a,b = self.members.get(left),self.members.get(right)
        if a is not None and b is not None and a.event == b.event:
            raise EvidenceStoreError("phase pair domain is invalid")
        key = None if a is None or b is None else (a.local_class,b.local_class,
            a.source_id == b.source_id,a.record_id == b.record_id,
            a.input_sha256 == b.input_sha256,a.contig == b.contig)
        if key is not None and key in self.hot:
            self.hits += 1
            self.hot.move_to_end(key)
            return self.hot[key][0]
        phase = self.facts.phase(left,right)
        value = phase_pattern(phase,left_sources=self.facts.allele(left).event_source_allele_ids,
            right_sources=self.facts.allele(right).event_source_allele_ids)
        self.evaluations += 1
        if key is not None:
            size = len(_canonical(value))
            if size <= MAX_CACHED_PATTERN_BYTES:
                while self.hot and (len(self.hot) >= HOT_PHASES or self.cache_bytes+size > MAX_CACHED_PATTERN_BYTES):
                    _,(_,old_size) = self.hot.popitem(last=False)
                    self.cache_bytes -= old_size
                self.hot[key] = (value,size)
                self.cache_bytes += size
        return value
