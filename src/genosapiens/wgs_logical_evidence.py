"""Bounded logical views of the single current normalized scientific artifact.

These values are reconstructed on explicit request. They are not alternate
saved scientific formats; their authority is the complete root closure.
"""
from contextlib import closing
from dataclasses import dataclass
import hashlib
import json
from typing import Literal

from genosapiens.evidence_model import (
    AlleleEvidence, EvidenceReference as PublicReference,
    Observation, ObservationContradiction, ScoreLedger,
)
from genosapiens.evidence_fusion import HypothesisEvidence
from genosapiens.evidence_serialization import _canonical
from genosapiens.genotype_hypotheses import ConstructedHypothesis
from genosapiens.track1_submission import SubmissionAlleleTrace, _allele_trace_key
from genosapiens.wgs_case_resources import bound_reference
from genosapiens.wgs_context_evidence import VerifiedContextEvidence
from genosapiens.wgs_context_population import PopulationAssessmentValues
from genosapiens.wgs_evidence_store import EvidenceReference, EvidenceStoreError
from genosapiens.wgs_provenance import PublicEvidenceRecord, public_record_id

MAX_LOGICAL_RECORD_BYTES = 128 * 1024**2


@dataclass(frozen=True, slots=True, repr=False)
class ScientificHypothesisDetail:
    context_ordinal: int
    context_id: str
    constructed: ConstructedHypothesis
    evidence: HypothesisEvidence
    ledger: ScoreLedger
    population_assessment: PopulationAssessmentValues
    population_policy: EvidenceReference
    eligibility_reasons: tuple[str, ...]
    allele_traces: tuple[SubmissionAlleleTrace, ...]

    def __post_init__(self):
        if (type(self.context_ordinal) is not int or self.context_ordinal < 1
                or self.context_id != self.constructed.admission_context_id
                or self.evidence.hypothesis != self.constructed.hypothesis
                or self.ledger.hypothesis_id != self.constructed.hypothesis.hypothesis_id
                or type(self.population_assessment) is not PopulationAssessmentValues
                or type(self.population_policy) is not EvidenceReference
                or self.population_policy.kind != "population-policy-evidence-v1"
                or self.eligibility_reasons != tuple(dict.fromkeys(self.eligibility_reasons))
                or self.allele_traces != tuple(sorted(self.allele_traces, key=_allele_trace_key))
                or (not self.eligibility_reasons) != bool(self.allele_traces)):
            raise EvidenceStoreError("invalid normalized hypothesis detail")


def hypothesis_detail(reader, record):
    policy = reader.reader.reference(reader.resources.population.name, kind="population-policy-evidence-v1")
    return ScientificHypothesisDetail(record["context_ordinal"], record["context_id"], record["constructed"],
        record["evidence"], record["ledger"], record["population_assessment"], policy,
        record["eligibility_reasons"], record["allele_traces"])


@dataclass(frozen=True, slots=True, repr=False)
class LogicalProvenanceRecord:
    kind: Literal["allele", "observation", "contradiction", "public"]
    record_id: str
    allele: AlleleEvidence | None = None
    observation: Observation | None = None
    contradiction: ObservationContradiction | None = None
    public: PublicEvidenceRecord | None = None

    def __post_init__(self):
        fields = {"allele": self.allele, "observation": self.observation,
            "contradiction": self.contradiction, "public": self.public}
        if self.kind not in fields or tuple(k for k, v in fields.items() if v is not None) != (self.kind,):
            raise EvidenceStoreError("invalid logical provenance record")
        value = fields[self.kind]
        identifier = (value.allele_id if self.kind == "allele" else value.scoring.observation_id
            if self.kind == "observation" else value.contradiction_id if self.kind == "contradiction"
            else public_record_id(value.reference))
        if identifier != self.record_id:
            raise EvidenceStoreError("logical provenance identity differs")


def iter_logical_provenance(reader):
    """Complete source rows, deduplicated on disk; no pair/hypothesis replay."""
    reader._require_open()
    db = reader.reader.connection
    table = "logical_provenance_seen"
    with reader.reader.read_scope():
        db.execute(f"CREATE TABLE {table} (kind TEXT,id TEXT,digest TEXT,PRIMARY KEY(kind,id))")
        db.execute("CREATE TABLE logical_provenance_public_sequences (digest TEXT PRIMARY KEY)")
        try:
            def unique(record):
                encoded = _canonical(record).encode("ascii")
                if len(encoded) > MAX_LOGICAL_RECORD_BYTES:
                    raise EvidenceStoreError("logical provenance row exceeds capacity", capacity=True)
                digest = hashlib.sha256(encoded).hexdigest()
                saved = db.execute(f"SELECT digest FROM {table} WHERE kind=? AND id=?", (record.kind, record.record_id)).fetchone()
                if saved is not None:
                    if saved != (digest,):
                        raise EvidenceStoreError("logical provenance row conflicts across contexts")
                    return False
                db.execute(f"INSERT INTO {table} VALUES (?,?,?)", (record.kind, record.record_id, digest))
                return True
            def public(owner, key, kind, binding, value, fields, assertions=(), gene_sha=None):
                ref = PublicReference(owner.evidence_id, key, owner.public_manifest_sha256)
                record = PublicEvidenceRecord(ref, kind, binding, json.loads(_canonical(value)), fields, assertions, gene_sha)
                return LogicalProvenanceRecord("public", public_record_id(ref), public=record)
            for value in reader.resources.phenotypes.observations:
                row = LogicalProvenanceRecord("observation", value.scoring.observation_id, observation=value)
                if unique(row):
                    yield row
            for value in reader.resources.phenotypes.view.contradictions:
                row = LogicalProvenanceRecord("contradiction", value.contradiction_id, contradiction=value)
                if unique(row):
                    yield row
            prepared = reader.resources.genes.prepared
            for owner in prepared.evidence:
                for key in sorted(owner.store.records):
                    value = owner.store.records[key]
                    assertions = prepared.assertions_by_record.get((owner.evidence_id, key), ())
                    row = public(owner, key, "disease_gene", owner.store.source_binding, value,
                        value.source_fields, assertions, owner.store.evidence_sha256)
                    if unique(row):
                        yield row
            owner = reader.resources.phenotypes.profile
            if owner is not None:
                for value in owner.store.annotations:
                    row = public(owner, value.annotation_id, "disease_profile", owner.store.source_binding, value, value.source_fields)
                    if unique(row):
                        yield row
            for context in reader.case.contexts():
                header = reader.reader.reconstruct(bound_reference(reader.reader, context["context"], "scientific-case-context-v1"))
                with VerifiedContextEvidence(reader.reader, bound_reference(reader.reader, header["members"], "context-evidence-v1")) as members:
                    with closing(members.iter_members()) as values:
                        for member in values:
                            row = LogicalProvenanceRecord("allele", member.allele.allele_id, allele=member.allele)
                            if unique(row):
                                yield row
                    with closing(members.iter_unique_public_records('logical_provenance_public_sequences')) as values:
                        for value in values:
                            if value.record is not None:
                                row = public(value, value.record.source_record_id, 'public', {}, value.record, value.record.source_fields)
                                if unique(row):
                                    yield row
        finally:
            db.execute("DROP TABLE logical_provenance_public_sequences")
            db.execute(f"DROP TABLE {table}")
