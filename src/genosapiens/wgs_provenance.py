"""Public source evidence helpers; persisted science uses normalized roots."""
from collections.abc import Mapping
from dataclasses import dataclass, fields, is_dataclass
import hashlib
from typing import Literal
from genosapiens.candidate_admission import AdmissionDecision, CandidateContextSeed
from genosapiens.bounded_wgs_config import bounded_wgs_execution_config_as_json
from genosapiens.evidence_model import (EvidenceReference, GeneAssertion, Observation,
    ObservationContradiction, PublicSourcedEvidence, DiseaseProfile, _sha, resolve_evidence_reference)
from genosapiens.evidence_serialization import _canonical, _plain
from genosapiens.wgs_failures import WgsFailure


class WgsProvenanceError(WgsFailure):
    """Static evidence closure failure."""

    def __init__(self, message, *, category="malformed_evidence"):
        super().__init__(
            "bounded_products", category,
            operation="bounded_product_validate", exception_kind="project_failure",
        )
        RuntimeError.__init__(self, message)



def _bytes(value):
    return (_canonical(value) + "\n").encode("ascii")



def execution_sha256(config):
    return hashlib.sha256(bounded_wgs_execution_config_as_json(config)).hexdigest()



def evidence_references(value):
    """Walk a bounded scientific value, retaining exact public identities."""
    if type(value) is EvidenceReference:
        yield value
    elif is_dataclass(value):
        for field in fields(value):
            yield from evidence_references(getattr(value, field.name))
    elif isinstance(value, Mapping):
        for child in value.values():
            yield from evidence_references(child)
    elif type(value) in (tuple, list):
        for child in value:
            yield from evidence_references(child)



def public_record_id(reference: EvidenceReference) -> str:
    return hashlib.sha256(_bytes(reference)).hexdigest()



@dataclass(frozen=True, slots=True, repr=False)
class PublicEvidenceRecord:
    reference: EvidenceReference
    owner_kind: Literal["public", "disease_profile", "disease_gene"]
    source_binding: Mapping[str, object]
    record: Mapping[str, object]
    source_fields: Mapping[str, object] | None
    assertions: tuple[GeneAssertion, ...]
    gene_store_sha256: str | None

    def __post_init__(self):
        if type(self.reference) is not EvidenceReference or self.owner_kind not in (
            "public", "disease_profile", "disease_gene",
        ):
            raise WgsProvenanceError("public provenance record is invalid")
        if (self.owner_kind != "disease_gene" and self.assertions
                or any(value.record_id != self.reference.record_id for value in self.assertions)):
            raise WgsProvenanceError("public assertion provenance is invalid")
        if (self.owner_kind == "disease_gene") != (self.gene_store_sha256 is not None):
            raise WgsProvenanceError("public store provenance is invalid")
        if self.gene_store_sha256 is not None:
            _sha(self.gene_store_sha256)



def public_records(document, values) -> tuple[PublicEvidenceRecord, ...]:
    """Project referenced public-store rows without embedding a whole store."""
    owners = {owner.evidence_id: owner for owner in (
        *document.disease_profiles, *document.disease_gene_evidence,
        *document.public_sourced_evidence,
    )}
    references = {public_record_id(ref): ref for ref in evidence_references(values)}
    # Preserve the source-equivalence component of each primary context. This
    # uses no disclosed diagnosis and cannot affect scoring, but later prior
    # proofs may cite an equivalence row absent from the primary score ledger.
    reachable = {value for hypothesis in document.hypotheses for value in (
        hypothesis.context.disease_id, hypothesis.context.native_disease_id,
    ) if value is not None}
    pending = True
    while pending:
        pending = False
        for owner in document.disease_gene_evidence:
            for relation in owner.store.relations:
                if (relation.kind != "equivalent" or not relation.subject_id
                        or not relation.object_id or not relation.record_ids
                        or not {relation.subject_id, relation.object_id} & reachable):
                    continue
                for value in (relation.subject_id, relation.object_id):
                    if value not in reachable:
                        reachable.add(value)
                        pending = True
                for record_id in relation.record_ids:
                    ref = EvidenceReference(owner.evidence_id, record_id,
                                            owner.public_manifest_sha256)
                    references[public_record_id(ref)] = ref
    # Retain all context-local profile annotations, including inheritance and
    # onset terms that are provenance but do not contribute phenotype scores.
    for owner in document.disease_profiles:
        for row in owner.store.annotations:
            if {row.source_fields.get("source_disease_id"),
                    row.source_fields.get("mapped_disease_id")} & reachable:
                ref = EvidenceReference(owner.evidence_id, row.annotation_id,
                                        owner.public_manifest_sha256)
                references[public_record_id(ref)] = ref
    # Context-local raw query evidence is retained even if it did not score.
    for owner in document.public_sourced_evidence:
        for record_id in owner.records:
            ref = EvidenceReference(owner.evidence_id, record_id,
                                    owner.public_manifest_sha256)
            references[public_record_id(ref)] = ref
    result = []
    for key in sorted(references):
        ref = references[key]
        raw = resolve_evidence_reference(document, ref)
        owner = owners[ref.evidence_id]
        if type(owner) is PublicSourcedEvidence:
            kind, binding = "public", {}
            record = owner.records[ref.record_id]
        elif type(owner) is DiseaseProfile:
            kind, binding = "disease_profile", owner.store.source_binding
            record = next(row for row in owner.store.annotations
                          if row.annotation_id == ref.record_id)
        else:
            kind, binding = "disease_gene", owner.store.source_binding
            record = owner.store.records[ref.record_id]
        assertions = (() if kind != "disease_gene" else tuple(
            value for value in owner.store.assertions if value.record_id == ref.record_id))
        result.append(PublicEvidenceRecord(ref, kind, binding, _plain(record), raw, assertions,
            owner.store.evidence_sha256 if kind == "disease_gene" else None))
    return tuple(result)



@dataclass(frozen=True, slots=True, repr=False)
class ContextProvenance:
    observations: tuple[Observation, ...]
    contradictions: tuple[ObservationContradiction, ...]
    admission_seed: CandidateContextSeed
    admission_decisions: tuple[AdmissionDecision, ...]
    public_records: tuple[PublicEvidenceRecord, ...]

    def __post_init__(self):
        if set(self.admission_seed.decision_ids) != {
            value.decision_id for value in self.admission_decisions
        }:
            raise WgsProvenanceError("context admission provenance is incomplete")
