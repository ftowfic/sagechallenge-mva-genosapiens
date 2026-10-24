"""Decision 6a: offline, diagnosis-blind context admission, BEFORE scoring.

Select original inputs and independent resource pins, then select the SHA-256 of
``admission_policy_bytes``. Pass that pin to ``generate_candidate_admission`` and
to ``validate_candidate_admission`` with the same independently held inputs.
The policy includes actual installed executable bytes and source bindings. It is
a development reproducibility pin, NOT task 9.6 acceptance or task 13.1 approval.

The public boundary accepts no association proofs, predicates, shortlist, scores,
diagnosis, narrative, acquisition function, model, or read executor. Annotation
is replayed from a complete publication and the ORIGINAL ORDERED VepRequests.
ENCODE, Open Targets and ClinVar use only their fixed offline file adapters.
The internal ``_evaluate`` kernel performs no I/O and accepts no callbacks.

None explicitly declares an unavailable source; it does not disable a rule.
Every required source/subject evaluation then remains unresolved. A declared
publication must include every derived allele, including native repeats. Present
link resources require their original reference files. Repeat catalogs retain
original bytes in DiseaseGeneEvidence and require independent contract/manifest
pins here; a missing repeat reference is an explicit unresolved assessment.
Profiles/genes/IC are the complete existing typed stores, not selected contexts;
their raw inputs are policy-bound and their interpretations are rebuilt.

Results contain HypothesisContext seeds, never GenotypeHypothesis or a partner
search certificate. The complete inventory (including nonadmitted alleles and
unresolved raw ALT slots) must reach tasks 9.2/9.3. Source technical obligations
are unchanged. Full audit bytes bind provenance; scientific_bytes deliberately
exclude physical row locators, run hashes and optional typed handoff enrichment.
"""

from __future__ import annotations

from collections import OrderedDict, defaultdict, deque
from collections.abc import Mapping
from contextlib import ExitStack, closing, contextmanager
from contextvars import ContextVar
from dataclasses import dataclass, fields, replace
import hashlib
import json
from pathlib import Path
import re
from types import MappingProxyType
from typing import Literal

from genosapiens import admission_inventory as inv, admission_links as links
from genosapiens import annotation_admission as annotation, repeat_catalog_admission as repeat
from genosapiens.disease_gene_evidence import (
    GeneRequirements, _requirements_disagree, assertion_native_disease_id,
    rebuild_disease_gene_evidence,
)
from genosapiens.disease_profiles import build_disease_profiles, disease_profiles_as_bytes
from genosapiens.evidence_model import (
    AlleleEvidence, AlleleReference, DiseaseGeneEvidence, DiseaseProfile,
    EvidenceReference, HypothesisContext, LocusBinding, PublicSourcedEvidence,
    ScoringView, SourcedRecord, _Controlled, _typed,
)
from genosapiens.evidence_serialization import _plain, from_json
from genosapiens.hpo_ic import HpoICStore, TRANSFORM, build_hpo_graph, build_hpo_ic, hpo_ic_as_bytes
from genosapiens.rare_disease import _normalise_curie
from genosapiens.vep_annotation import VepRequest
from genosapiens.wgs_admission import LedgerAdmissionInput
from genosapiens.allele_identity import ReferenceRole


SCHEMA = "genosapiens.candidate-admission"
EXECUTABLE_VERSION = "decision-6a-development"
Outcome = Literal["admitted", "evaluated_no_match", "unresolved"]
Support = Literal["present", "absent", "unknown"]
_SHA = re.compile(r"[0-9a-f]{64}")
_ENSG = re.compile(r"ENSG[0-9]{11}(?:\.[1-9][0-9]*)?")
_RULE_ADAPTERS = {
    inv.RULE_IDS[0]: ("annotation",), inv.RULE_IDS[1]: ("annotation",),
    inv.RULE_IDS[2]: ("encode",), inv.RULE_IDS[3]: ("repeat",),
    inv.RULE_IDS[4]: ("open_targets", "clinvar", "curated"),
}
_EXECUTABLE_FILES = (
    "candidate_admission.py", "wgs_admission_graph.py", "wgs_context_admission.py", "wgs_context_payloads.py", "admission_inventory.py", "annotation_admission.py", "longtr_admission.py",
    "wgs_association_projection.py", "wgs_association_codec.py", "wgs_admission_projection.py", "wgs_admission_inputs.py", "wgs_context_pipeline.py", "track1_association_scope.py",
    "wgs_analysis_eligibility.py", "wgs_allele_query_cache.py", "wgs_validation_audit.py",
    "admission_links.py", "repeat_catalog_admission.py", "disease_gene_evidence.py",
    "disease_profiles.py", "hpo_ic.py", "ontology_bundle.py", "evidence_model.py",
    "evidence_serialization.py", "allele_identity.py", "variant_inventory.py",
    "variant_normalization.py", "vcf_qc.py", "vep_annotation.py", "vep_native_annotation.py",
    "vep_bundle.py", "variant_class_resources.py", "open_targets.py", "clinvar.py",
    "reference_concordance.py", "grch38_reference.py", "submission_rendering.py", "manifest.py", "rare_disease.py",
)


def _bytes(value):
    return (json.dumps(_plain(value), sort_keys=True, separators=(",", ":"),
                       ensure_ascii=True, allow_nan=False) + "\n").encode("ascii")


def _hash(value):
    return hashlib.sha256(_bytes(value)).hexdigest()


_RECORD_BYTES = ContextVar('admission_record_bytes', default=None)
_RECORD_CACHE_BYTES = 8 * 1024**2
_RECORD_CACHE_ENTRIES = 8192


@contextmanager
def _record_byte_scope():
    """Only live worker-owned immutable records; never external serialization."""
    owner = [OrderedDict(), 0]
    token = _RECORD_BYTES.set(owner)
    try:
        yield owner
    finally:
        owner[0].clear()
        owner[1] = 0
        _RECORD_BYTES.reset(token)


def _record_bytes(value):
    owner = _RECORD_BYTES.get()
    if owner is None or not isinstance(value, _Controlled):
        return _bytes(value)
    cache, size = owner
    found = cache.get(id(value))
    if found is not None and found[0] is value:
        cache.move_to_end(id(value))
        return found[1]
    content = _bytes(value)
    if len(content) <= _RECORD_CACHE_BYTES:
        while cache and (size + len(content) > _RECORD_CACHE_BYTES
                         or len(cache) >= _RECORD_CACHE_ENTRIES):
            _, (_, previous) = cache.popitem(last=False)
            size -= len(previous)
        cache[id(value)] = (value, content)
        owner[1] = size + len(content)
    return content


POLICY_BYTES = _bytes({
    "schema": SCHEMA, "executable_version": EXECUTABLE_VERSION,
    "rules": _RULE_ADAPTERS, "inventory": inv.INVENTORY_POLICY_SHA256,
    "annotation": annotation.POLICY_SHA256, "links": links._POLICY,
    "repeat": repeat.POLICY, "synthetic_repeat": repeat.SYNTHETIC_ADAPTER_VERSION,
    "ic": TRANSFORM,
    "phenotype": "all-scopes; resolved-present x positive-P; Resnik/MICA > 0; no-top-K",
    "phenotype_unknown": "missing-store/profile/positive-observation-or-unresolved-comparison; not-negative",
    "disease": "source-native-anchor; explicit-equivalence-and-child-to-parent-proofs-only; mapped-view-not-equivalence",
    "clinical_condition": "existing-source-CURIE-spelling-normalization-for-join; retain-original-condition-and-source; no-concept-mapping",
    "target": "original-source-gene; Ensembl-stable-ID-with-version-retained; no-symbol/HGNC-crosswalk",
    "support": "rebuild-source-polarity-positive-or-contextual; no-numerical-aggregation",
    "specificity": "compatible-known-MOI-direct-curation-supersedes-umbrella; retain-equal-conflicts-unpenalized",
    "union": "native-disease,target-namespace,raw-gene/locus,exact-semantic-requirement-alternative",
    "unknown_requirements": "distinct-from-known; no-transitive-compatibility-merge",
    "dedup": "one-context; source-family/shared-PMID-dependence-groups; no-admission-score",
    "unavailable": "explicit-None-unresolved; present-invalid-or-incomplete-fails",
    "coverage": "all-raw-records/ALT-slots; all-required-source-evaluations; no-partner-certificate",
    "validation": "development-only; 9.6-and-13.1-pending",
})
POLICY_SHA256 = hashlib.sha256(POLICY_BYTES).hexdigest()


class CandidateAdmissionError(ValueError):
    """Static diagnostics, never protected input values or paths."""


_VERIFIED_HPO_STORES: dict[int, HpoICStore] = {}
_VERIFIED_PROFILE_STORES: dict[
    tuple[int, int], tuple[DiseaseProfile, HpoICStore]
] = {}
_VERIFIED_GENE_STORES: dict[int, object] = {}


def _remember_verified(cache: dict, key: object, value: object) -> None:
    if key not in cache and len(cache) >= 16:
        cache.pop(next(iter(cache)))
    cache[key] = value


def _verify_hpo_store(store: HpoICStore, authenticated_resources=None) -> None:
    if authenticated_resources is not None:
        from genosapiens.track1_sensitivity import _authenticated_ic_hash
        _authenticated_ic_hash(authenticated_resources, store)
        return
    cached = _VERIFIED_HPO_STORES.get(id(store))
    if cached is store:
        return
    graph = build_hpo_graph(store.graph.ontology)
    rebuilt = build_hpo_ic(
        graph,
        (dict(item.fields) for item in store.annotations),
        source_binding=store.source_binding,
    )
    if rebuilt != store or not rebuilt.ready:
        raise ValueError
    # Keep a strong reference so the object id cannot be recycled while cached.
    _remember_verified(_VERIFIED_HPO_STORES, id(store), store)


def _verify_profile_store(
    profiles: DiseaseProfile, store: HpoICStore
) -> None:
    key = id(profiles), id(store)
    cached = _VERIFIED_PROFILE_STORES.get(key)
    if cached is not None and cached[0] is profiles and cached[1] is store:
        return
    rebuilt = build_disease_profiles(
        store.graph,
        (item.source_fields for item in profiles.store.annotations),
        source_binding=profiles.store.source_binding,
    )
    if rebuilt != profiles.store:
        raise ValueError
    if not (
        store.source_binding.get("synthetic") is True
        and profiles.store.source_binding.get("synthetic") is True
    ):
        roles = {item["role"]: item for item in store.source_binding["sources"]}
        if len(roles) != 2 or set(roles) != {
            "hpo_ontology",
            "hpo_disease_annotations",
        }:
            raise ValueError
        for name, profile_name in (
            ("hpo_ontology", "hpo_ontology"),
            ("hpo_disease_annotations", "hpo_annotations"),
        ):
            if any(
                roles[name][field]
                != profiles.store.source_binding[profile_name][field]
                for field in ("role", "sha256", "release")
            ):
                raise ValueError
    _remember_verified(_VERIFIED_PROFILE_STORES, key, (profiles, store))


def _verify_gene_evidence(value: DiseaseGeneEvidence) -> None:
    store = value.store
    cached = _VERIFIED_GENE_STORES.get(id(store))
    if cached is store:
        return
    if rebuild_disease_gene_evidence(store) != store:
        raise ValueError
    _remember_verified(_VERIFIED_GENE_STORES, id(store), store)


@dataclass(frozen=True, slots=True, repr=False)
class AnnotationInput:
    artifacts: Mapping[str, bytes]
    requests: tuple[VepRequest, ...]
    expected_binding_sha256: str

    def __post_init__(self):
        if (not isinstance(self.artifacts, Mapping)
                or any(type(k) is not str or type(v) is not bytes for k, v in self.artifacts.items())
                or type(self.requests) is not tuple or any(type(r) is not VepRequest for r in self.requests)
                or type(self.expected_binding_sha256) is not str or not _SHA.fullmatch(self.expected_binding_sha256)):
            raise CandidateAdmissionError("invalid annotation input structure")
        object.__setattr__(self, "artifacts", MappingProxyType(dict(self.artifacts)))


@dataclass(frozen=True, slots=True, repr=False)
class RepeatCatalogInput:
    catalog: DiseaseGeneEvidence
    expected_contract_sha256: str
    expected_manifest_sha256: str
    reference: repeat.RepeatReferenceInput | None


@dataclass(frozen=True, slots=True, repr=False)
class AdmissionSources:
    """Closed resource selection; explicitly supply None for unavailable sources."""

    annotation: AnnotationInput | None
    encode: tuple[links.PinnedPublicSource, links.PinnedPublicSource] | None
    open_targets: links.PinnedPublicSource | None
    clinvar: links.PinnedPublicSource | None
    reference: links.PinnedPublicSource | None
    repeat_catalog: RepeatCatalogInput | None


@dataclass(frozen=True, slots=True, repr=False)
class ShardedAdmissionSources:
    """Closed resource selection for already verified annotation shards."""

    annotation_binding_sha256: str | None
    encode: tuple[links.PinnedPublicSource, links.PinnedPublicSource] | None
    open_targets: links.PinnedPublicSource | None
    clinvar: links.PinnedPublicSource | None
    reference: links.PinnedPublicSource | None
    repeat_catalog: RepeatCatalogInput | None
    regulatory_index_manifest_sha256: str | None = None

    def __post_init__(self):
        if (
            self.annotation_binding_sha256 is not None
            and (
                type(self.annotation_binding_sha256) is not str
                or _SHA.fullmatch(self.annotation_binding_sha256) is None
            )
        ):
            raise CandidateAdmissionError("invalid sharded annotation binding")
        if self.regulatory_index_manifest_sha256 is not None and (
            type(self.regulatory_index_manifest_sha256) is not str
            or _SHA.fullmatch(self.regulatory_index_manifest_sha256) is None
            or type(self.encode) is not tuple or len(self.encode) != 2
            or any(type(pin) is not links.PinnedPublicSource for pin in self.encode)
            or type(self.reference) is not links.PinnedPublicSource
        ):
            raise CandidateAdmissionError("invalid sharded regulatory index binding")


@dataclass(frozen=True, slots=True, repr=False)
class ShardedRegulatoryAssociations:
    """Complete source-labelled index reports for a bounded admission window."""

    index_manifest_sha256: str
    source_manifest_sha256s: tuple[str, str]
    reference_manifest_sha256: str
    proofs: tuple[links.AssociationProof, ...]

    def __post_init__(self):
        if (type(self.source_manifest_sha256s) is not tuple
                or len(self.source_manifest_sha256s) != 2
                or any(type(value) is not str or _SHA.fullmatch(value) is None
                    for value in (self.index_manifest_sha256, self.reference_manifest_sha256,
                                  *self.source_manifest_sha256s))
                or type(self.proofs) is not tuple
                or any(type(proof) is not links.AssociationProof
                       or proof.rule != 'regulatory_link' for proof in self.proofs)):
            raise CandidateAdmissionError("invalid sharded regulatory association binding")


@dataclass(frozen=True, slots=True, repr=False)
class ShardedClinicalProof:
    """One source-labelled proof emitted by a verified stage-scoped session."""

    adapter: Literal["open_targets", "clinvar"]
    proof: links.AssociationProof

    def __post_init__(self):
        if (
            self.adapter not in {"open_targets", "clinvar"}
            or type(self.proof) is not links.AssociationProof
            or self.proof.rule != "explicit_clinical_link"
        ):
            raise CandidateAdmissionError("invalid sharded clinical proof")


@dataclass(frozen=True, slots=True, repr=False)
class ShardedClinicalAssociations:
    """Authenticated bounded clinical reports for one admission shard."""

    index_manifest_sha256: str
    open_targets_manifest_sha256: str
    clinvar_manifest_sha256: str
    reference_manifest_sha256: str
    proofs: tuple[ShardedClinicalProof, ...]

    def __post_init__(self):
        if (
            any(_SHA.fullmatch(value) is None for value in (
                self.index_manifest_sha256,
                self.open_targets_manifest_sha256,
                self.clinvar_manifest_sha256,
                self.reference_manifest_sha256,
            ))
            or type(self.proofs) is not tuple
            or any(type(value) is not ShardedClinicalProof for value in self.proofs)
        ):
            raise CandidateAdmissionError(
                "invalid sharded clinical association binding"
            )


@dataclass(frozen=True, slots=True, repr=False)
class DiseaseSupportEntry(_Controlled):
    source_ref: EvidenceReference
    polarity: str
    applicability: str
    scope: str
    requirement_state: str
    positive: bool
    relation_refs: tuple[EvidenceReference, ...]
    reasons: tuple[str, ...]
    conflict_refs: tuple[EvidenceReference, ...]
    dependency_group: str


@dataclass(frozen=True, slots=True, repr=False)
class PhenotypeAdmission(_Controlled):
    proof_id: str
    native_disease_id: str
    observation_term: str
    profile_term: str
    mica_ids: tuple[str, ...]
    information_content: float
    observation_ids: tuple[str, ...]
    contradiction_ids: tuple[str, ...]
    evidence_refs: tuple[EvidenceReference, ...]
    outcome: Outcome
    context_ids: tuple[str, ...]
    reasons: tuple[str, ...]


@dataclass(frozen=True, slots=True, repr=False)
class AdmissionDecision(_Controlled):
    """One source outcome and optional disease context, not a trusted pass flag.

    decision_id is scientific: repeated physical source rows can share it. All
    those rows remain separately represented in this audit and the raw proofs.
    binding_sha256 addresses the complete input/policy/executable binding below.
    """

    decision_id: str
    subject_id: str
    subject_scientific_id: str
    rule_id: str
    adapter: str
    association_state: Literal["match", "no_match", "unresolved"]
    disease_support: Support
    phenotype_support: Support
    outcome: Outcome
    target_namespace: str | None
    target_id: str | None
    condition_id: str | None
    evaluated_context_id: str | None
    context_ids: tuple[str, ...]
    reasons: tuple[str, ...]
    evidence_refs: tuple[EvidenceReference, ...]
    source_refs: tuple[SourcedRecord, ...]
    binding_sha256: str
    rule_configuration_sha256: str
    support_ledger: tuple[DiseaseSupportEntry, ...]


@dataclass(frozen=True, slots=True, repr=False)
class CandidateContextSeed(_Controlled):
    context: HypothesisContext
    target_namespace: str
    routes: tuple[str, ...]
    allele_refs: tuple[AlleleReference, ...]
    phenotype_support: Support
    phenotype_proof_ids: tuple[str, ...]
    decision_ids: tuple[str, ...]
    support_ledger: tuple[DiseaseSupportEntry, ...]


@dataclass(frozen=True, slots=True, repr=False)
class CandidateAdmission(_Controlled):
    inventory: inv.AdmissionInventory
    view: ScoringView
    contexts: tuple[CandidateContextSeed, ...]
    decisions: tuple[AdmissionDecision, ...]
    phenotype_proofs: tuple[PhenotypeAdmission, ...]
    phenotype_limits: tuple[str, ...]
    disease_gene_evidence: tuple[DiseaseGeneEvidence, ...]
    profiles: DiseaseProfile | None
    public_sources: tuple[PublicSourcedEvidence, ...]
    annotation_proofs: tuple[annotation.AnnotationAssociationProof, ...]
    link_proofs: tuple[links.AssociationProof, ...]
    repeat_proofs: tuple[repeat.RepeatAssociationProof, ...]
    binding: Mapping[str, object]

    def as_bytes(self) -> bytes:
        return _bytes({"schema": SCHEMA, "classification": "controlled", "llm_policy": "no-llm", "result": self})

    def scientific_bytes(self) -> bytes:
        """Deduplicated unscored science; NOT a substitute for full audit validation."""
        decisions = {d.decision_id: {key: getattr(d, key) for key in (
            "decision_id", "subject_scientific_id", "rule_id", "adapter", "association_state",
            "disease_support", "phenotype_support", "outcome", "target_namespace", "target_id",
            "condition_id", "evaluated_context_id", "context_ids", "reasons", "rule_configuration_sha256",
        )} for d in self.decisions}
        return _bytes({
            "schema": SCHEMA, "policy_rules_sha256": POLICY_SHA256,
            "contexts": [{"context_id": s.context.context_id, "native_disease_id": s.context.native_disease_id,
                "target_namespace": s.target_namespace, "gene_id": s.context.gene_id, "locus_id": s.context.locus_id,
                "requirement": _requirement_key(s.context.requirements[0]), "routes": s.routes,
                "phenotype_support": s.phenotype_support, "decision_ids": s.decision_ids,
                "phenotype_proof_ids": s.phenotype_proof_ids} for s in self.contexts],
            "decisions": [decisions[key] for key in sorted(decisions)],
            "phenotype_proofs": [{key: getattr(p, key) for key in (
                "proof_id", "native_disease_id", "observation_term", "profile_term", "mica_ids",
                "information_content", "outcome", "context_ids", "reasons",
            )} for p in self.phenotype_proofs], "phenotype_limits": self.phenotype_limits,
        })


@dataclass(frozen=True, slots=True, repr=False)
class ShardedCandidateAdmission(_Controlled):
    """One bounded batch; global context membership is closed by task 7."""

    inputs: tuple[LedgerAdmissionInput, ...]
    contexts: tuple[CandidateContextSeed, ...]
    decisions: tuple[AdmissionDecision, ...]
    phenotype_proofs: tuple[PhenotypeAdmission, ...]
    phenotype_limits: tuple[str, ...]
    public_sources: tuple[PublicSourcedEvidence, ...]
    annotation_proofs: tuple[annotation.AnnotationAssociationProof, ...]
    link_proofs: tuple[links.AssociationProof, ...]
    repeat_proofs: tuple[repeat.RepeatAssociationProof, ...]
    binding: Mapping[str, object]

    def as_bytes(self) -> bytes:
        return _bytes({
            "schema": "genosapiens.sharded-candidate-admission",
            "classification": "controlled",
            "llm_policy": "no-llm",
            "result": self,
        })


def _pin(pin):
    if pin is None:
        return None
    if type(pin) is not links.PinnedPublicSource or any(
        type(value) is not str or not _SHA.fullmatch(value) for value in (pin.contract_sha256, pin.manifest_sha256)
    ) or any(not isinstance(path, Path) for path in (pin.contract_path, pin.manifest_path, pin.root)):
        raise CandidateAdmissionError("invalid public source pin")
    return {"contract_sha256": pin.contract_sha256, "manifest_sha256": pin.manifest_sha256}


def _source_bindings(sources):
    if type(sources) is not AdmissionSources:
        raise CandidateAdmissionError("closed admission sources required")
    a, r = sources.annotation, sources.repeat_catalog
    if a is not None and type(a) is not AnnotationInput or r is not None and type(r) is not RepeatCatalogInput:
        raise CandidateAdmissionError("invalid source adapter input")
    if sources.encode is not None and (type(sources.encode) is not tuple or len(sources.encode) != 2):
        raise CandidateAdmissionError("both ENCODE source pins required")
    result = {"annotation": a.expected_binding_sha256 if a else None,
        "encode": [_pin(p) for p in sources.encode] if sources.encode is not None else None,
        "open_targets": _pin(sources.open_targets), "clinvar": _pin(sources.clinvar),
        "reference": _pin(sources.reference), "repeat": None}
    if r is not None:
        if (type(r.catalog) is not DiseaseGeneEvidence or any(type(v) is not str or not _SHA.fullmatch(v)
                for v in (r.expected_contract_sha256, r.expected_manifest_sha256))
                or r.reference is not None and type(r.reference) is not repeat.RepeatReferenceInput):
            raise CandidateAdmissionError("invalid repeat source input")
        result["repeat"] = {"contract_sha256": r.expected_contract_sha256,
            "manifest_sha256": r.expected_manifest_sha256, "catalog": _hash(r.catalog.store),
            "reference": None if r.reference is None else {
                "contract_sha256": r.reference.expected_contract_sha256,
                "manifest_sha256": r.reference.expected_manifest_sha256,
                "bytes_sha256": _hash([hashlib.sha256(b).hexdigest() for b in (
                    r.reference.contract_bytes, r.reference.manifest_bytes, r.reference.aliases_bytes)])}}
    return result


def admission_policy_bytes(*, sources: AdmissionSources, profiles: DiseaseProfile | None,
                           genes: tuple[DiseaseGeneEvidence, ...] | None, ic_store: HpoICStore | None) -> bytes:
    """Offline pin selection, not validation or acceptance. No arbitrary opcodes.

    Select this digest independently of any emitted result. Immutable typed
    profile/gene/IC raw stores are permitted inputs; source interpretations are
    replayed by generation. Aggregation/score ledgers are intentionally not read.
    Executable hashes are actual fixed installed module files, not caller labels.
    """
    try:
        bindings = _source_bindings(sources)
        if (profiles is not None and type(profiles) is not DiseaseProfile
                or ic_store is not None and type(ic_store) is not HpoICStore
                or genes is not None and (type(genes) is not tuple or any(type(g) is not DiseaseGeneEvidence for g in genes))):
            raise ValueError
        return _bytes({"rules": json.loads(POLICY_BYTES), "sources": bindings,
            "profiles": None if profiles is None else [profiles.evidence_id, profiles.public_manifest_sha256,
                hashlib.sha256(disease_profiles_as_bytes(profiles.store)).hexdigest()],
            "genes": None if genes is None else sorted(
                ([g.evidence_id, g.public_manifest_sha256, _hash(g.store)] for g in genes), key=_bytes),
            "ic": None if ic_store is None else hashlib.sha256(hpo_ic_as_bytes(ic_store)).hexdigest(),
            "executables": {name: hashlib.sha256(Path(__file__).with_name(name).read_bytes()).hexdigest()
                            for name in _EXECUTABLE_FILES}})
    except (ValueError, TypeError, RuntimeError, OSError, KeyError, AttributeError):
        raise CandidateAdmissionError("invalid admission policy inputs") from None


def sharded_admission_policy_bytes(
    *,
    sources: ShardedAdmissionSources,
    profiles: DiseaseProfile | None,
    genes: tuple[DiseaseGeneEvidence, ...] | None,
    ic_store: HpoICStore | None,
    clinical_index_manifest_sha256: str | None = None,
) -> bytes:
    """Bind admission kernels and the current screened domain to closed sources."""

    try:
        if type(sources) is not ShardedAdmissionSources:
            raise ValueError
        resource_sources = AdmissionSources(
            None, sources.encode, sources.open_targets, sources.clinvar,
            sources.reference, sources.repeat_catalog,
        )
        value = json.loads(admission_policy_bytes(
            sources=resource_sources, profiles=profiles, genes=genes, ic_store=ic_store,
        ))
        from genosapiens.wgs_analysis_eligibility import METHOD_BYTES, METHOD_SHA256
        value['analysis_method'] = json.loads(METHOD_BYTES)
        value['analysis_method_sha256'] = METHOD_SHA256
        value["sources"]["annotation"] = sources.annotation_binding_sha256
        value["annotation_transport"] = "closed-bounded-shards"
        if sources.regulatory_index_manifest_sha256 is not None:
            value['regulatory_transport'] = {
                'kind': 'authenticated-public-index-bounded-report',
                'index_manifest_sha256': sources.regulatory_index_manifest_sha256,
                'raw_source_rows': 'replayed-by-pinned-locator',
            }
        if clinical_index_manifest_sha256 is not None:
            if _SHA.fullmatch(clinical_index_manifest_sha256) is None:
                raise ValueError
            value["clinical_transport"] = {
                "kind": "authenticated-public-index-stage-session",
                "index_manifest_sha256": clinical_index_manifest_sha256,
                "raw_source_rows": "replayed-by-pinned-locator",
            }
        return _bytes(value)
    except (ValueError, TypeError, RuntimeError, OSError, KeyError, AttributeError):
        raise CandidateAdmissionError("invalid sharded admission policy inputs") from None


def _native(allele):
    return next(i for i in allele.identities if i.identity_id == allele.native_identity_id)


def _requirement_key(requirement):
    return tuple(getattr(requirement, f.name) if f.name != "allelic_requirement" else (
        None if requirement.allelic_requirement_kind == "mode_of_inheritance" else requirement.allelic_requirement
    ) for f in fields(GeneRequirements) if f.name not in ("record_ids", "diagnostics"))


def _ref(public, key):
    return EvidenceReference(public.evidence_id, key, public.public_manifest_sha256)


def _unique(items):
    return tuple(value for _, value in sorted({_record_bytes(item): item for item in items}.items()))


def _target(public, assertion):
    if assertion.locus_id is not None:
        return assertion.source_family + ":locus", assertion.locus_id
    gene = assertion.gene_id
    if gene is None:
        return None, None
    if _ENSG.fullmatch(gene):
        return "Ensembl", gene
    # A native CURIE can be retained for phenotype seeds, never guessed from a symbol.
    if re.fullmatch(r"[A-Za-z][A-Za-z0-9_]*:[A-Za-z0-9_.-]+", gene):
        return gene.split(":", 1)[0], gene
    row = public.store.records[assertion.record_id].source_fields
    return "source:" + row["source_dataset"] + ":gene_id", gene


def _gene_join(namespace, left, right):
    if left == right:
        return True
    return bool(namespace == "Ensembl" and left and right and _ENSG.fullmatch(left)
                and _ENSG.fullmatch(right) and left.split(".")[0] == right.split(".")[0])


def _relation_edges(evidence):
    edges = defaultdict(list)
    for public in evidence:
        for relation in public.store.relations:
            a, b = relation.subject_id, relation.object_id
            if not a or not b or relation.kind not in ("equivalent", "narrower", "broader"):
                continue
            refs = tuple(_ref(public, key) for key in relation.record_ids)
            if relation.kind in ("equivalent", "narrower"):
                edges[a].append((b, int(relation.kind != "equivalent"), refs))
            if relation.kind in ("equivalent", "broader"):
                edges[b].append((a, int(relation.kind != "equivalent"), refs))
    return edges


def _reach(anchor, edges):
    """Typed child-to-parent reachability, keeping all shortest-path proof refs."""
    # Encode immutable proof references once per visited edge. Equal-distance
    # paths can revisit an edge many times while their complete proofs converge.
    # Byte sets preserve _unique's exact identity and final ordering without
    # repeatedly serializing every accumulated reference on every revisit.
    records, ordered_edges = {}, {}
    reached, queue = {anchor: (0, frozenset())}, deque([anchor])
    while queue:
        a = queue.popleft()
        distance, refs = reached[a]
        if a not in ordered_edges:
            prepared = []
            for b, step, proof in sorted(edges[a], key=_bytes):
                keys = set()
                for record in proof:
                    key = _record_bytes(record)
                    records[key] = record
                    keys.add(key)
                prepared.append((b, step, frozenset(keys)))
            ordered_edges[a] = prepared
        for b, step, proof in ordered_edges[a]:
            candidate = (distance + step, refs | proof)
            if b not in reached or candidate[0] < reached[b][0]:
                reached[b] = candidate
                queue.append(b)
            elif candidate[0] == reached[b][0]:
                combined = reached[b][1] | candidate[1]
                if combined != reached[b][1]:
                    reached[b] = (candidate[0], combined)
                    queue.append(b)
    return {node: (distance, tuple(records[key] for key in sorted(refs)))
            for node, (distance, refs) in reached.items()}


def _support(context, namespace, evidence, reached):
    """Task 7.4 categorical applicability only: no weights, caps or penalties."""
    entries, details = [], []
    selected = context.requirements[0]
    for public, assertions in evidence:
        for assertion in assertions:
            ns, target = _target(public, assertion)
            if ns != namespace or target != (context.locus_id or context.gene_id):
                continue
            native = assertion_native_disease_id(public.store, assertion)
            distance, relation_refs = reached.get(native, (None, ()))
            in_view = assertion.record_id in public.store.views.get(context.disease_id, ())
            if distance is None and not in_view:
                continue
            row = public.store.records[assertion.record_id].source_fields
            compatible = tuple(r for r in assertion.requirements if not _requirements_disagree(r, selected))
            matching = any(r.inheritance is not None and r.inheritance == selected.inheritance for r in compatible)
            # Unknown requirements do not contribute to a different known alternative.
            exact_model = any(_requirement_key(r) == _requirement_key(selected) for r in compatible)
            applicability = "unknown" if distance is None else "exact" if native == context.native_disease_id else (
                "equivalent" if distance == 0 else "broader")
            scope = "native_locus" if assertion.locus_id is not None else "indirect" if row.get("is_ontology_propagated") is True else (
                "broad" if distance and distance > 0 or row.get("disease_specificity") == "broad" else
                "direct" if row.get("association_scope") in ("source_direct", "ot_mapped_evidence")
                and row.get("is_ontology_propagated") is False and row.get("disease_specificity") in ("subtype", "disorder") else "contextual")
            direct = scope == "direct" and distance == 0 and assertion.evidence_kind == "curated"
            reasons = list(assertion.diagnostics)
            if not compatible:
                reasons.append("different_inheritance_or_requirements")
            elif not exact_model:
                reasons.append("separate_requirement_alternative")
            if distance is None:
                reasons.append("mapped_view_not_equivalence_applicability_unknown")
            entries.append(DiseaseSupportEntry(_ref(public, assertion.record_id), assertion.polarity, applicability,
                scope, "different" if not compatible else "exact_alternative" if exact_model else "unknown_or_other_alternative",
                distance is not None and exact_model and assertion.polarity in ("positive", "contextual"),
                relation_refs, tuple(reasons), (), ""))
            details.append((assertion, row, distance, direct, matching, compatible))
    specific = [i for i, (a, _, _, direct, matching, _) in enumerate(details)
                if direct and matching and a.polarity in ("positive", "negative")]
    for i, entry in enumerate(entries):
        a, row, distance, _, matching, compatible = details[i]
        reasons = list(entry.reasons)
        positive = entry.positive
        if specific and distance is not None and distance > 0:
            positive = False
            reasons.append("umbrella_evidence_superseded")
            if any(details[j][0].polarity != a.polarity for j in specific):
                reasons.append("umbrella_validity_conflict")
            if not compatible:
                reasons.append("umbrella_inheritance_conflict")
        conflicts = []
        for j, peer in enumerate(entries):
            b, other_row, other_distance, _, other_matching, other_requirements = details[j]
            if ({a.polarity, b.polarity} == {"positive", "negative"}
                    and distance is not None and distance == other_distance):
                known_scope = all(assertion.evidence_kind == "curated"
                    and source.get("is_ontology_propagated") is False
                    and source.get("association_scope") in ("source_direct", "ot_mapped_evidence")
                    and source.get("disease_specificity") in ("subtype", "disorder", "broad")
                    for assertion, source in ((a, row), (b, other_row)))
                if not compatible or not other_requirements:
                    reasons.append("opposing_assertion_has_incompatible_requirements")
                elif not known_scope or not matching or not other_matching:
                    reasons.append("opposing_assertion_scope_or_requirements_unknown")
                elif any(not _requirements_disagree(left, right) for left in compatible for right in other_requirements):
                    conflicts.append(peer.source_ref)
                    reasons.append("equal_specificity_validity_conflict_unpenalized")
                else:
                    reasons.append("opposing_assertion_has_incompatible_requirements")
        if a.evidence_kind == "aggregate" and any(
            j != i and details[j][0].source_family == a.source_family
            and details[j][0].evidence_kind == "curated" and entries[j].applicability in ("exact", "equivalent")
            and entries[j].positive for j in range(len(entries))
        ):
            reasons.append("aggregate_source_already_represented")
            # Keep an unknown alternative explicit; no numeric duplicate credit exists here.
        entries[i] = replace(entry, positive=positive, reasons=tuple(sorted(set(reasons))), conflict_refs=_unique(conflicts))
    # Unscored dependence groups retain every record. Unknown models never merge contexts.
    neighbours = defaultdict(set)
    for i, (a, row, *_rest) in enumerate(details):
        tokens = {"family:" + a.source_family}
        for publication in row.get("publications") or ():
            if publication and (match := re.fullmatch(r"(?:PMID:)?([1-9][0-9]*)", publication, re.I)):
                tokens.add("PMID:" + match[1])
        for token in tokens:
            neighbours[token].add(i)
    pending = set(range(len(entries)))
    while pending:
        members, queue = set(), [min(pending)]
        while queue:
            i = queue.pop()
            if i in members:
                continue
            members.add(i)
            queue.extend(j for group in neighbours.values() if i in group for j in group - members)
        group_id = "dependency-" + _hash(sorted(token for token, group in neighbours.items() if group & members))
        for i in members:
            entries[i] = replace(entries[i], dependency_group=group_id)
        pending -= members
    return tuple(sorted(entries, key=lambda e: _bytes(e.source_ref)))


@dataclass(frozen=True, slots=True)
class _Association:
    subject: str
    rule: str
    adapter: str
    state: str
    namespace: str | None
    target: str | None
    condition: str | None
    reasons: tuple[str, ...]
    source_refs: tuple[SourcedRecord, ...] = ()
    evidence_refs: tuple[EvidenceReference, ...] = ()
    context: HypothesisContext | None = None


@dataclass(frozen=True, slots=True)
class AssociationTarget:
    """Reached positive association geometry, without its per-input evidence."""
    namespace: str | None
    target: str | None
    condition: str | None
    context: HypothesisContext | None = None

    def __post_init__(self):
        if (self.target is None and self.condition is None
                or any(value is not None and (type(value) is not str or not value)
                    for value in (self.namespace,self.target,self.condition))
                or self.context is not None and type(self.context) is not HypothesisContext):
            raise CandidateAdmissionError('invalid reached association target')


def association_targets(associations):
    targets={}
    for value in associations:
        if type(value) is not _Association:
            raise CandidateAdmissionError('unvalidated association discovery')
        condition=_normalise_curie(value.condition)
        if value.state!='match' or value.target is None and condition is None:
            continue
        target=AssociationTarget(value.namespace,
            value.target.split('.')[0] if value.namespace=='Ensembl' and value.target else value.target,
            condition,value.context)
        targets[_bytes(target)]=target
    return tuple(targets[key] for key in sorted(targets))


def _association_reaches_context(value,context,namespace,reach_cache):
    """Shared target/condition/model applicability for discovery and evaluation."""
    if value.namespace is not None and (value.namespace!=namespace
            or not _gene_join(namespace,value.target,context.locus_id or context.gene_id)):
        return False
    condition=_normalise_curie(value.condition)
    if condition is not None and condition not in reach_cache[context.native_disease_id]:
        return False
    if value.target is None and value.condition is None:
        return False
    if value.context is not None and (context.native_disease_id!=value.context.native_disease_id
            or _requirement_key(context.requirements[0])!=_requirement_key(value.context.requirements[0])):
        return False
    return True


def _verify_repeat_input(value: RepeatCatalogInput, reference_manifest_sha256: str):
    if (
        value.catalog.store.source_binding.get("adapter_version")
        not in (repeat.ADAPTER_VERSION, repeat.SYNTHETIC_ADAPTER_VERSION)
        or value.catalog.store.source_binding["contract_sha256"]
        != value.expected_contract_sha256
        or value.catalog.public_manifest_sha256 != value.expected_manifest_sha256
        or rebuild_disease_gene_evidence(value.catalog.store) != value.catalog.store
    ):
        raise ValueError
    if value.reference is None:
        return
    reference = value.reference
    contract = repeat.ReferenceIdentityContract.from_mapping(
        repeat._json_bytes(
            reference.contract_bytes, reference.expected_contract_sha256,
        )
    )
    manifest = repeat.PublicResourceManifest.from_mapping(
        repeat._json_bytes(
            reference.manifest_bytes, reference.expected_manifest_sha256,
        )
    )
    repeat.verify_reference_manifest_contract(
        contract, manifest,
        contract_sha256=reference.expected_contract_sha256,
    )
    output = next(
        item for item in manifest.outputs
        if item.path == contract.expected_output("contig_aliases").path
    )
    aliases = repeat._json_bytes(reference.aliases_bytes, output.sha256)
    if (
        len(reference.aliases_bytes) != output.size_bytes
        or aliases["schema"] != "genosapiens.contig-aliases"
        or aliases["contract_id"] != contract.contract_id
        or reference.expected_manifest_sha256 != reference_manifest_sha256
    ):
        raise ValueError


def _load_associations(inputs, sources):
    """Only fixed offline adapters. Exhaust every stream before accepting results."""
    inventory = inv.build_admission_inventory(inputs)
    associations, annotations, link_proofs, repeat_proofs = [], (), [], []
    if sources.annotation is not None:
        a = sources.annotation
        inv.build_admission_inventory(inputs, alleles=tuple(r.allele for r in a.requests))
        if any(r.identity_id != r.allele.native_identity_id for r in a.requests):
            raise ValueError
        annotations = annotation.replay_annotation_admission(
            a.artifacts, expected_binding_sha256=a.expected_binding_sha256, requests=a.requests)
        published = tuple(from_json(_bytes(item), AlleleEvidence) for item in json.loads(a.artifacts["alleles.json"]))
        inventory = inv.build_admission_inventory(inputs, alleles=published)
        required = {key for key in inventory.required_rule_keys if key[0] in {a.allele_id for a in inventory.alleles}
                    and key[1] in (inv.RULE_IDS[0], inv.RULE_IDS[1])}
        if {p.rule_key for p in annotations} != required:
            raise ValueError
        for p in annotations:
            associations.append(_Association(p.allele_id, p.rule, "annotation",
                {"admitted": "match", "evaluated_no_match": "no_match", "unresolved": "unresolved"}[p.outcome],
                p.target_namespace, p.target_id, None, (*p.reasons, "geometry:" + str(p.geometry_role),
                    "overlap:" + p.overlap_certainty), (*p.source_refs, p.provenance)))
    if any(source is not None for source in (sources.encode, sources.open_targets, sources.clinvar)) and sources.reference is None:
        raise ValueError
    if sources.reference is not None:
        if sources.reference.manifest_sha256 != inventory.reference_manifest_sha256:
            raise ValueError
        with ExitStack() as stack:
            links._reference(links._Files(stack), sources.reference)
    for name, source, rule in (("encode", sources.encode, inv.RULE_IDS[2]),
                               ("open_targets", sources.open_targets, inv.RULE_IDS[4]),
                               ("clinvar", sources.clinvar, inv.RULE_IDS[4])):
        if source is None:
            continue
        selected = tuple(a for a in inventory.alleles if (a.allele_id, rule) in inventory.required_rule_keys)
        if name == "encode":
            proofs = tuple(links.encode_regulatory_report(selected, resources=source, reference=sources.reference))
        elif name == "open_targets":
            proofs = tuple(links.open_targets_clinical_report(selected, source=source, reference=sources.reference))
        else:
            proofs = tuple(links.clinvar_clinical_report(selected, source=source, reference=sources.reference))
        if {p.allele_id for p in proofs} != {a.allele_id for a in selected}:
            raise ValueError
        natives = {a.allele_id: _native(a) for a in selected}
        for p in proofs:
            if p.native != natives[p.allele_id] or p.rule.replace("_", "-") != rule:
                raise ValueError
            associations.append(_Association(p.allele_id, rule, name,
                {"admitted": "match", "evaluated_no_match": "no_match", "unresolved": "unresolved"}[p.outcome],
                p.target_namespace, p.target_id, p.condition_id, p.reasons, (*p.source_refs, p.provenance)))
        link_proofs.extend(proofs)
    if sources.repeat_catalog is not None:
        r = sources.repeat_catalog
        _verify_repeat_input(r, inventory.reference_manifest_sha256)
        for allele in inventory.alleles:
            if (allele.allele_id, inv.RULE_IDS[3]) not in inventory.required_rule_keys:
                continue
            proofs = repeat.evaluate_repeat_catalog(allele, r.catalog, reference=r.reference)
            if not proofs:
                raise ValueError
            repeat_proofs.extend(proofs)
            for p in proofs:
                namespace = (r.catalog.store.assertions[0].source_family + ":locus" if p.context is not None else
                             "Ensembl" if p.gene_id else None)
                associations.append(_Association(allele.allele_id, inv.RULE_IDS[3], "repeat", p.association_status,
                    namespace, p.locus_id if p.context is not None else p.gene_id,
                    p.context.native_disease_id if p.context else None, p.reasons, (), p.evidence_refs, p.context))
    present = {"annotation": sources.annotation, "encode": sources.encode,
               "open_targets": sources.open_targets, "clinvar": sources.clinvar,
               "curated": None, "repeat": sources.repeat_catalog}
    emitted = {(p.subject, p.rule, p.adapter) for p in associations}
    derived = {a.allele_id for a in inventory.alleles}
    for subject, rule in inventory.required_rule_keys:
        for adapter in _RULE_ADAPTERS[rule]:
            if (subject, rule, adapter) in emitted:
                continue
            if subject in derived and present[adapter] is not None:
                raise ValueError
            reason = "source-unavailable" if present[adapter] is None else "native-source-unresolved"
            associations.append(_Association(subject, rule, adapter, "unresolved", None, None, None, (reason,)))
    return inventory, tuple(associations), annotations, tuple(link_proofs), tuple(repeat_proofs)


_PREPARED_ADMISSION_SEAL = object()


@dataclass(frozen=True, slots=True, repr=False)
class _PreparedAdmissionPhenotype:
    view: ScoringView
    profiles: DiseaseProfile | None
    evidence: tuple[DiseaseGeneEvidence, ...]
    ic_store: HpoICStore | None
    observations: Mapping
    limits: tuple[str, ...]
    reach_cache: Mapping
    profile_status: Mapping
    matches: Mapping
    pheno: tuple[PhenotypeAdmission, ...]
    seal: object


def _prepare_admission_phenotype(view, profiles, evidence, ic_store):
    """Case-local producer; all members are immutable and bound to source owners."""
    observations = defaultdict(list)
    limits = []
    for o in view.observations:
        if o.status != "present":
            continue
        resolution = ic_store.graph.resolve(o.source_hpo_id) if ic_store is not None else None
        if resolution is not None and resolution.canonical_id:
            observations[resolution.canonical_id].append(o)
        else:
            limits.append("unresolved-positive-observation")
    if ic_store is None:
        limits.append("ic-source-unavailable")
    if profiles is None:
        limits.append("profile-source-unavailable")
    if not observations:
        limits.append("no-resolved-positive-observation")
    reach_cache, edges = {}, _relation_edges(evidence)
    anchors = {assertion_native_disease_id(g.store, a) for g in evidence for a in g.store.assertions}
    annotations = profiles.store.annotations if profiles else ()
    anchors.update(a.source_fields["source_disease_id"] for a in annotations)
    for anchor in sorted(anchors - {None}):
        reach_cache[anchor] = _reach(anchor, edges)
    pheno, profile_status, matches = [], defaultdict(list), {}
    grouped = defaultdict(list)
    for a in annotations:
        native = a.source_fields["source_disease_id"]
        if native and a.source_fields["aspect"] == "P":
            profile_status[native].append(a)
        if not native or a.polarity != "positive" or a.source_fields["aspect"] != "P" or not a.resolution or not a.resolution.canonical_id:
            continue
        for term in sorted(observations):
            pair = (term, a.resolution.canonical_id)
            if pair not in matches:
                matches[pair] = ic_store.resnik(*pair)
            match = matches[pair]
            if match.information_content is not None and match.information_content > 0:
                grouped[native, *pair].append(a)
    for (native, term, target), members in sorted(grouped.items()):
        match = matches[term, target]
        ids = tuple(sorted(o.observation_id for o in observations[term]))
        pheno.append(PhenotypeAdmission("phenotype-" + _hash([native, term, target, match.mica_ids, match.information_content]),
            native, term, target, match.mica_ids, match.information_content, ids,
            tuple(sorted(c.contradiction_id for c in view.contradictions if set(c.observation_ids) & set(ids))),
            _unique(_ref(profiles, a.annotation_id) for a in members), "unresolved", (), ()))
    return _PreparedAdmissionPhenotype(
        view, profiles, evidence, ic_store,
        MappingProxyType({key: tuple(value) for key, value in observations.items()}),
        tuple(limits),
        MappingProxyType({key: MappingProxyType(value) for key, value in reach_cache.items()}),
        MappingProxyType({key: tuple(value) for key, value in profile_status.items()}),
        MappingProxyType(matches), tuple(pheno), _PREPARED_ADMISSION_SEAL,
    )


@dataclass(frozen=True, slots=True, repr=False)
class _PreparedAdmissionGraph:
    phenotype: _PreparedAdmissionPhenotype
    candidates: Mapping
    ledgers: Mapping
    supported: Mapping
    phenotype_state: Mapping
    pheno: tuple[PhenotypeAdmission, ...]
    candidate_by_target: Mapping
    candidate_by_condition: Mapping
    seed_contexts: Mapping
    seal: object


def _prepare_admission_graph(prepared, *, reached_targets, execution_config=None, progress_path=None, progress=None):
    """Close only contexts reached by positive assessed source associations."""
    if type(reached_targets) is not tuple or any(type(v) is not AssociationTarget for v in reached_targets):
        raise CandidateAdmissionError('invalid graph target closure')
    from genosapiens.wgs_runtime import analysis_progress
    processed = 0
    analysis_progress(progress, 'context_assembly', 'candidate_construction', counts={'candidates_processed':0})
    observations,limits,reach_cache = prepared.observations,prepared.limits,prepared.reach_cache
    profile_status,matches,pheno = prepared.profile_status,prepared.matches,prepared.pheno
    evidence,ic_store = prepared.evidence,prepared.ic_store
    # Each candidate alternative is anchored in actual source requirements. No
    # union-find over partial models: an unknown cannot bridge AD and AR.
    candidates, candidate_origins = {}, defaultdict(set)
    phenotype_anchors, phenotype_by_anchor = defaultdict(set), defaultdict(list)
    for p in pheno:
        for native in reach_cache[p.native_disease_id]:
            phenotype_anchors[native].add(p.native_disease_id)
            phenotype_by_anchor[native].append(p)
    support_index = defaultdict(lambda: defaultdict(list))
    public_by_id = {g.evidence_id: g for g in evidence}
    for public in evidence:
        for a in public.store.assertions:
            native = assertion_native_disease_id(public.store, a)
            namespace, target = _target(public, a)
            if not native or not target:
                continue
            support_index[namespace, target][public.evidence_id].append(a)
            possible = {native, *phenotype_anchors[native]}
            for anchor in sorted(possible):
                relation_refs = reach_cache[anchor][native][1]
                for requirement in a.requirements:
                    key = "context-" + _hash([anchor, namespace, target, _requirement_key(requirement)])
                    refs = _unique((*(_ref(public, k) for k in (a.record_id, *requirement.record_ids)), *relation_refs))
                    locus_bindings = (LocusBinding(_ref(public, a.locus_record_id), ("LocusId",)),) if a.locus_id else ()
                    context = HypothesisContext(key, anchor, anchor, a.gene_id, a.locus_id, requirement.inheritance,
                        "source-requirement-context-only", (requirement,), refs, locus_bindings)
                    if key in candidates:
                        old, _ = candidates[key]
                        context = replace(old, requirements=_unique((*old.requirements, requirement)),
                            evidence_refs=_unique((*old.evidence_refs, *refs)),
                            locus_bindings=_unique((*old.locus_bindings, *locus_bindings)))
                    candidates[key] = (context, namespace)
                    processed += 1
                    if processed % 256 == 0:
                        analysis_progress(progress, 'context_assembly', 'candidate_construction',
                            counts={'candidates_processed':processed}, unit='candidates', completed=processed)
                    if anchor == native:
                        candidate_origins[key].add("variant")
    targets_by_gene,targets_by_condition=defaultdict(list),defaultdict(list)
    for target in reached_targets:
        if target.target is not None:
            targets_by_gene[target.namespace,target.target].append(target)
        else:
            targets_by_condition[target.condition].append(target)
    reached={}
    for key,(context,namespace) in candidates.items():
        if 'variant' not in candidate_origins[key]:
            continue
        target=context.locus_id or context.gene_id
        stable=target.split('.')[0] if namespace=='Ensembl' else target
        possible=(*targets_by_gene.get((namespace,stable),()),
            *(value for native in reach_cache[context.native_disease_id]
                for value in targets_by_condition.get(native,())))
        if any(_association_reaches_context(value,context,namespace,reach_cache) for value in possible):
            reached[key]=(context,namespace)
    candidates=reached
    analysis_progress(progress, 'context_assembly', 'candidate_construction',
        counts={'candidates_processed':processed, 'candidates_retained':len(candidates)},
        unit='candidates', completed=processed)
    ledgers, supported, phenotype_state, context_pheno = {}, {}, {}, {}
    progress_callback = progress
    progress = None
    if execution_config is not None:
        from genosapiens.bounded_wgs_config import BoundedWgsExecutionConfig
        from genosapiens.wgs_admission_graph import GraphProgress, parallel_support
        if type(execution_config) is not BoundedWgsExecutionConfig or execution_config.public_queries.concurrent_batches != 10:
            raise CandidateAdmissionError('invalid graph execution configuration')
        progress = GraphProgress(progress_path, len(candidates), progress=progress_callback)
        with closing(parallel_support(candidates, support_index, public_by_id, reach_cache, progress=progress)) as results:
            for key, ledger in results:
                if key not in candidates or key in ledgers:
                    raise CandidateAdmissionError('invalid graph support coverage')
                ledgers[key] = ledger
        if len(ledgers) != len(candidates):
            raise CandidateAdmissionError('incomplete graph support')
        ledgers = {key: ledgers[key] for key in sorted(candidates)}
    else:
        for key, (context, namespace) in sorted(candidates.items()):
            reached = reach_cache[context.native_disease_id]
            indexed = tuple((public_by_id[public_id], assertions) for public_id, assertions in
                            support_index[namespace, context.locus_id or context.gene_id].items())
            ledgers[key] = _support(context, namespace, indexed, reached)

    # Requirements and gene/locus alternatives do not affect phenotype state.
    # Cache the exact former expression once per native disease in this graph.
    state_by_native = {}
    for key, (context, namespace) in sorted(candidates.items()):
        supported[key] = any(e.positive for e in ledgers[key])
        # A native profile match can reach this native gene context by an actual
        # directed/equivalence proof, never by sharing an optional mapped view.
        native_id = context.native_disease_id
        if native_id not in state_by_native:
            reached = reach_cache[native_id]
            proofs = tuple(phenotype_by_anchor[native_id])
            rows = [a for native in reached for a in profile_status.get(native, ())]
            informative_rows = [a for a in rows if a.polarity == "positive" and a.resolution and a.resolution.canonical_id]
            unknown_rows = any(a.polarity in ("unknown", "unresolved") or not a.resolution or not a.resolution.canonical_id for a in rows)
            unknown_comparison = any(matches[term, a.resolution.canonical_id].information_content is None
                                     for term in observations for a in informative_rows)
            state = "present" if proofs else "unknown" if (
                not informative_rows or not observations or ic_store is None or unknown_rows or unknown_comparison
                or "unresolved-positive-observation" in limits) else "absent"
            state_by_native[native_id] = (proofs, state)
        context_pheno[key], phenotype_state[key] = state_by_native[native_id]
    phenotype_contexts = defaultdict(list)
    for key in candidates:
        if supported[key]:
            for p in context_pheno[key]:
                phenotype_contexts[p.proof_id].append(key)
    pheno = tuple(replace(p, context_ids=tuple(sorted(phenotype_contexts[p.proof_id])),
        evidence_refs=_unique((*p.evidence_refs, *(ref for key in phenotype_contexts[p.proof_id]
            for ref in reach_cache[p.native_disease_id][candidates[key][0].native_disease_id][1]))),
        outcome="admitted" if phenotype_contexts[p.proof_id] else "unresolved",
        reasons=("informative-positive-phenotype-association", *(() if phenotype_contexts[p.proof_id]
                  else ("no-applicable-positive-disease-support",)))) for p in pheno)
    candidate_by_target, candidate_by_condition = defaultdict(list), defaultdict(list)
    for key, (context, namespace) in candidates.items():
        if "variant" not in candidate_origins[key]:
            continue
        target = context.locus_id or context.gene_id
        stable = target.split(".")[0] if namespace == "Ensembl" else target
        candidate_by_target[namespace, stable].append(key)
        for native in reach_cache[context.native_disease_id]:
            candidate_by_condition[native].append(key)
    seed_contexts = {}
    by_context = defaultdict(list)
    for proof in pheno:
        for key in proof.context_ids:
            by_context[key].append(proof)
    for key,(context,namespace) in candidates.items():
        proofs = by_context.get(key,())
        refs = _unique((*context.evidence_refs, *(e.source_ref for e in ledgers[key]),
            *(r for e in ledgers[key] for r in e.relation_refs),
            *(r for p in proofs for r in p.evidence_refs)))
        seed_contexts[key] = (replace(context,evidence_refs=refs),namespace,
            phenotype_state[key],tuple(sorted(p.proof_id for p in proofs)),ledgers[key])
    graph = _PreparedAdmissionGraph(prepared,_freeze_policy(candidates),_freeze_policy(ledgers),
        _freeze_policy(supported),_freeze_policy(phenotype_state),pheno,
        _freeze_policy(dict(candidate_by_target)),_freeze_policy(dict(candidate_by_condition)),
        _freeze_policy(seed_contexts),_PREPARED_ADMISSION_SEAL)
    if progress is not None:
        progress.write('complete', force=True)
    return graph


def _evaluate(inventory, associations, view, profiles, evidence, ic_store, binding,
              prepared=None, graph=None):
    """Deterministic fixed union/evaluator kernel. No files, scores or callbacks."""
    if prepared is None:
        prepared = _prepare_admission_phenotype(view, profiles, evidence, ic_store)
    if (type(prepared) is not _PreparedAdmissionPhenotype
            or prepared.seal is not _PREPARED_ADMISSION_SEAL
            or prepared.view is not view or prepared.profiles is not profiles
            or prepared.evidence is not evidence or prepared.ic_store is not ic_store):
        raise ValueError
    observations, limits, reach_cache = prepared.observations, prepared.limits, prepared.reach_cache
    profile_status, matches, pheno = prepared.profile_status, prepared.matches, prepared.pheno

    if graph is None:
        graph = _prepare_admission_graph(prepared,reached_targets=association_targets(associations))
    if (type(graph) is not _PreparedAdmissionGraph or graph.seal is not _PREPARED_ADMISSION_SEAL
            or graph.phenotype is not prepared):
        raise ValueError
    return _evaluate_graph(inventory, associations, prepared, graph, binding)


def _context_subject_hash(*, allele=None, source=None, reference_role=None, reference_identity_id=None):
    """The existing scientific subject identity, shared with discovery."""
    if (allele is None)==(source is None):raise CandidateAdmissionError('invalid context subject')
    if allele is not None:
        return _hash([_native(allele).identity,sorted({call.sample.sample_name for call in allele.calls})])
    return _hash([reference_role,reference_identity_id,source.record.contig,
        source.record.position,source.record.reference,source.alternate])


def _evaluate_graph(inventory, associations, prepared, graph, binding, *,
                    include_phenotype=True, subject_hashes=None, binding_sha256=None):
    """Shared evaluator over an authenticated graph or its exact job projection.

    The WGS coordinator supplies all applicable candidates and already computed
    subject hashes. Only the once-per-case phenotype seed emission is separable
    from the unchanged association evaluation below.
    """
    reach_cache, limits = prepared.reach_cache, prepared.limits
    candidates,ledgers,supported = graph.candidates,graph.ledgers,graph.supported
    phenotype_state,pheno = graph.phenotype_state,graph.pheno
    candidate_by_target,candidate_by_condition = graph.candidate_by_target,graph.candidate_by_condition
    routes, allele_refs, decision_ids, decisions = defaultdict(set), defaultdict(list), defaultdict(set), []
    for p in pheno if include_phenotype else ():
        for key in p.context_ids:
            routes[key].add("phenotype_first")
    by_id = {a.allele_id: a for a in inventory.alleles}
    source_by_id = {a.source.source_allele_id: a for a in inventory.source_alleles}
    binding_sha256 = _hash(binding) if binding_sha256 is None else binding_sha256
    subject_hashes = {} if subject_hashes is None else dict(subject_hashes)
    for p in associations:
        allele = by_id.get(p.subject)
        if p.subject not in subject_hashes:
            if allele is not None:
                subject_hashes[p.subject] = _context_subject_hash(allele=allele)
            else:
                source = source_by_id[p.subject].source
                subject_hashes[p.subject] = _context_subject_hash(source=source,
                    reference_role=inventory.reference_role,reference_identity_id=inventory.reference_identity_id)
        subject_science = subject_hashes[p.subject]
        eligible = []
        stable = p.target.split(".")[0] if p.namespace == "Ensembl" and p.target else p.target
        condition = _normalise_curie(p.condition)
        possible = candidate_by_target.get((p.namespace,stable),()) if p.target else candidate_by_condition.get(condition,())
        for key in possible:
            context, namespace = candidates[key]
            if not _association_reaches_context(p,context,namespace,reach_cache):
                continue
            eligible.append(key)
        # A source-level no-match cannot be made disease-context-specific by
        # enumerating every disease assertion for its reported target.  Only a
        # positive target match is eligible to carry a context and its support
        # ledger; no-match/unresolved routes remain fully accounted once at the
        # source association level.
        evaluated_keys = sorted(eligible) if p.state == "match" else []
        for key in evaluated_keys or [None]:
            context = candidates[key][0] if key else None
            support = "present" if key and supported[key] else "absent" if key and any(
                e.polarity in ("negative", "no_positive_support") and e.applicability != "unknown"
                and e.requirement_state == "exact_alternative" for e in ledgers[key]
            ) else "unknown"
            outcome = "evaluated_no_match" if p.state == "no_match" else "admitted" if p.state == "match" and support == "present" else "unresolved"
            resulting = (key,) if outcome == "admitted" else ()
            reasons = tuple(sorted(set((*p.reasons, "context-only",
                *(("source-condition-spelling-normalized",) if condition != p.condition else ()), *( () if support == "present" else
                ("no-applicable-positive-disease-support",))))))
            phenotype = phenotype_state[key] if key else "unknown"
            scientific = [subject_science, p.rule, p.adapter, p.state, support, phenotype, outcome,
                          p.namespace, p.target, p.condition, key, reasons]
            decision_id = "decision-" + _hash(scientific)
            condition_refs = reach_cache[context.native_disease_id][condition][1] if context and condition else ()
            refs = _unique((*p.evidence_refs, *(context.evidence_refs if context else ()),
                *condition_refs,
                *(e.source_ref for e in ledgers.get(key, ())),
                *(r for e in ledgers.get(key, ()) for r in e.relation_refs)))
            decisions.append(AdmissionDecision(decision_id, p.subject, subject_science, p.rule, p.adapter, p.state,
                support, phenotype, outcome, p.namespace, p.target, p.condition, key, resulting, reasons, refs,
                p.source_refs, binding_sha256, _hash([POLICY_SHA256, p.rule, p.adapter]), ledgers.get(key, ())))
            if outcome == "admitted":
                routes[key].add("variant_first")
                allele_refs[key].append(AlleleReference(allele.allele_id, allele.native_identity_id))
                decision_ids[key].add(decision_id)
    seeds = []
    for key in sorted(routes):
        context,namespace,state,proof_ids,ledger = graph.seed_contexts[key]
        seeds.append(CandidateContextSeed(context,namespace,tuple(sorted(routes[key])),
            _unique(allele_refs[key]),state,proof_ids,tuple(sorted(decision_ids[key])),ledger))
    # Multiple source proofs can project to the same closed association (for
    # example, several transcript records for one allele/gene).  Retain every
    # proof upstream, but emit one content-identical decision so downstream
    # source-accounting has a set, not a duplicated multiset.
    return tuple(seeds), _unique(decisions), pheno, tuple(sorted(set(limits)))


def _prepare_admission_evidence(*, view, profiles, genes, ic_store, repeat_catalog, authenticated_resources=None):
    view = _typed(view, ScoringView)
    ids = {observation.observation_id for observation in view.observations}
    if any(not set(item.observation_ids) <= ids for item in view.contradictions):
        raise ValueError
    if ic_store is not None:
        _verify_hpo_store(ic_store, authenticated_resources)
    if profiles is not None:
        if type(profiles) is not DiseaseProfile:
            raise TypeError
        if ic_store is not None:
            _verify_profile_store(profiles, ic_store)
    evidence = list(genes or ())
    if repeat_catalog is not None:
        catalog = repeat_catalog.catalog
        if any(
            item.evidence_id == catalog.evidence_id and item.store != catalog.store
            for item in evidence
        ):
            raise ValueError
        if not any(item.evidence_id == catalog.evidence_id for item in evidence):
            evidence.append(catalog)
    if (
        len({item.evidence_id for item in evidence}) != len(evidence)
        or profiles is not None
        and profiles.evidence_id in {item.evidence_id for item in evidence}
    ):
        raise ValueError
    if any(type(item) is not DiseaseGeneEvidence for item in evidence):
        raise TypeError
    evidence = tuple(sorted((
        item if not item.aggregations else replace(item, aggregations=())
        for item in evidence
    ), key=lambda item: item.evidence_id))
    for item in evidence:
        _verify_gene_evidence(item)
    return view, evidence


def _attach_public_association_refs(seeds, decisions, associations):
    registry = defaultdict(dict)
    for association in associations:
        for record in association.source_refs:
            if record.public_manifest_sha256 is None:
                continue
            key = "admission-public-" + hashlib.sha256(_record_bytes(record)).hexdigest()
            registry[key][record.source_record_id] = record
    public_sources = tuple(
        PublicSourcedEvidence(
            key, next(iter(records.values())).public_manifest_sha256, records,
        )
        for key, records in sorted(registry.items())
    )
    refs_by_record = {
        _record_bytes(record): EvidenceReference(
            public.evidence_id, record.source_record_id,
            public.public_manifest_sha256,
        )
        for public in public_sources for record in public.records.values()
    }
    decisions = tuple(replace(
        decision,
        evidence_refs=_unique((
            *decision.evidence_refs,
            *(refs_by_record[content] for record in decision.source_refs
              if (content := _record_bytes(record)) in refs_by_record),
        )),
    ) for decision in decisions)
    context_refs = defaultdict(list)
    for decision in decisions:
        for key in decision.context_ids:
            context_refs[key].extend(decision.evidence_refs)
    seeds = tuple(replace(
        seed,
        context=replace(seed.context, evidence_refs=_unique((
            *seed.context.evidence_refs,
            *context_refs[seed.context.context_id],
        ))),
    ) if context_refs.get(seed.context.context_id) else seed for seed in seeds)
    return seeds, decisions, public_sources


@dataclass(frozen=True, slots=True)
class _LedgerSourceAllele:
    source: object


@dataclass(frozen=True, slots=True)
class _LedgerInventoryView:
    alleles: tuple[AlleleEvidence, ...]
    source_alleles: tuple[_LedgerSourceAllele, ...]
    reference_role: ReferenceRole
    reference_identity_id: str


def _load_sharded_associations(
    inputs: tuple[LedgerAdmissionInput, ...],
    proofs: tuple[annotation.AnnotationAssociationProof, ...],
    sources: ShardedAdmissionSources,
    *,
    reference_role: ReferenceRole,
    reference_identity_id: str,
    reference_manifest_sha256: str,
    clinical_associations: ShardedClinicalAssociations | None = None,
    regulatory_associations: ShardedRegulatoryAssociations | None = None,
    resource_preparation=None,
):
    if (
        type(inputs) is not tuple
        or not inputs
        or any(type(value) is not LedgerAdmissionInput for value in inputs)
        or type(proofs) is not tuple
        or any(type(value) is not annotation.AnnotationAssociationProof for value in proofs)
        or type(sources) is not ShardedAdmissionSources
        or type(reference_role) is not ReferenceRole
        or type(reference_identity_id) is not str
        or not reference_identity_id
        or type(reference_manifest_sha256) is not str
        or _SHA.fullmatch(reference_manifest_sha256) is None
        or clinical_associations is not None
        and type(clinical_associations) is not ShardedClinicalAssociations
        or regulatory_associations is not None
        and type(regulatory_associations) is not ShardedRegulatoryAssociations
    ):
        raise ValueError
    if any(
        right.retained_ordinal <= left.retained_ordinal
        for left, right in zip(inputs, inputs[1:])
    ):
        raise ValueError
    for name in (
        "parent_source_manifest_sha256",
        "parent_working_set_manifest_sha256",
        "native_source_contract_sha256",
    ):
        if len({getattr(value, name) for value in inputs}) != 1:
            raise ValueError
    alleles = tuple(value.allele for value in inputs if value.allele is not None)
    by_id = {value.allele_id: value for value in alleles}
    if len(by_id) != len(alleles):
        raise ValueError
    natives = {allele.allele_id: _native(allele) for allele in alleles}
    for allele in alleles:
        selected = natives[allele.allele_id]
        identity = selected.identity
        sequence = getattr(identity, "sequence", identity)
        identity_role = getattr(identity, "reference_role", None)
        identity_reference = getattr(identity, "reference_identity_id", None)
        if (
            selected.reference_manifest_sha256 != reference_manifest_sha256
            or identity_role is not reference_role
            and getattr(sequence, "reference_role", None) is not reference_role
            or identity_reference != reference_identity_id
            and getattr(sequence, "reference_identity_id", None) != reference_identity_id
        ):
            raise ValueError
    unresolved = tuple(
        _LedgerSourceAllele(value.unresolved_source)
        for value in inputs if value.unresolved_source is not None
    )
    inventory = _LedgerInventoryView(
        alleles, unresolved, reference_role, reference_identity_id,
    )
    required = {
        (value.subject_id, rule)
        for value in inputs for rule in value.required_rule_ids
    }
    associations = []
    required_annotation = {
        key for key in required
        if key[0] in by_id and key[1] in (inv.RULE_IDS[0], inv.RULE_IDS[1])
    }
    if sources.annotation_binding_sha256 is None:
        if proofs:
            raise ValueError
    elif {proof.rule_key for proof in proofs} != required_annotation:
        raise ValueError
    for proof in proofs:
        allele = by_id.get(proof.allele_id)
        if allele is None or proof.native != natives[allele.allele_id]:
            raise ValueError
        associations.append(_Association(
            proof.allele_id,
            proof.rule,
            "annotation",
            {
                "admitted": "match",
                "evaluated_no_match": "no_match",
                "unresolved": "unresolved",
            }[proof.outcome],
            proof.target_namespace,
            proof.target_id,
            None,
            (*proof.reasons, "geometry:" + str(proof.geometry_role),
             "overlap:" + proof.overlap_certainty),
            (*proof.source_refs, proof.provenance),
        ))

    if (
        any(source is not None for source in (
            sources.encode, sources.open_targets, sources.clinvar,
        ))
        and sources.reference is None
    ):
        raise ValueError
    if sources.reference is not None and sources.reference.manifest_sha256 != reference_manifest_sha256:
        raise ValueError
    if (sources.regulatory_index_manifest_sha256 is None) != (regulatory_associations is None):
        raise ValueError
    if regulatory_associations is not None:
        if (sources.encode is None or sources.reference is None
                or regulatory_associations.index_manifest_sha256 != sources.regulatory_index_manifest_sha256
                or regulatory_associations.source_manifest_sha256s != tuple(pin.manifest_sha256 for pin in sources.encode)
                or regulatory_associations.reference_manifest_sha256 != reference_manifest_sha256
                or any(proof.allele_id not in by_id
                       or (proof.allele_id, inv.RULE_IDS[2]) not in required
                       for proof in regulatory_associations.proofs)):
            raise ValueError
    if sources.reference is not None and (
        sources.encode is not None and regulatory_associations is None
        or clinical_associations is None and (
            regulatory_associations is None
            or sources.open_targets is not None or sources.clinvar is not None)
    ):
        with ExitStack() as stack:
            links._reference(links._Files(stack), sources.reference)
    if clinical_associations is not None:
        if (
            sources.open_targets is None
            or sources.clinvar is None
            or sources.reference is None
            or clinical_associations.open_targets_manifest_sha256
            != sources.open_targets.manifest_sha256
            or clinical_associations.clinvar_manifest_sha256
            != sources.clinvar.manifest_sha256
            or clinical_associations.reference_manifest_sha256
            != sources.reference.manifest_sha256
            or clinical_associations.reference_manifest_sha256
            != reference_manifest_sha256
        ):
            raise ValueError
    link_proofs = []
    for name, source, rule in (
        ("encode", sources.encode, inv.RULE_IDS[2]),
        ("open_targets", sources.open_targets, inv.RULE_IDS[4]),
        ("clinvar", sources.clinvar, inv.RULE_IDS[4]),
    ):
        if source is None:
            continue
        selected = tuple(
            allele for allele in alleles if (allele.allele_id, rule) in required
        )
        selected_ids = {allele.allele_id for allele in selected}
        if name == "encode":
            emitted = (tuple(proof for proof in regulatory_associations.proofs
                             if proof.allele_id in selected_ids)
                       if regulatory_associations is not None else
                       tuple(links.encode_regulatory_report(
                           selected, resources=source, reference=sources.reference)))
        elif clinical_associations is not None:
            emitted = tuple(
                value.proof for value in clinical_associations.proofs
                if value.adapter == name
                and value.proof.allele_id in selected_ids
            )
        elif name == "open_targets":
            emitted = tuple(links.open_targets_clinical_report(
                selected, source=source, reference=sources.reference,
            ))
        else:
            emitted = tuple(links.clinvar_clinical_report(
                selected, source=source, reference=sources.reference,
            ))
        if {proof.allele_id for proof in emitted} != selected_ids:
            raise ValueError
        for proof in emitted:
            if (
                proof.native != natives[proof.allele_id]
                or proof.rule.replace("_", "-") != rule
            ):
                raise ValueError
            associations.append(_Association(
                proof.allele_id, rule, name,
                {
                    "admitted": "match",
                    "evaluated_no_match": "no_match",
                    "unresolved": "unresolved",
                }[proof.outcome],
                proof.target_namespace, proof.target_id, proof.condition_id,
                proof.reasons, (*proof.source_refs, proof.provenance),
            ))
        link_proofs.extend(emitted)

    repeat_proofs = []
    if sources.repeat_catalog is not None:
        selected_source = sources.repeat_catalog
        if resource_preparation is None:
            _verify_repeat_input(selected_source, reference_manifest_sha256)
        elif (type(resource_preparation) is not PreparedShardedResources
                or resource_preparation.seal is not _PREPARED_SHARDED_SEAL or resource_preparation.sources is not sources
                or resource_preparation.repeat_reference_manifest_sha256 is not None
                and resource_preparation.repeat_reference_manifest_sha256 != reference_manifest_sha256):
            raise ValueError
        for allele in alleles:
            if (allele.allele_id, inv.RULE_IDS[3]) not in required:
                continue
            emitted = repeat.evaluate_repeat_catalog(
                allele, selected_source.catalog,
                reference=selected_source.reference,
            )
            if not emitted:
                raise ValueError
            repeat_proofs.extend(emitted)
            for proof in emitted:
                namespace = (
                    selected_source.catalog.store.assertions[0].source_family + ":locus"
                    if proof.context is not None else
                    "Ensembl" if proof.gene_id else None
                )
                associations.append(_Association(
                    allele.allele_id, inv.RULE_IDS[3], "repeat",
                    proof.association_status, namespace,
                    proof.locus_id if proof.context is not None else proof.gene_id,
                    proof.context.native_disease_id if proof.context else None,
                    proof.reasons, (), proof.evidence_refs, proof.context,
                ))

    present = {
        "annotation": sources.annotation_binding_sha256,
        "encode": sources.encode,
        "open_targets": sources.open_targets,
        "clinvar": sources.clinvar,
        "curated": None,
        "repeat": sources.repeat_catalog,
    }
    emitted_keys = {
        (value.subject, value.rule, value.adapter) for value in associations
    }
    for subject, rule in required:
        for adapter in _RULE_ADAPTERS[rule]:
            if (subject, rule, adapter) in emitted_keys:
                continue
            if subject in by_id and present[adapter] is not None:
                raise ValueError
            associations.append(_Association(
                subject, rule, adapter, "unresolved", None, None, None,
                ("source-unavailable" if present[adapter] is None
                 else "native-source-unresolved",),
            ))
    return (
        inventory, tuple(associations), tuple(link_proofs), tuple(repeat_proofs),
    )


_PREPARED_SHARDED_SEAL = object()


def _freeze_policy(value):
    if isinstance(value, dict):
        return MappingProxyType({key: _freeze_policy(child) for key, child in value.items()})
    if isinstance(value, list):
        return tuple(_freeze_policy(child) for child in value)
    return value


@dataclass(frozen=True, slots=True, repr=False)
class PreparedShardedResources:
    """Case-local invariant source validation, independent of graph reachability."""
    sources: ShardedAdmissionSources
    repeat_reference_manifest_sha256: str | None
    seal: object


def prepare_sharded_resources(sources):
    if type(sources) is not ShardedAdmissionSources:
        raise CandidateAdmissionError('invalid sharded association sources')
    repeat_reference=None
    if sources.repeat_catalog is not None:
        if sources.repeat_catalog.reference is not None:
            repeat_reference=sources.repeat_catalog.reference.expected_manifest_sha256
        _verify_repeat_input(sources.repeat_catalog,repeat_reference)
    return PreparedShardedResources(sources,repeat_reference,_PREPARED_SHARDED_SEAL)


@dataclass(frozen=True, slots=True, repr=False)
class _PreparedShardedAdmission:
    sources: ShardedAdmissionSources
    source_view: ScoringView
    view: ScoringView
    profiles: DiseaseProfile | None
    genes: tuple[DiseaseGeneEvidence, ...] | None
    ic_store: HpoICStore | None
    clinical_index_manifest_sha256: str | None
    policy: bytes
    policy_object: Mapping
    policy_sha256: str
    view_sha256: str
    evidence: tuple[DiseaseGeneEvidence, ...]
    phenotype: _PreparedAdmissionPhenotype
    graph: _PreparedAdmissionGraph
    resources: PreparedShardedResources
    reached_targets: frozenset[bytes]
    seal: object


def prepare_sharded_candidate_admission(
    *, sources: ShardedAdmissionSources, view: ScoringView,
    profiles: DiseaseProfile | None, genes: tuple[DiseaseGeneEvidence, ...] | None,
    ic_store: HpoICStore | None, reached_targets: tuple[AssociationTarget, ...],
    clinical_index_manifest_sha256: str | None = None,
    execution_config=None, graph_progress_path=None, resource_preparation=None, progress=None, authenticated_resources=None,
) -> _PreparedShardedAdmission:
    """Prepare fixed policy and phenotype relationships once for this case/view.

    The caller owns this bounded-lifetime object. It holds strong references to
    its immutable resource owners and is never shared through a global cache.
    """
    from genosapiens.wgs_runtime import analysis_progress
    from genosapiens.wgs_failures import WgsFailure
    analysis_progress(progress, 'context_assembly', 'resource_prepare')
    try:
        if resource_preparation is None:
            resource_preparation=prepare_sharded_resources(sources)
        if (type(resource_preparation) is not PreparedShardedResources
                or resource_preparation.seal is not _PREPARED_SHARDED_SEAL
                or resource_preparation.sources is not sources):
            raise ValueError
        source_view = view
        policy = sharded_admission_policy_bytes(
            sources=sources, profiles=profiles, genes=genes, ic_store=ic_store,
            clinical_index_manifest_sha256=clinical_index_manifest_sha256,
        )
        view, evidence = _prepare_admission_evidence(
            view=view, profiles=profiles, genes=genes, ic_store=ic_store,
            repeat_catalog=sources.repeat_catalog, authenticated_resources=authenticated_resources,
        )
        phenotype = _prepare_admission_phenotype(view,profiles,evidence,ic_store)
        return _PreparedShardedAdmission(
            sources, source_view, view, profiles, genes, ic_store,
            clinical_index_manifest_sha256, policy,
            _freeze_policy(json.loads(policy)), hashlib.sha256(policy).hexdigest(),
            _hash(view), evidence,
            phenotype,_prepare_admission_graph(phenotype, reached_targets=reached_targets,execution_config=execution_config,
                progress_path=graph_progress_path,progress=progress),resource_preparation,
            frozenset(_bytes(value) for value in reached_targets),
            _PREPARED_SHARDED_SEAL,
        )
    except WgsFailure:
        raise
    except (ValueError, TypeError, RuntimeError, OSError, KeyError, AttributeError):
        raise CandidateAdmissionError("invalid sharded admission preparation") from None


def _sharded_batch_binding(inputs,input_digests,prepared,*,reference_role,
        reference_identity_id,reference_manifest_sha256):
    return {
        "schema": "genosapiens.sharded-candidate-admission",
        "policy_sha256": prepared.policy_sha256,
        "policy": prepared.policy_object,
        "executable_version": EXECUTABLE_VERSION,
        "parent_source_manifest_sha256":
            inputs[0].parent_source_manifest_sha256,
        "parent_working_set_manifest_sha256":
            inputs[0].parent_working_set_manifest_sha256,
        "native_source_contract_sha256":
            inputs[0].native_source_contract_sha256,
        "first_retained_ordinal": inputs[0].retained_ordinal,
        "last_retained_ordinal": inputs[-1].retained_ordinal,
        # Preserve the existing canonical hash of ordered per-row hashes.
        # The private store separately seals concatenated original bytes.
        "input_rows_sha256": _hash(input_digests),
        "observation_sha256": prepared.view_sha256,
        "reference_role": reference_role.value,
        "reference_identity_id": reference_identity_id,
        "reference_manifest_sha256": reference_manifest_sha256,
        "genes_availability": "unavailable" if prepared.genes is None else "supplied",
        "scope": "bounded-context-admission-input; global-grouping-pending",
        }


def _generate_projected_candidate_admission(*,discovered,prepared,context_executor,expected_policy_sha256):
    """Internal bounded role: consume discovery's authenticated subjects."""
    from genosapiens.wgs_association_projection import AssociationBatch
    from genosapiens.wgs_context_admission import ContextAdmissionExecutor
    if (type(discovered) is not AssociationBatch or type(prepared) is not _PreparedShardedAdmission
            or prepared.seal is not _PREPARED_SHARDED_SEAL
            or type(context_executor) is not ContextAdmissionExecutor
            or prepared is not context_executor.prepared or prepared.policy_sha256!=expected_policy_sha256
            or prepared.clinical_index_manifest_sha256!=discovered.clinical_index_manifest_sha256):
        raise CandidateAdmissionError('invalid projected candidate admission owner')
    associations,_,_=discovered.require(prepared.resources)
    applicable=tuple(proof for row in discovered.rows for proof,eligible in zip(row.associations,
        (True,)*len(row.associations) if row.profile_applicability is None else row.profile_applicability,
        strict=True) if eligible)
    if not {_bytes(value) for value in association_targets(applicable)}<=prepared.reached_targets:
        raise CandidateAdmissionError('projected candidate targets are outside graph closure')
    manifest=discovered.manifest
    binding=_sharded_batch_binding(discovered.subjects,[value.admission_ref[4] for value in discovered.subjects],
        prepared,reference_role=manifest.reference_role,reference_identity_id=manifest.reference_identity_id,
        reference_manifest_sha256=manifest.reference_manifest_sha256)
    return context_executor.batch(subjects=discovered.subjects,associations=associations,prepared=prepared,
        binding=binding,association_owner=discovered)


def generate_sharded_candidate_admission(
    *,
    inputs: tuple[LedgerAdmissionInput, ...],
    annotation_proofs: tuple[annotation.AnnotationAssociationProof, ...],
    sources: ShardedAdmissionSources,
    reference_role: ReferenceRole,
    reference_identity_id: str,
    reference_manifest_sha256: str,
    view: ScoringView,
    profiles: DiseaseProfile | None,
    genes: tuple[DiseaseGeneEvidence, ...] | None,
    ic_store: HpoICStore | None,
    expected_policy_sha256: str,
    clinical_associations: ShardedClinicalAssociations | None = None,
    regulatory_associations: ShardedRegulatoryAssociations | None = None,
    prepared: _PreparedShardedAdmission | None = None,
) -> ShardedCandidateAdmission:
    """Evaluate a complete strict full-input batch for ordinary callers."""

    try:
        clinical_sha256 = (
            None if clinical_associations is None
            else clinical_associations.index_manifest_sha256
        )
        resources=prepared.resources if prepared is not None else prepare_sharded_resources(sources)
        inventory,associations,link_proofs,repeat_proofs=_load_sharded_associations(
            inputs,annotation_proofs,sources,reference_role=reference_role,
            reference_identity_id=reference_identity_id,reference_manifest_sha256=reference_manifest_sha256,
            clinical_associations=clinical_associations,regulatory_associations=regulatory_associations,
            resource_preparation=resources)
        if prepared is None:
            prepared = prepare_sharded_candidate_admission(
                sources=sources, view=view, profiles=profiles, genes=genes,
                ic_store=ic_store,
                clinical_index_manifest_sha256=clinical_sha256,
                reached_targets=association_targets(associations),resource_preparation=resources,
            )
        if (type(prepared) is not _PreparedShardedAdmission
                or prepared.seal is not _PREPARED_SHARDED_SEAL
                or prepared.sources is not sources or prepared.source_view is not view
                or prepared.profiles is not profiles or prepared.genes is not genes
                or prepared.ic_store is not ic_store
                or prepared.clinical_index_manifest_sha256 != clinical_sha256):
            raise ValueError
        if (
            type(expected_policy_sha256) is not str
            or prepared.policy_sha256 != expected_policy_sha256
        ):
            raise ValueError
        view, evidence = prepared.view, prepared.evidence
        if not {_bytes(value) for value in association_targets(associations)} <= prepared.reached_targets:
            raise ValueError
        binding=_sharded_batch_binding(inputs,[hashlib.sha256(value.as_bytes()).hexdigest() for value in inputs],
            prepared,reference_role=reference_role,reference_identity_id=reference_identity_id,
            reference_manifest_sha256=reference_manifest_sha256)
        contexts, decisions, phenotype, limits = _evaluate(
            inventory, associations, view, profiles, evidence, ic_store, binding,
            prepared.phenotype,prepared.graph,
        )
        contexts, decisions, public_sources = _attach_public_association_refs(
            contexts, decisions, associations,
        )
        if genes is None:
            limits = tuple(sorted((*limits, "disease-gene-source-unavailable")))
        expected = {
            (value.subject_id, rule, adapter)
            for value in inputs for rule in value.required_rule_ids
            for adapter in _RULE_ADAPTERS[rule]
        }
        if {
            (value.subject_id, value.rule_id, value.adapter)
            for value in decisions
        } != expected:
            raise ValueError
        return ShardedCandidateAdmission(
            inputs, contexts, decisions, phenotype, limits, public_sources,
            annotation_proofs, link_proofs, repeat_proofs,
            MappingProxyType(binding),
        )
    except (
        ValueError, TypeError, RuntimeError, OSError, KeyError, IndexError,
        AttributeError, StopIteration, OverflowError,
    ):
        raise CandidateAdmissionError(
            "invalid or incomplete sharded candidate admission inputs"
        ) from None


def generate_candidate_admission(*, inputs: inv.SuppliedInputs, sources: AdmissionSources,
                                 view: ScoringView, profiles: DiseaseProfile | None,
                                 genes: tuple[DiseaseGeneEvidence, ...] | None, ic_store: HpoICStore | None,
                                 expected_policy_sha256: str) -> CandidateAdmission:
    """Replay all original sources, reconcile coverage, then union context seeds."""
    try:
        policy = admission_policy_bytes(sources=sources, profiles=profiles, genes=genes, ic_store=ic_store)
        if type(expected_policy_sha256) is not str or hashlib.sha256(policy).hexdigest() != expected_policy_sha256:
            raise ValueError
        view, evidence = _prepare_admission_evidence(
            view=view,
            profiles=profiles,
            genes=genes,
            ic_store=ic_store,
            repeat_catalog=sources.repeat_catalog,
        )
        inventory, associations, annotations, link_proofs, repeat_proofs = _load_associations(inputs, sources)
        binding = {"schema": SCHEMA, "policy_sha256": expected_policy_sha256, "policy": json.loads(policy),
            "executable_version": EXECUTABLE_VERSION, "input_sha256": inputs.input_sha256,
            "inventory_sha256": hashlib.sha256(inv.admission_inventory_as_json(inventory).encode("ascii")).hexdigest(),
            "observation_sha256": _hash(view), "reference_role": inventory.reference_role.value,
            "reference_identity_id": inventory.reference_identity_id,
            "reference_manifest_sha256": inventory.reference_manifest_sha256,
            "reference_concordance_sha256": inventory.reference_concordance_sha256,
            "genes_availability": "unavailable" if genes is None else "supplied",
            "scope": "context-admission-only; no-genotype-or-partner-search-or-controlled-approval"}
        seeds, decisions, pheno, limits = _evaluate(inventory, associations, view, profiles, evidence, ic_store, binding)
        if genes is None:
            limits = tuple(sorted((*limits, "disease-gene-source-unavailable")))
        seeds, decisions, public_sources = _attach_public_association_refs(
            seeds, decisions, associations,
        )
        result = CandidateAdmission(inventory, view, seeds, decisions, pheno, limits, evidence, profiles,
                                    public_sources, annotations, link_proofs, repeat_proofs, binding)
        _account(result, inputs)
        return result
    except (ValueError, TypeError, RuntimeError, OSError, KeyError, IndexError, AttributeError, StopIteration, OverflowError):
        raise CandidateAdmissionError("invalid or incomplete candidate admission inputs") from None


def _account(result, inputs):
    inventory = result.inventory
    inv.validate_admission_accounting(inventory, inputs=inputs, alleles=inventory.alleles,
        source_record_ids=tuple(r.qc.source_record_id for r in inventory.records),
        allele_ids=tuple(a.allele_id for a in inventory.alleles),
        unresolved_entry_ids=tuple(u.entry_id for u in inventory.unresolved_entries),
        evaluated_rule_keys=tuple(sorted({(d.subject_id, d.rule_id) for d in result.decisions})))
    expected = {(subject, rule, adapter) for subject, rule in inventory.required_rule_keys for adapter in _RULE_ADAPTERS[rule]}
    if {(d.subject_id, d.rule_id, d.adapter) for d in result.decisions} != expected:
        raise CandidateAdmissionError("incomplete source adapter coverage")


def validate_candidate_admission(result: CandidateAdmission, *, inputs: inv.SuppliedInputs,
                                 sources: AdmissionSources, view: ScoringView, profiles: DiseaseProfile | None,
                                 genes: tuple[DiseaseGeneEvidence, ...] | None, ic_store: HpoICStore | None,
                                 expected_policy_sha256: str) -> None:
    """Recompute from independent inputs; a valid subset or self-hash cannot pass."""
    try:
        if type(result) is not CandidateAdmission:
            raise ValueError
        _account(result, inputs)
        expected = generate_candidate_admission(inputs=inputs, sources=sources, view=view, profiles=profiles,
            genes=genes, ic_store=ic_store, expected_policy_sha256=expected_policy_sha256)
        if result.as_bytes() != expected.as_bytes():
            raise ValueError
    except (ValueError, TypeError, RuntimeError, OSError, KeyError, AttributeError, OverflowError):
        raise CandidateAdmissionError("candidate admission differs from complete replayed inputs") from None
