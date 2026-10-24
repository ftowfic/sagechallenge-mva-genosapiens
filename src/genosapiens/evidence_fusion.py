"""Source-replayed exact evidence fusion for tasks 9.4 and 9.5.

This module supplies no production policy. The public ranker replays the existing
hypothesis/source adapters and derives a closed set of score facts; callers cannot
submit numerical charges or repeat lineage. Numerical policy selection, release
acceptance, and controlled-use freezing remain task 9.6 work.
"""

from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass, replace
from fractions import Fraction
import hashlib
from types import MappingProxyType, SimpleNamespace
from typing import Literal, get_args

from genosapiens.disease_gene_evidence import (
    AggregationPolicy, DiseaseGeneContext, aggregate_disease_gene_evidence,
)
from genosapiens.evidence_model import (
    AlleleEvidence, ComponentName, Completeness, DiseaseGeneEvidence,
    DiseaseProfile, EvidenceModelError, EvidenceReference, EvidenceDocument,
    GenotypeHypothesis, HypothesisPhase, IdentityStage, Observation,
    PublicSourcedEvidence, ScoringView, ScoreAdjustment, ScoreComponent,
    ScoreLedger, SourcedRecord, TieKey, _Controlled, _exact, _ids, _sha,
    _typed, resolve_evidence_reference,
)
from genosapiens.evidence_serialization import _canonical, _loads, _plain
from genosapiens.genotype_hypotheses import (
    ConstructedHypothesis, HypothesisGeneration, HypothesisReplay,
    replay_genotype_hypotheses,
)
from genosapiens.hpo_ic import HpoICStore, hpo_ic_as_bytes
from genosapiens.wgs_provenance import ContextProvenance, public_records
from genosapiens import longtr_admission as longtr
from genosapiens.phenotype_scoring import PhenotypeScore
from genosapiens.track1_sensitivity import (
    PhenotypeTransform, RankTie,
    _AuthenticatedSensitivityResources, _authenticated_ic_hash, order_ledgers,
    score_primary_phenotypes,
)


SCHEMA = "genosapiens.evidence-fusion"
RESULT_SCHEMA = "genosapiens.evidence-fusion-result"
EXECUTABLE_VERSION = "exact-affine-member-aware"
TRANSFORM = "exact-affine-clamp"
_PREPARED_GENERATION_SEAL = object()

Aggregation = Literal["shared", "member_mean", "shared_plus_member_mean"]
ChargeKind = Literal[
    "disease_gene", "molecular", "inheritance", "population", "clinical",
    "technical", "conflict",
]
SourceFactName = Literal[
    "molecular.admission.match",
    "molecular.admission.no_match",
    "molecular.admission.unresolved",
    "molecular.consequence_requirement.supported",
    "molecular.consequence_requirement.contradicted",
    "molecular.consequence_requirement.unresolved",
    "molecular.mechanism_direction.supported",
    "molecular.mechanism_direction.contradicted",
    "molecular.mechanism_direction.unresolved",
    "molecular.allele_origin.supported",
    "molecular.allele_origin.contradicted",
    "molecular.allele_origin.unresolved",
    "inheritance.model.supported",
    "inheritance.model.contradicted",
    "inheritance.model.unresolved",
    "clinical.assertion.positive",
    "clinical.assertion.negative",
    "clinical.assertion.unresolved",
    "clinical.assertion.conflict",
    "technical.read_validation_required",
    "technical.enumeration.truncated",
    "technical.enumeration.failed",
]
TieField = Literal[
    "phenotype", "disease_gene", "molecular", "inheritance", "population",
    "clinical", "technical_and_conflict_penalties", "completeness", "phase",
    "native_disease_id", "disease_id", "gene_id", "locus_id", "model",
    "hypothesis_id",
]

_COMPONENTS = get_args(ComponentName)
_NON_PHENOTYPE = _COMPONENTS[1:]
_OWNER = {
    "disease_gene": "disease_gene",
    "molecular": "molecular",
    "inheritance": "inheritance",
    "population": "population",
    "clinical": "clinical",
    "technical": "technical_and_conflict_penalties",
    "conflict": "technical_and_conflict_penalties",
}
_EXACT_TIES = {*_COMPONENTS, "completeness", "phase"}
_SUPPORT_TIES = set(_COMPONENTS) - {"technical_and_conflict_penalties"}
_COMPLETENESS = {
    Completeness.COMPLETE: "2", Completeness.UNKNOWN: "1", Completeness.INCOMPLETE: "0",
}
_PHASE = {
    HypothesisPhase.IN_TRANS: "4", HypothesisPhase.NOT_APPLICABLE: "3",
    HypothesisPhase.UNKNOWN: "2", HypothesisPhase.CONFLICTING: "1",
    HypothesisPhase.IN_CIS: "0",
}
_FACTS = get_args(SourceFactName)
_FACT_OWNER = {
    **dict.fromkeys(_FACTS[:12], "molecular"),
    **dict.fromkeys(_FACTS[12:15], "inheritance"),
    **dict.fromkeys(_FACTS[15:18], "clinical"),
    _FACTS[18]: "conflict",
    **dict.fromkeys(_FACTS[19:], "technical"),
}
_SOURCE_AGGREGATIONS = {
    "disease_gene": "shared",
    "molecular": "member_mean",
    "inheritance": "shared",
    "population": "shared",
    "clinical": "member_mean",
    "technical_and_conflict_penalties": "shared_plus_member_mean",
}
_ERRORS = (
    ValueError, TypeError, RuntimeError, AttributeError, KeyError, IndexError,
    OverflowError, RecursionError, StopIteration,
)
_PREPARED_FUSION_SEAL = object()


class EvidenceFusionError(EvidenceModelError):
    """Static diagnostics only; evidence and serialized results remain sensitive."""


def _bytes(value: object) -> bytes:
    return (_canonical(value) + "\n").encode("utf-8")


def _hash(value: object) -> str:
    return hashlib.sha256(_bytes(value)).hexdigest()


def _unique(values):
    by_value = {_canonical(value): value for value in values}
    return tuple(by_value[key] for key in sorted(by_value))


@dataclass(frozen=True, slots=True, repr=False)
class ExactAffineClamp(_Controlled):
    """Closed exact transform: clamp(raw * slope + intercept, lower, upper)."""

    slope: str
    intercept: str
    lower: str
    upper: str

    def __post_init__(self):
        _Controlled.__post_init__(self)
        slope, lower, upper = _exact(self.slope), _exact(self.lower), _exact(self.upper)
        _exact(self.intercept)
        if slope <= 0 or lower >= upper:
            raise EvidenceFusionError("invalid exact affine transform")

    def apply(self, raw: Fraction) -> Fraction:
        return min(_exact(self.upper), max(
            _exact(self.lower), raw * _exact(self.slope) + _exact(self.intercept),
        ))


@dataclass(frozen=True, slots=True, repr=False)
class ComponentPolicy(_Controlled):
    """One non-phenotype transform and candidate-independent aggregation rule."""

    name: ComponentName
    transform: ExactAffineClamp
    raw_neutral: str
    neutral_value: str
    weight: str
    aggregation: Aggregation

    def __post_init__(self):
        _Controlled.__post_init__(self)
        raw_neutral = _exact(self.raw_neutral)
        neutral, weight = _exact(self.neutral_value), _exact(self.weight)
        if self.name == "phenotype" or weight < 0 or self.transform.apply(raw_neutral) != neutral:
            raise EvidenceFusionError("invalid component policy")
        if self.name == "technical_and_conflict_penalties" and (
            _exact(self.transform.slope) != 1
            or _exact(self.transform.intercept) != 0
            or _exact(self.transform.lower) != 0
            or raw_neutral != 0
            or neutral != 0
        ):
            raise EvidenceFusionError("technical penalties require an exact identity clamp from zero")

    @property
    def parameters(self) -> Mapping[str, object]:
        pair_scale = {
            "shared": "hypothesis-level-shared-once",
            "member_mean": "arithmetic-mean-of-normalized-members-with-raw-neutral",
            "shared_plus_member_mean": "shared-once-plus-arithmetic-member-mean",
        }[self.aggregation]
        return MappingProxyType({
            "operation": TRANSFORM,
            "component": self.name,
            "aggregation": self.aggregation,
            "slope": self.transform.slope,
            "intercept": self.transform.intercept,
            "lower": self.transform.lower,
            "upper": self.transform.upper,
            "raw_neutral": self.raw_neutral,
            "neutral_value": self.neutral_value,
            "missing": "null-with-declared-neutral-no-weight-redistribution",
            "pair_scale": pair_scale,
        })

    @property
    def transform_sha256(self) -> str:
        return _hash(self.parameters)


@dataclass(frozen=True, slots=True, repr=False)
class PenaltyPolicy(_Controlled):
    """Explicit raw penalty amounts; all enter only the technical component."""

    weak_partner_floor: str
    weak_partner_scale: str
    weak_partner_cap: str
    incompleteness: str
    phase_unknown: str
    phase_conflicting: str
    phase_in_cis: str

    def __post_init__(self):
        _Controlled.__post_init__(self)
        floor = _exact(self.weak_partner_floor)
        scale, cap, incomplete, unknown, conflicting, cis = map(_exact, (
            self.weak_partner_scale, self.weak_partner_cap, self.incompleteness,
            self.phase_unknown, self.phase_conflicting, self.phase_in_cis,
        ))
        if (scale <= 0 or cap <= 0 or incomplete <= 0 or unknown <= 0
                or conflicting < unknown or cis < conflicting or floor < 0):
            raise EvidenceFusionError("invalid fusion penalty policy")


@dataclass(frozen=True, slots=True, repr=False)
class SourceFactMapping(_Controlled):
    """Policy value for one closed replay-derived fact, never a caller fact."""

    fact: SourceFactName
    value: str | None

    def __post_init__(self):
        _Controlled.__post_init__(self)
        if self.value is not None:
            value = _exact(self.value)
            if _FACT_OWNER[self.fact] in ("technical", "conflict") and value < 0:
                raise EvidenceFusionError("penalty fact values must be nonnegative")


@dataclass(frozen=True, slots=True, repr=False)
class TieRule(_Controlled):
    """Closed tie field. Numeric/text behavior is fixed by the field vocabulary."""

    field: TieField
    direction: Literal["ascending", "descending"]

    def __post_init__(self):
        _Controlled.__post_init__(self)
        if ((self.field in _SUPPORT_TIES or self.field in ("completeness", "phase"))
                and self.direction != "descending"):
            raise EvidenceFusionError("support and state tie fields must be descending")
        if (self.field == "technical_and_conflict_penalties" or self.field not in _EXACT_TIES
                and self.field not in _SUPPORT_TIES) and self.direction != "ascending":
            raise EvidenceFusionError("penalty and text tie fields must be ascending")

    @property
    def kind(self) -> Literal["exact", "text"]:
        return "exact" if self.field in _EXACT_TIES else "text"


@dataclass(frozen=True, slots=True, repr=False)
class FusionPolicy(_Controlled):
    """Caller-pinned policy input, not a task 9.6 freeze or production default."""

    policy_id: str
    phenotype: PhenotypeTransform
    components: tuple[ComponentPolicy, ...]
    source_facts: tuple[SourceFactMapping, ...]
    disease_gene_aggregation: AggregationPolicy
    penalties: PenaltyPolicy
    ties: tuple[TieRule, ...]
    derivation_sha256: str
    derivation_scope: Literal["public", "independently_synthetic"]

    def __post_init__(self):
        _Controlled.__post_init__(self)
        _ids((self.policy_id,))
        _sha(self.derivation_sha256)
        if tuple(item.name for item in self.components) != _NON_PHENOTYPE:
            raise EvidenceFusionError("fusion policy must declare every component in canonical order")
        if (tuple(item.aggregation for item in self.components)
                != tuple(_SOURCE_AGGREGATIONS[name] for name in _NON_PHENOTYPE)):
            raise EvidenceFusionError("component aggregation differs from the source adapter")
        if tuple(item.fact for item in self.source_facts) != _FACTS:
            raise EvidenceFusionError("fusion policy must map every source fact in canonical order")
        fields = tuple(item.field for item in self.ties)
        if (len(set(fields)) != len(fields) or not {"completeness", "phase"} <= set(fields)
                or not fields or fields[-1] != "hypothesis_id"):
            raise EvidenceFusionError("invalid deterministic tie contract")
        molecular = next(item for item in self.components if item.name == "molecular")
        floor = _exact(self.penalties.weak_partner_floor)
        if (molecular.aggregation == "shared"
                or not _exact(molecular.transform.lower) < floor <= _exact(molecular.transform.upper)):
            raise EvidenceFusionError("weak-partner policy requires member molecular evidence")


@dataclass(frozen=True, slots=True, repr=False)
class EvidenceCharge(_Controlled):
    """One atomic source claim with exactly one score-component owner.

    Parent/projection views reuse lineage_id + evidence_key. Presentation as
    shared or member evidence cannot create another charge. Distinct projected
    member evidence uses a distinct evidence key and therefore remains distinct.
    """

    kind: ChargeKind
    lineage_id: str
    evidence_key: str
    value: str | None
    member_id: str | None = None
    evidence_refs: tuple[EvidenceReference, ...] = ()
    source_values: tuple[SourcedRecord, ...] = ()

    def __post_init__(self):
        _Controlled.__post_init__(self)
        _ids((self.lineage_id, self.evidence_key))
        if self.member_id is not None:
            _ids((self.member_id,))
        if self.value is not None and self.kind in ("technical", "conflict") and _exact(self.value) < 0:
            raise EvidenceFusionError("penalty evidence must be nonnegative")
        if self.value is not None:
            _exact(self.value)

    @property
    def component(self) -> ComponentName:
        return _OWNER[self.kind]

    @property
    def ownership_key(self) -> tuple[str, str]:
        return self.lineage_id, self.evidence_key


@dataclass(frozen=True, slots=True, repr=False)
class HypothesisEvidence(_Controlled):
    """Closed ranker input; qualification labels are deliberately not inputs."""

    hypothesis: GenotypeHypothesis
    phenotype_score: PhenotypeScore
    charges: tuple[EvidenceCharge, ...]
    phenotype_evidence_refs: tuple[EvidenceReference, ...] = ()
    phenotype_source_values: tuple[SourcedRecord, ...] | None = None

    def __post_init__(self):
        _Controlled.__post_init__(self)
        members = set(self.hypothesis.evidence_subject_ids)
        if any(item.member_id is not None and item.member_id not in members for item in self.charges):
            raise EvidenceFusionError("evidence charge member is outside the hypothesis")


@dataclass(frozen=True, slots=True, repr=False)
class FusionRanking(_Controlled):
    """Deterministic controlled result; replay still requires independent inputs."""

    policy_id: str
    policy_sha256: str
    input_sha256: str
    ledgers: tuple[ScoreLedger, ...]
    order: tuple[str, ...]
    tie_rules: tuple[RankTie, ...]
    provenance: Mapping[str, object]

    def __post_init__(self):
        _Controlled.__post_init__(self)
        _ids((self.policy_id,))
        _sha(self.policy_sha256)
        _sha(self.input_sha256)
        identifiers = tuple(item.hypothesis_id for item in self.ledgers)
        _ids(identifiers)
        _ids(self.order)
        _ids(tuple(item.name for item in self.tie_rules))
        tie_contract = tuple((item.name, item.direction) for item in self.tie_rules)
        if (identifiers != tuple(sorted(identifiers))
                or len(self.order) != len(identifiers) or set(self.order) != set(identifiers)
                or any(item.configuration_id != self.policy_id
                       or item.configuration_sha256 != self.policy_sha256
                       or tuple((key.name, key.direction) for key in item.tie_keys) != tie_contract
                       for item in self.ledgers)
                or order_ledgers(self.ledgers, self.tie_rules) != self.order):
            raise EvidenceFusionError("invalid fusion ranking membership")

    def as_bytes(self) -> bytes:
        return _bytes({
            "schema": RESULT_SCHEMA,
            "classification": "controlled",
            "llm_policy": "no-llm",
            "controlled_use_ready": False,
            "competitive_output_allowed": False,
            "result": self,
        })


def policy_as_json(policy: FusionPolicy) -> bytes:
    """Canonical policy bytes. There is intentionally no module-level default."""
    try:
        policy = _typed(policy, FusionPolicy)
        return _bytes({"schema": SCHEMA, "executable_version": EXECUTABLE_VERSION, "policy": policy})
    except _ERRORS:
        raise EvidenceFusionError("cannot encode evidence-fusion policy") from None


def policy_from_json(content: bytes, *, expected_sha256: str) -> FusionPolicy:
    """Require the independently supplied hash and canonical policy bytes."""
    try:
        _sha(expected_sha256)
        if type(content) is not bytes or hashlib.sha256(content).hexdigest() != expected_sha256:
            raise ValueError
        data = _loads(content)
        if (type(data) is not dict or set(data) != {"schema", "executable_version", "policy"}
                or data["schema"] != SCHEMA or data["executable_version"] != EXECUTABLE_VERSION):
            raise ValueError
        policy = _typed(data["policy"], FusionPolicy, decode=True)
        if policy_as_json(policy) != content:
            raise ValueError
        return policy
    except _ERRORS:
        raise EvidenceFusionError("invalid evidence-fusion policy bytes or binding") from None


def ranking_from_json(content: bytes) -> FusionRanking:
    """Structural decode only; use rank_hypotheses(..., expected_result=...) to replay."""
    try:
        if type(content) is not bytes:
            raise ValueError
        data = _loads(content)
        if (type(data) is not dict
                or set(data) != {"schema", "classification", "llm_policy", "controlled_use_ready",
                                 "competitive_output_allowed", "result"}
                or (data["schema"], data["classification"], data["llm_policy"],
                    data["controlled_use_ready"], data["competitive_output_allowed"])
                != (RESULT_SCHEMA, "controlled", "no-llm", False, False)):
            raise ValueError
        result = _typed(data["result"], FusionRanking, decode=True)
        if result.as_bytes() != content:
            raise ValueError
        return result
    except _ERRORS:
        raise EvidenceFusionError("invalid evidence-fusion result JSON") from None


def ranking_ties(policy: FusionPolicy) -> tuple[RankTie, ...]:
    policy = _typed(policy, FusionPolicy)
    return tuple(RankTie(
        item.field,
        item.field if item.field in _COMPONENTS else "fixed",
        item.kind,
        item.direction,
    ) for item in policy.ties)


def _replay_generation(replay: HypothesisReplay) -> HypothesisGeneration:
    try:
        return replay_genotype_hypotheses(replay)
    except Exception:
        raise EvidenceFusionError("typed hypothesis replay required") from None


@dataclass(frozen=True, slots=True, repr=False, init=False)
class _PreparedFusion:
    """Run-local authenticated products; construction stays inside this module."""

    replay: HypothesisReplay
    generation: HypothesisGeneration
    document: EvidenceDocument
    ic_store: HpoICStore
    profile_evidence_id: str | None
    policy: FusionPolicy
    phenotype_scores: Mapping[str, PhenotypeScore]
    ranking: FusionRanking
    _seal: object
    _binding: tuple[object, ...]

    def __repr__(self):
        return "<controlled prepared fusion; no-llm>"


@dataclass(frozen=True, slots=True, repr=False, init=False)
class _PreparedGeneration:
    """Run-local source replay shared with evidence assembly and fusion."""

    replay: HypothesisReplay
    generation: HypothesisGeneration
    _seal: object
    _binding: tuple[object, ...]

    def __repr__(self):
        return "<controlled prepared generation; no-llm>"


def _prepare_generation(replay: HypothesisReplay) -> _PreparedGeneration:
    generation = _replay_generation(replay)
    prepared = object.__new__(_PreparedGeneration)
    binding = (replay, generation)
    for name, value in (
        ("replay", replay), ("generation", generation),
        ("_seal", _PREPARED_GENERATION_SEAL), ("_binding", binding),
    ):
        object.__setattr__(prepared, name, value)
    return prepared


def _is_prepared_generation(value: object, replay: HypothesisReplay) -> bool:
    if (
        type(value) is not _PreparedGeneration
        or getattr(value, "_seal", None) is not _PREPARED_GENERATION_SEAL
        or value.replay is not replay
    ):
        return False
    expected = (value.replay, value.generation)
    binding = getattr(value, "_binding", ())
    return len(binding) == len(expected) and all(
        actual is original for actual, original in zip(binding, expected, strict=True)
    )


def _prepared_fusion(
    document: EvidenceDocument, *, replay: HypothesisReplay, ic_store: HpoICStore,
    profile_evidence_id: str | None, policy_bytes: bytes, policy_sha256: str,
    derivation_bytes: bytes, expected_derivation_sha256: str,
    expected_result: FusionRanking | None = None,
    _generation_context: _PreparedGeneration | None = None,
    _authenticated_resources: _AuthenticatedSensitivityResources | None = None,
) -> _PreparedFusion:
    """Authenticate generation and ranking once for one in-process execution."""
    if _generation_context is None:
        generation = _replay_generation(replay)
    elif _is_prepared_generation(_generation_context, replay):
        generation = _generation_context.generation
    else:
        raise EvidenceFusionError("authenticated prepared generation required")
    ic_store = _bind_ic_store(replay, ic_store)
    ic_store_sha256 = None
    if _authenticated_resources is not None:
        _authenticated_ic_hash(_authenticated_resources, ic_store)
        ic_store_sha256 = _authenticated_resources.store_sha256
    policy = policy_from_json(policy_bytes, expected_sha256=policy_sha256)
    document = _typed(
        document,
        EvidenceDocument,
        _reuse_existing=_generation_context is not None,
    )
    evidence, scores, audit = _adapt_source_evidence(
        generation, document, policy=policy, ic_store=ic_store,
        profile_evidence_id=profile_evidence_id,
        _ic_store_sha256=ic_store_sha256,
    )
    ranking = _rank_adapted_hypotheses(
        evidence, policy_bytes=policy_bytes, policy_sha256=policy_sha256,
        derivation_bytes=derivation_bytes,
        expected_derivation_sha256=expected_derivation_sha256,
        source_replay=audit, expected_result=expected_result,
    )
    prepared = object.__new__(_PreparedFusion)
    bound = (
        replay, generation, document, ic_store, policy, scores, ranking,
    )
    for name, value in (
        ("replay", replay), ("generation", generation), ("document", document),
        ("ic_store", ic_store), ("profile_evidence_id", profile_evidence_id),
        ("policy", policy), ("phenotype_scores", scores), ("ranking", ranking),
        ("_seal", _PREPARED_FUSION_SEAL),
        ("_binding", bound),
    ):
        object.__setattr__(prepared, name, value)
    return prepared


def _is_prepared_fusion(value: object) -> bool:
    if (type(value) is not _PreparedFusion
            or getattr(value, "_seal", None) is not _PREPARED_FUSION_SEAL):
        return False
    expected = (
        value.replay, value.generation, value.document, value.ic_store,
        value.policy, value.phenotype_scores, value.ranking,
    )
    binding = getattr(value, "_binding", ())
    return len(binding) == len(expected) and all(
        actual is original for actual, original in zip(binding, expected, strict=True)
    )


def _bind_ic_store(replay: HypothesisReplay, supplied: HpoICStore) -> HpoICStore:
    replayed = replay.ic_store
    if (type(replayed) is not HpoICStore or type(supplied) is not HpoICStore
            or replayed is not supplied
            and hpo_ic_as_bytes(replayed) != hpo_ic_as_bytes(supplied)):
        raise EvidenceFusionError("phenotype IC store differs from hypothesis replay")
    return replayed


def _bind_document(generation: HypothesisGeneration, document: EvidenceDocument) -> None:
    generation = _typed(generation, HypothesisGeneration)
    document = _typed(document, EvidenceDocument)
    expected_hypotheses = tuple(item.hypothesis for item in generation.hypotheses)
    expected_profiles = () if generation.admission.profiles is None else (generation.admission.profiles,)
    required_public = {
        item.evidence_id: item for item in (
            *generation.admission.public_sources,
            *((generation.clinical.public_sourced_evidence) if generation.clinical is not None else ()),
        )
    }
    actual_public = {item.evidence_id: item for item in document.public_sourced_evidence}
    if (document.score_ledgers
            or document.scoring_view() != generation.admission.view
            or document.disease_profiles != expected_profiles
            or document.disease_gene_evidence != generation.admission.disease_gene_evidence
            or document.alleles != generation.alleles
            or document.hypotheses != expected_hypotheses
            or any(actual_public.get(key) != value for key, value in required_public.items())):
        raise EvidenceFusionError("evidence document differs from replayed generation")


def _longtr_parent_root(allele):
    data = longtr.metadata(allele)
    if not data:
        return None
    native = next(item for item in allele.identities
                  if item.identity_id == allele.native_identity_id)
    representation = data.get("repeat_representation")
    if representation == "allele_projection":
        if longtr.is_projection(allele) is not True:
            raise EvidenceFusionError("invalid LongTR projection lineage")
    elif representation == "locus_call":
        source = native.source
        source_mappings = tuple(item for item in allele.identities
                                if item.stage is IdentityStage.SOURCE)
        calls = allele.calls or ()
        profile = data.get("source_format") or {}
        expected_profile = {
            "repeat_identity": longtr.FORMAT,
            "input_sha256": source.source_sha256,
            "reference_role": native.identity.reference_role.value,
            "reference_identity_id": native.identity.reference_identity_id,
            "reference_manifest_sha256": native.reference_manifest_sha256,
        }
        expected_parent_id = "longtr-locus-" + hashlib.sha256(_canonical((
            source.source_sha256, source.source_record_id, data.get("sample_name"),
        )).encode("ascii")).hexdigest()
        expected_identity_sha256 = hashlib.sha256(
            _canonical(native.identity).encode("ascii")
        ).hexdigest()
        origin = (source.source_dataset, source.source_release,
                  source.source_record_id, source.source_sha256)
        if (data.get("protocol_sha256") != longtr.POLICY_SHA256
                or data.get("parent_call_id") != allele.allele_id
                or data.get("parent_call_id") != expected_parent_id
                or data.get("parent_identity_id") != allele.native_identity_id
                or data.get("parent_identity_sha256") != expected_identity_sha256
                or native.transformation_id != longtr.FORMAT
                or any(profile.get(key) != value for key, value in expected_profile.items())
                or not source_mappings or len(calls) != len(source_mappings)
                or set(native.derived_from_identity_ids)
                != {item.identity_id for item in source_mappings}
                or any(
                    item.identity != native.identity
                    or item.reference_manifest_sha256 != native.reference_manifest_sha256
                    or item.identity_id != (item.source.source_fields or {}).get("source_allele_id")
                    or (item.source.source_dataset, item.source.source_release,
                        item.source.source_record_id, item.source.source_sha256) != origin
                    or ((item.source.source_fields or {}).get("native_interpretation") or {}) != data
                    for item in source_mappings
                )
                or sorted(_canonical(item.source) for item in calls)
                != sorted(_canonical(item.source) for item in source_mappings)
                or any(item.sample is None or item.sample.sample_name != data.get("sample_name")
                       for item in calls)):
            raise EvidenceFusionError("invalid LongTR parent lineage")
    else:
        raise EvidenceFusionError("invalid LongTR representation lineage")
    source = native.source
    return (
        "longtr-parent", source.source_dataset, source.source_release,
        source.source_record_id, source.source_sha256, data["sample_name"],
        data["parent_call_id"], data["parent_identity_sha256"],
        native.reference_manifest_sha256, data["protocol_sha256"],
    )


def _lineage_for_allele(allele) -> str:
    longtr_root = _longtr_parent_root(allele)
    if longtr_root is not None:
        return "identity-lineage-" + _hash((longtr_root,))
    mappings = {item.identity_id: item for item in allele.identities}
    pending, seen, roots = [allele.native_identity_id], set(), []
    while pending:
        key = pending.pop()
        if key in seen:
            continue
        seen.add(key)
        mapping = mappings[key]
        if mapping.stage is IdentityStage.SOURCE:
            roots.append((mapping.identity_id, mapping.source.source_dataset,
                          mapping.source.source_release, mapping.source.source_record_id,
                          mapping.source.source_sha256))
        pending.extend(mapping.derived_from_identity_ids)
    if not roots:
        raise EvidenceFusionError("source identity lineage is unavailable")
    return "identity-lineage-" + _hash(tuple(sorted(roots)))


def _lineage_for_sources(values) -> str:
    sources = _unique(values)
    if not sources:
        raise EvidenceFusionError("source lineage is unavailable")
    return "source-lineage-" + _hash(tuple(
        (item.source_dataset, item.source_release, item.source_record_id, item.source_sha256)
        for item in sources
    ))


def _lineage_for_refs(values) -> str:
    refs = _unique(values)
    if not refs:
        raise EvidenceFusionError("public evidence lineage is unavailable")
    return "public-lineage-" + _hash(tuple(
        (item.evidence_id, item.record_id, item.public_manifest_sha256) for item in refs
    ))


def _mapped_charge(
    policy: FusionPolicy, fact: SourceFactName, *, lineage_id: str, evidence_key: str,
    member_id: str | None = None, evidence_refs=(), source_values=(),
) -> tuple[EvidenceCharge, ...]:
    mapping = next(item for item in policy.source_facts if item.fact == fact)
    if member_id is not None:
        evidence_key += "-member-" + _hash(member_id)[:16]
    return (EvidenceCharge(
        _FACT_OWNER[fact], lineage_id, evidence_key, mapping.value, member_id,
        _unique(evidence_refs), _unique(source_values),
    ),)


def _disease_context(hypothesis: GenotypeHypothesis) -> DiseaseGeneContext | None:
    return disease_gene_context(hypothesis.context)


def disease_gene_context(context) -> DiseaseGeneContext | None:
    names = (
        "allelic_requirement", "required_allele_count", "consequence_requirement",
        "mechanism_direction", "allele_origin",
    )
    alternatives = {tuple(getattr(item, name) for name in names)
                    for item in context.requirements}
    if len(alternatives) > 1:
        return None
    values = next(iter(alternatives), (None,) * len(names))
    return DiseaseGeneContext(
        context.disease_id, context.gene_id,
        context.inheritance, context.native_disease_id,
        *values, locus_id=context.locus_id,
    )


def _disease_gene_charges(generation, constructed, policy):
    return disease_context_charges(constructed.hypothesis.context, generation.admission.disease_gene_evidence, policy)


def disease_context_charges(hypothesis_context, gene_evidence, policy, *, prepared_gene_evidence=None):
    if prepared_gene_evidence is not None:
        from genosapiens.genotype_evaluation import PreparedGeneEvidence
        if type(prepared_gene_evidence) is not PreparedGeneEvidence:
            raise EvidenceFusionError("invalid prepared gene evidence")
        prepared_gene_evidence.require(gene_evidence)
    context = disease_gene_context(hypothesis_context)
    if context is None:
        return (), ()
    charges, audits = [], []
    referenced = {item.evidence_id for item in hypothesis_context.evidence_refs}
    for public in gene_evidence:
        if public.evidence_id not in referenced:
            continue
        aggregation = aggregate_disease_gene_evidence(
            public.store, context, policy=policy.disease_gene_aggregation,
            _prepared=None if prepared_gene_evidence is None else prepared_gene_evidence.stores[public.evidence_id],
        )
        audits.append(aggregation)
        record_ids = {item.record_id for item in aggregation.ledger}
        record_ids.update(key for item in aggregation.ledger
                          for key in (*item.relation_record_ids, *item.discounted_by_record_ids))
        record_ids.update(key for item in aggregation.penalties
                          for key in (*item.positive_record_ids, *item.negative_record_ids))
        refs = tuple(EvidenceReference(public.evidence_id, key, public.public_manifest_sha256)
                     for key in sorted(record_ids))
        lineage = _lineage_for_refs(refs) if refs else None
        if lineage is not None and any(item.raw_support is not None for item in aggregation.ledger):
            charges.append(EvidenceCharge(
                "disease_gene", lineage,
                "disease-gene-support-" + _hash((public.evidence_id, context)),
                str(Fraction(str(aggregation.positive_support))), None, refs,
            ))
        for penalty in aggregation.penalties:
            penalty_refs = tuple(EvidenceReference(
                public.evidence_id, key, public.public_manifest_sha256,
            ) for key in sorted((*penalty.positive_record_ids, *penalty.negative_record_ids)))
            charges.append(EvidenceCharge(
                "conflict", _lineage_for_refs(penalty_refs),
                "disease-gene-conflict-" + penalty.penalty_id,
                str(Fraction(str(penalty.amount))), None, penalty_refs,
            ))
    return tuple(charges), tuple(audits)


def _admission_charges(generation, constructed, policy, alleles):
    seeds = tuple(item for item in generation.admission.contexts
                  if item.context.context_id == constructed.admission_context_id)
    if len(seeds) != 1:
        raise EvidenceFusionError("replayed admission context is ambiguous")
    selected = set(seeds[0].decision_ids)
    members = set(constructed.hypothesis.evidence_subject_ids)
    return admission_decision_charges(generation.admission.decisions, selected, members, policy, alleles)


def admission_decision_charges(source_decisions, selected, members, policy, alleles):
    """Shared source adapter after the caller closes the selected decision set."""
    decisions = {}
    for item in source_decisions:
        if (item.decision_id not in selected or item.subject_id not in alleles
                or item.subject_id not in members):
            continue
        old = decisions.get(item.decision_id)
        if old is not None and (old.association_state != item.association_state
                                or old.subject_id != item.subject_id):
            raise EvidenceFusionError("replayed admission decision is inconsistent")
        decisions[item.decision_id] = item if old is None else replace(
            old,
            evidence_refs=_unique((*old.evidence_refs, *item.evidence_refs)),
            source_refs=_unique((*old.source_refs, *item.source_refs)),
        )
    charges = []
    for key in sorted(decisions):
        item = decisions[key]
        fact = "molecular.admission." + item.association_state
        charges.extend(_mapped_charge(
            policy, fact, lineage_id=_lineage_for_allele(alleles[item.subject_id]),
            evidence_key="admission-decision-" + item.decision_id,
            member_id=item.subject_id, evidence_refs=item.evidence_refs,
            source_values=item.source_refs,
        ))
    return tuple(charges)


def _requirement_charges(constructed, policy, alleles):
    return requirement_member_charges(constructed.members, policy, alleles)


def requirement_member_charges(members, policy, alleles):
    charges = []
    owned_fields = ("consequence_requirement", "mechanism_direction", "allele_origin")
    for member in members:
        for item in member.requirement_assessments:
            if item.field not in owned_fields or getattr(item.requirement, item.field) is None:
                continue
            fact = "molecular." + item.field + "." + item.state
            sources = tuple(source for basis in item.bases for source in basis.source_refs)
            refs = (*item.requirement_refs,
                    *(ref for basis in item.bases for ref in basis.evidence_refs))
            charges.extend(_mapped_charge(
                policy, fact, lineage_id=_lineage_for_allele(alleles[member.allele_id]),
                evidence_key="requirement-" + _hash((
                    item.requirement_index, item.field, item.requirement,
                    item.policy_sha256,
                )), member_id=member.allele_id, evidence_refs=refs,
                source_values=sources,
            ))
    return tuple(charges)


def _inheritance_charges(constructed, policy):
    model = constructed.model
    fact = "inheritance.model." + model.state
    refs = tuple(ref for item in model.requirement_facts for ref in item.source_refs)
    sources = (model.binding, *(source for member in model.members for source in member.locus.sources))
    return _mapped_charge(
        policy, fact, lineage_id=_lineage_for_sources(sources),
        evidence_key="genotype-model-" + model.binding.source_record_id,
        evidence_refs=refs, source_values=sources,
    )


def _population_charges(generation, constructed, alleles):
    matches = tuple(item for item in generation.population_assessments
                    if item.hypothesis == constructed.hypothesis)
    if not matches:
        return ()
    if len(matches) != 1:
        raise EvidenceFusionError("replayed population assessment is ambiguous")
    assessment, = matches
    return population_assessment_charges(assessment, alleles, assessment.policy_sha256)


def population_assessment_charges(assessment, alleles, policy_sha256):
    """Original source charge for already evaluated member population results."""
    if assessment.population_contribution is None:
        return ()
    sources = _unique(
        item.source
        for member in assessment.members
        for item in alleles[member.allele_id].population or ()
    )
    if not sources:
        raise EvidenceFusionError("population contribution has no source evidence")
    return (EvidenceCharge(
        "population", _lineage_for_sources(sources),
        "population-policy-" + _hash((policy_sha256, assessment.members)),
        assessment.population_contribution, None, (), sources,
    ),)


def _clinical_charges(constructed, policy):
    return clinical_member_charges(constructed.members, policy)


def clinical_member_charges(members, policy):
    charges = []
    for member in members:
        bases = tuple(member.clinical)
        conflicted = {item.source_ref for item in bases
                      if item.self_conflict or item.conflict_refs or item.positive and item.negative}
        pending = set(conflicted)
        by_ref = {item.source_ref: item for item in bases}
        while pending:
            component, queue = set(), [min(pending, key=_canonical)]
            while queue:
                ref = queue.pop()
                if ref in component:
                    continue
                component.add(ref)
                basis = by_ref[ref]
                neighbours = set(basis.conflict_refs) | {
                    peer.source_ref for peer in bases if ref in peer.conflict_refs
                }
                queue.extend((neighbours & conflicted) - component)
            pending -= component
            refs = tuple(sorted(component, key=_canonical))
            values = tuple(by_ref[ref].source for ref in refs)
            charges.extend(_mapped_charge(
                policy, "clinical.assertion.conflict",
                lineage_id=_lineage_for_refs(refs),
                evidence_key="clinical-conflict-" + _hash(refs),
                member_id=member.allele_id, evidence_refs=refs,
                source_values=values,
            ))
        groups = {}
        for basis in bases:
            if basis.source_ref in conflicted:
                continue
            applicable = basis.applicability in ("exact", "equivalent") and not basis.superseded_by
            state = ("positive" if applicable and basis.positive else
                     "negative" if applicable and basis.negative else "unresolved")
            groups.setdefault((basis.dependency_group, state), []).append(basis)
        for (dependency_group, state), grouped in sorted(groups.items()):
            source_refs = _unique(item.source_ref for item in grouped)
            refs = _unique(
                ref
                for item in grouped
                for ref in (item.source_ref, *item.relation_refs, *item.superseded_by)
            )
            charges.extend(_mapped_charge(
                policy, "clinical.assertion." + state,
                lineage_id=_lineage_for_refs(source_refs),
                evidence_key="clinical-assertion-group-" + _hash((
                    dependency_group, state, source_refs,
                )),
                member_id=member.allele_id,
                evidence_refs=refs,
                source_values=tuple(item.source for item in grouped),
            ))
    return tuple(charges)


def _technical_charges(constructed, policy, alleles):
    charges = list(technical_member_charges(constructed.model.members, policy, alleles))
    status = constructed.hypothesis.enumeration_status.value
    if status in ("truncated", "failed"):
        fact = "technical.enumeration." + status
        sources = constructed.hypothesis.enumeration_evidence
        charges.extend(_mapped_charge(
            policy, fact, lineage_id=_lineage_for_sources(sources),
            evidence_key="enumeration-" + status,
            source_values=sources,
        ))
    return tuple(charges)


def technical_member_charges(members, policy, alleles):
    charges = []
    for member in members:
        if not member.genotype.read_validation_required:
            continue
        allele = alleles[member.allele.allele_id]
        sources = tuple(item.source for item in allele.calls or ())
        charges.extend(_mapped_charge(
            policy, "technical.read_validation_required",
            lineage_id=_lineage_for_allele(allele),
            evidence_key="read-validation-required",
            member_id=allele.allele_id, evidence_refs=allele.evidence_refs,
            source_values=sources,
        ))
    return tuple(charges)


def _adapt_source_evidence(
    generation: HypothesisGeneration, document: EvidenceDocument, *, policy: FusionPolicy,
    ic_store: HpoICStore, profile_evidence_id: str | None,
    _ic_store_sha256: str | None = None,
):
    _bind_document(generation, document)
    if _ic_store_sha256 is None:
        ic_store_sha256 = hashlib.sha256(hpo_ic_as_bytes(ic_store)).hexdigest()
    else:
        _sha(_ic_store_sha256)
        ic_store_sha256 = _ic_store_sha256
    profile = generation.admission.profiles
    if profile is None:
        if profile_evidence_id is not None or document.disease_profiles:
            raise EvidenceFusionError("unavailable phenotype profile requires a null evidence ID")
        profile_store = None
    else:
        if (type(profile_evidence_id) is not str or profile_evidence_id != profile.evidence_id
                or document.disease_profiles != (profile,)):
            raise EvidenceFusionError("phenotype profile differs from replayed admission")
        profile_store = profile.store
    scores = score_primary_phenotypes(
        document.scoring_view(), document.hypotheses,
        ic_store=ic_store, profile_store=profile_store,
    )
    alleles = {item.allele_id: item for item in document.alleles}
    evidence, disease_audits = [], {}
    for constructed in generation.hypotheses:
        hypothesis = constructed.hypothesis
        score = scores[hypothesis.hypothesis_id]
        profile_refs = () if profile is None else tuple(EvidenceReference(
            profile.evidence_id, item["annotation_id"], profile.public_manifest_sha256,
        ) for item in score.ledger["profile_annotations"])
        disease, audits = _disease_gene_charges(generation, constructed, policy)
        disease_audits[hypothesis.hypothesis_id] = audits
        charges = (
            *disease,
            *_admission_charges(generation, constructed, policy, alleles),
            *_requirement_charges(constructed, policy, alleles),
            *_inheritance_charges(constructed, policy),
            *_population_charges(generation, constructed, alleles),
            *_clinical_charges(constructed, policy),
            *_technical_charges(constructed, policy, alleles),
        )
        for charge in charges:
            for ref in charge.evidence_refs:
                resolve_evidence_reference(document, ref)
        evidence.append(HypothesisEvidence(
            hypothesis, score, tuple(charges), _unique(profile_refs), None,
        ))
    audit = {
        "generation_binding": generation.binding,
        "generation_sha256": _hash(generation),
        "document_sha256": _hash(document),
        "phenotype_scores_sha256": _hash(scores),
        "hpo_ic_store_sha256": ic_store_sha256,
        "hpo_ic_source_binding": ic_store.source_binding,
        "profile_evidence_id": profile_evidence_id,
        "phenotype_profile": "unavailable" if profile is None else "source_bound",
        "disease_gene_aggregations": disease_audits,
        "population_assessments": {
            item.hypothesis.hypothesis_id: item for item in generation.population_assessments
        },
        "qualification_assessments": {
            item.hypothesis.hypothesis_id: item.members for item in generation.hypotheses
        },
        "hypothesis_limits": {
            item.hypothesis.hypothesis_id: item.limits for item in generation.hypotheses
        },
        "source_fact_mappings": policy.source_facts,
    }
    return tuple(evidence), scores, audit


@dataclass(frozen=True, slots=True, repr=False)
class AdaptedContextEvidence:
    evidence: tuple[HypothesisEvidence, ...]
    provenance: ContextProvenance


def adapt_bounded_context_evidence(
    *,
    run_manifest_sha256: str,
    observations: tuple[Observation, ...],
    scoring_view: ScoringView,
    profile: DiseaseProfile | None,
    disease_gene_evidence: tuple[DiseaseGeneEvidence, ...],
    alleles: tuple[AlleleEvidence, ...],
    hypotheses: tuple[ConstructedHypothesis, ...],
    population_assessments: tuple,
    admission_seed,
    admission_decisions: tuple,
    public_sources: tuple[PublicSourcedEvidence, ...],
    clinical_public_sources: tuple[PublicSourcedEvidence, ...],
    policy: FusionPolicy,
    ic_store: HpoICStore,
    profile_evidence_id: str | None,
    prepared_scores=None,
) -> AdaptedContextEvidence:
    """Adapt one complete bounded context through the shared score kernels."""

    try:
        _sha(run_manifest_sha256)
        observations = _typed(observations, tuple[Observation, ...])
        scoring_view = _typed(scoring_view, ScoringView)
        alleles = _typed(alleles, tuple[AlleleEvidence, ...])
        hypotheses = _typed(
            hypotheses, tuple[ConstructedHypothesis, ...],
        )
        genes = _typed(
            disease_gene_evidence, tuple[DiseaseGeneEvidence, ...],
        )
        public_sources = _typed(
            public_sources, tuple[PublicSourcedEvidence, ...],
        )
        clinical_public_sources = _typed(
            clinical_public_sources, tuple[PublicSourcedEvidence, ...],
        )
        if (
            tuple(value.scoring for value in observations)
            != scoring_view.observations
            or {
                value.scoring.observation_id: value.context_provenance
                for value in observations
            }
            != dict(scoring_view.context_provenance)
        ):
            raise ValueError
        input_profile_store = None if profile is None else profile.store
        if profile is None:
            if profile_evidence_id is not None:
                raise ValueError
            profiles = ()
            profile_store = None
        else:
            profile = _typed(profile, DiseaseProfile)
            if profile_evidence_id != profile.evidence_id:
                raise ValueError
            profiles = (profile,)
            profile_store = profile.store
        if prepared_scores is not None:
            if prepared_scores.profile_store is not input_profile_store:
                raise ValueError
            # _typed authenticates a canonical copy. Score against the exact
            # already-bound owner retained by the case-local preparation.
            profile_store = input_profile_store
        decisions = {
            value.decision_id: value for value in admission_decisions
        }
        if (
            len(decisions) != len(admission_decisions)
            or set(decisions) != set(admission_seed.decision_ids)
        ):
            raise ValueError
        public = {}
        for value in (*public_sources, *clinical_public_sources):
            previous = public.get(value.evidence_id)
            if previous is None:
                public[value.evidence_id] = value
            elif previous.public_manifest_sha256 != (
                value.public_manifest_sha256
            ):
                raise ValueError
            else:
                records = dict(previous.records)
                for key, record in value.records.items():
                    if key in records and records[key] != record:
                        raise ValueError
                    records[key] = record
                public[value.evidence_id] = PublicSourcedEvidence(
                    value.evidence_id,
                    value.public_manifest_sha256,
                    MappingProxyType(records),
                )
        document = EvidenceDocument(
            run_manifest_sha256,
            observations=observations,
            contradictions=scoring_view.contradictions,
            disease_profiles=profiles,
            disease_gene_evidence=genes,
            alleles=alleles,
            hypotheses=tuple(value.hypothesis for value in hypotheses),
            public_sourced_evidence=tuple(
                public[key] for key in sorted(public)
            ),
        )
        if document.scoring_view() != scoring_view:
            raise ValueError
        generation = SimpleNamespace(
            admission=SimpleNamespace(
                contexts=(admission_seed,),
                decisions=tuple(decisions[key] for key in sorted(decisions)),
                disease_gene_evidence=genes,
                profiles=profile,
            ),
            population_assessments=population_assessments,
        )
        scores = score_primary_phenotypes(
            scoring_view,
            document.hypotheses,
            ic_store=ic_store,
            profile_store=profile_store,
            prepared_scores=prepared_scores,
        )
        by_id = {value.allele_id: value for value in alleles}
        evidence = []
        for constructed in hypotheses:
            hypothesis = constructed.hypothesis
            score = scores[hypothesis.hypothesis_id]
            profile_refs = () if profile is None else tuple(
                EvidenceReference(
                    profile.evidence_id,
                    item["annotation_id"],
                    profile.public_manifest_sha256,
                )
                for item in score.ledger["profile_annotations"]
            )
            disease, _audits = _disease_gene_charges(
                generation, constructed, policy,
            )
            charges = (
                *disease,
                *_admission_charges(
                    generation, constructed, policy, by_id,
                ),
                *_requirement_charges(constructed, policy, by_id),
                *_inheritance_charges(constructed, policy),
                *_population_charges(generation, constructed, by_id),
                *_clinical_charges(constructed, policy),
                *_technical_charges(constructed, policy, by_id),
            )
            for charge in charges:
                for reference in charge.evidence_refs:
                    resolve_evidence_reference(document, reference)
            evidence.append(HypothesisEvidence(
                hypothesis,
                score,
                tuple(charges),
                _unique(profile_refs),
                None,
            ))
        return AdaptedContextEvidence(tuple(evidence), ContextProvenance(
            observations, scoring_view.contradictions, admission_seed,
            tuple(decisions[key] for key in sorted(decisions)),
            public_records(document, (alleles, hypotheses, evidence,
                                      admission_seed, admission_decisions)),
        ))
    except EvidenceFusionError:
        raise
    except Exception:
        raise EvidenceFusionError(
            "bounded context source evidence could not be adapted"
        ) from None


def _canonical_evidence(value: HypothesisEvidence) -> HypothesisEvidence:
    value = _typed(value, HypothesisEvidence)
    return replace(
        value,
        charges=canonical_evidence_charges(value.charges),
        phenotype_evidence_refs=_unique(value.phenotype_evidence_refs),
        phenotype_source_values=(None if value.phenotype_source_values is None
                                 else _unique(value.phenotype_source_values)),
    )


def canonical_evidence_charges(values):
    """The single atomic ownership/conflict closure for score source adapters."""
    values = _typed(values, tuple[EvidenceCharge, ...])
    charges = {}
    for original in values:
        charge = replace(
            original,
            evidence_refs=_unique(original.evidence_refs),
            source_values=_unique(original.source_values),
        )
        old = charges.get(charge.ownership_key)
        if old is not None and old != charge:
            raise EvidenceFusionError("one atomic evidence item has multiple score owners or values")
        charges[charge.ownership_key] = charge
    return tuple(charges[key] for key in sorted(charges))


def _sum(values: tuple[EvidenceCharge, ...]) -> Fraction | None:
    available = tuple(_exact(item.value) for item in values if item.value is not None)
    return sum(available, Fraction(0)) if available else None


def _raw_component(
    settings: ComponentPolicy, charges: tuple[EvidenceCharge, ...], members: tuple[str, ...],
) -> tuple[Fraction | None, Mapping[str, Fraction | None], tuple[EvidenceCharge, ...]]:
    selected = tuple(item for item in charges if item.component == settings.name)
    shared = tuple(item for item in selected if item.member_id is None)
    member_charges = {member: tuple(item for item in selected if item.member_id == member) for member in members}
    if (settings.aggregation == "shared" and any(item.member_id is not None for item in selected)
            or settings.aggregation == "member_mean" and shared):
        raise EvidenceFusionError("evidence charge scope differs from component aggregation")
    shared_raw = _sum(shared)
    member_raw = MappingProxyType({member: _sum(member_charges[member]) for member in members})
    available = shared_raw is not None or any(value is not None for value in member_raw.values())
    if not available:
        return None, member_raw, selected
    if settings.aggregation == "shared":
        return shared_raw, member_raw, selected
    mean = sum((value if value is not None else _exact(settings.raw_neutral)
                for value in member_raw.values()), Fraction(0)) / len(members)
    if settings.aggregation == "member_mean":
        return mean, member_raw, selected
    return (shared_raw or Fraction(0)) + mean, member_raw, selected


def component_numeric_value(settings: ComponentPolicy, raw: Fraction | None, *,
                            member_raw: Mapping[str, Fraction | None] | None = None) -> Fraction | None:
    """The shared exact transform, including neutral-aware member means."""
    if raw is None:
        normalized = None
    elif settings.aggregation == "member_mean":
        if not member_raw:
            raise EvidenceFusionError("member aggregation requires member values")
        normalized = sum((
            settings.transform.apply(
                value if value is not None else _exact(settings.raw_neutral)
            )
            for value in member_raw.values()
        ), Fraction(0)) / len(member_raw)
    else:
        normalized = settings.transform.apply(raw)
    return normalized


def _component(
    settings: ComponentPolicy, raw: Fraction | None, charges: tuple[EvidenceCharge, ...], *,
    member_raw: Mapping[str, Fraction | None] | None = None,
    adjustments: tuple[ScoreAdjustment, ...] = (), extra_sources: tuple[SourcedRecord, ...] = (),
    extra_refs: tuple[EvidenceReference, ...] = (),
) -> ScoreComponent:
    normalized = component_numeric_value(settings,raw,member_raw=member_raw)
    sources = _unique((*extra_sources, *(source for item in charges for source in item.source_values)))
    refs = _unique((*extra_refs, *(ref for item in charges for ref in item.evidence_refs)))
    return ScoreComponent(
        settings.name,
        None if raw is None else str(raw),
        None if normalized is None else str(normalized),
        settings.neutral_value,
        settings.weight,
        TRANSFORM + ":" + settings.name,
        settings.transform_sha256,
        settings.parameters,
        refs,
        None if not charges and not extra_sources else sources,
        adjustments,
    )


def generated_penalty_values(*, completeness, phase, paired, policy, molecular, member_raw):
    """Evaluate the original penalty amounts and identity-aware weakest member.

    Source proofs are bound separately after this exact numeric decision.
    Missing member values remain neutral and cannot imply observed weakness.
    """
    generated = []
    if completeness is Completeness.INCOMPLETE:
        generated.append(("incompleteness",_exact(policy.penalties.incompleteness),None,None))
    if paired:
        phase_amount = {
            HypothesisPhase.UNKNOWN: policy.penalties.phase_unknown,
            HypothesisPhase.CONFLICTING: policy.penalties.phase_conflicting,
            HypothesisPhase.IN_CIS: policy.penalties.phase_in_cis,
        }.get(phase)
        if phase_amount is not None:
            generated.append(("phase-"+phase.value,_exact(phase_amount),None,None))
        # Missing molecular evidence is neutral/unknown, not an observed weak
        # partner. The member mean still exposes the missing neutral contribution.
        if all(raw is not None for raw in member_raw.values()):
            weakest_id, weakest_value = min((
                (member, molecular.transform.apply(raw))
                for member, raw in member_raw.items()
            ), key=lambda item: (item[1], item[0]))
            shortfall = max(Fraction(0), _exact(policy.penalties.weak_partner_floor) - weakest_value)
            amount = min(_exact(policy.penalties.weak_partner_cap),
                         shortfall * _exact(policy.penalties.weak_partner_scale))
            if amount:
                generated.append(("weak-partner",amount,weakest_id,weakest_value))
    return tuple(generated)


def _penalty_adjustments(
    value: HypothesisEvidence, policy: FusionPolicy, molecular: ComponentPolicy,
    member_raw: Mapping[str, Fraction | None], molecular_charges: tuple[EvidenceCharge, ...],
) -> tuple[tuple[str, Fraction, str, tuple[EvidenceReference, ...], tuple[SourcedRecord, ...]], ...]:
    hypothesis = value.hypothesis
    generated = []
    for kind,amount,weakest_id,weakest_value in generated_penalty_values(completeness=hypothesis.completeness,
            phase=hypothesis.phase,paired=len(hypothesis.allele_refs)==2,policy=policy,molecular=molecular,member_raw=member_raw):
        if kind == "incompleteness":
            generated.append((kind,amount,"incomplete genotype under the declared hypothesis requirements",
                              (),hypothesis.completeness_evidence))
        elif kind.startswith("phase-"):
            generated.append((kind,amount,"pair phase is "+hypothesis.phase.value,(),hypothesis.phase_evidence))
        else:
            member_evidence = tuple(item for item in molecular_charges if item.member_id == weakest_id)
            generated.append((
                    "weak-partner-" + _hash((hypothesis.hypothesis_id, weakest_id))[:16], amount,
                    "weakest molecular member=" + weakest_id + ";normalized=" + str(weakest_value)
                    + ";floor=" + policy.penalties.weak_partner_floor,
                    _unique(ref for item in member_evidence for ref in item.evidence_refs),
                    _unique(source for item in member_evidence for source in item.source_values),
            ))
    return tuple(generated)


def _tie_key(rule: TieRule, hypothesis: GenotypeHypothesis, components: Mapping[str, ScoreComponent]) -> TieKey:
    return scientific_tie_key(rule,context=hypothesis.context,completeness=hypothesis.completeness,
                              phase=hypothesis.phase,hypothesis_id=hypothesis.hypothesis_id,components=components)


def scientific_tie_key(rule, *, context, completeness, phase, hypothesis_id, components):
    """Shared original tie values over already evaluated scientific outcomes."""
    if rule.field in _COMPONENTS:
        component = components[rule.field]
        value = component.normalized if component.normalized is not None else component.neutral_value
    elif rule.field == "completeness":
        value = _COMPLETENESS[completeness]
    elif rule.field == "phase":
        value = _PHASE[phase]
    else:
        value = {
            "native_disease_id": context.native_disease_id,
            "disease_id": context.disease_id,
            "gene_id": context.gene_id,
            "locus_id": context.locus_id,
            "model": context.model,
            "hypothesis_id": hypothesis_id,
        }[rule.field]
    return TieKey(rule.field, value, rule.direction)


def technical_penalty_values(settings, members, base_raw, charges, generated):
    """Exact technical additions, including each original saturation refusal."""
    running_raw = _exact(settings.raw_neutral)
    requests = []
    for charge in charges:
        if charge.value is None:
            requests.append(None)
            continue
        requested = _exact(charge.value)
        if charge.member_id is not None:
            requested /= len(members)
        before = settings.transform.apply(running_raw)
        running_raw += requested
        after = settings.transform.apply(running_raw)
        if after - before != requested:
            raise EvidenceFusionError("technical penalty would saturate its declared transform")
        requests.append(requested)
    if base_raw is not None and running_raw != base_raw:
        raise EvidenceFusionError("technical penalty aggregation is not replayable")
    for requested in generated:
        before = settings.transform.apply(running_raw)
        running_raw += requested
        after = settings.transform.apply(running_raw)
        if after - before != requested:
            raise EvidenceFusionError("generated penalty would saturate its declared transform")
    return running_raw if base_raw is not None or generated else None, tuple(requests)


def total_score_value(components):
    """Shared exact weighted sum of (component name, weight, effective value)."""
    return sum(((-1 if name == "technical_and_conflict_penalties" else 1)
                * _exact(weight) * _exact(value) for name, weight, value in components), Fraction(0))


def _score(value: HypothesisEvidence, policy: FusionPolicy, policy_sha256: str) -> ScoreLedger:
    value = _canonical_evidence(value)
    hypothesis = value.hypothesis
    members = tuple(sorted(hypothesis.evidence_subject_ids))
    phenotype = replace(
        policy.phenotype.component(value.phenotype_score),
        evidence_refs=value.phenotype_evidence_refs,
        source_values=value.phenotype_source_values,
    )
    components = [phenotype]
    member_values = {}
    selected_charges = {}
    for settings in policy.components:
        if settings.name == "technical_and_conflict_penalties":
            continue
        raw, member_raw, charges = _raw_component(settings, value.charges, members)
        components.append(_component(settings, raw, charges, member_raw=member_raw))
        member_values[settings.name] = member_raw
        selected_charges[settings.name] = charges

    molecular = next(item for item in policy.components if item.name == "molecular")
    technical = next(item for item in policy.components if item.name == "technical_and_conflict_penalties")
    base_raw, _, technical_charges = _raw_component(technical, value.charges, members)
    generated = _penalty_adjustments(
        value, policy, molecular, member_values["molecular"], selected_charges["molecular"],
    )
    penalty_raw, requests = technical_penalty_values(technical, members, base_raw, technical_charges,
                                                    tuple(requested for _, requested, *_ in generated))
    adjustments = []
    for charge, requested in zip(technical_charges, requests, strict=True):
        if requested:
            adjustments.append(ScoreAdjustment(
                "fusion-evidence-" + _hash(charge.ownership_key)[:16],
                str(requested),
                charge.kind + " evidence charged once from " + charge.lineage_id + "/" + charge.evidence_key,
                charge.evidence_refs,
            ))
    for key, requested, reason, refs, _ in generated:
        if requested:
            adjustments.append(ScoreAdjustment(
                "fusion-" + key,
                str(requested),
                reason,
                refs,
            ))
    components.append(_component(
        technical,
        penalty_raw,
        technical_charges,
        adjustments=tuple(adjustments),
        extra_sources=_unique(source for *_, sources in generated for source in sources),
        extra_refs=_unique(ref for *_, refs, _sources in generated for ref in refs),
    ))
    by_name = {item.name: item for item in components}
    ordered_components = tuple(by_name[name] for name in _COMPONENTS)
    total = total_score_value((item.name, item.weight,
        item.normalized if item.normalized is not None else item.neutral_value) for item in ordered_components)
    provenance = {
        "schema": SCHEMA,
        "executable_version": EXECUTABLE_VERSION,
        "policy_sha256": policy_sha256,
        "derivation_sha256": policy.derivation_sha256,
        "derivation_scope": policy.derivation_scope,
        "policy": _plain(policy),
        "release_status": "development-policy-input; task-9.6-freeze-required",
        "evidence_ownership": "lineage+evidence-key; one-component-and-scope-only",
    }
    return ScoreLedger(
        hypothesis.hypothesis_id,
        ordered_components,
        str(total),
        tuple(_tie_key(rule, hypothesis, by_name) for rule in policy.ties),
        policy.policy_id,
        policy_sha256,
        provenance,
    )


def _rank_adapted_hypotheses(
    evidence: tuple[HypothesisEvidence, ...], *, policy_bytes: bytes, policy_sha256: str,
    derivation_bytes: bytes, expected_derivation_sha256: str,
    source_replay: Mapping[str, object], expected_result: FusionRanking | None = None,
) -> FusionRanking:
    """Internal exact kernel; public callers enter through source replay below."""
    try:
        policy = policy_from_json(policy_bytes, expected_sha256=policy_sha256)
        _sha(expected_derivation_sha256)
        if (type(derivation_bytes) is not bytes
                or hashlib.sha256(derivation_bytes).hexdigest() != expected_derivation_sha256
                or policy.derivation_sha256 != expected_derivation_sha256):
            raise ValueError
        if not isinstance(source_replay, Mapping) or not source_replay:
            raise ValueError
        evidence = tuple(sorted(
            (_canonical_evidence(item) for item in _typed(evidence, tuple[HypothesisEvidence, ...])),
            key=lambda item: item.hypothesis.hypothesis_id,
        ))
        _ids(tuple(item.hypothesis.hypothesis_id for item in evidence))
        ledgers = tuple(_score(item, policy, policy_sha256) for item in evidence)
        ties = ranking_ties(policy)
        order = order_ledgers(ledgers, ties)
        source_replay = _plain(source_replay)
        input_sha256 = _hash({"evidence": evidence, "source_replay": source_replay})
        result = FusionRanking(
            policy.policy_id,
            policy_sha256,
            input_sha256,
            ledgers,
            order,
            ties,
            {
                "schema": SCHEMA,
                "executable_version": EXECUTABLE_VERSION,
                "input_sha256": input_sha256,
                "policy_sha256": policy_sha256,
                "derivation_sha256": expected_derivation_sha256,
                "derivation_bytes_verified": True,
                "candidate_membership": "source-replayed-fixed-generation",
                "candidate_set_normalization": "forbidden",
                "source_replay": source_replay,
                "atomic_evidence_charges": {
                    item.hypothesis.hypothesis_id: _plain(item.charges) for item in evidence
                },
                "controlled_use_ready": False,
            },
        )
        if expected_result is not None and result.as_bytes() != _typed(expected_result, FusionRanking).as_bytes():
            raise ValueError
        return result
    except _ERRORS:
        raise EvidenceFusionError("evidence fusion replay failed") from None


def rank_hypotheses(
    document: EvidenceDocument, *, replay: HypothesisReplay, ic_store: HpoICStore,
    profile_evidence_id: str | None, policy_bytes: bytes, policy_sha256: str,
    derivation_bytes: bytes, expected_derivation_sha256: str,
    expected_result: FusionRanking | None = None,
) -> FusionRanking:
    """Replay sources, derive the fixed score facts, and order every hypothesis.

    HypothesisEvidence/EvidenceCharge are internal adapter records. This boundary
    does not accept either type from callers and independently regenerates all
    candidate membership, score facts, phenotype values, and identity lineage.
    """
    try:
        return _prepared_fusion(
            document, replay=replay, ic_store=ic_store,
            profile_evidence_id=profile_evidence_id,
            policy_bytes=policy_bytes, policy_sha256=policy_sha256,
            derivation_bytes=derivation_bytes,
            expected_derivation_sha256=expected_derivation_sha256,
            expected_result=expected_result,
        ).ranking
    except _ERRORS:
        raise EvidenceFusionError("source-bound evidence fusion replay failed") from None
