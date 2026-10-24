"""Closed hypothesis population policy, not annotation, admission, or readiness.

The pure kernel authenticates no source. consume_population_policy replays both
annotation adapters against independently supplied pins before invoking it.
Public exception curation is configuration, never inferred biological evidence.
"""

from __future__ import annotations

from collections.abc import Mapping
from dataclasses import InitVar, dataclass, fields, replace
from enum import StrEnum
from fractions import Fraction
import hashlib
from types import MappingProxyType
from typing import Literal

from genosapiens import clinical_annotation as clinical, population_annotation as population
from genosapiens.admission_links import PinnedPublicSource
from genosapiens.allele_identity import AlleleIdentity
from genosapiens.evidence_model import (
    AlleleEvidence, GenotypeHypothesis, HypothesisContext, SourcedRecord,
    _Controlled, _exact, _hints, _sha, _typed, _validation_session,
)
from genosapiens.evidence_serialization import _canonical, _loads
from genosapiens.gnomad_projection import EXPECTED_REGION_FLAGS
from genosapiens.manifest import PublicResourceManifest


SCHEMA = "genosapiens.population-policy"
RESULT_SCHEMA = "genosapiens.population-policy-result"
CURATION_SCHEMA = "genosapiens.population-exception-curation"
Category = Literal["curated_pathogenic", "founder", "reduced_penetrance", "technical_uncertainty"]
Signal = Literal[
    "counts_missing", "inconsistent_ac_an", "inconsistent_genotype_counts", "zero_an", "low_an", "source_filtered",
    "filters_unknown", "discordant", "discordance_unknown", "not_callable",
    "callability_unknown", "region_flagged", "region_unknown",
]
_CATEGORIES = ("curated_pathogenic", "founder", "reduced_penetrance", "technical_uncertainty")
_REGIONS = tuple("REGION_" + flag.upper() for flag in EXPECTED_REGION_FLAGS)
# Full source vocabularies from gnomad-joint-projection.json. Maxima are
# publisher-native selections, not restricted to the six retained ancestry rows.
_FREQUENCY_ANCESTRIES = {
    subset: (*population._ANCESTRIES, "asj", "fin", "remaining", *(() if subset == "exomes" else ("ami",)))
    for subset in population.GNOMAD_SUBSETS
}
_FAF_ANCESTRIES = {
    subset: tuple(a for a in population._ANCESTRIES if subset != "genomes" or a != "mid")
    for subset in population.GNOMAD_SUBSETS
}


class PopulationPolicyError(ValueError):
    """Static diagnostics only; no paths, alleles, source rows, or chained errors."""


class PopulationOutcome(StrEnum):
    NOT_POPULATION_EXCLUDED = "not_population_excluded"
    POPULATION_EXCLUDED = "population_excluded"


@dataclass(frozen=True, slots=True, repr=False)
class ContextKey(_Controlled):
    native_disease_id: str
    gene_id: str | None
    locus_id: str | None
    inheritance: str
    model: str

    def __post_init__(self):
        _Controlled.__post_init__(self)
        if (not all(v.strip() for v in (self.native_disease_id, self.inheritance, self.model))
                or not (self.gene_id or self.locus_id)
                or any(v is not None and not v.strip() for v in (self.gene_id, self.locus_id))):
            raise PopulationPolicyError("invalid population context key")

    @classmethod
    def from_context(cls, context: HypothesisContext) -> ContextKey | None:
        if context.native_disease_id is None or context.inheritance is None:
            return None
        return cls(context.native_disease_id, context.gene_id, context.locus_id,
                   context.inheritance, context.model)


@dataclass(frozen=True, slots=True, repr=False)
class FrequencyRule(_Controlled):
    rule_id: str
    subset: Literal["exomes", "genomes", "joint"]
    stratum: Literal["adjusted", "ancestry", "grpmax", "fafmax95", "fafmax99"]
    ancestry: str | None
    metric: Literal["AF", "FAF95", "FAF99"]
    threshold: str
    minimum_an: int
    allowed_filters: tuple[str, ...]
    region_flags: tuple[str, ...]
    unreliable_signals: tuple[Signal, ...]

    def __post_init__(self):
        _Controlled.__post_init__(self)
        ancestries = (_FREQUENCY_ANCESTRIES[self.subset] if self.stratum == "grpmax" else
                      _FAF_ANCESTRIES[self.subset] if self.stratum in ("fafmax95", "fafmax99") else population._ANCESTRIES)
        if (not self.rule_id.strip() or not 0 < _exact(self.threshold) <= 1
                or self.minimum_an < 1
                or self.stratum == "adjusted" and self.ancestry is not None
                or self.stratum == "ancestry" and self.ancestry is None
                or self.ancestry is not None and self.ancestry not in ancestries
                or not {"inconsistent_ac_an", "inconsistent_genotype_counts", "zero_an", "low_an"} <= set(self.unreliable_signals)
                or not set(self.region_flags) <= set(_REGIONS)
                or any(len(set(v)) != len(v) for v in (
                    self.allowed_filters, self.region_flags, self.unreliable_signals,
                )) or any(not v.strip() for v in self.allowed_filters)):
            raise PopulationPolicyError("invalid population frequency rule")


@dataclass(frozen=True, slots=True, repr=False)
class PenaltyTransform(_Controlled):
    kind: Literal["capped_linear_excess", "categorical_excess"]
    scale: str
    cap: str

    def __post_init__(self):
        _Controlled.__post_init__(self)
        if _exact(self.scale) <= 0 or not 0 < _exact(self.cap) <= 1:
            raise PopulationPolicyError("invalid population penalty transform")


@dataclass(frozen=True, slots=True, repr=False)
class RescueRules(_Controlled):
    pathogenic_significances: tuple[str, ...]
    pathogenic_review_statuses: tuple[str, ...]
    pathogenic_conflicts: Literal["block_assertion_rescue", "retain_without_adjudication"]
    pathogenic_match: Literal["exact_native_condition_and_gene_target"]
    technical_signals: tuple[Signal, ...]
    technical_scope: Literal["every_reliable_excess_record_or_every_unreliable_excess_record"]
    founder: Literal["exact_public_curated_exception"]
    reduced_penetrance: Literal["exact_public_curated_exception"]
    combination: Literal["any_applies"]

    def __post_init__(self):
        _Controlled.__post_init__(self)
        for values in (self.pathogenic_significances, self.pathogenic_review_statuses, self.technical_signals):
            if not values or len(set(values)) != len(values) or any(not v.strip() for v in values):
                raise PopulationPolicyError("invalid population rescue rules")
        if any(v.lower().replace("_", " ") not in (
            "pathogenic", "likely pathogenic", "pathogenic/likely pathogenic",
        ) for v in self.pathogenic_significances):
            raise PopulationPolicyError("nonpathogenic significance in pathogenic rescue rule")


@dataclass(frozen=True, slots=True, repr=False)
class ContextPolicy(_Controlled):
    context: ContextKey
    rules: tuple[FrequencyRule, ...]
    transform: PenaltyTransform
    rescues: RescueRules
    record_aggregation: Literal["maximum"]
    member_aggregation: Literal["maximum"]
    comparison: Literal["strict_greater_than"]
    frequency_arithmetic: Literal["AF_exact_AC_over_AN;FAF_decimal_typed_value"]
    unreliable_frequency: Literal["neutral"]
    unknown_frequency: Literal["neutral_not_excluded"]
    neutral_penalty: Literal["0"]
    exclusion: Literal["any_reliable_excess_member_without_rescue"]

    def __post_init__(self):
        _Controlled.__post_init__(self)
        selectors = [(r.subset, r.stratum, r.ancestry, r.metric) for r in self.rules]
        if (not self.rules or len({r.rule_id for r in self.rules}) != len(self.rules)
                or len(set(selectors)) != len(selectors)):
            raise PopulationPolicyError("duplicate or missing population rules")


@dataclass(frozen=True, slots=True, repr=False)
class CuratedException(_Controlled):
    """Compiled PUBLIC policy entry, reproduced from an independent curation record."""

    category: Literal["founder", "reduced_penetrance"]
    native_allele: AlleleIdentity
    reference_manifest_sha256: str
    context: ContextKey
    source: SourcedRecord

    def __post_init__(self):
        _Controlled.__post_init__(self)
        _sha(self.reference_manifest_sha256)
        if self.source.public_manifest_sha256 is None or self.source.source_fields is None:
            raise PopulationPolicyError("public exception source binding required")


@dataclass(frozen=True, slots=True, repr=False)
class PopulationPolicy(_Controlled):
    policy_id: str
    contexts: tuple[ContextPolicy, ...]
    exceptions: tuple[CuratedException, ...]

    def __post_init__(self):
        _Controlled.__post_init__(self)
        if (not self.policy_id.strip()
                or len({c.context for c in self.contexts}) != len(self.contexts)
                or len({_canonical(e) for e in self.exceptions}) != len(self.exceptions)
                or any(e.context not in {c.context for c in self.contexts} for e in self.exceptions)):
            raise PopulationPolicyError("invalid population policy")


@dataclass(frozen=True, slots=True, repr=False)
class PreparedPopulationPolicy:
    """One owner's immutable policy/index; never a source-authentication token."""

    policy: PopulationPolicy
    policy_sha256: str
    contexts_by_key: Mapping[ContextKey, ContextPolicy]

    def __post_init__(self):
        try:
            # Keep the validated immutable copy. All results made by this owner
            # share it; validating and then retaining the caller's object would
            # not establish that ownership boundary.
            policy = _typed(self.policy, PopulationPolicy)
            _sha(self.policy_sha256)
            if not isinstance(self.contexts_by_key, Mapping):
                raise ValueError
            index = dict(self.contexts_by_key)
            if index != {context.context: context for context in policy.contexts}:
                raise ValueError
            object.__setattr__(self, "policy", policy)
            object.__setattr__(self, "contexts_by_key", MappingProxyType(
                {context.context: context for context in policy.contexts}))
        except Exception:
            raise PopulationPolicyError("invalid prepared population policy") from None


def prepare_population_policy(policy_bytes: bytes, *, expected_sha256: str) -> PreparedPopulationPolicy:
    """Authenticate actual policy bytes once; the owner separately validates curation.

    Complete embedded-policy result validation/serialization remains mandatory.
    """
    policy = policy_from_json(policy_bytes, expected_sha256=expected_sha256)
    return PreparedPopulationPolicy(policy, expected_sha256,
        {context.context: context for context in policy.contexts})


@dataclass(frozen=True, slots=True, repr=False)
class FrequencyAssessment(_Controlled):
    rule_id: str
    population_index: int | None
    frequency: str | None
    quality_signals: tuple[Signal, ...]
    reliable: bool
    threshold: str
    excess_ratio: str | None
    transformed_penalty: str | None
    applied_penalty: str | None
    status: Literal["unknown", "unreliable", "within_threshold", "excess"]


@dataclass(frozen=True, slots=True, repr=False)
class ClinicalAssessment(_Controlled):
    clinical_index: int
    status: Literal["qualifying", "no_match", "unknown"]
    reason: str


@dataclass(frozen=True, slots=True, repr=False)
class RescueAssessment(_Controlled):
    category: Category
    status: Literal["applies", "no_match", "unknown", "not_assessed"]
    clinical_indices: tuple[int, ...]
    exception_indices: tuple[int, ...]
    frequency_indices: tuple[int, ...]
    reason: str


@dataclass(frozen=True, slots=True, repr=False)
class MemberAssessment(_Controlled):
    allele_id: str
    native_identity_id: str
    status: Literal["not_assessed", "unknown", "partial", "assessed"]
    frequencies: tuple[FrequencyAssessment, ...]
    clinical: tuple[ClinicalAssessment, ...]
    rescues: tuple[RescueAssessment, ...]
    penalty: str | None
    outcome: PopulationOutcome


@dataclass(frozen=True, slots=True, repr=False)
class PopulationPolicyResult(_Controlled):
    hypothesis: GenotypeHypothesis
    alleles: tuple[AlleleEvidence, ...]
    policy: PopulationPolicy
    policy_sha256: str
    selected_context: ContextKey | None
    members: tuple[MemberAssessment, ...]
    population_penalty: str | None
    population_contribution: str | None
    outcome: PopulationOutcome
    reason: str
    _prepared: InitVar[PreparedPopulationPolicy | None] = None

    def __post_init__(self, _prepared):
        if _prepared is None:
            # External construction and decoding still validate the entire
            # logical result. The owner is never part of the wire contract.
            _Controlled.__post_init__(self)
            return
        if (type(_prepared) is not PreparedPopulationPolicy
                or self.policy is not _prepared.policy
                or self.policy_sha256 != _prepared.policy_sha256):
            raise PopulationPolicyError("population result policy owner changed")
        with _validation_session():
            for field in fields(self):
                if field.name != "policy":
                    object.__setattr__(self, field.name, _typed(
                        getattr(self, field.name), _hints(type(self))[field.name]))


def policy_as_json(policy: PopulationPolicy) -> bytes:
    try:
        return (_canonical({"schema": SCHEMA, "policy": _typed(policy, PopulationPolicy)}) + "\n").encode()
    except Exception:
        raise PopulationPolicyError("cannot encode population policy") from None


def policy_from_json(content: bytes, *, expected_sha256: str) -> PopulationPolicy:
    """Hash actual supplied bytes. The independent pin is not acceptance or curation."""
    try:
        _sha(expected_sha256)
        if type(content) is not bytes or hashlib.sha256(content).hexdigest() != expected_sha256:
            raise ValueError
        data = _loads(content)
        if set(data) != {"schema", "policy"} or data["schema"] != SCHEMA:
            raise ValueError
        return _typed(data["policy"], PopulationPolicy, decode=True)
    except Exception:
        raise PopulationPolicyError("invalid population policy bytes or binding") from None


def result_as_json(result: PopulationPolicyResult) -> bytes:
    try:
        return (_canonical({"schema": RESULT_SCHEMA, "classification": "controlled", "llm_policy": "no-llm",
                            "result": _typed(result, PopulationPolicyResult)}) + "\n").encode()
    except Exception:
        raise PopulationPolicyError("cannot encode population policy result") from None


def result_from_json(content: bytes) -> PopulationPolicyResult:
    """Structural decoding only, never a source-authentication token."""
    try:
        data = _loads(content)
        if set(data) != {"schema", "classification", "llm_policy", "result"} or (
            data["schema"], data["classification"], data["llm_policy"]
        ) != (RESULT_SCHEMA, "controlled", "no-llm"):
            raise ValueError
        return _typed(data["result"], PopulationPolicyResult, decode=True)
    except Exception:
        raise PopulationPolicyError("invalid population policy result JSON") from None


def evaluate_population_members(key: ContextKey | None, alleles: tuple[AlleleEvidence, ...], *,
                                prepared_policy: PreparedPopulationPolicy) -> tuple[MemberAssessment, ...]:
    """Evaluate complete member population inputs through the shared kernel.

    Context/model policy, raw source array order and exception indices are the
    actual dependencies. Phase, completeness and other members cannot alter an
    individual assessment. Combining their penalties/outcomes is separate.
    """
    if type(prepared_policy) is not PreparedPopulationPolicy:
        raise PopulationPolicyError("population policy owner is invalid")
    policy = prepared_policy.policy
    settings = prepared_policy.contexts_by_key.get(key)
    return evaluate_population_member_inputs(key, alleles, settings=settings,
        indexed_exceptions=tuple(enumerate(policy.exceptions)))


def evaluate_population_member_inputs(key, alleles, *, settings, indexed_exceptions):
    """Original member kernel with explicit context and indexed exception inputs.

    Exception ordinals refer to the complete policy, including when an owner
    supplies only this context's relevant entries. They are never renumbered.
    Source authentication and complete policy closure belong to the caller.
    """
    key = _typed(key, ContextKey | None)
    alleles = _typed(alleles, tuple[AlleleEvidence, ...])
    settings = _typed(settings, ContextPolicy | None)
    indexed_exceptions = _typed(indexed_exceptions, tuple[tuple[int, CuratedException], ...])
    if (settings is not None and settings.context != key
            or any(i < 0 for i, _ in indexed_exceptions)
            or tuple(i for i, _ in indexed_exceptions) != tuple(sorted({i for i, _ in indexed_exceptions}))):
        raise PopulationPolicyError("population member policy inputs are invalid")
    members = []
    for allele in alleles:
        if settings is None:
            members.append(MemberAssessment(allele.allele_id, allele.native_identity_id, "not_assessed", (), (),
                tuple(RescueAssessment(c, "not_assessed", (), (), (), "no_exact_context_policy") for c in _CATEGORIES),
                None, PopulationOutcome.NOT_POPULATION_EXCLUDED))
            continue
        rows = []
        for rule in sorted(settings.rules, key=lambda r: r.rule_id):
            matches = [(i, p) for i, p in enumerate(allele.population or ()) if (
                p.subset, p.stratum, p.ancestry
            ) == (rule.subset, rule.stratum, rule.ancestry)]
            if not matches:
                rows.append(FrequencyAssessment(rule.rule_id, None, None, (), False, rule.threshold,
                                                None, None, None, "unknown"))
            for index, record in matches:
                signals = set()
                ac, an = record.allele_count, record.allele_number
                if ac is None or an is None:
                    signals.add("counts_missing")
                if ac is not None and an is not None and ac > an:
                    signals.add("inconsistent_ac_an")
                # Sum only known contributions. Missing raw counts stay null;
                # even this lower bound cannot exceed either known AC or AN.
                genotype_lower_bound = sum(weight * count for weight, count in (
                    (2, record.homozygote_count), (1, record.hemizygote_count),
                ) if count is not None)
                if any(genotype_lower_bound > count for count in (ac, an) if count is not None):
                    signals.add("inconsistent_genotype_counts")
                if an == 0:
                    signals.add("zero_an")
                if an is not None and an < rule.minimum_an:
                    signals.add("low_an")
                if record.filters is None or None in record.filters:
                    signals.add("filters_unknown")
                if any(f is not None and f not in rule.allowed_filters for f in record.filters or ()):
                    signals.add("source_filtered")
                if record.discordant is None:
                    signals.add("discordance_unknown")
                elif record.discordant:
                    signals.add("discordant")
                if record.callable is None:
                    signals.add("callability_unknown")
                elif not record.callable:
                    signals.add("not_callable")
                raw_context = (record.source.source_fields or {}).get("raw_context", {})
                for flag in rule.region_flags:
                    value = raw_context.get(flag)
                    if value == "1":
                        signals.add("region_flagged")
                    elif value != "0":
                        signals.add("region_unknown")
                frequency = None
                if rule.metric == "AF":
                    if record.allele_frequency is not None and ac is not None and an is not None and 0 <= ac <= an and an > 0:
                        frequency = _exact(str(ac)) / an
                else:
                    value = (record.filtering_allele_frequency if rule.metric == "FAF95"
                             else record.filtering_allele_frequency_99)
                    if value is not None:
                        # str(float) may be scientific notation; Fraction accepts it.
                        frequency = Fraction(str(value))
                reliable = frequency is not None and not signals.intersection(rule.unreliable_signals)
                excess = None if frequency is None else max(0, frequency / _exact(rule.threshold) - 1)
                penalty = None if excess is None else min(_exact(settings.transform.cap), _exact(settings.transform.scale) * (
                    excess if settings.transform.kind == "capped_linear_excess" else int(excess > 0)
                ))
                rows.append(FrequencyAssessment(rule.rule_id, index, None if frequency is None else str(frequency),
                    tuple(sorted(signals)), reliable, rule.threshold, None if excess is None else str(excess),
                    None if penalty is None else str(penalty), str(penalty) if reliable else None,
                    "unknown" if frequency is None else "unreliable" if not reliable else "excess" if excess > 0 else "within_threshold"))

        clinical_rows = []
        rescue = settings.rescues
        for index, record in enumerate(allele.clinical or ()):
            if (record.explicit_condition_link is not True or not record.assertion_id
                    or not record.condition_ids or record.target_id is None or key.gene_id is None
                    or key.locus_id is not None):
                status, reason = "unknown", "individual_native_condition_and_gene_target_required"
            elif record.condition_ids != (key.native_disease_id,) or record.target_id != key.gene_id:
                status, reason = "no_match", "different_native_condition_or_target"
            elif not record.significance or record.review_status is None or any(s is None for s in record.significance):
                status, reason = "unknown", "significance_or_review_unavailable"
            elif (not set(record.significance) <= set(rescue.pathogenic_significances)
                    or record.review_status not in rescue.pathogenic_review_statuses):
                status, reason = "no_match", "configured_significance_or_review_not_met"
            elif record.conflict_refs and rescue.pathogenic_conflicts == "block_assertion_rescue":
                status, reason = "no_match", "conflict_blocks_this_rescue_without_adjudication"
            else:
                status, reason = "qualifying", "individual_exact_native_assertion_meets_policy"
            clinical_rows.append(ClinicalAssessment(index, status, reason))
        qualifying = tuple(c.clinical_index for c in clinical_rows if c.status == "qualifying")
        rescues = [RescueAssessment("curated_pathogenic", "applies" if qualifying else
            "unknown" if not clinical_rows or any(c.status == "unknown" for c in clinical_rows) else "no_match",
            qualifying, (), (), "assertions_are_boolean_rescue_only;no_clinical_or_conflict_score")]
        native = next(i for i in allele.identities if i.identity_id == allele.native_identity_id)
        for category in ("founder", "reduced_penetrance"):
            matched = tuple(i for i, entry in indexed_exceptions if entry.category == category
                and entry.context == key and entry.native_allele == native.identity
                and entry.reference_manifest_sha256 == native.reference_manifest_sha256)
            rescues.append(RescueAssessment(category, "applies" if matched else "unknown", (), matched, (),
                "exact_public_curated_entry" if matched else "no_exact_public_curated_entry;not_negative_evidence"))
        excess_rows = tuple(i for i, row in enumerate(rows) if row.reliable and _exact(row.excess_ratio) > 0)
        technical_rows = excess_rows or tuple(i for i, row in enumerate(rows)
                                             if row.excess_ratio is not None and _exact(row.excess_ratio) > 0)
        technical_matches = tuple(i for i in technical_rows if set(rows[i].quality_signals).intersection(rescue.technical_signals))
        rescues.append(RescueAssessment("technical_uncertainty",
            "applies" if technical_rows and technical_matches == technical_rows else
            "unknown" if not any(r.frequency is not None for r in rows) else "no_match",
            (), (), technical_matches, "source_uncertainty_only;mandatory_read_checks_unchanged"))
        penalties = [_exact(row.applied_penalty) for row in rows if row.applied_penalty is not None]
        excluded = bool(excess_rows) and not any(r.status == "applies" for r in rescues)
        members.append(MemberAssessment(allele.allele_id, allele.native_identity_id,
            "unknown" if not penalties else "assessed" if all(r.reliable for r in rows) else "partial",
            tuple(rows), tuple(clinical_rows), tuple(rescues), str(max(penalties)) if penalties else None,
            PopulationOutcome.POPULATION_EXCLUDED if excluded else PopulationOutcome.NOT_POPULATION_EXCLUDED))
    return tuple(members)


def summarize_population_members(members: tuple[MemberAssessment, ...]):
    """Combine existing member assessments with the original maximum/any rules."""
    return summarize_population_values(tuple((m.penalty,m.outcome) for m in members))


def summarize_population_values(values):
    penalties = [_exact(penalty) for penalty,_ in values if penalty is not None]
    penalty = max(penalties) if penalties else None
    excluded = any(outcome is PopulationOutcome.POPULATION_EXCLUDED for _,outcome in values)
    return (None if penalty is None else str(penalty),None if penalty is None else str(-penalty),
            PopulationOutcome.POPULATION_EXCLUDED if excluded else PopulationOutcome.NOT_POPULATION_EXCLUDED)


def evaluate_population_policy(
    hypothesis: GenotypeHypothesis, alleles: tuple[AlleleEvidence, ...], *,
    prepared_policy: PreparedPopulationPolicy,
) -> PopulationPolicyResult:
    """Pure deterministic kernel. Inputs/curation are NOT authenticated here.

    Every raw record survives in alleles. No call, consequence, optional mapping,
    phase, or read-validation label changes the population-only outcome.
    """
    try:
        if type(prepared_policy) is not PreparedPopulationPolicy:
            raise ValueError
        policy, policy_sha256 = prepared_policy.policy, prepared_policy.policy_sha256
        hypothesis = _typed(hypothesis, GenotypeHypothesis)
        alleles = _typed(alleles, tuple[AlleleEvidence, ...])
        by_id = {a.allele_id: a for a in alleles}
        if (len(by_id) != len(alleles) or set(by_id) != set(hypothesis.evidence_subject_ids)
                or any(by_id[r.allele_id].native_identity_id != r.identity_id for r in hypothesis.evidence_subject_refs)):
            raise ValueError
        # Canonical unordered member order; source array order is never changed.
        hypothesis = replace(hypothesis, allele_refs=tuple(sorted(hypothesis.allele_refs, key=lambda r: r.allele_id)))
        alleles = tuple(by_id[k] for k in sorted(by_id))
        key = ContextKey.from_context(hypothesis.context)
        settings = prepared_policy.contexts_by_key.get(key)
        members = evaluate_population_members(key, alleles, prepared_policy=prepared_policy)
        penalty,contribution,outcome = summarize_population_members(members)
        return PopulationPolicyResult(hypothesis, alleles, policy, policy_sha256, key if settings else None,
            tuple(members), penalty,contribution,outcome,
            "no_exact_context_policy" if settings is None else "population_only;not_all_checks_or_scientific_acceptance",
            _prepared=prepared_policy)
    except Exception:
        raise PopulationPolicyError("population policy evaluation failed") from None


def validate_population_curation(
    policy: PopulationPolicy, public_curation_bytes: Mapping[str, bytes],
) -> None:
    """Replay every selected exception against its independent public source bytes.

    Shared by single-hypothesis consumption and whole-universe construction so
    source validation need not be repeated for every member combination. This
    verifies curation bindings, not empirical acceptance of the selected policy.
    """
    try:
        policy = _typed(policy, PopulationPolicy)
        public_curation_bytes = dict(public_curation_bytes)
        required_bytes = set()
        for entry in policy.exceptions:
            source = entry.source
            for digest in (source.source_sha256, source.public_manifest_sha256):
                content = public_curation_bytes[digest]
                if type(content) is not bytes or hashlib.sha256(content).hexdigest() != digest:
                    raise ValueError
                required_bytes.add(digest)
            raw = _loads(public_curation_bytes[source.source_sha256])
            manifest = _typed(_loads(public_curation_bytes[source.public_manifest_sha256]), PublicResourceManifest, decode=True)
            if (set(raw) != {"schema", "classification", "release", "record_id", "claim", "citations"}
                    or raw["schema"] != CURATION_SCHEMA or raw["classification"] != "public"
                    or raw["release"] != source.source_release or raw["record_id"] != source.source_record_id
                    or (manifest.resource_id, manifest.release) != (source.source_dataset, source.source_release)
                    or manifest.assembly != "GRCh38" or manifest.transformation.schema != CURATION_SCHEMA
                    or _canonical(raw) != _canonical(source.source_fields)
                    or _canonical(raw["claim"]) != _canonical({
                        "category": entry.category, "native_allele": entry.native_allele,
                        "reference_manifest_sha256": entry.reference_manifest_sha256, "context": entry.context,
                    }) or type(raw["citations"]) is not list or not raw["citations"]
                    or any(type(c) is not str or not c.startswith("https://") for c in raw["citations"])
                    or not any(o.sha256 == source.source_sha256 and o.size_bytes == len(public_curation_bytes[source.source_sha256])
                               and o.retention.policy == "durable" and o.retention.delete_after is None for o in manifest.outputs)):
                raise ValueError
        if set(public_curation_bytes) != required_bytes:
            raise ValueError
    except Exception:
        raise PopulationPolicyError("population exception curation replay failed") from None


def consume_population_policy(
    hypothesis: GenotypeHypothesis, *, policy_bytes: bytes, policy_sha256: str,
    population_result: population.PopulationAnnotationResult,
    population_requests: tuple[population.PopulationRequest, ...],
    joint: population.GnomadJointReference, all_sites: population.GnomadAllSitesReference,
    population_reference: population.PopulationReferenceCompatibility,
    runtime_image: str, limits: population.PopulationQueryLimits,
    clinical_result: clinical.ClinicalAnnotationResult,
    clinvar: PinnedPublicSource | None, open_targets: PinnedPublicSource | None,
    clinical_reference: PinnedPublicSource,
    public_curation_bytes: Mapping[str, bytes],
    expected_result: PopulationPolicyResult | None = None,
) -> PopulationPolicyResult:
    """Replay actual annotations and public curation, then recompute the kernel.

    The supported chain is population -> clinical. Original population must be
    empty: the append-only adapter cannot authenticate caller-carried records.
    All independent resource/policy pins are supplied by the caller's frozen stack,
    not recovered from reports. expected_result enables serialized-result replay.
    No result constructor, source field, or codec can replace this boundary.
    """
    try:
        prepared = prepare_population_policy(policy_bytes, expected_sha256=policy_sha256)
        if any(r.allele.population for r in population_requests):
            raise ValueError
        population.validate_population_annotation(population_result, population_requests, joint=joint,
            all_sites=all_sites, reference=population_reference, runtime_image=runtime_image, limits=limits)
        clinical.verify_clinical_annotation(clinical_result, population_result.alleles,
            clinvar=clinvar, open_targets=open_targets, reference=clinical_reference)
        validate_population_curation(prepared.policy, public_curation_bytes)
        # Replay covers the entire batch; evaluate precisely the requested members.
        selected = tuple(a for a in clinical_result.alleles if a.allele_id in hypothesis.evidence_subject_ids)
        result = evaluate_population_policy(hypothesis, selected, prepared_policy=prepared)
        if expected_result is not None and result != _typed(expected_result, PopulationPolicyResult):
            raise ValueError
        return result
    except Exception:
        raise PopulationPolicyError("population policy source replay failed") from None
