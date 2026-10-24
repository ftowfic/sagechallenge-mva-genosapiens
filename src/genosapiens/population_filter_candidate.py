"""Conservative, inactive candidate for early population work elimination.

This pure evaluator consumes already authenticated source projections. It does
not authenticate resources, replace the active population score/penalty kernel,
or authorize removal from a published scientific domain. The experimental
driver keeps frequency opportunity separate from complete exclusion decisions.
"""
from dataclasses import dataclass
from enum import Enum
from fractions import Fraction
from typing import Mapping, Sequence

from genosapiens.common_allele_filter import (
    ClinicalProtectionFacts, CommonAlleleFilterPolicy, PrefilterOccurrence,
    _protection_reasons, prefilter_bypass_reasons,
)
from genosapiens.gnomad_faf95_adapter import float32_interval
from genosapiens.population_policy import ContextKey

POLICY_ID = "track1-conservative-population-candidate-v1"
THRESHOLD = Fraction(1, 20)
MINIMUM_AN = 2000
ANCESTRIES = ("afr", "amr", "eas", "nfe", "sas")
WORK_POLICY_ID = "track1-protected-population-work-filter-v1"
WORK_POLICY_PATH = "config/track1/protected-population-work-filter-v1.json"
WORK_POLICY_SHA256 = "3c0d96e7621507498eeda8fc9ddd964fa2626b7f27a8b1e7eef4434d6c5d719e"


@dataclass(frozen=True, slots=True)
class PopulationWorkPolicy:
    """Pinned integration policy, distinct from experimental decision authority.

    Loading the policy does not qualify a release or establish occurrence,
    model, resource or protection coverage. The context owner must supply those
    proofs before using the unchanged candidate predicate to omit work.
    """

    content: bytes

    def __post_init__(self):
        import hashlib
        if (type(self.content) is not bytes
                or hashlib.sha256(self.content).hexdigest() != WORK_POLICY_SHA256):
            raise ValueError("protected population work policy binding changed")

    @property
    def policy_id(self):
        return WORK_POLICY_ID

    @property
    def sha256(self):
        return WORK_POLICY_SHA256


def load_population_work_policy(content: bytes, *, expected_sha256: str) -> PopulationWorkPolicy:
    """Accept only the canonical fixed policy; no caller-selected threshold."""
    if expected_sha256 != WORK_POLICY_SHA256:
        raise ValueError("protected population work policy authority changed")
    return PopulationWorkPolicy(content)


@dataclass(frozen=True, slots=True)
class FrequencySupport:
    qualified: bool
    supporting_ancestries: tuple[str, ...]
    genome_only_ancestries: tuple[str, ...]
    reasons: tuple[str, ...]


@dataclass(frozen=True, slots=True)
class CandidateDecision:
    disposition: str
    reasons: tuple[str, ...]
    frequency: FrequencySupport
    policy_id: str = POLICY_ID

    @property
    def excluded(self):
        return self.disposition == "exclude_common"

    @property
    def outcome(self):
        if self.excluded:
            return OptionalFilterOutcome.COMPLETE_EXCLUSION
        if self.disposition == "retain":
            return OptionalFilterOutcome.PROTECTED
        return OptionalFilterOutcome.UNSUPPORTED


class OptionalFilterOutcome(str, Enum):
    COMPLETE_EXCLUSION = "complete_protected_exclusion"
    PROTECTED = "protected_route_retained"
    UNSUPPORTED = "unsupported_optimization_retained"


class PopulationCandidateBindingError(ValueError):
    """Required identity failure. Optional-filter recovery must not catch this."""

    operation = "population_work_filter"
    failure_code = "required_population_binding_mismatch"


class OptionalPopulationFilterError(RuntimeError):
    """Local predicate failure after required source authentication."""

    operation = "population_work_filter"
    failure_code = "optional_population_predicate_failed"


def optional_frequency_support(records):
    # Only the optional pure predicate runs inside this boundary. Resource
    # queries, ownership checks, and source authentication stay outside it.
    try:
        return frequency_support(records)
    except (ValueError, TypeError, ArithmeticError) as error:
        raise OptionalPopulationFilterError("optional population predicate failed") from error


def _counts_usable(record):
    if record is None:
        return False
    ac, an = record.get("allele_count"), record.get("allele_number")
    if type(ac) is not int or type(an) is not int or not 0 <= ac <= an or an <= 0:
        return False
    lower = 0
    for field, weight in (("homozygote_count", 2), ("hemizygote_count", 1)):
        count = record.get(field)
        if count is not None:
            if type(count) is not int or count < 0:
                return False
            lower += weight * count
    return lower <= ac


def _quality_usable(record):
    if not _counts_usable(record):
        return False
    filters = record.get("filters")
    flags = record.get("raw_context", {})
    if (not isinstance(filters, (list, tuple)) or any(f != "PASS" for f in filters)
            or record.get("callable") is not True or not isinstance(flags, Mapping)
            or flags.get("REGION_NOT_CALLED_IN_" + record["subset"].upper()) != "0"):
        return False
    # This is exome interval QC, not a genome QC flag. Missing is not PASS.
    return record["subset"] != "exomes" or flags.get("REGION_FAIL_INTERVAL_QC") == "0"


def _bound_above(record):
    if not _quality_usable(record) or record["allele_number"] < MINIMUM_AN:
        return False
    value = record.get("filtering_allele_frequency_99")
    if type(value) not in (str, float, int):
        return False
    try:
        bound = Fraction(str(value))
        if not 0 <= bound <= Fraction(record["allele_count"], record["allele_number"]):
            return False
        # Reuse the established float32/text rounding guard. This is a storage
        # interval around source FAF99, not an additional statistical interval.
        return bound > THRESHOLD and float32_interval(str(value)).lower_bound > THRESHOLD
    except (ValueError, ArithmeticError):
        return False


def _outside_exome_calling_regions(record):
    if record is None:
        return False
    flags = record.get("raw_context", {})
    return (isinstance(flags, Mapping)
            and all(flags.get(key) == "1" for key in (
                "REGION_NOT_CALLED_IN_EXOMES", "REGION_OUTSIDE_BROAD_CALLING_REGION",
                "REGION_OUTSIDE_UKB_CALLING_REGION", "REGION_OUTSIDE_BROAD_CAPTURE_REGION",
                "REGION_OUTSIDE_UKB_CAPTURE_REGION"))
            and record.get("allele_count") is None and record.get("allele_number") is None
            and record.get("filters") is None and record.get("callable") is None)


def frequency_support(records: Sequence[Mapping]) -> FrequencySupport:
    """Source FAF99 per continental stratum; no pooling or count imputation.

FAF99 is the publisher's per-stratum Poisson confidence bound, not the prior
experiment's multiplicity-adjusted exact-binomial bound or a cohort error rate.
Different exome/genome frequencies can support the same >5% decision. Any
one-sided support with an unresolved comparison vetoes this candidate. Exome
absence is allowed only outside both declared calling/capture regions; this
does not assert zero physical coverage or supply an all-sites witness.
"""
    by = {}
    for row in records:
        if not isinstance(row, Mapping):
            raise ValueError("invalid population candidate projection")
        if row.get("subset") not in ("genomes", "exomes") or row.get("stratum") != "ancestry":
            continue
        ancestry = row.get("ancestry")
        if ancestry not in ANCESTRIES:
            continue
        key = (row["subset"], ancestry)
        if key in by:
            return FrequencySupport(False, (), (), ("duplicate_frequency_selector",))
        by[key] = row
    supported, genome_only, reasons = [], [], []
    for ancestry in ANCESTRIES:
        genome, exome = by.get(("genomes", ancestry)), by.get(("exomes", ancestry))
        gh, eh = _bound_above(genome), _bound_above(exome)
        if gh and eh:
            supported.append(ancestry)
        elif gh and _outside_exome_calling_regions(exome) and genome.get("discordant") is not True:
            supported.append(ancestry)
            genome_only.append(ancestry)
        elif gh or eh:
            reasons.append(ancestry + ":comparison_unresolved_or_contradictory")
    qualified = bool(supported) and not reasons
    if not supported:
        reasons.append("no_quality_qualified_continental_FAF99_above_5pct")
    return FrequencySupport(qualified, tuple(supported), tuple(genome_only), tuple(reasons))


def evaluate_candidate(*, context: ContextKey | None, occurrence: PrefilterOccurrence | None,
                       records: Sequence[Mapping], protection: ClinicalProtectionFacts | None,
                       protection_policy: CommonAlleleFilterPolicy) -> CandidateDecision:
    """Apply existing representation and exception protections before exclusion.

The caller binds occurrence, source records and protection-query facts to the
same exact allele/reference. Empty exception arrays or absent assertions are
never substituted for a completed authenticated protection query.
"""
    if type(protection_policy) is not CommonAlleleFilterPolicy:
        raise ValueError("invalid candidate protection policy")
    support = optional_frequency_support(records)
    if context is None:
        return CandidateDecision("retain_unknown", ("context_applicability_unknown",), support)
    if type(context) is not ContextKey:
        raise ValueError("invalid candidate context")
    if occurrence is None:
        return CandidateDecision("retain_unknown", ("occurrence_protections_unavailable",), support)
    bypass = prefilter_bypass_reasons(occurrence)
    if bypass:
        return CandidateDecision("retain", bypass, support)
    if occurrence.reference_equivalence != "exact":
        return CandidateDecision("retain_unknown", ("reference_equivalence_uncertain",), support)
    if protection is None:
        return CandidateDecision("retain_unknown", ("exception_review_unavailable",), support)
    if type(protection) is not ClinicalProtectionFacts:
        raise ValueError("invalid candidate protection facts")
    status, reasons = _protection_reasons(protection, protection_policy)
    if status != "continue":
        return CandidateDecision(status, tuple(reasons), support)
    if not support.qualified:
        return CandidateDecision("retain_unknown", support.reasons, support)
    return CandidateDecision("exclude_common", ("FAF99_lower_bound_above_5pct",
        "complete_pinned_protection_query_clear"), support)


def removable_shared_classes(consumers, excluded_hypotheses):
    """A shared class is removable only when every recorded consumer is excluded.

Consumers is a complete, caller-authenticated mapping from class IDs to sets
of hypothesis IDs. Empty/partial domains cannot establish removal authority.
"""
    excluded = frozenset(excluded_hypotheses)
    return frozenset(key for key, uses in consumers.items() if uses and frozenset(uses) <= excluded)


def context_rescue_vetoes(assessment, *, protection, protection_policy, clinical_evidence=None):
    """Scheduling interpretation of existing, unchanged rescue assessments.

    The context owner authenticates the assessment's exact model/member and
    the protection query's occurrence/reference/resource bindings. This pure
    helper cannot grant that authority. D41.4a/b discharge only missing curated
    entries and absent/aggregate clinical uncertainty after complete clearance;
    the latter also requires the exact typed clinical inputs. It never
    changes the scientific assessment or its population outcome/penalty.
    """
    from genosapiens.population_policy import MemberAssessment
    if type(assessment) is not MemberAssessment or type(protection_policy) is not CommonAlleleFilterPolicy:
        raise ValueError('invalid context rescue inputs')
    if protection is not None and type(protection) is not ClinicalProtectionFacts:
        raise ValueError('invalid context rescue protection')
    expected = {'curated_pathogenic', 'founder', 'reduced_penetrance', 'technical_uncertainty'}
    if len(assessment.rescues) != len(expected) or {r.category for r in assessment.rescues} != expected:
        raise ValueError('incomplete context rescue coverage')
    clear = protection is not None and _protection_reasons(protection, protection_policy)[0] == 'continue'
    aggregate_only = False
    if clinical_evidence is not None:
        from genosapiens.evidence_model import ClinicalEvidence
        if (type(clinical_evidence) is not tuple
                or any(type(row) is not ClinicalEvidence for row in clinical_evidence)
                or tuple(row.clinical_index for row in assessment.clinical) != tuple(range(len(clinical_evidence)))):
            raise ValueError('incomplete clinical rescue input coverage')
        aggregate_only = all(row.status == 'no_match' or (
            row.status == 'unknown' and row.reason == 'individual_native_condition_and_gene_target_required'
            and clinical_evidence[row.clinical_index].assertion_id is None
            and clinical_evidence[row.clinical_index].explicit_condition_link is not True)
            for row in assessment.clinical)
    reasons = []
    for rescue in assessment.rescues:
        if rescue.status == 'no_match':
            continue
        if (clear and rescue.category in {'founder', 'reduced_penetrance'}
                and rescue.status == 'unknown'
                and rescue.reason == 'no_exact_public_curated_entry;not_negative_evidence'):
            continue
        if (clear and aggregate_only and rescue.category == 'curated_pathogenic'
                and rescue.status == 'unknown'
                and rescue.reason == 'assertions_are_boolean_rescue_only;no_clinical_or_conflict_score'):
            continue
        reasons.append('context_rescue:' + rescue.category + ':' + rescue.status)
    return tuple(reasons)


def evaluate_source_candidate(*, value, routing, context, records,
                              population_allele_id, population_reference_sha256,
                              protection_session, protection_policy,
                              clinical_manifest_sha256, exception_manifest_sha256):
    """Evaluate one authenticated source occurrence without publishing a filter.

    The caller authenticates population records and binds their allele/reference
    arguments, and opens the protection session against the pinned complete
    clinical/exception index. Source routing and selected-call protections use
    the preparation owners. A protection result is consumed only for the query
    issued here; absent source routing is never inferred from annotation fields.
    All source occurrences must exclude before removing a normalized member.
    """
    from genosapiens.wgs_source_routing import (
        SourceOccurrenceFacts, SourceRoutingContext, classify_source_occurrence,
    )
    from genosapiens.wgs_public_filter import (
        _prefilter_occurrence, _clinical_requests, clinical_protection_facts,
    )
    if (type(value) is not SourceOccurrenceFacts or type(routing) is not SourceRoutingContext
            or value.reference_manifest_sha256 != population_reference_sha256):
        raise PopulationCandidateBindingError("candidate source/population reference binding changed")
    if (type(protection_policy) is not CommonAlleleFilterPolicy
            or context is not None and type(context) is not ContextKey):
        raise PopulationCandidateBindingError("candidate policy or context binding changed")
    route = classify_source_occurrence(value, contract=routing.contract, renderer=routing.renderer)
    # Native routing deliberately hides a normalized identifier. Authenticate
    # supplied population evidence against the independent normalization fact.
    normalized = value.normalization.normalized_allele
    expected_allele_id = None if normalized is None else normalized.normalized_allele_id
    if (population_allele_id != expected_allele_id
            and not (population_allele_id is None and route.normalized_allele_id is None)):
        raise PopulationCandidateBindingError("candidate source/population allele binding changed")
    occurrence = _prefilter_occurrence(value, routing.renderer, routing.contract, route=route)
    bypass = prefilter_bypass_reasons(occurrence)
    if bypass:
        return (CandidateDecision("retain", bypass,
            FrequencySupport(False, (), (), ("protected_route_frequency_not_evaluated",))),
            occurrence, None)
    support = optional_frequency_support(records)
    protection = None
    # Preserve source/native/phase/technical guards without issuing a lookup.
    if context is not None and not prefilter_bypass_reasons(occurrence) and support.qualified:
        query, = _clinical_requests(((value, occurrence),), routing.renderer)
        result, = protection_session.query_many((query,))
        protection = clinical_protection_facts(result.status, result.records,
            clinical_manifest_sha256=clinical_manifest_sha256,
            exception_manifest_sha256=exception_manifest_sha256)
    decision = evaluate_candidate(context=context, occurrence=occurrence, records=records,
        protection=protection, protection_policy=protection_policy)
    if (context is not None and not prefilter_bypass_reasons(occurrence)
            and not support.qualified):
        decision = CandidateDecision("retain_unknown", support.reasons, support)
    return decision, occurrence, protection


def evaluate_saved_candidate(*, saved, context, records, population_allele_id,
                             source_manifest_sha256, protection_policy,
                             clinical_manifest_sha256, exception_manifest_sha256):
    """Reuse a caller-authenticated scoped filter receipt and its source binding.

    The completed filter shard and parent source/reference manifest must already
    be authenticated. Only the explicit successful protection-check reason
    restores negative coverage: a bare no-match status does not prove it.
    This does not reuse the old frequency decision or activate domain pruning.
    """
    import hashlib
    import json
    from genosapiens.common_allele_filter import (
        CLINICAL_ROLES, PrefilterDecision, common_allele_filter_policy_as_json,
    )
    from genosapiens.evidence_model import _typed
    from genosapiens.wgs_filter import FILTER_EVIDENCE_SCHEMA, ledger_filter_decision
    from genosapiens.wgs_ledger import FilterDecision
    from genosapiens.track1_early_scope import EarlyScopeDecision
    if (type(saved) is not FilterDecision
            or saved.parent_source_manifest_sha256 != source_manifest_sha256
            or saved.normalized_allele_id != population_allele_id
            or saved.policy_sha256 != hashlib.sha256(
                common_allele_filter_policy_as_json(protection_policy)).hexdigest()):
        raise ValueError("candidate saved protection binding changed")
    evidence = json.loads(saved.evidence_json)
    if (set(evidence) != {"schema", "version", "classification", "llm_policy", "occurrence", "decision",
            "scope_decision", "resource_sha256", "population_lookup", "clinical_lookup", "request_counts"}
            or evidence["schema"] != FILTER_EVIDENCE_SCHEMA or evidence["version"] != 2
            or evidence["classification"] != "controlled" or evidence["llm_policy"] != "no-llm"):
        raise ValueError("candidate saved protection envelope changed")
    scope = _typed(evidence["scope_decision"], EarlyScopeDecision, decode=True)
    if scope.status != 'in_scope':
        raise ValueError("candidate saved protection is outside the admitted scope")
    occurrence = _typed(evidence["occurrence"], PrefilterOccurrence, decode=True)
    old = _typed(evidence["decision"], PrefilterDecision, decode=True)
    if (old.clinical_manifest_sha256 != clinical_manifest_sha256
            or old.exception_manifest_sha256 != exception_manifest_sha256
            or ledger_filter_decision(occurrence, old,
                scope_decision=scope, **evidence["request_counts"],
                parent_source_manifest_sha256=source_manifest_sha256) != saved):
        raise ValueError("candidate saved protection source tuple changed")
    proof = old.proof
    clear = "no_protection_found_in_pinned_sources" in old.reason_codes
    protection = ClinicalProtectionFacts(proof.protection_query_status,
        proof.protected_significance_tokens, proof.conflict_tokens,
        proof.curated_exact_match, proof.curated_interval_match,
        "complete" if clear else "uncertain", CLINICAL_ROLES,
        clinical_manifest_sha256, exception_manifest_sha256, proof.assertion_scopes)
    if clear and _protection_reasons(protection, protection_policy)[0] != "continue":
        raise ValueError("candidate saved protection clearance contradicts proof")
    return evaluate_candidate(context=context, occurrence=occurrence, records=records,
        protection=protection, protection_policy=protection_policy), occurrence, protection
