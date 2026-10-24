"""Public Track 1 experiment, evaluation, and pre-held-out transition.

Ranking execution receives only the hash-bound runner view.  Evaluator truth is
joined afterwards by this coordinator and never enters the bounded scientific worker.
The records here are validation-only: even an accepted public candidate is not
the later controlled-use base release.
"""

from __future__ import annotations

from collections.abc import Mapping
from dataclasses import asdict, dataclass
from fractions import Fraction
import hashlib
import json
import re
from typing import Literal

from genosapiens.evidence_model import EnumerationStatus
from genosapiens.wgs_public_execution import BoundedPublicExecution
from genosapiens.wgs_validation_projection import ValidationRecoveryBinding, TargetRecovery
from genosapiens.track1_policy_bindings import PolicyBindings
from genosapiens.track1_roster import REQUIRED_MODEL_COVERAGE, RosterCase, Track1Roster
from genosapiens.track1_validation_plan import (
    AcceptanceCriterion,
    InputKind,
    RankerCaseView,
    RankerInput,
    RankerResource,
    ValidationPlan,
    validation_plan_as_json,
)


SCHEMA = "genosapiens.track1-public-validation-experiment"
CANDIDATE_SCHEMA = "genosapiens.track1-pre-held-out-candidate"
ACCEPTANCE_SCHEMA = "genosapiens.track1-consolidated-public-acceptance"
SHA256 = re.compile(r"^[0-9a-f]{64}$")
IDENTIFIER = re.compile(r"^[A-Za-z0-9][A-Za-z0-9._:-]*$")
AttemptStatus = Literal["succeeded", "failed", "timed_out", "unsupported", "truncated"]
REQUIRED_ACCEPTANCE_RESULTS = (
    "authorized-evidence-and-publication-contracts",
    "bounded-real-tool-vcf-first",
    "current-contract-engineering",
    "evaluator-and-numerical-policy",
    "giab-genomic-representation",
    "public-caller-preparation",
    "public-resource-license-and-usability",
    "security-and-data-boundaries",
)
MODEL_COVERAGE_KEYS = (
    "clinical-annotation-completed",
    "clinical-unasserted-exploratory-pair",
    "completed-population-lookup-unknown-rescue",
    "independent-non-population-contradiction",
    "population-annotation-completed",
    "population-excluded-unresolved-partner",
    "unknown-count-hypothesis",
    "unresolved-partner-singleton-violation",
)
POPULATION_LOOKUP_KEYS = tuple(sorted((
    "requested_alleles", "supported_alleles", "unique_loci",
    "bcftools_invocations", "tabix_invocations", "matched", "no_match",
    "unsupported_native_class", "no_authoritative_normalized_identity",
    "ambiguous_sequence", "mitochondrial_outside_joint_product",
    "contig_outside_joint_product",
)))


class Track1ValidationRunnerError(ValueError):
    """Fail-closed validation diagnostics without case-derived payloads."""


def _bytes(value: object) -> bytes:
    return (
        json.dumps(
            value,
            sort_keys=True,
            separators=(",", ":"),
            ensure_ascii=True,
            allow_nan=False,
        )
        + "\n"
    ).encode("ascii")


def _hash(value: bytes) -> str:
    return hashlib.sha256(value).hexdigest()


def _fraction_text(value: Fraction) -> str:
    return (
        str(value.numerator)
        if value.denominator == 1
        else f"{value.numerator}/{value.denominator}"
    )


def _sha(value: object) -> None:
    if type(value) is not str or SHA256.fullmatch(value) is None:
        raise Track1ValidationRunnerError("invalid validation binding")


def _identifier(value: object) -> None:
    if type(value) is not str or IDENTIFIER.fullmatch(value) is None:
        raise Track1ValidationRunnerError("invalid validation identifier")


@dataclass(frozen=True, slots=True)
class PolicyCandidate:
    policy_id: str
    bindings: PolicyBindings
    fixed_settings: tuple[str, ...]
    variable_settings: tuple[str, ...]

    def __post_init__(self) -> None:
        _identifier(self.policy_id)
        if (
            type(self.bindings) is not PolicyBindings
            or type(self.fixed_settings) is not tuple
            or type(self.variable_settings) is not tuple
            or not self.fixed_settings
            or any(type(item) is not str or not item for item in (*self.fixed_settings, *self.variable_settings))
            or self.fixed_settings != tuple(sorted(set(self.fixed_settings)))
            or self.variable_settings != tuple(sorted(set(self.variable_settings)))
            or set(self.fixed_settings) & set(self.variable_settings)
        ):
            raise Track1ValidationRunnerError("invalid finite policy candidate")


@dataclass(frozen=True, slots=True)
class ExperimentProtocol:
    experiment_id: str
    policy_family_id: str
    objective_metrics: tuple[str, ...]
    deterministic_selection_ties: tuple[str, ...]
    maximum_policy_candidates: int
    maximum_attempts_per_case: int
    stopping_condition: Literal["complete-declared-family-or-budget-exhausted"]
    retry_aggregation: Literal["last-attempt-is-terminal-all-attempts-retained"]

    def __post_init__(self) -> None:
        _identifier(self.experiment_id)
        _identifier(self.policy_family_id)
        if (
            self.objective_metrics
            != (
                "top_10_full_genotype_recovery",
                "top_1_full_genotype_recovery",
                "exact_genotype_recovery",
                "candidate_burden",
            )
            or self.deterministic_selection_ties
            != ("lower_candidate_burden", "lexical_policy_id")
            or type(self.maximum_policy_candidates) is not int
            or self.maximum_policy_candidates < 1
            or type(self.maximum_attempts_per_case) is not int
            or self.maximum_attempts_per_case < 1
            or self.stopping_condition
            != "complete-declared-family-or-budget-exhausted"
            or self.retry_aggregation
            != "last-attempt-is-terminal-all-attempts-retained"
        ):
            raise Track1ValidationRunnerError("invalid finite experiment protocol")


@dataclass(frozen=True, slots=True)
class ResourceMeasurements:
    wall_seconds: str | None
    cpu_seconds: str | None
    peak_memory_bytes: int | None
    peak_disk_bytes: int | None

    def __post_init__(self) -> None:
        for value in (self.wall_seconds, self.cpu_seconds):
            if value is not None:
                try:
                    if type(value) is not str or Fraction(value) < 0:
                        raise ValueError
                except (ValueError, ZeroDivisionError, TypeError):
                    raise Track1ValidationRunnerError("invalid attempt measurement") from None
        for value in (self.peak_memory_bytes, self.peak_disk_bytes):
            if value is not None and (type(value) is not int or value < 0):
                raise Track1ValidationRunnerError("invalid attempt measurement")


PREPARATION_STAGE_PREFIX = "preparation:"
PREPARATION_STAGES = frozenset({
    "preparation", "materialization", "input-descriptors", "call-preparation",
    "inventory", "annotation", "annotation-replay", "phenotype-ingestion",
    "variant-first-resource-query", "scientific-policy-binding",
    "hypothesis-replay-binding", "resource-topology", "manifest-qualification",
    "candidate-admission", "population-policy-binding",
})
EXECUTION_STAGE_PREFIX = "execution:"
EXECUTION_STAGES = frozenset({
    "execution", "candidate-authentication", "hypothesis-generation",
    "evidence-assembly", "ranking-sensitivity-and-draft", "bounded-scientific-execution",
})


@dataclass(frozen=True, slots=True)
class ValidationAttempt:
    case_id: str
    split: Literal["development", "held_out"]
    scope: str
    policy_id: str
    policy_bindings_sha256: str
    attempt_number: int
    status: AttemptStatus
    execution_sha256: str | None
    input_sha256s: tuple[str, ...]
    resource_descriptor_sha256s: tuple[str, ...]
    stage_statuses: tuple[tuple[str, str], ...]
    recovery_binding: ValidationRecoveryBinding | None
    eligible_member_groups: tuple[tuple[str, ...], ...]
    candidate_burden: int | None
    hypothesis_count: int | None
    required_enumerations: int | None
    completed_enumerations: int | None
    measurements: ResourceMeasurements
    eligible_epcr: tuple[str, ...] = ()
    association_state_counts: tuple[tuple[str, int], ...] = ()
    admission_outcome_counts: tuple[tuple[str, int], ...] = ()
    model_coverage_counts: tuple[tuple[str, int], ...] = ()
    clinical_only_eligible_member_groups: tuple[tuple[str, ...], ...] = ()
    source_inventory_status: Literal["complete", "not_assessed"] = "not_assessed"
    broader_unresolved_source_count: int | None = None
    broader_unresolved_relevance_source_count: int | None = None
    population_lookup_counts: tuple[tuple[str, int], ...] = ()

    @property
    def preparation_stage(self) -> str | None:
        return next((name.removeprefix(PREPARATION_STAGE_PREFIX)
                     for name, _status in self.stage_statuses
                     if name.startswith(PREPARATION_STAGE_PREFIX)), None)

    @property
    def execution_stage(self) -> str | None:
        return next((name.removeprefix(EXECUTION_STAGE_PREFIX)
                     for name, _status in self.stage_statuses
                     if name.startswith(EXECUTION_STAGE_PREFIX)), None)

    def __post_init__(self) -> None:
        if self.source_inventory_status not in {"complete", "not_assessed"} or any(
            value is not None and (type(value) is not int or value < 0)
            for value in (self.broader_unresolved_source_count, self.broader_unresolved_relevance_source_count)
        ):
            raise Track1ValidationRunnerError("invalid validation source coverage")
        for value in (self.case_id, self.scope, self.policy_id):
            _identifier(value)
        _sha(self.policy_bindings_sha256)
        if self.execution_sha256 is not None:
            _sha(self.execution_sha256)
        for value in (*self.input_sha256s, *self.resource_descriptor_sha256s):
            _sha(value)
        try:
            epcr = tuple(Fraction(value) for value in self.eligible_epcr)
        except (ValueError, ZeroDivisionError, TypeError):
            raise Track1ValidationRunnerError(
                "invalid validation EPCR values"
            ) from None
        groups = (
            self.eligible_member_groups,
            self.clinical_only_eligible_member_groups,
        )
        if (
            self.split not in {"development", "held_out"}
            or self.status not in {
                "succeeded", "failed", "timed_out", "unsupported", "truncated",
            }
            or type(self.attempt_number) is not int
            or self.attempt_number < 1
            or type(self.measurements) is not ResourceMeasurements
            or type(self.stage_statuses) is not tuple
            or not self.stage_statuses
            or any(
                type(item) is not tuple
                or len(item) != 2
                or type(item[0]) is not str
                or not item[0]
                or item[1] not in {"succeeded", "succeeded-empty", "failed"}
                for item in self.stage_statuses
            )
            or len({item[0] for item in self.stage_statuses}) != len(self.stage_statuses)
            or sum(name.startswith(PREPARATION_STAGE_PREFIX) for name, _ in self.stage_statuses) > 1
            or sum(name.startswith(EXECUTION_STAGE_PREFIX) for name, _ in self.stage_statuses) > 1
            or (any(name.startswith(PREPARATION_STAGE_PREFIX) for name, _ in self.stage_statuses)
                and any(name.startswith(EXECUTION_STAGE_PREFIX) for name, _ in self.stage_statuses))
            or any(
                name.startswith(PREPARATION_STAGE_PREFIX) and (
                    name.removeprefix(PREPARATION_STAGE_PREFIX) not in PREPARATION_STAGES
                    or status != "failed" or self.status != "failed"
                    or self.execution_sha256 is not None
                )
                for name, status in self.stage_statuses
            )
            or any(
                name.startswith(EXECUTION_STAGE_PREFIX) and (
                    name.removeprefix(EXECUTION_STAGE_PREFIX) not in EXECUTION_STAGES
                    or status != "failed" or self.status != "failed"
                    or self.execution_sha256 is not None
                )
                for name, status in self.stage_statuses
            )
            or len(self.eligible_member_groups) > 10
            or self.recovery_binding is not None and type(self.recovery_binding) is not ValidationRecoveryBinding
            or (self.status in {"succeeded", "truncated"}) != (self.recovery_binding is not None)
            or len(self.eligible_epcr) != len(self.eligible_member_groups)
            or any(type(value) is not str for value in self.eligible_epcr)
            or any(value <= 0 for value in epcr)
            or any(
                left <= right
                for left, right in zip(epcr, epcr[1:])
            )
            or any(
                values != tuple(sorted(values))
                or any(not name or type(count) is not int or count < 0 for name, count in values)
                or len({name for name, _count in values}) != len(values)
                for values in (
                    self.association_state_counts,
                    self.admission_outcome_counts,
                    self.model_coverage_counts,
                    self.population_lookup_counts,
                )
            )
            or self.model_coverage_counts and tuple(name for name, _ in self.model_coverage_counts) != MODEL_COVERAGE_KEYS
            or self.population_lookup_counts and tuple(name for name, _ in self.population_lookup_counts) != POPULATION_LOOKUP_KEYS
            or not set(self.clinical_only_eligible_member_groups) <= set(self.eligible_member_groups)
            or any(group != tuple(sorted(set(group))) or not group for values in groups for group in values)
            or len(set(self.eligible_member_groups)) != len(self.eligible_member_groups)
            or any(
                value is not None and (type(value) is not int or value < 0)
                for value in (
                    self.candidate_burden,
                    self.hypothesis_count,
                    self.required_enumerations,
                    self.completed_enumerations,
                )
            )
            or self.required_enumerations is not None
            and self.completed_enumerations is not None
            and self.completed_enumerations > self.required_enumerations
            or self.status in {"succeeded", "truncated"} and self.execution_sha256 is None
            or self.status in {"failed", "timed_out", "unsupported"} and self.execution_sha256 is not None
            or self.status in {"succeeded", "truncated"} and any(
                value is None
                for value in (
                    self.candidate_burden,
                    self.hypothesis_count,
                    self.required_enumerations,
                    self.completed_enumerations,
                )
            )
            or self.status in {"failed", "timed_out", "unsupported"} and any(
                value is not None
                for value in (
                    self.candidate_burden,
                    self.hypothesis_count,
                    self.required_enumerations,
                    self.completed_enumerations,
                )
            )
        ):
            raise Track1ValidationRunnerError("invalid validation attempt")


@dataclass(frozen=True, slots=True)
class CaseEvaluation:
    case_id: str
    split: Literal["development", "held_out"]
    scope: str
    policy_id: str
    status: AttemptStatus
    target_kind: str
    target_member_count: int
    admitted_target_members: int
    generated_target_members: int
    exact_target_generated: bool
    exact_target_rank: int | None
    exact_target_eligible_rank: int | None
    required_enumerations: int | None
    completed_enumerations: int | None
    candidate_burden: int | None
    route_admitted_target_members: tuple[tuple[str, int], ...]
    association_state_counts: tuple[tuple[str, int], ...]
    admission_outcome_counts: tuple[tuple[str, int], ...]
    challenge_rank_points: str
    challenge_f_max: str
    model_coverage_counts: tuple[tuple[str, int], ...] = ()
    required_model_coverage: tuple[str, ...] = ()
    clinical_only_target_generated: bool = False
    clinical_only_target_rank: int | None = None
    clinical_only_target_eligible_rank: int | None = None
    target_population_excluded: bool = False
    source_inventory_status: Literal["complete", "not_assessed"] = "not_assessed"
    broader_unresolved_source_count: int | None = None
    broader_unresolved_relevance_source_count: int | None = None

    def __post_init__(self) -> None:
        if self.source_inventory_status not in {"complete", "not_assessed"} or any(
            value is not None and (type(value) is not int or value < 0)
            for value in (self.broader_unresolved_source_count, self.broader_unresolved_relevance_source_count)
        ):
            raise Track1ValidationRunnerError("invalid evaluation source coverage")
        for value in (self.case_id, self.scope, self.policy_id, self.target_kind):
            _identifier(value)
        try:
            rank_points = Fraction(self.challenge_rank_points)
            f_max = Fraction(self.challenge_f_max)
        except (ValueError, ZeroDivisionError, TypeError):
            raise Track1ValidationRunnerError(
                "invalid case evaluation metric"
            ) from None
        counts = (
            self.target_member_count,
            self.admitted_target_members,
            self.generated_target_members,
        )
        optional_counts = (
            self.required_enumerations,
            self.completed_enumerations,
            self.candidate_burden,
        )
        if (
            self.split not in {"development", "held_out"}
            or self.status not in {
                "succeeded", "failed", "timed_out", "unsupported", "truncated",
            }
            or self.target_kind not in {"singleton", "full_pair", "locus"}
            or any(type(value) is not int or value < 0 for value in counts)
            or self.target_member_count < 1
            or self.admitted_target_members > self.target_member_count
            or self.generated_target_members > self.target_member_count
            or type(self.exact_target_generated) is not bool
            or self.exact_target_generated
            and self.generated_target_members != self.target_member_count
            or any(
                value is not None and (type(value) is not int or value < 0)
                for value in optional_counts
            )
            or self.required_enumerations is not None
            and self.completed_enumerations is not None
            and self.completed_enumerations > self.required_enumerations
            or any(
                value is not None and (type(value) is not int or value < 1)
                for value in (self.exact_target_rank, self.exact_target_eligible_rank)
            )
            or rank_points < 0
            or rank_points > 100
            or f_max < 0
            or f_max > 1
            or any(
                values != tuple(sorted(values))
                or len({name for name, _count in values}) != len(values)
                or any(
                    type(name) is not str
                    or not name
                    or type(count) is not int
                    or count < 0
                    for name, count in values
                )
                for values in (
                    self.route_admitted_target_members,
                    self.association_state_counts,
                    self.admission_outcome_counts,
                    self.model_coverage_counts,
                )
            )
            or self.required_model_coverage != tuple(sorted(set(self.required_model_coverage)))
            or not set(self.required_model_coverage) <= set(MODEL_COVERAGE_KEYS)
            or self.model_coverage_counts and tuple(name for name, _ in self.model_coverage_counts) != MODEL_COVERAGE_KEYS
            or any(type(value) is not bool for value in (self.clinical_only_target_generated, self.target_population_excluded))
            or any(value is not None and (type(value) is not int or value < 1) for value in (self.clinical_only_target_rank, self.clinical_only_target_eligible_rank))
            or any(
                count > self.target_member_count
                for _route, count in self.route_admitted_target_members
            )
        ):
            raise Track1ValidationRunnerError("invalid case evaluation")


@dataclass(frozen=True, slots=True)
class DevelopmentExperiment:
    plan_sha256: str
    roster_sha256: str
    protocol: ExperimentProtocol
    candidates: tuple[PolicyCandidate, ...]
    attempts: tuple[ValidationAttempt, ...]
    evaluations: tuple[CaseEvaluation, ...]
    selected_policy_id: str | None
    selection_status: Literal["selected", "no-acceptable-candidate"]
    schema: Literal["genosapiens.track1-public-validation-experiment"] = SCHEMA
    classification: Literal["public"] = "public"
    validation_only: Literal[True] = True
    controlled_use_authorized: Literal[False] = False

    def __post_init__(self) -> None:
        for value in (self.plan_sha256, self.roster_sha256):
            _sha(value)
        if (
            self.schema != SCHEMA
            or self.classification != "public"
            or self.validation_only is not True
            or self.controlled_use_authorized is not False
            or self.selection_status not in {
                "selected", "no-acceptable-candidate",
            }
            or type(self.protocol) is not ExperimentProtocol
            or not self.candidates
            or any(type(item) is not PolicyCandidate for item in self.candidates)
            or tuple(item.policy_id for item in self.candidates)
            != tuple(sorted({item.policy_id for item in self.candidates}))
            or any(type(item) is not ValidationAttempt for item in self.attempts)
            or any(type(item) is not CaseEvaluation for item in self.evaluations)
            or (self.selected_policy_id is None) != (self.selection_status == "no-acceptable-candidate")
            or self.selected_policy_id is not None
            and self.selected_policy_id not in {item.policy_id for item in self.candidates}
        ):
            raise Track1ValidationRunnerError("invalid development experiment")

    def as_bytes(self) -> bytes:
        return _bytes(asdict(self))


@dataclass(frozen=True, slots=True)
class PreHeldOutCandidate:
    candidate_id: str
    plan_sha256: str
    roster_sha256: str
    development_experiment_sha256: str
    development_derivation_sha256: str
    selected_policy_id: str
    policy_bindings: PolicyBindings
    numerical_policy_evidence_sha256: str
    evaluator_sha256: str
    executable_stack_sha256: str
    runtime_image_digest: str
    resource_manifest_sha256s: tuple[str, ...]
    acceptance_criteria: tuple[AcceptanceCriterion, ...]
    schema: Literal["genosapiens.track1-pre-held-out-candidate"] = CANDIDATE_SCHEMA
    classification: Literal["public"] = "public"
    validation_only: Literal[True] = True
    controlled_use_authorized: Literal[False] = False
    held_out_independence: Literal["prespecified"] = "prespecified"

    def __post_init__(self) -> None:
        _identifier(self.candidate_id)
        _identifier(self.selected_policy_id)
        for value in (
            self.plan_sha256,
            self.roster_sha256,
            self.development_experiment_sha256,
            self.development_derivation_sha256,
            self.numerical_policy_evidence_sha256,
            self.evaluator_sha256,
            self.executable_stack_sha256,
            *self.resource_manifest_sha256s,
        ):
            _sha(value)
        if (
            self.schema != CANDIDATE_SCHEMA
            or self.classification != "public"
            or self.validation_only is not True
            or self.controlled_use_authorized is not False
            or self.held_out_independence != "prespecified"
            or type(self.policy_bindings) is not PolicyBindings
            or not re.fullmatch(r"sha256:[0-9a-f]{64}", self.runtime_image_digest)
            or not self.resource_manifest_sha256s
            or self.resource_manifest_sha256s
            != tuple(sorted(set(self.resource_manifest_sha256s)))
            or not self.acceptance_criteria
            or any(
                type(item) is not AcceptanceCriterion
                or item.comparator not in {">=", "<="}
                or not re.fullmatch(r"(?:0|[1-9][0-9]*)(?:\.[0-9]+)?", item.threshold)
                for item in self.acceptance_criteria
            )
            or tuple(item.metric_id for item in self.acceptance_criteria)
            != tuple(sorted({item.metric_id for item in self.acceptance_criteria}))
            or any(
                item.derivation_sha256 != self.development_derivation_sha256
                for item in self.acceptance_criteria
            )
        ):
            raise Track1ValidationRunnerError("invalid pre-held-out candidate")

    def as_bytes(self) -> bytes:
        return _bytes(asdict(self))

    @property
    def sha256(self) -> str:
        return _hash(self.as_bytes())


@dataclass(frozen=True, slots=True)
class AcceptanceResultReference:
    result_id: str
    status: Literal["accepted", "rejected"]
    candidate_sha256: str
    dependency_set_sha256: str
    assessment_scope_sha256: str
    artifact_sha256: str

    def __post_init__(self) -> None:
        _identifier(self.result_id)
        if self.status not in {"accepted", "rejected"}:
            raise Track1ValidationRunnerError("invalid acceptance result status")
        for value in (
            self.candidate_sha256,
            self.dependency_set_sha256,
            self.assessment_scope_sha256,
            self.artifact_sha256,
        ):
            _sha(value)


@dataclass(frozen=True, slots=True)
class ConsolidatedAcceptance:
    candidate_sha256: str
    held_out_attempts_sha256: str
    held_out_metric_values: tuple[tuple[str, str], ...]
    criteria_status: Literal["accepted", "rejected"]
    result_references: tuple[AcceptanceResultReference, ...]
    status: Literal["accepted", "rejected"]
    limitations: tuple[str, ...]
    schema: Literal["genosapiens.track1-consolidated-public-acceptance"] = ACCEPTANCE_SCHEMA
    classification: Literal["public"] = "public"
    controlled_use_authorized: Literal[False] = False

    def __post_init__(self) -> None:
        _sha(self.candidate_sha256)
        _sha(self.held_out_attempts_sha256)
        try:
            metric_values = tuple(
                (name, Fraction(value))
                for name, value in self.held_out_metric_values
            )
        except (ValueError, ZeroDivisionError, TypeError):
            raise Track1ValidationRunnerError(
                "invalid consolidated acceptance metric"
            ) from None
        if (
            self.schema != ACCEPTANCE_SCHEMA
            or self.classification != "public"
            or self.controlled_use_authorized is not False
            or self.criteria_status not in {"accepted", "rejected"}
            or self.status not in {"accepted", "rejected"}
            or any(
                type(item) is not AcceptanceResultReference
                for item in self.result_references
            )
            or tuple(item.result_id for item in self.result_references)
            != REQUIRED_ACCEPTANCE_RESULTS
            or any(
                item.candidate_sha256 != self.candidate_sha256
                for item in self.result_references
            )
            or self.held_out_metric_values
            != tuple(sorted(dict(self.held_out_metric_values).items()))
            or any(
                type(name) is not str
                or not name
                or type(value) is not str
                or fraction < 0
                for (name, value), (_metric_name, fraction) in zip(
                    self.held_out_metric_values, metric_values
                )
            )
            or not self.limitations
            or any(type(item) is not str or not item for item in self.limitations)
            or (self.status == "accepted")
            != (
                self.criteria_status == "accepted"
                and all(item.status == "accepted" for item in self.result_references)
            )
        ):
            raise Track1ValidationRunnerError("invalid consolidated acceptance")

    def as_bytes(self) -> bytes:
        return _bytes(asdict(self))

    @property
    def sha256(self) -> str:
        return _hash(self.as_bytes())


def _member_keys(allele) -> tuple[str, ...]:
    """Recover evaluator-side public source identities without decoding IDs."""
    keys = set()
    for call in allele.calls or ():
        fields = call.source.source_fields or {}
        columns = fields.get("columns")
        number = fields.get("allele_number")
        if (
            isinstance(columns, (tuple, list))
            and len(columns) >= 5
            and type(number) is int
            and number >= 1
        ):
            alternates = str(columns[4]).split(",")
            if number <= len(alternates):
                keys.add(
                    ":".join(
                        (
                            str(columns[0]),
                            str(columns[1]),
                            str(columns[3]),
                            alternates[number - 1],
                        )
                    )
                )
        interpretation = fields.get("native_interpretation") or {}
        if isinstance(interpretation, Mapping):
            locus = interpretation.get("locus_id")
            if type(locus) is str and locus:
                keys.add(locus)
            parent = interpretation.get("parent_call_id")
            if type(parent) is str and parent.startswith("longtr-region:"):
                keys.add(parent)
    return tuple(sorted(keys))


def record_execution_attempt(
    execution: BoundedPublicExecution, *, roster_case: RosterCase,
    policy: PolicyCandidate, attempt_number: int, measurements: ResourceMeasurements,
    population_selector: tuple[str, str, str | None] | None = None,
) -> ValidationAttempt:
    """Reconstruct finite acceptance metrics from complete authenticated shards."""
    from genosapiens.wgs_validation_metrics import bounded_validation_metrics
    if (type(execution) is not BoundedPublicExecution
            or execution.manifest.case_id != roster_case.case_id
            or execution.manifest.candidate_stack_sha256 != policy.bindings.sha256):
        raise Track1ValidationRunnerError("execution does not match roster case")
    manifest = execution.manifest
    metrics = bounded_validation_metrics(execution.products, execution.audit,
        audit_directory=execution.audit_directory,
        scratch_path=execution.audit_directory / "acceptance-metrics.sqlite",
        epcr_assignment=execution.epcr_assignment, population_selector=population_selector)
    return ValidationAttempt(case_id=roster_case.case_id, split=roster_case.split,
        scope=roster_case.scope, policy_id=policy.policy_id, policy_bindings_sha256=policy.bindings.sha256,
        attempt_number=attempt_number, execution_sha256=manifest.sha256,
        input_sha256s=manifest.runner_input_sha256s,
        resource_descriptor_sha256s=manifest.runner_resource_descriptor_sha256s,
        stage_statuses=tuple((key, "succeeded") for key, _ in manifest.output_bindings),
        measurements=measurements, **metrics)


def _population_frequency_observations(alleles, selector):
    """Use the declared rule's exact source stratum, never a hardcoded subset."""
    return tuple(sorted({
        (_hash(member.encode("ascii")), str(Fraction(row.allele_count, row.allele_number)), row.allele_number)
        for allele in alleles for member in _member_keys(allele)
        for row in allele.population or ()
        if (row.subset, row.stratum, row.ancestry) == selector
        and row.allele_count is not None and row.allele_number is not None and row.allele_number > 0
        and 0 <= row.allele_count <= row.allele_number
    }))


def _population_lookup_counts(annotation):
    """Expose only fixed operational counters from the authenticated annotation."""
    if annotation is None:
        return ()
    fields = annotation.report.source_fields
    counts = dict.fromkeys(POPULATION_LOOKUP_KEYS, 0)
    for name in ("requested_alleles", "supported_alleles", "unique_loci",
                 "bcftools_invocations", "tabix_invocations"):
        counts[name] = fields["query_summary"][name]
    for item in fields["assessments"]:
        name = item["reason"] if item["status"] == "unsupported" else item["status"]
        if name not in POPULATION_LOOKUP_KEYS:
            raise Track1ValidationRunnerError("invalid population lookup diagnostic")
        counts[name] += 1
    return tuple(sorted(counts.items()))


def _model_coverage(generation, group, *, excluded_allele_ids=()):
    """Audit replayed model states; rank recovery is evaluated independently."""
    coverage = {name: 0 for name in MODEL_COVERAGE_KEYS}
    coverage["population-annotation-completed"] = int(generation.population_annotation is not None)
    coverage["clinical-annotation-completed"] = int(generation.clinical is not None)
    excluded_ids = set(excluded_allele_ids)
    excluded_groups = set()
    for assessment in generation.population_assessments:
        if assessment.outcome.value == "population_excluded":
            excluded_groups.add(group(assessment.hypothesis))
        for member in assessment.members:
            if member.outcome.value == "population_excluded":
                excluded_ids.add(member.allele_id)
            if member.status == "assessed" and any(row.reliable for row in member.frequencies) and any(rescue.status == "unknown" for rescue in member.rescues):
                coverage["completed-population-lookup-unknown-rescue"] += 1
    clinical_ids = set()
    for item in generation.hypotheses:
        hypothesis = item.hypothesis
        if item.members and all(member.state == "qualifying" and member.clinical_support_refs for member in item.members):
            clinical_ids.add(hypothesis.hypothesis_id)
        if len(hypothesis.evidence_subject_ids) == 2 and item.model.state != "contradicted" and any(member.state == "unresolved" for member in item.members) and not any(member.clinical_support_refs for member in item.members):
            coverage["clinical-unasserted-exploratory-pair"] += 1
        if hypothesis.completeness.value == "unknown" and item.model.state != "contradicted" and all(member.genotype.alternate_present is True for member in item.members):
            coverage["unknown-count-hypothesis"] += 1
    for search in generation.searches:
        unresolved = {member.allele_id for member in search.member_assessments if member.state == "unresolved"}
        coverage["population-excluded-unresolved-partner"] += len(unresolved & excluded_ids)
        for members, model in zip(search.pair_assessments, search.pair_models, strict=True):
            if model.state == "contradicted" and any(member.allele_id in excluded_ids for member in members) and any(member.state == "nonqualifying" for member in members):
                coverage["independent-non-population-contradiction"] += 1
        for item in generation.hypotheses:
            if item.admission_context_id != search.context_id or not item.conditional_singleton_search_satisfied:
                continue
            own = set(item.hypothesis.evidence_subject_ids)
            if search.enumeration_status is not EnumerationStatus.COMPLETE or any(any(member.allele_id in own for member in pair) and any(member.allele_id not in own and member.state != "nonqualifying" for member in pair) for pair in search.pair_assessments):
                coverage["unresolved-partner-singleton-violation"] += 1
    return coverage, excluded_groups, clinical_ids


def _terminal_attempts(
    attempts: tuple[ValidationAttempt, ...],
    *,
    split: str,
    policy_id: str,
) -> dict[str, ValidationAttempt]:
    selected = [
        item for item in attempts
        if item.split == split and item.policy_id == policy_id
    ]
    grouped: dict[str, list[ValidationAttempt]] = {}
    for item in selected:
        grouped.setdefault(item.case_id, []).append(item)
    terminal = {}
    for case_id, values in grouped.items():
        values.sort(key=lambda item: item.attempt_number)
        if tuple(item.attempt_number for item in values) != tuple(range(1, len(values) + 1)):
            raise Track1ValidationRunnerError("validation retry accounting is incomplete")
        terminal[case_id] = values[-1]
    return terminal


def evaluate_attempts(
    roster: Track1Roster,
    attempts: tuple[ValidationAttempt, ...],
    *,
    split: Literal["development", "held_out"],
    policy_id: str,
    recovery_resolver=None,
) -> tuple[CaseEvaluation, ...]:
    """Join isolated evaluator truth only after terminal runner attempts exist."""
    cases = tuple(item for item in roster.cases if item.split == split)
    terminal = _terminal_attempts(attempts, split=split, policy_id=policy_id)
    if set(terminal) != {item.case_id for item in cases}:
        raise Track1ValidationRunnerError("validation split attempt accounting is incomplete")
    evaluated = []
    for case in cases:
        attempt = terminal[case.case_id]
        target = tuple(sorted(case.truth_members))
        if attempt.recovery_binding is None:
            recovery = TargetRecovery(0, 0, False, None, False, None, False, ())
        else:
            if recovery_resolver is None:
                raise Track1ValidationRunnerError("authenticated normalized recovery resolver required")
            recovery = recovery_resolver(attempt, target)
            if type(recovery) is not TargetRecovery:
                raise Track1ValidationRunnerError("invalid normalized target recovery")
        generated_rank = recovery.generated_rank
        eligible_rank = next(
            (index for index, group in enumerate(attempt.eligible_member_groups, 1) if group == target), None)
        route_recall = recovery.route_recall
        full_rank = next(
            (
                index
                for index, group in enumerate(attempt.eligible_member_groups, 1)
                if group == target
            ),
            None,
        )
        partial_rank = None
        if full_rank is None and len(target) == 2:
            partial_rank = next(
                (
                    index
                    for index, group in enumerate(attempt.eligible_member_groups, 1)
                    if set(group) & set(target)
                ),
                None,
            )

        def tier(rank: int | None) -> Fraction:
            if rank is None or rank > 10:
                return Fraction(0)
            value = 100 if rank == 1 else 50 if rank <= 3 else 25 if rank <= 5 else 10
            return Fraction(value)

        points = tier(full_rank)
        if full_rank is None and partial_rank is not None:
            points = tier(partial_rank) / 2
        predictions: set[str] = set()
        best_f = Fraction(0)
        for group in attempt.eligible_member_groups:
            predictions.update(group)
            true_positive = len(predictions & set(target))
            if true_positive:
                best_f = max(
                    best_f,
                    Fraction(2 * true_positive, len(predictions) + len(target)),
                )
        evaluated.append(CaseEvaluation(
            case.case_id,
            split,
            case.scope,
            policy_id,
            attempt.status,
            case.target_kind,
            len(target),
            recovery.admitted_member_count,
            recovery.generated_member_count,
            recovery.exact_generated,
            generated_rank,
            eligible_rank,
            attempt.required_enumerations,
            attempt.completed_enumerations,
            attempt.candidate_burden,
            route_recall,
            attempt.association_state_counts,
            attempt.admission_outcome_counts,
            _fraction_text(points),
            _fraction_text(best_f),
            attempt.model_coverage_counts,
            case.required_model_coverage,
            recovery.clinical_generated,
            recovery.clinical_rank,
            next((index for index, value in enumerate(attempt.clinical_only_eligible_member_groups, 1) if value == target), None),
            recovery.population_excluded,
            attempt.source_inventory_status,
            attempt.broader_unresolved_source_count,
            attempt.broader_unresolved_relevance_source_count,
        ))
    return tuple(evaluated)


def complete_development_experiment(
    *,
    plan: ValidationPlan,
    roster: Track1Roster,
    protocol: ExperimentProtocol,
    candidates: tuple[PolicyCandidate, ...],
    attempts: tuple[ValidationAttempt, ...],
    recovery_resolver=None,
) -> DevelopmentExperiment:
    """Complete the declared family and deterministically select one candidate."""
    if (
        type(plan) is not ValidationPlan
        or type(roster) is not Track1Roster
        or type(protocol) is not ExperimentProtocol
        or not candidates
        or len(candidates) > protocol.maximum_policy_candidates
        or tuple(item.policy_id for item in candidates)
        != tuple(sorted({item.policy_id for item in candidates}))
        or any(type(item) is not ValidationAttempt for item in attempts)
    ):
        raise Track1ValidationRunnerError("invalid development experiment inputs")
    evaluations = tuple(
        item
        for candidate in candidates
        for item in evaluate_attempts(
            roster, attempts, split="development", policy_id=candidate.policy_id, recovery_resolver=recovery_resolver
        )
    )
    for candidate in candidates:
        per_case = [item for item in attempts if item.policy_id == candidate.policy_id]
        if any(
            len(values) > protocol.maximum_attempts_per_case
            for case_id in {item.case_id for item in per_case}
            if (values := [item for item in per_case if item.case_id == case_id])
        ):
            raise Track1ValidationRunnerError("development attempt budget exceeded")
    by_policy = {
        candidate.policy_id: tuple(
            item for item in evaluations if item.policy_id == candidate.policy_id
        )
        for candidate in candidates
    }

    def objective(candidate: PolicyCandidate):
        values = by_policy[candidate.policy_id]
        top10 = sum(item.exact_target_rank is not None and item.exact_target_rank <= 10 for item in values)
        top1 = sum(item.exact_target_rank == 1 for item in values)
        exact = sum(item.exact_target_generated for item in values)
        burden = sum(item.candidate_burden or 0 for item in values)
        return (-top10, -top1, -exact, burden, candidate.policy_id)

    successful = [
        candidate for candidate in candidates
        if all(item.status == "succeeded" for item in by_policy[candidate.policy_id])
        and model_coverage_complete(by_policy[candidate.policy_id])
    ]
    selected = min(successful, key=objective) if successful else None
    return DevelopmentExperiment(
        _hash(validation_plan_as_json(plan)),
        roster.roster_sha256,
        protocol,
        candidates,
        attempts,
        evaluations,
        None if selected is None else selected.policy_id,
        "no-acceptable-candidate" if selected is None else "selected",
    )


def freeze_preheldout_candidate(
    experiment: DevelopmentExperiment,
    *,
    derivation_bytes: bytes,
    numerical_policy_evidence_sha256: str,
    evaluator_sha256: str,
    executable_stack_sha256: str,
    runtime_image_digest: str,
    resource_manifest_sha256s: tuple[str, ...],
    acceptance_criteria: tuple[AcceptanceCriterion, ...],
) -> PreHeldOutCandidate:
    """Authenticate the completed development record before held-out access."""
    if (
        type(experiment) is not DevelopmentExperiment
        or experiment.selection_status != "selected"
        or experiment.selected_policy_id is None
        or type(derivation_bytes) is not bytes
        or not derivation_bytes
        or not model_coverage_complete(tuple(
            value for value in experiment.evaluations if value.policy_id == experiment.selected_policy_id
        ))
    ):
        raise Track1ValidationRunnerError("development produced no candidate")
    selected = next(
        item for item in experiment.candidates
        if item.policy_id == experiment.selected_policy_id
    )
    derivation_sha256 = _hash(derivation_bytes)
    candidate_id = "track1-base-candidate-" + _hash(_bytes({
        "experiment": _hash(experiment.as_bytes()),
        "derivation": derivation_sha256,
        "policy": selected.bindings.sha256,
    }))[:16]
    return PreHeldOutCandidate(
        candidate_id,
        experiment.plan_sha256,
        experiment.roster_sha256,
        _hash(experiment.as_bytes()),
        derivation_sha256,
        selected.policy_id,
        selected.bindings,
        numerical_policy_evidence_sha256,
        evaluator_sha256,
        executable_stack_sha256,
        runtime_image_digest,
        tuple(sorted(resource_manifest_sha256s)),
        acceptance_criteria,
    )


def development_ranker_view(
    roster: Track1Roster,
    case_id: str,
    *,
    runner_resources: tuple[RankerResource, ...],
) -> RankerCaseView:
    """Build the truth-free development view from the authenticated roster."""
    if (
        type(roster) is not Track1Roster
        or not runner_resources
        or any(type(item) is not RankerResource for item in runner_resources)
    ):
        raise Track1ValidationRunnerError("development runner inputs are invalid")
    cases = tuple(item for item in roster.cases if item.case_id == case_id)
    if len(cases) != 1 or cases[0].split != "development":
        raise Track1ValidationRunnerError("requested case is not development")
    view = roster.runner_case(case_id)
    return _ranker_view(view, runner_resources)


def _ranker_view(view, runner_resources: tuple[RankerResource, ...]) -> RankerCaseView:
    """Project the common roster bindings without case purpose or truth."""
    inputs = tuple(
        RankerInput(
            InputKind.HPO_TEXT if item.role == "hpo_text" else InputKind.VCF,
            item.sha256,
        )
        for item in view.inputs
    )
    resources = (*runner_resources, *(
        RankerResource("runner-source-format", value)
        for value in view.source_format_sha256s
    ), RankerResource("runner-public-background", view.background_sha256))
    rendered = RankerCaseView(inputs, tuple(resources))
    lowered = _bytes(asdict(rendered)).lower()
    if any(value in lowered for value in (b"truth", b"causal", b"oracle", b"recipe", b"seed")):
        raise Track1ValidationRunnerError("runner view crosses truth boundary")
    return rendered


def authenticated_held_out_ranker_view(
    candidate: PreHeldOutCandidate,
    experiment: DevelopmentExperiment,
    roster: Track1Roster,
    case_id: str,
    *,
    runner_resources: tuple[RankerResource, ...],
) -> RankerCaseView:
    """Release one truth-free held-out view only after exact candidate replay."""
    if (
        type(candidate) is not PreHeldOutCandidate
        or type(experiment) is not DevelopmentExperiment
        or type(roster) is not Track1Roster
        or candidate.development_experiment_sha256 != _hash(experiment.as_bytes())
        or candidate.plan_sha256 != experiment.plan_sha256
        or candidate.roster_sha256 != roster.roster_sha256
        or candidate.selected_policy_id != experiment.selected_policy_id
        or candidate.policy_bindings != next(
            item.bindings for item in experiment.candidates
            if item.policy_id == experiment.selected_policy_id
        )
        or not runner_resources
        or any(type(item) is not RankerResource for item in runner_resources)
    ):
        raise Track1ValidationRunnerError("pre-held-out candidate authentication failed")
    cases = tuple(item for item in roster.cases if item.case_id == case_id)
    if len(cases) != 1 or cases[0].split != "held_out":
        raise Track1ValidationRunnerError("requested case is not held out")
    view = roster.runner_case(case_id, allow_held_out=True)
    return _ranker_view(view, runner_resources)


def evaluation_metrics(values: tuple[CaseEvaluation, ...]) -> dict[str, Fraction]:
    """Aggregate exact evaluator-side numerators over a complete split."""
    if not values or any(type(item) is not CaseEvaluation for item in values):
        raise Track1ValidationRunnerError("evaluation metric inputs are incomplete")

    def fraction(numerator: int, denominator: int) -> Fraction:
        return Fraction(numerator, denominator) if denominator else Fraction(0)

    succeeded = tuple(item for item in values if item.status == "succeeded")
    singleton = tuple(item for item in values if item.target_kind == "singleton")
    pairs = tuple(item for item in values if item.target_kind == "full_pair")
    targets = sum(item.target_member_count for item in values)
    required = sum(item.required_enumerations or 0 for item in values)
    completed = sum(item.completed_enumerations or 0 for item in values)
    burden_values = tuple(
        item.candidate_burden
        for item in values
        if item.candidate_burden is not None
    )
    route_names = tuple(sorted({
        route
        for item in values
        for route, _count in item.route_admitted_target_members
    }))
    association_names = tuple(sorted({
        name for item in values for name, _count in item.association_state_counts
    }))
    outcome_names = tuple(sorted({
        name for item in values for name, _count in item.admission_outcome_counts
    }))
    metrics = {
        "admission_causal_allele_recall": fraction(
            sum(item.admitted_target_members for item in values), targets,
        ),
        "candidate_generation_recall": fraction(
            sum(item.generated_target_members for item in values), targets,
        ),
        "enumeration_coverage": fraction(completed, required),
        "exact_full_pair_recovery": fraction(
            sum(item.exact_target_generated for item in pairs), len(pairs),
        ),
        "exact_singleton_recovery": fraction(
            sum(item.exact_target_generated for item in singleton), len(singleton),
        ),
        "full_pair_member_availability": fraction(
            sum(item.generated_target_members for item in pairs),
            sum(item.target_member_count for item in pairs),
        ),
        "top_10_full_genotype_recovery": fraction(
            sum(item.exact_target_rank is not None and item.exact_target_rank <= 10 for item in values),
            len(values),
        ),
        "top_1_full_genotype_recovery": fraction(
            sum(item.exact_target_rank == 1 for item in values), len(values),
        ),
        "successful_case_fraction": fraction(len(succeeded), len(values)),
        "context_recovery": fraction(
            sum(item.exact_target_generated for item in values), len(values),
        ),
        "mean_candidate_burden": fraction(sum(burden_values), len(burden_values)),
        "maximum_candidate_burden": Fraction(max(burden_values, default=0)),
        "challenge_rank_points": fraction(
            sum(Fraction(item.challenge_rank_points) for item in values), len(values),
        ),
        "challenge_f_max": fraction(
            sum(Fraction(item.challenge_f_max) for item in values), len(values),
        ),
    }
    for route in route_names:
        metrics[f"route_{route}_causal_allele_recall"] = fraction(
            sum(dict(item.route_admitted_target_members).get(route, 0) for item in values),
            targets,
        )
    for name in association_names:
        metrics[f"association_state_{name}_count"] = Fraction(sum(
            dict(item.association_state_counts).get(name, 0) for item in values
        ))
    for name in outcome_names:
        metrics[f"admission_outcome_{name}_count"] = Fraction(sum(
            dict(item.admission_outcome_counts).get(name, 0) for item in values
        ))
    metrics["required_model_coverage_fraction"] = fraction(sum(all(dict(item.model_coverage_counts).get(name, 0) > 0 for name in item.required_model_coverage) for item in values), len(values))
    metrics["unresolved_partner_singleton_violation_count"] = Fraction(sum(dict(item.model_coverage_counts).get("unresolved-partner-singleton-violation", 0) for item in values))
    metrics["clinical_only_exact_genotype_recovery"] = fraction(sum(item.clinical_only_target_generated for item in values), len(values))
    metrics["clinical_only_top_10_genotype_recovery"] = fraction(sum(item.clinical_only_target_rank is not None and item.clinical_only_target_rank <= 10 for item in values), len(values))
    metrics["clinical_only_eligible_top_10_genotype_recovery"] = fraction(sum(item.clinical_only_target_eligible_rank is not None and item.clinical_only_target_eligible_rank <= 10 for item in values), len(values))
    metrics["generated_target_population_excluded_fraction"] = fraction(sum(item.target_population_excluded for item in values), len(values))
    for name in MODEL_COVERAGE_KEYS:
        metrics["model_coverage_" + name.replace("-", "_") + "_count"] = Fraction(sum(dict(item.model_coverage_counts).get(name, 0) for item in values))
    return metrics


def model_coverage_complete(values: tuple[CaseEvaluation, ...]) -> bool:
    """Required scenario presence is a structural gate, not a recovery floor."""
    declared = {name for value in values for name in value.required_model_coverage}
    return bool(values) and bool(declared) and declared <= set(REQUIRED_MODEL_COVERAGE) and all(
        item.status == "succeeded"
        and dict(item.model_coverage_counts).get("clinical-annotation-completed", 0) == 1
        and dict(item.model_coverage_counts).get("population-annotation-completed", 0) == 1
        and all(dict(item.model_coverage_counts).get(name, 0) > 0 for name in item.required_model_coverage)
        and dict(item.model_coverage_counts).get("unresolved-partner-singleton-violation", 0) == 0
        for item in values
    )


def consolidate_public_acceptance(
    candidate: PreHeldOutCandidate,
    roster: Track1Roster,
    held_out_attempts: tuple[ValidationAttempt, ...],
    result_references: tuple[AcceptanceResultReference, ...],
    *,
    limitations: tuple[str, ...],
    recovery_resolver=None,
) -> ConsolidatedAcceptance:
    """Evaluate fixed held-out criteria and reference matching 12.x evidence."""
    if (
        type(candidate) is not PreHeldOutCandidate
        or type(roster) is not Track1Roster
        or candidate.roster_sha256 != roster.roster_sha256
        or any(type(item) is not ValidationAttempt for item in held_out_attempts)
        or any(type(item) is not AcceptanceResultReference for item in result_references)
        or tuple(item.result_id for item in result_references)
        != REQUIRED_ACCEPTANCE_RESULTS
        or any(item.candidate_sha256 != candidate.sha256 for item in result_references)
    ):
        raise Track1ValidationRunnerError("acceptance evidence does not match candidate")
    evaluations = evaluate_attempts(
        roster,
        held_out_attempts,
        split="held_out",
        policy_id=candidate.selected_policy_id,
        recovery_resolver=recovery_resolver,
    )
    metrics = evaluation_metrics(evaluations)
    criteria_pass = model_coverage_complete(evaluations)
    for criterion in candidate.acceptance_criteria:
        actual = metrics.get(criterion.metric_id)
        if actual is None:
            raise Track1ValidationRunnerError("acceptance criterion has no held-out metric")
        threshold = Fraction(criterion.threshold)
        criteria_pass &= (
            actual >= threshold
            if criterion.comparator == ">="
            else actual <= threshold
        )
    criteria_status = "accepted" if criteria_pass else "rejected"
    status = "accepted" if criteria_pass and all(
        item.status == "accepted" for item in result_references
    ) else "rejected"
    metric_values = tuple(sorted(
        (key, str(value.numerator) if value.denominator == 1 else f"{value.numerator}/{value.denominator}")
        for key, value in metrics.items()
    ))
    attempt_bytes = _bytes(tuple(asdict(item) for item in held_out_attempts))
    return ConsolidatedAcceptance(
        candidate.sha256,
        _hash(attempt_bytes),
        metric_values,
        criteria_status,
        result_references,
        status,
        limitations,
    )
