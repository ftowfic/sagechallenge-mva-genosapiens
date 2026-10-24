"""Structural G7 review for staged, lesion-qualified candidate validation plans."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Mapping

from genosapiens.track2_gates import COVERAGE_STATES, GateResult


REQUIRED_STAGES = {
    "longer_term_clonal",
    "proximal_molecular_rescue",
    "subtype_functional_rescue",
    "target_engagement",
    "toxicity",
}
REQUIRED_READOUTS = {
    "apoptosis",
    "cell_cycle",
    "genetic_toxicity",
    "longer_term_clonal",
    "proliferation",
    "proximal_rescue",
    "subtype_relevant_rescue",
    "target_engagement",
    "viability",
}


class ValidationPlanError(ValueError):
    """Raised when G7 audit metadata or a plan identity is undefined."""


@dataclass(frozen=True, slots=True)
class ModelQualificationStage:
    genotype_matched: bool
    subtype_matched: bool
    discriminating_baseline_lesion: str
    baseline_assay: str
    correction_method: str
    correction_normalization_lesion: str
    go_no_go_before_drug_testing: bool
    status: str
    evidence_ids: tuple[str, ...]


@dataclass(frozen=True, slots=True)
class ValidationControls:
    untreated: bool
    vehicle: bool
    wild_type_or_corrected: bool
    target_engagement_positive: bool
    target_engagement_negative: bool
    orthogonal_mechanism: bool
    harmful_comparator_state: str
    harmful_comparator_rationale: str


@dataclass(frozen=True, slots=True)
class DoseResponsePlan:
    concentrations: tuple[float, ...]
    unit: str
    human_exposure_upper_bound: float
    exposure_evidence_id: str
    schedule: str
    duration: str
    target_occupancy_or_engagement_measure: str


@dataclass(frozen=True, slots=True)
class ValidationPlan:
    plan_id: str
    active_ingredient_id: str
    model_qualification: ModelQualificationStage
    stages: tuple[str, ...]
    controls: ValidationControls
    dose_response: DoseResponsePlan
    readouts: tuple[str, ...]
    primary_endpoint: str
    meaningful_effect_threshold: str
    biological_replicate_rationale: str
    independent_clones_or_donors_state: str
    independent_clones_or_donors_rationale: str
    assay_quality_controls: tuple[str, ...]
    blinded_quantification_state: str
    blinded_quantification_rationale: str
    statistical_analysis: str
    multiplicity_handling: str
    stopping_criteria: tuple[str, ...]
    aneuploidy_rescue_requires_no_arrest_or_selective_cell_loss: bool
    candidate_experiment_status: str
    existing_result_interpretation: str
    evidence_ids: tuple[str, ...]

    def __post_init__(self) -> None:
        if not self.plan_id or not self.active_ingredient_id or not self.evidence_ids:
            raise ValidationPlanError("validation plan identity and evidence are required")


@dataclass(frozen=True, slots=True)
class StructuralG7Review:
    gate_result: GateResult
    plan_structurally_complete: bool
    drug_testing_authorized_in_model: bool
    empirical_rescue_evidence_created: bool
    efficacy_evidence_created: bool
    failed_requirements: tuple[str, ...]


def review_structural_g7(
    plan: ValidationPlan,
    *,
    coverage_state: str,
    review_state: str,
    decision_policy: Mapping[str, object],
) -> StructuralG7Review:
    if coverage_state not in COVERAGE_STATES:
        raise ValidationPlanError("G7 coverage state is undefined")
    review = decision_policy.get("review")
    allowed_reviews = set(review.get("independence_states", ())) if isinstance(review, Mapping) else set()
    if review_state not in allowed_reviews:
        raise ValidationPlanError("G7 review state is undefined")
    model = plan.model_qualification
    model_failures: list[str] = []
    if not model.genotype_matched:
        model_failures.append("genotype_matched_model")
    if not model.subtype_matched:
        model_failures.append("subtype_matched_model")
    if not model.discriminating_baseline_lesion or not model.baseline_assay:
        model_failures.append("discriminating_baseline_lesion")
    if model.correction_method not in {"complementation", "isogenic_correction"}:
        model_failures.append("isogenic_correction_or_complementation")
    if (
        not model.correction_normalization_lesion
        or model.correction_normalization_lesion != model.discriminating_baseline_lesion
    ):
        model_failures.append("correction_normalizes_discriminating_lesion")
    if not model.go_no_go_before_drug_testing:
        model_failures.append("model_qualification_precedes_drug_testing")
    if model.status not in {"failed", "passed", "planned"} or not model.evidence_ids:
        model_failures.append("model_qualification_status")
    if model.status == "failed":
        model_failures.append("existing_model_qualification_failure")

    design_failures: list[str] = []
    if set(plan.stages) != REQUIRED_STAGES:
        design_failures.append("complete_staged_plan")
    controls = plan.controls
    if not all(
        (
            controls.untreated,
            controls.vehicle,
            controls.wild_type_or_corrected,
            controls.target_engagement_positive,
            controls.target_engagement_negative,
            controls.orthogonal_mechanism,
        )
    ):
        design_failures.append("required_controls")
    if controls.harmful_comparator_state not in {
        "included",
        "not_available_with_rationale",
    } or not controls.harmful_comparator_rationale:
        design_failures.append("harmful_or_discriminating_comparator")
    dose = plan.dose_response
    if (
        len(dose.concentrations) < 3
        or any(value < 0 for value in dose.concentrations)
        or tuple(sorted(set(dose.concentrations))) != dose.concentrations
        or not dose.unit
        or dose.human_exposure_upper_bound <= 0
        or max(dose.concentrations) > dose.human_exposure_upper_bound
        or not dose.exposure_evidence_id
        or not dose.schedule
        or not dose.duration
        or not dose.target_occupancy_or_engagement_measure
    ):
        design_failures.append("exposure_bounded_dose_response")
    if set(plan.readouts) != REQUIRED_READOUTS:
        design_failures.append("mechanism_rescue_and_toxicity_readouts")
    if plan.primary_endpoint not in REQUIRED_READOUTS or not plan.meaningful_effect_threshold:
        design_failures.append("primary_endpoint_and_effect_threshold")
    if not plan.biological_replicate_rationale:
        design_failures.append("biological_replicate_rationale")
    if plan.independent_clones_or_donors_state not in {
        "included",
        "not_feasible_with_rationale",
    } or not plan.independent_clones_or_donors_rationale:
        design_failures.append("independent_clones_or_donors")
    if not plan.assay_quality_controls:
        design_failures.append("assay_quality_controls")
    if plan.blinded_quantification_state not in {
        "included",
        "not_feasible_with_rationale",
    } or not plan.blinded_quantification_rationale:
        design_failures.append("blinded_quantification")
    if not plan.statistical_analysis or not plan.multiplicity_handling:
        design_failures.append("statistical_analysis_and_multiplicity")
    if not plan.stopping_criteria:
        design_failures.append("stopping_criteria")
    if not plan.aneuploidy_rescue_requires_no_arrest_or_selective_cell_loss:
        design_failures.append("cytostasis_and_selective_cell_loss_exclusion")
    if plan.candidate_experiment_status not in {"executed", "not_executed"}:
        design_failures.append("candidate_experiment_status")
    if plan.existing_result_interpretation not in {
        "aneuploid_reduction_with_arrest_or_death",
        "no_existing_result",
        "qualified_unconfounded_result",
        "viability_only",
    }:
        design_failures.append("existing_result_interpretation")
    if plan.existing_result_interpretation in {
        "aneuploid_reduction_with_arrest_or_death",
        "viability_only",
    }:
        design_failures.append("existing_result_is_not_mechanism_rescue")

    if model_failures:
        outcome, reasons = "fail", ("unqualified_model",)
    elif design_failures:
        outcome, reasons = "fail", ("validation_design_insufficient",)
    else:
        outcome, reasons = "pass", ("pass",)
    complete = not model_failures and not design_failures
    testing_authorized = complete and model.status == "passed"
    empirical_created = (
        complete
        and plan.candidate_experiment_status == "executed"
        and plan.existing_result_interpretation == "qualified_unconfounded_result"
    )
    evidence_ids = tuple(sorted(set(plan.evidence_ids) | set(model.evidence_ids)))
    return StructuralG7Review(
        gate_result=GateResult(
            "G7", outcome, reasons, evidence_ids, coverage_state,
            decision_policy["policy_version"], review_state,
        ),
        plan_structurally_complete=complete,
        drug_testing_authorized_in_model=testing_authorized,
        empirical_rescue_evidence_created=empirical_created,
        efficacy_evidence_created=False,
        failed_requirements=tuple(model_failures + design_failures),
    )
