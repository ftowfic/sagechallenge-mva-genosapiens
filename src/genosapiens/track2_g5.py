"""Versioned G5 target-site exposure and effect-magnitude evaluation."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Mapping

from genosapiens.track2_exposure import ActivityExposureComparison
from genosapiens.track2_gates import COVERAGE_STATES, GateResult


class G5EvaluationError(ValueError):
    """Raised when exposure uncertainty or G5 audit metadata is incomplete."""


@dataclass(frozen=True, slots=True)
class ExposureUncertainty:
    name: str
    material: bool
    bounded: bool
    evidence_ids: tuple[str, ...]

    def __post_init__(self) -> None:
        if not self.name or not self.evidence_ids:
            raise G5EvaluationError("exposure uncertainty is incomplete")


def evaluate_g5(
    *,
    comparison: ActivityExposureComparison | None,
    uncertainties: tuple[ExposureUncertainty, ...],
    time_above_threshold_supported: bool | None,
    evidence_ids: tuple[str, ...],
    coverage_state: str,
    review_state: str,
    decision_policy: Mapping[str, object],
) -> GateResult:
    policy = decision_policy.get("g5_exposure_policy")
    review = decision_policy.get("review")
    allowed_reviews = set(review.get("independence_states", ())) if isinstance(review, Mapping) else set()
    if (
        not isinstance(policy, Mapping)
        or not evidence_ids
        or coverage_state not in COVERAGE_STATES
        or review_state not in allowed_reviews
    ):
        raise G5EvaluationError("G5 policy or audit metadata is incomplete")
    if len({item.name for item in uncertainties}) != len(uncertainties):
        raise G5EvaluationError("exposure uncertainties must be unique")
    if comparison is None:
        outcome = "not_evaluable"
        reasons = ("target_site_exposure_unknown", "rescue_magnitude_unknown")
    else:
        evidence_ids = tuple(
            sorted(
                set(evidence_ids)
                | {comparison.activity_evidence_id, comparison.exposure_evidence_id}
                | {evidence for item in uncertainties for evidence in item.evidence_ids}
            )
        )
        ratio = comparison.potency_exposure_ratio
        effect_state = comparison.effect_magnitude_state
        if (
            (ratio is not None and ratio < policy["minimum_potency_exposure_ratio"])
            or time_above_threshold_supported is False
            or effect_state
            in {"below_endpoint_relevant_magnitude", "wrong_direction"}
        ):
            outcome = "fail"
            reasons = ("exposure_infeasible",)
        else:
            unknown_reasons: list[str] = []
            if ratio is None or time_above_threshold_supported is None:
                unknown_reasons.append("target_site_exposure_unknown")
            if effect_state == "magnitude_unknown":
                unknown_reasons.append("rescue_magnitude_unknown")
            allowed = set(policy["bounded_secondary_uncertainties"])
            unbounded = tuple(
                item
                for item in uncertainties
                if item.material or not item.bounded or item.name not in allowed
            )
            if unbounded:
                unknown_reasons.append("target_site_exposure_unknown")
            if unknown_reasons:
                outcome = "not_evaluable"
                reasons = tuple(dict.fromkeys(unknown_reasons))
            elif uncertainties:
                outcome = "bounded_uncertainty"
                reasons = ("bounded_secondary_uncertainty",)
            else:
                outcome = "pass"
                reasons = ("pass",)
    allowed_reasons = set(decision_policy["gates"]["G5"]["reason_codes"])
    if set(reasons) - allowed_reasons:
        raise G5EvaluationError("G5 emitted an undefined reason code")
    return GateResult(
        gate_id="G5",
        outcome=outcome,
        reason_codes=reasons,
        evidence_ids=evidence_ids,
        coverage_state=coverage_state,
        policy_version=decision_policy["policy_version"],
        review_state=review_state,
    )
