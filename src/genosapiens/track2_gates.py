"""Versioned, fail-closed G0-G4 candidate truth table."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Mapping


GATE_IDS = tuple(f"G{number}" for number in range(5))
REQUIRED_INPUTS = {
    "input_integrity",
    "ingredient_identity",
    "current_market_approval",
    "named_jurisdiction",
    "exact_subtype_applicability",
    "mechanism_readiness",
    "lesion_qualification",
    "credible_mechanism_alternatives",
    "drug_action_direction",
    "direct_contradictions",
    "independent_pharmacology",
    "independent_rescue_lineage",
    "translational_maturity",
}
ALLOWED_VALUES = {
    "input_integrity": {"invalid", "valid"},
    "ingredient_identity": {"ambiguous", "non_medication", "resolved", "unresolved"},
    "current_market_approval": {
        "current_us",
        "not_current",
        "reviewed_named_jurisdiction",
        "unknown",
    },
    "exact_subtype_applicability": {
        "exact_subtype",
        "other_subtype",
        "shared_mechanism_supported",
        "umbrella_only",
        "unresolved",
    },
    "mechanism_readiness": {"bounded_uncertainty", "not_evaluable", "not_ready", "ready", "prediction_accepted"},
    "lesion_qualification": {"qualified", "unknown", "unqualified"},
    "credible_mechanism_alternatives": {"aligned", "bounded", "harmful", "none", "unknown"},
    "drug_action_direction": {"aligned", "binding_only", "context_unresolved", "opposed", "unknown"},
    "direct_contradictions": {"absent", "present", "unknown"},
    "independent_pharmacology": {"copied_only", "independent", "missing", "pathway_proximity"},
    "independent_rescue_lineage": {"copied_only", "independent", "missing", "pathway_proximity"},
}
COVERAGE_STATES = {"complete", "incomplete", "no_result", "success", "truncated", "unavailable"}


class GateEvaluationError(ValueError):
    """Raised when a gate input or audit field is omitted or undefined."""


@dataclass(frozen=True, slots=True)
class GateResult:
    gate_id: str
    outcome: str
    reason_codes: tuple[str, ...]
    evidence_ids: tuple[str, ...]
    coverage_state: str
    policy_version: str
    review_state: str


@dataclass(frozen=True, slots=True)
class G0G4Evaluation:
    results: tuple[GateResult, ...]
    stopped_after_g0: bool


def _decision(outcome: str, *reasons: str) -> tuple[str, tuple[str, ...]]:
    return outcome, tuple(reasons)


def _g0(values: Mapping[str, str]) -> tuple[str, tuple[str, ...]]:
    return (
        _decision("pass", "pass")
        if values["input_integrity"] == "valid"
        else _decision("fail", "integrity_failure")
    )


def _g1(values: Mapping[str, str]) -> tuple[str, tuple[str, ...]]:
    if values["ingredient_identity"] != "resolved":
        return _decision("fail", "identity_ambiguous")
    approval = values["current_market_approval"]
    jurisdiction = values["named_jurisdiction"]
    if approval == "current_us" and jurisdiction == "US":
        return _decision("pass", "pass")
    if approval == "reviewed_named_jurisdiction" and jurisdiction != "US":
        return _decision("pass", "non_us_reviewed_exception", "pass")
    return _decision("fail", "not_currently_marketed")


def _g2(values: Mapping[str, str]) -> tuple[str, tuple[str, ...]]:
    applicability = values["exact_subtype_applicability"]
    if applicability in {"umbrella_only", "other_subtype", "unresolved"}:
        return _decision("not_evaluable", "umbrella_only")
    if values["credible_mechanism_alternatives"] == "harmful":
        return _decision("fail", "direct_contradiction")
    if values["lesion_qualification"] == "unqualified":
        return _decision("fail", "unqualified_model")
    if values["lesion_qualification"] == "unknown":
        return _decision("not_evaluable", "unqualified_model")
    readiness = values["mechanism_readiness"]
    alternatives = values["credible_mechanism_alternatives"]
    if readiness in {"not_evaluable", "not_ready", "prediction_accepted"} or alternatives == "unknown":
        return _decision("not_evaluable", "direction_unknown")
    if readiness == "bounded_uncertainty" or alternatives == "bounded":
        return _decision("bounded_uncertainty", "bounded_secondary_uncertainty")
    return _decision("pass", "pass")


def _g3(values: Mapping[str, str]) -> tuple[str, tuple[str, ...]]:
    if values["drug_action_direction"] == "opposed" or values["direct_contradictions"] == "present":
        return _decision("fail", "direct_contradiction")
    if (
        values["drug_action_direction"] in {"binding_only", "context_unresolved", "unknown"}
        or values["direct_contradictions"] == "unknown"
    ):
        return _decision("not_evaluable", "direction_unknown")
    return _decision("pass", "pass")


def _g4(values: Mapping[str, str]) -> tuple[str, tuple[str, ...]]:
    support = {
        values["independent_pharmacology"],
        values["independent_rescue_lineage"],
    }
    if "copied_only" in support:
        return _decision("fail", "copied_evidence")
    if "pathway_proximity" in support:
        return _decision("fail", "insufficient_directness")
    if "missing" in support or values["translational_maturity"] == "computational_inference":
        return _decision("not_evaluable", "insufficient_directness")
    return _decision("pass", "pass")


TRUTH_TABLE = {"G0": _g0, "G1": _g1, "G2": _g2, "G3": _g3, "G4": _g4}


def evaluate_g0_g4(
    *,
    input_values: Mapping[str, str],
    evidence_ids_by_gate: Mapping[str, tuple[str, ...]],
    coverage_by_gate: Mapping[str, str],
    review_state_by_gate: Mapping[str, str],
    decision_policy: Mapping[str, object],
) -> G0G4Evaluation:
    if set(input_values) != REQUIRED_INPUTS:
        raise GateEvaluationError("G0-G4 inputs are incomplete; no field defaults to pass")
    for name, allowed in ALLOWED_VALUES.items():
        if input_values[name] not in allowed:
            raise GateEvaluationError(f"undefined gate input state: {name}")
    maturity = decision_policy.get("evidence_maturity")
    if not isinstance(maturity, list) or input_values["translational_maturity"] not in maturity:
        raise GateEvaluationError("undefined translational maturity")
    if not input_values["named_jurisdiction"]:
        raise GateEvaluationError("named jurisdiction is required")
    review = decision_policy.get("review")
    review_states = set(review.get("independence_states", ())) if isinstance(review, Mapping) else set()
    if (
        set(evidence_ids_by_gate) != set(GATE_IDS)
        or set(coverage_by_gate) != set(GATE_IDS)
        or set(review_state_by_gate) != set(GATE_IDS)
    ):
        raise GateEvaluationError("G0-G4 audit metadata is incomplete")
    policy_version = decision_policy.get("policy_version")
    if not isinstance(policy_version, str) or not policy_version:
        raise GateEvaluationError("decision policy version is missing")
    results: list[GateResult] = []
    for gate_id in GATE_IDS:
        evidence_ids = evidence_ids_by_gate[gate_id]
        coverage = coverage_by_gate[gate_id]
        review_state = review_state_by_gate[gate_id]
        if (
            not evidence_ids
            or coverage not in COVERAGE_STATES
            or review_state not in review_states
        ):
            raise GateEvaluationError(f"{gate_id} audit metadata is incomplete")
        outcome, reasons = TRUTH_TABLE[gate_id](input_values)
        allowed_reasons = set(decision_policy["gates"][gate_id]["reason_codes"])
        if set(reasons) - allowed_reasons:
            raise GateEvaluationError(f"{gate_id} emitted an undefined reason code")
        results.append(
            GateResult(
                gate_id=gate_id,
                outcome=outcome,
                reason_codes=reasons,
                evidence_ids=evidence_ids,
                coverage_state=coverage,
                policy_version=policy_version,
                review_state=review_state,
            )
        )
        if gate_id == "G0" and outcome == "fail":
            return G0G4Evaluation(tuple(results), stopped_after_g0=True)
    return G0G4Evaluation(tuple(results), stopped_after_g0=False)
