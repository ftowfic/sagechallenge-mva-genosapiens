"""Reviewed intervention objectives and signed target hypotheses."""

from __future__ import annotations

from dataclasses import dataclass

from genosapiens.track2_endpoint import EndpointSelection
from genosapiens.track2_mechanism import MechanismModel


class TargetHypothesisError(ValueError):
    """Raised when a target lacks a reviewed, signed path to rescue."""


@dataclass(frozen=True, slots=True)
class InterventionObjective:
    objective_id: str
    mechanism_scenario_id: str
    endpoint_id: str
    desired_change: str
    direction: str
    compartment: str
    tissue: str
    stage: str
    assumptions: tuple[str, ...]
    harm_conditions: tuple[str, ...]
    mechanism_review_id: str
    endpoint_review_id: str
    evidence_ids: tuple[str, ...]

    def as_track2_payload(self) -> dict:
        return {
            "assumptions": list(self.assumptions),
            "compartment": self.compartment,
            "desired_change": self.desired_change,
            "direction": self.direction,
            "harm_conditions": list(self.harm_conditions),
            "stage": self.stage,
            "tissue": self.tissue,
        }


@dataclass(frozen=True, slots=True)
class SignedPathEdge:
    from_id: str
    to_id: str
    sign: int
    evidence_ids: tuple[str, ...]

    def __post_init__(self) -> None:
        if (
            not self.from_id
            or not self.to_id
            or self.from_id == self.to_id
            or self.sign not in {-1, 1}
            or not self.evidence_ids
        ):
            raise TargetHypothesisError("rescue path edges must be signed and evidenced")


@dataclass(frozen=True, slots=True)
class TargetHypothesis:
    target_id: str
    causal_role: str
    required_modulation: str
    signed_rescue_path: tuple[SignedPathEdge, ...]
    tissue: str
    stage: str
    alternatives: tuple[str, ...]
    conflicts: tuple[str, ...]
    harm_conditions: tuple[str, ...]
    evidence_ids: tuple[str, ...]
    objective_id: str
    ready_for_drug_retrieval: bool

    def as_track2_payload(self) -> dict:
        return {
            "alternatives": list(self.alternatives),
            "causal_role": self.causal_role,
            "conflicts": list(self.conflicts),
            "harm_conditions": list(self.harm_conditions),
            "required_modulation": self.required_modulation,
            "signed_rescue_path": [
                {
                    "from_id": edge.from_id,
                    "to_id": edge.to_id,
                    "sign": edge.sign,
                    "evidence_ids": list(edge.evidence_ids),
                }
                for edge in self.signed_rescue_path
            ],
            "stage": self.stage,
            "target_id": self.target_id,
            "tissue": self.tissue,
        }


def derive_intervention_objective(
    *,
    objective_id: str,
    mechanism: MechanismModel,
    endpoint: EndpointSelection,
    desired_change: str,
    compartment: str,
    assumptions: tuple[str, ...],
    harm_conditions: tuple[str, ...],
    mechanism_review_id: str,
    endpoint_review_id: str,
    evidence_ids: tuple[str, ...],
) -> InterventionObjective:
    if (
        not objective_id
        or not desired_change
        or endpoint.expected_direction == "unknown"
        or not compartment
        or not endpoint.tissue
        or not endpoint.stage
        or not assumptions
        or not harm_conditions
        or not mechanism_review_id
        or not endpoint_review_id
        or not evidence_ids
        or endpoint.relevance != "supported"
    ):
        raise TargetHypothesisError(
            "intervention objective requires reviewed mechanism and endpoint evidence"
        )
    return InterventionObjective(
        objective_id=objective_id,
        mechanism_scenario_id=mechanism.scenario_id,
        endpoint_id=endpoint.endpoint_id,
        desired_change=desired_change,
        direction=endpoint.expected_direction,
        compartment=compartment,
        tissue=endpoint.tissue,
        stage=endpoint.stage,
        assumptions=assumptions,
        harm_conditions=harm_conditions,
        mechanism_review_id=mechanism_review_id,
        endpoint_review_id=endpoint_review_id,
        evidence_ids=evidence_ids,
    )


def build_target_hypothesis(
    *,
    objective: InterventionObjective,
    target_id: str,
    causal_role: str,
    required_modulation: str,
    signed_rescue_path: tuple[SignedPathEdge, ...],
    tissue: str,
    stage: str,
    alternatives: tuple[str, ...],
    conflicts: tuple[str, ...],
    harm_conditions: tuple[str, ...],
    evidence_ids: tuple[str, ...],
) -> TargetHypothesis:
    if (
        not target_id
        or causal_role not in {"direct", "indirect"}
        or required_modulation not in {"decrease", "increase", "partial_decrease", "partial_increase"}
        or not signed_rescue_path
        or tissue != objective.tissue
        or stage != objective.stage
        or not evidence_ids
    ):
        raise TargetHypothesisError("target hypothesis is incomplete")
    if signed_rescue_path[0].from_id != target_id:
        raise TargetHypothesisError("signed rescue path does not begin at the target")
    if signed_rescue_path[-1].to_id != objective.endpoint_id:
        raise TargetHypothesisError("signed rescue path does not reach the selected endpoint")
    if any(
        left.to_id != right.from_id
        for left, right in zip(signed_rescue_path, signed_rescue_path[1:])
    ):
        raise TargetHypothesisError("signed rescue path is not contiguous")
    sign = -1 if required_modulation in {"decrease", "partial_decrease"} else 1
    for edge in signed_rescue_path:
        sign *= edge.sign
    if objective.direction not in {"increase", "decrease"} or sign != (1 if objective.direction == "increase" else -1):
        raise TargetHypothesisError("target modulation does not produce the required endpoint direction")
    if causal_role == "indirect" and (
        len(signed_rescue_path) < 2 or not alternatives or not harm_conditions
    ):
        raise TargetHypothesisError(
            "indirect targets require a signed path, alternatives, and harm conditions"
        )
    return TargetHypothesis(
        target_id=target_id,
        causal_role=causal_role,
        required_modulation=required_modulation,
        signed_rescue_path=signed_rescue_path,
        tissue=tissue,
        stage=stage,
        alternatives=alternatives,
        conflicts=conflicts,
        harm_conditions=harm_conditions,
        evidence_ids=evidence_ids,
        objective_id=objective.objective_id,
        ready_for_drug_retrieval=True,
    )


def require_drug_retrieval_ready(hypothesis: TargetHypothesis) -> None:
    if not hypothesis.ready_for_drug_retrieval:
        raise TargetHypothesisError("target hypothesis is not ready for drug retrieval")
