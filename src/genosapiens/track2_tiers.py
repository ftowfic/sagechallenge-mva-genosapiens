"""Exact candidate-tier transitions over complete G0-G7 gate results."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Iterable, Mapping

from genosapiens.track2_gates import GateResult
from genosapiens.track2_policy import GATE_OUTCOMES, GATES


class TierDecisionError(ValueError):
    """Raised when incomplete or undefined gates could otherwise imply a tier."""


@dataclass(frozen=True, slots=True)
class CandidateTierDecision:
    tier: str
    primary_shortlist_member: bool
    report_label: str
    decisive_gate_id: str | None
    decision_reasons: tuple[str, ...]
    policy_version: str


def _blocking_result(
    by_gate: Mapping[str, GateResult], precedence: Iterable[str]
) -> GateResult | None:
    for entry in precedence:
        gate_id, outcome = entry.split(":", 1)
        result = by_gate[gate_id]
        if result.outcome == outcome:
            return result
    return None


def decide_candidate_tier(
    gate_results: Iterable[GateResult], decision_policy: Mapping[str, object]
) -> CandidateTierDecision:
    results = tuple(gate_results)
    by_gate = {result.gate_id: result for result in results}
    if len(results) != len(GATES) or set(by_gate) != GATES:
        raise TierDecisionError("all G0-G7 results are required; no missing gate defaults to pass")
    policy_version = decision_policy.get("policy_version")
    if not isinstance(policy_version, str) or not policy_version:
        raise TierDecisionError("decision policy version is missing")
    for result in results:
        if (
            result.outcome not in GATE_OUTCOMES
            or result.policy_version != policy_version
            or not result.reason_codes
            or not result.evidence_ids
            or not result.coverage_state
            or not result.review_state
        ):
            raise TierDecisionError("gate result is incomplete or incompatible")
        allowed_reasons = set(decision_policy["gates"][result.gate_id]["reason_codes"])
        if set(result.reason_codes) - allowed_reasons:
            raise TierDecisionError("gate result uses an undefined reason code")
    blocking = _blocking_result(by_gate, decision_policy["blocking_precedence"])
    all_reasons = tuple(
        reason
        for gate_id in sorted(by_gate)
        for reason in by_gate[gate_id].reason_codes
        if reason != "pass"
    )
    if blocking is not None:
        tier = "contraindicated" if blocking.outcome == "contraindicating" else "excluded"
        return CandidateTierDecision(
            tier=tier,
            primary_shortlist_member=False,
            report_label=tier,
            decisive_gate_id=blocking.gate_id,
            decision_reasons=all_reasons,
            policy_version=policy_version,
        )
    bounded = {gate_id for gate_id, result in by_gate.items() if result.outcome == "bounded_uncertainty"}
    not_evaluable = {
        gate_id for gate_id, result in by_gate.items() if result.outcome == "not_evaluable"
    }
    if bounded:
        conditional = decision_policy["tier_transitions"]["conditional"]
        allowed_bounded = set(conditional["bounded_uncertainty_only"])
        required_pass = set(conditional["required_pass"])
        forbidden = set(conditional["forbidden_reasons"])
        bounded_reasons = {
            reason for gate_id in bounded for reason in by_gate[gate_id].reason_codes
        }
        if (
            not_evaluable
            or not bounded <= allowed_bounded
            or any(by_gate[gate_id].outcome != "pass" for gate_id in required_pass)
            or bounded_reasons != {"bounded_secondary_uncertainty"}
            or set(all_reasons) & forbidden
        ):
            raise TierDecisionError("bounded uncertainty does not satisfy conditional policy")
        return CandidateTierDecision(
            tier="conditional",
            primary_shortlist_member=True,
            report_label="conditional mechanistically supported research hypothesis; not evidence of efficacy",
            decisive_gate_id=None,
            decision_reasons=all_reasons,
            policy_version=policy_version,
        )
    if not_evaluable:
        return CandidateTierDecision(
            tier="exploratory",
            primary_shortlist_member=False,
            report_label="exploratory research hypothesis; not evidence of efficacy",
            decisive_gate_id=sorted(not_evaluable)[0],
            decision_reasons=all_reasons,
            policy_version=policy_version,
        )
    if any(result.outcome != "pass" for result in results):
        raise TierDecisionError("gate outcomes do not match a configured tier transition")
    return CandidateTierDecision(
        tier="supported",
        primary_shortlist_member=True,
        report_label="mechanistically supported research hypothesis; not evidence of efficacy",
        decisive_gate_id=None,
        decision_reasons=(),
        policy_version=policy_version,
    )
