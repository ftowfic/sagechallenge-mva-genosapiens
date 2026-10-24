"""Route-completeness milestone and deferred computational-method boundary."""

from __future__ import annotations

from dataclasses import asdict, dataclass
from datetime import date, datetime
import hashlib
import json
from pathlib import Path
from typing import Mapping

from genosapiens.track2_routes import RoutePortfolio
from genosapiens.track2_gates import GateResult
from genosapiens.track2_records import validate_track2_record, load_record_registry
from genosapiens.track2_tiers import decide_candidate_tier


class MilestoneError(ValueError):
    """Raised when route coverage or computational provenance is incomplete."""


RESEARCH_DOMAINS = {
    "approval", "mechanism", "drug_action", "current_warnings", "negative_evidence",
    "harmful_alternatives", "target_liabilities", "claims", "validation",
}
TRANSLATIONAL_GAPS = {"target_site_exposure", "rescue_magnitude", "pediatric_applicability"}


@dataclass(frozen=True, slots=True)
class ResearchGap:
    field: str
    state: str
    resolving_observation: str
    stopping_condition: str

    def __post_init__(self) -> None:
        if (self.field not in TRANSLATIONAL_GAPS
                or self.state not in {"unknown", "not_assessed"}
                or not self.resolving_observation.strip() or not self.stopping_condition.strip()):
            raise MilestoneError("research gap needs an explicit resolution and stopping condition")


@dataclass(frozen=True, slots=True)
class ResearchSubmissionReview:
    review_id: str
    reviewer_role: str
    independence_state: str
    reviewed_at: str
    evidence_by_domain: Mapping[str, tuple[str, ...]]
    evidence_records: Mapping[str, Mapping[str, object]]
    gaps: tuple[ResearchGap, ...]


@dataclass(frozen=True, slots=True)
class ResearchEligibilityDecision:
    candidate_id: str
    eligible: bool
    scientific_tier: str
    reason_codes: tuple[str, ...]
    gaps: tuple[ResearchGap, ...]
    evidence_ids: tuple[str, ...]
    review_id: str | None
    review_state: str
    policy_version: str
    input_sha256: str


def public_research_decision(decision: ResearchEligibilityDecision | Mapping[str, object]) -> dict[str, object]:
    """Projection excludes private review identities and input fingerprints."""
    value = asdict(decision) if isinstance(decision, ResearchEligibilityDecision) else decision
    return {key: value[key] for key in ("eligible", "scientific_tier", "reason_codes", "gaps")}


def assess_research_eligibility(
    candidate_id: str,
    gate_results: tuple[GateResult, ...],
    review: ResearchSubmissionReview | None,
    decision_policy: Mapping[str, object],
) -> ResearchEligibilityDecision:
    """Replay research membership from gates and reviewed source-bound records.

    This adds report membership; it never overrides any scientific gate or tier.
    The frozen synthesis must also verify record/source bytes and their manifests.
    """
    tier = decide_candidate_tier(gate_results, decision_policy)
    policy = decision_policy["research_submission"]
    gates = {result.gate_id: result for result in gate_results}
    reviewed_states = set(decision_policy["review"]["independence_states"]) - {"not_reviewed"}
    reasons = []
    if not candidate_id:
        raise MilestoneError("research candidate identity is required")
    for gate_id in policy["required_pass"]:
        gate = gates[gate_id]
        if (gate.outcome != "pass" or gate.review_state not in reviewed_states
                or gate.coverage_state not in {"complete", "success"}):
            reasons.append("core_gate_not_satisfied:" + gate_id)
    if gates["G2"].outcome not in {"pass", "bounded_uncertainty"}:
        reasons.append("core_gate_not_satisfied:G2")
    if (gates["G2"].review_state not in reviewed_states
            or gates["G2"].coverage_state not in {"complete", "success"}):
        reasons.append("required_mechanism_review_missing")
    if tier.tier in {"excluded", "contraindicated"}:
        reasons.append("blocking_scientific_finding")
    if gates["G5"].outcome == "not_evaluable" and (
        not set(gates["G5"].reason_codes) <= set(policy["allowed_exploratory_reasons"])
    ):
        reasons.append("unsupported_exploratory_gap")
    evidence_ids = set()
    gaps = ()
    if review is None:
        reasons.append("required_research_review_missing")
    else:
        try:
            reviewed_at = datetime.fromisoformat(review.reviewed_at.replace("Z", "+00:00"))
            if (not review.review_id or reviewed_at.tzinfo is None
                    or review.reviewer_role not in decision_policy["review"]["roles"]
                    or review.independence_state not in {"independent", "self_review_disclosed"}):
                raise MilestoneError("research review audit is incomplete")
            if set(review.evidence_by_domain) != RESEARCH_DOMAINS:
                reasons.append("required_research_evidence_missing")
            root = Path(__file__).resolve().parents[2]
            registry = load_record_registry(root / "config/track2/record-registry-v1.json")
            schema = json.loads((root / "schemas/track2/record-envelope-v1.schema.json").read_text())
            for domain, ids in review.evidence_by_domain.items():
                if not ids or len(ids) != len(set(ids)):
                    reasons.append("required_research_evidence_missing:" + domain)
                for evidence_id in ids:
                    record = review.evidence_records.get(evidence_id)
                    if record is None:
                        reasons.append("unresolved_research_evidence:" + domain)
                        continue
                    validate_track2_record(record, schema, registry)
                    if (record["record_type"] != "evidence_record"
                            or record["record_id"] != evidence_id
                            or record["missing_state"] != "present"
                            or not record["references"]
                            or any(reference["manifest_sha256"] is None
                                   or reference["role"] != "source"
                                   or reference["target_type"] != "source_snapshot"
                                   for reference in record["references"])):
                        reasons.append("unusable_research_evidence:" + domain)
                    evidence_ids.add(evidence_id)
            role_gates = {"approval": "G1", "mechanism": "G2", "drug_action": "G3",
                          "current_warnings": "G6", "validation": "G7"}
            for domain, gate_id in role_gates.items():
                if not set(review.evidence_by_domain.get(domain, ())) & set(gates[gate_id].evidence_ids):
                    reasons.append("review_gate_evidence_mismatch:" + domain)
            sources = json.loads((root / "config/track2/source-domain-matrix-v1.json").read_text())
            for domain, source_domain in (("current_warnings", "current_label"), ("approval", "approval_marketing")):
                requirement = sources["domains"][source_domain]
                for evidence_id in review.evidence_by_domain.get(domain, ()):
                    record = review.evidence_records.get(evidence_id)
                    if record is None:
                        continue
                    context = record["payload"]["context"]
                    try:
                        accessed = date.fromisoformat(context["accessed_on"])
                    except (KeyError, TypeError, ValueError):
                        reasons.append("current_authoritative_source_missing:" + domain)
                        continue
                    if (not context.get("document_version")
                            or context.get("authority") not in requirement["primary"] + requirement["alternate"]
                            or not 0 <= (reviewed_at.date() - accessed).days <= requirement["freshness_days"]):
                        reasons.append("current_authoritative_source_missing:" + domain)
            for evidence_id in review.evidence_by_domain.get("mechanism", ()):
                record = review.evidence_records.get(evidence_id)
                if record is not None and not (
                    record["payload"]["applicability"] == "exact_subtype" or (
                        record["payload"]["applicability"] == "shared_mechanism"
                        and record["payload"]["context"].get("transferability_rationale")
                        and record["payload"]["context"].get("transferability_evidence_ids")
                    )
                ):
                    reasons.append("exact_or_transferable_mechanism_missing")
            def lineages(domain):
                return {
                    review.evidence_records[evidence_id]["payload"]["lineage_key"]
                    for evidence_id in review.evidence_by_domain.get(domain, ())
                    if evidence_id in review.evidence_records
                    and review.evidence_records[evidence_id]["payload"]["independence"] == "independent"
                    and review.evidence_records[evidence_id]["payload"]["lineage_key"]
                }
            if not any(action != mechanism for action in lineages("drug_action")
                       for mechanism in lineages("mechanism")):
                reasons.append("independent_directional_support_missing")
            for domain in ("mechanism", "drug_action"):
                for evidence_id in review.evidence_by_domain.get(domain, ()):
                    record = review.evidence_records.get(evidence_id)
                    if record is not None and record["payload"]["result_direction"] not in {"increase", "decrease"}:
                        reasons.append("signed_support_missing:" + domain)
            gaps = tuple(sorted(review.gaps, key=lambda gap: gap.field))
            if len({gap.field for gap in gaps}) != len(gaps):
                raise MilestoneError("research gaps must be unique")
            required_gaps = {reason.removesuffix("_unknown")
                for reason in gates["G5"].reason_codes if reason.endswith("_unknown")}
            if required_gaps - {gap.field for gap in gaps}:
                reasons.append("translational_gap_followup_missing")
        except (KeyError, TypeError, ValueError) as error:
            raise MilestoneError("research review or source records are invalid") from error
    semantic = {"candidate_id": candidate_id, "gates": [asdict(gates[key]) for key in sorted(gates)],
                "review": asdict(review) if review else None, "policy": decision_policy,
                "source_requirements": sources if review else None}
    digest = hashlib.sha256(json.dumps(semantic, sort_keys=True, separators=(",", ":"),
                            allow_nan=False).encode()).hexdigest()
    return ResearchEligibilityDecision(
        candidate_id, not reasons, tier.tier, tuple(sorted(set(reasons))) or ("research_evidence_contract_satisfied",),
        gaps, tuple(sorted(evidence_ids)), review.review_id if review else None,
        review.independence_state if review else "not_reviewed", tier.policy_version, digest,
    )


@dataclass(frozen=True, slots=True)
class RouteMilestoneDecision:
    final_prioritization_allowed: bool
    organizer_milestone_allowed: bool
    route_coverage_claim: str
    no_candidate_claim: str
    blocking_routes: tuple[tuple[str, str, str], ...]
    coverage_limited_routes: tuple[tuple[str, str, str], ...]
    candidate_tier_ceilings: tuple[tuple[str, str], ...]
    deferred_routes: tuple[tuple[str, str, str], ...] = ()
    public_scope_sha256: str | None = None
    selected_routes: tuple[str, ...] = ()
    analysis_complete: bool = False


@dataclass(frozen=True, slots=True)
class ExternalComputationalHypothesis:
    method: str
    metadata: tuple[tuple[str, object], ...]
    tier_ceiling: str
    sole_primary_support_allowed: bool
    score_substitution_used: bool


def assess_route_milestone(
    portfolio: RoutePortfolio,
    *,
    candidate_route_dependencies: Mapping[str, tuple[tuple[str, str], ...]],
    decision_policy: Mapping[str, object],
) -> RouteMilestoneDecision:
    if not portfolio.all_routes_assessed:
        raise MilestoneError("every route requires an applicability decision")
    policy = decision_policy.get("route_completeness")
    if not isinstance(policy, Mapping):
        raise MilestoneError("route completeness policy is missing")
    states = {
        (item.objective_id, item.route): (
            "not_applicable" if item.attempt is None else item.attempt.state
        )
        for item in portfolio.assessments
        if item.route in portfolio.selected_routes
    }
    blocking_states = set(policy["milestone_blocking_states"])
    limited_states = set(policy["coverage_limiting_states"])
    blocking = tuple(
        sorted(
            (objective_id, route, state)
            for (objective_id, route), state in states.items()
            if state in blocking_states
        )
    )
    limited = tuple(
        sorted(
            (objective_id, route, state)
            for (objective_id, route), state in states.items()
            if state in limited_states
        )
    )
    ceilings: list[tuple[str, str]] = []
    for candidate_id, dependencies in candidate_route_dependencies.items():
        if not candidate_id or not dependencies or any(key not in states for key in dependencies):
            raise MilestoneError("candidate route dependency is unresolved")
        if all(states[key] in limited_states for key in dependencies):
            ceilings.append((candidate_id, policy["affected_candidate_tier_ceiling"]))
    allowed = not blocking and portfolio.all_applicable_routes_attempted
    return RouteMilestoneDecision(
        final_prioritization_allowed=allowed,
        organizer_milestone_allowed=allowed,
        route_coverage_claim=(policy["limited_coverage_claim"] if limited else "complete"),
        no_candidate_claim=(
            "not_exhaustive" if limited or portfolio.deferred_routes else "bounded_exhaustive"
        ),
        blocking_routes=blocking,
        coverage_limited_routes=limited,
        candidate_tier_ceilings=tuple(sorted(ceilings)),
        deferred_routes=portfolio.deferred_routes,
        public_scope_sha256=portfolio.public_scope_sha256,
        selected_routes=portfolio.selected_routes,
        analysis_complete=portfolio.analysis_complete,
    )


def assess_computational_candidate_request(
    *,
    method: str,
    externally_supplied: bool,
    metadata: Mapping[str, object],
    score_substitution: float | None,
    decision_policy: Mapping[str, object],
) -> ExternalComputationalHypothesis:
    constraints = decision_policy.get("candidate_generation_constraints")
    if not isinstance(constraints, Mapping):
        raise MilestoneError("candidate-generation constraints are missing")
    prohibited = set(constraints["prohibited_internal_methods"])
    if method not in prohibited:
        raise MilestoneError("computational candidate method is not defined")
    if score_substitution is not None or constraints["score_substitution_allowed"] is not False:
        raise MilestoneError("computational score substitution is prohibited")
    if not externally_supplied:
        raise MilestoneError("deferred computational candidate generation is rejected")
    required = set(constraints["required_external_metadata"])
    if set(metadata) != required or any(
        value is None or value == "" for value in metadata.values()
    ):
        raise MilestoneError("external computational hypothesis metadata is incomplete")
    return ExternalComputationalHypothesis(
        method=method,
        metadata=tuple(sorted(metadata.items())),
        tier_ceiling=constraints["external_computational_tier_ceiling"],
        sole_primary_support_allowed=not constraints["cannot_be_sole_primary_support"],
        score_substitution_used=False,
    )
