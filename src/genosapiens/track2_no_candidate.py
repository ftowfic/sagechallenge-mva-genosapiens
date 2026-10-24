"""Complete exclusions and honest no-primary-candidate objective outcomes."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Mapping

from genosapiens.track2_gates import GateResult
from genosapiens.track2_tiers import CandidateTierDecision
from genosapiens.track2_milestone import (
    ResearchSubmissionReview, ResearchEligibilityDecision, assess_research_eligibility,
)


class NoCandidateError(ValueError):
    """Raised when failed candidates or no-candidate evidence are omitted."""


def native_objective_outcomes(plan,result,assessment):
    statuses={r['objective_id']:r['status'] for r in assessment.objectives}
    return [dict(objective_id=o['objective_id'],assessment_status=statuses[o['objective_id']],
        discovery_candidate_ids=sorted({r['canonical_active_ingredient_id'] for r in result['concordance'] if r['objective_id']==o['objective_id']}),
        paired_research_candidate_ids=sorted({r['canonical_active_ingredient_id'] for r in result['concordance']
            if r['objective_id']==o['objective_id'] and r['paired_research_eligible']}),
        clinical_tier='unassessed',submission_ready=False) for o in plan['objectives']]


@dataclass(frozen=True, slots=True)
class CandidateExclusion:
    active_ingredient_id: str
    tier: str
    decisive_gate_ids: tuple[str, ...]
    reason_codes: tuple[str, ...]
    evidence_gaps: tuple[str, ...]
    resolution_experiments: tuple[str, ...]


@dataclass(frozen=True, slots=True)
class ObjectiveCandidateOutcome:
    objective_id: str
    primary_candidate_ids: tuple[str, ...]
    exclusions: tuple[CandidateExclusion, ...]
    track2_objective_status: str
    submission_ready: bool
    research_candidate_ids: tuple[str, ...] = ()
    research_decisions: tuple[ResearchEligibilityDecision, ...] = ()
    scientific_outcome: str = "not_assessed"
    assessment_complete: bool = False
    category_research_decisions: tuple[Mapping[str, object], ...] = ()
    assessment_status: str = 'not_performed'
    assessment_reasons: tuple[Mapping[str, object], ...] = ()


def include_lincs_research_membership(outcomes, result):
    """Join research membership; retain scientific exclusions, gates and tiers."""
    from dataclasses import replace
    updated = []
    for outcome in outcomes:
        decisions = tuple(row for row in result['candidate_decisions']
                          if outcome.objective_id in row['objective_ids'])
        ids = tuple(sorted(set(outcome.research_candidate_ids) |
                           {row['candidate_id'] for row in decisions if row['research_eligible']}))
        updated.append(replace(outcome, research_candidate_ids=ids, category_research_decisions=decisions,
            track2_objective_status='research_hypothesis_identified' if ids else outcome.track2_objective_status))
    return tuple(updated)


def build_objective_candidate_outcome(
    *,
    objective_id: str,
    decisions: Mapping[str, CandidateTierDecision],
    gate_results: Mapping[str, tuple[GateResult, ...]],
    evidence_gaps: Mapping[str, tuple[str, ...]],
    resolution_experiments: Mapping[str, tuple[str, ...]],
    other_submission_requirements_met: bool,
    research_reviews: Mapping[str, ResearchSubmissionReview | None] | None = None,
    decision_policy: Mapping[str, object] | None = None,
    assessment: Mapping[str, object],
) -> ObjectiveCandidateOutcome:
    from genosapiens.track2_assessment import assess_objective
    replayed = assess_objective(assessment['objective_id'], assessment['scenario_id'],
        assessment['supported'], assessment['routes'])
    if replayed != assessment or objective_id != assessment['objective_id']:
        raise NoCandidateError('Objective assessment disagrees with retained search evidence.')
    assessment_complete = assessment['status'] == 'complete'
    candidate_ids = set(decisions)
    if (
        not objective_id
        or set(gate_results) != candidate_ids
        or set(evidence_gaps) != candidate_ids
        or set(resolution_experiments) != candidate_ids
    ):
        raise NoCandidateError("candidate disposition inputs are incomplete")
    primary: list[str] = []
    exclusions: list[CandidateExclusion] = []
    for candidate_id in sorted(candidate_ids):
        decision = decisions[candidate_id]
        if decision.primary_shortlist_member:
            if decision.tier not in {"supported", "conditional"}:
                raise NoCandidateError("primary shortlist contains a non-primary tier")
            primary.append(candidate_id)
            continue
        if decision.tier not in {"contraindicated", "excluded", "exploratory"}:
            raise NoCandidateError("candidate has an undefined non-primary disposition")
        nonpassing = tuple(
            sorted(
                result.gate_id
                for result in gate_results[candidate_id]
                if result.outcome not in {"pass", "bounded_uncertainty"}
            )
        )
        reasons = tuple(
            reason
            for result in sorted(gate_results[candidate_id], key=lambda item: item.gate_id)
            for reason in result.reason_codes
            if reason != "pass"
        )
        if (
            not nonpassing
            or not reasons
            or not evidence_gaps[candidate_id]
            or not resolution_experiments[candidate_id]
        ):
            raise NoCandidateError(
                "every non-primary candidate requires decisive gates, gaps, and experiments"
            )
        exclusions.append(
            CandidateExclusion(
                active_ingredient_id=candidate_id,
                tier=decision.tier,
                decisive_gate_ids=nonpassing,
                reason_codes=reasons,
                evidence_gaps=evidence_gaps[candidate_id],
                resolution_experiments=resolution_experiments[candidate_id],
            )
        )
    research_decisions = ()
    if research_reviews is not None:
        if decision_policy is None or set(research_reviews) != candidate_ids:
            raise NoCandidateError("research eligibility requires every candidate and the decision policy")
        research_decisions = tuple(assess_research_eligibility(
            candidate_id, gate_results[candidate_id], research_reviews[candidate_id], decision_policy,
        ) for candidate_id in sorted(candidate_ids))
        if any(decision.scientific_tier != decisions[decision.candidate_id].tier
               for decision in research_decisions):
            raise NoCandidateError("research eligibility disagrees with the scientific decision replay")
    research_ids = tuple(decision.candidate_id for decision in research_decisions if decision.eligible)
    return ObjectiveCandidateOutcome(
        objective_id=objective_id,
        primary_candidate_ids=tuple(primary),
        exclusions=tuple(exclusions),
        track2_objective_status=(
            "research_hypothesis_identified" if research_ids else (
                "no_plausible_candidate_identified" if assessment_complete else 'assessment_' + assessment['status']
            )
        ),
        submission_ready=(assessment_complete and bool(research_ids) and other_submission_requirements_met),
        research_candidate_ids=research_ids,
        research_decisions=research_decisions,
        scientific_outcome=("supported_or_conditional_candidate_identified" if primary else
                            "no_supported_or_conditional_candidate" if assessment_complete else
                            'assessment_' + assessment['status']),
        assessment_complete=assessment_complete,
        assessment_status=assessment['status'],
        assessment_reasons=tuple(assessment['reasons']),
    )
