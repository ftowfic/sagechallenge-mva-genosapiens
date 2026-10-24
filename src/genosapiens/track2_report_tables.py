"""Deterministic machine-readable and Markdown Track 2 evidence tables."""

from __future__ import annotations

from dataclasses import asdict, dataclass, fields, is_dataclass
import json
import re
from typing import Any, Mapping

from genosapiens.track2_candidate_union import CanonicalCandidate
from genosapiens.track2_coverage import CoverageRecord
from genosapiens.track2_endpoint import EndpointSelection
from genosapiens.track2_exposure import ActivityExposureComparison
from genosapiens.track2_g6 import G6SafetyReview
from genosapiens.track2_gates import GateResult
from genosapiens.track2_input import AnalysisOrigin
from genosapiens.track2_mechanism import MechanismReadiness
from genosapiens.track2_no_candidate import CandidateExclusion
from genosapiens.track2_patient_context import OptionalPatientContext
from genosapiens.track2_population_applicability import PopulationApplicability
from genosapiens.track2_policy import GATES
from genosapiens.track2_target_hypothesis import InterventionObjective, TargetHypothesis
from genosapiens.track2_tiers import CandidateTierDecision, decide_candidate_tier
from genosapiens.track2_milestone import (
    ResearchEligibilityDecision, ResearchSubmissionReview, ResearchGap, assess_research_eligibility,
)


REPORT_TABLE_SCHEMA = "genosapiens.track2-report-tables/v1"


def exploratory_tables(plan,result):
    """Keep native nomination flags independent of clinical gates and tiers."""
    from genosapiens.track2_mechanism import assess_mechanism_readiness
    return dict(schema='genosapiens.track2-exploratory-tables/v1',
        scenarios=[dict(scenario_id=r['scenario_id'],gene_id=r['gene_id'],
            mechanism_acceptance=asdict(assess_mechanism_readiness(r['mechanism_inputs'])),
            eligibility=r['discovery_eligibility']) for r in plan['mechanism_assessment']['scenarios']],
        objectives=plan['exploratory_discovery']['objectives'],routes=result['route_outcomes'],
        candidates=[dict(**row,gates={f'G{i}':'unassessed' for i in range(8)},clinical_tier='unassessed')
            for row in result['concordance']])
TABLE_NAMES = (
    "mechanism_readiness",
    "endpoints",
    "target_hypotheses",
    "candidates",
    "exclusions",
    "gates",
)
SUBTYPE_APPLICABILITY_STATES = {
    "exact_subtype",
    "other_subtype",
    "shared_mechanism",
    "umbrella",
    "unresolved",
}
EVIDENCE_DIMENSIONS = {
    "applicability",
    "directness",
    "exact_subtype",
    "exposure_feasibility",
    "mechanism_readiness",
    "pediatric_practicality",
    "risk_of_bias",
    "robustness",
    "safety",
    "source_independence",
    "target_liability",
    "translational_maturity",
}


class ReportTableError(ValueError):
    """Raised when a report table would contain an unresolved or drifting row."""


def _require_text(value: object, label: str) -> str:
    if not isinstance(value, str) or not value:
        raise ReportTableError(f"{label} is required")
    return value


def _require_text_tuple(values: tuple[str, ...], label: str) -> None:
    if (
        not values
        or len(values) != len(set(values))
        or any(not isinstance(value, str) or not value for value in values)
    ):
        raise ReportTableError(f"{label} must be a non-empty unique list")


@dataclass(frozen=True, slots=True)
class ReportEvidenceReference:
    evidence_id: str
    source_identity: str
    citation: str
    provenance: str

    def __post_init__(self) -> None:
        if ":" not in self.evidence_id or any(
            not isinstance(value, str) or not value
            for value in (
                self.evidence_id,
                self.source_identity,
                self.citation,
                self.provenance,
            )
        ):
            raise ReportTableError("report evidence reference is incomplete")


@dataclass(frozen=True, slots=True)
class ReportAssessment:
    state: str
    summary: str
    evidence_ids: tuple[str, ...]

    def __post_init__(self) -> None:
        _require_text(self.state, "assessment state")
        _require_text(self.summary, "assessment summary")
        _require_text_tuple(self.evidence_ids, "assessment evidence IDs")


@dataclass(frozen=True, slots=True)
class EvidenceLinkedStatement:
    statement: str
    evidence_ids: tuple[str, ...]

    def __post_init__(self) -> None:
        _require_text(self.statement, "evidence-linked statement")
        _require_text_tuple(self.evidence_ids, "statement evidence IDs")


@dataclass(frozen=True, slots=True)
class ApprovalContext:
    jurisdiction: str
    state: str
    product_ids: tuple[str, ...]
    evidence_ids: tuple[str, ...]

    def __post_init__(self) -> None:
        _require_text(self.jurisdiction, "approval jurisdiction")
        _require_text(self.state, "approval state")
        _require_text_tuple(self.product_ids, "approval product IDs")
        _require_text_tuple(self.evidence_ids, "approval evidence IDs")


@dataclass(frozen=True, slots=True)
class ResolvedReportContext:
    origin: AnalysisOrigin
    subtype_applicability: str
    translational_maturity: str
    effect_magnitude: str
    child_applicability: str
    general_mva_applicability: str
    evidence_ids: tuple[str, ...]
    policy_version: str
    exposure: ReportAssessment
    safety: ReportAssessment
    decision_reasons: tuple[str, ...]

    def __post_init__(self) -> None:
        for value, label in (
            (self.origin.scenario_id, "scenario ID"),
            (self.origin.input_origin, "Input origin"),
            (self.origin.source_sha256, "Input source hash"),
            (self.translational_maturity, "translational maturity"),
            (self.effect_magnitude, "effect magnitude"),
            (self.child_applicability, "child applicability"),
            (self.general_mva_applicability, "general MVA applicability"),
            (self.policy_version, "policy version"),
        ):
            _require_text(value, label)
        if self.subtype_applicability not in SUBTYPE_APPLICABILITY_STATES:
            raise ReportTableError("subtype applicability is undefined")
        _require_text_tuple(self.evidence_ids, "row-context evidence IDs")
        _require_text_tuple(self.decision_reasons, "row decision reasons")


@dataclass(frozen=True, slots=True)
class MechanismTableRecord:
    context: ResolvedReportContext
    readiness: MechanismReadiness


@dataclass(frozen=True, slots=True)
class EndpointTableRecord:
    context: ResolvedReportContext
    scenario_id: str
    endpoint: EndpointSelection


@dataclass(frozen=True, slots=True)
class TargetHypothesisTableRecord:
    context: ResolvedReportContext
    objective: InterventionObjective
    target: TargetHypothesis


@dataclass(frozen=True, slots=True)
class CandidateTableRecord:
    context: ResolvedReportContext
    scenario_id: str
    candidate: CanonicalCandidate
    selected_endpoint_id: str
    target_hypothesis_ids: tuple[str, ...]
    approval_contexts: tuple[ApprovalContext, ...]
    proposed_actions: tuple[EvidenceLinkedStatement, ...]
    route_coverage: tuple[CoverageRecord, ...]
    evidence_dimensions: Mapping[str, str]
    mechanism_alternative_robustness: EvidenceLinkedStatement
    target_liabilities: tuple[EvidenceLinkedStatement, ...]
    reachable_polypharmacology: tuple[EvidenceLinkedStatement, ...]
    current_management: tuple[EvidenceLinkedStatement, ...]
    negative_evidence: tuple[EvidenceLinkedStatement, ...]
    exposure_comparison: ActivityExposureComparison | None
    safety_review: G6SafetyReview
    population_applicability: PopulationApplicability
    patient_context: OptionalPatientContext
    conflicts: tuple[str, ...]
    uncertainties: tuple[str, ...]
    gate_results: tuple[GateResult, ...]
    tier_decision: CandidateTierDecision
    exclusion: CandidateExclusion | None
    display_order: int
    research_submission_review: ResearchSubmissionReview | None = None
    research_eligibility: ResearchEligibilityDecision | None = None


@dataclass(frozen=True, slots=True)
class ReportTableInputs:
    evidence_catalog: tuple[ReportEvidenceReference, ...]
    mechanism_records: tuple[MechanismTableRecord, ...]
    endpoint_records: tuple[EndpointTableRecord, ...]
    target_records: tuple[TargetHypothesisTableRecord, ...]
    candidate_records: tuple[CandidateTableRecord, ...]


@dataclass(frozen=True, slots=True)
class ReportTableRow:
    row_id: str
    origin: AnalysisOrigin
    subtype_applicability: str
    translational_maturity: str
    effect_magnitude: str
    child_applicability: str
    general_mva_applicability: str
    evidence_ids: tuple[str, ...]
    policy_version: str
    exposure: ReportAssessment
    safety: ReportAssessment
    decision_reasons: tuple[str, ...]
    details: Mapping[str, Any]

    def as_dict(self) -> dict[str, Any]:
        return {
            "row_id": self.row_id,
            "analysis_origin": _plain(self.origin),
            "subtype_applicability": self.subtype_applicability,
            "translational_maturity": self.translational_maturity,
            "effect_magnitude": self.effect_magnitude,
            "child_applicability": self.child_applicability,
            "general_mva_applicability": self.general_mva_applicability,
            "evidence_ids": list(self.evidence_ids),
            "policy_version": self.policy_version,
            "exposure": _plain(self.exposure),
            "safety": _plain(self.safety),
            "decision_reasons": list(self.decision_reasons),
            "details": _plain(self.details),
        }


@dataclass(frozen=True, slots=True)
class Track2ReportTables:
    evidence_catalog: tuple[ReportEvidenceReference, ...]
    policy_versions: tuple[str, ...]
    mechanism_readiness: tuple[ReportTableRow, ...]
    endpoints: tuple[ReportTableRow, ...]
    target_hypotheses: tuple[ReportTableRow, ...]
    candidates: tuple[ReportTableRow, ...]
    exclusions: tuple[ReportTableRow, ...]
    gates: tuple[ReportTableRow, ...]

    def as_dict(self) -> dict[str, Any]:
        return {
            "schema": REPORT_TABLE_SCHEMA,
            "evidence_catalog": [_plain(item) for item in self.evidence_catalog],
            "policy_versions": list(self.policy_versions),
            "tables": {
                name: [row.as_dict() for row in getattr(self, name)]
                for name in TABLE_NAMES
            },
        }


def _plain(value: Any) -> Any:
    if is_dataclass(value):
        return {field.name: _plain(getattr(value, field.name)) for field in fields(value)}
    if isinstance(value, Mapping):
        return {str(key): _plain(value[key]) for key in sorted(value)}
    if isinstance(value, (list, tuple)):
        return [_plain(item) for item in value]
    return value


def _evidence_from_statements(
    statements: tuple[EvidenceLinkedStatement, ...],
) -> set[str]:
    return {
        evidence_id
        for statement in statements
        for evidence_id in statement.evidence_ids
    }


def _context_evidence(context: ResolvedReportContext) -> set[str]:
    return (
        set(context.evidence_ids)
        | set(context.exposure.evidence_ids)
        | set(context.safety.evidence_ids)
    )


def _make_row(
    *,
    row_id: str,
    context: ResolvedReportContext,
    details: Mapping[str, Any],
    evidence_ids: set[str],
    evidence_catalog: set[str],
) -> ReportTableRow:
    _require_text(row_id, "report row ID")
    resolved = tuple(sorted(_context_evidence(context) | evidence_ids))
    if not resolved or set(resolved) - evidence_catalog:
        raise ReportTableError(f"report row has unresolved evidence: {row_id}")
    return ReportTableRow(
        row_id=row_id,
        origin=context.origin,
        subtype_applicability=context.subtype_applicability,
        translational_maturity=context.translational_maturity,
        effect_magnitude=context.effect_magnitude,
        child_applicability=context.child_applicability,
        general_mva_applicability=context.general_mva_applicability,
        evidence_ids=resolved,
        policy_version=context.policy_version,
        exposure=context.exposure,
        safety=context.safety,
        decision_reasons=context.decision_reasons,
        details=dict(details),
    )


def _validate_context(
    context: ResolvedReportContext,
    decision_policies: Mapping[str, Mapping[str, object]],
) -> None:
    policy = decision_policies.get(context.policy_version)
    if policy is None or policy.get("policy_version") != context.policy_version:
        raise ReportTableError("report row references an unresolved decision policy")
    maturity = policy.get("evidence_maturity")
    if not isinstance(maturity, list) or context.translational_maturity not in maturity:
        raise ReportTableError("report row uses undefined translational maturity")


def _mechanism_details(readiness: MechanismReadiness) -> dict[str, Any]:
    return {
        "allele_evidence": tuple(
            sorted(readiness.allele_evidence, key=lambda item: item.allele_id)
        ),
        "confidence_inputs": readiness.confidence_inputs,
        "conflicts": tuple(sorted(readiness.conflicts)),
        "direction": readiness.direction,
        "genotype_completeness": readiness.genotype_completeness,
        "inheritance_consistency": readiness.inheritance_consistency,
        "joint_evidence": readiness.joint_evidence,
        "mechanism_confidence": readiness.mechanism_confidence,
        "missing_direction": readiness.missing_direction,
        "phase": readiness.phase,
        "readiness_state": readiness.readiness_state,
        "technical_support": readiness.technical_support,
        "transcript_relevance": readiness.transcript_relevance,
        "prediction_acceptance": readiness.prediction_acceptance,
    }


def _candidate_details(record: CandidateTableRecord) -> dict[str, Any]:
    candidate = record.candidate

    def statement_key(item: EvidenceLinkedStatement) -> tuple[str, tuple[str, ...]]:
        return item.statement, item.evidence_ids

    return {
        "active_ingredient_id": candidate.active_ingredient_id,
        "approval_jurisdictions": tuple(
            sorted(
                record.approval_contexts,
                key=lambda item: (item.jurisdiction, item.state, item.product_ids),
            )
        ),
        "causal_paths": tuple(
            sorted(rationale.causal_path for rationale in candidate.rationales)
        ),
        "conflicts": tuple(sorted(record.conflicts)),
        "current_management": tuple(sorted(record.current_management, key=statement_key)),
        "evidence_dimensions": record.evidence_dimensions,
        "g0_g7": tuple(
            sorted(record.gate_results, key=lambda item: int(item.gate_id[1:]))
        ),
        "generation_rationales": tuple(
            sorted(candidate.rationales, key=lambda item: item.rationale_id)
        ),
        "generation_routes": tuple(sorted(candidate.generation_routes)),
        "independent_rationale_count": candidate.independent_rationale_count,
        "mechanism_alternative_robustness": record.mechanism_alternative_robustness,
        "negative_evidence": tuple(sorted(record.negative_evidence, key=statement_key)),
        "patient_context": record.patient_context,
        "priority_tier": record.tier_decision,
        "proposed_actions": tuple(sorted(record.proposed_actions, key=statement_key)),
        "reachable_polypharmacology": tuple(
            sorted(record.reachable_polypharmacology, key=statement_key)
        ),
        "resolved_drug_identity": {
            "exact_substance_ids": tuple(sorted(candidate.exact_substance_ids)),
            "identity_source_ids": tuple(sorted(candidate.identity_source_ids)),
            "product_ids": tuple(sorted(candidate.product_ids)),
        },
        "route_coverage": tuple(
            sorted(record.route_coverage, key=lambda item: item.coverage_id)
        ),
        "selected_endpoint_id": record.selected_endpoint_id,
        "target_hypothesis_ids": tuple(sorted(record.target_hypothesis_ids)),
        "target_liabilities": tuple(
            sorted(record.target_liabilities, key=statement_key)
        ),
        "target_site_exposure": record.exposure_comparison,
        "chromosome_instability_and_other_safety": record.safety_review,
        "child_and_general_applicability": record.population_applicability,
        "uncertainties": tuple(sorted(record.uncertainties)),
        "research_submission_review": record.research_submission_review,
        "research_eligibility": record.research_eligibility,
    }


def _candidate_evidence(record: CandidateTableRecord) -> set[str]:
    evidence = {
        evidence_id
        for rationale in record.candidate.rationales
        for evidence_id in rationale.evidence_ids
    }
    evidence.update(
        evidence_id
        for approval in record.approval_contexts
        for evidence_id in approval.evidence_ids
    )
    evidence.update(_evidence_from_statements(record.proposed_actions))
    evidence.update(record.mechanism_alternative_robustness.evidence_ids)
    for statements in (
        record.target_liabilities,
        record.reachable_polypharmacology,
        record.current_management,
        record.negative_evidence,
    ):
        evidence.update(_evidence_from_statements(statements))
    if record.exposure_comparison is not None:
        evidence.update({record.exposure_comparison.activity_evidence_id,
                         record.exposure_comparison.exposure_evidence_id})
    evidence.update(record.safety_review.gate_result.evidence_ids)
    evidence.update(record.safety_review.candidate_benefit_evidence_ids)
    evidence.update(record.population_applicability.child_applicability.evidence_ids)
    evidence.update(record.population_applicability.general_mva_applicability.evidence_ids)
    evidence.update(
        dimension.evidence_id
        for dimension in record.patient_context.dimensions
        if dimension.evidence_id is not None
    )
    evidence.update(
        evidence_id
        for gate in record.gate_results
        for evidence_id in gate.evidence_ids
    )
    if record.research_submission_review is not None:
        evidence.update(evidence_id for ids in record.research_submission_review.evidence_by_domain.values()
                        for evidence_id in ids)
    return evidence


def _validate_candidate(
    record: CandidateTableRecord,
    *,
    mechanisms: Mapping[str, MechanismTableRecord],
    endpoints: Mapping[str, EndpointTableRecord],
    targets: Mapping[str, TargetHypothesisTableRecord],
    decision_policies: Mapping[str, Mapping[str, object]],
) -> None:
    candidate_id = record.candidate.active_ingredient_id
    if (
        not candidate_id
        or record.context.origin.scenario_id != record.scenario_id
        or record.scenario_id not in mechanisms
        or record.selected_endpoint_id not in endpoints
        or not record.target_hypothesis_ids
        or set(record.target_hypothesis_ids) - set(targets)
        or not record.candidate.rationales
        or not record.approval_contexts
        or not record.proposed_actions
        or not record.route_coverage
        or set(record.evidence_dimensions) != EVIDENCE_DIMENSIONS
        or any(not isinstance(value, str) or not value for value in record.evidence_dimensions.values())
        or not record.target_liabilities
        or not record.reachable_polypharmacology
        or not record.current_management
        or not record.negative_evidence
        or not record.uncertainties
        or not isinstance(record.display_order, int)
        or isinstance(record.display_order, bool)
        or record.display_order < 0
    ):
        raise ReportTableError("candidate report record is incomplete")
    if (
        record.context.origin != mechanisms[record.scenario_id].context.origin
        or record.context.origin != endpoints[record.selected_endpoint_id].context.origin
    ):
        raise ReportTableError("candidate crosses a source scenario")
    if any(
        target_record.context.origin != record.context.origin
        or target_record.objective.endpoint_id != record.selected_endpoint_id
        for target_id, target_record in targets.items()
        if target_id in record.target_hypothesis_ids
    ):
        raise ReportTableError("candidate target or endpoint crosses a source scenario")
    objective_ids = {
        targets[target_id].objective.objective_id
        for target_id in record.target_hypothesis_ids
    }
    if any(
        rationale.objective_id not in objective_ids
        for rationale in record.candidate.rationales
    ):
        raise ReportTableError("candidate rationale does not resolve to a target objective")
    gates = {gate.gate_id: gate for gate in record.gate_results}
    if len(record.gate_results) != len(GATES) or set(gates) != GATES:
        raise ReportTableError("candidate report requires complete G0-G7 results")
    if gates["G2"].outcome in {"pass", "bounded_uncertainty"} and (
        endpoints[record.selected_endpoint_id].context.subtype_applicability
        != record.context.subtype_applicability
        or any(targets[target_id].context.subtype_applicability
               != record.context.subtype_applicability
               for target_id in record.target_hypothesis_ids)
    ):
        raise ReportTableError("candidate subtype support disagrees with its endpoint or target")
    policy = decision_policies[record.context.policy_version]
    expected_decision = decide_candidate_tier(record.gate_results, policy)
    if expected_decision != record.tier_decision:
        raise ReportTableError("candidate tier does not reproduce the decision policy")
    expected_research = assess_research_eligibility(
        candidate_id, record.gate_results, record.research_submission_review, policy,
    )
    if record.research_eligibility is not None and record.research_eligibility != expected_research:
        raise ReportTableError("research eligibility does not reproduce source and decision inputs")
    if (
        record.tier_decision.policy_version != record.context.policy_version
        or any(gate.policy_version != record.context.policy_version for gate in record.gate_results)
        or gates["G6"] != record.safety_review.gate_result
        or record.context.exposure.state != gates["G5"].outcome
        or record.context.safety.state != record.safety_review.subtype_safety_conclusion
        or record.context.effect_magnitude
        != (record.exposure_comparison.effect_magnitude_state
            if record.exposure_comparison is not None else "magnitude_unknown")
        or record.context.child_applicability
        != record.population_applicability.child_applicability.state
        or record.context.general_mva_applicability
        != record.population_applicability.general_mva_applicability.state
        or (record.exposure_comparison is not None
            and record.exposure_comparison.active_ingredient_id != candidate_id)
        or record.context.decision_reasons
        != (record.tier_decision.decision_reasons or ("pass",))
        or (record.exposure_comparison is not None and not {
            record.exposure_comparison.activity_evidence_id,
            record.exposure_comparison.exposure_evidence_id,
        }.issubset(record.context.exposure.evidence_ids))
        or not set(gates["G6"].evidence_ids).issubset(
            record.context.safety.evidence_ids
        )
    ):
        raise ReportTableError("candidate report context disagrees with structured records")
    expected_dimensions = {
        "exact_subtype": record.context.subtype_applicability,
        "exposure_feasibility": gates["G5"].outcome,
        "mechanism_readiness": mechanisms[
            record.scenario_id
        ].readiness.readiness_state,
        "safety": gates["G6"].outcome,
        "translational_maturity": record.context.translational_maturity,
        "pediatric_practicality": record.population_applicability.child_applicability.state,
    }
    if any(
        record.evidence_dimensions[name] != expected
        for name, expected in expected_dimensions.items()
    ):
        raise ReportTableError("candidate evidence dimensions disagree with source records")
    if (
        {item.route for item in record.route_coverage}
        != set(record.candidate.generation_routes)
        or any(
            set(approval.product_ids) - set(record.candidate.product_ids)
            for approval in record.approval_contexts
        )
    ):
        raise ReportTableError("candidate route coverage or approval identity is unresolved")
    if record.candidate.independent_rationale_count != len(record.candidate.rationales):
        raise ReportTableError("candidate rationale count is inconsistent")
    if record.tier_decision.primary_shortlist_member:
        if record.exclusion is not None:
            raise ReportTableError("primary candidate cannot also be an exclusion")
    else:
        nonpassing = tuple(
            gate_id
            for gate_id in sorted(gates, key=lambda value: int(value[1:]))
            if gates[gate_id].outcome not in {"pass", "bounded_uncertainty"}
        )
        reasons = tuple(
            reason
            for gate_id in sorted(gates, key=lambda value: int(value[1:]))
            for reason in gates[gate_id].reason_codes
            if reason != "pass"
        )
        if (
            record.exclusion is None
            or record.exclusion.active_ingredient_id != candidate_id
            or record.exclusion.tier != record.tier_decision.tier
            or record.exclusion.decisive_gate_ids != nonpassing
            or record.exclusion.reason_codes != reasons
            or not record.exclusion.evidence_gaps
            or not record.exclusion.resolution_experiments
        ):
            raise ReportTableError("non-primary candidate requires a complete exclusion")


def build_report_tables(
    inputs: ReportTableInputs,
    *,
    decision_policies: Mapping[str, Mapping[str, object]],
    allow_unqualified: bool = False,
) -> Track2ReportTables:
    """Build all six tables with strict cross-record referential integrity."""
    if (
        not inputs.evidence_catalog
        or not allow_unqualified and (not inputs.mechanism_records
        or not inputs.endpoint_records
        or not inputs.target_records)
    ):
        raise ReportTableError("all core report table inputs are required")
    catalog = {item.evidence_id: item for item in inputs.evidence_catalog}
    if len(catalog) != len(inputs.evidence_catalog):
        raise ReportTableError("report evidence catalog contains duplicate IDs")
    policy_versions = tuple(sorted(decision_policies))
    if not policy_versions:
        raise ReportTableError("at least one decision policy is required")
    catalog_ids = set(catalog)

    mechanism_by_hypothesis: dict[str, MechanismTableRecord] = {}
    mechanism_rows: list[ReportTableRow] = []
    for record in sorted(
        inputs.mechanism_records, key=lambda item: item.readiness.scenario_id
    ):
        context = record.context
        scenario_id = record.readiness.scenario_id
        _validate_context(context, decision_policies)
        if context.origin.scenario_id != scenario_id or scenario_id in mechanism_by_hypothesis:
            raise ReportTableError("mechanism row has an unresolved or duplicate source origin")
        mechanism_by_hypothesis[scenario_id] = record
        evidence = {
            evidence_id
            for allele in record.readiness.allele_evidence
            for evidence_id in allele.evidence_ids
        }
        mechanism_rows.append(
            _make_row(
                row_id=f"mechanism-readiness:{scenario_id}",
                context=context,
                details=_mechanism_details(record.readiness),
                evidence_ids=evidence,
                evidence_catalog=catalog_ids,
            )
        )

    endpoint_by_id: dict[str, EndpointTableRecord] = {}
    endpoint_rows: list[ReportTableRow] = []
    for record in sorted(inputs.endpoint_records, key=lambda item: item.endpoint.endpoint_id):
        context = record.context
        endpoint_id = record.endpoint.endpoint_id
        _validate_context(context, decision_policies)
        if (
            record.scenario_id not in mechanism_by_hypothesis
            or context.origin.scenario_id != record.scenario_id
            or context.origin
            != mechanism_by_hypothesis[record.scenario_id].context.origin
            or endpoint_id in endpoint_by_id
        ):
            raise ReportTableError("endpoint row has an unresolved hypothesis or duplicate ID")
        endpoint_by_id[endpoint_id] = record
        endpoint_rows.append(
            _make_row(
                row_id=endpoint_id,
                context=context,
                details={
                    **record.endpoint.as_track2_payload(),
                    "relevance": record.endpoint.relevance,
                    "support": tuple(
                        sorted(record.endpoint.support, key=lambda item: item.evidence_id)
                    ),
                },
                evidence_ids={item.evidence_id for item in record.endpoint.support},
                evidence_catalog=catalog_ids,
            )
        )

    target_by_id: dict[str, TargetHypothesisTableRecord] = {}
    target_rows: list[ReportTableRow] = []
    for record in sorted(inputs.target_records, key=lambda item: item.target.target_id):
        context = record.context
        target_id = record.target.target_id
        _validate_context(context, decision_policies)
        if (
            context.origin.scenario_id != record.objective.mechanism_scenario_id
            or record.objective.mechanism_scenario_id not in mechanism_by_hypothesis
            or record.objective.endpoint_id not in endpoint_by_id
            or context.origin
            != mechanism_by_hypothesis[
                record.objective.mechanism_scenario_id
            ].context.origin
            or context.origin
            != endpoint_by_id[record.objective.endpoint_id].context.origin
            or record.target.objective_id != record.objective.objective_id
            or target_id in target_by_id
        ):
            raise ReportTableError("target row has an unresolved mechanism, endpoint, or objective")
        target_by_id[target_id] = record
        evidence = set(record.objective.evidence_ids) | set(record.target.evidence_ids)
        evidence.update(
            evidence_id
            for edge in record.target.signed_rescue_path
            for evidence_id in edge.evidence_ids
        )
        target_rows.append(
            _make_row(
                row_id=target_id,
                context=context,
                details={
                    "intervention_objective": record.objective,
                    **{
                        **record.target.as_track2_payload(),
                        "alternatives": sorted(record.target.alternatives),
                        "conflicts": sorted(record.target.conflicts),
                        "harm_conditions": sorted(record.target.harm_conditions),
                    },
                    "ready_for_drug_retrieval": record.target.ready_for_drug_retrieval,
                },
                evidence_ids=evidence,
                evidence_catalog=catalog_ids,
            )
        )

    candidate_rows: list[ReportTableRow] = []
    exclusion_rows: list[ReportTableRow] = []
    gate_rows: list[tuple[tuple[int, str, int], ReportTableRow]] = []
    seen_candidates: set[str] = set()
    for record in sorted(
        inputs.candidate_records,
        key=lambda item: (item.display_order, item.candidate.active_ingredient_id),
    ):
        context = record.context
        candidate_id = record.candidate.active_ingredient_id
        _validate_context(context, decision_policies)
        if candidate_id in seen_candidates:
            raise ReportTableError("candidate table contains duplicate active ingredients")
        seen_candidates.add(candidate_id)
        _validate_candidate(
            record,
            mechanisms=mechanism_by_hypothesis,
            endpoints=endpoint_by_id,
            targets=target_by_id,
            decision_policies=decision_policies,
        )
        evidence = _candidate_evidence(record)
        mechanism = mechanism_by_hypothesis[record.scenario_id].readiness
        evidence.update(
            evidence_id
            for allele in mechanism.allele_evidence
            for evidence_id in allele.evidence_ids
        )
        evidence.update(
            item.evidence_id
            for item in endpoint_by_id[record.selected_endpoint_id].endpoint.support
        )
        for target_id in record.target_hypothesis_ids:
            target_record = target_by_id[target_id]
            evidence.update(target_record.objective.evidence_ids)
            evidence.update(target_record.target.evidence_ids)
            evidence.update(
                evidence_id
                for edge in target_record.target.signed_rescue_path
                for evidence_id in edge.evidence_ids
            )
        details = _candidate_details(record)
        details["research_eligibility"] = assess_research_eligibility(
            candidate_id, record.gate_results, record.research_submission_review,
            decision_policies[context.policy_version],
        )
        candidate_rows.append(
            _make_row(
                row_id=candidate_id,
                context=context,
                details=details,
                evidence_ids=evidence,
                evidence_catalog=catalog_ids,
            )
        )
        if record.exclusion is not None:
            exclusion_rows.append(
                _make_row(
                    row_id=f"candidate-exclusion:{candidate_id}",
                    context=context,
                    details={**details, "exclusion": record.exclusion},
                    evidence_ids=evidence,
                    evidence_catalog=catalog_ids,
                )
            )
        for gate in record.gate_results:
            gate_rows.append(
                (
                    (record.display_order, candidate_id, int(gate.gate_id[1:])),
                    _make_row(
                        row_id=f"decision-gate:{candidate_id}:{gate.gate_id}",
                        context=context,
                        details={
                            "active_ingredient_id": candidate_id,
                            "coverage_state": gate.coverage_state,
                            "gate_id": gate.gate_id,
                            "gate_outcome": gate.outcome,
                            "gate_reason_codes": gate.reason_codes,
                            "gate_review_state": gate.review_state,
                            "priority_tier": record.tier_decision.tier,
                        },
                        evidence_ids=set(gate.evidence_ids) | evidence,
                        evidence_catalog=catalog_ids,
                    ),
                )
            )

    return Track2ReportTables(
        evidence_catalog=tuple(sorted(catalog.values(), key=lambda item: item.evidence_id)),
        policy_versions=policy_versions,
        mechanism_readiness=tuple(mechanism_rows),
        endpoints=tuple(endpoint_rows),
        target_hypotheses=tuple(target_rows),
        candidates=tuple(candidate_rows),
        exclusions=tuple(exclusion_rows),
        gates=tuple(row for _, row in sorted(gate_rows, key=lambda item: item[0])),
    )


def report_tables_as_json(tables: Track2ReportTables) -> str:
    """Serialize report tables deterministically for package checksums."""
    return json.dumps(tables.as_dict(), indent=2, sort_keys=True, ensure_ascii=True) + "\n"


def replay_research_membership(
    row: ReportTableRow, decision_policy: Mapping[str, object],
) -> ResearchEligibilityDecision:
    expected = assess_research_eligibility(
        row.row_id, row.details["g0_g7"], row.details.get("research_submission_review"), decision_policy,
    )
    if row.details.get("research_eligibility") != expected:
        raise ReportTableError("research eligibility does not reproduce source and decision inputs")
    if set(expected.evidence_ids) - set(row.evidence_ids):
        raise ReportTableError("research evidence does not resolve to the candidate row")
    return expected


def replay_serialized_research_membership(
    row: Mapping[str, Any], decision_policy: Mapping[str, object],
) -> ResearchEligibilityDecision:
    """Recheck exported tables before public projection or release preparation."""
    if row.get('details', {}).get('record_kind') == 'discovery_review':
        raise ReportTableError('discovery review rows have no scientific membership assessment')
    try:
        details = row["details"]
        gates = tuple(GateResult(**{
            **gate, "reason_codes": tuple(gate["reason_codes"]),
            "evidence_ids": tuple(gate["evidence_ids"]),
        }) for gate in details["g0_g7"])
        raw_review = details.get("research_submission_review")
        review = None if raw_review is None else ResearchSubmissionReview(**{
            **raw_review,
            "evidence_by_domain": {domain: tuple(ids) for domain, ids in raw_review["evidence_by_domain"].items()},
            "gaps": tuple(ResearchGap(**gap) for gap in raw_review["gaps"]),
        })
        expected = assess_research_eligibility(row["row_id"], gates, review, decision_policy)
        stored = details["research_eligibility"]
        if (type(stored["eligible"]) is not bool
                or stored != _plain(expected)
                or details["priority_tier"] != _plain(decide_candidate_tier(gates, decision_policy))
                or set(expected.evidence_ids) - set(row["evidence_ids"])):
            raise ReportTableError("serialized research eligibility failed source and decision replay")
        return expected
    except (KeyError, TypeError, ValueError, AttributeError) as error:
        raise ReportTableError("serialized research eligibility failed source and decision replay") from error


def _markdown_value(value: Any) -> str:
    plain = _plain(value)
    if isinstance(plain, str):
        rendered = plain
    elif plain is None:
        rendered = "null"
    elif isinstance(plain, (dict, list)):
        rendered = json.dumps(plain, sort_keys=True, separators=(",", ":"))
    else:
        rendered = str(plain).lower() if isinstance(plain, bool) else str(plain)
    if isinstance(plain, (dict, list)):
        # Preserve source records as literal data, distinct from report prose.
        fence = '`' * (max((len(m.group()) for m in re.finditer(r'`+', rendered)), default=0) + 1)
        rendered = fence + ' ' + rendered + ' ' + fence
    return rendered.replace("\\", "\\\\").replace("|", "\\|").replace("\r\n", "\n").replace("\n", "<br>")


def _render_markdown_table(rows: tuple[ReportTableRow, ...]) -> list[str]:
    if not rows:
        return ["_No rows._"]
    detail_keys = tuple(sorted({key for row in rows for key in row.details}))
    common_headers = (
        "Row ID",
        "Scenario",
        "Input Origin",
        "Source SHA-256",
        "Subtype Applicability",
        "Translational Maturity",
        "Effect Magnitude",
        "Child Applicability",
        "General MVA Applicability",
        "Evidence IDs",
        "Policy Version",
        "Exposure",
        "Safety",
        "Decision Reasons",
    )
    detail_headers = tuple(key.replace("_", " ").title() for key in detail_keys)
    headers = common_headers + detail_headers
    rendered = [
        "| " + " | ".join(headers) + " |",
        "|" + "|".join("---" for _ in headers) + "|",
    ]
    for row in rows:
        values: tuple[Any, ...] = (
            row.row_id,
            row.origin.scenario_id,
            row.origin.input_origin,
            row.origin.source_sha256,
            row.subtype_applicability,
            row.translational_maturity,
            row.effect_magnitude,
            row.child_applicability,
            row.general_mva_applicability,
            row.evidence_ids,
            row.policy_version,
            row.exposure,
            row.safety,
            row.decision_reasons,
        ) + tuple(row.details.get(key) for key in detail_keys)
        rendered.append("| " + " | ".join(_markdown_value(value) for value in values) + " |")
    return rendered


def _render_discovery_table(rows, name):
    from genosapiens.track2_lincs_priority import GROUPS, md
    if name == 'gates':
        grouped = {}
        for row in rows:
            grouped.setdefault(row.details['active_ingredient_id'], {})[row.details['gate_id']] = row
        lines = ['| Medication | G0 | G1 | G2 | G3 | G4 | G5 | G6 | G7 |',
                 '|---|' + '---|' * 8]
        for gates in grouped.values():
            first = next(iter(gates.values()))
            values = [gates[g].details['gate_outcome'].replace('_', ' ') if g in gates else 'not assessed'
                      for g in (f'G{i}' for i in range(8))]
            lines.append('| ' + md(first.details['name']) + ' | ' + ' | '.join(values) + ' |')
        return lines
    lines = ['| Medication | Review rank | Local agreement | Research status | Reason / next step |',
             '|---|---:|---|---|---|']
    for row in rows:
        d = row.details
        priority = d['review_priority']
        color, group = GROUPS[priority['agreement_group']]
        reason = d.get('exclusion_summary',
            ('Eligible by repeated LINCS agreement alone. ' if 'lincs_repeated_consistent' in d.get('eligibility_bases', []) else '')
            + 'Independent drug action, rescue, exposure, safety and validation review remain separate.')
        # Reformat older retained summaries without changing their source records.
        reason = reason.replace('Expression hit; independent drug review still missing.',
                                'The expression checks found a match. Independent medication review remains missing.')
        reason = reason.replace('No assessed cell line passed both directional checks; no LINCS nomination.',
                                'No assessed cell line passed both directional checks. There is no LINCS nomination.')
        lines.append(f"| {md(d['name'])} | {priority['review_rank']} | {color} {group} | "
                     + ('Eligible' if d['research_eligible'] else 'Not eligible on current evidence')
                     + f" | {md(reason)} |")
    return lines


def report_tables_as_markdown(tables: Track2ReportTables, *, heading_level=1) -> str:
    """Render the six human-readable tables from the same immutable rows."""
    titles = {
        "mechanism_readiness": "Mechanism Readiness",
        "endpoints": "Endpoints",
        "target_hypotheses": "Target Hypotheses",
        "candidates": "Candidates",
        "exclusions": "Exclusions",
        "gates": "G0-G7 Gates",
    }
    lines = ["#" * heading_level + " Track 2 Structured Evidence Tables", ""]
    for name in TABLE_NAMES:
        lines.extend(("#" * (heading_level + 1) + f" {titles[name]}", ""))
        rows = getattr(tables, name)
        assessed = tuple(r for r in rows if r.details.get('record_kind') != 'discovery_review')
        discovery = tuple(r for r in rows if r.details.get('record_kind') == 'discovery_review')
        if assessed or not discovery:
            lines.extend(_render_markdown_table(assessed))
        if discovery:
            from genosapiens.track2_lincs_priority import CAVEAT
            lines += ['', '**Cell-specific LINCS**', '', CAVEAT, '']
            if name == 'gates':
                lines += ['G0 records public input that passed the input checks. G1–G7 remain unassessed for these medications. '
                          'An FDA snapshot link alone does not pass G1. Full source bindings and reasons are in tables.json.', '']
            lines.extend(_render_discovery_table(discovery, name))
        lines.append("")
    return "\n".join(lines)
