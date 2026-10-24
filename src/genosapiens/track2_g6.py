"""MVA-specific G6 chromosome-instability and cumulative safety review."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Iterable, Mapping

from genosapiens.track2_gates import COVERAGE_STATES, GateResult


MVA_HAZARD_DOMAINS = {
    "aneugenicity",
    "carcinogenicity",
    "chromosome_missegregation",
    "clastogenicity",
    "cumulative_exposure",
    "gonadal_toxicity",
    "immune_reserve",
    "marrow_reserve",
    "micronuclei",
    "mutagenicity",
    "secondary_malignancy",
}
FINDINGS = {"negative_test", "no_harm_observed", "positive_hazard", "unknown"}
APPLICABILITY = {
    "exact_subtype",
    "external_nontransferable",
    "general_human",
    "transferable",
    "unresolved",
}
EXPOSURE_RELEVANCE = {"not_supported", "supported", "unknown"}
HAZARD_DECISIONS = {"blocking", "nonblocking", "unresolved"}
POPULATION_STATES = {"conditional", "not_assessed", "not_evaluable", "not_supported", "supported"}


class G6EvaluationError(ValueError):
    """Raised when a safety claim omits an MVA hazard or overstates reassurance."""


@dataclass(frozen=True, slots=True)
class MVAHazardAssessment:
    domain: str
    finding: str
    applicability: str
    exposure_relevance: str
    hazard_decision: str
    rationale: str
    evidence_ids: tuple[str, ...]

    def __post_init__(self) -> None:
        if (
            self.domain not in MVA_HAZARD_DOMAINS
            or self.finding not in FINDINGS
            or self.applicability not in APPLICABILITY
            or self.exposure_relevance not in EXPOSURE_RELEVANCE
            or self.hazard_decision not in HAZARD_DECISIONS
            or not self.rationale
            or not self.evidence_ids
        ):
            raise G6EvaluationError("MVA hazard assessment is incomplete")
        if self.hazard_decision == "blocking" and not (
            self.finding == "positive_hazard"
            and self.applicability in {"exact_subtype", "general_human", "transferable"}
            and self.exposure_relevance == "supported"
        ):
            raise G6EvaluationError("blocking hazard is not transferable at relevant exposure")
        if (
            self.finding in {"negative_test", "no_harm_observed"}
            and self.applicability == "external_nontransferable"
            and self.hazard_decision != "unresolved"
        ):
            raise G6EvaluationError(
                "external absence of harm cannot establish selected-subtype safety"
            )


@dataclass(frozen=True, slots=True)
class G6SafetyReview:
    gate_result: GateResult
    blocking_hazard_domains: tuple[str, ...]
    cancer_surveillance_requirements: tuple[str, ...]
    cancer_prevention_claim_supported: bool
    subtype_safety_conclusion: str
    candidate_benefit_evidence_ids: tuple[str, ...]


def evaluate_mva_g6(
    assessments: Iterable[MVAHazardAssessment],
    *,
    child_applicability_state: str,
    general_mva_applicability_state: str,
    cancer_surveillance_requirements: tuple[str, ...],
    claimed_cancer_prevention_from_proximal_biomarker: bool,
    candidate_benefit_evidence_ids: tuple[str, ...],
    coverage_state: str,
    review_state: str,
    decision_policy: Mapping[str, object],
) -> G6SafetyReview:
    records = tuple(sorted(assessments, key=lambda item: item.domain))
    if (
        len(records) != len(MVA_HAZARD_DOMAINS)
        or {item.domain for item in records} != MVA_HAZARD_DOMAINS
    ):
        raise G6EvaluationError("all MVA hazard domains require explicit review")
    if (
        child_applicability_state not in POPULATION_STATES
        or general_mva_applicability_state not in POPULATION_STATES
        or not cancer_surveillance_requirements
        or not candidate_benefit_evidence_ids
        or coverage_state not in COVERAGE_STATES
    ):
        raise G6EvaluationError("G6 population, surveillance, or audit context is incomplete")
    if claimed_cancer_prevention_from_proximal_biomarker:
        raise G6EvaluationError(
            "proximal chromosome-instability biomarkers cannot establish cancer prevention"
        )
    review = decision_policy.get("review")
    allowed_reviews = set(review.get("independence_states", ())) if isinstance(review, Mapping) else set()
    if review_state not in allowed_reviews:
        raise G6EvaluationError("G6 review state is undefined")
    blocking = tuple(item for item in records if item.hazard_decision == "blocking")
    external_absence = any(
        item.finding in {"negative_test", "no_harm_observed"}
        and item.applicability == "external_nontransferable"
        for item in records
    )
    unresolved = any(
        item.finding == "unknown"
        or item.hazard_decision == "unresolved"
        or item.applicability == "unresolved"
        or (
            item.finding == "positive_hazard"
            and item.exposure_relevance == "unknown"
        )
        for item in records
    )
    population_incomplete = (
        child_applicability_state in {"not_assessed", "not_evaluable", "not_supported"}
        or general_mva_applicability_state in {"not_assessed", "not_evaluable", "not_supported"}
    )
    if blocking:
        outcome, reasons = "contraindicating", ("blocking_hazard",)
        conclusion = "blocking_hazard_identified"
    elif unresolved or external_absence or population_incomplete:
        outcome, reasons = "not_evaluable", ("safety_evidence_incomplete",)
        conclusion = "not_established"
    else:
        outcome, reasons = "pass", ("pass",)
        conclusion = "no_blocking_hazard_identified_not_proof_of_safety"
    evidence_ids = tuple(
        sorted(
            set(candidate_benefit_evidence_ids)
            | {evidence for item in records for evidence in item.evidence_ids}
        )
    )
    gate = GateResult(
        "G6",
        outcome,
        reasons,
        evidence_ids,
        coverage_state,
        decision_policy["policy_version"],
        review_state,
    )
    return G6SafetyReview(
        gate_result=gate,
        blocking_hazard_domains=tuple(item.domain for item in blocking),
        cancer_surveillance_requirements=cancer_surveillance_requirements,
        cancer_prevention_claim_supported=False,
        subtype_safety_conclusion=conclusion,
        candidate_benefit_evidence_ids=candidate_benefit_evidence_ids,
    )
