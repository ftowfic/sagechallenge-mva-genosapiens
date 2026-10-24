"""Mechanism-linked endpoint selection with surrogate-claim boundaries."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Iterable


ENDPOINT_CLASSES = {
    "cancer_risk_outcome",
    "fixed_developmental_injury",
    "ongoing_pathology",
    "proximal_rescue",
    "subtype_functional_rescue",
    "symptomatic_complication",
    "target_engagement",
    "toxicity",
}
SUPPORT_ROLES = {"modifiability", "relevance"}
EVIDENCE_TYPES = {
    "clinical_outcome",
    "fixed_injury_evidence",
    "intervention_response",
    "longitudinal_reversibility",
    "mechanistic_evidence",
    "natural_history",
    "track1_phenotype_contribution",
}
MODIFIABILITY_TYPES = {
    "clinical_outcome",
    "intervention_response",
    "longitudinal_reversibility",
}
CELLULAR_CANCER_SURROGATES = {
    "aneuploid_cell_fraction",
    "micronuclei",
    "segregation_error_rate",
}


class EndpointSelectionError(ValueError):
    """Raised when an endpoint overstates relevance, modifiability, or meaning."""


@dataclass(frozen=True, slots=True)
class EndpointSupport:
    evidence_id: str
    role: str
    evidence_type: str

    def __post_init__(self) -> None:
        if (
            not self.evidence_id
            or self.role not in SUPPORT_ROLES
            or self.evidence_type not in EVIDENCE_TYPES
        ):
            raise EndpointSelectionError("endpoint support is incomplete")
        if (
            self.role == "modifiability"
            and self.evidence_type == "track1_phenotype_contribution"
        ):
            raise EndpointSelectionError(
                "Track 1 phenotype contribution establishes relevance, not modifiability"
            )


@dataclass(frozen=True, slots=True)
class EndpointSelection:
    endpoint_id: str
    endpoint_class: str
    measurement: str
    expected_direction: str
    causal_proximity: str
    tissue: str
    stage: str
    meaningful_interpretation: str
    support: tuple[EndpointSupport, ...]
    relevance: str
    modifiability: str

    def as_track2_payload(self) -> dict:
        return {
            "causal_proximity": self.causal_proximity,
            "endpoint_class": self.endpoint_class,
            "expected_direction": self.expected_direction,
            "meaningful_interpretation": self.meaningful_interpretation,
            "measurement": self.measurement,
            "modifiability": self.modifiability,
            "stage": self.stage,
            "tissue": self.tissue,
        }


def select_endpoint(
    *,
    endpoint_id: str,
    endpoint_class: str,
    measurement: str,
    expected_direction: str,
    causal_proximity: str,
    tissue: str,
    stage: str,
    meaningful_interpretation: str,
    support: Iterable[EndpointSupport],
) -> EndpointSelection:
    support_records = tuple(support)
    if (
        not endpoint_id
        or endpoint_class not in ENDPOINT_CLASSES
        or not measurement
        or expected_direction not in {"decrease", "increase", "maintain", "unknown"}
        or causal_proximity
        not in {"clinical", "distal", "proximal", "target", "unknown"}
        or not tissue
        or not stage
        or not meaningful_interpretation
        or not support_records
        or len({item.evidence_id for item in support_records}) != len(support_records)
    ):
        raise EndpointSelectionError("endpoint selection is incomplete")
    normalized_measurement = measurement.lower().replace(" ", "_")
    interpretation = meaningful_interpretation.lower().replace("-", " ")
    if normalized_measurement in CELLULAR_CANCER_SURROGATES and (
        endpoint_class == "cancer_risk_outcome"
        or "cancer incidence" in interpretation
        or "tumor incidence" in interpretation
    ):
        raise EndpointSelectionError(
            "segregation, micronuclei, or aneuploid-cell endpoints cannot establish reduced cancer incidence"
        )
    if endpoint_class == "cancer_risk_outcome" and not any(
        item.evidence_type == "clinical_outcome" and item.role == "relevance"
        for item in support_records
    ):
        raise EndpointSelectionError(
            "cancer-risk outcomes require direct clinical outcome relevance"
        )
    relevance = (
        "supported"
        if any(item.role == "relevance" for item in support_records)
        else "unknown"
    )
    modifiability_records = [
        item for item in support_records if item.role == "modifiability"
    ]
    if any(item.evidence_type in MODIFIABILITY_TYPES for item in modifiability_records):
        modifiability = "supported"
    elif endpoint_class == "fixed_developmental_injury" or any(
        item.evidence_type == "fixed_injury_evidence" for item in support_records
    ):
        modifiability = "not_supported"
    else:
        modifiability = "unknown"
    return EndpointSelection(
        endpoint_id=endpoint_id,
        endpoint_class=endpoint_class,
        measurement=measurement,
        expected_direction=expected_direction,
        causal_proximity=causal_proximity,
        tissue=tissue,
        stage=stage,
        meaningful_interpretation=meaningful_interpretation,
        support=support_records,
        relevance=relevance,
        modifiability=modifiability,
    )
