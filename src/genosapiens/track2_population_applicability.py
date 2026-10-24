"""Separate child-specific and general-MVA applicability conclusions."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Iterable


DIMENSIONS = {
    "age_context",
    "approved_regimen",
    "developmental_safety",
    "exposure",
    "formulation",
    "organ_function",
    "preferences",
    "treatment_history",
}
STATES = {"concern", "not_applicable", "supported", "unknown"}


class PopulationApplicabilityError(ValueError):
    """Raised when child or general applicability silently omits a material factor."""


@dataclass(frozen=True, slots=True)
class ApplicabilityDimension:
    name: str
    state: str
    material: bool
    evidence_ids: tuple[str, ...]
    rationale: str

    def __post_init__(self) -> None:
        if (
            self.name not in DIMENSIONS
            or self.state not in STATES
            or not self.rationale
            or (
                self.state in {"supported", "concern"}
                and not self.evidence_ids
            )
            or (
                self.state in {"unknown", "not_applicable"}
                and self.evidence_ids
            )
        ):
            raise PopulationApplicabilityError("applicability dimension is incomplete")


@dataclass(frozen=True, slots=True)
class ApplicabilityConclusion:
    population: str
    state: str
    dimensions: tuple[ApplicabilityDimension, ...]
    limiting_dimensions: tuple[str, ...]
    evidence_ids: tuple[str, ...]


@dataclass(frozen=True, slots=True)
class PopulationApplicability:
    child_applicability: ApplicabilityConclusion
    general_mva_applicability: ApplicabilityConclusion


def _conclude(
    population: str, dimensions: Iterable[ApplicabilityDimension]
) -> ApplicabilityConclusion:
    records = tuple(sorted(dimensions, key=lambda item: item.name))
    if {item.name for item in records} != DIMENSIONS or len(records) != len(DIMENSIONS):
        raise PopulationApplicabilityError(
            f"{population} applicability omits a required material dimension"
        )
    material = tuple(item for item in records if item.material)
    concerns = tuple(item.name for item in material if item.state == "concern")
    unknowns = tuple(item.name for item in material if item.state == "unknown")
    if concerns:
        state = "not_supported"
        limiting = concerns
    elif unknowns:
        # Evidence in other dimensions cannot bound absent exposure, age or
        # developmental-safety evidence for this population.
        state = "not_evaluable"
        limiting = unknowns
    else:
        state = "supported"
        limiting = ()
    return ApplicabilityConclusion(
        population=population,
        state=state,
        dimensions=records,
        limiting_dimensions=limiting,
        evidence_ids=tuple(
            sorted({evidence_id for item in records for evidence_id in item.evidence_ids})
        ),
    )


def assess_population_applicability(
    child_dimensions: Iterable[ApplicabilityDimension],
    general_mva_dimensions: Iterable[ApplicabilityDimension],
) -> PopulationApplicability:
    """Compute each population independently; neither conclusion can alter the other."""
    return PopulationApplicability(
        child_applicability=_conclude("child", child_dimensions),
        general_mva_applicability=_conclude("general_MVA", general_mva_dimensions),
    )
