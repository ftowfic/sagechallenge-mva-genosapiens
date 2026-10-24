"""Context-preserving directional pharmacology records and synthesis."""

from __future__ import annotations

from dataclasses import dataclass
from itertools import combinations
from typing import Iterable


SOURCES = {"chembl", "gtopdb", "primary-study", "regulatory"}
ASSAY_TYPES = {"binding", "clinical", "functional", "regulatory"}
DIRECTIONS = {"decrease", "increase", "neutral", "unknown"}
QUALIFIERS = {"<", "<=", "=", ">", ">=", "approximately", "not_reported"}
CONFIDENCE = {"high", "low", "moderate", "not_reported"}


class PharmacologyError(ValueError):
    """Raised when an assay record loses material pharmacology context."""


@dataclass(frozen=True, slots=True)
class ConcentrationResponseRange:
    minimum: float
    maximum: float
    unit: str

    def __post_init__(self) -> None:
        if self.minimum < 0 or self.maximum < self.minimum or not self.unit:
            raise PharmacologyError("concentration-response range is invalid")


@dataclass(frozen=True, slots=True)
class PharmacologyRecord:
    evidence_id: str
    source: str
    drug_id: str
    action: str
    result_direction: str
    target_id: str
    assay_type: str
    assay: str
    endpoint: str
    value: float | None
    unit: str | None
    qualifier: str
    species: str
    system: str
    concentration_response_range: ConcentrationResponseRange | None
    maximum_effect_value: float | None
    maximum_effect_unit: str | None
    partial_modulation: bool | None
    selectivity: str
    uncertainty: str
    confidence: str
    document_lineage: tuple[str, ...]

    def __post_init__(self) -> None:
        if (
            not self.evidence_id
            or self.source not in SOURCES
            or not self.drug_id
            or not self.action
            or self.result_direction not in DIRECTIONS
            or not self.target_id
            or self.assay_type not in ASSAY_TYPES
            or not self.assay
            or not self.endpoint
            or self.qualifier not in QUALIFIERS
            or not self.species
            or not self.system
            or not self.selectivity
            or not self.uncertainty
            or self.confidence not in CONFIDENCE
            or not self.document_lineage
            or len(self.document_lineage) != len(set(self.document_lineage))
        ):
            raise PharmacologyError("pharmacology record is incomplete")
        if (self.value is None) != (self.unit is None):
            raise PharmacologyError("pharmacology value and unit must appear together")
        if (self.maximum_effect_value is None) != (self.maximum_effect_unit is None):
            raise PharmacologyError("maximum effect value and unit must appear together")
        if self.value is None and self.qualifier != "not_reported":
            raise PharmacologyError("missing values require a not_reported qualifier")
        if self.value is not None and self.qualifier == "not_reported":
            raise PharmacologyError("reported values require a measurement qualifier")
        if self.assay_type == "binding" and (
            self.partial_modulation is not None or self.maximum_effect_value is not None
        ):
            raise PharmacologyError(
                "binding-only records cannot claim functional modulation or maximum effect"
            )

    @property
    def interpretation(self) -> str:
        if self.assay_type == "binding":
            return "binding_only"
        if self.value is None and self.maximum_effect_value is None:
            return "magnitude_unknown"
        if self.qualifier in {"<", "<=", ">", ">="}:
            return "censored"
        return "directional_measurement"

    @property
    def comparability_key(self) -> tuple[str, ...]:
        return (
            self.drug_id,
            self.target_id,
            self.assay_type,
            self.assay,
            self.endpoint,
            self.unit or "not_reported",
            self.species,
            self.system,
            self.maximum_effect_unit or "not_reported",
        )


@dataclass(frozen=True, slots=True)
class PharmacologySynthesis:
    records: tuple[PharmacologyRecord, ...]
    binding_only_ids: tuple[str, ...]
    magnitude_unknown_ids: tuple[str, ...]
    censored_ids: tuple[str, ...]
    contradiction_groups: tuple[tuple[str, ...], ...]
    incomparable_pairs: tuple[tuple[str, str], ...]


def synthesize_pharmacology(
    records: Iterable[PharmacologyRecord],
) -> PharmacologySynthesis:
    ordered = tuple(sorted(records, key=lambda item: item.evidence_id))
    if len(ordered) != len({item.evidence_id for item in ordered}):
        raise PharmacologyError("pharmacology evidence IDs must be unique")
    by_context: dict[tuple[str, ...], list[PharmacologyRecord]] = {}
    for record in ordered:
        by_context.setdefault(record.comparability_key, []).append(record)
    contradictions = tuple(
        tuple(item.evidence_id for item in context_records)
        for _key, context_records in sorted(by_context.items())
        if len(
            {
                item.result_direction
                for item in context_records
                if item.result_direction != "unknown"
            }
        )
        > 1
    )
    incomparable = tuple(
        (left.evidence_id, right.evidence_id)
        for left, right in combinations(ordered, 2)
        if left.drug_id == right.drug_id
        and left.target_id == right.target_id
        and left.comparability_key != right.comparability_key
    )
    return PharmacologySynthesis(
        records=ordered,
        binding_only_ids=tuple(
            item.evidence_id for item in ordered if item.interpretation == "binding_only"
        ),
        magnitude_unknown_ids=tuple(
            item.evidence_id
            for item in ordered
            if item.interpretation == "magnitude_unknown"
        ),
        censored_ids=tuple(
            item.evidence_id for item in ordered if item.interpretation == "censored"
        ),
        contradiction_groups=contradictions,
        incomparable_pairs=incomparable,
    )
