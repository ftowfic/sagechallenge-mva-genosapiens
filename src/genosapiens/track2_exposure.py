"""Context-preserving normalization of pharmacology and human exposure."""

from __future__ import annotations

from dataclasses import dataclass


MOLAR_TO_NM = {"M": 1_000_000_000.0, "mM": 1_000_000.0, "uM": 1_000.0, "nM": 1.0}
MASS_G_PER_L = {"ug/mL": 1e-3, "mg/L": 1e-3, "ng/mL": 1e-6, "ug/L": 1e-6}
CONCENTRATION_UNITS = set(MOLAR_TO_NM) | set(MASS_G_PER_L)
QUALIFIERS = {"<", "<=", "=", ">", ">=", "approximately"}
CONCENTRATION_BASES = {"free", "nominal", "total"}
ANALYTE_TYPES = {"active_metabolite", "parent"}
PK_METRICS = {"cmax", "steady_state", "trough"}
LOCATIONS = {"extracellular", "intracellular", "plasma", "tissue"}
DIRECTIONS = {"decrease", "increase"}


class ExposureNormalizationError(ValueError):
    """Raised when incompatible potency, exposure, or effect values are combined."""


@dataclass(frozen=True, slots=True)
class ConcentrationMeasure:
    value: float
    unit: str
    qualifier: str
    molecular_weight_g_mol: float | None

    def __post_init__(self) -> None:
        if (
            self.value < 0
            or self.unit not in CONCENTRATION_UNITS
            or self.qualifier not in QUALIFIERS
            or (
                self.unit in MASS_G_PER_L
                and (
                    self.molecular_weight_g_mol is None
                    or self.molecular_weight_g_mol <= 0
                )
            )
        ):
            raise ExposureNormalizationError("concentration measure is not normalizable")

    @property
    def normalized_nm(self) -> float:
        if self.unit in MOLAR_TO_NM:
            return self.value * MOLAR_TO_NM[self.unit]
        return (
            self.value
            * MASS_G_PER_L[self.unit]
            / self.molecular_weight_g_mol
            * 1_000_000_000.0
        )


@dataclass(frozen=True, slots=True)
class EffectMeasure:
    value: float
    unit: str

    def __post_init__(self) -> None:
        if self.value < 0 or not self.unit:
            raise ExposureNormalizationError("effect measure is invalid")


@dataclass(frozen=True, slots=True)
class ActivityEvidence:
    evidence_id: str
    active_ingredient_id: str
    analyte_id: str
    analyte_type: str
    target_id: str
    potency_metric: str
    potency: ConcentrationMeasure
    concentration_response_range: tuple[ConcentrationMeasure, ConcentrationMeasure]
    maximum_effect: EffectMeasure | None
    partial_modulation: bool | None
    minimum_endpoint_relevant_modulation: EffectMeasure | None
    result_direction: str
    required_endpoint_direction: str
    concentration_basis: str
    tissue: str
    location: str
    target_occupancy: EffectMeasure | None

    def __post_init__(self) -> None:
        lower, upper = self.concentration_response_range
        if (
            not self.evidence_id
            or not self.active_ingredient_id
            or not self.analyte_id
            or self.analyte_type not in ANALYTE_TYPES
            or not self.target_id
            or not self.potency_metric
            or lower.normalized_nm > upper.normalized_nm
            or self.result_direction not in DIRECTIONS
            or self.required_endpoint_direction not in DIRECTIONS
            or self.concentration_basis not in CONCENTRATION_BASES
            or not self.tissue
            or self.location not in LOCATIONS
        ):
            raise ExposureNormalizationError("activity evidence is incomplete or inconsistent")


@dataclass(frozen=True, slots=True)
class ExposureEvidence:
    evidence_id: str
    active_ingredient_id: str
    analyte_id: str
    analyte_type: str
    concentration: ConcentrationMeasure
    concentration_basis: str
    pk_metric: str
    time_above_threshold: str
    schedule: str
    duration: str
    tissue: str
    location: str
    active_metabolite_ids: tuple[str, ...]
    target_occupancy: EffectMeasure | None

    def __post_init__(self) -> None:
        if (
            not self.evidence_id
            or not self.active_ingredient_id
            or not self.analyte_id
            or self.analyte_type not in ANALYTE_TYPES
            or self.concentration_basis not in CONCENTRATION_BASES
            or self.pk_metric not in PK_METRICS
            or not self.time_above_threshold
            or not self.schedule
            or not self.duration
            or not self.tissue
            or self.location not in LOCATIONS
            or len(self.active_metabolite_ids) != len(set(self.active_metabolite_ids))
        ):
            raise ExposureNormalizationError("exposure evidence is incomplete")


@dataclass(frozen=True, slots=True)
class ActivityExposureComparison:
    activity_evidence_id: str
    exposure_evidence_id: str
    active_ingredient_id: str
    analyte_id: str
    potency_nm: float
    exposure_nm: float
    potency_qualifier: str
    exposure_qualifier: str
    potency_exposure_ratio: float | None
    normalized_response_range_nm: tuple[float, float]
    pk_metric: str
    concentration_basis: str
    time_above_threshold: str
    schedule: str
    duration: str
    tissue: str
    location: str
    active_metabolite_ids: tuple[str, ...]
    target_occupancy: EffectMeasure | None
    maximum_effect: EffectMeasure | None
    partial_modulation: bool | None
    minimum_endpoint_relevant_modulation: EffectMeasure | None
    effect_magnitude_state: str


def compare_activity_to_exposure(
    activity: ActivityEvidence, exposure: ExposureEvidence
) -> ActivityExposureComparison:
    if (
        activity.active_ingredient_id != exposure.active_ingredient_id
        or activity.analyte_id != exposure.analyte_id
        or activity.analyte_type != exposure.analyte_type
        or activity.concentration_basis != exposure.concentration_basis
        or activity.tissue != exposure.tissue
        or activity.location != exposure.location
    ):
        raise ExposureNormalizationError(
            "incompatible analyte, basis, tissue, or intracellular context"
        )
    exact_enough = {activity.potency.qualifier, exposure.concentration.qualifier} <= {
        "=",
        "approximately",
    }
    ratio = (
        exposure.concentration.normalized_nm / activity.potency.normalized_nm
        if exact_enough and activity.potency.normalized_nm > 0
        else None
    )
    maximum = activity.maximum_effect
    required = activity.minimum_endpoint_relevant_modulation
    if maximum is None or required is None:
        effect_state = "magnitude_unknown"
    elif maximum.unit != required.unit:
        raise ExposureNormalizationError(
            "potency and incompatible effect-magnitude values cannot be combined"
        )
    elif activity.result_direction != activity.required_endpoint_direction:
        effect_state = "wrong_direction"
    elif maximum.value < required.value:
        effect_state = "below_endpoint_relevant_magnitude"
    else:
        effect_state = "endpoint_relevant_magnitude_supported"
    occupancy = exposure.target_occupancy or activity.target_occupancy
    if (
        exposure.target_occupancy is not None
        and activity.target_occupancy is not None
        and exposure.target_occupancy.unit != activity.target_occupancy.unit
    ):
        raise ExposureNormalizationError("target occupancy units are incompatible")
    return ActivityExposureComparison(
        activity_evidence_id=activity.evidence_id,
        exposure_evidence_id=exposure.evidence_id,
        active_ingredient_id=activity.active_ingredient_id,
        analyte_id=activity.analyte_id,
        potency_nm=activity.potency.normalized_nm,
        exposure_nm=exposure.concentration.normalized_nm,
        potency_qualifier=activity.potency.qualifier,
        exposure_qualifier=exposure.concentration.qualifier,
        potency_exposure_ratio=ratio,
        normalized_response_range_nm=tuple(
            item.normalized_nm for item in activity.concentration_response_range
        ),
        pk_metric=exposure.pk_metric,
        concentration_basis=exposure.concentration_basis,
        time_above_threshold=exposure.time_above_threshold,
        schedule=exposure.schedule,
        duration=exposure.duration,
        tissue=exposure.tissue,
        location=exposure.location,
        active_metabolite_ids=exposure.active_metabolite_ids,
        target_occupancy=occupancy,
        maximum_effect=maximum,
        partial_modulation=activity.partial_modulation,
        minimum_endpoint_relevant_modulation=required,
        effect_magnitude_state=effect_state,
    )
