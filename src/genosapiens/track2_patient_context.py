"""Explicitly authorized structured patient context, never raw VCF-derived PGx."""

from __future__ import annotations

from dataclasses import dataclass
from datetime import datetime, timezone
from typing import Callable, Mapping


CONTEXT_DIMENSIONS = {
    "comorbidities",
    "current_medications",
    "organ_function",
    "pharmacogenomics",
    "preferences",
}


class PatientContextError(ValueError):
    """Raised when optional patient context is unauthorized, raw, or ambiguous."""


@dataclass(frozen=True, slots=True)
class PatientContextAuthorization:
    authorization_id: str
    authorized_dimensions: tuple[str, ...]
    purpose: str
    reviewer_role: str
    reviewed_at: str

    def __post_init__(self) -> None:
        try:
            timestamp = datetime.fromisoformat(self.reviewed_at.replace("Z", "+00:00"))
        except (TypeError, ValueError) as error:
            raise PatientContextError("authorization review timestamp is invalid") from error
        if (
            not self.authorization_id
            or not self.authorized_dimensions
            or set(self.authorized_dimensions) - CONTEXT_DIMENSIONS
            or self.purpose != "track2_candidate_context"
            or not self.reviewer_role
            or timestamp.tzinfo is None
            or timestamp.utcoffset() != timezone.utc.utcoffset(timestamp)
        ):
            raise PatientContextError("patient-context authorization is incomplete")


@dataclass(frozen=True, slots=True)
class PatientContextDimension:
    name: str
    status: str
    value: object | None
    evidence_id: str | None


@dataclass(frozen=True, slots=True)
class OptionalPatientContext:
    dimensions: tuple[PatientContextDimension, ...]
    patient_specific_suitability: str
    low_risk_inferred_from_absence: bool
    wild_type_inferred_from_absence: bool
    authorization_id: str | None


def not_assessed_patient_context() -> OptionalPatientContext:
    return OptionalPatientContext(
        dimensions=tuple(
            PatientContextDimension(name, "not_assessed", None, None)
            for name in sorted(CONTEXT_DIMENSIONS)
        ),
        patient_specific_suitability="not_assessed",
        low_risk_inferred_from_absence=False,
        wild_type_inferred_from_absence=False,
        authorization_id=None,
    )


def _contains_vcf_reference(value: object) -> bool:
    if isinstance(value, Mapping):
        return any(
            str(key).lower() in {"vcf", "vcf_path", "raw_vcf"}
            or _contains_vcf_reference(item)
            for key, item in value.items()
        )
    if isinstance(value, (list, tuple)):
        return any(_contains_vcf_reference(item) for item in value)
    return isinstance(value, str) and value.lower().endswith((".vcf", ".vcf.gz"))


def ingest_authorized_patient_context(
    *,
    source_kind: str,
    payload_loader: Callable[[], Mapping[str, object]],
    authorization: PatientContextAuthorization,
) -> OptionalPatientContext:
    if source_kind != "reviewed_structured_input":
        raise PatientContextError("raw VCF or unreviewed patient context is prohibited")
    payload = payload_loader()
    if not isinstance(payload, Mapping) or set(payload) - CONTEXT_DIMENSIONS:
        raise PatientContextError("structured patient-context fields are undefined")
    if _contains_vcf_reference(payload):
        raise PatientContextError("patient context cannot reference or read a VCF")
    unauthorized = set(payload) - set(authorization.authorized_dimensions)
    if unauthorized:
        raise PatientContextError("patient-context dimension is not authorized")
    dimensions: list[PatientContextDimension] = []
    for name in sorted(CONTEXT_DIMENSIONS):
        item = payload.get(name)
        if item is None:
            dimensions.append(PatientContextDimension(name, "not_assessed", None, None))
            continue
        if (
            not isinstance(item, Mapping)
            or set(item) != {"evidence_id", "value"}
            or not item["evidence_id"]
            or item["value"] is None
            or item["value"] == ""
        ):
            raise PatientContextError("patient-context dimension must be structured and evidenced")
        dimensions.append(
            PatientContextDimension(name, "assessed", item["value"], str(item["evidence_id"]))
        )
    return OptionalPatientContext(
        dimensions=tuple(dimensions),
        patient_specific_suitability="assessed_not_decided",
        low_risk_inferred_from_absence=False,
        wild_type_inferred_from_absence=False,
        authorization_id=authorization.authorization_id,
    )
