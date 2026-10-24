"""Reviewed manual imports for evidence without reliable structured interfaces."""

from __future__ import annotations

from dataclasses import asdict, dataclass
from datetime import date, datetime
import hashlib
import json
import re
from typing import Any


SHA256 = re.compile(r"^[0-9a-f]{64}$")
DOCUMENT_TYPES = {"label", "pharmacokinetic_study", "regulator_document", "safety_record", "primary_publication", "dataset_record"}


class ManualImportError(ValueError):
    """Raised when a reviewed manual extraction is not reproducible."""


@dataclass(frozen=True, slots=True)
class DocumentLocator:
    page: str | None
    section: str | None

    def __post_init__(self) -> None:
        if not self.page and not self.section:
            raise ManualImportError("manual import requires a page or section locator")


@dataclass(frozen=True, slots=True)
class ExtractedField:
    name: str
    raw_text: str
    normalized_value: Any
    unit: str | None
    transcription_uncertainty: str
    interpretation_uncertainty: str

    def __post_init__(self) -> None:
        if (
            not self.name
            or not self.raw_text
            or self.normalized_value is None
            or not self.transcription_uncertainty
            or not self.interpretation_uncertainty
        ):
            raise ManualImportError(
                "manual extracted fields require transcription and interpretation uncertainty"
            )
        try:
            json.dumps(self.normalized_value, sort_keys=True)
        except (TypeError, ValueError) as error:
            raise ManualImportError("normalized manual value must be JSON-compatible") from error


@dataclass(frozen=True, slots=True)
class ManualEvidenceImport:
    import_id: str
    document_type: str
    document_identity: str
    document_version: str
    document_checksum: str
    document_url: str
    accessed_on: str
    locator: DocumentLocator
    extracted_fields: tuple[ExtractedField, ...]
    reviewer_role: str
    review_record_id: str
    review_time: str
    structured_interface_status: str
    overall_transcription_uncertainty: str
    overall_interpretation_uncertainty: str

    def __post_init__(self) -> None:
        try:
            date.fromisoformat(self.accessed_on)
            reviewed = datetime.fromisoformat(self.review_time.replace("Z", "+00:00"))
        except (TypeError, ValueError) as error:
            raise ManualImportError("manual import dates must use ISO formats") from error
        if (
            self.document_type not in DOCUMENT_TYPES
            or not self.document_identity
            or not self.document_version
            or not SHA256.fullmatch(self.document_checksum)
            or not self.document_url.startswith("https://")
            or not self.extracted_fields
            or len({field.name for field in self.extracted_fields})
            != len(self.extracted_fields)
            or not self.reviewer_role
            or not self.review_record_id
            or reviewed.tzinfo is None
            or self.structured_interface_status not in {"unavailable_or_unreliable", "reviewed_manual_extraction"}
            or not self.overall_transcription_uncertainty
            or not self.overall_interpretation_uncertainty
        ):
            raise ManualImportError("manual evidence import is incomplete")
        expected_id = make_manual_import_id(
            self.document_type,
            self.document_identity,
            self.document_version,
            self.document_checksum,
            self.locator,
        )
        if self.import_id != expected_id:
            raise ManualImportError("manual import ID is not stable")

    def as_track2_payload(self) -> dict[str, Any]:
        return {
            "document_checksum": self.document_checksum,
            "document_identity": self.document_identity,
            "document_version": self.document_version,
            "extracted_fields": [asdict(field) for field in self.extracted_fields],
            "locator": asdict(self.locator),
            "review_time": self.review_time,
            "reviewer_role": self.reviewer_role,
            "uncertainty": {
                "transcription": self.overall_transcription_uncertainty,
                "interpretation": self.overall_interpretation_uncertainty,
            },
        }


def make_manual_import_id(
    document_type: str,
    document_identity: str,
    document_version: str,
    document_checksum: str,
    locator: DocumentLocator,
) -> str:
    content = json.dumps(
        [
            document_type,
            document_identity,
            document_version,
            document_checksum,
            asdict(locator),
        ],
        sort_keys=True,
        separators=(",", ":"),
    )
    return "manual-import:" + hashlib.sha256(content.encode()).hexdigest()[:24]


def make_manual_evidence_import(
    *,
    document_type: str,
    document_identity: str,
    document_version: str,
    document_checksum: str,
    document_url: str,
    accessed_on: str,
    locator: DocumentLocator,
    extracted_fields: tuple[ExtractedField, ...],
    reviewer_role: str,
    review_record_id: str,
    review_time: str,
    structured_interface_status: str,
    overall_transcription_uncertainty: str,
    overall_interpretation_uncertainty: str,
) -> ManualEvidenceImport:
    return ManualEvidenceImport(
        import_id=make_manual_import_id(
            document_type,
            document_identity,
            document_version,
            document_checksum,
            locator,
        ),
        document_type=document_type,
        document_identity=document_identity,
        document_version=document_version,
        document_checksum=document_checksum,
        document_url=document_url,
        accessed_on=accessed_on,
        locator=locator,
        extracted_fields=extracted_fields,
        reviewer_role=reviewer_role,
        review_record_id=review_record_id,
        review_time=review_time,
        structured_interface_status=structured_interface_status,
        overall_transcription_uncertainty=overall_transcription_uncertainty,
        overall_interpretation_uncertainty=overall_interpretation_uncertainty,
    )
