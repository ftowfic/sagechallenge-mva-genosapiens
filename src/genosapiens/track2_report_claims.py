"""Typed, evidence-resolved claims for generated Track 2 reports."""

from __future__ import annotations

from dataclasses import dataclass, replace
import json
import math
import re
from typing import Any, Mapping

from genosapiens.track2_records import make_record_id
from genosapiens.track2_report_tables import Track2ReportTables


REPORT_CLAIM_SCHEMA = "genosapiens.track2-report-claims/v1"
CLAIM_KINDS = {
    "assumption",
    "computational_hypothesis",
    "computed_result",
    "inference",
    "manual_judgment",
    "sourced_fact",
}
PUBLICATION_CLASSIFICATIONS = {"controlled", "public"}
CLAIM_ID = re.compile(r"^report-claim:[0-9a-f]{24}$")


class ReportClaimError(ValueError):
    """Raised when a material report claim is unlabeled or unsupported."""


def _validate_json_value(value: Any, path: str = "value") -> None:
    if value is None or isinstance(value, (str, bool, int)):
        return
    if isinstance(value, float):
        if not math.isfinite(value):
            raise ReportClaimError(f"claim {path} contains a non-finite number")
        return
    if isinstance(value, list):
        for index, item in enumerate(value):
            _validate_json_value(item, f"{path}[{index}]")
        return
    if isinstance(value, Mapping):
        if any(not isinstance(key, str) for key in value):
            raise ReportClaimError(f"claim {path} contains a non-string object key")
        for key in sorted(value):
            _validate_json_value(value[key], f"{path}.{key}")
        return
    raise ReportClaimError(f"claim {path} is not JSON-compatible")


def _plain_json(value: Any) -> Any:
    if isinstance(value, Mapping):
        return {key: _plain_json(value[key]) for key in sorted(value)}
    if isinstance(value, list):
        return [_plain_json(item) for item in value]
    return value


@dataclass(frozen=True, slots=True)
class ReportClaim:
    claim_id: str
    claim_kind: str
    text: str | None
    value: Any | None
    evidence_ids: tuple[str, ...]
    publication_classification: str

    def __post_init__(self) -> None:
        if CLAIM_ID.fullmatch(self.claim_id) is None:
            raise ReportClaimError("report claim ID is not stable or typed")
        if self.claim_kind not in CLAIM_KINDS:
            raise ReportClaimError("report claim kind is undefined")
        if self.publication_classification not in PUBLICATION_CLASSIFICATIONS:
            raise ReportClaimError("report claim publication classification is undefined")
        has_text = isinstance(self.text, str) and bool(self.text.strip())
        has_value = self.value is not None
        if has_text == has_value:
            raise ReportClaimError("report claim requires exactly one text or value")
        if self.text is not None and not has_text:
            raise ReportClaimError("report claim text cannot be blank")
        if (
            not self.evidence_ids
            or len(self.evidence_ids) != len(set(self.evidence_ids))
            or any(not isinstance(item, str) or not item for item in self.evidence_ids)
        ):
            raise ReportClaimError("every report claim requires unique evidence IDs")
        if has_value:
            _validate_json_value(self.value)

    def as_track2_payload(self) -> dict[str, Any]:
        """Return the exact payload declared by the Track 2 record registry."""
        return {
            "claim_kind": self.claim_kind,
            "publication_classification": self.publication_classification,
            "text": self.text,
            "value": _plain_json(self.value),
        }

    def as_dict(self) -> dict[str, Any]:
        return {
            "claim_id": self.claim_id,
            **self.as_track2_payload(),
            "evidence_ids": list(self.evidence_ids),
        }


@dataclass(frozen=True, slots=True)
class ReportClaimBundle:
    claims: tuple[ReportClaim, ...]
    evidence_catalog_ids: tuple[str, ...]

    def as_dict(self) -> dict[str, Any]:
        return {
            "schema": REPORT_CLAIM_SCHEMA,
            "evidence_catalog_ids": list(self.evidence_catalog_ids),
            "claims": [claim.as_dict() for claim in self.claims],
        }


def make_report_claim(
    *,
    source_identity: str,
    role: str,
    claim_kind: str,
    evidence_ids: tuple[str, ...],
    publication_classification: str,
    text: str | None = None,
    value: Any | None = None,
) -> ReportClaim:
    """Create a claim whose identity does not depend on prose or row order."""
    return ReportClaim(
        claim_id=make_record_id("report_claim", source_identity, role),
        claim_kind=claim_kind,
        text=text,
        value=_plain_json(value),
        evidence_ids=tuple(sorted(evidence_ids)),
        publication_classification=publication_classification,
    )


def build_report_claim_bundle(
    claims: tuple[ReportClaim, ...],
    *,
    report_tables: Track2ReportTables,
    excluded_evidence_ids: tuple[str, ...] = (),
) -> ReportClaimBundle:
    """Resolve every material claim to usable report evidence and sort it."""
    if not claims:
        raise ReportClaimError("at least one typed report claim is required")
    catalog_ids = tuple(
        sorted(reference.evidence_id for reference in report_tables.evidence_catalog)
    )
    if not catalog_ids or len(catalog_ids) != len(set(catalog_ids)):
        raise ReportClaimError("report evidence catalog is empty or ambiguous")
    excluded = set(excluded_evidence_ids)
    normalized: list[ReportClaim] = []
    seen_claim_ids: set[str] = set()
    for claim in claims:
        if claim.claim_id in seen_claim_ids:
            raise ReportClaimError("report claim IDs must be unique")
        seen_claim_ids.add(claim.claim_id)
        unresolved = set(claim.evidence_ids) - set(catalog_ids)
        if unresolved:
            raise ReportClaimError("report claim has unresolved evidence")
        if set(claim.evidence_ids) & excluded:
            raise ReportClaimError("report claim relies on excluded evidence")
        normalized.append(
            replace(
                claim,
                value=_plain_json(claim.value),
                evidence_ids=tuple(sorted(claim.evidence_ids)),
            )
        )
    return ReportClaimBundle(
        claims=tuple(sorted(normalized, key=lambda item: item.claim_id)),
        evidence_catalog_ids=catalog_ids,
    )


def report_claims_as_json(bundle: ReportClaimBundle) -> str:
    """Serialize typed claims deterministically for report generation."""
    return json.dumps(bundle.as_dict(), indent=2, sort_keys=True, ensure_ascii=True) + "\n"
