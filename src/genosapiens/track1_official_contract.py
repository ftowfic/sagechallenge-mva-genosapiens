"""Closed validation of the final public Track 1 instructions and evaluator review."""

from __future__ import annotations

import hashlib
from pathlib import Path
import re
from types import MappingProxyType
from typing import Mapping

from genosapiens.bounded_json import BoundedJSONError, load_bounded_json


SCHEMA = "genosapiens.track1-official-contract-review"
REPOSITORY = "SageBio/rare-disease-real-kid-mva-hackathon-2026"
EXPECTED_COLUMNS = (
    "proband_id", "chrom_1", "pos_1", "ref_1", "alt_1",
    "chrom_2", "pos_2", "ref_2", "alt_2", "epcr", "finding_type", "notes",
)
COMPETITIVE_CLASSES = (
    "mitochondrial_sequence", "mosaic_sequence", "regulatory_sequence", "sequence",
)
AUDIT_ONLY_CLASSES = ("copy_number", "repeat_locus", "structural")


class OfficialContractError(ValueError):
    """The reviewed public contract is missing, malformed, or stale."""


def load_official_contract_review(
    review_path: str | Path,
    *,
    evaluator_pin_path: str | Path,
) -> Mapping[str, object]:
    """Load one exact review and prove its evaluator identity against the final pin."""

    try:
        value, raw = load_bounded_json(review_path, maximum_bytes=65536)
        pin, _pin_raw = load_bounded_json(evaluator_pin_path, maximum_bytes=65536)
        if type(value) is not dict or set(value) != {
            "schema", "review_id", "reviewed_at", "space", "assets",
            "submission_contract", "evaluator_contract", "class_representation",
            "operator_facts", "organizer_clarification",
        }:
            raise ValueError
        if value["schema"] != SCHEMA:
            raise ValueError
        space = value["space"]
        if type(space) is not dict or set(space) != {"repository", "revision", "last_modified"}:
            raise ValueError
        revision = space["revision"]
        if space["repository"] != REPOSITORY or re.fullmatch(r"[0-9a-f]{40}", revision) is None:
            raise ValueError

        assets = value["assets"]
        if type(assets) is not list or not assets:
            raise ValueError
        roles: dict[str, Mapping[str, object]] = {}
        paths: set[str] = set()
        for asset in assets:
            if type(asset) is not dict or set(asset) != {
                "role", "path", "sha256", "git_blob_sha1", "size_bytes",
            }:
                raise ValueError
            if (
                asset["role"] in roles
                or asset["path"] in paths
                or re.fullmatch(r"[0-9a-f]{64}", asset["sha256"]) is None
                or re.fullmatch(r"[0-9a-f]{40}", asset["git_blob_sha1"]) is None
                or type(asset["size_bytes"]) is not int
                or asset["size_bytes"] <= 0
            ):
                raise ValueError
            roles[asset["role"]] = asset
            paths.add(asset["path"])

        evaluator = roles["evaluator"]
        evaluator_contract = value["evaluator_contract"]
        if (
            pin["schema"] != "genosapiens.challenge-evaluator-pin"
            or pin["source"]["repository"] != REPOSITORY
            or pin["source"]["revision"] != revision
            or pin["source"]["path"] != evaluator["path"]
            or pin["source"]["sha256"] != evaluator["sha256"]
            or pin["source"]["git_blob_sha1"] != evaluator["git_blob_sha1"]
            or pin["source"]["size_bytes"] != evaluator["size_bytes"]
            or evaluator_contract["source_revision"] != revision
            or evaluator_contract["source_sha256"] != evaluator["sha256"]
            or evaluator_contract["selected_submission_policy"] != "primary-only"
        ):
            raise ValueError

        submission = value["submission_contract"]
        classes = value["class_representation"]
        clarification = value["organizer_clarification"]
        operator = value["operator_facts"]
        if (
            tuple(submission["columns"]) != EXPECTED_COLUMNS
            or submission["candidate_row_limit"] != 10
            or submission["coordinate_assembly"] != "GRCh38"
            or submission["epcr_interval"] != "(0,1]"
            or submission["proband_id"] != "PROBAND01"
            or set(submission["finding_types"]) != {"primary", "secondary"}
            or tuple(classes["competitive_exact_sequence_classes"]) != COMPETITIVE_CLASSES
            or tuple(classes["audit_only_without_exact_sequence_tuple"]) != AUDIT_ONLY_CLASSES
            or classes["public_native_ref_check"] != "required-additional-check-without-input-reference-repair"
            or classes["unsupported_proxy_policy"] != "never-render-a-proxy-sequence-allele"
            or clarification["contact_attempted"] is not False
            or clarification["status"] != "not-requested-no-contact-authorized"
            or operator["retention_due_date"] != "unresolved-authoritative-calendar-date"
            or operator["video_status"] != "operator-confirmation-required"
        ):
            raise ValueError
        digest = hashlib.sha256(raw).hexdigest()
        return MappingProxyType({"review_sha256": digest, **value})
    except (BoundedJSONError, KeyError, TypeError, ValueError):
        raise OfficialContractError("official Track 1 contract review is invalid") from None
