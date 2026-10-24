"""Canonical Track 1 mechanism-ready handoff validation."""

from __future__ import annotations

from dataclasses import dataclass
import hashlib
import json
from pathlib import Path
import re
from typing import Any, Mapping

from genosapiens.json_schema import JsonSchemaError, validate_json_schema


from genosapiens.wgs_handoff import (
    HANDOFF_SCHEMA, MODES, ShardedHandoffManifest, handoff_manifest_from_json,
    validate_handoff_bundle,
)
CONTRACT_SCHEMA = "genosapiens.track1-handoff-contract"
SHA256 = re.compile(r"^[0-9a-f]{64}$")


class HandoffError(ValueError):
    """Raised when the canonical Track 1 handoff is invalid."""


@dataclass(frozen=True, slots=True)
class PublishedHandoff:
    path: Path
    sha256: str
    size_bytes: int


def _json_object(path: str | Path, label: str) -> dict[str, Any]:
    source = Path(path)
    try:
        value = json.loads(source.read_text(encoding="utf-8"))
    except OSError as error:
        raise HandoffError(f"cannot read {label} {source}: {error}") from error
    except json.JSONDecodeError as error:
        raise HandoffError(f"invalid JSON in {label} {source}: {error}") from error
    if not isinstance(value, dict):
        raise HandoffError(f"{label} must contain a JSON object")
    return value




def public_evidence_reference_id(
    *, evidence_id: str, source_record_id: str, manifest_sha256: str,
) -> str:
    """Project-owned reference to one unchanged native public record tuple."""
    if (not isinstance(evidence_id, str) or not evidence_id.strip()
            or not isinstance(source_record_id, str) or not source_record_id.strip()
            or not isinstance(manifest_sha256, str) or SHA256.fullmatch(manifest_sha256) is None):
        raise HandoffError("native public evidence reference tuple is incomplete")
    content = (json.dumps(
        [evidence_id, source_record_id, manifest_sha256],
        separators=(",", ":"), ensure_ascii=True, allow_nan=False,
    ) + "\n").encode("ascii")
    return "track1-evidence:" + hashlib.sha256(content).hexdigest()


@dataclass(frozen=True, slots=True)
class HandoffContract:
    canonical_schema: str
    canonical_artifact_path: str
    allowed_analysis_modes: tuple[str, ...]
    forbidden_substitutes: tuple[str, ...]

    @classmethod
    def from_mapping(cls, data: Mapping[str, Any]) -> "HandoffContract":
        required = {
            "allowed_analysis_modes",
            "canonical_artifact_path",
            "canonical_schema",
            "contract_owner",
            "forbidden_substitutes",
            "schema",
        }
        if set(data) != required:
            raise HandoffError("handoff contract fields are incomplete or unknown")
        if data.get("schema") != CONTRACT_SCHEMA:
            raise HandoffError("unsupported handoff contract schema")
        if data.get("contract_owner") != "track1":
            raise HandoffError("the handoff compatibility policy must be Track 1-owned")
        if data.get("canonical_schema") != HANDOFF_SCHEMA:
            raise HandoffError("handoff contract does not identify the canonical schema")
        modes = data.get("allowed_analysis_modes")
        substitutes = data.get("forbidden_substitutes")
        if not isinstance(modes, list) or not modes or len(modes) != len(set(modes)):
            raise HandoffError("allowed analysis modes must be a non-empty unique list")
        if not isinstance(substitutes, list) or set(substitutes) != {
            "case_vcf",
            "challenge_csv",
            "manually_authored_handoff",
            "track2_owned_schema",
        }:
            raise HandoffError("forbidden handoff substitutes are incomplete")
        artifact_path = data.get("canonical_artifact_path")
        if (
            not isinstance(artifact_path, str)
            or artifact_path
            != "controlled-work/track1/handoff/mechanism-ready-handoff.json"
        ):
            raise HandoffError("canonical handoff path must be relative to the selected run root")
        return cls(
            canonical_schema=HANDOFF_SCHEMA,
            canonical_artifact_path=artifact_path,
            allowed_analysis_modes=tuple(modes),
            forbidden_substitutes=tuple(substitutes),
        )


def load_handoff_contract(path: str | Path) -> HandoffContract:
    return HandoffContract.from_mapping(_json_object(path, "handoff contract"))


def validate_track1_handoff(manifest, schema, contract):
    """Validate current manifest metadata; integrity separately replays shards."""
    if (type(manifest) is not ShardedHandoffManifest
            or type(contract) is not HandoffContract
            or manifest.schema != contract.canonical_schema
            or manifest.canonical_artifact_path != contract.canonical_artifact_path
            or set(contract.allowed_analysis_modes) != set(MODES)):
        raise HandoffError("the current sharded handoff contract is required")
    try:
        validate_json_schema(json.loads(manifest.as_bytes()), schema)
    except JsonSchemaError:
        raise HandoffError("sharded handoff metadata violates its pinned schema") from None


def validate_handoff_bytes(content: bytes, *, schema_path: Path, contract_path: Path):
    """Validate metadata against the exact current Track 1 schema and contract."""
    try:
        if len(content) > 8 * 1024 * 1024:
            raise ValueError("handoff manifest exceeds capacity")
        manifest = handoff_manifest_from_json(content,
            expected_sha256=hashlib.sha256(content).hexdigest())
        if (manifest.schema_sha256 != hashlib.sha256(schema_path.read_bytes()).hexdigest()
                or manifest.contract_sha256 != hashlib.sha256(contract_path.read_bytes()).hexdigest()):
            raise ValueError("handoff schema or contract binding changed")
        validate_track1_handoff(manifest, _json_object(schema_path, "Track 1 handoff schema"),
            load_handoff_contract(contract_path))
        return manifest
    except (ValueError, RuntimeError, TypeError):
        raise HandoffError("current sharded handoff binding is invalid") from None






