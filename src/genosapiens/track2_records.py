"""Canonical Track 2 record envelopes, hashes, and reference integrity."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path, PurePosixPath
import re
from typing import Any, Iterable, Mapping

from genosapiens.json_schema import JsonSchemaError, validate_json_schema


REGISTRY_SCHEMA = "genosapiens.track2-record-registry/v1"
RECORD_SCHEMA = "genosapiens.track2-record/v1"
SHA256 = re.compile(r"^[0-9a-f]{64}$")


class Track2RecordError(ValueError):
    """Raised when a record is unstable, ambiguous, or unresolved."""


def _json_object(path: str | Path, label: str) -> dict[str, Any]:
    try:
        data = json.loads(Path(path).read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as error:
        raise Track2RecordError(f"cannot load {label}: {error}") from error
    if not isinstance(data, dict):
        raise Track2RecordError(f"{label} must be a JSON object")
    return data


def load_record_registry(path: str | Path) -> dict[str, Any]:
    registry = _json_object(path, "Track 2 record registry")
    required = {
        "artifact_roots",
        "missing_states",
        "public_manifest_forbidden_fields",
        "record_schema",
        "record_types",
        "schema",
    }
    if set(registry) != required or registry.get("schema") != REGISTRY_SCHEMA:
        raise Track2RecordError("record registry is incomplete or incompatible")
    if registry.get("record_schema") != RECORD_SCHEMA:
        raise Track2RecordError("record registry identifies the wrong envelope")
    record_types = registry.get("record_types")
    if not isinstance(record_types, dict) or not record_types:
        raise Track2RecordError("record type definitions are required")
    for name, definition in record_types.items():
        if not isinstance(definition, Mapping) or set(definition) != {
            "classification",
            "id_prefix",
            "nullable_payload_fields",
            "reference_roles",
            "required_payload_fields",
        }:
            raise Track2RecordError(f"record type {name} is incomplete")
        if not isinstance(definition["required_payload_fields"], list):
            raise Track2RecordError(f"record type {name} payload fields must be a list")
        if not set(definition["nullable_payload_fields"]).issubset(
            definition["required_payload_fields"]
        ):
            raise Track2RecordError(f"record type {name} has undeclared nullable fields")
    return registry


def make_record_id(record_type: str, source_identity: str, role: str) -> str:
    if not record_type or not source_identity or not role:
        raise Track2RecordError("stable record IDs require type, source identity, and role")
    digest = hashlib.sha256(
        f"{record_type}\0{source_identity}\0{role}".encode()
    ).hexdigest()[:24]
    prefix = record_type.replace("_", "-")
    return f"{prefix}:{digest}"


def _canonicalize(value: Any) -> Any:
    if isinstance(value, Mapping):
        return {key: _canonicalize(value[key]) for key in sorted(value)}
    if isinstance(value, list):
        normalized = [_canonicalize(item) for item in value]
        return sorted(
            normalized,
            key=lambda item: json.dumps(item, sort_keys=True, separators=(",", ":")),
        )
    if value is None or isinstance(value, (str, int, float, bool)):
        return value
    raise Track2RecordError(f"record contains non-JSON value: {type(value).__name__}")


def canonical_record_bytes(record: Mapping[str, Any]) -> bytes:
    semantic = {key: value for key, value in record.items() if key != "content_sha256"}
    return (
        json.dumps(_canonicalize(semantic), sort_keys=True, separators=(",", ":")) + "\n"
    ).encode()


def record_content_sha256(record: Mapping[str, Any]) -> str:
    return hashlib.sha256(canonical_record_bytes(record)).hexdigest()


def add_record_hash(record: Mapping[str, Any]) -> dict[str, Any]:
    hashed = dict(record)
    hashed["content_sha256"] = record_content_sha256(hashed)
    return hashed


def canonical_hashed_record_bytes(record: Mapping[str, Any]) -> bytes:
    return (
        json.dumps(_canonicalize(record), sort_keys=True, separators=(",", ":")) + "\n"
    ).encode()


def _walk_public_values(value: Any, forbidden_fields: set[str]) -> None:
    if isinstance(value, Mapping):
        normalized_forbidden = {
            re.sub(r"[^a-z0-9]", "", item.casefold()) for item in forbidden_fields
        }
        leaked = sorted(
            key
            for key in value
            if re.sub(r"[^a-z0-9]", "", str(key).casefold())
            in normalized_forbidden
        )
        if leaked:
            raise Track2RecordError(
                f"public manifest contains controlled identity field: {leaked[0]}"
            )
        for item in value.values():
            _walk_public_values(item, forbidden_fields)
    elif isinstance(value, list):
        for item in value:
            _walk_public_values(item, forbidden_fields)
    elif isinstance(value, str) and (
        "controlled-work/" in value or value.startswith("/handoff/")
    ):
        raise Track2RecordError("public manifest contains a controlled path or identity")


def validate_track2_record(
    record: Mapping[str, Any],
    envelope_schema: Mapping[str, Any],
    registry: Mapping[str, Any],
) -> None:
    if record.get("schema") != RECORD_SCHEMA:
        raise Track2RecordError("unsupported Track 2 record contract")
    try:
        validate_json_schema(record, envelope_schema)
    except JsonSchemaError as error:
        raise Track2RecordError(f"record envelope validation failed: {error}") from error
    record_type = record["record_type"]
    definition = registry["record_types"].get(record_type)
    if not isinstance(definition, Mapping):
        raise Track2RecordError(f"unregistered Track 2 record type: {record_type}")
    expected_id = rf"^{re.escape(definition['id_prefix'])}:[0-9a-f]{{24}}$"
    if re.fullmatch(expected_id, record["record_id"]) is None:
        raise Track2RecordError("record ID is unstable or has the wrong typed prefix")
    if record["classification"] not in definition["classification"]:
        raise Track2RecordError("record classification is not allowed for its type")
    roots = registry["artifact_roots"][record["classification"]]
    artifact_path = PurePosixPath(record["artifact_path"])
    if (
        artifact_path.is_absolute()
        or ".." in artifact_path.parts
        or "\\" in record["artifact_path"]
        or not any(record["artifact_path"].startswith(root) for root in roots)
    ):
        raise Track2RecordError("record artifact path is outside its classification root")
    if record["missing_state"] not in registry["missing_states"]:
        raise Track2RecordError("record uses an ambiguous missing state")
    payload = record["payload"]
    required_fields = set(definition["required_payload_fields"])
    if set(payload) != required_fields:
        raise Track2RecordError("record payload has omitted or unknown fields")
    nullable = set(definition["nullable_payload_fields"])
    null_fields = {key for key, value in payload.items() if value is None}
    if not null_fields.issubset(nullable):
        raise Track2RecordError("record payload contains an undeclared null")
    if record_content_sha256(record) != record["content_sha256"]:
        raise Track2RecordError("record content is noncanonical or its hash is invalid")
    if record["record_type"] == "research_run_manifest":
        _walk_public_values(
            payload, set(registry["public_manifest_forbidden_fields"])
        )


def validate_record_bundle(
    records: Iterable[Mapping[str, Any]],
    envelope_schema: Mapping[str, Any],
    registry: Mapping[str, Any],
) -> None:
    indexed: dict[str, Mapping[str, Any]] = {}
    records_tuple = tuple(records)
    for record in records_tuple:
        validate_track2_record(record, envelope_schema, registry)
        if record["record_id"] in indexed:
            raise Track2RecordError(f"duplicate record ID: {record['record_id']}")
        indexed[record["record_id"]] = record
    for record in records_tuple:
        definition = registry["record_types"][record["record_type"]]
        references = record["references"]
        roles: dict[str, list[Mapping[str, Any]]] = {}
        for reference in references:
            roles.setdefault(reference["role"], []).append(reference)
            if reference["external"]:
                if not isinstance(reference["manifest_sha256"], str) or not SHA256.fullmatch(
                    reference["manifest_sha256"]
                ):
                    raise Track2RecordError("external references require a manifest SHA-256")
            else:
                target = indexed.get(reference["target_id"])
                if target is None:
                    raise Track2RecordError(
                        f"unresolved record reference: {reference['target_id']}"
                    )
                if target["record_type"] != reference["target_type"]:
                    raise Track2RecordError("record reference target type does not match")
        declared_roles = definition["reference_roles"]
        if set(roles) - set(declared_roles):
            raise Track2RecordError("record contains an undeclared reference role")
        for role, constraint in declared_roles.items():
            role_references = roles.get(role, [])
            minimum = constraint["minimum"]
            maximum = constraint["maximum"]
            if len(role_references) < minimum or (
                maximum is not None and len(role_references) > maximum
            ):
                raise Track2RecordError(f"reference cardinality failed for role {role}")
            if any(
                reference["target_type"] not in constraint["targets"]
                for reference in role_references
            ):
                raise Track2RecordError(f"reference role {role} has an invalid target type")


def load_canonical_record(
    path: str | Path,
    envelope_schema: Mapping[str, Any],
    registry: Mapping[str, Any],
) -> dict[str, Any]:
    source = Path(path)
    try:
        content = source.read_bytes()
        record = json.loads(content)
    except (OSError, json.JSONDecodeError) as error:
        raise Track2RecordError(f"cannot load Track 2 record: {error}") from error
    if not isinstance(record, dict):
        raise Track2RecordError("Track 2 record must be an object")
    validate_track2_record(record, envelope_schema, registry)
    if content != canonical_hashed_record_bytes(record):
        raise Track2RecordError("record file is not canonical JSON")
    return record
