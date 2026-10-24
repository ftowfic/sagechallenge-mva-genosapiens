"""Explicit public/controlled artifact classification and case inventories."""

from __future__ import annotations

from dataclasses import asdict, dataclass, replace
from datetime import datetime, timezone
from enum import StrEnum
import hashlib
import json
import os
from pathlib import Path, PurePosixPath
import re
from typing import Any, Mapping
import uuid


INVENTORY_SCHEMA = "genosapiens.case-inventory"
OPERATOR_INVENTORY_SCHEMA = "genosapiens.case-inventory/v2"
IDENTIFIER = re.compile(r"^[A-Za-z0-9][A-Za-z0-9._-]*$")
SHA256 = re.compile(r"^[0-9a-f]{64}$")
PUBLIC_PREFIXES = (
    ("manifests",),
    ("releases",),
    ("derived",),
    ("scratch", "public-build"),
)
CONTROLLED_PREFIXES = (
    ("controlled-work",),
    ("scratch", "controlled"),
)


class GovernanceError(RuntimeError):
    """Raised when an artifact violates an explicit data boundary."""


class ArtifactClassification(StrEnum):
    PUBLIC = "public"
    CONTROLLED = "controlled"


def _has_prefix(path: PurePosixPath, prefixes: tuple[tuple[str, ...], ...]) -> bool:
    return any(path.parts[: len(prefix)] == prefix for prefix in prefixes)


def validate_relative_classification(
    relative_path: str,
    classification: ArtifactClassification,
) -> None:
    path = PurePosixPath(relative_path)
    if path.is_absolute() or ".." in path.parts or path == PurePosixPath("."):
        raise GovernanceError("artifact path must be a safe relative path")
    public = _has_prefix(path, PUBLIC_PREFIXES)
    controlled = _has_prefix(path, CONTROLLED_PREFIXES)
    if classification is ArtifactClassification.PUBLIC and not public:
        raise GovernanceError(f"public artifact is outside public roots: {relative_path}")
    if classification is ArtifactClassification.CONTROLLED and not controlled:
        raise GovernanceError(
            f"controlled artifact is outside controlled roots: {relative_path}"
        )
    if public and controlled:
        raise GovernanceError(f"artifact path has ambiguous classification: {relative_path}")


def _file_identity(path: Path) -> tuple[str, int]:
    digest = hashlib.sha256()
    size = 0
    with path.open("rb") as source:
        while chunk := source.read(1024 * 1024):
            digest.update(chunk)
            size += len(chunk)
    return digest.hexdigest(), size


@dataclass(frozen=True, slots=True)
class ArtifactRecord:
    artifact_id: str
    relative_path: str
    classification: ArtifactClassification
    kind: str
    sha256: str
    size_bytes: int
    retention: str
    status: str = "present"
    deleted_on: str | None = None
    operator_path: str | None = None

    def __post_init__(self) -> None:
        if not IDENTIFIER.fullmatch(self.artifact_id):
            raise GovernanceError("artifact_id is invalid")
        validate_relative_classification(self.relative_path, self.classification)
        if self.operator_path is not None:
            path = Path(self.operator_path)
            expected = 'controlled-work/operator-paths/' + hashlib.sha256(self.operator_path.encode()).hexdigest()
            if (self.classification is not ArtifactClassification.CONTROLLED
                    or not path.is_absolute() or '..' in path.parts
                    or str(path) != self.operator_path or self.relative_path != expected):
                raise GovernanceError('explicit controlled path binding is invalid')
        if not self.kind:
            raise GovernanceError("artifact kind is required")
        if not SHA256.fullmatch(self.sha256):
            raise GovernanceError("artifact SHA-256 is invalid")
        if not isinstance(self.size_bytes, int) or self.size_bytes < 0:
            raise GovernanceError("artifact size_bytes must be non-negative")
        expected_retention = (
            "retain" if self.classification is ArtifactClassification.PUBLIC else "delete"
        )
        if self.retention != expected_retention:
            raise GovernanceError(
                f"{self.classification.value} artifacts require {expected_retention!r} retention"
            )
        if self.status not in {"present", "deleted"}:
            raise GovernanceError("artifact status must be present or deleted")
        if (self.status == "deleted") != (self.deleted_on is not None):
            raise GovernanceError("deleted status and deleted_on must be set together")

    @classmethod
    def from_file(
        cls,
        reference_root: str | Path,
        path: str | Path,
        classification: ArtifactClassification,
        kind: str,
    ) -> "ArtifactRecord":
        root = Path(reference_root).resolve()
        artifact = Path(path).resolve()
        if not artifact.is_relative_to(root) or not artifact.is_file():
            raise GovernanceError("artifact must be a file under reference_root")
        relative = artifact.relative_to(root).as_posix()
        validate_relative_classification(relative, classification)
        checksum, size = _file_identity(artifact)
        artifact_id = hashlib.sha256(
            f"{classification.value}\0{relative}\0{checksum}".encode()
        ).hexdigest()[:24]
        return cls(
            artifact_id=artifact_id,
            relative_path=relative,
            classification=classification,
            kind=kind,
            sha256=checksum,
            size_bytes=size,
            retention="retain"
            if classification is ArtifactClassification.PUBLIC
            else "delete",
        )

    def mark_deleted(self, deleted_on: str | None = None) -> "ArtifactRecord":
        if self.classification is ArtifactClassification.PUBLIC:
            raise GovernanceError("public artifacts cannot be marked for controlled deletion")
        timestamp = deleted_on or datetime.now(timezone.utc).isoformat()
        return replace(self, status="deleted", deleted_on=timestamp)


@dataclass(frozen=True, slots=True)
class CaseInventory:
    case_id: str
    artifacts: tuple[ArtifactRecord, ...]
    schema: str = INVENTORY_SCHEMA

    def __post_init__(self) -> None:
        if self.schema not in {INVENTORY_SCHEMA, OPERATOR_INVENTORY_SCHEMA}:
            raise GovernanceError('inventory schema is unsupported')
        if any(record.operator_path is not None for record in self.artifacts):
            # Older readers must refuse these locations, rather than report a
            # nonexistent legacy path as successfully deleted.
            object.__setattr__(self, 'schema', OPERATOR_INVENTORY_SCHEMA)
        if not IDENTIFIER.fullmatch(self.case_id):
            raise GovernanceError("case_id is invalid")
        identifiers = [artifact.artifact_id for artifact in self.artifacts]
        paths = [artifact.relative_path for artifact in self.artifacts]
        if len(identifiers) != len(set(identifiers)):
            raise GovernanceError("inventory artifact IDs must be unique")
        if len(paths) != len(set(paths)):
            raise GovernanceError("inventory artifact paths must be unique")

    @classmethod
    def from_mapping(cls, data: Mapping[str, Any]) -> "CaseInventory":
        try:
            if (any(item.get('operator_path') is not None for item in data['artifacts'])
                    and data['schema'] != OPERATOR_INVENTORY_SCHEMA):
                raise GovernanceError('explicit operator paths require inventory schema v2')
            artifacts = tuple(
                ArtifactRecord(
                    artifact_id=item["artifact_id"],
                    relative_path=item["relative_path"],
                    classification=ArtifactClassification(item["classification"]),
                    kind=item["kind"],
                    sha256=item["sha256"],
                    size_bytes=item["size_bytes"],
                    retention=item["retention"],
                    status=item["status"],
                    deleted_on=item.get("deleted_on"),
                    operator_path=item.get("operator_path"),
                )
                for item in data["artifacts"]
            )
            return cls(
                schema=data["schema"],
                case_id=data["case_id"],
                artifacts=artifacts,
            )
        except (KeyError, TypeError, ValueError) as error:
            raise GovernanceError(f"invalid case inventory: {error}") from error


def inventory_as_bytes(inventory: CaseInventory) -> bytes:
    data = asdict(inventory)
    # Preserve the bytes of existing inventories when no explicit path is used.
    for record in data['artifacts']:
        if record['operator_path'] is None:
            del record['operator_path']
    return (
        json.dumps(data, sort_keys=True, separators=(",", ":")) + "\n"
    ).encode()


def write_inventory(inventory: CaseInventory, path: str | Path) -> None:
    destination = Path(path)
    destination.parent.mkdir(parents=True, exist_ok=True, mode=0o700)
    partial = destination.parent / f".{destination.name}.{uuid.uuid4().hex}.partial"
    try:
        with partial.open("xb") as output:
            os.chmod(partial, 0o600)
            output.write(inventory_as_bytes(inventory))
            output.flush()
            os.fsync(output.fileno())
        os.replace(partial, destination)
        os.chmod(destination, 0o600)
    finally:
        partial.unlink(missing_ok=True)


def load_inventory(path: str | Path) -> CaseInventory:
    source = Path(path)
    try:
        data = json.loads(source.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as error:
        raise GovernanceError(f"cannot load case inventory: {error}") from error
    if not isinstance(data, Mapping):
        raise GovernanceError("case inventory root must be an object")
    return CaseInventory.from_mapping(data)


@dataclass(frozen=True, slots=True)
class CleanupReport:
    case_id: str
    completed_on: str
    deleted_artifacts: tuple[str, ...]
    already_absent_artifacts: tuple[str, ...]
    identity_mismatches_deleted: tuple[str, ...]
    preserved_public_artifacts: tuple[str, ...]
    verified: bool
    schema: str = "genosapiens.deletion-verification"


def _path_for_record(root: Path, record: ArtifactRecord) -> Path:
    if record.operator_path is not None:
        path = Path(record.operator_path)
        if path.resolve() != path or any(p.is_symlink() for p in (path, *path.parents)):
            raise GovernanceError('explicit controlled artifact path changed or uses a symbolic link')
        return path
    path = root / record.relative_path
    prefixes = (
        CONTROLLED_PREFIXES
        if record.classification is ArtifactClassification.CONTROLLED
        else PUBLIC_PREFIXES
    )
    allowed_roots = [(root.joinpath(*prefix)).resolve() for prefix in prefixes]
    resolved_parent = path.parent.resolve()
    if not any(
        resolved_parent == allowed or resolved_parent.is_relative_to(allowed)
        for allowed in allowed_roots
    ):
        raise GovernanceError(
            f"artifact path escapes its {record.classification.value} root: "
            f"{record.relative_path}"
        )
    return path


def _matches_record(path: Path, record: ArtifactRecord) -> bool:
    if not path.is_file() or path.is_symlink():
        return False
    checksum, size = _file_identity(path)
    return checksum == record.sha256 and size == record.size_bytes


def cleanup_controlled_artifacts(
    inventory: CaseInventory,
    reference_root: str | Path,
    *,
    completed_on: str | None = None,
) -> tuple[CaseInventory, CleanupReport]:
    root = Path(reference_root).resolve()
    timestamp = completed_on or datetime.now(timezone.utc).isoformat()
    deleted: list[str] = []
    already_absent: list[str] = []
    identity_mismatches: list[str] = []
    preserved_public: list[str] = []

    public_paths = {_path_for_record(root, record) for record in inventory.artifacts
                    if record.classification is ArtifactClassification.PUBLIC}
    if any(_path_for_record(root, record) in public_paths for record in inventory.artifacts
           if record.classification is ArtifactClassification.CONTROLLED):
        raise GovernanceError('a controlled artifact aliases an inventoried public artifact')

    for record in inventory.artifacts:
        if record.classification is not ArtifactClassification.PUBLIC:
            continue
        path = _path_for_record(root, record)
        if not _matches_record(path, record):
            raise GovernanceError(
                f"public artifact changed before cleanup: {record.relative_path}"
            )
        preserved_public.append(record.artifact_id)

    updated: list[ArtifactRecord] = []
    for record in inventory.artifacts:
        if record.classification is ArtifactClassification.PUBLIC:
            updated.append(record)
            continue
        path = _path_for_record(root, record)
        if path.is_symlink():
            raise GovernanceError(
                f"refusing to follow controlled artifact symlink: {record.relative_path}"
            )
        if path.exists():
            if not path.is_file():
                raise GovernanceError(
                    f"controlled artifact is not a file: {record.relative_path}"
                )
            if not _matches_record(path, record):
                if record.operator_path is not None:
                    raise GovernanceError('operator-selected artifact changed before cleanup')
                identity_mismatches.append(record.artifact_id)
            path.unlink()
            deleted.append(record.artifact_id)
        else:
            already_absent.append(record.artifact_id)
        if path.exists() or path.is_symlink():
            raise GovernanceError(
                f"controlled artifact deletion could not be verified: {record.relative_path}"
            )
        updated.append(record if record.status == "deleted" else record.mark_deleted(timestamp))

    for record in inventory.artifacts:
        if record.classification is ArtifactClassification.PUBLIC:
            path = _path_for_record(root, record)
            if not _matches_record(path, record):
                raise GovernanceError(
                    f"public artifact changed during cleanup: {record.relative_path}"
                )

    updated_inventory = CaseInventory(inventory.case_id, tuple(updated), inventory.schema)
    report = CleanupReport(
        case_id=inventory.case_id,
        completed_on=timestamp,
        deleted_artifacts=tuple(deleted),
        already_absent_artifacts=tuple(already_absent),
        identity_mismatches_deleted=tuple(identity_mismatches),
        preserved_public_artifacts=tuple(preserved_public),
        verified=True,
    )
    return updated_inventory, report


def write_cleanup_report(report: CleanupReport, path: str | Path) -> None:
    destination = Path(path)
    destination.parent.mkdir(parents=True, exist_ok=True, mode=0o700)
    content = (
        json.dumps(asdict(report), sort_keys=True, separators=(",", ":")) + "\n"
    ).encode()
    partial = destination.parent / f".{destination.name}.{uuid.uuid4().hex}.partial"
    try:
        with partial.open("xb") as output:
            os.chmod(partial, 0o600)
            output.write(content)
            output.flush()
            os.fsync(output.fileno())
        os.replace(partial, destination)
        os.chmod(destination, 0o600)
    finally:
        partial.unlink(missing_ok=True)
