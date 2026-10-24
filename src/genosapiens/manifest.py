"""Immutable public-resource manifests with output verification."""

from __future__ import annotations

from dataclasses import asdict, dataclass
from datetime import date
import hashlib
import json
import os
from pathlib import Path, PurePosixPath
import re
from typing import Any, Mapping
import uuid

from genosapiens.governance import (
    ArtifactClassification,
    GovernanceError,
    validate_relative_classification,
)


MANIFEST_SCHEMA = "genosapiens.public-resource"
SHA256 = re.compile(r"^[0-9a-f]{64}$")
DIGEST = re.compile(r"^sha256:[0-9a-f]{64}$")
IDENTIFIER = re.compile(r"^[a-z0-9][a-z0-9._-]*$")


class ManifestError(RuntimeError):
    """Raised when resource provenance or output identity is invalid."""


class ManifestMismatch(ManifestError):
    """Raised when an immutable manifest or completed output differs."""


@dataclass(frozen=True, slots=True)
class SourceIdentity:
    url: str
    sha256: str | None = None
    object_generation: str | None = None

    def __post_init__(self) -> None:
        if not self.url.startswith("https://"):
            raise ManifestError("resource source must use HTTPS")
        if self.sha256 is None and self.object_generation is None:
            raise ManifestError("source requires a SHA-256 or object generation")
        if self.sha256 is not None and not SHA256.fullmatch(self.sha256):
            raise ManifestError("source SHA-256 is invalid")
        if self.object_generation is not None and not self.object_generation:
            raise ManifestError("source object generation cannot be empty")


@dataclass(frozen=True, slots=True)
class Transformation:
    schema: str
    parameters_sha256: str

    def __post_init__(self) -> None:
        if not self.schema:
            raise ManifestError("transformation schema is required")
        if not SHA256.fullmatch(self.parameters_sha256):
            raise ManifestError("transformation parameters SHA-256 is invalid")


@dataclass(frozen=True, slots=True)
class RetentionPolicy:
    policy: str
    delete_after: str | None

    def __post_init__(self) -> None:
        if self.policy not in {"durable", "transient-until-validation"}:
            raise ManifestError("unsupported retention policy")
        if self.delete_after is not None:
            try:
                date.fromisoformat(self.delete_after)
            except ValueError as error:
                raise ManifestError("delete_after must be an ISO date") from error


@dataclass(frozen=True, slots=True)
class ResourceOutput:
    path: str
    size_bytes: int
    sha256: str
    retention: RetentionPolicy

    def __post_init__(self) -> None:
        path = PurePosixPath(self.path)
        if path.is_absolute() or ".." in path.parts or path == PurePosixPath("."):
            raise ManifestError("output path must be a safe relative path")
        try:
            validate_relative_classification(
                self.path, ArtifactClassification.PUBLIC
            )
        except GovernanceError as error:
            raise ManifestError(str(error)) from error
        if not isinstance(self.size_bytes, int) or self.size_bytes < 0:
            raise ManifestError("output size_bytes must be a non-negative integer")
        if not SHA256.fullmatch(self.sha256):
            raise ManifestError("output SHA-256 is invalid")


@dataclass(frozen=True, slots=True)
class PublicResourceManifest:
    resource_id: str
    release: str
    source: SourceIdentity
    license: str
    assembly: str
    retrieved_on: str
    transformation: Transformation
    image_digest: str
    outputs: tuple[ResourceOutput, ...]
    schema: str = MANIFEST_SCHEMA
    classification: str = "public"

    def __post_init__(self) -> None:
        if self.schema != MANIFEST_SCHEMA:
            raise ManifestError(f"manifest schema must be {MANIFEST_SCHEMA!r}")
        if self.classification != "public":
            raise ManifestError("public resource manifests require public classification")
        if not IDENTIFIER.fullmatch(self.resource_id):
            raise ManifestError("resource_id is invalid")
        for name, value in (
            ("release", self.release),
            ("license", self.license),
            ("assembly", self.assembly),
        ):
            if not value:
                raise ManifestError(f"{name} is required")
        try:
            date.fromisoformat(self.retrieved_on)
        except ValueError as error:
            raise ManifestError("retrieved_on must be an ISO date") from error
        if not DIGEST.fullmatch(self.image_digest):
            raise ManifestError("image_digest must be a SHA-256 digest")
        if not self.outputs:
            raise ManifestError("manifest requires at least one output")
        paths = [output.path for output in self.outputs]
        if len(paths) != len(set(paths)):
            raise ManifestError("manifest output paths must be unique")

    @classmethod
    def from_mapping(cls, data: Mapping[str, Any]) -> "PublicResourceManifest":
        try:
            source = SourceIdentity(**data["source"])
            transformation = Transformation(**data["transformation"])
            outputs = tuple(
                ResourceOutput(
                    path=item["path"],
                    size_bytes=item["size_bytes"],
                    sha256=item["sha256"],
                    retention=RetentionPolicy(**item["retention"]),
                )
                for item in data["outputs"]
            )
            return cls(
                schema=data["schema"],
                classification=data["classification"],
                resource_id=data["resource_id"],
                release=data["release"],
                source=source,
                license=data["license"],
                assembly=data["assembly"],
                retrieved_on=data["retrieved_on"],
                transformation=transformation,
                image_digest=data["image_digest"],
                outputs=outputs,
            )
        except (KeyError, TypeError, ValueError) as error:
            raise ManifestError(f"invalid resource manifest: {error}") from error


@dataclass(frozen=True, slots=True)
class ManifestPublication:
    path: Path
    manifest_sha256: str
    reused: bool


def load_public_resource_manifest(path: str | Path) -> PublicResourceManifest:
    """Load the shared public-resource contract without converting its identity."""
    source = Path(path)
    try:
        data = json.loads(source.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as error:
        raise ManifestError(f"cannot load public resource manifest: {error}") from error
    if not isinstance(data, Mapping):
        raise ManifestError("public resource manifest root must be an object")
    return PublicResourceManifest.from_mapping(data)


def _canonical_bytes(manifest: PublicResourceManifest) -> bytes:
    return (
        json.dumps(asdict(manifest), sort_keys=True, separators=(",", ":")) + "\n"
    ).encode()


def public_resource_manifest_from_bytes(content: bytes) -> PublicResourceManifest:
    """Validate exact native provenance metadata without opening its outputs."""
    from genosapiens.evidence_model import _typed
    from genosapiens.evidence_serialization import _loads

    try:
        value = _typed(_loads(content), PublicResourceManifest, decode=True)
        if _canonical_bytes(value) != content:
            raise ValueError
        return value
    except (ValueError, TypeError, RuntimeError, AttributeError):
        raise ManifestError("invalid canonical public resource metadata") from None


def validate_public_evidence_manifest_bytes(content: bytes) -> None:
    """Check current native resource/bundle metadata, never raw source rows."""
    from genosapiens.evidence_serialization import _loads
    from genosapiens.rare_disease import BUNDLE_MANIFEST_SCHEMA, validate_bundle_manifest_metadata

    try:
        data = _loads(content)
        if isinstance(data, Mapping) and data.get("schema") == BUNDLE_MANIFEST_SCHEMA:
            validate_bundle_manifest_metadata(data)
        else:
            public_resource_manifest_from_bytes(content)
    except (ValueError, TypeError, RuntimeError, AttributeError):
        raise ManifestError("invalid current public evidence metadata") from None


def _file_identity(path: Path) -> tuple[str, int]:
    digest = hashlib.sha256()
    size = 0
    with path.open("rb") as source:
        while chunk := source.read(1024 * 1024):
            digest.update(chunk)
            size += len(chunk)
    return digest.hexdigest(), size


def verify_outputs(manifest: PublicResourceManifest, reference_root: str | Path) -> None:
    root = Path(reference_root).resolve()
    for output in manifest.outputs:
        path = root / output.path
        if not path.resolve().is_relative_to(root) or not path.is_file():
            raise ManifestMismatch(f"manifest output is missing or invalid: {output.path}")
        digest, size = _file_identity(path)
        if size != output.size_bytes or digest != output.sha256:
            raise ManifestMismatch(f"manifest output checksum mismatch: {output.path}")


def _reuse_manifest(
    path: Path,
    content: bytes,
    manifest_digest: str,
    root: Path,
) -> ManifestPublication:
    if path.is_symlink() or not path.is_file():
        raise ManifestMismatch(f"manifest path is not a regular file: {path}")
    try:
        existing_data = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as error:
        raise ManifestMismatch(f"cannot read completed manifest: {error}") from error
    existing = PublicResourceManifest.from_mapping(existing_data)
    if _canonical_bytes(existing) != content:
        raise ManifestMismatch("completed resource manifest is immutable")
    verify_outputs(existing, root)
    return ManifestPublication(path, manifest_digest, reused=True)


def publish_manifest(
    manifest: PublicResourceManifest,
    manifest_path: str | Path,
    reference_root: str | Path,
) -> ManifestPublication:
    path = Path(manifest_path)
    root = Path(reference_root).resolve()
    if not path.is_absolute() or not path.parent.resolve().is_relative_to(root):
        raise ManifestError("manifest path must be under reference_root")
    content = _canonical_bytes(manifest)
    manifest_digest = hashlib.sha256(content).hexdigest()
    if path.exists() or path.is_symlink():
        return _reuse_manifest(path, content, manifest_digest, root)

    verify_outputs(manifest, root)
    path.parent.mkdir(parents=True, exist_ok=True)
    partial_root = path.parent / ".partial"
    if partial_root.is_symlink():
        raise ManifestMismatch(f"manifest partial path is a symlink: {partial_root}")
    partial_root.mkdir(mode=0o700, exist_ok=True)
    if not partial_root.resolve().is_relative_to(root):
        raise ManifestMismatch("manifest partial path must be under reference_root")
    partial = partial_root / f"{path.name}.{uuid.uuid4().hex}.partial"
    try:
        with partial.open("xb") as output:
            os.chmod(partial, 0o600)
            output.write(content)
            output.flush()
            os.fsync(output.fileno())
        os.chmod(partial, 0o444)
        try:
            os.link(partial, path)
        except FileExistsError:
            return _reuse_manifest(path, content, manifest_digest, root)
        except OSError as error:
            raise ManifestMismatch(f"cannot publish manifest: {error}") from error
        directory_fd = os.open(path.parent, os.O_RDONLY)
        try:
            os.fsync(directory_fd)
        finally:
            os.close(directory_fd)
    finally:
        partial.unlink(missing_ok=True)
    return ManifestPublication(path, manifest_digest, reused=False)
