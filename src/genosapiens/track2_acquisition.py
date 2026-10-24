"""Atomic, immutable acquisition of provenance-approved Track 2 public queries."""

from __future__ import annotations

from dataclasses import dataclass
from datetime import date
import hashlib
import json
import os
from pathlib import Path
import re
from typing import Any, Callable, Mapping
import uuid

from genosapiens.download import DownloadRequest, DownloadResult, download_atomic
from genosapiens.manifest import (
    ManifestPublication,
    PublicResourceManifest,
    ResourceOutput,
    RetentionPolicy,
    SourceIdentity,
    Transformation,
    publish_manifest,
)
from genosapiens.track2_public_query import (
    AuditLogger,
    PublicSeedRegistry,
    validate_public_query,
)


SHA256 = re.compile(r"^[0-9a-f]{64}$")
DIGEST = re.compile(r"^sha256:[0-9a-f]{64}$")


class Track2AcquisitionError(RuntimeError):
    """Raised when a public snapshot cannot be acquired or published safely."""


@dataclass(frozen=True, slots=True)
class PublicAcquisitionPlan:
    response_url: str
    expected_size_bytes: int
    expected_sha256: str
    license: str
    redistribution: str
    transformation_schema: str
    transformation_parameters_sha256: str
    accessed_on: str
    image_digest: str

    def __post_init__(self) -> None:
        try:
            date.fromisoformat(self.accessed_on)
        except (TypeError, ValueError) as error:
            raise Track2AcquisitionError("accessed_on must be an ISO date") from error
        if (
            not self.response_url.startswith("https://")
            or not isinstance(self.expected_size_bytes, int)
            or isinstance(self.expected_size_bytes, bool)
            or self.expected_size_bytes <= 0
            or not SHA256.fullmatch(self.expected_sha256)
            or not self.license
            or not self.transformation_schema
            or not SHA256.fullmatch(self.transformation_parameters_sha256)
            or not DIGEST.fullmatch(self.image_digest)
        ):
            raise Track2AcquisitionError("public acquisition plan is incomplete")
        if self.redistribution != "raw_publication_allowed":
            raise Track2AcquisitionError(
                "non-redistributable raw content cannot be published"
            )


@dataclass(frozen=True, slots=True)
class PublicAcquisitionResult:
    query_id: str
    response: DownloadResult
    metadata_path: Path
    metadata_sha256: str
    manifest: ManifestPublication


def _canonical_bytes(value: Mapping[str, Any]) -> bytes:
    return (
        json.dumps(value, sort_keys=True, separators=(",", ":"), ensure_ascii=True)
        + "\n"
    ).encode()


def publish_immutable_public_bytes(path: Path, content: bytes) -> str:
    expected = hashlib.sha256(content).hexdigest()
    if path.exists() or path.is_symlink():
        if path.is_symlink() or not path.is_file() or path.read_bytes() != content:
            raise Track2AcquisitionError("completed acquisition metadata is immutable")
        return expected
    path.parent.mkdir(parents=True, exist_ok=True)
    partial_root = path.parent / ".partial"
    if partial_root.is_symlink():
        raise Track2AcquisitionError("acquisition partial path cannot be a symlink")
    partial_root.mkdir(mode=0o700, exist_ok=True)
    partial = partial_root / f"{path.name}.{uuid.uuid4().hex}.partial"
    try:
        with partial.open("xb") as output:
            os.chmod(partial, 0o600)
            output.write(content)
            output.flush()
            os.fsync(output.fileno())
        try:
            os.link(partial, path)
        except FileExistsError:
            if path.is_symlink() or not path.is_file() or path.read_bytes() != content:
                raise Track2AcquisitionError("completed acquisition metadata is immutable")
        os.chmod(path, 0o444)
        directory_fd = os.open(path.parent, os.O_RDONLY)
        try:
            os.fsync(directory_fd)
        finally:
            os.close(directory_fd)
    finally:
        partial.unlink(missing_ok=True)
    return expected


def acquire_public_query_snapshot(
    query: Mapping[str, Any],
    registry: PublicSeedRegistry,
    expansion_policy: Mapping[str, Any],
    plan: PublicAcquisitionPlan,
    reference_root: str | Path,
    *,
    logger: AuditLogger,
    opener: Callable[..., Any] | None = None,
    sleeper: Callable[[float], None] | None = None,
    public_scope: Mapping[str, Any] | None = None,
    route: str | None = None,
    domain: str | None = None,
) -> PublicAcquisitionResult:
    """Publish response and metadata only after boundary and checksum validation."""
    validate_public_query(query, registry, expansion_policy, logger,
                          public_scope=public_scope, route=route, domain=domain)
    root = Path(reference_root).resolve()
    query_id = str(query["query_id"])
    suffix = query_id.removeprefix("public-query:")
    if not re.fullmatch(r"[0-9a-f]{24}", suffix):
        raise Track2AcquisitionError("query ID is not a stable public-query ID")
    response_path = root / f"releases/track2/acquisition/{suffix}/response.bin"
    download_kwargs: dict[str, Any] = {}
    if opener is not None:
        download_kwargs["opener"] = opener
    if sleeper is not None:
        download_kwargs["sleeper"] = sleeper
    response = download_atomic(
        DownloadRequest(
            url=plan.response_url,
            destination=response_path,
            expected_size_bytes=plan.expected_size_bytes,
            expected_sha256=plan.expected_sha256,
        ),
        root,
        **download_kwargs,
    )
    provenance = query["provenance"]
    metadata = {
        "schema": "genosapiens.track2-public-acquisition/v1",
        "query": {
            "query_id": query_id,
            "adapter": query["adapter"],
            "parameters": query["parameters"],
        },
        "seed": {
            "entity_id": query["seed_id"],
            "entity_type": query["entity_type"],
            "manifest_sha256": provenance["manifest_sha256"],
            "provenance_token": provenance["token"],
            "source_reference_id": provenance["source_reference_id"],
        },
        "expansion": query["parameters"],
        "response": {
            "url": plan.response_url,
            "sha256": response.sha256,
            "size_bytes": response.size_bytes,
            "source_generation": response.source_generation,
        },
        "license": plan.license,
        "redistribution": plan.redistribution,
        "transformation": {
            "schema": plan.transformation_schema,
            "parameters_sha256": plan.transformation_parameters_sha256,
        },
        "accessed_on": plan.accessed_on,
    }
    metadata_content = _canonical_bytes(metadata)
    metadata_path = root / f"derived/track2/acquisition/{suffix}/metadata.json"
    metadata_sha256 = publish_immutable_public_bytes(metadata_path, metadata_content)
    manifest = PublicResourceManifest(
        resource_id="track2-public-acquisition",
        release=f"{plan.accessed_on}-{suffix}",
        source=SourceIdentity(plan.response_url, sha256=response.sha256),
        license=plan.license,
        assembly="not-applicable",
        retrieved_on=plan.accessed_on,
        transformation=Transformation(
            plan.transformation_schema, plan.transformation_parameters_sha256
        ),
        image_digest=plan.image_digest,
        outputs=(
            ResourceOutput(
                path=response_path.relative_to(root).as_posix(),
                size_bytes=response.size_bytes,
                sha256=response.sha256,
                retention=RetentionPolicy("durable", None),
            ),
            ResourceOutput(
                path=metadata_path.relative_to(root).as_posix(),
                size_bytes=len(metadata_content),
                sha256=metadata_sha256,
                retention=RetentionPolicy("durable", None),
            ),
        ),
    )
    publication = publish_manifest(
        manifest,
        root / f"manifests/track2/acquisition/{suffix}.json",
        root,
    )
    return PublicAcquisitionResult(
        query_id=query_id,
        response=response,
        metadata_path=metadata_path,
        metadata_sha256=metadata_sha256,
        manifest=publication,
    )
