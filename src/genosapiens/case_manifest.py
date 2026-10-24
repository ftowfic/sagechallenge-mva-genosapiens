"""Deterministic controlled-case manifests and preparation reuse tracking."""

from __future__ import annotations

from dataclasses import asdict, dataclass, replace
from enum import StrEnum
import hashlib
import json
import os
from pathlib import Path
import re
import stat
from typing import Any, Mapping
import uuid

from genosapiens.bounded_json import BoundedJSONError, load_bounded_json
from genosapiens.governance import (
    ArtifactClassification,
    ArtifactRecord,
    GovernanceError,
    validate_relative_classification,
)
from genosapiens.llm_guard import LLMPolicy


CASE_MANIFEST_SCHEMA = "genosapiens.case-manifest"
SHA256 = re.compile(r"^[0-9a-f]{64}$")
IDENTIFIER = re.compile(r"^[A-Za-z0-9][A-Za-z0-9._:-]*$")
SAMPLE_IDENTIFIER = re.compile(r"^[^\x00-\x20\x7f]+$")


class CaseManifestError(RuntimeError):
    """Raised when a controlled case manifest is incomplete or inconsistent."""


class CaseFileClass(StrEnum):
    HPO_LIST = "hpo-list"
    VCF = "vcf"
    VCF_INDEX = "vcf-index"
    FASTQ_R1 = "fastq-r1"
    FASTQ_R2 = "fastq-r2"


class RetentionStatus(StrEnum):
    ACTIVE = "active"
    DELETION_PENDING = "deletion-pending"
    DELETED = "deleted"


@dataclass(frozen=True, slots=True)
class CaseInput:
    input_id: str
    inventory_artifact_id: str
    role: str
    file_class: CaseFileClass
    sha256: str
    size_bytes: int
    classification: ArtifactClassification = ArtifactClassification.CONTROLLED
    llm_policy: LLMPolicy = LLMPolicy.NO_LLM
    retention: str = "delete"

    def __post_init__(self) -> None:
        if (
            IDENTIFIER.fullmatch(self.input_id) is None
            or IDENTIFIER.fullmatch(self.inventory_artifact_id) is None
            or IDENTIFIER.fullmatch(self.role) is None
            or SHA256.fullmatch(self.sha256) is None
            or not isinstance(self.size_bytes, int)
            or isinstance(self.size_bytes, bool)
            or self.size_bytes < 0
            or self.classification is not ArtifactClassification.CONTROLLED
            or self.llm_policy is not LLMPolicy.NO_LLM
            or self.retention != "delete"
            or self.input_id != _input_id(self.role, self.file_class, self.sha256)
        ):
            raise CaseManifestError("case input binding is invalid")


@dataclass(frozen=True, slots=True)
class PhenotypeSource:
    input_id: str

    def __post_init__(self) -> None:
        if IDENTIFIER.fullmatch(self.input_id) is None:
            raise CaseManifestError("phenotype source binding is invalid")


@dataclass(frozen=True, slots=True)
class ConfigurationBinding:
    schema: str
    sha256: str

    def __post_init__(self) -> None:
        if (
            not isinstance(self.schema, str)
            or not self.schema
            or any(ord(character) < 32 for character in self.schema)
            or SHA256.fullmatch(self.sha256) is None
        ):
            raise CaseManifestError("configuration binding is invalid")


@dataclass(frozen=True, slots=True)
class ReferenceBundleBinding:
    bundle_id: str
    manifest_sha256: str

    def __post_init__(self) -> None:
        if (
            IDENTIFIER.fullmatch(self.bundle_id) is None
            or SHA256.fullmatch(self.manifest_sha256) is None
        ):
            raise CaseManifestError("reference bundle binding is invalid")


@dataclass(frozen=True, slots=True)
class CaseManifest:
    case_id: str
    sample_id: str
    inputs: tuple[CaseInput, ...]
    available_file_classes: tuple[CaseFileClass, ...]
    phenotype_source: PhenotypeSource
    configuration: ConfigurationBinding
    reference_bundle: ReferenceBundleBinding
    retention_status: RetentionStatus
    reuse_count: int
    reused: bool
    classification: ArtifactClassification
    llm_policy: LLMPolicy
    content_sha256: str
    schema: str = CASE_MANIFEST_SCHEMA

    def __post_init__(self) -> None:
        input_ids = [item.input_id for item in self.inputs]
        artifact_ids = [item.inventory_artifact_id for item in self.inputs]
        roles = [item.role for item in self.inputs]
        expected_classes = tuple(sorted({item.file_class for item in self.inputs}))
        phenotype_inputs = {
            item.input_id: item.file_class
            for item in self.inputs
            if item.file_class is CaseFileClass.HPO_LIST
        }
        if (
            self.schema != CASE_MANIFEST_SCHEMA
            or re.fullmatch(r"case-[0-9a-f]{24}", self.case_id) is None
            or SAMPLE_IDENTIFIER.fullmatch(self.sample_id) is None
            or not self.inputs
            or self.inputs != tuple(sorted(self.inputs, key=lambda item: item.role))
            or len(input_ids) != len(set(input_ids))
            or len(artifact_ids) != len(set(artifact_ids))
            or len(roles) != len(set(roles))
            or CaseFileClass.VCF not in expected_classes
            or self.available_file_classes != expected_classes
            or len(phenotype_inputs) != 1
            or phenotype_inputs.get(self.phenotype_source.input_id)
            is not CaseFileClass.HPO_LIST
            or not isinstance(self.reuse_count, int)
            or isinstance(self.reuse_count, bool)
            or self.reuse_count < 0
            or self.reused != (self.reuse_count > 0)
            or self.classification is not ArtifactClassification.CONTROLLED
            or self.llm_policy is not LLMPolicy.NO_LLM
            or SHA256.fullmatch(self.content_sha256) is None
            or self.case_id != _case_id(_identity_semantic(self))
            or self.content_sha256 != _manifest_sha256(self)
        ):
            raise CaseManifestError("case manifest is invalid")


@dataclass(frozen=True, slots=True)
class CasePreparationResult:
    manifest: CaseManifest
    reused: bool


def _canonical_bytes(value: Mapping[str, Any]) -> bytes:
    return (
        json.dumps(value, sort_keys=True, separators=(",", ":"), ensure_ascii=True)
        + "\n"
    ).encode()


def _input_id(role: str, file_class: CaseFileClass, sha256: str) -> str:
    digest = hashlib.sha256(
        f"{role}\0{file_class.value}\0{sha256}".encode()
    ).hexdigest()[:24]
    return f"input-{digest}"


def case_input_from_artifact(
    artifact: ArtifactRecord,
    *,
    role: str,
    file_class: CaseFileClass,
) -> CaseInput:
    if (
        artifact.classification is not ArtifactClassification.CONTROLLED
        or artifact.retention != "delete"
        or artifact.status != "present"
    ):
        raise CaseManifestError("case inputs must be present controlled artifacts")
    return CaseInput(
        input_id=_input_id(role, file_class, artifact.sha256),
        inventory_artifact_id=artifact.artifact_id,
        role=role,
        file_class=file_class,
        sha256=artifact.sha256,
        size_bytes=artifact.size_bytes,
    )


def _identity_semantic(manifest: CaseManifest) -> dict[str, Any]:
    return {
        "sample_id": manifest.sample_id,
        "inputs": [
            {
                key: value
                for key, value in asdict(item).items()
                if key != "inventory_artifact_id"
            }
            for item in manifest.inputs
        ],
        "available_file_classes": list(manifest.available_file_classes),
        "phenotype_source": asdict(manifest.phenotype_source),
        "configuration": asdict(manifest.configuration),
        "reference_bundle": asdict(manifest.reference_bundle),
        "retention_status": manifest.retention_status,
        "classification": manifest.classification,
        "llm_policy": manifest.llm_policy,
        "schema": manifest.schema,
    }


def _case_id(semantic: Mapping[str, Any]) -> str:
    return f"case-{hashlib.sha256(_canonical_bytes(semantic)).hexdigest()[:24]}"


def _manifest_sha256(manifest: CaseManifest) -> str:
    value = asdict(manifest)
    value.pop("content_sha256")
    return hashlib.sha256(_canonical_bytes(value)).hexdigest()


def create_case_manifest(
    *,
    sample_id: str,
    inputs: tuple[CaseInput, ...],
    configuration: ConfigurationBinding,
    reference_bundle: ReferenceBundleBinding,
    retention_status: RetentionStatus = RetentionStatus.ACTIVE,
    reuse_count: int = 0,
) -> CaseManifest:
    ordered_inputs = tuple(sorted(inputs, key=lambda item: item.role))
    phenotype_candidates = tuple(
        item
        for item in ordered_inputs
        if item.file_class is CaseFileClass.HPO_LIST
    )
    if len(phenotype_candidates) != 1:
        raise CaseManifestError("case manifest requires exactly one phenotype source")
    phenotype_source = PhenotypeSource(input_id=phenotype_candidates[0].input_id)
    values: dict[str, Any] = {
        "sample_id": sample_id,
        "inputs": ordered_inputs,
        "available_file_classes": tuple(
            sorted({item.file_class for item in ordered_inputs})
        ),
        "phenotype_source": phenotype_source,
        "configuration": configuration,
        "reference_bundle": reference_bundle,
        "retention_status": retention_status,
        "reuse_count": reuse_count,
        "reused": reuse_count > 0,
        "classification": ArtifactClassification.CONTROLLED,
        "llm_policy": LLMPolicy.NO_LLM,
        "schema": CASE_MANIFEST_SCHEMA,
    }
    semantic = {
        key: value
        for key, value in values.items()
        if key not in {"reuse_count", "reused"}
    }
    semantic["inputs"] = [
        {
            key: value
            for key, value in asdict(item).items()
            if key != "inventory_artifact_id"
        }
        for item in ordered_inputs
    ]
    case_id = _case_id(_jsonable(semantic))
    content_sha256 = hashlib.sha256(
        _canonical_bytes(_jsonable({**values, "case_id": case_id}))
    ).hexdigest()
    return CaseManifest(case_id=case_id, content_sha256=content_sha256, **values)


def _jsonable(value: Any) -> Any:
    if hasattr(value, "__dataclass_fields__"):
        return {key: _jsonable(child) for key, child in asdict(value).items()}
    if isinstance(value, Mapping):
        return {key: _jsonable(child) for key, child in value.items()}
    if isinstance(value, (tuple, list)):
        return [_jsonable(child) for child in value]
    return value


def case_manifest_as_bytes(manifest: CaseManifest) -> bytes:
    return _canonical_bytes(asdict(manifest))


def load_case_manifest(path: str | Path) -> CaseManifest:
    try:
        data, content = load_bounded_json(path)
        expected = {
            "available_file_classes",
            "case_id",
            "classification",
            "configuration",
            "content_sha256",
            "inputs",
            "llm_policy",
            "phenotype_source",
            "reference_bundle",
            "retention_status",
            "reuse_count",
            "reused",
            "sample_id",
            "schema",
        }
        if (
            not isinstance(data, Mapping)
            or set(data) != expected
            or not isinstance(data["inputs"], list)
            or not isinstance(data["available_file_classes"], list)
            or not isinstance(data["phenotype_source"], Mapping)
            or not isinstance(data["configuration"], Mapping)
            or not isinstance(data["reference_bundle"], Mapping)
            or any(
                not isinstance(item, Mapping)
                or set(item)
                != {
                    "classification",
                    "file_class",
                    "input_id",
                    "inventory_artifact_id",
                    "llm_policy",
                    "retention",
                    "role",
                    "sha256",
                    "size_bytes",
                }
                for item in data["inputs"]
            )
            or set(data["phenotype_source"]) != {"input_id"}
            or set(data["configuration"]) != {"schema", "sha256"}
            or set(data["reference_bundle"]) != {"bundle_id", "manifest_sha256"}
        ):
            raise CaseManifestError("case manifest is invalid")
        manifest = CaseManifest(
            case_id=data["case_id"],
            sample_id=data["sample_id"],
            inputs=tuple(
                CaseInput(
                    input_id=item["input_id"],
                    inventory_artifact_id=item["inventory_artifact_id"],
                    role=item["role"],
                    file_class=CaseFileClass(item["file_class"]),
                    sha256=item["sha256"],
                    size_bytes=item["size_bytes"],
                    classification=ArtifactClassification(item["classification"]),
                    llm_policy=LLMPolicy(item["llm_policy"]),
                    retention=item["retention"],
                )
                for item in data["inputs"]
            ),
            available_file_classes=tuple(
                CaseFileClass(item) for item in data["available_file_classes"]
            ),
            phenotype_source=PhenotypeSource(
                input_id=data["phenotype_source"]["input_id"],
            ),
            configuration=ConfigurationBinding(**data["configuration"]),
            reference_bundle=ReferenceBundleBinding(**data["reference_bundle"]),
            retention_status=RetentionStatus(data["retention_status"]),
            reuse_count=data["reuse_count"],
            reused=data["reused"],
            classification=ArtifactClassification(data["classification"]),
            llm_policy=LLMPolicy(data["llm_policy"]),
            content_sha256=data["content_sha256"],
            schema=data["schema"],
        )
    except (
        BoundedJSONError,
        KeyError,
        TypeError,
        ValueError,
        GovernanceError,
    ) as error:
        raise CaseManifestError("case manifest is invalid") from error
    if case_manifest_as_bytes(manifest) != content:
        raise CaseManifestError("case manifest is not canonical")
    return manifest


def _write_case_manifest(manifest: CaseManifest, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True, mode=0o700)
    partial = path.parent / f".{path.name}.{uuid.uuid4().hex}.partial"
    try:
        with partial.open("xb") as output:
            os.chmod(partial, 0o600)
            output.write(case_manifest_as_bytes(manifest))
            output.flush()
            os.fsync(output.fileno())
        os.replace(partial, path)
        os.chmod(path, 0o600)
    finally:
        partial.unlink(missing_ok=True)


def prepare_case_manifest(
    manifest: CaseManifest,
    *,
    reference_root: str | Path,
    output_path: str | Path,
) -> CasePreparationResult:
    root = Path(reference_root).resolve(strict=True)
    output = Path(output_path)
    if output.is_symlink():
        raise CaseManifestError("case manifest output is invalid")
    try:
        relative = output.resolve(strict=False).relative_to(root).as_posix()
        validate_relative_classification(relative, ArtifactClassification.CONTROLLED)
    except (OSError, ValueError, GovernanceError) as error:
        raise CaseManifestError("case manifest output is invalid") from error
    if output.exists():
        try:
            if stat.S_IMODE(output.stat().st_mode) != 0o600:
                raise CaseManifestError("published case manifest is immutable")
        except OSError as error:
            raise CaseManifestError("case manifest output is invalid") from error
        existing = load_case_manifest(output)
        if _identity_semantic(existing) != _identity_semantic(manifest):
            raise CaseManifestError("existing case manifest has a different identity")
        updated = create_case_manifest(
            sample_id=existing.sample_id,
            inputs=manifest.inputs,
            configuration=existing.configuration,
            reference_bundle=existing.reference_bundle,
            retention_status=existing.retention_status,
            reuse_count=existing.reuse_count + 1,
        )
        _write_case_manifest(updated, output)
        return CasePreparationResult(updated, True)
    _write_case_manifest(manifest, output)
    return CasePreparationResult(manifest, False)
