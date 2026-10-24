"""Bootstrap open references for non-small-variant analysis lanes."""

from __future__ import annotations

from collections import Counter
from collections.abc import Callable, Mapping
from dataclasses import asdict, dataclass
from datetime import date
import gzip
import hashlib
import json
import math
import os
from pathlib import Path, PurePosixPath
import re
import shutil
from typing import Any, TextIO
import uuid

from genosapiens.download import DownloadRequest, DownloadResult, download_atomic
from genosapiens.manifest import (
    ManifestError,
    PublicResourceManifest,
    ResourceOutput,
    RetentionPolicy,
    SourceIdentity,
    Transformation,
    publish_manifest,
    verify_outputs,
)
from genosapiens.storage import DEFAULT_RESERVE_BYTES, CapacityPlan, require_capacity
from genosapiens.tabix import TabixError, TabixIndex


SCHEMA = "genosapiens.variant-class-resource-bundle"
TRANSFORMATION_SCHEMA = "genosapiens.source-preserving-class-resource-bootstrap"
BUNDLE_MANIFEST_SCHEMA = "genosapiens.variant-class-resource-bundle-manifest"
EXPECTED_LANES = {
    "sv_cnv",
    "mitochondrial",
    "mosaic",
    "regulatory",
    "repeat_expansion",
}
EXPECTED_LANE_COUNTS = {
    "sv_cnv": 2,
    "mitochondrial": 1,
    "mosaic": 1,
    "regulatory": 2,
    "repeat_expansion": 1,
}
VALIDATION_KINDS = {"vcf", "bed", "tsv", "repeat_catalog"}
IDENTIFIER = re.compile(r"^[a-z0-9][a-z0-9._-]*$")
SHA256 = re.compile(r"^[0-9a-f]{64}$")
MD5 = re.compile(r"^[0-9a-f]{32}$")
DIGEST = re.compile(r"^sha256:[0-9a-f]{64}$")
INFO_HEADER = re.compile(r"^##INFO=<ID=([^,>]+),")


class VariantClassResourceError(RuntimeError):
    """Raised when class-specific resource acquisition or validation fails."""


def _required_string(data: Mapping[str, Any], name: str, context: str) -> str:
    value = data.get(name)
    if not isinstance(value, str) or not value:
        raise VariantClassResourceError(f"{context}.{name} must be a non-empty string")
    return value


def _optional_string(
    data: Mapping[str, Any], name: str, context: str
) -> str | None:
    value = data.get(name)
    if value is not None and (not isinstance(value, str) or not value):
        raise VariantClassResourceError(
            f"{context}.{name} must be null or a non-empty string"
        )
    return value


def _positive_int(data: Mapping[str, Any], name: str, context: str) -> int:
    value = data.get(name)
    if not isinstance(value, int) or isinstance(value, bool) or value <= 0:
        raise VariantClassResourceError(f"{context}.{name} must be a positive integer")
    return value


def _nonnegative_int(data: Mapping[str, Any], name: str, context: str) -> int:
    value = data.get(name)
    if not isinstance(value, int) or isinstance(value, bool) or value < 0:
        raise VariantClassResourceError(
            f"{context}.{name} must be a non-negative integer"
        )
    return value


def _safe_path(value: str, context: str) -> PurePosixPath:
    path = PurePosixPath(value)
    if path.is_absolute() or path == PurePosixPath(".") or ".." in path.parts:
        raise VariantClassResourceError(f"{context} must be a safe relative path")
    return path


def _string_tuple(value: Any, context: str, *, allow_empty: bool = False) -> tuple[str, ...]:
    if (
        not isinstance(value, list)
        or (not allow_empty and not value)
        or any(not isinstance(item, str) or not item for item in value)
        or len(value) != len(set(value))
    ):
        raise VariantClassResourceError(f"{context} must be a unique string array")
    return tuple(value)


@dataclass(frozen=True, slots=True)
class ResourceArtifact:
    role: str
    url: str
    path: str
    size_bytes: int
    sha256: str | None
    expected_md5: str | None
    object_generation: str | None

    @classmethod
    def from_mapping(cls, data: Mapping[str, Any]) -> "ResourceArtifact":
        expected = {
            "role",
            "url",
            "path",
            "size_bytes",
            "sha256",
            "expected_md5",
            "object_generation",
        }
        if set(data) != expected:
            raise VariantClassResourceError("class resource artifact fields are invalid")
        artifact = cls(
            role=_required_string(data, "role", "artifact"),
            url=_required_string(data, "url", "artifact"),
            path=_required_string(data, "path", "artifact"),
            size_bytes=_positive_int(data, "size_bytes", "artifact"),
            sha256=_optional_string(data, "sha256", "artifact"),
            expected_md5=_optional_string(data, "expected_md5", "artifact"),
            object_generation=_optional_string(data, "object_generation", "artifact"),
        )
        if (
            not IDENTIFIER.fullmatch(artifact.role)
            or not artifact.url.startswith("https://")
            or (artifact.sha256 is not None and not SHA256.fullmatch(artifact.sha256))
            or (
                artifact.expected_md5 is not None
                and not MD5.fullmatch(artifact.expected_md5)
            )
            or (
                artifact.sha256 is None
                and artifact.expected_md5 is None
                and artifact.object_generation is None
            )
        ):
            raise VariantClassResourceError(
                f"invalid class resource artifact identity: {artifact.role}"
            )
        path = _safe_path(artifact.path, "class resource artifact path")
        if len(path.parts) < 3 or path.parts[0] != "releases":
            raise VariantClassResourceError(
                "class resource artifact path must be inside a release directory"
            )
        return artifact

    def source_identity(self) -> SourceIdentity:
        if self.sha256 is None and self.object_generation is None:
            raise VariantClassResourceError(
                "primary artifact requires a SHA-256 or object generation"
            )
        return SourceIdentity(
            self.url,
            sha256=self.sha256,
            object_generation=self.object_generation,
        )


def _validate_sentinel_mapping(
    validation: Mapping[str, Any], kind: str
) -> Mapping[str, Any]:
    sentinel = validation.get("sentinel")
    if not isinstance(sentinel, Mapping):
        raise VariantClassResourceError(f"{kind} validation sentinel must be an object")
    if kind == "vcf":
        expected = {
            "chromosome",
            "position",
            "id",
            "reference",
            "alternate",
            "filter",
            "info",
        }
        if set(sentinel) != expected or not isinstance(sentinel.get("info"), Mapping):
            raise VariantClassResourceError("VCF sentinel fields are invalid")
        _positive_int(sentinel, "position", "VCF sentinel")
        for name in expected - {"position", "info"}:
            _required_string(sentinel, name, "VCF sentinel")
        info = sentinel["info"]
        if not info or any(
            not isinstance(key, str)
            or not key
            or not isinstance(value, str)
            or not value
            for key, value in info.items()
        ):
            raise VariantClassResourceError("VCF sentinel INFO values are invalid")
    elif kind in {"bed", "tsv"}:
        fields = sentinel.get("fields")
        if set(sentinel) != {"fields"} or not isinstance(fields, list) or not fields:
            raise VariantClassResourceError(f"{kind} sentinel fields are invalid")
        if any(not isinstance(field, str) for field in fields):
            raise VariantClassResourceError(f"{kind} sentinel fields must be strings")
    else:
        expected = {
            "locus_id",
            "main_reference_region",
            "reference_regions",
            "repeat_unit",
            "locus_structure",
            "variant_types",
            "inheritance",
            "disease_symbol",
            "normal_max",
            "pathogenic_min",
            "intermediate_range",
        }
        if set(sentinel) != expected:
            raise VariantClassResourceError("repeat catalog sentinel fields are invalid")
        for name in expected - {
            "normal_max",
            "pathogenic_min",
            "reference_regions",
            "variant_types",
        }:
            _required_string(sentinel, name, "repeat sentinel")
        reference_regions = _string_tuple(
            sentinel.get("reference_regions"), "repeat sentinel reference_regions"
        )
        variant_types_value = sentinel.get("variant_types")
        if (
            not isinstance(variant_types_value, list)
            or not variant_types_value
            or any(
                not isinstance(item, str) or not item for item in variant_types_value
            )
        ):
            raise VariantClassResourceError(
                "repeat sentinel variant_types must be a string array"
            )
        variant_types = tuple(variant_types_value)
        if (
            len(reference_regions) != len(variant_types)
            or sentinel["main_reference_region"] not in reference_regions
        ):
            raise VariantClassResourceError(
                "repeat sentinel regions and variant types are inconsistent"
            )
        _positive_int(sentinel, "normal_max", "repeat sentinel")
        _positive_int(sentinel, "pathogenic_min", "repeat sentinel")
    return sentinel


def _validate_validation_mapping(
    data: Mapping[str, Any], artifact_roles: set[str]
) -> Mapping[str, Any]:
    kind = _required_string(data, "kind", "validation")
    if kind not in VALIDATION_KINDS:
        raise VariantClassResourceError(f"unsupported validation kind: {kind}")
    if kind == "vcf":
        expected = {
            "kind",
            "artifact_role",
            "index_artifact_role",
            "required_info_fields",
            "sentinel",
        }
        index_role = data.get("index_artifact_role")
        if index_role is not None and (
            not isinstance(index_role, str) or index_role not in artifact_roles
        ):
            raise VariantClassResourceError("VCF index artifact role is invalid")
        required = _string_tuple(data.get("required_info_fields"), "VCF INFO fields")
        sentinel = _validate_sentinel_mapping(data, kind)
        if set(sentinel["info"]) - set(required):
            raise VariantClassResourceError("VCF sentinel INFO fields are not required")
    elif kind in {"bed", "tsv"}:
        expected = {"kind", "artifact_role", "column_count", "sentinel"}
        if "interval_convention" in data or "interpretation" in data:
            expected.update({"interval_convention", "interpretation"})
            _required_string(data, "interval_convention", f"{kind} validation")
            _required_string(data, "interpretation", f"{kind} validation")
        count = _positive_int(data, "column_count", f"{kind} validation")
        sentinel = _validate_sentinel_mapping(data, kind)
        if len(sentinel["fields"]) != count:
            raise VariantClassResourceError(f"{kind} sentinel width is invalid")
    else:
        expected = {"kind", "artifact_role", "expected_record_count", "sentinel"}
        _positive_int(data, "expected_record_count", "repeat catalog validation")
        _validate_sentinel_mapping(data, kind)
    if set(data) != expected:
        raise VariantClassResourceError(f"{kind} validation fields are invalid")
    artifact_role = _required_string(data, "artifact_role", "validation")
    if artifact_role not in artifact_roles:
        raise VariantClassResourceError("validation artifact role is not declared")
    return dict(data)


@dataclass(frozen=True, slots=True)
class ClassResource:
    lane: str
    resource_id: str
    release: str
    assembly: str
    sequence_identity: str
    license: str
    license_url: str
    primary_artifact_role: str
    artifacts: tuple[ResourceArtifact, ...]
    validation: Mapping[str, Any]

    @classmethod
    def from_mapping(cls, data: Mapping[str, Any]) -> "ClassResource":
        expected = {
            "lane",
            "resource_id",
            "release",
            "assembly",
            "sequence_identity",
            "license",
            "license_url",
            "primary_artifact_role",
            "artifacts",
            "validation",
        }
        artifacts_data = data.get("artifacts")
        validation_data = data.get("validation")
        if (
            set(data) != expected
            or not isinstance(artifacts_data, list)
            or not artifacts_data
            or not all(isinstance(item, Mapping) for item in artifacts_data)
            or not isinstance(validation_data, Mapping)
        ):
            raise VariantClassResourceError("class resource fields are invalid")
        artifacts = tuple(ResourceArtifact.from_mapping(item) for item in artifacts_data)
        roles = [artifact.role for artifact in artifacts]
        paths = [artifact.path for artifact in artifacts]
        primary_role = _required_string(
            data, "primary_artifact_role", "class resource"
        )
        resource = cls(
            lane=_required_string(data, "lane", "class resource"),
            resource_id=_required_string(data, "resource_id", "class resource"),
            release=_required_string(data, "release", "class resource"),
            assembly=_required_string(data, "assembly", "class resource"),
            sequence_identity=_required_string(
                data, "sequence_identity", "class resource"
            ),
            license=_required_string(data, "license", "class resource"),
            license_url=_required_string(data, "license_url", "class resource"),
            primary_artifact_role=primary_role,
            artifacts=artifacts,
            validation=_validate_validation_mapping(validation_data, set(roles)),
        )
        parents = {str(PurePosixPath(path).parent) for path in paths}
        if (
            resource.lane not in EXPECTED_LANES
            or not IDENTIFIER.fullmatch(resource.resource_id)
            or not resource.license_url.startswith("https://")
            or len(roles) != len(set(roles))
            or len(paths) != len(set(paths))
            or len(parents) != 1
            or resource.assembly != "GRCh38"
            or primary_role not in roles
            or not _valid_lane_contract(resource)
        ):
            raise VariantClassResourceError(
                f"class resource identities are inconsistent: {resource.resource_id}"
            )
        resource.artifact(primary_role).source_identity()
        return resource

    def artifact(self, role: str) -> ResourceArtifact:
        for artifact in self.artifacts:
            if artifact.role == role:
                return artifact
        raise VariantClassResourceError(
            f"resource {self.resource_id} lacks artifact role {role}"
        )

    @property
    def release_root(self) -> PurePosixPath:
        return PurePosixPath(self.artifacts[0].path).parent


def _valid_lane_contract(resource: ClassResource) -> bool:
    kind = resource.validation["kind"]
    required = set(resource.validation.get("required_info_fields", ()))
    if resource.lane == "sv_cnv":
        sv_fields = {"END", "SVTYPE", "AC", "AN", "AF"}
        cnv_fields = {"END", "SVTYPE", "Genes", "SC", "SN", "SF"}
        is_sv = sv_fields.issubset(required)
        is_cnv = cnv_fields.issubset(required)
        return (
            kind == "vcf"
            and is_sv != is_cnv
            and (not is_sv or resource.validation["index_artifact_role"] is not None)
        )
    if resource.lane == "mitochondrial":
        mitochondrial_fields = {
            "AN",
            "AC_hom",
            "AC_het",
            "AF_hom",
            "AF_het",
            "max_hl",
        }
        return (
            kind == "vcf"
            and mitochondrial_fields.issubset(required)
            and resource.validation["index_artifact_role"] is not None
        )
    if resource.lane == "mosaic":
        return (
            kind == "bed"
            and resource.validation.get("interval_convention")
            == "zero-based-half-open"
            and resource.validation.get("interpretation") == "artifact_warning"
        )
    if resource.lane == "regulatory":
        return kind in {"bed", "tsv"}
    return kind == "repeat_catalog"


@dataclass(frozen=True, slots=True)
class RegulatoryJoin:
    registry_resource_id: str
    links_resource_id: str
    registry_legacy_id_column: int
    links_legacy_id_column: int

    @classmethod
    def from_mapping(cls, data: Mapping[str, Any]) -> "RegulatoryJoin":
        expected = {
            "registry_resource_id",
            "links_resource_id",
            "registry_legacy_id_column",
            "links_legacy_id_column",
        }
        if set(data) != expected:
            raise VariantClassResourceError("regulatory join fields are invalid")
        value = cls(
            registry_resource_id=_required_string(
                data, "registry_resource_id", "regulatory join"
            ),
            links_resource_id=_required_string(
                data, "links_resource_id", "regulatory join"
            ),
            registry_legacy_id_column=_nonnegative_int(
                data, "registry_legacy_id_column", "regulatory join"
            ),
            links_legacy_id_column=_nonnegative_int(
                data, "links_legacy_id_column", "regulatory join"
            ),
        )
        return value


@dataclass(frozen=True, slots=True)
class VariantClassBundleContract:
    bundle_id: str
    resources: tuple[ClassResource, ...]
    regulatory_join: RegulatoryJoin
    schema: str = SCHEMA

    @classmethod
    def from_mapping(cls, data: Mapping[str, Any]) -> "VariantClassBundleContract":
        expected = {"schema", "bundle_id", "resources", "regulatory_join"}
        resources_data = data.get("resources")
        join_data = data.get("regulatory_join")
        if (
            set(data) != expected
            or data.get("schema") != SCHEMA
            or not isinstance(resources_data, list)
            or not resources_data
            or not all(isinstance(item, Mapping) for item in resources_data)
            or not isinstance(join_data, Mapping)
        ):
            raise VariantClassResourceError("class resource bundle fields are invalid")
        resources = tuple(ClassResource.from_mapping(item) for item in resources_data)
        contract = cls(
            schema=SCHEMA,
            bundle_id=_required_string(data, "bundle_id", "bundle"),
            resources=resources,
            regulatory_join=RegulatoryJoin.from_mapping(join_data),
        )
        resource_ids = [resource.resource_id for resource in resources]
        paths = [artifact.path for resource in resources for artifact in resource.artifacts]
        release_roots = [resource.release_root for resource in resources]
        lane_counts = Counter(resource.lane for resource in resources)
        by_id = {resource.resource_id: resource for resource in resources}
        join = contract.regulatory_join
        if (
            not IDENTIFIER.fullmatch(contract.bundle_id)
            or len(resource_ids) != len(set(resource_ids))
            or len(paths) != len(set(paths))
            or len(release_roots) != len(set(release_roots))
            or set(lane_counts) != EXPECTED_LANES
            or lane_counts != EXPECTED_LANE_COUNTS
            or join.registry_resource_id not in by_id
            or join.links_resource_id not in by_id
            or by_id[join.registry_resource_id].validation["kind"] != "bed"
            or by_id[join.links_resource_id].validation["kind"] != "tsv"
            or by_id[join.registry_resource_id].lane != "regulatory"
            or by_id[join.links_resource_id].lane != "regulatory"
        ):
            raise VariantClassResourceError("class resource bundle identities are inconsistent")
        registry_width = by_id[join.registry_resource_id].validation["column_count"]
        links_width = by_id[join.links_resource_id].validation["column_count"]
        if (
            join.registry_legacy_id_column >= registry_width
            or join.links_legacy_id_column >= links_width
        ):
            raise VariantClassResourceError("regulatory join column is out of range")
        registry_fields = by_id[join.registry_resource_id].validation["sentinel"][
            "fields"
        ]
        links_fields = by_id[join.links_resource_id].validation["sentinel"]["fields"]
        if (
            registry_fields[join.registry_legacy_id_column]
            != links_fields[join.links_legacy_id_column]
        ):
            raise VariantClassResourceError(
                "regulatory sentinel join differs in the bundle contract"
            )
        sv_profiles = Counter(
            "sv"
            if {"AC", "AN", "AF"}.issubset(
                resource.validation["required_info_fields"]
            )
            else "cnv"
            for resource in resources
            if resource.lane == "sv_cnv"
        )
        if sv_profiles != {"sv": 1, "cnv": 1}:
            raise VariantClassResourceError(
                "bundle requires separate SV and CNV population resources"
            )
        return contract

    def resource(self, resource_id: str) -> ClassResource:
        for resource in self.resources:
            if resource.resource_id == resource_id:
                return resource
        raise VariantClassResourceError(f"bundle lacks resource {resource_id}")


def load_variant_class_bundle_contract(
    path: str | Path,
) -> VariantClassBundleContract:
    return _read_variant_class_bundle_contract(Path(path))[0]


def _read_variant_class_bundle_contract(
    source: Path,
) -> tuple[VariantClassBundleContract, str]:
    try:
        content = source.read_bytes()
        data = json.loads(content)
    except (OSError, UnicodeDecodeError, json.JSONDecodeError) as error:
        raise VariantClassResourceError(
            f"cannot read class resource contract {source}: {error}"
        ) from error
    if not isinstance(data, Mapping):
        raise VariantClassResourceError("class resource contract root must be an object")
    return VariantClassBundleContract.from_mapping(data), hashlib.sha256(content).hexdigest()


@dataclass(frozen=True, slots=True)
class ResourceValidation:
    resource_id: str
    lane: str
    kind: str
    sentinel: str
    details: Mapping[str, Any]


@dataclass(frozen=True, slots=True)
class VariantClassBootstrapResult:
    bundle_id: str
    bundle_manifest_path: str
    bundle_manifest_sha256: str
    manifest_paths: tuple[str, ...]
    validation: tuple[ResourceValidation, ...]
    reused: bool


Downloader = Callable[[DownloadRequest, str | Path], DownloadResult]
CapacityChecker = Callable[[str | Path, CapacityPlan], Any]


def _file_identity(path: Path) -> tuple[str, str, int]:
    sha256 = hashlib.sha256()
    md5 = hashlib.md5(usedforsecurity=False)
    size = 0
    try:
        with path.open("rb") as source:
            while chunk := source.read(1024 * 1024):
                sha256.update(chunk)
                md5.update(chunk)
                size += len(chunk)
    except OSError as error:
        raise VariantClassResourceError(f"cannot read class resource {path}: {error}") from error
    return sha256.hexdigest(), md5.hexdigest(), size


def _verify_artifact(path: Path, artifact: ResourceArtifact) -> tuple[str, int]:
    sha256, md5, size = _file_identity(path)
    if (
        size != artifact.size_bytes
        or artifact.sha256 is not None
        and sha256 != artifact.sha256
        or artifact.expected_md5 is not None
        and md5 != artifact.expected_md5
    ):
        raise VariantClassResourceError(
            f"class resource checksum mismatch: {artifact.role}"
        )
    return sha256, size


def _text_source(path: Path) -> TextIO:
    try:
        if path.suffix in {".gz", ".bgz"}:
            return gzip.open(path, "rt", encoding="utf-8")
        return path.open("rt", encoding="utf-8")
    except OSError as error:
        raise VariantClassResourceError(f"cannot open class resource {path}: {error}") from error


def _parse_info(value: str) -> Mapping[str, str | None]:
    parsed: dict[str, str | None] = {}
    for item in value.split(";"):
        key, separator, field_value = item.partition("=")
        if not key or key in parsed:
            raise VariantClassResourceError("VCF record has malformed or duplicate INFO")
        parsed[key] = field_value if separator else None
    return parsed


def _validated_count_frequency(
    info: Mapping[str, str | None],
    count_field: str,
    denominator_field: str,
    frequency_field: str,
    context: str,
) -> tuple[int, int, float]:
    count = int(str(info[count_field]).split(",", 1)[0])
    denominator = int(str(info[denominator_field]).split(",", 1)[0])
    frequency = float(str(info[frequency_field]).split(",", 1)[0])
    if (
        denominator <= 0
        or count < 0
        or count > denominator
        or not math.isfinite(frequency)
        or not 0 <= frequency <= 1
    ):
        raise VariantClassResourceError(f"{context} count or frequency is invalid")
    derived = count / denominator
    if abs(frequency - derived) > max(1e-8, derived * 0.02):
        raise VariantClassResourceError(
            f"{context} frequency disagrees with count/denominator"
        )
    return count, denominator, frequency


def _read_vcf_header(path: Path, required_fields: tuple[str, ...]) -> tuple[str, ...]:
    info_fields: list[str] = []
    try:
        with _text_source(path) as source:
            for line in source:
                value = line.rstrip("\r\n")
                match = INFO_HEADER.match(value)
                if match:
                    info_fields.append(match.group(1))
                if value.startswith("#CHROM\t"):
                    columns = tuple(value.split("\t"))
                    if columns != (
                        "#CHROM",
                        "POS",
                        "ID",
                        "REF",
                        "ALT",
                        "QUAL",
                        "FILTER",
                        "INFO",
                    ):
                        raise VariantClassResourceError(
                            "class resource VCF must contain eight site-only columns"
                        )
                    break
                if not value.startswith("#"):
                    raise VariantClassResourceError("class resource VCF header is incomplete")
            else:
                raise VariantClassResourceError("class resource VCF lacks column declarations")
    except (OSError, EOFError, gzip.BadGzipFile, UnicodeDecodeError) as error:
        raise VariantClassResourceError(f"cannot read class resource VCF: {error}") from error
    if len(info_fields) != len(set(info_fields)):
        raise VariantClassResourceError("class resource VCF has duplicate INFO declarations")
    missing = sorted(set(required_fields) - set(info_fields))
    if missing:
        raise VariantClassResourceError(
            f"class resource VCF lacks INFO fields: {', '.join(missing)}"
        )
    return tuple(info_fields)


def _scan_vcf(path: Path, sentinel: Mapping[str, Any]) -> list[str]:
    try:
        with _text_source(path) as source:
            for line in source:
                if line.startswith("#"):
                    continue
                fields = line.rstrip("\r\n").split("\t")
                if len(fields) != 8:
                    raise VariantClassResourceError("class resource VCF row width is invalid")
                if fields[:5] == [
                    sentinel["chromosome"],
                    str(sentinel["position"]),
                    sentinel["id"],
                    sentinel["reference"],
                    sentinel["alternate"],
                ]:
                    return fields
    except (OSError, EOFError, gzip.BadGzipFile, UnicodeDecodeError) as error:
        raise VariantClassResourceError(f"cannot scan class resource VCF: {error}") from error
    raise VariantClassResourceError("class resource VCF sentinel was not found")


def _validate_vcf(
    resource: ClassResource, paths: Mapping[str, Path]
) -> ResourceValidation:
    validation = resource.validation
    required = tuple(validation["required_info_fields"])
    primary = paths[validation["artifact_role"]]
    _read_vcf_header(primary, required)
    sentinel = validation["sentinel"]
    index_role = validation["index_artifact_role"]
    if index_role is None:
        fields = _scan_vcf(primary, sentinel)
    else:
        try:
            index = TabixIndex.load(paths[index_role])
            index.require_vcf_layout()
            records = index.query(
                primary,
                sentinel["chromosome"],
                sentinel["position"],
                sentinel["position"],
            )
        except (KeyError, OSError, TabixError) as error:
            raise VariantClassResourceError(
                f"class resource tabix query failed: {error}"
            ) from error
        matches = [
            record.rstrip("\r\n").split("\t")
            for record in records
            if record.rstrip("\r\n").split("\t")[:5]
            == [
                sentinel["chromosome"],
                str(sentinel["position"]),
                sentinel["id"],
                sentinel["reference"],
                sentinel["alternate"],
            ]
        ]
        if len(matches) != 1:
            raise VariantClassResourceError(
                "class resource index did not recover exactly one sentinel"
            )
        fields = matches[0]
    if len(fields) != 8 or fields[6] != sentinel["filter"]:
        raise VariantClassResourceError("class resource VCF sentinel filter is invalid")
    info = _parse_info(fields[7])
    expected_info = sentinel["info"]
    if any(info.get(key) != value for key, value in expected_info.items()):
        raise VariantClassResourceError("class resource VCF sentinel INFO differs")
    if any(info.get(key) is None for key in required):
        raise VariantClassResourceError("class resource VCF sentinel lacks INFO values")
    try:
        if {"AC", "AN", "AF"}.issubset(info):
            _validated_count_frequency(info, "AC", "AN", "AF", "SV allele")
        if {"SC", "SN", "SF"}.issubset(info):
            _validated_count_frequency(info, "SC", "SN", "SF", "CNV site")
        if resource.lane == "mitochondrial":
            hom_count, allele_number, _ = _validated_count_frequency(
                info, "AC_hom", "AN", "AF_hom", "mitochondrial homoplasmic"
            )
            het_count, het_allele_number, _ = _validated_count_frequency(
                info, "AC_het", "AN", "AF_het", "mitochondrial heteroplasmic"
            )
            maximum = float(str(info["max_hl"]))
            if (
                het_allele_number != allele_number
                or hom_count + het_count > allele_number
                or not math.isfinite(maximum)
                or not 0 <= maximum <= 1
            ):
                raise VariantClassResourceError(
                    "mitochondrial sentinel counts are invalid"
                )
    except (TypeError, ValueError, ZeroDivisionError) as error:
        raise VariantClassResourceError(
            "class resource VCF sentinel frequency values are invalid"
        ) from error
    return ResourceValidation(
        resource_id=resource.resource_id,
        lane=resource.lane,
        kind="vcf",
        sentinel=f"{fields[0]}:{fields[1]}:{fields[3]}:{fields[4]}",
        details={key: info[key] for key in required},
    )


def _validate_delimited(
    resource: ClassResource, paths: Mapping[str, Path]
) -> ResourceValidation:
    validation = resource.validation
    primary = paths[validation["artifact_role"]]
    expected = list(validation["sentinel"]["fields"])
    count = validation["column_count"]
    found: list[str] | None = None
    try:
        with _text_source(primary) as source:
            for line in source:
                if not line.strip() or line.startswith("#"):
                    continue
                fields = line.rstrip("\r\n").split("\t")
                if len(fields) != count:
                    raise VariantClassResourceError(
                        f"{validation['kind']} resource row width is invalid"
                    )
                if fields == expected:
                    found = fields
                    break
    except (OSError, EOFError, gzip.BadGzipFile, UnicodeDecodeError) as error:
        raise VariantClassResourceError(
            f"cannot scan {validation['kind']} class resource: {error}"
        ) from error
    if found is None:
        raise VariantClassResourceError(
            f"{validation['kind']} class resource sentinel was not found"
        )
    details: dict[str, Any] = {"fields": tuple(found)}
    if validation["kind"] == "bed":
        try:
            start, end = int(found[1]), int(found[2])
        except (IndexError, ValueError) as error:
            raise VariantClassResourceError("BED sentinel coordinates are invalid") from error
        if start < 0 or end <= start:
            raise VariantClassResourceError("BED sentinel interval is invalid")
        if resource.lane == "mosaic":
            details.update(
                {
                    "interval_convention": validation["interval_convention"],
                    "interpretation": validation["interpretation"],
                }
            )
    return ResourceValidation(
        resource_id=resource.resource_id,
        lane=resource.lane,
        kind=validation["kind"],
        sentinel="\t".join(found),
        details=details,
    )


def _validate_repeat_catalog(
    resource: ClassResource, paths: Mapping[str, Path]
) -> ResourceValidation:
    validation = resource.validation
    try:
        records = json.loads(paths[validation["artifact_role"]].read_text(encoding="utf-8"))
    except (OSError, UnicodeDecodeError, json.JSONDecodeError) as error:
        raise VariantClassResourceError(f"cannot read repeat catalog: {error}") from error
    if not isinstance(records, list) or len(records) != validation["expected_record_count"]:
        raise VariantClassResourceError("repeat catalog record count differs")
    if any(
        not isinstance(record, Mapping)
        or not all(
            field in record
            for field in ("LocusId", "ReferenceRegion", "LocusStructure", "VariantType")
        )
        for record in records
    ):
        raise VariantClassResourceError("repeat catalog contains an invalid record")
    sentinel = validation["sentinel"]
    matches = [record for record in records if record.get("LocusId") == sentinel["locus_id"]]
    if len(matches) != 1:
        raise VariantClassResourceError("repeat catalog sentinel is not unique")
    record = matches[0]
    diseases = record.get("Diseases")
    if not isinstance(diseases, list):
        raise VariantClassResourceError("repeat catalog sentinel diseases are invalid")
    disease_matches = [
        disease
        for disease in diseases
        if isinstance(disease, Mapping)
        and disease.get("Symbol") == sentinel["disease_symbol"]
    ]
    expected_record = {
        "MainReferenceRegion": sentinel["main_reference_region"],
        "ReferenceRegion": sentinel["reference_regions"],
        "RepeatUnit": sentinel["repeat_unit"],
        "LocusStructure": sentinel["locus_structure"],
        "VariantType": sentinel["variant_types"],
        "Inheritance": sentinel["inheritance"],
    }
    expected_disease = {
        "Inheritance": sentinel["inheritance"],
        "NormalMax": sentinel["normal_max"],
        "PathogenicMin": sentinel["pathogenic_min"],
        "IntermediateRange": sentinel["intermediate_range"],
    }
    if (
        any(record.get(key) != value for key, value in expected_record.items())
        or len(disease_matches) != 1
        or any(
            disease_matches[0].get(key) != value
            for key, value in expected_disease.items()
        )
    ):
        raise VariantClassResourceError("repeat catalog sentinel fields differ")
    return ResourceValidation(
        resource_id=resource.resource_id,
        lane=resource.lane,
        kind="repeat_catalog",
        sentinel=sentinel["locus_id"],
        details={**expected_record, **expected_disease},
    )


def validate_class_resource(
    resource: ClassResource, paths: Mapping[str, Path]
) -> ResourceValidation:
    if set(paths) != {artifact.role for artifact in resource.artifacts}:
        raise VariantClassResourceError(
            f"resource paths are incomplete: {resource.resource_id}"
        )
    for artifact in resource.artifacts:
        _verify_artifact(paths[artifact.role], artifact)
    kind = resource.validation["kind"]
    if kind == "vcf":
        return _validate_vcf(resource, paths)
    if kind in {"bed", "tsv"}:
        return _validate_delimited(resource, paths)
    return _validate_repeat_catalog(resource, paths)


def _validate_regulatory_join(
    contract: VariantClassBundleContract,
    validation: Mapping[str, ResourceValidation],
) -> str:
    join = contract.regulatory_join
    registry = validation[join.registry_resource_id].details["fields"]
    links = validation[join.links_resource_id].details["fields"]
    registry_id = registry[join.registry_legacy_id_column]
    links_id = links[join.links_legacy_id_column]
    if registry_id != links_id:
        raise VariantClassResourceError("regulatory cCRE-to-gene sentinel join failed")
    return str(registry_id)


def _require_within_root(root: Path, path: Path, context: str) -> None:
    if not path.is_absolute() or not path.is_relative_to(root):
        raise VariantClassResourceError(f"{context} escapes reference_root: {path}")
    relative = path.relative_to(root)
    if ".." in relative.parts:
        raise VariantClassResourceError(f"{context} escapes reference_root: {path}")
    current = root
    for part in relative.parts:
        current /= part
        if current.is_symlink():
            raise VariantClassResourceError(f"{context} contains a symlink: {current}")


def _verify_release_tree(root: Path, resource: ClassResource) -> None:
    release = root / str(resource.release_root)
    _require_within_root(root, release, "class resource release")
    if release.is_symlink() or not release.is_dir():
        raise VariantClassResourceError(
            f"class resource release is not a directory: {resource.resource_id}"
        )
    expected = {PurePosixPath(artifact.path).name for artifact in resource.artifacts}
    try:
        entries = tuple(release.iterdir())
    except OSError as error:
        raise VariantClassResourceError(
            f"cannot inspect class resource release: {error}"
        ) from error
    if {path.name for path in entries} != expected or any(
        path.is_symlink() or not path.is_file() for path in entries
    ):
        raise VariantClassResourceError(
            f"class resource release is incomplete or contaminated: {resource.resource_id}"
        )
    for artifact in resource.artifacts:
        _verify_artifact(root / artifact.path, artifact)


def _copy_file(source: Path, destination: Path) -> None:
    try:
        with source.open("rb") as input_file, destination.open("xb") as output_file:
            shutil.copyfileobj(input_file, output_file, length=1024 * 1024)
            output_file.flush()
            os.fsync(output_file.fileno())
        os.chmod(destination, 0o444)
    except OSError as error:
        raise VariantClassResourceError(f"cannot stage class resource: {error}") from error


def _publish_release(
    root: Path, resource: ClassResource, sources: Mapping[str, Path]
) -> None:
    release = root / str(resource.release_root)
    _require_within_root(root, release, "class resource release")
    if release.exists() or release.is_symlink():
        _verify_release_tree(root, resource)
        return
    release.parent.mkdir(parents=True, exist_ok=True)
    partial_root = release.parent / ".partial"
    _require_within_root(root, partial_root, "class resource release partial path")
    partial_root.mkdir(mode=0o700, exist_ok=True)
    staging = partial_root / f"{release.name}.{uuid.uuid4().hex}.partial"
    staging.mkdir(mode=0o700)
    try:
        for artifact in resource.artifacts:
            _copy_file(
                sources[artifact.role], staging / PurePosixPath(artifact.path).name
            )
        os.chmod(staging, 0o755)
        staging_fd = os.open(staging, os.O_RDONLY)
        try:
            os.fsync(staging_fd)
        finally:
            os.close(staging_fd)
        try:
            os.replace(staging, release)
        except OSError:
            if not release.exists():
                raise
            shutil.rmtree(staging, ignore_errors=True)
            _verify_release_tree(root, resource)
        directory_fd = os.open(release.parent, os.O_RDONLY)
        try:
            os.fsync(directory_fd)
        finally:
            os.close(directory_fd)
    except Exception:
        shutil.rmtree(staging, ignore_errors=True)
        raise
    _verify_release_tree(root, resource)


def _manifest_path(root: Path, resource: ClassResource) -> Path:
    return root / f"manifests/variant-class-resources/{resource.resource_id}.json"


def _load_manifest(path: Path) -> PublicResourceManifest:
    try:
        data = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as error:
        raise VariantClassResourceError(
            f"cannot read class resource manifest {path}: {error}"
        ) from error
    if not isinstance(data, Mapping):
        raise VariantClassResourceError("class resource manifest root must be an object")
    try:
        return PublicResourceManifest.from_mapping(data)
    except ManifestError as error:
        raise VariantClassResourceError(f"invalid class resource manifest: {error}") from error


def _resource_manifest(
    resource: ClassResource,
    contract_sha256: str,
    image_digest: str,
    retrieved_on: str,
    root: Path,
) -> PublicResourceManifest:
    outputs = []
    for artifact in resource.artifacts:
        sha256, _, size = _file_identity(root / artifact.path)
        outputs.append(
            ResourceOutput(
                path=artifact.path,
                size_bytes=size,
                sha256=sha256,
                retention=RetentionPolicy("durable", None),
            )
        )
    return PublicResourceManifest(
        resource_id=resource.resource_id,
        release=resource.release,
        source=resource.artifact(resource.primary_artifact_role).source_identity(),
        license=resource.license,
        assembly=resource.assembly,
        retrieved_on=retrieved_on,
        transformation=Transformation(TRANSFORMATION_SCHEMA, contract_sha256),
        image_digest=image_digest,
        outputs=tuple(outputs),
    )


def _verify_manifest_contract(
    manifest: PublicResourceManifest,
    resource: ClassResource,
    contract_sha256: str,
) -> None:
    outputs = {output.path: output for output in manifest.outputs}
    if (
        manifest.resource_id != resource.resource_id
        or manifest.release != resource.release
        or manifest.source
        != resource.artifact(resource.primary_artifact_role).source_identity()
        or manifest.license != resource.license
        or manifest.assembly != resource.assembly
        or manifest.transformation
        != Transformation(TRANSFORMATION_SCHEMA, contract_sha256)
        or set(outputs) != {artifact.path for artifact in resource.artifacts}
        or any(
            outputs[artifact.path].size_bytes != artifact.size_bytes
            or artifact.sha256 is not None
            and outputs[artifact.path].sha256 != artifact.sha256
            for artifact in resource.artifacts
        )
    ):
        raise VariantClassResourceError(
            f"class resource manifest differs from contract: {resource.resource_id}"
        )


def _bundle_manifest_path(root: Path, contract: VariantClassBundleContract) -> Path:
    return root / f"manifests/variant-class-resources/{contract.bundle_id}.bundle.json"


def _bundle_manifest_data(
    contract: VariantClassBundleContract,
    contract_sha256: str,
    manifests: Mapping[str, tuple[Path, str]],
    validations: tuple[ResourceValidation, ...],
    regulatory_join_id: str,
    root: Path,
) -> Mapping[str, Any]:
    resources = []
    for resource in contract.resources:
        manifest_path, manifest_sha256 = manifests[resource.resource_id]
        resources.append(
            {
                "lane": resource.lane,
                "resource_id": resource.resource_id,
                "release": resource.release,
                "assembly": resource.assembly,
                "sequence_identity": resource.sequence_identity,
                "license": resource.license,
                "license_url": resource.license_url,
                "sources": [
                    {
                        "role": artifact.role,
                        "url": artifact.url,
                        "sha256": artifact.sha256,
                        "expected_md5": artifact.expected_md5,
                        "object_generation": artifact.object_generation,
                    }
                    for artifact in resource.artifacts
                ],
                "manifest_path": str(manifest_path.relative_to(root)),
                "manifest_sha256": manifest_sha256,
            }
        )
    return {
        "schema": BUNDLE_MANIFEST_SCHEMA,
        "classification": "public",
        "bundle_id": contract.bundle_id,
        "contract_sha256": contract_sha256,
        "resources": resources,
        "validation": [asdict(item) for item in validations],
        "regulatory_join_sentinel": regulatory_join_id,
    }


def _publish_json_immutable(data: Mapping[str, Any], path: Path, root: Path) -> str:
    _require_within_root(root, path, "class resource bundle manifest")
    content = (json.dumps(data, sort_keys=True, separators=(",", ":")) + "\n").encode()
    digest = hashlib.sha256(content).hexdigest()
    if path.exists() or path.is_symlink():
        if path.is_symlink() or not path.is_file():
            raise VariantClassResourceError(
                "class resource bundle manifest is not a regular file"
            )
        if path.read_bytes() != content:
            raise VariantClassResourceError(
                "completed class resource bundle manifest is immutable"
            )
        return digest
    path.parent.mkdir(parents=True, exist_ok=True)
    partial_root = path.parent / ".partial"
    _require_within_root(root, partial_root, "class resource manifest partial path")
    partial_root.mkdir(mode=0o700, exist_ok=True)
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
            if path.is_symlink() or not path.is_file() or path.read_bytes() != content:
                raise VariantClassResourceError(
                    "concurrent class resource bundle manifest differs"
                )
        directory_fd = os.open(path.parent, os.O_RDONLY)
        try:
            os.fsync(directory_fd)
        finally:
            os.close(directory_fd)
    finally:
        partial.unlink(missing_ok=True)
    return digest


def bootstrap_variant_class_resources(
    *,
    contract_path: str | Path,
    reference_root: str | Path,
    image_digest: str,
    retrieved_on: str | None = None,
    scratch_root: str | Path | None = None,
    reserve_bytes: int = DEFAULT_RESERVE_BYTES,
    downloader: Downloader = download_atomic,
    capacity_checker: CapacityChecker = require_capacity,
) -> VariantClassBootstrapResult:
    contract_file = Path(contract_path)
    try:
        root = Path(reference_root).resolve(strict=True)
    except OSError as error:
        raise VariantClassResourceError(
            f"cannot resolve class resource reference_root: {error}"
        ) from error
    contract, contract_sha256 = _read_variant_class_bundle_contract(contract_file)
    if not DIGEST.fullmatch(image_digest):
        raise VariantClassResourceError(
            "class resource bootstrap image_digest must be a SHA-256 digest"
        )
    bundle_manifest_path = _bundle_manifest_path(root, contract)
    _require_within_root(root, bundle_manifest_path, "class resource bundle manifest")
    bundle_existed = bundle_manifest_path.exists() or bundle_manifest_path.is_symlink()
    if bundle_existed and (
        bundle_manifest_path.is_symlink() or not bundle_manifest_path.is_file()
    ):
        raise VariantClassResourceError(
            "class resource bundle manifest is not a regular file"
        )
    scratch_base = (
        Path(scratch_root) if scratch_root is not None else root / "scratch/public-build"
    )
    if not scratch_base.is_absolute():
        raise VariantClassResourceError("class resource scratch_root must be absolute")
    _require_within_root(root, scratch_base, "class resource scratch_root")
    retrieval_date = retrieved_on or date.today().isoformat()
    publications: dict[str, tuple[Path, str]] = {}
    validations: list[ResourceValidation] = []

    for resource in contract.resources:
        manifest_path = _manifest_path(root, resource)
        _require_within_root(root, manifest_path, "class resource manifest")
        release = root / str(resource.release_root)
        _require_within_root(root, release, "class resource release")
        if manifest_path.exists() or manifest_path.is_symlink():
            manifest = _load_manifest(manifest_path)
            _verify_manifest_contract(manifest, resource, contract_sha256)
            try:
                verify_outputs(manifest, root)
            except ManifestError as error:
                raise VariantClassResourceError(
                    f"completed class resource is invalid: {resource.resource_id}: {error}"
                ) from error
            _verify_release_tree(root, resource)
        else:
            if bundle_existed:
                raise VariantClassResourceError(
                    "completed class resource bundle is missing a resource manifest: "
                    f"{resource.resource_id}"
                )
            if not release.exists():
                scratch = scratch_base / "variant-class-resources" / resource.resource_id
                _require_within_root(root, scratch, "class resource scratch")
                scratch.mkdir(parents=True, exist_ok=True, mode=0o700)
                _require_within_root(root, scratch, "class resource scratch")
                source_paths = {
                    artifact.role: scratch / PurePosixPath(artifact.path).name
                    for artifact in resource.artifacts
                }
                for source_path in source_paths.values():
                    _require_within_root(
                        root, source_path, "class resource scratch artifact"
                    )
                    _require_within_root(
                        root,
                        source_path.parent / ".partial",
                        "class resource download partial path",
                    )
                missing = [
                    artifact
                    for artifact in resource.artifacts
                    if not source_paths[artifact.role].exists()
                ]
                capacity_checker(
                    root,
                    CapacityPlan(
                        operation=f"bootstrap {resource.resource_id}",
                        source_bytes=sum(item.size_bytes for item in missing),
                        output_bytes=sum(item.size_bytes for item in resource.artifacts),
                        reserve_bytes=reserve_bytes,
                    ),
                )
                for artifact in resource.artifacts:
                    downloader(
                        DownloadRequest(
                            url=artifact.url,
                            destination=source_paths[artifact.role],
                            expected_size_bytes=artifact.size_bytes,
                            expected_sha256=artifact.sha256,
                            expected_md5=artifact.expected_md5,
                            expected_generation=artifact.object_generation,
                            attempts=3,
                            timeout_seconds=600,
                        ),
                        root,
                    )
                validate_class_resource(resource, source_paths)
                _publish_release(root, resource, source_paths)
                for path in source_paths.values():
                    path.unlink(missing_ok=True)
            else:
                if any(
                    artifact.sha256 is None and artifact.expected_md5 is None
                    for artifact in resource.artifacts
                ):
                    raise VariantClassResourceError(
                        "cannot adopt a generation-only class resource release "
                        f"without its manifest: {resource.resource_id}"
                    )
                _verify_release_tree(root, resource)
            manifest = _resource_manifest(
                resource, contract_sha256, image_digest, retrieval_date, root
            )
            publication = publish_manifest(manifest, manifest_path, root)
            manifest_path = publication.path
        final_paths = {
            artifact.role: root / artifact.path for artifact in resource.artifacts
        }
        validations.append(validate_class_resource(resource, final_paths))
        publications[resource.resource_id] = (
            manifest_path,
            _file_identity(manifest_path)[0],
        )

    validation_tuple = tuple(validations)
    by_resource = {item.resource_id: item for item in validation_tuple}
    join_id = _validate_regulatory_join(contract, by_resource)
    bundle_data = _bundle_manifest_data(
        contract,
        contract_sha256,
        publications,
        validation_tuple,
        join_id,
        root,
    )
    bundle_sha256 = _publish_json_immutable(
        bundle_data, bundle_manifest_path, root
    )
    return VariantClassBootstrapResult(
        bundle_id=contract.bundle_id,
        bundle_manifest_path=str(bundle_manifest_path),
        bundle_manifest_sha256=bundle_sha256,
        manifest_paths=tuple(
            str(publications[resource.resource_id][0])
            for resource in contract.resources
        ),
        validation=validation_tuple,
        reused=bundle_existed,
    )


def variant_class_result_as_json(result: VariantClassBootstrapResult) -> str:
    return json.dumps(asdict(result), indent=2, sort_keys=True) + "\n"
