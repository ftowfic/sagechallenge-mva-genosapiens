"""Bootstrap and validate a pinned Open Targets Platform Parquet subset."""

from __future__ import annotations

from collections.abc import Callable, Mapping, Sequence
from dataclasses import asdict, dataclass
from datetime import date
from fnmatch import fnmatchcase
import hashlib
import json
import os
from pathlib import Path, PurePosixPath
import re
import shutil
from typing import Any

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


SCHEMA = "genosapiens.open-targets-contract"
TRANSFORMATION_SCHEMA = "genosapiens.open-targets-parquet-bootstrap"
DIGEST = re.compile(r"^sha256:[0-9a-f]{64}$")
SHA256 = re.compile(r"^[0-9a-f]{64}$")
SHA1 = re.compile(r"^[0-9a-f]{40}$")
IDENTIFIER = re.compile(r"^[A-Za-z0-9][A-Za-z0-9._-]*$")
INTEGRITY_LINE = re.compile(r"^([0-9a-f]{40})  (\./[^\r\n]+)$")
METADATA_NAMES = {
    "croissant.json",
    "release_data_integrity",
    "release_data_integrity.sha1",
}
TYPE_NAMES = {
    "cr:Float32": "FLOAT",
    "cr:Float64": "DOUBLE",
    "cr:Int32": "INTEGER",
    "cr:Int64": "BIGINT",
    "sc:Boolean": "BOOLEAN",
    "sc:Text": "VARCHAR",
}


class OpenTargetsError(RuntimeError):
    """Raised when Open Targets acquisition or validation fails."""


def _required_string(data: Mapping[str, Any], name: str, context: str) -> str:
    value = data.get(name)
    if not isinstance(value, str) or not value:
        raise OpenTargetsError(f"{context}.{name} must be a non-empty string")
    return value


def _positive_int(data: Mapping[str, Any], name: str, context: str) -> int:
    value = data.get(name)
    if not isinstance(value, int) or isinstance(value, bool) or value <= 0:
        raise OpenTargetsError(f"{context}.{name} must be a positive integer")
    return value


def _safe_path(value: str, context: str) -> PurePosixPath:
    path = PurePosixPath(value)
    if path.is_absolute() or path == PurePosixPath(".") or ".." in path.parts:
        raise OpenTargetsError(f"{context} must be a safe relative path")
    return path


@dataclass(frozen=True, slots=True)
class MetadataDocument:
    name: str
    size_bytes: int
    sha1: str
    sha256: str

    @classmethod
    def from_mapping(cls, data: Mapping[str, Any]) -> "MetadataDocument":
        if set(data) != {"name", "size_bytes", "sha1", "sha256"}:
            raise OpenTargetsError("Open Targets metadata fields are invalid")
        document = cls(
            name=_required_string(data, "name", "metadata"),
            size_bytes=_positive_int(data, "size_bytes", "metadata"),
            sha1=_required_string(data, "sha1", "metadata"),
            sha256=_required_string(data, "sha256", "metadata"),
        )
        if (
            document.name not in METADATA_NAMES
            or not SHA1.fullmatch(document.sha1)
            or not SHA256.fullmatch(document.sha256)
        ):
            raise OpenTargetsError(f"invalid metadata identity: {document.name}")
        return document


@dataclass(frozen=True, slots=True)
class DatasetContract:
    name: str
    partition_count: int
    size_bytes: int
    row_count: int
    partition_size_bytes: tuple[int, ...]

    @classmethod
    def from_mapping(cls, data: Mapping[str, Any]) -> "DatasetContract":
        expected = {
            "name",
            "partition_count",
            "size_bytes",
            "row_count",
            "partition_size_bytes",
        }
        sizes = data.get("partition_size_bytes")
        if (
            set(data) != expected
            or not isinstance(sizes, list)
            or any(not isinstance(value, int) or isinstance(value, bool) or value <= 0 for value in sizes)
        ):
            raise OpenTargetsError("Open Targets dataset fields are invalid")
        dataset = cls(
            name=_required_string(data, "name", "dataset"),
            partition_count=_positive_int(data, "partition_count", "dataset"),
            size_bytes=_positive_int(data, "size_bytes", "dataset"),
            row_count=_positive_int(data, "row_count", "dataset"),
            partition_size_bytes=tuple(sizes),
        )
        if (
            not IDENTIFIER.fullmatch(dataset.name)
            or len(dataset.partition_size_bytes) != dataset.partition_count
            or sum(dataset.partition_size_bytes) != dataset.size_bytes
        ):
            raise OpenTargetsError(f"dataset partition identity is invalid: {dataset.name}")
        return dataset


@dataclass(frozen=True, slots=True)
class SemanticJoin:
    source: str
    target: str
    reason: str
    source_prefix: str | None = None

    @classmethod
    def from_mapping(cls, data: Mapping[str, Any]) -> "SemanticJoin":
        if set(data) not in (
            {"source", "target", "reason"},
            {"source", "target", "reason", "source_prefix"},
        ):
            raise OpenTargetsError("Open Targets semantic join fields are invalid")
        source_prefix = data.get("source_prefix")
        if source_prefix is not None and (
            not isinstance(source_prefix, str) or not source_prefix
        ):
            raise OpenTargetsError("semantic_join.source_prefix must be a non-empty string")
        join = cls(
            source=_required_string(data, "source", "semantic_join"),
            target=_required_string(data, "target", "semantic_join"),
            reason=_required_string(data, "reason", "semantic_join"),
            source_prefix=source_prefix,
        )
        if "/" not in join.source or "/" not in join.target:
            raise OpenTargetsError("semantic join fields must be Croissant field IDs")
        return join


@dataclass(frozen=True, slots=True)
class OpenTargetsContract:
    resource: str
    release: str
    published_date: str
    croissant_date_published: str
    base_url: str
    output_prefix: str
    assembly: str
    license: str
    croissant_license: str
    metadata: tuple[MetadataDocument, ...]
    datasets: tuple[DatasetContract, ...]
    semantic_joins: tuple[SemanticJoin, ...]
    schema: str = SCHEMA

    @classmethod
    def from_mapping(cls, data: Mapping[str, Any]) -> "OpenTargetsContract":
        expected = {
            "schema",
            "resource",
            "release",
            "published_date",
            "croissant_date_published",
            "base_url",
            "output_prefix",
            "assembly",
            "license",
            "croissant_license",
            "metadata",
            "datasets",
            "semantic_joins",
        }
        metadata = data.get("metadata")
        datasets = data.get("datasets")
        joins = data.get("semantic_joins")
        if (
            set(data) != expected
            or data.get("schema") != SCHEMA
            or not isinstance(metadata, list)
            or not all(isinstance(item, Mapping) for item in metadata)
            or not isinstance(datasets, list)
            or not datasets
            or not all(isinstance(item, Mapping) for item in datasets)
            or not isinstance(joins, list)
            or not all(isinstance(item, Mapping) for item in joins)
        ):
            raise OpenTargetsError("Open Targets contract fields or schema are invalid")
        contract = cls(
            resource=_required_string(data, "resource", "contract"),
            release=_required_string(data, "release", "contract"),
            published_date=_required_string(data, "published_date", "contract"),
            croissant_date_published=_required_string(
                data, "croissant_date_published", "contract"
            ),
            base_url=_required_string(data, "base_url", "contract"),
            output_prefix=_required_string(data, "output_prefix", "contract"),
            assembly=_required_string(data, "assembly", "contract"),
            license=_required_string(data, "license", "contract"),
            croissant_license=_required_string(data, "croissant_license", "contract"),
            metadata=tuple(MetadataDocument.from_mapping(item) for item in metadata),
            datasets=tuple(DatasetContract.from_mapping(item) for item in datasets),
            semantic_joins=tuple(SemanticJoin.from_mapping(item) for item in joins),
        )
        try:
            date.fromisoformat(contract.published_date)
            date.fromisoformat(contract.croissant_date_published)
        except ValueError as error:
            raise OpenTargetsError("Open Targets publication dates must be ISO dates") from error
        metadata_names = [item.name for item in contract.metadata]
        dataset_names = [item.name for item in contract.datasets]
        if (
            not IDENTIFIER.fullmatch(contract.resource)
            or not re.fullmatch(r"[0-9]{2}\.[0-9]{2}", contract.release)
            or contract.base_url.endswith("/")
            or not contract.base_url.startswith("https://")
            or not contract.croissant_license.startswith("https://")
            or str(_safe_path(contract.output_prefix, "output_prefix"))
            != contract.output_prefix
            or set(metadata_names) != METADATA_NAMES
            or len(metadata_names) != len(set(metadata_names))
            or dataset_names != sorted(dataset_names)
            or len(dataset_names) != len(set(dataset_names))
            or not contract.output_prefix.endswith(f"/{contract.release}")
        ):
            raise OpenTargetsError("Open Targets contract identities are inconsistent")
        return contract

    def metadata_document(self, name: str) -> MetadataDocument:
        for document in self.metadata:
            if document.name == name:
                return document
        raise OpenTargetsError(f"Open Targets metadata document is missing: {name}")


def load_open_targets_contract(path: str | Path) -> OpenTargetsContract:
    source = Path(path)
    try:
        value = json.loads(source.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as error:
        raise OpenTargetsError(f"cannot read Open Targets contract {source}: {error}") from error
    if not isinstance(value, Mapping):
        raise OpenTargetsError("Open Targets contract root must be an object")
    return OpenTargetsContract.from_mapping(value)


@dataclass(frozen=True, slots=True)
class ReleaseArtifact:
    dataset: str | None
    relative_path: str
    size_bytes: int
    sha1: str
    sha256: str | None = None


@dataclass(frozen=True, slots=True)
class FieldModel:
    field_id: str
    expression: str
    data_type: str | None
    repeated_depth: int
    nested: bool


@dataclass(frozen=True, slots=True)
class DatasetModel:
    name: str
    fields: tuple[FieldModel, ...]
    top_level_fields: tuple[str, ...]
    key_fields: tuple[str, ...]


@dataclass(frozen=True, slots=True)
class DatasetValidation:
    name: str
    partitions: int
    size_bytes: int
    row_count: int
    key_fields: tuple[str, ...]
    nested_fields: int


@dataclass(frozen=True, slots=True)
class OpenTargetsValidation:
    release: str
    datasets: tuple[DatasetValidation, ...]
    partition_count: int
    size_bytes: int
    row_count: int
    declared_joins: int
    semantic_joins: int


@dataclass(frozen=True, slots=True)
class OpenTargetsBootstrapResult:
    resource: str
    release: str
    manifest_path: str
    manifest_sha256: str
    validation: OpenTargetsValidation
    reused: bool


Downloader = Callable[[DownloadRequest, str | Path], DownloadResult]
CapacityChecker = Callable[[str | Path, CapacityPlan], Any]


def _identity(path: Path) -> tuple[str, str, int]:
    sha256 = hashlib.sha256()
    sha1 = hashlib.sha1(usedforsecurity=False)
    size = 0
    try:
        with path.open("rb") as source:
            while chunk := source.read(1024 * 1024):
                sha256.update(chunk)
                sha1.update(chunk)
                size += len(chunk)
    except OSError as error:
        raise OpenTargetsError(f"cannot read Open Targets artifact {path}: {error}") from error
    return sha256.hexdigest(), sha1.hexdigest(), size


def _json_object(path: Path, context: str) -> Mapping[str, Any]:
    try:
        value = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as error:
        raise OpenTargetsError(f"cannot read {context} {path}: {error}") from error
    if not isinstance(value, Mapping):
        raise OpenTargetsError(f"{context} root must be an object")
    return value


def _parse_integrity(path: Path) -> Mapping[str, str]:
    checksums: dict[str, str] = {}
    try:
        lines = path.read_text(encoding="utf-8").splitlines()
    except (OSError, UnicodeDecodeError) as error:
        raise OpenTargetsError(f"cannot read Open Targets integrity index: {error}") from error
    for number, line in enumerate(lines, start=1):
        match = INTEGRITY_LINE.fullmatch(line)
        if match is None:
            raise OpenTargetsError(f"malformed integrity index line {number}")
        checksum, raw_path = match.groups()
        relative = raw_path.removeprefix("./")
        _safe_path(relative, f"integrity index line {number}")
        if relative in checksums:
            raise OpenTargetsError(f"duplicate integrity index path: {relative}")
        checksums[relative] = checksum
    if not checksums:
        raise OpenTargetsError("Open Targets integrity index is empty")
    return checksums


def _quote_identifier(value: str) -> str:
    if not IDENTIFIER.fullmatch(value):
        raise OpenTargetsError(f"invalid Croissant field name: {value}")
    return f'"{value}"'


def _nested_expression(
    parent: str,
    remaining: Sequence[tuple[str, bool]],
    parent_repeated: bool,
    depth: int,
) -> str:
    if not remaining:
        return parent
    name, repeated = remaining[0]
    if parent_repeated:
        variable = f"x{depth}"
        child = f"{variable}.{_quote_identifier(name)}"
        body = _nested_expression(child, remaining[1:], repeated, depth + 1)
        return f"list_transform({parent}, {variable} -> {body})"
    child = f"({parent}).{_quote_identifier(name)}"
    return _nested_expression(child, remaining[1:], repeated, depth + 1)


def _field_models(dataset: str, fields: Any) -> tuple[FieldModel, ...]:
    if not isinstance(fields, list) or not fields:
        raise OpenTargetsError(f"Croissant record set {dataset} has no fields")
    models: list[FieldModel] = []

    def visit(field: Any, parents: tuple[tuple[str, bool], ...]) -> None:
        if not isinstance(field, Mapping):
            raise OpenTargetsError(f"Croissant field in {dataset} must be an object")
        name = _required_string(field, "name", f"Croissant field in {dataset}")
        field_id = _required_string(field, "@id", f"Croissant field in {dataset}")
        repeated = field.get("repeated", False)
        if not isinstance(repeated, bool):
            raise OpenTargetsError(f"Croissant repeated flag is invalid: {field_id}")
        segments = parents + ((name, repeated),)
        expected_id = f"{dataset}/" + "/".join(item[0] for item in segments)
        if field_id != expected_id:
            raise OpenTargetsError(f"Croissant field ID is inconsistent: {field_id}")
        expression = _quote_identifier(segments[0][0])
        if len(segments) > 1:
            expression = _nested_expression(
                expression, segments[1:], segments[0][1], 1
            )
        data_type = field.get("dataType")
        if data_type is not None and data_type not in TYPE_NAMES:
            raise OpenTargetsError(f"unsupported Croissant type for {field_id}: {data_type}")
        models.append(
            FieldModel(
                field_id=field_id,
                expression=expression,
                data_type=data_type,
                repeated_depth=sum(item[1] for item in segments),
                nested=len(segments) > 1,
            )
        )
        children = field.get("subField", [])
        if not isinstance(children, list):
            raise OpenTargetsError(f"Croissant subField is invalid: {field_id}")
        for child in children:
            visit(child, segments)

    for item in fields:
        visit(item, ())
    ids = [model.field_id for model in models]
    if len(ids) != len(set(ids)):
        raise OpenTargetsError(f"Croissant record set has duplicate fields: {dataset}")
    return tuple(models)


def _reference_target(field: Mapping[str, Any]) -> str | None:
    reference = field.get("references")
    if reference is None:
        return None
    if not isinstance(reference, Mapping):
        raise OpenTargetsError("Croissant field reference must be an object")
    target = reference.get("field")
    if not isinstance(target, Mapping):
        raise OpenTargetsError("Croissant field reference lacks field")
    return _required_string(target, "@id", "Croissant reference field")


def _references(fields: Any) -> tuple[tuple[str, str], ...]:
    joins: list[tuple[str, str]] = []

    def visit(field: Any) -> None:
        if not isinstance(field, Mapping):
            raise OpenTargetsError("Croissant field must be an object")
        target = _reference_target(field)
        if target is not None:
            joins.append((_required_string(field, "@id", "Croissant field"), target))
        children = field.get("subField", [])
        if not isinstance(children, list):
            raise OpenTargetsError("Croissant subField must be an array")
        for child in children:
            visit(child)

    if not isinstance(fields, list):
        raise OpenTargetsError("Croissant fields must be an array")
    for item in fields:
        visit(item)
    return tuple(joins)


def _release_model(
    release: Path, contract: OpenTargetsContract
) -> tuple[tuple[ReleaseArtifact, ...], tuple[DatasetModel, ...], tuple[tuple[str, str], ...]]:
    croissant = _json_object(release / "croissant.json", "Open Targets Croissant metadata")
    if (
        croissant.get("version") != contract.release
        or croissant.get("datePublished") != contract.croissant_date_published
        or croissant.get("license") != contract.croissant_license
    ):
        raise OpenTargetsError("Croissant release identity differs from the contract")
    record_sets_value = croissant.get("recordSet")
    distributions_value = croissant.get("distribution")
    if not isinstance(record_sets_value, list) or not isinstance(distributions_value, list):
        raise OpenTargetsError("Croissant recordSet and distribution must be arrays")
    record_sets = {
        item.get("name"): item
        for item in record_sets_value
        if isinstance(item, Mapping) and isinstance(item.get("name"), str)
    }
    distributions = {
        item.get("@id"): item
        for item in distributions_value
        if isinstance(item, Mapping) and isinstance(item.get("@id"), str)
    }
    selected = {dataset.name for dataset in contract.datasets}
    if not selected.issubset(record_sets):
        missing = sorted(selected - set(record_sets))
        raise OpenTargetsError(f"Croissant lacks selected record sets: {', '.join(missing)}")

    integrity_document = contract.metadata_document("release_data_integrity")
    checksum_text = (release / "release_data_integrity.sha1").read_text(
        encoding="ascii"
    )
    expected_checksum_text = (
        f"{integrity_document.sha1}  release_data_integrity\n"
    )
    if checksum_text != expected_checksum_text:
        raise OpenTargetsError("publisher integrity checksum file has unexpected contents")
    integrity = _parse_integrity(release / "release_data_integrity")

    artifacts: list[ReleaseArtifact] = [
        ReleaseArtifact(
            dataset=None,
            relative_path=document.name,
            size_bytes=document.size_bytes,
            sha1=document.sha1,
            sha256=document.sha256,
        )
        for document in contract.metadata
    ]
    models: list[DatasetModel] = []
    joins: list[tuple[str, str]] = []
    all_field_ids: set[str] = set()
    for dataset in contract.datasets:
        record_set = record_sets[dataset.name]
        if record_set.get("@id") != dataset.name:
            raise OpenTargetsError(f"Croissant record set ID is invalid: {dataset.name}")
        fields = record_set.get("field")
        field_models = _field_models(dataset.name, fields)
        field_ids = {field.field_id for field in field_models}
        all_field_ids.update(field_ids)
        key_value = record_set.get("key")
        if isinstance(key_value, Mapping):
            key_items = [key_value]
        elif isinstance(key_value, list) and key_value:
            key_items = key_value
        else:
            raise OpenTargetsError(f"Croissant record set lacks a key: {dataset.name}")
        key_fields = tuple(
            _required_string(item, "@id", f"Croissant key in {dataset.name}")
            for item in key_items
            if isinstance(item, Mapping)
        )
        if len(key_fields) != len(key_items) or not set(key_fields).issubset(field_ids):
            raise OpenTargetsError(f"Croissant key is invalid: {dataset.name}")
        models.append(
            DatasetModel(
                name=dataset.name,
                fields=field_models,
                top_level_fields=tuple(
                    _required_string(item, "name", f"Croissant field in {dataset.name}")
                    for item in fields
                    if isinstance(item, Mapping)
                ),
                key_fields=key_fields,
            )
        )
        joins.extend(_references(fields))

        fileset = distributions.get(f"{dataset.name}-fileset")
        if not isinstance(fileset, Mapping):
            raise OpenTargetsError(f"Croissant lacks FileSet for {dataset.name}")
        includes = _required_string(fileset, "includes", f"FileSet {dataset.name}")
        expected_pattern = f"{dataset.name}/*.parquet"
        if (
            fileset.get("encodingFormat") != "application/vnd.apache.parquet"
            or includes != expected_pattern
        ):
            raise OpenTargetsError(f"Croissant FileSet is invalid: {dataset.name}")
        paths = sorted(
            path
            for path in integrity
            if path.startswith("output/")
            and fnmatchcase(path.removeprefix("output/"), includes)
        )
        if len(paths) != dataset.partition_count:
            raise OpenTargetsError(f"partition count differs for {dataset.name}")
        for path, size in zip(paths, dataset.partition_size_bytes, strict=True):
            artifacts.append(
                ReleaseArtifact(
                    dataset=dataset.name,
                    relative_path=path,
                    size_bytes=size,
                    sha1=integrity[path],
                )
            )

    declared = set(joins)
    if len(declared) != len(joins):
        raise OpenTargetsError("Croissant contains duplicate selected joins")
    for source, target in joins:
        if source not in all_field_ids or target not in all_field_ids:
            raise OpenTargetsError(f"Croissant join leaves the allowlist: {source} -> {target}")
    for join in contract.semantic_joins:
        if join.source not in all_field_ids or join.target not in all_field_ids:
            raise OpenTargetsError(
                f"semantic join references an unavailable field: {join.source} -> {join.target}"
            )
        if (join.source, join.target) in declared:
            raise OpenTargetsError(f"semantic join duplicates Croissant: {join.source}")
    return tuple(artifacts), tuple(models), tuple(sorted(joins))


def _verify_artifacts(release: Path, artifacts: Sequence[ReleaseArtifact]) -> tuple[ReleaseArtifact, ...]:
    verified: list[ReleaseArtifact] = []
    for artifact in artifacts:
        path = release / artifact.relative_path
        sha256, sha1, size = _identity(path)
        if (
            size != artifact.size_bytes
            or sha1 != artifact.sha1
            or (artifact.sha256 is not None and sha256 != artifact.sha256)
        ):
            raise OpenTargetsError(f"artifact identity mismatch: {artifact.relative_path}")
        verified.append(
            ReleaseArtifact(
                dataset=artifact.dataset,
                relative_path=artifact.relative_path,
                size_bytes=size,
                sha1=sha1,
                sha256=sha256,
            )
        )
    return tuple(verified)


def _type_shape(value: str) -> tuple[str, int]:
    result = value
    depth = 0
    while result.endswith("[]"):
        result = result[:-2]
        depth += 1
    return result, depth


def _validate_parquet(
    release: Path,
    contract: OpenTargetsContract,
    artifacts: Sequence[ReleaseArtifact],
    models: Sequence[DatasetModel],
) -> tuple[DatasetValidation, ...]:
    try:
        import duckdb
    except ImportError as error:
        raise OpenTargetsError("DuckDB is required to validate Open Targets Parquet") from error

    by_dataset: dict[str, list[Path]] = {dataset.name: [] for dataset in contract.datasets}
    for artifact in artifacts:
        if artifact.dataset is not None:
            by_dataset[artifact.dataset].append(release / artifact.relative_path)
    model_by_name = {model.name: model for model in models}
    validations: list[DatasetValidation] = []
    connection = duckdb.connect(database=":memory:")
    try:
        for dataset in contract.datasets:
            paths = sorted(by_dataset[dataset.name])
            model = model_by_name[dataset.name]
            baseline: tuple[tuple[str, ...], tuple[str, ...]] | None = None
            rows = 0
            for path in paths:
                relation = connection.from_parquet(str(path))
                schema = (
                    tuple(relation.columns),
                    tuple(str(value) for value in relation.types),
                )
                if baseline is None:
                    baseline = schema
                elif schema != baseline:
                    raise OpenTargetsError(
                        f"Parquet partition schema drift in {dataset.name}: {path.name}"
                    )
                rows += int(relation.aggregate("count(*)").fetchone()[0])
            if baseline is None or set(baseline[0]) != set(model.top_level_fields):
                raise OpenTargetsError(f"Parquet top-level schema differs for {dataset.name}")
            relation = connection.from_parquet(
                [str(path) for path in paths], union_by_name=True
            )
            for field in model.fields:
                try:
                    observed = str(relation.project(field.expression).types[0])
                except Exception as error:
                    raise OpenTargetsError(
                        f"Parquet lacks Croissant field {field.field_id}: {error}"
                    ) from error
                observed_type, observed_depth = _type_shape(observed)
                if observed_depth != field.repeated_depth:
                    raise OpenTargetsError(
                        f"Parquet cardinality differs for {field.field_id}: "
                        f"expected list depth {field.repeated_depth}, observed {observed_depth}"
                    )
                if field.data_type is None:
                    continue
                expected = TYPE_NAMES[field.data_type]
                if observed_type != expected:
                    raise OpenTargetsError(
                        f"Parquet type differs for {field.field_id}: "
                        f"expected {expected}, observed {observed}"
                    )
            if rows != dataset.row_count:
                raise OpenTargetsError(
                    f"row count differs for {dataset.name}: expected {dataset.row_count}, observed {rows}"
                )
            observed_size = sum(path.stat().st_size for path in paths)
            if observed_size != dataset.size_bytes:
                raise OpenTargetsError(f"dataset byte size differs for {dataset.name}")
            validations.append(
                DatasetValidation(
                    name=dataset.name,
                    partitions=len(paths),
                    size_bytes=observed_size,
                    row_count=rows,
                    key_fields=model.key_fields,
                    nested_fields=sum(field.nested for field in model.fields),
                )
            )
    finally:
        connection.close()
    return tuple(validations)


def _expected_files(artifacts: Sequence[ReleaseArtifact]) -> set[str]:
    return {artifact.relative_path for artifact in artifacts}


def _verify_release_tree(release: Path, artifacts: Sequence[ReleaseArtifact]) -> None:
    if release.is_symlink() or not release.is_dir():
        raise OpenTargetsError(f"Open Targets release is not a regular directory: {release}")
    observed: set[str] = set()
    expected = _expected_files(artifacts)
    expected_directories = {
        parent.as_posix()
        for value in expected
        for parent in PurePosixPath(value).parents
        if parent != PurePosixPath(".")
    }
    try:
        for path in release.rglob("*"):
            if path.is_symlink():
                raise OpenTargetsError(f"Open Targets release contains a symlink: {path}")
            if path.is_file():
                observed.add(path.relative_to(release).as_posix())
            elif path.is_dir():
                relative = path.relative_to(release).as_posix()
                if relative not in expected_directories:
                    raise OpenTargetsError(
                        f"Open Targets release contains an unexpected directory: {path}"
                    )
            else:
                raise OpenTargetsError(
                    f"Open Targets release contains a special filesystem entry: {path}"
                )
    except OSError as error:
        raise OpenTargetsError(f"cannot inspect Open Targets release: {error}") from error
    if observed != expected:
        raise OpenTargetsError("Open Targets release is incomplete or contaminated")


def validate_open_targets_release(
    contract: OpenTargetsContract, release: Path
) -> tuple[OpenTargetsValidation, tuple[ReleaseArtifact, ...]]:
    artifacts, models, joins = _release_model(release, contract)
    _verify_release_tree(release, artifacts)
    verified = _verify_artifacts(release, artifacts)
    datasets = _validate_parquet(release, contract, verified, models)
    validation = OpenTargetsValidation(
        release=contract.release,
        datasets=datasets,
        partition_count=sum(item.partitions for item in datasets),
        size_bytes=sum(item.size_bytes for item in datasets),
        row_count=sum(item.row_count for item in datasets),
        declared_joins=len(joins),
        semantic_joins=len(contract.semantic_joins),
    )
    return validation, verified


def _manifest_path(
    root: Path,
    contract: OpenTargetsContract,
    contract_sha256: str | None = None,
) -> Path:
    legacy = root / f"manifests/open-targets/{contract.resource}-{contract.release}.json"
    if contract_sha256 is None or not legacy.exists():
        return legacy
    try:
        document = json.loads(legacy.read_text(encoding="utf-8"))
        if (
            document.get("schema") == "genosapiens.public-resource"
            and document.get("transformation", {}).get("parameters_sha256")
            == contract_sha256
        ):
            return legacy
    except (OSError, json.JSONDecodeError, AttributeError):
        pass
    return legacy.with_name(
        f"{contract.resource}-{contract.release}.{contract_sha256[:12]}.json"
    )


def _load_manifest(path: Path) -> PublicResourceManifest:
    value = _json_object(path, "Open Targets manifest")
    try:
        return PublicResourceManifest.from_mapping(value)
    except ManifestError as error:
        raise OpenTargetsError(f"invalid Open Targets manifest: {error}") from error


def _manifest(
    contract: OpenTargetsContract,
    contract_sha256: str,
    image_digest: str,
    retrieved_on: str,
    artifacts: Sequence[ReleaseArtifact],
) -> PublicResourceManifest:
    croissant = contract.metadata_document("croissant.json")
    return PublicResourceManifest(
        resource_id=contract.resource,
        release=contract.release,
        source=SourceIdentity(
            url=f"{contract.base_url}/croissant.json",
            sha256=croissant.sha256,
        ),
        license=contract.license,
        assembly=contract.assembly,
        retrieved_on=retrieved_on,
        transformation=Transformation(TRANSFORMATION_SCHEMA, contract_sha256),
        image_digest=image_digest,
        outputs=tuple(
            ResourceOutput(
                path=f"{contract.output_prefix}/{artifact.relative_path}",
                size_bytes=artifact.size_bytes,
                sha256=artifact.sha256 or "",
                retention=RetentionPolicy("durable", None),
            )
            for artifact in artifacts
        ),
    )


def _verify_manifest_contract(
    manifest: PublicResourceManifest,
    contract: OpenTargetsContract,
    contract_sha256: str,
    artifacts: Sequence[ReleaseArtifact],
) -> None:
    croissant = contract.metadata_document("croissant.json")
    expected_paths = {
        f"{contract.output_prefix}/{artifact.relative_path}" for artifact in artifacts
    }
    if (
        manifest.resource_id != contract.resource
        or manifest.release != contract.release
        or manifest.assembly != contract.assembly
        or manifest.license != contract.license
        or manifest.source.url != f"{contract.base_url}/croissant.json"
        or manifest.source.sha256 != croissant.sha256
        or manifest.transformation
        != Transformation(TRANSFORMATION_SCHEMA, contract_sha256)
        or {output.path for output in manifest.outputs} != expected_paths
    ):
        raise OpenTargetsError("Open Targets manifest identity differs from its contract")


def _require_within_root(root: Path, path: Path, context: str) -> None:
    if not path.is_absolute() or not path.is_relative_to(root):
        raise OpenTargetsError(f"{context} escapes reference_root: {path}")
    current = root
    for part in path.relative_to(root).parts:
        current /= part
        if current.is_symlink():
            raise OpenTargetsError(f"{context} contains a symlink: {current}")


def _download_metadata(
    contract: OpenTargetsContract,
    destination: Path,
    root: Path,
    downloader: Downloader,
) -> None:
    for document in contract.metadata:
        downloader(
            DownloadRequest(
                url=f"{contract.base_url}/{document.name}",
                destination=destination / document.name,
                expected_size_bytes=document.size_bytes,
                expected_sha256=document.sha256,
                expected_sha1=document.sha1,
                attempts=3,
                timeout_seconds=300,
            ),
            root,
        )


def _download_partitions(
    contract: OpenTargetsContract,
    destination: Path,
    root: Path,
    downloader: Downloader,
    artifacts: Sequence[ReleaseArtifact],
) -> None:
    for artifact in artifacts:
        if artifact.dataset is None:
            continue
        downloader(
            DownloadRequest(
                url=f"{contract.base_url}/{artifact.relative_path}",
                destination=destination / artifact.relative_path,
                expected_size_bytes=artifact.size_bytes,
                expected_sha1=artifact.sha1,
                attempts=3,
                timeout_seconds=300,
            ),
            root,
        )


def _remove_partial_directories(root: Path) -> None:
    for path in sorted(root.rglob(".partial"), reverse=True):
        try:
            path.rmdir()
        except OSError as error:
            raise OpenTargetsError(f"Open Targets partial directory is not empty: {path}") from error


def _make_read_only(root: Path) -> None:
    for path in root.rglob("*"):
        os.chmod(path, 0o755 if path.is_dir() else 0o444)
    os.chmod(root, 0o755)


def bootstrap_open_targets(
    *,
    contract_path: str | Path,
    reference_root: str | Path,
    image_digest: str,
    retrieved_on: str | None = None,
    scratch_root: str | Path | None = None,
    reserve_bytes: int = DEFAULT_RESERVE_BYTES,
    downloader: Downloader = download_atomic,
    capacity_checker: CapacityChecker = require_capacity,
) -> OpenTargetsBootstrapResult:
    contract_file = Path(contract_path)
    root = Path(reference_root).resolve(strict=True)
    contract = load_open_targets_contract(contract_file)
    if not DIGEST.fullmatch(image_digest):
        raise OpenTargetsError("Open Targets image_digest must be a SHA-256 digest")
    contract_sha256 = _identity(contract_file)[0]
    release = root / contract.output_prefix
    manifest_path = _manifest_path(root, contract, contract_sha256)
    _require_within_root(root, release, "Open Targets release")
    _require_within_root(root, manifest_path, "Open Targets manifest")

    if manifest_path.exists() or manifest_path.is_symlink():
        artifacts, _models, _joins = _release_model(release, contract)
        manifest = _load_manifest(manifest_path)
        _verify_manifest_contract(manifest, contract, contract_sha256, artifacts)
        try:
            verify_outputs(manifest, root)
        except ManifestError as error:
            raise OpenTargetsError(f"completed Open Targets bundle is invalid: {error}") from error
        validation, _verified = validate_open_targets_release(contract, release)
        return OpenTargetsBootstrapResult(
            resource=contract.resource,
            release=contract.release,
            manifest_path=str(manifest_path),
            manifest_sha256=_identity(manifest_path)[0],
            validation=validation,
            reused=True,
        )

    if release.exists() or release.is_symlink():
        validation, artifacts = validate_open_targets_release(contract, release)
        publication = publish_manifest(
            _manifest(
                contract,
                contract_sha256,
                image_digest,
                retrieved_on or date.today().isoformat(),
                artifacts,
            ),
            manifest_path,
            root,
        )
        return OpenTargetsBootstrapResult(
            resource=contract.resource,
            release=contract.release,
            manifest_path=str(publication.path),
            manifest_sha256=publication.manifest_sha256,
            validation=validation,
            reused=False,
        )

    scratch_base = (
        Path(scratch_root) if scratch_root is not None else root / "scratch/public-build"
    )
    if not scratch_base.is_absolute():
        raise OpenTargetsError("Open Targets scratch_root must be absolute")
    _require_within_root(root, scratch_base, "Open Targets scratch_root")
    scratch = scratch_base / "open-targets" / contract.release
    _require_within_root(root, scratch, "Open Targets scratch")
    scratch.mkdir(parents=True, exist_ok=True, mode=0o700)
    total_source_bytes = sum(item.size_bytes for item in contract.metadata) + sum(
        item.size_bytes for item in contract.datasets
    )
    existing_bytes = sum(
        path.stat().st_size for path in scratch.rglob("*") if path.is_file()
    )
    capacity_checker(
        root,
        CapacityPlan(
            operation=f"bootstrap Open Targets Platform {contract.release}",
            source_bytes=max(total_source_bytes - existing_bytes, 0),
            reserve_bytes=reserve_bytes,
        ),
    )
    _download_metadata(contract, scratch, root, downloader)
    artifacts, _models, _joins = _release_model(scratch, contract)
    _download_partitions(contract, scratch, root, downloader, artifacts)
    _remove_partial_directories(scratch)
    validation, verified = validate_open_targets_release(contract, scratch)
    _make_read_only(scratch)
    release.parent.mkdir(parents=True, exist_ok=True)
    try:
        os.replace(scratch, release)
    except OSError:
        if not release.exists():
            raise
        shutil.rmtree(scratch, ignore_errors=True)
        validation, verified = validate_open_targets_release(contract, release)
    directory_fd = os.open(release.parent, os.O_RDONLY)
    try:
        os.fsync(directory_fd)
    finally:
        os.close(directory_fd)
    publication = publish_manifest(
        _manifest(
            contract,
            contract_sha256,
            image_digest,
            retrieved_on or date.today().isoformat(),
            verified,
        ),
        manifest_path,
        root,
    )
    return OpenTargetsBootstrapResult(
        resource=contract.resource,
        release=contract.release,
        manifest_path=str(publication.path),
        manifest_sha256=publication.manifest_sha256,
        validation=validation,
        reused=publication.reused,
    )


def open_targets_result_as_json(result: OpenTargetsBootstrapResult) -> str:
    return json.dumps(asdict(result), indent=2, sort_keys=True) + "\n"
