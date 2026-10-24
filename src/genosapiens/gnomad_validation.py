"""Official-source concordance validation for the staged gnomAD pilot."""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from dataclasses import asdict, dataclass
import json
import math
import os
from pathlib import Path
import re
import struct
from typing import Any

from genosapiens.gnomad_hail import (
    GnomadHailError,
    HAIL_VERSION,
    PUBLIC_GNOMAD_SOURCE,
    PYSPARK_VERSION,
    RUN_ID,
    _validate_runtime_confinement,
    initialize_local_hail,
    prepare_hail_scratch,
)
from genosapiens.gnomad_pilot import (
    BCFTOOLS_VERSION,
    COMPLETION_FILE,
    PILOT_CHROMOSOMES,
    PILOT_COMPLETION_SCHEMA,
    PILOT_SCHEMA,
    SHA256_DIGEST,
    _fsync_directory,
    _pilot_run_lock,
    _remove_failed_run,
    _require_bcftools,
    _run_bcftools,
    _sha256_file,
    _source_inventory,
    _verify_source_root_metadata,
)
from genosapiens.gnomad_projection import (
    DATASET_VCF_PREFIXES,
    GnomadProjectionContract,
    GnomadProjectionLayout,
    gnomad_vcf_metadata,
    load_gnomad_projection_contract,
    project_gnomad_joint_row,
    resolve_gnomad_projection_layout,
)


CONTRACT_SCHEMA = "genosapiens.gnomad-pilot-validation-contract"
VALIDATION_SCHEMA = "genosapiens.gnomad-pilot-validation"
VALIDATION_COMPLETION_SCHEMA = "genosapiens.gnomad-pilot-validation-completion"
VALIDATION_REPORT_FILE = "gnomad-pilot-validation-report.json"
VALIDATION_COMPLETION_FILE = "gnomad-pilot-validation-complete.json"
INFO_HEADER = re.compile(
    r'^##INFO=<ID=([^,]+),Number=([^,]+),Type=([^,]+),Description="(.*)">$'
)
ALLELE = re.compile(r"^[ACGT]+$")
SUPPORTED_EXPECTATIONS = frozenset(
    {
        "ac-an-derived-frequency",
        "ancestry-values",
        "genome-mid-faf-missing",
        "exome-data-missing",
        "exome-filters-missing",
        "exome-filters-empty",
        "exome-filters-nonempty",
        "genome-filters-empty",
        "genome-filters-nonempty",
        "discordance-missing",
        "discordance-false",
        "discordance-true",
    }
)
REQUIRED_EXPECTATIONS = SUPPORTED_EXPECTATIONS
PILOT_REPORT_FIELDS = {
    "schema",
    "mode",
    "input_classification",
    "source",
    "projection_id",
    "projection_contract_sha256",
    "source_root_metadata_url",
    "source_root_metadata_object_generation",
    "source_root_metadata_md5_base64",
    "source_root_metadata_verified",
    "source_object_inventory_sha256",
    "source_object_count",
    "source_object_bytes",
    "source_object_inventory_verified",
    "source_written_hail_version",
    "source_partitions",
    "chromosomes",
    "measurements",
    "measured",
    "extrapolation",
    "gate",
    "hail_version",
    "pyspark_version",
    "bcftools_version",
    "image_digest",
    "workers",
    "spark_master",
    "output_format",
    "remote_bytes_method",
    "memory_method",
    "scratch_method",
    "scratch_root",
    "run_root",
    "report_path",
}
MEASUREMENT_FIELDS = {
    "chromosome",
    "rows",
    "remote_bytes",
    "runtime_ns",
    "bcf_bytes",
    "csi_bytes",
    "output_bytes",
    "bcf_path",
    "csi_path",
    "bcf_sha256",
    "csi_sha256",
}
SEMANTIC_INFO_FIELDS = (
    "EX_AC",
    "EX_AN",
    "EX_NHOMALT",
    "EX_NHEMI",
    "EX_AC_ANC",
    "EX_AN_ANC",
    "EX_NHOMALT_ANC",
    "EX_NHEMI_ANC",
    "EX_FILTERS",
    "EX_FILTERS_DEFINED",
    "GN_FILTERS",
    "GN_FILTERS_DEFINED",
    "GN_FAF95_ANC",
    "GN_FAF99_ANC",
    "REGION_NOT_CALLED_IN_EXOMES",
    "REGION_NOT_CALLED_IN_GENOMES",
    "DISC_UNION_P",
    "DISC_UNION_TEST",
    "DISC_UNION_ANC",
    "DISC",
)


class GnomadValidationError(GnomadHailError):
    """Raised when pilot evidence differs from its pinned official source."""


def _mapping(value: Any, context: str) -> Mapping[str, Any]:
    if isinstance(value, Mapping):
        return value
    try:
        converted = dict(value)
    except (TypeError, ValueError) as error:
        raise GnomadValidationError(f"{context} must be an object") from error
    if not all(isinstance(key, str) for key in converted):
        raise GnomadValidationError(f"{context} keys must be strings")
    return converted


def _exact_fields(data: Mapping[str, Any], expected: set[str], context: str) -> None:
    if set(data) != expected:
        raise GnomadValidationError(f"{context} fields are invalid")


def _string(data: Mapping[str, Any], name: str, context: str) -> str:
    value = data.get(name)
    if not isinstance(value, str) or not value:
        raise GnomadValidationError(f"{context}.{name} must be a non-empty string")
    return value


def _sha256(data: Mapping[str, Any], name: str, context: str) -> str:
    value = _string(data, name, context)
    if len(value) != 64 or any(
        character not in "0123456789abcdef" for character in value
    ):
        raise GnomadValidationError(f"{context}.{name} must be a SHA-256 digest")
    return value


def _strings(value: Any, context: str) -> tuple[str, ...]:
    if (
        not isinstance(value, list)
        or not value
        or any(not isinstance(item, str) or not item for item in value)
        or len(value) != len(set(value))
    ):
        raise GnomadValidationError(f"{context} must be a unique string array")
    return tuple(value)


@dataclass(frozen=True, slots=True)
class GnomadSentinel:
    identifier: str
    chromosome: str
    position: int
    reference: str
    alternate: str
    expectations: tuple[str, ...]

    @classmethod
    def from_mapping(cls, value: Mapping[str, Any]) -> "GnomadSentinel":
        _exact_fields(
            value,
            {"id", "chromosome", "position", "reference", "alternate", "expectations"},
            "sentinel",
        )
        position = value.get("position")
        if not isinstance(position, int) or isinstance(position, bool) or position <= 0:
            raise GnomadValidationError("sentinel.position must be a positive integer")
        sentinel = cls(
            identifier=_string(value, "id", "sentinel"),
            chromosome=_string(value, "chromosome", "sentinel"),
            position=position,
            reference=_string(value, "reference", "sentinel"),
            alternate=_string(value, "alternate", "sentinel"),
            expectations=_strings(value.get("expectations"), "sentinel.expectations"),
        )
        if (
            not RUN_ID.fullmatch(sentinel.identifier)
            or sentinel.chromosome not in PILOT_CHROMOSOMES
            or not ALLELE.fullmatch(sentinel.reference)
            or not ALLELE.fullmatch(sentinel.alternate)
            or sentinel.reference == sentinel.alternate
            or not set(sentinel.expectations).issubset(SUPPORTED_EXPECTATIONS)
        ):
            raise GnomadValidationError("sentinel identity or expectations are invalid")
        return sentinel


@dataclass(frozen=True, slots=True)
class GnomadValidationContract:
    validation_id: str
    pilot_run_id: str
    projection_id: str
    projection_contract_sha256: str
    pilot_report_sha256: str
    source_object_inventory_sha256: str
    ancestry_order: tuple[str, ...]
    sentinels: tuple[GnomadSentinel, ...]


def load_gnomad_validation_contract(path: str | Path) -> GnomadValidationContract:
    try:
        value = json.loads(Path(path).read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as error:
        raise GnomadValidationError(
            f"cannot load gnomAD validation contract: {error}"
        ) from error
    data = _mapping(value, "validation contract")
    _exact_fields(
        data,
        {
            "schema",
            "validation_id",
            "pilot_run_id",
            "projection_id",
            "projection_contract_sha256",
            "pilot_report_sha256",
            "source_object_inventory_sha256",
            "ancestry_order",
            "sentinels",
        },
        "validation contract",
    )
    if data.get("schema") != CONTRACT_SCHEMA:
        raise GnomadValidationError(
            f"validation contract schema must be {CONTRACT_SCHEMA!r}"
        )
    sentinel_values = data.get("sentinels")
    if not isinstance(sentinel_values, list) or not sentinel_values:
        raise GnomadValidationError("validation contract sentinels must be an array")
    sentinels = tuple(
        GnomadSentinel.from_mapping(_mapping(item, "sentinel"))
        for item in sentinel_values
    )
    keys = [
        (item.chromosome, item.position, item.reference, item.alternate)
        for item in sentinels
    ]
    identifiers = [item.identifier for item in sentinels]
    coverage = set().union(*(set(item.expectations) for item in sentinels))
    if (
        len(keys) != len(set(keys))
        or len(identifiers) != len(set(identifiers))
        or set(item.chromosome for item in sentinels) != set(PILOT_CHROMOSOMES)
        or coverage != REQUIRED_EXPECTATIONS
    ):
        raise GnomadValidationError(
            "validation sentinels must be unique and cover both shards and all semantics"
        )
    contract = GnomadValidationContract(
        validation_id=_string(data, "validation_id", "validation contract"),
        pilot_run_id=_string(data, "pilot_run_id", "validation contract"),
        projection_id=_string(data, "projection_id", "validation contract"),
        projection_contract_sha256=_sha256(
            data, "projection_contract_sha256", "validation contract"
        ),
        pilot_report_sha256=_sha256(
            data, "pilot_report_sha256", "validation contract"
        ),
        source_object_inventory_sha256=_sha256(
            data, "source_object_inventory_sha256", "validation contract"
        ),
        ancestry_order=_strings(data.get("ancestry_order"), "ancestry_order"),
        sentinels=sentinels,
    )
    if not RUN_ID.fullmatch(contract.validation_id) or not RUN_ID.fullmatch(
        contract.pilot_run_id
    ):
        raise GnomadValidationError("validation identifiers are not path-safe")
    return contract


@dataclass(frozen=True, slots=True)
class GnomadShardValidation:
    chromosome: str
    rows: int
    bcf_path: str
    bcf_bytes: int
    bcf_sha256: str
    csi_path: str
    csi_bytes: int
    csi_sha256: str
    header_info_fields: int
    header_verified: bool
    index_row_count_verified: bool


@dataclass(frozen=True, slots=True)
class GnomadFrequencyValidation:
    dataset: str
    ancestry: str | None
    ac: int | None
    an: int | None
    derived_af: float | None
    source_af: float | None
    comparable: bool


@dataclass(frozen=True, slots=True)
class GnomadSentinelValidation:
    identifier: str
    chromosome: str
    position: int
    reference: str
    alternate: str
    expectations: tuple[str, ...]
    expectation_checks: Mapping[str, bool]
    bcf_filter: str
    csi_query_matches: int
    info_fields_checked: int
    semantic_observations: Mapping[str, Any]
    frequencies: tuple[GnomadFrequencyValidation, ...]


@dataclass(frozen=True, slots=True)
class GnomadPilotValidationReport:
    schema: str
    mode: str
    input_classification: str
    accepted: bool
    mismatch_count: int
    validation_id: str
    validation_contract_sha256: str
    projection_id: str
    projection_contract_sha256: str
    source: str
    source_root_metadata_object_generation: str
    source_root_metadata_md5_base64: str
    source_object_inventory_sha256: str
    source_object_count: int
    source_object_bytes: int
    pilot_report_path: str
    pilot_report_sha256: str
    pilot_image_digest: str
    validator_image_digest: str
    ancestry_order: tuple[str, ...]
    shards: tuple[GnomadShardValidation, ...]
    sentinels: tuple[GnomadSentinelValidation, ...]
    hail_version: str
    pyspark_version: str
    bcftools_version: str
    workers: int
    spark_master: str
    scratch_root: str
    report_path: str


def _json_safe(value: Any) -> Any:
    if isinstance(value, float) and not math.isfinite(value):
        if math.isnan(value):
            return "NaN"
        return "Infinity" if value > 0 else "-Infinity"
    if isinstance(value, Mapping):
        return {str(key): _json_safe(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [_json_safe(item) for item in value]
    return value


def validation_report_as_json(report: GnomadPilotValidationReport) -> str:
    return json.dumps(_json_safe(asdict(report)), indent=2, sort_keys=True) + "\n"


def _load_json_object(path: Path, context: str) -> Mapping[str, Any]:
    try:
        return _mapping(json.loads(path.read_text(encoding="utf-8")), context)
    except (OSError, json.JSONDecodeError) as error:
        raise GnomadValidationError(f"cannot load {context}: {error}") from error


def _validate_pilot_report(
    pilot_root: Path,
    contract: GnomadValidationContract,
    projection: GnomadProjectionContract,
    projection_sha256: str,
) -> tuple[Mapping[str, Any], tuple[Mapping[str, Any], ...], Path]:
    completion_path = pilot_root / COMPLETION_FILE
    report_path = pilot_root / "gnomad-pilot-report.json"
    if pilot_root.is_symlink() or not pilot_root.is_dir():
        raise GnomadValidationError("completed gnomAD pilot directory is unavailable")
    completion = _load_json_object(completion_path, "pilot completion marker")
    if completion_path.is_symlink() or report_path.is_symlink():
        raise GnomadValidationError("pilot report paths cannot be symbolic links")
    _exact_fields(
        completion,
        {"schema", "report_path", "report_sha256"},
        "pilot completion",
    )
    if (
        completion.get("schema") != PILOT_COMPLETION_SCHEMA
        or completion.get("report_path") != report_path.name
        or completion.get("report_sha256") != contract.pilot_report_sha256
        or _sha256_file(report_path) != contract.pilot_report_sha256
    ):
        raise GnomadValidationError("pilot completion marker or report checksum differs")
    report = _load_json_object(report_path, "pilot report")
    _exact_fields(report, PILOT_REPORT_FIELDS, "pilot report")
    if (
        report.get("schema") != PILOT_SCHEMA
        or report.get("mode") != "chromosome-21-22-pilot"
        or report.get("input_classification") != "public"
        or report.get("source") != PUBLIC_GNOMAD_SOURCE
        or report.get("projection_id") != contract.projection_id
        or report.get("projection_id") != projection.projection_id
        or report.get("projection_contract_sha256") != projection_sha256
        or projection_sha256 != contract.projection_contract_sha256
        or report.get("source_object_inventory_sha256")
        != contract.source_object_inventory_sha256
        or report.get("source_root_metadata_object_generation")
        != projection.source.root_metadata_object_generation
        or report.get("source_root_metadata_md5_base64")
        != projection.source.root_metadata_md5_base64
        or report.get("source_root_metadata_verified") is not True
        or report.get("source_object_inventory_verified") is not True
        or report.get("source_partitions") != projection.source.expected_partitions
        or tuple(report.get("chromosomes", ())) != PILOT_CHROMOSOMES
        or report.get("hail_version") != HAIL_VERSION
        or report.get("pyspark_version") != PYSPARK_VERSION
        or report.get("bcftools_version") != BCFTOOLS_VERSION
        or report.get("run_root") != str(pilot_root)
        or report.get("report_path") != str(report_path)
        or _mapping(report.get("gate"), "pilot gate").get("allowed") is not True
    ):
        raise GnomadValidationError("pilot report identity or continuation gate differs")
    measurements_value = report.get("measurements")
    if not isinstance(measurements_value, list) or len(measurements_value) != 2:
        raise GnomadValidationError("pilot report measurements are invalid")
    measurements = tuple(
        _mapping(value, "pilot measurement") for value in measurements_value
    )
    if any(set(value) != MEASUREMENT_FIELDS for value in measurements):
        raise GnomadValidationError("pilot measurement fields are invalid")
    return report, measurements, report_path


def _validate_header(header: str, projection: GnomadProjectionContract) -> int:
    observed: dict[str, tuple[str, str, str]] = {}
    columns: tuple[str, ...] | None = None
    for line in header.splitlines():
        match = INFO_HEADER.fullmatch(line)
        if match:
            identifier, number, field_type, description = match.groups()
            if identifier in observed:
                raise GnomadValidationError(f"duplicate BCF INFO declaration: {identifier}")
            observed[identifier] = (number, field_type, description)
        elif line.startswith("#CHROM"):
            columns = tuple(line.split("\t"))
    expected = gnomad_vcf_metadata(projection)["info"]
    expected_values = {
        identifier: (
            metadata["Number"],
            metadata["Type"],
            metadata["Description"],
        )
        for identifier, metadata in expected.items()
    }
    if observed != expected_values:
        raise GnomadValidationError("BCF INFO header differs from the projection contract")
    if columns != ("#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO"):
        raise GnomadValidationError("pilot BCF must be site-only with eight VCF columns")
    return len(observed)


def _positive_integer(value: Any, context: str) -> int:
    if not isinstance(value, int) or isinstance(value, bool) or value <= 0:
        raise GnomadValidationError(f"{context} must be a positive integer")
    return value


def _validate_shards(
    pilot_root: Path,
    measurements: tuple[Mapping[str, Any], ...],
    projection: GnomadProjectionContract,
) -> tuple[GnomadShardValidation, ...]:
    output_root = pilot_root / "pilot-output"
    expected_names = {
        f"gnomad.joint.v4.1.{chromosome}.bcf"
        for chromosome in PILOT_CHROMOSOMES
    }
    expected_names.update(f"{name}.csi" for name in tuple(expected_names))
    if output_root.is_symlink() or not output_root.is_dir():
        raise GnomadValidationError("pilot output directory is unavailable")
    try:
        observed_names = {path.name for path in output_root.iterdir()}
    except OSError as error:
        raise GnomadValidationError(f"cannot inspect pilot output directory: {error}") from error
    if observed_names != expected_names:
        raise GnomadValidationError("pilot output directory contains unexpected files")

    results = []
    for chromosome, measurement in zip(PILOT_CHROMOSOMES, measurements, strict=True):
        if measurement.get("chromosome") != chromosome:
            raise GnomadValidationError("pilot measurements are not chr21 then chr22")
        bcf_path = output_root / f"gnomad.joint.v4.1.{chromosome}.bcf"
        csi_path = Path(f"{bcf_path}.csi")
        if (
            bcf_path.is_symlink()
            or csi_path.is_symlink()
            or not bcf_path.is_file()
            or not csi_path.is_file()
            or measurement.get("bcf_path") != str(bcf_path)
            or measurement.get("csi_path") != str(csi_path)
            or measurement.get("bcf_bytes") != bcf_path.stat().st_size
            or measurement.get("csi_bytes") != csi_path.stat().st_size
            or measurement.get("bcf_sha256") != _sha256_file(bcf_path)
            or measurement.get("csi_sha256") != _sha256_file(csi_path)
        ):
            raise GnomadValidationError(f"{chromosome} BCF or CSI identity differs")
        rows = _positive_integer(measurement.get("rows"), f"{chromosome} rows")
        try:
            indexed_rows = int(_run_bcftools(("index", "--nrecords", str(bcf_path))))
        except ValueError as error:
            raise GnomadValidationError(
                f"bcftools returned an invalid row count for {chromosome}"
            ) from error
        if indexed_rows != rows:
            raise GnomadValidationError(f"{chromosome} CSI row count differs")
        header_fields = _validate_header(
            _run_bcftools(("view", "--header-only", str(bcf_path))), projection
        )
        results.append(
            GnomadShardValidation(
                chromosome=chromosome,
                rows=rows,
                bcf_path=str(bcf_path),
                bcf_bytes=bcf_path.stat().st_size,
                bcf_sha256=str(measurement["bcf_sha256"]),
                csi_path=str(csi_path),
                csi_bytes=csi_path.stat().st_size,
                csi_sha256=str(measurement["csi_sha256"]),
                header_info_fields=header_fields,
                header_verified=True,
                index_row_count_verified=True,
            )
        )
    return tuple(results)


def _parse_scalar(value: str, field_type: str, context: str) -> Any:
    if value == ".":
        return None
    try:
        if field_type == "Integer":
            return int(value)
        if field_type == "Float":
            return float(value)
        if field_type == "String":
            return value
    except ValueError as error:
        raise GnomadValidationError(f"invalid {field_type} value for {context}") from error
    raise GnomadValidationError(f"unsupported BCF INFO type for {context}: {field_type}")


def _parse_info(value: str, projection: GnomadProjectionContract) -> dict[str, Any]:
    metadata = gnomad_vcf_metadata(projection)["info"]
    result: dict[str, Any] = {identifier: None for identifier in metadata}
    seen: set[str] = set()
    if value != ".":
        for item in value.split(";"):
            if "=" not in item:
                raise GnomadValidationError(f"BCF INFO entry lacks a value: {item}")
            identifier, raw = item.split("=", 1)
            if identifier not in metadata or identifier in seen:
                raise GnomadValidationError(f"unknown or duplicate BCF INFO field: {identifier}")
            seen.add(identifier)
            number = metadata[identifier]["Number"]
            field_type = metadata[identifier]["Type"]
            if number == "1":
                result[identifier] = _parse_scalar(raw, field_type, identifier)
            else:
                values = [
                    _parse_scalar(element, field_type, identifier)
                    for element in raw.split(",")
                ]
                if number != "." and len(values) != int(number):
                    raise GnomadValidationError(
                        f"BCF INFO array has the wrong length: {identifier}"
                    )
                result[identifier] = values
    for prefix in ("EX", "GN"):
        defined = result[f"{prefix}_FILTERS_DEFINED"]
        filters = result[f"{prefix}_FILTERS"]
        if defined not in (0, 1) or (defined == 0 and filters is not None):
            raise GnomadValidationError(f"{prefix} filter availability is inconsistent")
        if defined == 1 and filters is None:
            result[f"{prefix}_FILTERS"] = []
    return result


def _query_sentinel(
    bcf_path: Path,
    sentinel: GnomadSentinel,
    projection: GnomadProjectionContract,
) -> tuple[str, dict[str, Any], int]:
    region = f"{sentinel.chromosome}:{sentinel.position}-{sentinel.position}"
    output = _run_bcftools(
        (
            "query",
            "--regions",
            region,
            "--format",
            "%CHROM\\t%POS\\t%REF\\t%ALT\\t%FILTER\\t%INFO\\n",
            str(bcf_path),
        )
    )
    matches = []
    for line in output.splitlines():
        fields = line.split("\t")
        if len(fields) != 6:
            raise GnomadValidationError("bcftools sentinel query returned malformed output")
        chromosome, position, reference, alternate, filter_value, info = fields
        if (
            chromosome == sentinel.chromosome
            and position == str(sentinel.position)
            and reference == sentinel.reference
            and alternate == sentinel.alternate
        ):
            matches.append((filter_value, _parse_info(info, projection)))
    if len(matches) != 1:
        raise GnomadValidationError(
            f"CSI query recovered {len(matches)} exact records for {sentinel.identifier}"
        )
    filter_value, info = matches[0]
    if filter_value != "PASS":
        raise GnomadValidationError(
            f"pilot projection FILTER is not PASS for {sentinel.identifier}"
        )
    return filter_value, info, len(matches)


def _projected_info(
    projected: Mapping[str, Any], projection: GnomadProjectionContract
) -> dict[str, Any]:
    info: dict[str, Any] = {}
    ancestries = projection.major_genetic_ancestries
    for dataset_contract in projection.datasets:
        dataset = _mapping(projected[dataset_contract.name], dataset_contract.name)
        adjusted = _mapping(dataset["adjusted"], f"{dataset_contract.name}.adjusted")
        ancestry_values = _mapping(
            dataset["ancestries"], f"{dataset_contract.name}.ancestries"
        )
        faf = _mapping(dataset["faf"], f"{dataset_contract.name}.faf")
        grpmax = _mapping(dataset["grpmax"], f"{dataset_contract.name}.grpmax")
        fafmax = _mapping(dataset["fafmax"], f"{dataset_contract.name}.fafmax")
        prefix = DATASET_VCF_PREFIXES[dataset_contract.name]

        def ancestry_array(field: str) -> list[Any]:
            return [
                _mapping(
                    ancestry_values[item], f"{dataset_contract.name}.{item}"
                ).get(field)
                for item in ancestries
            ]

        info.update(
            {
                f"{prefix}_AC": adjusted.get("ac"),
                f"{prefix}_AN": adjusted.get("an"),
                f"{prefix}_NHOMALT": adjusted.get("homozygote_count"),
                f"{prefix}_NHEMI": adjusted.get("hemizygote_count"),
                f"{prefix}_AC_ANC": ancestry_array("ac"),
                f"{prefix}_AN_ANC": ancestry_array("an"),
                f"{prefix}_NHOMALT_ANC": ancestry_array("homozygote_count"),
                f"{prefix}_NHEMI_ANC": ancestry_array("hemizygote_count"),
                f"{prefix}_FAF95": faf.get("faf95"),
                f"{prefix}_FAF99": faf.get("faf99"),
                f"{prefix}_FAF95_ANC": ancestry_array("faf95"),
                f"{prefix}_FAF99_ANC": ancestry_array("faf99"),
                f"{prefix}_GRPMAX_AC": grpmax.get("AC"),
                f"{prefix}_GRPMAX_AF": grpmax.get("AF"),
                f"{prefix}_GRPMAX_AN": grpmax.get("AN"),
                f"{prefix}_GRPMAX_NHOMALT": grpmax.get("homozygote_count"),
                f"{prefix}_GRPMAX_ANC": grpmax.get("gen_anc"),
                f"{prefix}_FAFMAX95": fafmax.get("faf95_max"),
                f"{prefix}_FAFMAX95_ANC": fafmax.get("faf95_max_gen_anc"),
                f"{prefix}_FAFMAX99": fafmax.get("faf99_max"),
                f"{prefix}_FAFMAX99_ANC": fafmax.get("faf99_max_gen_anc"),
            }
        )
        if dataset_contract.has_filters:
            filters = dataset.get("filters")
            info[f"{prefix}_FILTERS"] = filters
            info[f"{prefix}_FILTERS_DEFINED"] = 0 if filters is None else 1

    region = _mapping(projected["region_flags"], "region_flags")
    for flag in projection.region_flags:
        value = region.get(flag)
        info[f"REGION_{flag.upper()}"] = None if value is None else int(value)
    discordance = _mapping(projected["discordance"], "discordance")
    overall = _mapping(discordance["overall"], "discordance.overall")
    ancestry_discordance = _mapping(
        discordance["ancestries"], "discordance.ancestries"
    )
    cmh = _mapping(
        discordance["cochran_mantel_haenszel_test"],
        "discordance.cochran_mantel_haenszel_test",
    )
    union = _mapping(discordance["stat_union"], "discordance.stat_union")
    is_discrepant = discordance.get("is_discrepant")
    info.update(
        {
            "DISC_OVERALL_P": overall.get("p_value"),
            "DISC_OVERALL_OR": overall.get("odds_ratio"),
            "DISC_ANC_P": [
                _mapping(
                    ancestry_discordance[item], f"discordance.{item}"
                ).get("p_value")
                for item in ancestries
            ],
            "DISC_ANC_OR": [
                _mapping(
                    ancestry_discordance[item], f"discordance.{item}"
                ).get("odds_ratio")
                for item in ancestries
            ],
            "DISC_CMH_P": cmh.get("p_value"),
            "DISC_CMH_CHISQ": cmh.get("chisq"),
            "DISC_UNION_P": union.get("p_value"),
            "DISC_UNION_TEST": union.get("stat_test_name"),
            "DISC_UNION_ANC": (
                None if union.get("gen_ancs") is None else sorted(union["gen_ancs"])
            ),
            "DISC": None if is_discrepant is None else int(is_discrepant),
        }
    )
    expected_fields = set(gnomad_vcf_metadata(projection)["info"])
    if set(info) != expected_fields:
        raise GnomadValidationError("projected source INFO fields are incomplete")
    return info


def _compare_value(expected: Any, observed: Any, context: str) -> None:
    if expected is None or observed is None:
        if expected is not observed:
            raise GnomadValidationError(f"source/BCF missing state differs for {context}")
        return
    if isinstance(expected, Sequence) and not isinstance(expected, (str, bytes)):
        if not isinstance(observed, list) or len(expected) != len(observed):
            raise GnomadValidationError(f"source/BCF array differs for {context}")
        for index, (expected_item, observed_item) in enumerate(
            zip(expected, observed, strict=True)
        ):
            _compare_value(expected_item, observed_item, f"{context}[{index}]")
        return
    if isinstance(expected, float):
        if not isinstance(observed, (int, float)) or isinstance(observed, bool):
            raise GnomadValidationError(f"source/BCF type differs for {context}")
        if math.isnan(expected) or math.isnan(float(observed)):
            raise GnomadValidationError(f"NaN is not valid concordance evidence for {context}")
        if math.isinf(expected) or math.isinf(float(observed)):
            if expected != float(observed):
                raise GnomadValidationError(f"source/BCF value differs for {context}")
            return
        try:
            bcf_expected = struct.unpack("!f", struct.pack("!f", expected))[0]
        except OverflowError:
            bcf_expected = math.copysign(math.inf, expected)
        if not math.isclose(
            bcf_expected, float(observed), rel_tol=5e-6, abs_tol=1e-45
        ):
            raise GnomadValidationError(f"source/BCF value differs for {context}")
        return
    if expected != observed:
        raise GnomadValidationError(f"source/BCF value differs for {context}")


def _frequency_validations(
    source_row: Mapping[str, Any],
    layout: GnomadProjectionLayout,
    info: Mapping[str, Any],
    projection: GnomadProjectionContract,
) -> tuple[GnomadFrequencyValidation, ...]:
    results = []
    for dataset_contract in projection.datasets:
        indices = layout.dataset(dataset_contract.name)
        source_value = source_row.get(dataset_contract.row_field)
        source = (
            None
            if source_value is None
            else _mapping(source_value, dataset_contract.name)
        )
        prefix = DATASET_VCF_PREFIXES[dataset_contract.name]
        checks: list[tuple[str | None, int, int | None, str, str]] = [
            (None, indices.adjusted_frequency, None, f"{prefix}_AC", f"{prefix}_AN")
        ]
        checks.extend(
            (
                ancestry.genetic_ancestry,
                ancestry.adjusted_frequency,
                position,
                f"{prefix}_AC_ANC",
                f"{prefix}_AN_ANC",
            )
            for position, ancestry in enumerate(indices.ancestries)
        )
        for ancestry, source_index, array_index, ac_field, an_field in checks:
            source_frequency = None
            if source is not None:
                frequencies = source.get("freq")
                if isinstance(frequencies, Sequence) and source_index < len(frequencies):
                    source_frequency = frequencies[source_index]
            source_af = None
            if source_frequency is not None:
                source_af = _mapping(source_frequency, "source frequency").get("AF")
            ac_value = info[ac_field]
            an_value = info[an_field]
            ac = ac_value if array_index is None else ac_value[array_index]
            an = an_value if array_index is None else an_value[array_index]
            comparable = ac is not None and an is not None and an > 0
            derived = ac / an if comparable else None
            if comparable:
                if source_af is None or not math.isclose(
                    derived, float(source_af), rel_tol=1e-12, abs_tol=1e-15
                ):
                    raise GnomadValidationError(
                        f"AC/AN-derived AF differs from source for {dataset_contract.name}"
                        f"/{ancestry or 'adjusted'}"
                    )
            elif ac is None or an is None:
                if source_af is not None:
                    raise GnomadValidationError(
                        f"missing AC/AN has a source AF for {dataset_contract.name}"
                        f"/{ancestry or 'adjusted'}"
                    )
            results.append(
                GnomadFrequencyValidation(
                    dataset=dataset_contract.name,
                    ancestry=ancestry,
                    ac=ac,
                    an=an,
                    derived_af=derived,
                    source_af=source_af,
                    comparable=comparable,
                )
            )
    return tuple(results)


def _expectation_checks(
    sentinel: GnomadSentinel,
    info: Mapping[str, Any],
    frequencies: tuple[GnomadFrequencyValidation, ...],
) -> Mapping[str, bool]:
    ancestry_counts = [
        item.ac
        for item in frequencies
        if item.ancestry is not None and item.ac is not None
    ]
    missing_exome_scalar_fields = (
        "EX_AC",
        "EX_AN",
        "EX_NHOMALT",
        "EX_NHEMI",
        "EX_FAF95",
        "EX_FAF99",
        "EX_GRPMAX_AC",
        "EX_GRPMAX_AF",
        "EX_GRPMAX_AN",
        "EX_GRPMAX_NHOMALT",
        "EX_GRPMAX_ANC",
        "EX_FAFMAX95",
        "EX_FAFMAX95_ANC",
        "EX_FAFMAX99",
        "EX_FAFMAX99_ANC",
    )
    missing_exome_array_fields = (
        "EX_AC_ANC",
        "EX_AN_ANC",
        "EX_NHOMALT_ANC",
        "EX_NHEMI_ANC",
        "EX_FAF95_ANC",
        "EX_FAF99_ANC",
    )
    checks = {
        "ac-an-derived-frequency": any(item.comparable for item in frequencies),
        "ancestry-values": len(ancestry_counts) >= 2 and len(set(ancestry_counts)) >= 2,
        "genome-mid-faf-missing": (
            info["GN_FAF95_ANC"][3] is None and info["GN_FAF99_ANC"][3] is None
        ),
        "exome-data-missing": (
            all(info[field] is None for field in missing_exome_scalar_fields)
            and all(
                info[field] == [None] * len(info[field])
                for field in missing_exome_array_fields
            )
        ),
        "exome-filters-missing": (
            info["EX_FILTERS_DEFINED"] == 0 and info["EX_FILTERS"] is None
        ),
        "exome-filters-empty": (
            info["EX_FILTERS_DEFINED"] == 1 and info["EX_FILTERS"] == []
        ),
        "exome-filters-nonempty": (
            info["EX_FILTERS_DEFINED"] == 1
            and isinstance(info["EX_FILTERS"], list)
            and len(info["EX_FILTERS"]) > 0
        ),
        "genome-filters-empty": (
            info["GN_FILTERS_DEFINED"] == 1 and info["GN_FILTERS"] == []
        ),
        "genome-filters-nonempty": (
            info["GN_FILTERS_DEFINED"] == 1
            and isinstance(info["GN_FILTERS"], list)
            and len(info["GN_FILTERS"]) > 0
        ),
        "discordance-missing": info["DISC"] is None,
        "discordance-false": info["DISC"] == 0,
        "discordance-true": info["DISC"] == 1,
    }
    selected = {name: checks[name] for name in sentinel.expectations}
    failed = [name for name, passed in selected.items() if not passed]
    if failed:
        raise GnomadValidationError(
            f"sentinel {sentinel.identifier} did not satisfy: {', '.join(failed)}"
        )
    return selected


def _plain_source_value(value: Any) -> Any:
    if isinstance(value, Mapping):
        return {key: _plain_source_value(item) for key, item in value.items()}
    if isinstance(value, (set, frozenset)):
        return {_plain_source_value(item) for item in value}
    if isinstance(value, (str, bytes, int, float, bool)) or value is None:
        return value
    items = getattr(value, "items", None)
    if callable(items):
        return {key: _plain_source_value(item) for key, item in items()}
    try:
        return [_plain_source_value(item) for item in value]
    except TypeError:
        return value


def _collect_source_rows(
    table: Any,
    hail_module: Any,
    projection: GnomadProjectionContract,
    sentinels: tuple[GnomadSentinel, ...],
) -> Mapping[tuple[str, int, str, str], Mapping[str, Any]]:
    intervals = [
        hail_module.parse_locus_interval(
            f"{item.chromosome}:{item.position}-{item.position + 1}",
            reference_genome=projection.source.assembly,
        )
        for item in sentinels
    ]
    rows = hail_module.filter_intervals(table, intervals).collect()
    result: dict[tuple[str, int, str, str], Mapping[str, Any]] = {}
    for row_value in rows:
        row = _mapping(row_value, "source row")
        locus = row.get("locus")
        alleles = row.get("alleles")
        if (
            locus is None
            or not isinstance(alleles, Sequence)
            or isinstance(alleles, (str, bytes))
            or len(alleles) != 2
        ):
            continue
        key = (str(locus.contig), int(locus.position), str(alleles[0]), str(alleles[1]))
        if key in result:
            raise GnomadValidationError(f"official source duplicated sentinel allele: {key}")
        result[key] = _plain_source_value(row)
    required = {
        (item.chromosome, item.position, item.reference, item.alternate)
        for item in sentinels
    }
    if set(result).intersection(required) != required:
        missing = sorted(required - set(result))
        raise GnomadValidationError(f"official source lacks sentinel alleles: {missing}")
    return {key: result[key] for key in required}


def _write_validation_report(report: GnomadPilotValidationReport) -> None:
    path = Path(report.report_path)
    completion = path.parent / VALIDATION_COMPLETION_FILE
    temporary = path.with_suffix(".json.partial")
    completion_temporary = completion.with_suffix(".json.partial")
    try:
        with temporary.open("x", encoding="utf-8") as handle:
            handle.write(validation_report_as_json(report))
            handle.flush()
            os.fsync(handle.fileno())
        os.replace(temporary, path)
        _fsync_directory(path.parent)
        document = {
            "report_path": path.name,
            "report_sha256": _sha256_file(path),
            "schema": VALIDATION_COMPLETION_SCHEMA,
        }
        with completion_temporary.open("x", encoding="utf-8") as handle:
            json.dump(document, handle, indent=2, sort_keys=True)
            handle.write("\n")
            handle.flush()
            os.fsync(handle.fileno())
        os.replace(completion_temporary, completion)
        _fsync_directory(path.parent)
    except OSError as error:
        raise GnomadValidationError(f"cannot publish validation report {path}: {error}") from error


def _prepare_validation_publication(path: Path, completion: Path) -> None:
    if completion.exists() or completion.is_symlink():
        raise GnomadValidationError("completed gnomAD pilot validation already exists")
    if path.exists() or path.is_symlink():
        if path.is_symlink() or not path.is_file():
            raise GnomadValidationError("incomplete validation report path is invalid")
        try:
            path.unlink()
            _fsync_directory(path.parent)
        except OSError as error:
            raise GnomadValidationError(
                f"cannot remove incomplete validation report {path}: {error}"
            ) from error
    for partial in (
        path.with_suffix(".json.partial"),
        completion.with_suffix(".json.partial"),
    ):
        try:
            partial.unlink(missing_ok=True)
        except OSError as error:
            raise GnomadValidationError(
                f"cannot remove incomplete validation artifact {partial}: {error}"
            ) from error


def _validate_image_digest(image_digest: str) -> None:
    if (
        not isinstance(image_digest, str)
        or len(image_digest) != len(SHA256_DIGEST) + 64
        or not image_digest.startswith(SHA256_DIGEST)
        or any(character not in "0123456789abcdef" for character in image_digest[7:])
    ):
        raise GnomadValidationError("image_digest must be a SHA-256 digest")


def validate_gnomad_pilot(
    *,
    validation_contract_path: str | Path,
    projection_contract_path: str | Path,
    reference_root: str | Path,
    scratch_root: str | Path,
    workers: int,
    image_digest: str,
) -> GnomadPilotValidationReport:
    contract = load_gnomad_validation_contract(validation_contract_path)
    projection = load_gnomad_projection_contract(projection_contract_path)
    _validate_image_digest(image_digest)
    if projection.source.hail_path != PUBLIC_GNOMAD_SOURCE:
        raise GnomadValidationError(
            "gnomAD Hail source must match the pinned public table"
        )
    if not isinstance(workers, int) or isinstance(workers, bool) or not 1 <= workers <= 16:
        raise GnomadValidationError("workers must be an integer between 1 and 16")
    if tuple(projection.major_genetic_ancestries) != contract.ancestry_order:
        raise GnomadValidationError("validation ancestry order differs from projection")
    projection_sha256 = _sha256_file(projection_contract_path)
    validation_sha256 = _sha256_file(validation_contract_path)
    pilot_root = Path(scratch_root) / "gnomad-hail" / contract.pilot_run_id
    report_path = pilot_root / VALIDATION_REPORT_FILE
    completion_path = pilot_root / VALIDATION_COMPLETION_FILE
    validation_run_id = f"{contract.pilot_run_id}-validation"
    validation_root = Path(scratch_root) / "gnomad-hail" / validation_run_id

    _require_bcftools()
    with _pilot_run_lock(reference_root, scratch_root, contract.pilot_run_id):
        _prepare_validation_publication(report_path, completion_path)
        pilot_report, measurements, pilot_report_path = _validate_pilot_report(
            pilot_root, contract, projection, projection_sha256
        )
        shards = _validate_shards(pilot_root, measurements, projection)
        _verify_source_root_metadata(projection)
        inventory_before = _source_inventory(projection)
        if (
            inventory_before.sha256 != contract.source_object_inventory_sha256
            or inventory_before.object_count != pilot_report.get("source_object_count")
            or inventory_before.total_bytes != pilot_report.get("source_object_bytes")
        ):
            raise GnomadValidationError("live gnomAD source inventory differs from pilot")

        if validation_root.exists() or validation_root.is_symlink():
            _remove_failed_run(validation_root, scratch_root)
        layout = prepare_hail_scratch(
            reference_root, scratch_root, validation_run_id
        )
        hail_module: Any | None = None
        try:
            hail_module = initialize_local_hail(layout, workers)
            table = hail_module.read_table(projection.source.hail_path)
            if table.n_partitions() != projection.source.expected_partitions:
                raise GnomadValidationError("official source partition count differs")
            source_globals = hail_module.eval(table.index_globals())
            projection_layout = resolve_gnomad_projection_layout(
                projection, source_globals
            )
            source_rows = _collect_source_rows(
                table, hail_module, projection, contract.sentinels
            )
            _validate_runtime_confinement(hail_module, layout)
        finally:
            if hail_module is not None:
                hail_module.stop()
            _remove_failed_run(layout.run_root, scratch_root)

        inventory_after = _source_inventory(projection)
        if inventory_after != inventory_before:
            raise GnomadValidationError(
                "gnomAD source object inventory changed during validation"
            )

        shard_paths = {item.chromosome: Path(item.bcf_path) for item in shards}
        sentinel_results = []
        for sentinel in contract.sentinels:
            filter_value, bcf_info, matches = _query_sentinel(
                shard_paths[sentinel.chromosome], sentinel, projection
            )
            key = (
                sentinel.chromosome,
                sentinel.position,
                sentinel.reference,
                sentinel.alternate,
            )
            source_row = source_rows[key]
            projected = project_gnomad_joint_row(
                source_row, projection_layout, is_autosome_or_par=True
            )
            source_info = _projected_info(projected, projection)
            for identifier in gnomad_vcf_metadata(projection)["info"]:
                _compare_value(
                    source_info[identifier],
                    bcf_info[identifier],
                    f"{sentinel.identifier}.{identifier}",
                )
            frequencies = _frequency_validations(
                source_row, projection_layout, bcf_info, projection
            )
            expectation_checks = _expectation_checks(
                sentinel, bcf_info, frequencies
            )
            sentinel_results.append(
                GnomadSentinelValidation(
                    identifier=sentinel.identifier,
                    chromosome=sentinel.chromosome,
                    position=sentinel.position,
                    reference=sentinel.reference,
                    alternate=sentinel.alternate,
                    expectations=sentinel.expectations,
                    expectation_checks=expectation_checks,
                    bcf_filter=filter_value,
                    csi_query_matches=matches,
                    info_fields_checked=len(bcf_info),
                    semantic_observations=_json_safe(
                        {
                            identifier: bcf_info[identifier]
                            for identifier in SEMANTIC_INFO_FIELDS
                        }
                    ),
                    frequencies=frequencies,
                )
            )

        report = GnomadPilotValidationReport(
            schema=VALIDATION_SCHEMA,
            mode="official-source-chr21-chr22-sentinel-validation",
            input_classification="public",
            accepted=True,
            mismatch_count=0,
            validation_id=contract.validation_id,
            validation_contract_sha256=validation_sha256,
            projection_id=projection.projection_id,
            projection_contract_sha256=projection_sha256,
            source=projection.source.hail_path,
            source_root_metadata_object_generation=(
                projection.source.root_metadata_object_generation
            ),
            source_root_metadata_md5_base64=projection.source.root_metadata_md5_base64,
            source_object_inventory_sha256=inventory_after.sha256,
            source_object_count=inventory_after.object_count,
            source_object_bytes=inventory_after.total_bytes,
            pilot_report_path=str(pilot_report_path),
            pilot_report_sha256=contract.pilot_report_sha256,
            pilot_image_digest=str(pilot_report["image_digest"]),
            validator_image_digest=image_digest,
            ancestry_order=contract.ancestry_order,
            shards=shards,
            sentinels=tuple(sentinel_results),
            hail_version=HAIL_VERSION,
            pyspark_version=PYSPARK_VERSION,
            bcftools_version=BCFTOOLS_VERSION,
            workers=workers,
            spark_master=f"local[{workers}]",
            scratch_root=str(Path(scratch_root)),
            report_path=str(report_path),
        )
        _write_validation_report(report)
        return report
