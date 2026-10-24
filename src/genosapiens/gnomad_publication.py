"""Serial, fail-closed publication of the reduced public gnomAD BCF shards."""

from __future__ import annotations

from collections.abc import Iterator, Mapping, Sequence
from dataclasses import asdict, dataclass, replace
from datetime import date
import hashlib
import json
import os
from pathlib import Path, PurePosixPath
import re
import shutil
import subprocess
from typing import Any
import uuid

from genosapiens.gnomad_hail import (
    GnomadHailError,
    PUBLIC_GNOMAD_SOURCE,
    _validate_runtime_confinement,
    initialize_local_hail,
    prepare_hail_scratch,
)
from genosapiens.gnomad_pilot import (
    BCFTOOLS_VERSION,
    OUTPUT_CAP_BYTES,
    PILOT_CHROMOSOMES,
    SHA256_DIGEST,
    GnomadSourceInventory,
    _ResourceSampler,
    _fsync_directory,
    _fsync_file,
    _pilot_run_lock,
    _remove_failed_run,
    _require_bcftools,
    _run_bcftools,
    _sha256_file,
    _source_inventory,
    _verify_source_root_metadata,
    materialize_gnomad_shard,
)
from genosapiens.gnomad_projection import (
    GnomadProjectionContract,
    load_gnomad_projection_contract,
    resolve_gnomad_projection_layout,
)
from genosapiens.gnomad_validation import (
    VALIDATION_COMPLETION_FILE,
    VALIDATION_COMPLETION_SCHEMA,
    VALIDATION_REPORT_FILE,
    VALIDATION_SCHEMA,
    GnomadPilotValidationReport,
    GnomadValidationContract,
    _validate_header,
    _validate_pilot_report,
    load_gnomad_validation_contract,
)
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
from genosapiens.storage import CapacityPlan, DEFAULT_RESERVE_BYTES, require_capacity


PUBLICATION_SCHEMA = "genosapiens.gnomad-joint-publication"
PUBLICATION_TRANSFORMATION_SCHEMA = "genosapiens.gnomad-joint-bcf-publication"
PUBLICATION_RESULT_SCHEMA = "genosapiens.gnomad-joint-publication-result"
CANONICAL_CONTIGS = tuple(f"chr{number}" for number in range(1, 23)) + (
    "chrX",
    "chrY",
)
RECORD_LINE = re.compile(r"^([^\t]+)\t([1-9][0-9]*)\t([^\t]+)\t([^\t]+)$")
IDENTIFIER = re.compile(r"^[a-z0-9][a-z0-9._-]*$")
SHA256 = re.compile(r"^[0-9a-f]{64}$")


class GnomadPublicationError(GnomadHailError):
    """Raised when the full gnomAD publication cannot be safely accepted."""


def _mapping(value: Any, context: str) -> Mapping[str, Any]:
    if not isinstance(value, Mapping):
        raise GnomadPublicationError(f"{context} must be an object")
    if not all(isinstance(key, str) for key in value):
        raise GnomadPublicationError(f"{context} keys must be strings")
    return value


def _string(data: Mapping[str, Any], name: str, context: str) -> str:
    value = data.get(name)
    if not isinstance(value, str) or not value:
        raise GnomadPublicationError(f"{context}.{name} must be a non-empty string")
    return value


def _safe_relative(value: str, context: str) -> str:
    path = PurePosixPath(value)
    if path.is_absolute() or path == PurePosixPath(".") or ".." in path.parts:
        raise GnomadPublicationError(f"{context} must be a safe relative path")
    return str(path)


@dataclass(frozen=True, slots=True)
class GnomadPublicationContract:
    resource_id: str
    release: str
    assembly: str
    projection_id: str
    pilot_validation_id: str
    canonical_contigs: tuple[str, ...]
    output_directory: str
    manifest_path: str
    aggregate_output_cap_bytes: int
    schema: str = PUBLICATION_SCHEMA

    @classmethod
    def from_mapping(cls, data: Mapping[str, Any]) -> "GnomadPublicationContract":
        expected = {
            "schema",
            "resource_id",
            "release",
            "assembly",
            "projection_id",
            "pilot_validation_id",
            "canonical_contigs",
            "output_directory",
            "manifest_path",
            "aggregate_output_cap_bytes",
        }
        contigs = data.get("canonical_contigs")
        cap = data.get("aggregate_output_cap_bytes")
        if (
            set(data) != expected
            or data.get("schema") != PUBLICATION_SCHEMA
            or not isinstance(contigs, list)
            or not isinstance(cap, int)
            or isinstance(cap, bool)
        ):
            raise GnomadPublicationError("gnomAD publication contract fields are invalid")
        contract = cls(
            resource_id=_string(data, "resource_id", "publication contract"),
            release=_string(data, "release", "publication contract"),
            assembly=_string(data, "assembly", "publication contract"),
            projection_id=_string(data, "projection_id", "publication contract"),
            pilot_validation_id=_string(
                data, "pilot_validation_id", "publication contract"
            ),
            canonical_contigs=tuple(contigs),
            output_directory=_safe_relative(
                _string(data, "output_directory", "publication contract"),
                "publication output_directory",
            ),
            manifest_path=_safe_relative(
                _string(data, "manifest_path", "publication contract"),
                "publication manifest_path",
            ),
            aggregate_output_cap_bytes=cap,
        )
        if (
            not IDENTIFIER.fullmatch(contract.resource_id)
            or not IDENTIFIER.fullmatch(contract.projection_id)
            or not IDENTIFIER.fullmatch(contract.pilot_validation_id)
            or contract.canonical_contigs != CANONICAL_CONTIGS
            or contract.aggregate_output_cap_bytes != OUTPUT_CAP_BYTES
            or not contract.output_directory.startswith("derived/gnomad/")
            or not contract.manifest_path.startswith("manifests/gnomad/")
        ):
            raise GnomadPublicationError("gnomAD publication scope is invalid")
        return contract


def load_gnomad_publication_contract(path: str | Path) -> GnomadPublicationContract:
    try:
        data = json.loads(Path(path).read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as error:
        raise GnomadPublicationError(
            f"cannot load gnomAD publication contract: {error}"
        ) from error
    return GnomadPublicationContract.from_mapping(_mapping(data, "publication contract"))


@dataclass(frozen=True, slots=True)
class GnomadPublishedShard:
    chromosome: str
    rows: int
    bcf_path: str
    bcf_bytes: int
    bcf_sha256: str
    csi_path: str
    csi_bytes: int
    csi_sha256: str
    header_info_fields: int
    single_contig_verified: bool
    coordinate_allele_order_verified: bool
    index_row_count_verified: bool
    csi_queryable: bool

    @property
    def output_bytes(self) -> int:
        return self.bcf_bytes + self.csi_bytes


@dataclass(frozen=True, slots=True)
class GnomadFullPublicationResult:
    resource_id: str
    release: str
    output_directory: str
    aggregate_output_bytes: int
    aggregate_output_cap_bytes: int
    reserve_bytes: int
    manifest_path: str
    manifest_sha256: str
    shards: tuple[GnomadPublishedShard, ...]
    reused: bool
    schema: str = PUBLICATION_RESULT_SCHEMA


@dataclass(frozen=True, slots=True)
class _AcceptedPilot:
    pilot_root: Path
    pilot_report: Mapping[str, Any]
    measurements: tuple[Mapping[str, Any], ...]
    validation_report_sha256: str
    source_inventory: GnomadSourceInventory


def _file_is_regular(path: Path, context: str) -> None:
    if path.is_symlink() or not path.is_file():
        raise GnomadPublicationError(f"{context} is not a regular file: {path}")


def _expected_shard_name(chromosome: str) -> str:
    return f"gnomad.joint.v4.1.{chromosome}.bcf"


def _stream_bcf_records(path: Path) -> Iterator[tuple[str, int, str, str]]:
    command = [
        "bcftools",
        "query",
        "--format",
        "%CHROM\\t%POS\\t%REF\\t%ALT\\n",
        str(path),
    ]
    try:
        process = subprocess.Popen(
            command,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            bufsize=1,
        )
    except OSError as error:
        raise GnomadPublicationError(f"cannot execute bcftools: {error}") from error
    assert process.stdout is not None
    assert process.stderr is not None
    try:
        for raw_line in process.stdout:
            match = RECORD_LINE.fullmatch(raw_line.rstrip("\n"))
            if match is None:
                raise GnomadPublicationError("bcftools record stream is malformed")
            chromosome, position, reference, alternate = match.groups()
            yield chromosome, int(position), reference, alternate
        stderr = process.stderr.read().strip()
        if process.wait() != 0:
            raise GnomadPublicationError(
                f"bcftools record stream failed: {stderr or 'unknown error'}"
            )
    finally:
        process.stdout.close()
        process.stderr.close()
        if process.poll() is None:
            process.kill()
            process.wait()


def _parse_index_count(value: str, chromosome: str) -> int:
    try:
        rows = int(value)
    except ValueError as error:
        raise GnomadPublicationError(
            f"bcftools returned an invalid row count for {chromosome}: {value!r}"
        ) from error
    if rows <= 0:
        raise GnomadPublicationError(f"{chromosome} shard has no records")
    return rows


def _validate_index_membership(path: Path, chromosome: str, rows: int) -> None:
    lines = [
        line.split("\t")
        for line in _run_bcftools(("index", "--stats", str(path))).splitlines()
    ]
    if len(lines) != 1 or len(lines[0]) != 3:
        raise GnomadPublicationError(
            f"{chromosome} CSI must contain exactly one contig"
        )
    observed_chromosome, _length, observed_rows = lines[0]
    if observed_chromosome != chromosome:
        raise GnomadPublicationError(
            f"{chromosome} CSI contains a different contig: {observed_chromosome}"
        )
    try:
        if int(observed_rows) != rows:
            raise GnomadPublicationError(f"{chromosome} CSI statistics row count differs")
    except ValueError as error:
        raise GnomadPublicationError(
            f"{chromosome} CSI statistics row count is invalid"
        ) from error


def _validate_queryability(
    path: Path, record: tuple[str, int, str, str]
) -> None:
    chromosome, position, reference, alternate = record
    output = _run_bcftools(
        (
            "query",
            "--regions",
            f"{chromosome}:{position}-{position}",
            "--format",
            "%CHROM\\t%POS\\t%REF\\t%ALT\\n",
            str(path),
        )
    )
    matches = []
    for line in output.splitlines():
        match = RECORD_LINE.fullmatch(line)
        if match is None:
            raise GnomadPublicationError("CSI query returned malformed records")
        values = match.groups()
        if values == (chromosome, str(position), reference, alternate):
            matches.append(values)
    if len(matches) != 1:
        raise GnomadPublicationError(
            f"{chromosome} CSI query did not recover exactly one first record"
        )


def _validate_shard(
    path: Path,
    chromosome: str,
    projection: GnomadProjectionContract,
) -> GnomadPublishedShard:
    csi_path = Path(f"{path}.csi")
    _file_is_regular(path, f"{chromosome} BCF")
    _file_is_regular(csi_path, f"{chromosome} CSI")
    rows = _parse_index_count(
        _run_bcftools(("index", "--nrecords", str(path))), chromosome
    )
    _validate_index_membership(path, chromosome, rows)
    header_fields = _validate_header(
        _run_bcftools(("view", "--header-only", str(path))), projection
    )

    previous: tuple[int, str, str] | None = None
    first: tuple[str, int, str, str] | None = None
    streamed_rows = 0
    for observed_chromosome, position, reference, alternate in _stream_bcf_records(path):
        if observed_chromosome != chromosome:
            raise GnomadPublicationError(
                f"{chromosome} shard contains {observed_chromosome} records"
            )
        key = (position, reference, alternate)
        if previous is not None and key <= previous:
            raise GnomadPublicationError(
                f"{chromosome} shard is not in coordinate and allele order"
            )
        previous = key
        first = first or (observed_chromosome, position, reference, alternate)
        streamed_rows += 1
    if streamed_rows != rows or first is None:
        raise GnomadPublicationError(f"{chromosome} BCF record count differs from CSI")
    _validate_queryability(path, first)

    return GnomadPublishedShard(
        chromosome=chromosome,
        rows=rows,
        bcf_path=str(path),
        bcf_bytes=path.stat().st_size,
        bcf_sha256=_sha256_file(path),
        csi_path=str(csi_path),
        csi_bytes=csi_path.stat().st_size,
        csi_sha256=_sha256_file(csi_path),
        header_info_fields=header_fields,
        single_contig_verified=True,
        coordinate_allele_order_verified=True,
        index_row_count_verified=True,
        csi_queryable=True,
    )


def _validate_bundle(
    release: Path,
    contract: GnomadPublicationContract,
    projection: GnomadProjectionContract,
) -> tuple[GnomadPublishedShard, ...]:
    if release.is_symlink() or not release.is_dir():
        raise GnomadPublicationError("gnomAD publication directory is unavailable")
    expected_names = {
        name
        for chromosome in contract.canonical_contigs
        for name in (_expected_shard_name(chromosome), f"{_expected_shard_name(chromosome)}.csi")
    }
    try:
        observed_names = {path.name for path in release.iterdir()}
    except OSError as error:
        raise GnomadPublicationError(
            f"cannot inspect gnomAD publication directory: {error}"
        ) from error
    if observed_names != expected_names:
        raise GnomadPublicationError("gnomAD publication contains unexpected shards")
    return tuple(
        _validate_shard(release / _expected_shard_name(chromosome), chromosome, projection)
        for chromosome in contract.canonical_contigs
    )


def _accepted_pilot(
    *,
    validation_contract: GnomadValidationContract,
    validation_contract_path: Path,
    projection: GnomadProjectionContract,
    projection_contract_path: Path,
    scratch_root: Path,
) -> _AcceptedPilot:
    projection_sha256 = _sha256_file(projection_contract_path)
    pilot_root = scratch_root / "gnomad-hail" / validation_contract.pilot_run_id
    pilot_report, measurements, pilot_report_path = _validate_pilot_report(
        pilot_root, validation_contract, projection, projection_sha256
    )
    report_path = pilot_root / VALIDATION_REPORT_FILE
    completion_path = pilot_root / VALIDATION_COMPLETION_FILE
    _file_is_regular(report_path, "gnomAD pilot validation report")
    _file_is_regular(completion_path, "gnomAD pilot validation completion marker")
    try:
        completion = _mapping(
            json.loads(completion_path.read_text(encoding="utf-8")),
            "gnomAD pilot validation completion marker",
        )
        report = _mapping(
            json.loads(report_path.read_text(encoding="utf-8")),
            "gnomAD pilot validation report",
        )
    except (OSError, json.JSONDecodeError) as error:
        raise GnomadPublicationError(
            f"cannot load accepted gnomAD pilot validation: {error}"
        ) from error
    if set(completion) != {"schema", "report_path", "report_sha256"} or (
        completion.get("schema") != VALIDATION_COMPLETION_SCHEMA
        or completion.get("report_path") != report_path.name
        or completion.get("report_sha256") != _sha256_file(report_path)
    ):
        raise GnomadPublicationError("gnomAD pilot validation completion differs")
    if set(report) != set(GnomadPilotValidationReport.__dataclass_fields__) or (
        report.get("schema") != VALIDATION_SCHEMA
        or report.get("accepted") is not True
        or report.get("mismatch_count") != 0
        or report.get("validation_id") != validation_contract.validation_id
        or report.get("validation_contract_sha256")
        != _sha256_file(validation_contract_path)
        or report.get("projection_id") != projection.projection_id
        or report.get("projection_contract_sha256") != projection_sha256
        or report.get("source") != PUBLIC_GNOMAD_SOURCE
        or report.get("pilot_report_path") != str(pilot_report_path)
        or report.get("pilot_report_sha256") != validation_contract.pilot_report_sha256
        or tuple(report.get("ancestry_order", ()))
        != validation_contract.ancestry_order
        or report.get("bcftools_version") != BCFTOOLS_VERSION
    ):
        raise GnomadPublicationError("gnomAD pilot validation is not accepted")
    shards = report.get("shards")
    sentinels = report.get("sentinels")
    if not isinstance(shards, list) or not isinstance(sentinels, list):
        raise GnomadPublicationError("gnomAD pilot validation evidence is invalid")
    if tuple(
        item.get("chromosome") for item in shards if isinstance(item, Mapping)
    ) != PILOT_CHROMOSOMES:
        raise GnomadPublicationError("gnomAD pilot validation shard order differs")
    if {item.get("identifier") for item in sentinels if isinstance(item, Mapping)} != {
        item.identifier for item in validation_contract.sentinels
    }:
        raise GnomadPublicationError("gnomAD pilot validation sentinels differ")
    for measurement, shard in zip(measurements, shards, strict=True):
        shard_data = _mapping(shard, "gnomAD pilot validation shard")
        if (
            shard_data.get("rows") != measurement.get("rows")
            or shard_data.get("bcf_path") != measurement.get("bcf_path")
            or shard_data.get("bcf_bytes") != measurement.get("bcf_bytes")
            or shard_data.get("bcf_sha256") != measurement.get("bcf_sha256")
            or shard_data.get("csi_path") != measurement.get("csi_path")
            or shard_data.get("csi_bytes") != measurement.get("csi_bytes")
            or shard_data.get("csi_sha256") != measurement.get("csi_sha256")
            or shard_data.get("header_verified") is not True
            or shard_data.get("index_row_count_verified") is not True
        ):
            raise GnomadPublicationError("gnomAD pilot validation shard identity differs")
    for sentinel in sentinels:
        sentinel_data = _mapping(sentinel, "gnomAD pilot validation sentinel")
        checks = sentinel_data.get("expectation_checks")
        if (
            sentinel_data.get("csi_query_matches") != 1
            or not isinstance(checks, Mapping)
            or not checks
            or not all(value is True for value in checks.values())
        ):
            raise GnomadPublicationError("gnomAD pilot validation sentinel failed")
    return _AcceptedPilot(
        pilot_root=pilot_root,
        pilot_report=pilot_report,
        measurements=measurements,
        validation_report_sha256=_sha256_file(report_path),
        source_inventory=GnomadSourceInventory(
            object_count=int(report["source_object_count"]),
            total_bytes=int(report["source_object_bytes"]),
            sha256=str(report["source_object_inventory_sha256"]),
        ),
    )


def _transformation_parameters_sha256(
    publication_contract_path: Path,
    projection_contract_path: Path,
    validation_contract_path: Path,
    validation_contract: GnomadValidationContract,
) -> str:
    values = {
        "pilot_report_sha256": validation_contract.pilot_report_sha256,
        "projection_contract_sha256": _sha256_file(projection_contract_path),
        "publication_contract_sha256": _sha256_file(publication_contract_path),
        "source_object_inventory_sha256": validation_contract.source_object_inventory_sha256,
        "validation_contract_sha256": _sha256_file(validation_contract_path),
    }
    return hashlib.sha256(
        json.dumps(values, sort_keys=True, separators=(",", ":")).encode("ascii")
    ).hexdigest()


def _manifest(
    *,
    contract: GnomadPublicationContract,
    projection: GnomadProjectionContract,
    retrieved_on: str,
    image_digest: str,
    parameters_sha256: str,
    shards: Sequence[GnomadPublishedShard],
    reference_root: Path,
) -> PublicResourceManifest:
    outputs = []
    for shard in shards:
        for path, size, digest in (
            (Path(shard.bcf_path), shard.bcf_bytes, shard.bcf_sha256),
            (Path(shard.csi_path), shard.csi_bytes, shard.csi_sha256),
        ):
            outputs.append(
                ResourceOutput(
                    path=str(path.relative_to(reference_root)),
                    size_bytes=size,
                    sha256=digest,
                    retention=RetentionPolicy("durable", None),
                )
            )
    return PublicResourceManifest(
        resource_id=contract.resource_id,
        release=contract.release,
        source=SourceIdentity(
            projection.source.root_metadata_url,
            object_generation=projection.source.root_metadata_object_generation,
        ),
        license=projection.source.license,
        assembly=contract.assembly,
        retrieved_on=retrieved_on,
        transformation=Transformation(PUBLICATION_TRANSFORMATION_SCHEMA, parameters_sha256),
        image_digest=image_digest,
        outputs=tuple(outputs),
    )


def _expected_output_paths(contract: GnomadPublicationContract) -> set[str]:
    return {
        str(PurePosixPath(contract.output_directory) / name)
        for chromosome in contract.canonical_contigs
        for name in (
            _expected_shard_name(chromosome),
            f"{_expected_shard_name(chromosome)}.csi",
        )
    }


def _load_manifest(path: Path) -> PublicResourceManifest:
    _file_is_regular(path, "gnomAD publication manifest")
    try:
        return PublicResourceManifest.from_mapping(
            _mapping(json.loads(path.read_text(encoding="utf-8")), "gnomAD manifest")
        )
    except (OSError, json.JSONDecodeError, ManifestError) as error:
        raise GnomadPublicationError(
            f"cannot load gnomAD publication manifest: {error}"
        ) from error


def _verify_manifest_contract(
    manifest: PublicResourceManifest,
    contract: GnomadPublicationContract,
    projection: GnomadProjectionContract,
    parameters_sha256: str,
) -> None:
    if (
        manifest.resource_id != contract.resource_id
        or manifest.release != contract.release
        or manifest.license != projection.source.license
        or manifest.assembly != contract.assembly
        or manifest.source.url != projection.source.root_metadata_url
        or manifest.source.sha256 is not None
        or manifest.source.object_generation
        != projection.source.root_metadata_object_generation
        or manifest.transformation
        != Transformation(PUBLICATION_TRANSFORMATION_SCHEMA, parameters_sha256)
        or {output.path for output in manifest.outputs} != _expected_output_paths(contract)
    ):
        raise GnomadPublicationError("completed gnomAD manifest differs from contract")


def _require_within_root(root: Path, path: Path, context: str) -> None:
    if not path.resolve(strict=False).is_relative_to(root):
        raise GnomadPublicationError(f"{context} escapes reference_root")


def _link_or_copy(source: Path, destination: Path) -> None:
    _file_is_regular(source, "accepted pilot shard")
    try:
        os.link(source, destination, follow_symlinks=False)
    except OSError:
        try:
            with (
                source.open("rb") as input_handle,
                destination.open("xb") as output_handle,
            ):
                shutil.copyfileobj(input_handle, output_handle, length=1024 * 1024)
                output_handle.flush()
                os.fsync(output_handle.fileno())
        except OSError as error:
            raise GnomadPublicationError(
                f"cannot stage accepted pilot shard {source}: {error}"
            ) from error
    _fsync_file(destination)


def _stage_pilot_shards(
    staging: Path,
    accepted: _AcceptedPilot,
) -> None:
    source_root = accepted.pilot_root / "pilot-output"
    for chromosome, measurement in zip(
        PILOT_CHROMOSOMES, accepted.measurements, strict=True
    ):
        bcf_source = source_root / _expected_shard_name(chromosome)
        csi_source = Path(f"{bcf_source}.csi")
        bcf_destination = staging / bcf_source.name
        csi_destination = staging / csi_source.name
        if (
            measurement.get("bcf_path") != str(bcf_source)
            or measurement.get("csi_path") != str(csi_source)
        ):
            raise GnomadPublicationError("accepted pilot shard path differs")
        _link_or_copy(bcf_source, bcf_destination)
        _link_or_copy(csi_source, csi_destination)


def _pilot_output_matches(
    shards: Sequence[GnomadPublishedShard], accepted: _AcceptedPilot
) -> None:
    by_chromosome = {shard.chromosome: shard for shard in shards}
    for chromosome, measurement in zip(
        PILOT_CHROMOSOMES, accepted.measurements, strict=True
    ):
        shard = by_chromosome.get(chromosome)
        if shard is None or (
            shard.rows != measurement.get("rows")
            or shard.bcf_bytes != measurement.get("bcf_bytes")
            or shard.bcf_sha256 != measurement.get("bcf_sha256")
            or shard.csi_bytes != measurement.get("csi_bytes")
            or shard.csi_sha256 != measurement.get("csi_sha256")
        ):
            raise GnomadPublicationError(
                f"published {chromosome} pilot shard differs from accepted validation"
            )


def _projected_values(report: Mapping[str, Any]) -> tuple[int, int]:
    extrapolation = _mapping(report.get("extrapolation"), "pilot extrapolation")
    output_bytes = extrapolation.get("projected_output_bytes")
    scratch_bytes = extrapolation.get("projected_peak_scratch_bytes")
    if (
        not isinstance(output_bytes, int)
        or isinstance(output_bytes, bool)
        or output_bytes <= 0
        or not isinstance(scratch_bytes, int)
        or isinstance(scratch_bytes, bool)
        or scratch_bytes < 0
    ):
        raise GnomadPublicationError("accepted pilot extrapolation is invalid")
    return output_bytes, scratch_bytes


def _make_release_read_only(release: Path, contract: GnomadPublicationContract) -> None:
    for chromosome in contract.canonical_contigs:
        bcf_path = release / _expected_shard_name(chromosome)
        csi_path = Path(f"{bcf_path}.csi")
        os.chmod(bcf_path, 0o444)
        os.chmod(csi_path, 0o444)
    os.chmod(release, 0o755)
    _fsync_directory(release)


def _result(
    *,
    contract: GnomadPublicationContract,
    shards: tuple[GnomadPublishedShard, ...],
    reserve_bytes: int,
    manifest_path: Path,
    manifest_sha256: str,
    reused: bool,
) -> GnomadFullPublicationResult:
    return GnomadFullPublicationResult(
        resource_id=contract.resource_id,
        release=contract.release,
        output_directory=contract.output_directory,
        aggregate_output_bytes=sum(shard.output_bytes for shard in shards),
        aggregate_output_cap_bytes=contract.aggregate_output_cap_bytes,
        reserve_bytes=reserve_bytes,
        manifest_path=str(manifest_path),
        manifest_sha256=manifest_sha256,
        shards=shards,
        reused=reused,
    )


def _validate_image_digest(image_digest: str) -> None:
    if (
        not isinstance(image_digest, str)
        or len(image_digest) != len(SHA256_DIGEST) + 64
        or not image_digest.startswith(SHA256_DIGEST)
        or SHA256.fullmatch(image_digest[7:]) is None
    ):
        raise GnomadPublicationError("image_digest must be a SHA-256 digest")


def publish_gnomad_joint_frequency(
    *,
    publication_contract_path: str | Path,
    projection_contract_path: str | Path,
    validation_contract_path: str | Path,
    reference_root: str | Path,
    scratch_root: str | Path,
    reserve_bytes: int,
    image_digest: str,
    retrieved_on: str | None = None,
) -> GnomadFullPublicationResult:
    """Publish the configured canonical shards only after the accepted pilot gate."""
    publication_path = Path(publication_contract_path)
    projection_path = Path(projection_contract_path)
    validation_path = Path(validation_contract_path)
    root = Path(reference_root).resolve(strict=True)
    scratch = Path(scratch_root)
    contract = load_gnomad_publication_contract(publication_path)
    projection = load_gnomad_projection_contract(projection_path)
    validation_contract = load_gnomad_validation_contract(validation_path)
    _validate_image_digest(image_digest)
    if (
        projection.source.hail_path != PUBLIC_GNOMAD_SOURCE
        or projection.projection_id != contract.projection_id
        or validation_contract.projection_id != contract.projection_id
        or validation_contract.validation_id != contract.pilot_validation_id
    ):
        raise GnomadPublicationError("gnomAD publication contracts differ")
    if not isinstance(reserve_bytes, int) or isinstance(reserve_bytes, bool):
        raise GnomadPublicationError("reserve_bytes must be an integer")
    if reserve_bytes < DEFAULT_RESERVE_BYTES:
        raise GnomadPublicationError("publication must preserve at least 100 GB")
    _require_within_root(root, scratch, "gnomAD scratch_root")
    release = root / contract.output_directory
    manifest_path = root / contract.manifest_path
    _require_within_root(root, release, "gnomAD publication directory")
    _require_within_root(root, manifest_path, "gnomAD publication manifest")
    parameters_sha256 = _transformation_parameters_sha256(
        publication_path, projection_path, validation_path, validation_contract
    )

    def reuse_completed() -> GnomadFullPublicationResult | None:
        if not (manifest_path.exists() or manifest_path.is_symlink()):
            return None
        manifest = _load_manifest(manifest_path)
        _verify_manifest_contract(manifest, contract, projection, parameters_sha256)
        try:
            verify_outputs(manifest, root)
        except ManifestError as error:
            raise GnomadPublicationError(
                f"completed gnomAD publication is invalid: {error}"
            ) from error
        shards = _validate_bundle(release, contract, projection)
        if sum(shard.output_bytes for shard in shards) > contract.aggregate_output_cap_bytes:
            raise GnomadPublicationError("completed gnomAD publication exceeds 150 GB")
        return _result(
            contract=contract,
            shards=shards,
            reserve_bytes=reserve_bytes,
            manifest_path=manifest_path,
            manifest_sha256=_sha256_file(manifest_path),
            reused=True,
        )

    _require_bcftools()
    completed = reuse_completed()
    if completed is not None:
        return completed

    with _pilot_run_lock(reference_root, scratch_root, "gnomad-joint-full-publication"):
        completed = reuse_completed()
        if completed is not None:
            return completed
        accepted = _accepted_pilot(
            validation_contract=validation_contract,
            validation_contract_path=validation_path,
            projection=projection,
            projection_contract_path=projection_path,
            scratch_root=scratch,
        )
        projected_output, projected_scratch = _projected_values(accepted.pilot_report)
        if projected_output > contract.aggregate_output_cap_bytes:
            raise GnomadPublicationError("accepted pilot projects output above 150 GB")
        require_capacity(
            root,
            CapacityPlan(
                operation="publish reduced gnomAD v4.1 joint BCF shards",
                output_bytes=projected_output,
                scratch_bytes=projected_scratch,
                reserve_bytes=reserve_bytes,
            ),
        )
        _verify_source_root_metadata(projection)
        inventory_before = _source_inventory(projection)
        if inventory_before != accepted.source_inventory:
            raise GnomadPublicationError("live gnomAD source differs from accepted pilot")

        if release.exists() or release.is_symlink():
            shards = _validate_bundle(release, contract, projection)
            if sum(shard.output_bytes for shard in shards) > contract.aggregate_output_cap_bytes:
                raise GnomadPublicationError("gnomAD publication exceeds 150 GB")
            publication = publish_manifest(
                _manifest(
                    contract=contract,
                    projection=projection,
                    retrieved_on=retrieved_on or date.today().isoformat(),
                    image_digest=image_digest,
                    parameters_sha256=parameters_sha256,
                    shards=shards,
                    reference_root=root,
                ),
                manifest_path,
                root,
            )
            _remove_failed_run(accepted.pilot_root, scratch)
            return _result(
                contract=contract,
                shards=shards,
                reserve_bytes=reserve_bytes,
                manifest_path=publication.path,
                manifest_sha256=publication.manifest_sha256,
                reused=publication.reused,
            )

        release.parent.mkdir(parents=True, exist_ok=True)
        partial_root = release.parent / ".partial"
        _require_within_root(root, partial_root, "gnomAD publication partial directory")
        if partial_root.is_symlink():
            raise GnomadPublicationError("gnomAD publication partial directory is a symlink")
        partial_root.mkdir(mode=0o700, exist_ok=True)
        staging = partial_root / f"{release.name}.{uuid.uuid4().hex}.partial"
        staging.mkdir(mode=0o700)
        _stage_pilot_shards(staging, accepted)
        staged_pilot = tuple(
            _validate_shard(
                staging / _expected_shard_name(chromosome), chromosome, projection
            )
            for chromosome in PILOT_CHROMOSOMES
        )
        _pilot_output_matches(staged_pilot, accepted)

        run_id = f"gnomad-joint-full-{uuid.uuid4().hex}"
        layout = prepare_hail_scratch(reference_root, scratch_root, run_id)
        export_root = layout.run_root / "export"
        export_root.mkdir()
        sampler = _ResourceSampler(
            (
                layout.hail_tmp,
                layout.spark_local,
                layout.spark_warehouse,
                layout.checkpoints,
                layout.logs,
                layout.home,
                layout.cache,
                layout.ivy,
                layout.java_tmp,
                layout.python_tmp,
                export_root,
            )
        )
        hail_module: Any | None = None
        sampler_started = False
        try:
            sampler.start()
            sampler_started = True
            hail_module = initialize_local_hail(layout, workers=16)
            table = hail_module.read_table(projection.source.hail_path)
            if table.n_partitions() != projection.source.expected_partitions:
                raise GnomadPublicationError("gnomAD source partition count differs")
            projection_layout = resolve_gnomad_projection_layout(
                projection, hail_module.eval(table.index_globals())
            )
            for chromosome in contract.canonical_contigs:
                if chromosome in PILOT_CHROMOSOMES:
                    continue
                require_capacity(
                    root,
                    CapacityPlan(
                        operation=f"materialize reduced gnomAD {chromosome}",
                        output_bytes=max(
                            projected_output
                            - sum(
                                path.stat().st_size
                                for path in staging.iterdir()
                                if path.is_file()
                            ),
                            0,
                        ),
                        scratch_bytes=projected_scratch,
                        reserve_bytes=reserve_bytes,
                    ),
                )
                materialize_gnomad_shard(
                    hail_module=hail_module,
                    table=table,
                    contract=projection,
                    projection_layout=projection_layout,
                    chromosome=chromosome,
                    export_root=export_root,
                    output_root=staging,
                    sampler=sampler,
                )
            _validate_runtime_confinement(hail_module, layout)
        finally:
            if hail_module is not None:
                hail_module.stop()
            if sampler_started:
                sampler.close()

        inventory_after = _source_inventory(projection)
        if inventory_after != inventory_before:
            raise GnomadPublicationError("gnomAD source changed during full publication")
        shards = _validate_bundle(staging, contract, projection)
        _pilot_output_matches(shards, accepted)
        aggregate_output = sum(shard.output_bytes for shard in shards)
        if aggregate_output > contract.aggregate_output_cap_bytes:
            raise GnomadPublicationError("gnomAD publication exceeds 150 GB")
        if shutil.disk_usage(root).free < reserve_bytes:
            raise GnomadPublicationError("gnomAD publication consumed the free-space reserve")
        _make_release_read_only(staging, contract)
        try:
            os.replace(staging, release)
        except OSError as error:
            raise GnomadPublicationError(
                f"cannot atomically publish gnomAD shards: {error}"
            ) from error
        _fsync_directory(release.parent)
        publication = publish_manifest(
            _manifest(
                contract=contract,
                projection=projection,
                retrieved_on=retrieved_on or date.today().isoformat(),
                image_digest=image_digest,
                parameters_sha256=parameters_sha256,
                shards=tuple(
                    replace(
                        shard,
                        bcf_path=str(release / Path(shard.bcf_path).name),
                        csi_path=str(release / Path(shard.csi_path).name),
                    )
                    for shard in shards
                ),
                reference_root=root,
            ),
            manifest_path,
            root,
        )
        published_shards = tuple(
            replace(
                shard,
                bcf_path=str(release / Path(shard.bcf_path).name),
                csi_path=str(release / Path(shard.csi_path).name),
            )
            for shard in shards
        )
        _remove_failed_run(layout.run_root, scratch)
        _remove_failed_run(accepted.pilot_root, scratch)
        return _result(
            contract=contract,
            shards=published_shards,
            reserve_bytes=reserve_bytes,
            manifest_path=publication.path,
            manifest_sha256=publication.manifest_sha256,
            reused=publication.reused,
        )


def full_publication_result_as_json(result: GnomadFullPublicationResult) -> str:
    return json.dumps(asdict(result), indent=2, sort_keys=True) + "\n"
