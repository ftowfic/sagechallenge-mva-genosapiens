"""Pinned local gnomAD v4.1 all-sites allele-number lookups.

The all-sites tables are intentionally separate from the reduced joint BCF:
they establish locus callability when a variant record is absent, but do not
turn an absent locus or a missing AN into an alternate count of zero.
"""

from __future__ import annotations

from collections.abc import Callable, Iterator, Mapping, Sequence
from contextlib import contextmanager
from dataclasses import asdict, dataclass
from datetime import date
from enum import StrEnum
import fcntl
import gzip
import hashlib
import json
import os
from pathlib import Path, PurePosixPath
import re
import shutil
import subprocess
from typing import Any, BinaryIO
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
from genosapiens.storage import CapacityPlan, DEFAULT_RESERVE_BYTES, require_capacity


SCHEMA = "genosapiens.gnomad-all-sites"
SOURCE_MANIFEST_TRANSFORMATION_SCHEMA = "genosapiens.gnomad-all-sites-transform"
BUNDLE_SCHEMA = "genosapiens.gnomad-all-sites-bundle"
BUNDLE_TRANSFORMATION_SCHEMA = "genosapiens.gnomad-all-sites-bundle"
SOURCE_NAMES = ("exomes", "genomes")
EXOMES_SOURCE_HEADER = (
    "locus\tAN\toutside_broad_capture_region\toutside_ukb_capture_region\t"
    "outside_broad_calling_region\toutside_ukb_calling_region"
)
GENOMES_SOURCE_HEADER = "locus\tAN"
DERIVED_HEADER = (
    "#contig\tposition\tAN\toutside_broad_capture_region\t"
    "outside_ukb_capture_region\toutside_broad_calling_region\t"
    "outside_ukb_calling_region"
)
FLAG_NAMES = (
    "outside_broad_capture_region",
    "outside_ukb_capture_region",
    "outside_broad_calling_region",
    "outside_ukb_calling_region",
)
EXPECTED_SOURCES = {
    "exomes": {
        "source_path": "release/4.1/tsv/exomes/gnomad.exomes.v4.1.allele_number_all_sites.tsv.bgz",
        "url": (
            "https://storage.googleapis.com/download/storage/v1/b/"
            "gcp-public-data--gnomad/o/release%2F4.1%2Ftsv%2Fexomes%2F"
            "gnomad.exomes.v4.1.allele_number_all_sites.tsv.bgz?"
            "generation=1713312294386520&alt=media"
        ),
        "object_generation": "1713312294386520",
        "size_bytes": 1_149_738_480,
        "md5": "ee71ce1ccfa5c1d9dd86e1ee1b1d11e2",
        "source_header": EXOMES_SOURCE_HEADER,
    },
    "genomes": {
        "source_path": "release/4.1/tsv/genomes/gnomad.genomes.v4.1.allele_number_all_sites.tsv.bgz",
        "url": (
            "https://storage.googleapis.com/download/storage/v1/b/"
            "gcp-public-data--gnomad/o/release%2F4.1%2Ftsv%2Fgenomes%2F"
            "gnomad.genomes.v4.1.allele_number_all_sites.tsv.bgz?"
            "generation=1713312294896637&alt=media"
        ),
        "object_generation": "1713312294896637",
        "size_bytes": 12_012_641_029,
        "md5": "7101516dd79d48d10c28fa548b22884a",
        "source_header": GENOMES_SOURCE_HEADER,
    },
}
EXPECTED_OUTPUT_DIRECTORY = "derived/gnomad/4.1/all-sites/gnomad-v4.1-all-sites-an-v1"
EXPECTED_BUNDLE_MANIFEST_PATH = "manifests/gnomad/4.1/all-sites-an-v1.json"
MAXIMUM_DERIVED_BYTES = 16_000_000_000
SHA256 = re.compile(r"^[0-9a-f]{64}$")
MD5 = re.compile(r"^[0-9a-f]{32}$")
DIGEST = re.compile(r"^sha256:[0-9a-f]{64}$")
LOCUS = re.compile(r"^(chr[^:\t]+):([1-9][0-9]*)$")


class GnomadAllSitesError(RuntimeError):
    """Raised when the all-sites resource cannot be safely used."""


def _mapping(value: Any, context: str) -> Mapping[str, Any]:
    if not isinstance(value, Mapping) or not all(isinstance(key, str) for key in value):
        raise GnomadAllSitesError(f"{context} must be an object")
    return value


def _string(data: Mapping[str, Any], field: str, context: str) -> str:
    value = data.get(field)
    if not isinstance(value, str) or not value:
        raise GnomadAllSitesError(f"{context}.{field} must be a non-empty string")
    return value


def _positive_int(data: Mapping[str, Any], field: str, context: str) -> int:
    value = data.get(field)
    if not isinstance(value, int) or isinstance(value, bool) or value <= 0:
        raise GnomadAllSitesError(f"{context}.{field} must be a positive integer")
    return value


def _safe_relative(value: str, context: str) -> str:
    path = PurePosixPath(value)
    if path.is_absolute() or path == PurePosixPath(".") or ".." in path.parts:
        raise GnomadAllSitesError(f"{context} must be a safe relative path")
    return str(path)


@dataclass(frozen=True, slots=True)
class GnomadAllSitesSource:
    name: str
    source_path: str
    url: str
    object_generation: str
    size_bytes: int
    md5: str
    source_header: str

    @classmethod
    def from_mapping(cls, data: Mapping[str, Any]) -> "GnomadAllSitesSource":
        expected = {
            "name",
            "source_path",
            "url",
            "object_generation",
            "size_bytes",
            "md5",
            "source_header",
        }
        if set(data) != expected:
            raise GnomadAllSitesError("all-sites source fields are invalid")
        source = cls(
            name=_string(data, "name", "all-sites source"),
            source_path=_string(data, "source_path", "all-sites source"),
            url=_string(data, "url", "all-sites source"),
            object_generation=_string(data, "object_generation", "all-sites source"),
            size_bytes=_positive_int(data, "size_bytes", "all-sites source"),
            md5=_string(data, "md5", "all-sites source"),
            source_header=_string(data, "source_header", "all-sites source"),
        )
        if (
            source.name not in SOURCE_NAMES
            or not source.url.startswith("https://")
            or not source.object_generation.isdecimal()
            or not MD5.fullmatch(source.md5)
            or source.source_path != EXPECTED_SOURCES[source.name]["source_path"]
            or source.url != EXPECTED_SOURCES[source.name]["url"]
            or source.object_generation
            != EXPECTED_SOURCES[source.name]["object_generation"]
            or source.size_bytes != EXPECTED_SOURCES[source.name]["size_bytes"]
            or source.md5 != EXPECTED_SOURCES[source.name]["md5"]
            or source.source_header
            != EXPECTED_SOURCES[source.name]["source_header"]
        ):
            raise GnomadAllSitesError("all-sites source is not the pinned v4.1 object")
        return source

    @property
    def output_name(self) -> str:
        return f"gnomad.{self.name}.v4.1.allele_number_all_sites.lookup.tsv.bgz"


@dataclass(frozen=True, slots=True)
class GnomadAllSitesContract:
    resource_id: str
    release: str
    assembly: str
    license: str
    license_url: str
    derived_header: str
    output_directory: str
    bundle_manifest_path: str
    maximum_derived_bytes: int
    sources: tuple[GnomadAllSitesSource, ...]
    schema: str = SCHEMA

    @classmethod
    def from_mapping(cls, data: Mapping[str, Any]) -> "GnomadAllSitesContract":
        expected = {
            "schema",
            "resource_id",
            "release",
            "assembly",
            "license",
            "license_url",
            "derived_header",
            "output_directory",
            "bundle_manifest_path",
            "maximum_derived_bytes",
            "sources",
        }
        sources = data.get("sources")
        if (
            set(data) != expected
            or data.get("schema") != SCHEMA
            or not isinstance(sources, list)
            or len(sources) != len(SOURCE_NAMES)
            or not all(isinstance(item, Mapping) for item in sources)
        ):
            raise GnomadAllSitesError("all-sites contract fields are invalid")
        contract = cls(
            resource_id=_string(data, "resource_id", "all-sites contract"),
            release=_string(data, "release", "all-sites contract"),
            assembly=_string(data, "assembly", "all-sites contract"),
            license=_string(data, "license", "all-sites contract"),
            license_url=_string(data, "license_url", "all-sites contract"),
            derived_header=_string(data, "derived_header", "all-sites contract"),
            output_directory=_safe_relative(
                _string(data, "output_directory", "all-sites contract"),
                "all-sites output_directory",
            ),
            bundle_manifest_path=_safe_relative(
                _string(data, "bundle_manifest_path", "all-sites contract"),
                "all-sites bundle_manifest_path",
            ),
            maximum_derived_bytes=_positive_int(
                data, "maximum_derived_bytes", "all-sites contract"
            ),
            sources=tuple(GnomadAllSitesSource.from_mapping(item) for item in sources),
        )
        if (
            contract.resource_id != "gnomad-v4.1-all-sites-an"
            or contract.release != "v4.1"
            or contract.assembly != "GRCh38"
            or contract.license != "ODC-ODbL-1.0"
            or contract.license_url
            != "https://opendatacommons.org/licenses/odbl/1-0/"
            or contract.derived_header != DERIVED_HEADER
            or contract.output_directory != EXPECTED_OUTPUT_DIRECTORY
            or PurePosixPath(contract.bundle_manifest_path).parent
            != PurePosixPath(EXPECTED_BUNDLE_MANIFEST_PATH).parent
            or PurePosixPath(contract.bundle_manifest_path).suffix != ".json"
            or contract.maximum_derived_bytes != MAXIMUM_DERIVED_BYTES
            or tuple(source.name for source in contract.sources) != SOURCE_NAMES
        ):
            raise GnomadAllSitesError("all-sites contract is outside the fixed v4.1 scope")
        return contract

    def source(self, name: str) -> GnomadAllSitesSource:
        for source in self.sources:
            if source.name == name:
                return source
        raise GnomadAllSitesError(f"unknown all-sites source: {name}")


def load_gnomad_all_sites_contract(path: str | Path) -> GnomadAllSitesContract:
    try:
        data = json.loads(Path(path).read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as error:
        raise GnomadAllSitesError(f"cannot load all-sites contract: {error}") from error
    return GnomadAllSitesContract.from_mapping(_mapping(data, "all-sites contract"))


@dataclass(frozen=True, slots=True)
class AllSitesRegionFlags:
    outside_broad_capture_region: bool | None
    outside_ukb_capture_region: bool | None
    outside_broad_calling_region: bool | None
    outside_ukb_calling_region: bool | None


@dataclass(frozen=True, slots=True)
class AllSitesDatasetCoverage:
    an: int | None
    flags: AllSitesRegionFlags

    @property
    def callable(self) -> bool:
        return self.an is not None and self.an > 0


@dataclass(frozen=True, slots=True)
class AllSitesLocusCoverage:
    locus: str
    exomes: AllSitesDatasetCoverage | None
    genomes: AllSitesDatasetCoverage | None

    @property
    def joint_an(self) -> int | None:
        """Return an exact joint AN only when both source ANs are defined."""
        if (
            self.exomes is None
            or self.genomes is None
            or self.exomes.an is None
            or self.genomes.an is None
        ):
            return None
        return self.exomes.an + self.genomes.an

    @property
    def callable(self) -> bool:
        return any(
            value is not None and value.callable
            for value in (self.exomes, self.genomes)
        )


class GnomadAlleleState(StrEnum):
    OBSERVED_ALTERNATE = "observed_alternate"
    AC0_FILTERED_OBSERVATION = "ac0_filtered_observation"
    CALLABLE_ZERO_ALTERNATE = "callable_zero_alternate"
    UNKNOWN_NO_CALLABLE_DATA = "unknown_no_callable_data"


def classify_gnomad_allele(
    *,
    exact_observed_alternate: bool,
    ac0_filtered_observation: bool,
    coverage: AllSitesLocusCoverage,
) -> GnomadAlleleState:
    """Classify population evidence without treating an absent row as AC=0."""
    if not isinstance(exact_observed_alternate, bool) or not isinstance(
        ac0_filtered_observation, bool
    ):
        raise GnomadAllSitesError("gnomAD observation indicators must be Boolean")
    if exact_observed_alternate and ac0_filtered_observation:
        raise GnomadAllSitesError("an allele cannot be observed and AC0-filtered")
    if exact_observed_alternate:
        return GnomadAlleleState.OBSERVED_ALTERNATE
    if ac0_filtered_observation:
        return GnomadAlleleState.AC0_FILTERED_OBSERVATION
    if coverage.callable:
        return GnomadAlleleState.CALLABLE_ZERO_ALTERNATE
    return GnomadAlleleState.UNKNOWN_NO_CALLABLE_DATA


@dataclass(frozen=True, slots=True)
class AllSitesTransform:
    source_name: str
    source_sha256: str
    source_records: int
    derived_records: int
    locus_order_sha256: str
    contigs: tuple[str, ...]

    @classmethod
    def from_mapping(cls, data: Mapping[str, Any]) -> "AllSitesTransform":
        expected = {
            "source_name",
            "source_sha256",
            "source_records",
            "derived_records",
            "locus_order_sha256",
            "contigs",
        }
        contigs = data.get("contigs")
        if set(data) != expected or not isinstance(contigs, list):
            raise GnomadAllSitesError("all-sites transform record is invalid")
        value = cls(
            source_name=_string(data, "source_name", "all-sites transform"),
            source_sha256=_string(data, "source_sha256", "all-sites transform"),
            source_records=_positive_int(data, "source_records", "all-sites transform"),
            derived_records=_positive_int(data, "derived_records", "all-sites transform"),
            locus_order_sha256=_string(
                data, "locus_order_sha256", "all-sites transform"
            ),
            contigs=tuple(contigs),
        )
        if (
            value.source_name not in SOURCE_NAMES
            or not SHA256.fullmatch(value.source_sha256)
            or not SHA256.fullmatch(value.locus_order_sha256)
            or value.source_records != value.derived_records
            or not value.contigs
            or any(not isinstance(contig, str) or not contig.startswith("chr") for contig in value.contigs)
            or len(value.contigs) != len(set(value.contigs))
        ):
            raise GnomadAllSitesError("all-sites transform values are invalid")
        return value


@dataclass(frozen=True, slots=True)
class GnomadAllSitesBootstrapResult:
    resource_id: str
    release: str
    output_directory: str
    bundle_manifest_path: str
    bundle_manifest_sha256: str
    transforms: tuple[AllSitesTransform, ...]
    reused: bool


Downloader = Callable[[DownloadRequest, str | Path], DownloadResult]
CapacityChecker = Callable[[str | Path, CapacityPlan], Any]
CommandRunner = Callable[[Sequence[str]], str]


def _sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    try:
        with path.open("rb") as handle:
            for chunk in iter(lambda: handle.read(1024 * 1024), b""):
                digest.update(chunk)
    except OSError as error:
        raise GnomadAllSitesError(f"cannot hash {path}: {error}") from error
    return digest.hexdigest()


def _fsync_directory(path: Path) -> None:
    try:
        descriptor = os.open(path, os.O_RDONLY)
        try:
            os.fsync(descriptor)
        finally:
            os.close(descriptor)
    except OSError as error:
        raise GnomadAllSitesError(f"cannot sync directory {path}: {error}") from error


def _require_within_root(root: Path, path: Path, context: str) -> None:
    if not path.is_absolute() or not path.resolve(strict=False).is_relative_to(root):
        raise GnomadAllSitesError(f"{context} escapes reference_root")
    current = root
    for part in path.relative_to(root).parts:
        current /= part
        if current.is_symlink():
            raise GnomadAllSitesError(f"{context} contains a symlink: {current}")


def _parse_locus(value: str) -> tuple[str, int]:
    match = LOCUS.fullmatch(value)
    if match is None:
        raise GnomadAllSitesError(f"invalid all-sites locus: {value!r}")
    return match.group(1), int(match.group(2))


def _parse_an(value: str, context: str) -> int | None:
    if value in {"", "."}:
        return None
    try:
        parsed = int(value)
    except ValueError as error:
        raise GnomadAllSitesError(f"{context} AN is invalid") from error
    if parsed < 0:
        raise GnomadAllSitesError(f"{context} AN is negative")
    return parsed


def _parse_flag(value: str, context: str) -> bool | None:
    if value in {"", "."}:
        return None
    if value in {"0", "false", "False"}:
        return False
    if value in {"1", "true", "True"}:
        return True
    raise GnomadAllSitesError(f"{context} flag is invalid")


def _format_flag(value: bool | None) -> str:
    if value is None:
        return "."
    return "1" if value else "0"


def _parse_source_record(
    line: str, source: GnomadAllSitesSource
) -> tuple[str, int, int | None, AllSitesRegionFlags]:
    fields = line.split("\t")
    header_fields = source.source_header.split("\t")
    if len(fields) != len(header_fields):
        raise GnomadAllSitesError("all-sites source record does not match its header")
    contig, position = _parse_locus(fields[0])
    flags = AllSitesRegionFlags(
        **{
            name: (
                _parse_flag(fields[header_fields.index(name)], f"all-sites {contig}:{position}")
                if name in header_fields
                else None
            )
            for name in FLAG_NAMES
        }
    )
    return contig, position, _parse_an(fields[1], f"all-sites {contig}:{position}"), flags


def _parse_derived_record(line: str) -> tuple[str, int, AllSitesDatasetCoverage]:
    fields = line.split("\t")
    if len(fields) != 7:
        raise GnomadAllSitesError("derived all-sites record does not have seven columns")
    try:
        position = int(fields[1])
    except ValueError as error:
        raise GnomadAllSitesError("derived all-sites position is invalid") from error
    if not fields[0].startswith("chr") or position <= 0:
        raise GnomadAllSitesError("derived all-sites locus is invalid")
    flags = AllSitesRegionFlags(
        **{
            name: _parse_flag(value, f"derived all-sites {fields[0]}:{position}")
            for name, value in zip(FLAG_NAMES, fields[3:], strict=True)
        }
    )
    return fields[0], position, AllSitesDatasetCoverage(
        an=_parse_an(fields[2], f"derived all-sites {fields[0]}:{position}"),
        flags=flags,
    )


def _write_derived_record(
    output: BinaryIO,
    contig: str,
    position: int,
    an: int | None,
    flags: AllSitesRegionFlags,
) -> None:
    output.write(
        (
            "\t".join(
                (
                    contig,
                    str(position),
                    "." if an is None else str(an),
                    *(_format_flag(getattr(flags, name)) for name in FLAG_NAMES),
                )
            )
            + "\n"
        ).encode("ascii")
    )


def _run_command(arguments: Sequence[str]) -> str:
    try:
        completed = subprocess.run(
            list(arguments),
            check=False,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
        )
    except OSError as error:
        raise GnomadAllSitesError(f"cannot execute {arguments[0]}: {error}") from error
    if completed.returncode != 0:
        detail = completed.stderr.strip().splitlines()
        raise GnomadAllSitesError(
            f"{arguments[0]} failed: {detail[-1] if detail else completed.returncode}"
        )
    return completed.stdout


def _require_index_tools() -> None:
    _run_command(("bgzip", "--version"))
    _run_command(("tabix", "--version"))


def _stop_process(process: subprocess.Popen[bytes] | None) -> None:
    if process is None or process.poll() is not None:
        return
    process.kill()
    process.wait()


def _process_error(process: subprocess.Popen[bytes], tool: str) -> None:
    stderr = process.stderr.read().decode("utf-8", errors="replace").strip()
    if process.wait() != 0:
        raise GnomadAllSitesError(f"{tool} failed: {stderr or 'unknown error'}")


def _transform_source(
    source_path: Path,
    output_path: Path,
    source: GnomadAllSitesSource,
    source_sha256: str,
    *,
    command_runner: CommandRunner = _run_command,
) -> AllSitesTransform:
    """Stream one source through bgzip without materializing an uncompressed TSV."""
    del command_runner  # Kept injectable for a uniform bootstrap test seam.
    decompressor: subprocess.Popen[bytes] | None = None
    compressor: subprocess.Popen[bytes] | None = None
    source_records = 0
    derived_records = 0
    order = hashlib.sha256()
    contigs: list[str] = []
    seen_contigs: set[str] = set()
    current_contig: str | None = None
    previous_position = 0
    try:
        output_path.parent.mkdir(parents=True, exist_ok=True)
        with output_path.open("xb") as compressed_output:
            decompressor = subprocess.Popen(
                ("bgzip", "--decompress", "--stdout", str(source_path)),
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
            )
            compressor = subprocess.Popen(
                ("bgzip", "--stdout", "--threads", "1"),
                stdin=subprocess.PIPE,
                stdout=compressed_output,
                stderr=subprocess.PIPE,
            )
            assert decompressor.stdout is not None
            assert compressor.stdin is not None
            reader = decompressor.stdout
            header = reader.readline()
            try:
                header_text = header.rstrip(b"\r\n").decode("ascii")
            except UnicodeDecodeError as error:
                raise GnomadAllSitesError("all-sites source header is not ASCII") from error
            if header_text != source.source_header:
                raise GnomadAllSitesError("all-sites source header differs from the contract")
            compressor.stdin.write((DERIVED_HEADER + "\n").encode("ascii"))
            for raw_line in reader:
                try:
                    line = raw_line.rstrip(b"\r\n").decode("ascii")
                except UnicodeDecodeError as error:
                    raise GnomadAllSitesError("all-sites source record is not ASCII") from error
                if not line:
                    raise GnomadAllSitesError("all-sites source contains a blank record")
                contig, position, an, flags = _parse_source_record(line, source)
                if contig != current_contig:
                    if contig in seen_contigs:
                        raise GnomadAllSitesError("all-sites source contigs are not grouped")
                    seen_contigs.add(contig)
                    contigs.append(contig)
                    current_contig = contig
                    previous_position = 0
                if position <= previous_position:
                    raise GnomadAllSitesError("all-sites source is not position sorted")
                previous_position = position
                source_records += 1
                order.update(f"{contig}:{position}\n".encode("ascii"))
                _write_derived_record(
                    compressor.stdin,
                    contig,
                    position,
                    an,
                    flags,
                )
                derived_records += 1
            compressor.stdin.close()
            _process_error(compressor, "bgzip compression")
            _process_error(decompressor, "bgzip decompression")
            compressed_output.flush()
            os.fsync(compressed_output.fileno())
        if not source_records or source_records != derived_records:
            raise GnomadAllSitesError("all-sites source-to-derived record count differs")
        _run_command(
            (
                "tabix",
                "--force",
                "--sequence",
                "1",
                "--begin",
                "2",
                "--end",
                "2",
                "--comment",
                "#",
                str(output_path),
            )
        )
        index = Path(f"{output_path}.tbi")
        if index.is_symlink() or not index.is_file():
            raise GnomadAllSitesError("tabix did not create an all-sites index")
        return AllSitesTransform(
            source_name=source.name,
            source_sha256=source_sha256,
            source_records=source_records,
            derived_records=derived_records,
            locus_order_sha256=order.hexdigest(),
            contigs=tuple(contigs),
        )
    except BaseException:
        _stop_process(compressor)
        _stop_process(decompressor)
        output_path.unlink(missing_ok=True)
        Path(f"{output_path}.tbi").unlink(missing_ok=True)
        raise


def _scan_derived_lookup(path: Path) -> tuple[int, str, tuple[str, ...], str]:
    count = 0
    order = hashlib.sha256()
    contigs: list[str] = []
    seen_contigs: set[str] = set()
    current_contig: str | None = None
    previous_position = 0
    first_line: str | None = None
    try:
        with gzip.open(path, "rt", encoding="ascii", newline="") as source:
            header = source.readline().rstrip("\r\n")
            if header != DERIVED_HEADER:
                raise GnomadAllSitesError("derived all-sites header differs from the contract")
            for line in source:
                stripped = line.rstrip("\r\n")
                if not stripped:
                    raise GnomadAllSitesError("derived all-sites lookup contains a blank record")
                contig, position, _coverage = _parse_derived_record(stripped)
                if contig != current_contig:
                    if contig in seen_contigs:
                        raise GnomadAllSitesError("derived all-sites contigs are not grouped")
                    seen_contigs.add(contig)
                    contigs.append(contig)
                    current_contig = contig
                    previous_position = 0
                if position <= previous_position:
                    raise GnomadAllSitesError("derived all-sites lookup is not position sorted")
                previous_position = position
                count += 1
                order.update(f"{contig}:{position}\n".encode("ascii"))
                first_line = first_line or stripped
    except (OSError, EOFError, gzip.BadGzipFile, UnicodeDecodeError) as error:
        raise GnomadAllSitesError(f"cannot read derived all-sites lookup: {error}") from error
    if not count or first_line is None:
        raise GnomadAllSitesError("derived all-sites lookup has no records")
    return count, order.hexdigest(), tuple(contigs), first_line


def _validate_derived_lookup(
    path: Path,
    transform: AllSitesTransform,
    *,
    command_runner: CommandRunner = _run_command,
) -> None:
    if path.is_symlink() or not path.is_file():
        raise GnomadAllSitesError(f"derived all-sites lookup is missing: {path}")
    index = Path(f"{path}.tbi")
    if index.is_symlink() or not index.is_file():
        raise GnomadAllSitesError(f"derived all-sites index is missing: {index}")
    count, order, contigs, first_line = _scan_derived_lookup(path)
    if (
        count != transform.source_records
        or count != transform.derived_records
        or order != transform.locus_order_sha256
        or contigs != transform.contigs
    ):
        raise GnomadAllSitesError("source-to-derived all-sites record order differs")
    indexed_contigs = tuple(
        line for line in command_runner(("tabix", "--list-chroms", str(path))).splitlines() if line
    )
    if indexed_contigs != contigs:
        raise GnomadAllSitesError("all-sites tabix contigs differ from the lookup")
    first_contig, first_position, _coverage = _parse_derived_record(first_line)
    query = command_runner(
        ("tabix", str(path), f"{first_contig}:{first_position}-{first_position}")
    )
    if query.splitlines() != [first_line]:
        raise GnomadAllSitesError("all-sites tabix query did not recover the first record")


def _source_manifest_path(contract: GnomadAllSitesContract, name: str) -> str:
    bundle = PurePosixPath(contract.bundle_manifest_path)
    return str(bundle.with_name(f"{bundle.stem}.{name}.json"))


def _source_manifest(
    *,
    contract: GnomadAllSitesContract,
    source: GnomadAllSitesSource,
    output_path: Path,
    root: Path,
    contract_sha256: str,
    image_digest: str,
    retrieved_on: str,
) -> PublicResourceManifest:
    paths = (output_path, Path(f"{output_path}.tbi"))
    return PublicResourceManifest(
        resource_id=f"gnomad-v4.1-all-sites-an-{source.name}",
        release=contract.release,
        source=SourceIdentity(source.url, object_generation=source.object_generation),
        license=contract.license,
        assembly=contract.assembly,
        retrieved_on=retrieved_on,
        transformation=Transformation(SOURCE_MANIFEST_TRANSFORMATION_SCHEMA, contract_sha256),
        image_digest=image_digest,
        outputs=tuple(
            ResourceOutput(
                path=str(path.relative_to(root)),
                size_bytes=path.stat().st_size,
                sha256=_sha256_file(path),
                retention=RetentionPolicy("durable", None),
            )
            for path in paths
        ),
    )


def _verify_source_manifest_contract(manifest, contract, source, contract_sha256):
    output_path = PurePosixPath(contract.output_directory) / source.output_name
    expected_outputs = {str(output_path), str(output_path) + ".tbi"}
    if (
        manifest.resource_id != f"gnomad-v4.1-all-sites-an-{source.name}"
        or manifest.release != contract.release
        or manifest.assembly != contract.assembly
        or manifest.license != contract.license
        or manifest.source != SourceIdentity(source.url, object_generation=source.object_generation)
        or manifest.transformation != Transformation(SOURCE_MANIFEST_TRANSFORMATION_SCHEMA, contract_sha256)
        or {item.path for item in manifest.outputs} != expected_outputs
        or any(item.retention != RetentionPolicy("durable", None) for item in manifest.outputs)
    ):
        raise GnomadAllSitesError("all-sites source manifest differs from the contract")


def _validate_source_manifest(
    *,
    manifest_path: Path,
    contract: GnomadAllSitesContract,
    source: GnomadAllSitesSource,
    root: Path,
    contract_sha256: str,
    transform: AllSitesTransform,
    command_runner: CommandRunner,
) -> None:
    try:
        manifest = PublicResourceManifest.from_mapping(
            _mapping(json.loads(manifest_path.read_text(encoding="utf-8")), "source manifest")
        )
    except (OSError, json.JSONDecodeError, ManifestError) as error:
        raise GnomadAllSitesError(f"cannot load all-sites source manifest: {error}") from error
    output_path = root / contract.output_directory / source.output_name
    _verify_source_manifest_contract(manifest, contract, source, contract_sha256)
    try:
        verify_outputs(manifest, root)
    except ManifestError as error:
        raise GnomadAllSitesError(f"all-sites source outputs are invalid: {error}") from error
    _validate_derived_lookup(output_path, transform, command_runner=command_runner)


def _canonical_json(value: Mapping[str, Any]) -> bytes:
    return (json.dumps(value, sort_keys=True, separators=(",", ":")) + "\n").encode("utf-8")


def _bundle_document(
    *,
    contract: GnomadAllSitesContract,
    contract_sha256: str,
    image_digest: str,
    retrieved_on: str,
    transforms: Sequence[AllSitesTransform],
    source_manifest_sha256s: Mapping[str, str],
) -> dict[str, Any]:
    by_name = {transform.source_name: transform for transform in transforms}
    return {
        "schema": BUNDLE_SCHEMA,
        "resource_id": contract.resource_id,
        "release": contract.release,
        "assembly": contract.assembly,
        "license": contract.license,
        "license_url": contract.license_url,
        "retrieved_on": retrieved_on,
        "image_digest": image_digest,
        "transformation": {
            "schema": BUNDLE_TRANSFORMATION_SCHEMA,
            "parameters_sha256": contract_sha256,
        },
        "output_directory": contract.output_directory,
        "sources": [
            {
                "name": source.name,
                "source_path": source.source_path,
                "url": source.url,
                "object_generation": source.object_generation,
                "size_bytes": source.size_bytes,
                "md5": source.md5,
                "source_manifest_path": _source_manifest_path(contract, source.name),
                "source_manifest_sha256": source_manifest_sha256s[source.name],
                "transform": asdict(by_name[source.name]),
            }
            for source in contract.sources
        ],
    }


def _load_bundle_document(path: Path) -> Mapping[str, Any]:
    try:
        document = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as error:
        raise GnomadAllSitesError(f"cannot load all-sites bundle manifest: {error}") from error
    return _mapping(document, "all-sites bundle manifest")


def validate_bundle_metadata(document, contract, contract_sha256):
    """Validate the current bundle and every declared source before data I/O."""
    expected = {
        "schema",
        "resource_id",
        "release",
        "assembly",
        "license",
        "license_url",
        "retrieved_on",
        "image_digest",
        "transformation",
        "output_directory",
        "sources",
    }
    sources = document.get("sources")
    transformation = document.get("transformation")
    if (
        set(document) != expected
        or document.get("schema") != BUNDLE_SCHEMA
        or document.get("resource_id") != contract.resource_id
        or document.get("release") != contract.release
        or document.get("assembly") != contract.assembly
        or document.get("license") != contract.license
        or document.get("license_url") != contract.license_url
        or document.get("output_directory") != contract.output_directory
        or not isinstance(document.get("retrieved_on"), str)
        or not DIGEST.fullmatch(str(document.get("image_digest")))
        or not isinstance(transformation, Mapping)
        or transformation
        != {
            "schema": BUNDLE_TRANSFORMATION_SCHEMA,
            "parameters_sha256": contract_sha256,
        }
        or not isinstance(sources, list)
        or len(sources) != len(SOURCE_NAMES)
    ):
        raise GnomadAllSitesError("all-sites bundle manifest differs from the contract")
    entries = []
    expected_source_fields = {
        "name",
        "source_path",
        "url",
        "object_generation",
        "size_bytes",
        "md5",
        "source_manifest_path",
        "source_manifest_sha256",
        "transform",
    }
    for expected_source, value in zip(contract.sources, sources, strict=True):
        source_data = _mapping(value, "all-sites bundle source")
        if (
            set(source_data) != expected_source_fields
            or any(
                source_data.get(field) != getattr(expected_source, field)
                for field in ("name", "source_path", "url", "object_generation", "size_bytes", "md5")
            )
            or source_data.get("source_manifest_path")
            != _source_manifest_path(contract, expected_source.name)
            or not isinstance(source_data.get("source_manifest_sha256"), str)
            or not SHA256.fullmatch(source_data["source_manifest_sha256"])
            or not isinstance(source_data.get("transform"), Mapping)
        ):
            raise GnomadAllSitesError("all-sites bundle source identity is invalid")
        transform = AllSitesTransform.from_mapping(source_data["transform"])
        if transform.source_name != expected_source.name:
            raise GnomadAllSitesError("all-sites bundle transform source differs")
        entries.append((source_data, transform))
    return tuple(entries)


def _validate_bundle(
    *,
    bundle_path: Path,
    contract: GnomadAllSitesContract,
    root: Path,
    contract_sha256: str,
    command_runner: CommandRunner,
) -> tuple[AllSitesTransform, ...]:
    entries = validate_bundle_metadata(_load_bundle_document(bundle_path), contract, contract_sha256)
    transforms = []
    for expected_source, (source_data, transform) in zip(contract.sources, entries, strict=True):
        source_manifest_path = root / str(source_data["source_manifest_path"])
        _require_within_root(root, source_manifest_path, "all-sites source manifest")
        if _sha256_file(source_manifest_path) != source_data["source_manifest_sha256"]:
            raise GnomadAllSitesError("all-sites source manifest checksum differs")
        _validate_source_manifest(
            manifest_path=source_manifest_path,
            contract=contract,
            source=expected_source,
            root=root,
            contract_sha256=contract_sha256,
            transform=transform,
            command_runner=command_runner,
        )
        transforms.append(transform)
    return tuple(transforms)


def _publish_bundle(path: Path, document: Mapping[str, Any], root: Path) -> str:
    content = _canonical_json(document)
    digest = hashlib.sha256(content).hexdigest()
    if path.exists() or path.is_symlink():
        if path.is_symlink() or not path.is_file() or path.read_bytes() != content:
            raise GnomadAllSitesError("completed all-sites bundle manifest is immutable")
        return digest
    path.parent.mkdir(parents=True, exist_ok=True)
    partial_root = path.parent / ".partial"
    if partial_root.is_symlink():
        raise GnomadAllSitesError("all-sites bundle manifest partial path is a symlink")
    partial_root.mkdir(mode=0o700, exist_ok=True)
    partial = partial_root / f"{path.name}.{uuid.uuid4().hex}.partial"
    try:
        with partial.open("xb") as handle:
            os.chmod(partial, 0o600)
            handle.write(content)
            handle.flush()
            os.fsync(handle.fileno())
        os.chmod(partial, 0o444)
        try:
            os.link(partial, path)
        except FileExistsError:
            if path.is_symlink() or not path.is_file() or path.read_bytes() != content:
                raise GnomadAllSitesError("completed all-sites bundle manifest is immutable")
        _fsync_directory(path.parent)
    except OSError as error:
        raise GnomadAllSitesError(f"cannot publish all-sites bundle manifest: {error}") from error
    finally:
        partial.unlink(missing_ok=True)
    return digest


@contextmanager
def _bootstrap_lock(scratch_root: Path) -> Iterator[None]:
    scratch_root.mkdir(parents=True, exist_ok=True, mode=0o700)
    lock_path = scratch_root / ".gnomad-all-sites.lock"
    try:
        descriptor = os.open(lock_path, os.O_RDWR | os.O_CREAT | os.O_NOFOLLOW, 0o600)
    except OSError as error:
        raise GnomadAllSitesError(f"cannot claim all-sites bootstrap lock: {error}") from error
    try:
        try:
            fcntl.flock(descriptor, fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError as error:
            raise GnomadAllSitesError("all-sites bootstrap is already active") from error
        yield
    finally:
        os.close(descriptor)


def _make_release_read_only(release: Path, contract: GnomadAllSitesContract) -> None:
    for source in contract.sources:
        output = release / source.output_name
        os.chmod(output, 0o444)
        os.chmod(Path(f"{output}.tbi"), 0o444)
    os.chmod(release, 0o755)
    _fsync_directory(release)


def bootstrap_gnomad_all_sites(
    *,
    contract_path: str | Path,
    reference_root: str | Path,
    scratch_root: str | Path,
    reserve_bytes: int,
    image_digest: str,
    retrieved_on: str | None = None,
    downloader: Downloader = download_atomic,
    capacity_checker: CapacityChecker = require_capacity,
    command_runner: CommandRunner = _run_command,
) -> GnomadAllSitesBootstrapResult:
    """Download, transform, index, and atomically publish the two v4.1 lookups."""
    contract_file = Path(contract_path)
    root = Path(reference_root).resolve(strict=True)
    scratch = Path(scratch_root)
    contract = load_gnomad_all_sites_contract(contract_file)
    if not DIGEST.fullmatch(image_digest):
        raise GnomadAllSitesError("all-sites bootstrap image_digest must be a SHA-256 digest")
    if not isinstance(reserve_bytes, int) or isinstance(reserve_bytes, bool):
        raise GnomadAllSitesError("all-sites reserve_bytes must be an integer")
    if reserve_bytes < DEFAULT_RESERVE_BYTES:
        raise GnomadAllSitesError("all-sites bootstrap must preserve at least 100 GB")
    _require_within_root(root, scratch, "all-sites scratch_root")
    release = root / contract.output_directory
    bundle_path = root / contract.bundle_manifest_path
    _require_within_root(root, release, "all-sites output directory")
    _require_within_root(root, bundle_path, "all-sites bundle manifest")
    contract_sha256 = _sha256_file(contract_file)

    if bundle_path.exists() or bundle_path.is_symlink():
        transforms = _validate_bundle(
            bundle_path=bundle_path,
            contract=contract,
            root=root,
            contract_sha256=contract_sha256,
            command_runner=command_runner,
        )
        return GnomadAllSitesBootstrapResult(
            resource_id=contract.resource_id,
            release=contract.release,
            output_directory=str(release),
            bundle_manifest_path=str(bundle_path),
            bundle_manifest_sha256=_sha256_file(bundle_path),
            transforms=transforms,
            reused=True,
        )
    if release.exists() or release.is_symlink():
        raise GnomadAllSitesError("all-sites output exists without its immutable bundle manifest")

    _require_index_tools()
    capacity_checker(
        root,
        CapacityPlan(
            operation="bootstrap gnomAD v4.1 all-sites allele-number lookups",
            source_bytes=max(source.size_bytes for source in contract.sources),
            output_bytes=contract.maximum_derived_bytes,
            reserve_bytes=reserve_bytes,
        ),
    )
    retrieval_date = retrieved_on or date.today().isoformat()
    run_root = scratch / "gnomad-all-sites" / uuid.uuid4().hex
    _require_within_root(root, run_root, "all-sites transient directory")
    staging: Path | None = None
    published_release = False
    try:
        with _bootstrap_lock(scratch):
            if bundle_path.exists() or bundle_path.is_symlink():
                transforms = _validate_bundle(
                    bundle_path=bundle_path,
                    contract=contract,
                    root=root,
                    contract_sha256=contract_sha256,
                    command_runner=command_runner,
                )
                return GnomadAllSitesBootstrapResult(
                    resource_id=contract.resource_id,
                    release=contract.release,
                    output_directory=str(release),
                    bundle_manifest_path=str(bundle_path),
                    bundle_manifest_sha256=_sha256_file(bundle_path),
                    transforms=transforms,
                    reused=True,
                )
            if release.exists() or release.is_symlink():
                raise GnomadAllSitesError(
                    "all-sites output exists without its immutable bundle manifest"
                )
            run_root.mkdir(parents=True, mode=0o700)
            release.parent.mkdir(parents=True, exist_ok=True)
            partial_root = release.parent / ".partial"
            _require_within_root(root, partial_root, "all-sites publication staging")
            partial_root.mkdir(mode=0o700, exist_ok=True)
            staging = partial_root / f"{release.name}.{uuid.uuid4().hex}.partial"
            staging.mkdir(mode=0o700)
            transforms = []
            for source in contract.sources:
                raw_path = run_root / f"{source.name}.source.tsv.gz"
                output_path = staging / source.output_name
                try:
                    downloaded = downloader(
                        DownloadRequest(
                            url=source.url,
                            destination=raw_path,
                            expected_size_bytes=source.size_bytes,
                            expected_md5=source.md5,
                            expected_generation=source.object_generation,
                            attempts=3,
                            timeout_seconds=300,
                        ),
                        root,
                    )
                    if (
                        downloaded.destination != raw_path
                        or downloaded.size_bytes != source.size_bytes
                        or downloaded.md5 != source.md5
                        or downloaded.source_generation != source.object_generation
                        or not SHA256.fullmatch(downloaded.sha256)
                    ):
                        raise GnomadAllSitesError("all-sites download identity differs")
                    transform = _transform_source(
                        raw_path,
                        output_path,
                        source,
                        downloaded.sha256,
                        command_runner=command_runner,
                    )
                    _validate_derived_lookup(
                        output_path, transform, command_runner=command_runner
                    )
                    transforms.append(transform)
                finally:
                    raw_path.unlink(missing_ok=True)
            if sum(
                path.stat().st_size for path in staging.iterdir() if path.is_file()
            ) > contract.maximum_derived_bytes:
                raise GnomadAllSitesError("derived all-sites lookups exceed their capacity cap")
            _make_release_read_only(staging, contract)
            os.replace(staging, release)
            staging = None
            published_release = True
            _fsync_directory(release.parent)

            source_manifest_sha256s: dict[str, str] = {}
            for source, transform in zip(contract.sources, transforms, strict=True):
                output_path = release / source.output_name
                publication = publish_manifest(
                    _source_manifest(
                        contract=contract,
                        source=source,
                        output_path=output_path,
                        root=root,
                        contract_sha256=contract_sha256,
                        image_digest=image_digest,
                        retrieved_on=retrieval_date,
                    ),
                    root / _source_manifest_path(contract, source.name),
                    root,
                )
                source_manifest_sha256s[source.name] = publication.manifest_sha256
                _validate_source_manifest(
                    manifest_path=publication.path,
                    contract=contract,
                    source=source,
                    root=root,
                    contract_sha256=contract_sha256,
                    transform=transform,
                    command_runner=command_runner,
                )
            bundle_digest = _publish_bundle(
                bundle_path,
                _bundle_document(
                    contract=contract,
                    contract_sha256=contract_sha256,
                    image_digest=image_digest,
                    retrieved_on=retrieval_date,
                    transforms=transforms,
                    source_manifest_sha256s=source_manifest_sha256s,
                ),
                root,
            )
            return GnomadAllSitesBootstrapResult(
                resource_id=contract.resource_id,
                release=contract.release,
                output_directory=str(release),
                bundle_manifest_path=str(bundle_path),
                bundle_manifest_sha256=bundle_digest,
                transforms=tuple(transforms),
                reused=False,
            )
    except ManifestError as error:
        raise GnomadAllSitesError(f"cannot publish all-sites source manifest: {error}") from error
    finally:
        shutil.rmtree(run_root, ignore_errors=True)
        try:
            run_root.parent.rmdir()
        except OSError:
            pass
        if staging is not None and not published_release:
            shutil.rmtree(staging, ignore_errors=True)


class GnomadAllSitesLookup:
    """Local tabix API for separate all-sites exome and genome callability."""

    def __init__(
        self,
        *,
        contract_path: str | Path,
        reference_root: str | Path,
        command_runner: CommandRunner = _run_command,
    ) -> None:
        self.contract = load_gnomad_all_sites_contract(contract_path)
        self.root = Path(reference_root).resolve(strict=True)
        self.command_runner = command_runner
        self.paths = {
            source.name: self.root / self.contract.output_directory / source.output_name
            for source in self.contract.sources
        }
        for path in self.paths.values():
            _require_within_root(self.root, path, "all-sites lookup path")

    def _dataset_coverage(
        self, source: GnomadAllSitesSource, contig: str, position: int
    ) -> AllSitesDatasetCoverage | None:
        path = self.paths[source.name]
        output = self.command_runner(
            ("tabix", str(path), f"{contig}:{position}-{position}")
        )
        matches = []
        for line in output.splitlines():
            observed_contig, observed_position, coverage = _parse_derived_record(line)
            if observed_contig == contig and observed_position == position:
                matches.append(coverage)
        if len(matches) > 1:
            raise GnomadAllSitesError("all-sites lookup returned duplicate loci")
        return matches[0] if matches else None

    def lookup(self, locus: str) -> AllSitesLocusCoverage:
        contig, position = _parse_locus(locus)
        return AllSitesLocusCoverage(
            locus=f"{contig}:{position}",
            exomes=self._dataset_coverage(self.contract.source("exomes"), contig, position),
            genomes=self._dataset_coverage(self.contract.source("genomes"), contig, position),
        )


def noncontribution_query_basis(text, *, contig, start, end, status):
    """Interpret a successful complete-span query without changing lookup AN.

    The caller must authenticate the complete pinned resource and query scope.
    A failed/unindexed query never proves absence; malformed successful output
    is an integrity failure. Missing published AN remains unknown.
    """
    if status not in ('matched', 'no_match'):
        return None
    if (type(text) is not str or type(start) is not int or type(end) is not int
            or start < 1 or end < start or type(contig) is not str):
        raise GnomadAllSitesError('invalid noncontribution query geometry')
    seen = set()
    unknown = False
    callable_locus = False
    for line in text.splitlines():
        c, position, coverage = _parse_derived_record(line)
        if c != contig or not start <= position <= end or position in seen:
            raise GnomadAllSitesError('invalid noncontribution query coverage')
        seen.add(position)
        unknown |= coverage.an is None
        callable_locus |= coverage.an is not None and coverage.an > 0
    if bool(seen) != (status == 'matched'):
        raise GnomadAllSitesError('noncontribution query status conflicts')
    if unknown or callable_locus:
        return None
    return 'explicit_zero_callable' if seen else 'absent_from_complete_all_sites'


def gnomad_all_sites_result_as_json(result: GnomadAllSitesBootstrapResult) -> str:
    return json.dumps(asdict(result), indent=2, sort_keys=True) + "\n"
