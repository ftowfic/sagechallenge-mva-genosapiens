"""Bootstrap and verify the pinned Ensembl VEP offline reference bundle."""

from __future__ import annotations

from collections.abc import Callable, Mapping, Sequence
from dataclasses import asdict, dataclass
from datetime import date
import gzip
import hashlib
import json
import os
from pathlib import Path, PurePosixPath
import re
import shutil
import subprocess
import tarfile
import tempfile
from typing import Any
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


SCHEMA = "genosapiens.ensembl-vep-bundle"
CACHE_TRANSFORMATION = "genosapiens.ensembl-vep-cache-extraction"
FASTA_TRANSFORMATION = "genosapiens.ensembl-vep-fasta-preparation"
INVENTORY_SCHEMA = "genosapiens.ensembl-vep-cache-inventory"
SHA256 = re.compile(r"^[0-9a-f]{64}$")
DIGEST = re.compile(r"^sha256:[0-9a-f]{64}$")
IDENTIFIER = re.compile(r"^[a-z0-9][a-z0-9._-]*$")
VEP_VERSION = re.compile(r"ensembl-vep\s*:\s*([0-9]+(?:\.[0-9]+)?)")


class VepBundleError(RuntimeError):
    """Raised when VEP acquisition, provenance, or offline validation fails."""


def _require_string(data: Mapping[str, Any], name: str) -> str:
    value = data.get(name)
    if not isinstance(value, str) or not value:
        raise VepBundleError(f"VEP contract field {name} must be a non-empty string")
    return value


def _require_positive_int(data: Mapping[str, Any], name: str) -> int:
    value = data.get(name)
    if not isinstance(value, int) or isinstance(value, bool) or value <= 0:
        raise VepBundleError(f"VEP contract field {name} must be a positive integer")
    return value


def _safe_relative_path(value: str, description: str) -> PurePosixPath:
    path = PurePosixPath(value)
    if path.is_absolute() or path == PurePosixPath(".") or ".." in path.parts:
        raise VepBundleError(f"{description} must be a safe relative path")
    return path


@dataclass(frozen=True, slots=True)
class VepSource:
    resource_id: str
    release: str
    url: str
    filename: str
    size_bytes: int
    sha256: str
    publisher_checksum_url: str
    publisher_bsd_sum: int
    publisher_blocks: int
    license: str
    license_url: str

    @classmethod
    def from_mapping(cls, data: Mapping[str, Any]) -> "VepSource":
        expected = {
            "resource_id",
            "release",
            "url",
            "filename",
            "size_bytes",
            "sha256",
            "publisher_checksum_url",
            "publisher_bsd_sum",
            "publisher_blocks",
            "license",
            "license_url",
        }
        if set(data) != expected:
            raise VepBundleError("VEP source fields are invalid")
        source = cls(
            resource_id=_require_string(data, "resource_id"),
            release=_require_string(data, "release"),
            url=_require_string(data, "url"),
            filename=_require_string(data, "filename"),
            size_bytes=_require_positive_int(data, "size_bytes"),
            sha256=_require_string(data, "sha256"),
            publisher_checksum_url=_require_string(data, "publisher_checksum_url"),
            publisher_bsd_sum=_require_positive_int(data, "publisher_bsd_sum"),
            publisher_blocks=_require_positive_int(data, "publisher_blocks"),
            license=_require_string(data, "license"),
            license_url=_require_string(data, "license_url"),
        )
        if (
            not IDENTIFIER.fullmatch(source.resource_id)
            or not source.url.startswith("https://")
            or not source.publisher_checksum_url.startswith("https://")
            or not source.license_url.startswith("https://")
            or PurePosixPath(source.filename).name != source.filename
            or not SHA256.fullmatch(source.sha256)
        ):
            raise VepBundleError(f"invalid VEP source: {source.resource_id}")
        return source


@dataclass(frozen=True, slots=True)
class VepCacheContract:
    source: VepSource
    output_root: str
    cache_root: str
    info_path: str
    info_sha256: str
    member_count: int
    regular_file_count: int
    directory_count: int
    materialized_directory_count: int
    total_file_bytes: int
    maximum_file_bytes: int

    @classmethod
    def from_mapping(cls, data: Mapping[str, Any]) -> "VepCacheContract":
        expected = {
            "source",
            "output_root",
            "cache_root",
            "info_path",
            "info_sha256",
            "member_count",
            "regular_file_count",
            "directory_count",
            "materialized_directory_count",
            "total_file_bytes",
            "maximum_file_bytes",
        }
        if set(data) != expected or not isinstance(data.get("source"), Mapping):
            raise VepBundleError("VEP cache fields are invalid")
        value = cls(
            source=VepSource.from_mapping(data["source"]),
            output_root=_require_string(data, "output_root"),
            cache_root=_require_string(data, "cache_root"),
            info_path=_require_string(data, "info_path"),
            info_sha256=_require_string(data, "info_sha256"),
            member_count=_require_positive_int(data, "member_count"),
            regular_file_count=_require_positive_int(data, "regular_file_count"),
            directory_count=_require_positive_int(data, "directory_count"),
            materialized_directory_count=_require_positive_int(
                data, "materialized_directory_count"
            ),
            total_file_bytes=_require_positive_int(data, "total_file_bytes"),
            maximum_file_bytes=_require_positive_int(data, "maximum_file_bytes"),
        )
        output_root = _safe_relative_path(value.output_root, "VEP output_root")
        cache_root = _safe_relative_path(value.cache_root, "VEP cache_root")
        info_path = _safe_relative_path(value.info_path, "VEP info_path")
        if (
            not cache_root.is_relative_to(output_root)
            or cache_root == output_root
            or not info_path.is_relative_to(cache_root)
            or info_path == cache_root
            or value.member_count
            != value.regular_file_count + value.directory_count
            or value.materialized_directory_count < value.directory_count
            or not SHA256.fullmatch(value.info_sha256)
        ):
            raise VepBundleError("VEP cache archive contract is inconsistent")
        return value


@dataclass(frozen=True, slots=True)
class VepFastaContract:
    source: VepSource
    output_path: str
    output_size_bytes: int
    output_sha256: str
    fai_size_bytes: int
    fai_sha256: str
    sequence_count: int
    required_sequence_lengths: Mapping[str, int]

    @classmethod
    def from_mapping(cls, data: Mapping[str, Any]) -> "VepFastaContract":
        expected = {
            "source",
            "output_path",
            "output_size_bytes",
            "output_sha256",
            "fai_size_bytes",
            "fai_sha256",
            "sequence_count",
            "required_sequence_lengths",
        }
        lengths = data.get("required_sequence_lengths")
        if (
            set(data) != expected
            or not isinstance(data.get("source"), Mapping)
            or not isinstance(lengths, Mapping)
            or not lengths
            or any(
                not isinstance(name, str)
                or not name
                or not isinstance(length, int)
                or isinstance(length, bool)
                or length <= 0
                for name, length in lengths.items()
            )
        ):
            raise VepBundleError("VEP FASTA fields are invalid")
        value = cls(
            source=VepSource.from_mapping(data["source"]),
            output_path=_require_string(data, "output_path"),
            output_size_bytes=_require_positive_int(data, "output_size_bytes"),
            output_sha256=_require_string(data, "output_sha256"),
            fai_size_bytes=_require_positive_int(data, "fai_size_bytes"),
            fai_sha256=_require_string(data, "fai_sha256"),
            sequence_count=_require_positive_int(data, "sequence_count"),
            required_sequence_lengths=dict(lengths),
        )
        if (
            not SHA256.fullmatch(value.output_sha256)
            or not SHA256.fullmatch(value.fai_sha256)
        ):
            raise VepBundleError("VEP FASTA output checksums are invalid")
        _safe_relative_path(value.output_path, "VEP FASTA output_path")
        return value

    @property
    def fai_path(self) -> str:
        return f"{self.output_path}.fai"


@dataclass(frozen=True, slots=True)
class VepVersions:
    cache: str
    assembly: str
    source_assembly: str
    mane: str
    mane_source_url: str
    genebuild: str
    gencode: str
    clinvar: str
    gnomad_exomes: str
    gnomad_genomes: str

    @classmethod
    def from_mapping(cls, data: Mapping[str, Any]) -> "VepVersions":
        expected = {
            "cache",
            "assembly",
            "source_assembly",
            "mane",
            "mane_source_url",
            "genebuild",
            "gencode",
            "clinvar",
            "gnomad_exomes",
            "gnomad_genomes",
        }
        if set(data) != expected:
            raise VepBundleError("VEP version fields are invalid")
        value = cls(**{name: _require_string(data, name) for name in expected})
        if not value.mane_source_url.startswith("https://"):
            raise VepBundleError("VEP MANE source URL must use HTTPS")
        return value


@dataclass(frozen=True, slots=True)
class VepSentinel:
    chromosome: str
    position: int
    identifier: str
    reference: str
    alternate: str
    gene_symbol: str
    gene_id: str
    transcript_id: str
    consequence: str
    existing_variation: str
    mane_select: str

    @classmethod
    def from_mapping(cls, data: Mapping[str, Any]) -> "VepSentinel":
        expected = {
            "chromosome",
            "position",
            "identifier",
            "reference",
            "alternate",
            "gene_symbol",
            "gene_id",
            "transcript_id",
            "consequence",
            "existing_variation",
            "mane_select",
        }
        if set(data) != expected:
            raise VepBundleError("VEP sentinel fields are invalid")
        position = _require_positive_int(data, "position")
        values = {name: _require_string(data, name) for name in expected - {"position"}}
        value = cls(position=position, **values)
        if (
            not re.fullmatch(r"[ACGTN]+", value.reference)
            or not re.fullmatch(r"[ACGTN]+", value.alternate)
            or value.reference == value.alternate
            or not re.fullmatch(r"ENSG[0-9]{11}", value.gene_id)
            or not re.fullmatch(r"ENST[0-9]{11}", value.transcript_id)
        ):
            raise VepBundleError("VEP sentinel identity is invalid")
        return value


@dataclass(frozen=True, slots=True)
class VepBundleContract:
    bundle_id: str
    tool_version: str
    tool_image: str
    tool_image_digest: str
    cache: VepCacheContract
    fasta: VepFastaContract
    versions: VepVersions
    sentinel: VepSentinel
    schema: str = SCHEMA

    @classmethod
    def from_mapping(cls, data: Mapping[str, Any]) -> "VepBundleContract":
        expected = {
            "schema",
            "bundle_id",
            "tool_version",
            "tool_image",
            "tool_image_digest",
            "cache",
            "fasta",
            "versions",
            "sentinel",
        }
        nested = (data.get("cache"), data.get("fasta"), data.get("versions"), data.get("sentinel"))
        if (
            set(data) != expected
            or data.get("schema") != SCHEMA
            or not all(isinstance(value, Mapping) for value in nested)
        ):
            raise VepBundleError("VEP bundle fields or schema are invalid")
        value = cls(
            schema=SCHEMA,
            bundle_id=_require_string(data, "bundle_id"),
            tool_version=_require_string(data, "tool_version"),
            tool_image=_require_string(data, "tool_image"),
            tool_image_digest=_require_string(data, "tool_image_digest"),
            cache=VepCacheContract.from_mapping(data["cache"]),
            fasta=VepFastaContract.from_mapping(data["fasta"]),
            versions=VepVersions.from_mapping(data["versions"]),
            sentinel=VepSentinel.from_mapping(data["sentinel"]),
        )
        output_root = _safe_relative_path(value.cache.output_root, "VEP output_root")
        if (
            not IDENTIFIER.fullmatch(value.bundle_id)
            or not re.fullmatch(r"[0-9]+\.[0-9]+", value.tool_version)
            or not value.tool_image.startswith("ensemblorg/ensembl-vep:")
            or not DIGEST.fullmatch(value.tool_image_digest)
            or not _safe_relative_path(
                value.fasta.output_path, "VEP FASTA output_path"
            ).is_relative_to(output_root)
            or value.cache.source.release != value.versions.cache
            or value.fasta.source.release != value.versions.cache
        ):
            raise VepBundleError("VEP bundle identities are inconsistent")
        return value


def load_vep_bundle_contract(path: str | Path) -> VepBundleContract:
    contract_path = Path(path)
    try:
        data = json.loads(contract_path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as error:
        raise VepBundleError(f"cannot read VEP contract {contract_path}: {error}") from error
    if not isinstance(data, Mapping):
        raise VepBundleError("VEP contract root must be an object")
    return VepBundleContract.from_mapping(data)


@dataclass(frozen=True, slots=True)
class VepBundleValidation:
    cache_members: int
    cache_files: int
    cache_directories: int
    cache_file_bytes: int
    cache_inventory_sha256: str
    fasta_sequences: int
    cache_version: str
    assembly: str
    mane_version: str
    genebuild: str
    gencode_version: str
    clinvar_version: str
    gnomad_exomes_version: str
    gnomad_genomes_version: str


@dataclass(frozen=True, slots=True)
class VepSentinelResult:
    tool_version: str
    variant: str
    gene_symbol: str
    gene_id: str
    transcript_id: str
    consequence: str
    existing_variation: str
    mane_select: str
    gnomad_exomes_af: str
    gnomad_genomes_af: str
    clinical_significance: str


@dataclass(frozen=True, slots=True)
class VepBootstrapResult:
    bundle_id: str
    bundle_manifest_path: str
    bundle_manifest_sha256: str
    cache_manifest_path: str
    cache_manifest_sha256: str
    fasta_manifest_path: str
    fasta_manifest_sha256: str
    validation: VepBundleValidation
    sentinel: VepSentinelResult
    reused: bool


Downloader = Callable[[DownloadRequest, str | Path], DownloadResult]
CapacityChecker = Callable[[str | Path, CapacityPlan], Any]
SentinelRunner = Callable[[VepBundleContract, Path, Path], VepSentinelResult]


def _file_identity(path: Path) -> tuple[str, int]:
    digest = hashlib.sha256()
    size = 0
    try:
        with path.open("rb") as source:
            while chunk := source.read(1024 * 1024):
                digest.update(chunk)
                size += len(chunk)
    except OSError as error:
        raise VepBundleError(f"cannot read VEP artifact {path}: {error}") from error
    return digest.hexdigest(), size


def _require_within_root(root: Path, path: Path, description: str) -> None:
    if not path.is_absolute() or not path.is_relative_to(root):
        raise VepBundleError(f"{description} escapes reference_root: {path}")
    relative = path.relative_to(root)
    current = root
    for part in relative.parts:
        current /= part
        if current.is_symlink():
            raise VepBundleError(f"{description} contains a symlink: {current}")
    try:
        resolved = path.resolve()
    except OSError as error:
        raise VepBundleError(f"cannot resolve {description} {path}: {error}") from error
    if not resolved.is_relative_to(root):
        raise VepBundleError(f"{description} escapes reference_root: {path}")


def _canonical_json(data: Mapping[str, Any]) -> bytes:
    return (json.dumps(data, sort_keys=True, separators=(",", ":")) + "\n").encode()


def _write_new_file(path: Path, content: bytes) -> tuple[str, int]:
    try:
        with path.open("xb") as output:
            output.write(content)
            output.flush()
            os.fsync(output.fileno())
        os.chmod(path, 0o444)
    except OSError as error:
        raise VepBundleError(f"cannot write VEP artifact {path}: {error}") from error
    return hashlib.sha256(content).hexdigest(), len(content)


def _parse_info(path: Path, contract: VepBundleContract) -> Mapping[str, str]:
    digest, _size = _file_identity(path)
    if digest != contract.cache.info_sha256:
        raise VepBundleError("VEP cache info.txt checksum mismatch")
    metadata: dict[str, str] = {}
    try:
        for raw_line in path.read_text(encoding="utf-8").splitlines():
            fields = raw_line.split("\t")
            if len(fields) != 2 or not all(fields):
                raise VepBundleError("VEP cache info.txt row is invalid")
            key, value = fields
            if key in metadata:
                raise VepBundleError(f"duplicate VEP cache metadata field: {key}")
            metadata[key] = value
    except (OSError, UnicodeDecodeError) as error:
        raise VepBundleError(f"cannot parse VEP cache info.txt: {error}") from error
    versions = contract.versions
    expected = {
        "species": "homo_sapiens",
        "assembly": versions.assembly,
        "source_assembly": versions.source_assembly,
        "source_genebuild": versions.genebuild,
        "source_gencode": versions.gencode,
        "source_ClinVar": versions.clinvar,
        "source_gnomADe": versions.gnomad_exomes,
        "source_gnomADg": versions.gnomad_genomes,
        "var_type": "tabix",
    }
    for key, value in expected.items():
        if metadata.get(key) != value:
            raise VepBundleError(f"VEP cache metadata mismatch for {key}")
    return metadata


def _extract_cache(
    archive: Path, staging: Path, contract: VepBundleContract
) -> tuple[str, int]:
    cache = contract.cache
    expected_root = PurePosixPath(cache.cache_root).relative_to(cache.output_root)
    seen: set[str] = set()
    entries: list[dict[str, Any]] = []
    members = files = directories = total_bytes = 0
    try:
        with tarfile.open(archive, mode="r|gz") as source:
            for member in source:
                members += 1
                path = PurePosixPath(member.name)
                archive_top = PurePosixPath(expected_root.parts[0])
                if (
                    path.is_absolute()
                    or path == PurePosixPath(".")
                    or ".." in path.parts
                    or (path != archive_top and not path.is_relative_to(expected_root))
                ):
                    raise VepBundleError(f"unsafe VEP cache archive path: {member.name}")
                normalized = str(path)
                if normalized in seen:
                    raise VepBundleError(f"duplicate VEP cache archive path: {member.name}")
                seen.add(normalized)
                destination = staging / normalized
                if member.isdir():
                    directories += 1
                    destination.mkdir(parents=True, exist_ok=True, mode=0o755)
                    continue
                if not member.isreg():
                    raise VepBundleError(f"unsupported VEP cache archive member: {member.name}")
                if member.size < 0 or member.size > cache.maximum_file_bytes:
                    raise VepBundleError(f"oversized VEP cache archive member: {member.name}")
                files += 1
                total_bytes += member.size
                if total_bytes > cache.total_file_bytes:
                    raise VepBundleError("VEP cache archive exceeds the byte contract")
                destination.parent.mkdir(parents=True, exist_ok=True, mode=0o755)
                extracted = source.extractfile(member)
                if extracted is None:
                    raise VepBundleError(f"cannot read VEP cache member: {member.name}")
                digest = hashlib.sha256()
                written = 0
                with destination.open("xb") as output:
                    while chunk := extracted.read(1024 * 1024):
                        output.write(chunk)
                        digest.update(chunk)
                        written += len(chunk)
                    output.flush()
                    os.fsync(output.fileno())
                if written != member.size:
                    raise VepBundleError(f"truncated VEP cache member: {member.name}")
                os.chmod(destination, 0o444)
                entries.append(
                    {"path": normalized, "size_bytes": written, "sha256": digest.hexdigest()}
                )
    except (OSError, tarfile.TarError) as error:
        raise VepBundleError(f"cannot extract VEP cache: {error}") from error
    if (
        members != cache.member_count
        or files != cache.regular_file_count
        or directories != cache.directory_count
        or total_bytes != cache.total_file_bytes
    ):
        raise VepBundleError("VEP cache archive member totals do not match the contract")
    _parse_info(staging / str(PurePosixPath(cache.info_path).relative_to(cache.output_root)), contract)
    inventory = {
        "schema": INVENTORY_SCHEMA,
        "cache_root": str(expected_root),
        "entries": sorted(entries, key=lambda item: item["path"]),
    }
    return _write_new_file(staging / "cache.inventory.json", _canonical_json(inventory))


def _prepare_fasta(
    compressed: Path, staging: Path, contract: VepBundleContract
) -> None:
    fasta = contract.fasta
    relative = PurePosixPath(fasta.output_path).relative_to(contract.cache.output_root)
    destination = staging / str(relative)
    destination.parent.mkdir(parents=True, exist_ok=True, mode=0o755)
    digest = hashlib.sha256()
    entries: list[tuple[str, int, int, int, int]] = []
    names: set[str] = set()
    offset = 0
    name: str | None = None
    sequence_offset = length = line_bases = line_width = 0
    try:
        with gzip.open(compressed, "rb") as source, destination.open("xb") as output:
            for line in source:
                output.write(line)
                digest.update(line)
                if line.startswith(b">"):
                    if name is not None:
                        entries.append((name, length, sequence_offset, line_bases, line_width))
                    name = line[1:].split(None, 1)[0].decode("ascii")
                    if not name or name in names:
                        raise VepBundleError("VEP FASTA contains a duplicate sequence name")
                    names.add(name)
                    sequence_offset = offset + len(line)
                    length = line_bases = line_width = 0
                else:
                    if name is None:
                        raise VepBundleError("VEP FASTA sequence precedes its header")
                    sequence = line.rstrip(b"\r\n")
                    if not sequence or not re.fullmatch(rb"[ACGTNacgtn]+", sequence):
                        raise VepBundleError(f"VEP FASTA sequence is invalid: {name}")
                    if line_bases == 0:
                        line_bases = len(sequence)
                        line_width = len(line)
                    elif len(sequence) > line_bases:
                        raise VepBundleError(f"VEP FASTA line width is invalid: {name}")
                    length += len(sequence)
                offset += len(line)
            if name is not None:
                entries.append((name, length, sequence_offset, line_bases, line_width))
            output.flush()
            os.fsync(output.fileno())
    except (OSError, EOFError, UnicodeDecodeError) as error:
        raise VepBundleError(f"cannot prepare VEP FASTA: {error}") from error
    if (
        offset != fasta.output_size_bytes
        or digest.hexdigest() != fasta.output_sha256
        or len(entries) != fasta.sequence_count
    ):
        raise VepBundleError("VEP FASTA identity does not match the contract")
    lengths = {entry[0]: entry[1] for entry in entries}
    for required, expected_length in fasta.required_sequence_lengths.items():
        if lengths.get(required) != expected_length:
            raise VepBundleError(f"VEP FASTA sequence length mismatch: {required}")
    os.chmod(destination, 0o444)
    fai = b"".join(
        f"{entry_name}\t{entry_length}\t{entry_offset}\t{bases}\t{width}\n".encode()
        for entry_name, entry_length, entry_offset, bases, width in entries
    )
    if len(fai) != fasta.fai_size_bytes or hashlib.sha256(fai).hexdigest() != fasta.fai_sha256:
        raise VepBundleError("generated VEP FASTA index identity does not match the contract")
    _write_new_file(Path(f"{destination}.fai"), fai)


def _load_inventory(path: Path, contract: VepBundleContract) -> list[Mapping[str, Any]]:
    try:
        data = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as error:
        raise VepBundleError(f"cannot read VEP cache inventory: {error}") from error
    expected_root = str(
        PurePosixPath(contract.cache.cache_root).relative_to(contract.cache.output_root)
    )
    if (
        not isinstance(data, Mapping)
        or set(data) != {"schema", "cache_root", "entries"}
        or data.get("schema") not in {INVENTORY_SCHEMA, f"{INVENTORY_SCHEMA}/v1"}
        or data.get("cache_root") != expected_root
        or not isinstance(data.get("entries"), list)
    ):
        raise VepBundleError("VEP cache inventory structure is invalid")
    entries = data["entries"]
    if len(entries) != contract.cache.regular_file_count:
        raise VepBundleError("VEP cache inventory file count mismatch")
    return entries


def _validate_release(
    release: Path, contract: VepBundleContract
) -> VepBundleValidation:
    if not release.is_dir() or release.is_symlink():
        raise VepBundleError(f"VEP release is not a directory: {release}")
    inventory_path = release / "cache.inventory.json"
    entries = _load_inventory(inventory_path, contract)
    expected_paths: set[str] = {"cache.inventory.json"}
    total_bytes = 0
    seen: set[str] = set()
    for entry in entries:
        if not isinstance(entry, Mapping) or set(entry) != {"path", "size_bytes", "sha256"}:
            raise VepBundleError("VEP cache inventory entry is invalid")
        path_value = entry.get("path")
        size = entry.get("size_bytes")
        digest = entry.get("sha256")
        if (
            not isinstance(path_value, str)
            or str(_safe_relative_path(path_value, "VEP inventory path")) in seen
            or not isinstance(size, int)
            or isinstance(size, bool)
            or size < 0
            or not isinstance(digest, str)
            or not SHA256.fullmatch(digest)
        ):
            raise VepBundleError("VEP cache inventory entry identity is invalid")
        normalized = str(_safe_relative_path(path_value, "VEP inventory path"))
        cache_root = PurePosixPath(contract.cache.cache_root).relative_to(
            contract.cache.output_root
        )
        if not PurePosixPath(normalized).is_relative_to(cache_root):
            raise VepBundleError("VEP inventory path escapes the cache root")
        seen.add(normalized)
        expected_paths.add(normalized)
        path = release / normalized
        if path.is_symlink() or not path.is_file():
            raise VepBundleError(f"VEP cache output is missing or invalid: {normalized}")
        actual_digest, actual_size = _file_identity(path)
        if actual_size != size or actual_digest != digest:
            raise VepBundleError(f"VEP cache output checksum mismatch: {normalized}")
        total_bytes += size
    if total_bytes != contract.cache.total_file_bytes:
        raise VepBundleError("VEP cache inventory byte count mismatch")
    expected_paths.update(
        {
            str(PurePosixPath(contract.fasta.output_path).relative_to(contract.cache.output_root)),
            str(PurePosixPath(contract.fasta.fai_path).relative_to(contract.cache.output_root)),
        }
    )
    actual_paths: set[str] = set()
    directory_count = 0
    cache_top = PurePosixPath(
        PurePosixPath(contract.cache.cache_root)
        .relative_to(contract.cache.output_root)
        .parts[0]
    )
    for path in release.rglob("*"):
        if path.is_symlink():
            raise VepBundleError(f"VEP release contains a symlink: {path}")
        relative = str(path.relative_to(release))
        if path.is_file():
            actual_paths.add(relative)
        elif path.is_dir():
            if PurePosixPath(relative).is_relative_to(cache_top):
                directory_count += 1
        else:
            raise VepBundleError(f"VEP release contains an invalid entry: {path}")
    if (
        actual_paths != expected_paths
        or directory_count != contract.cache.materialized_directory_count
    ):
        missing = sorted(expected_paths - actual_paths)[:5]
        unexpected = sorted(actual_paths - expected_paths)[:5]
        raise VepBundleError(
            "VEP release directory is incomplete or contaminated: "
            f"files={len(actual_paths)}/{len(expected_paths)}, "
            "cache_directories="
            f"{directory_count}/{contract.cache.materialized_directory_count}, "
            f"missing={missing}, unexpected={unexpected}"
        )
    fasta_path = release / str(
        PurePosixPath(contract.fasta.output_path).relative_to(contract.cache.output_root)
    )
    fasta_identity = _file_identity(fasta_path)
    fai_identity = _file_identity(Path(f"{fasta_path}.fai"))
    if fasta_identity != (contract.fasta.output_sha256, contract.fasta.output_size_bytes):
        raise VepBundleError("VEP FASTA checksum mismatch")
    if fai_identity != (contract.fasta.fai_sha256, contract.fasta.fai_size_bytes):
        raise VepBundleError("VEP FASTA index checksum mismatch")
    _parse_info(
        release
        / str(PurePosixPath(contract.cache.info_path).relative_to(contract.cache.output_root)),
        contract,
    )
    inventory_sha256 = _file_identity(inventory_path)[0]
    versions = contract.versions
    return VepBundleValidation(
        cache_members=contract.cache.member_count,
        cache_files=contract.cache.regular_file_count,
        cache_directories=contract.cache.materialized_directory_count,
        cache_file_bytes=total_bytes,
        cache_inventory_sha256=inventory_sha256,
        fasta_sequences=contract.fasta.sequence_count,
        cache_version=versions.cache,
        assembly=versions.source_assembly,
        mane_version=versions.mane,
        genebuild=versions.genebuild,
        gencode_version=versions.gencode,
        clinvar_version=versions.clinvar,
        gnomad_exomes_version=versions.gnomad_exomes,
        gnomad_genomes_version=versions.gnomad_genomes,
    )


def _parse_sentinel_vcf(
    vcf: str, tool_output: str, contract: VepBundleContract
) -> VepSentinelResult:
    version_match = VEP_VERSION.search(tool_output)
    if version_match is None or version_match.group(1) != contract.tool_version:
        raise VepBundleError("VEP executable version does not match the contract")
    format_fields: list[str] | None = None
    records: list[list[str]] = []
    for line in vcf.splitlines():
        if line.startswith("##INFO=<ID=CSQ"):
            match = re.search(r"Format: ([^\"]+)", line)
            if match is not None:
                format_fields = match.group(1).rstrip(">").split("|")
        elif line and not line.startswith("#"):
            records.append(line.split("\t"))
    if format_fields is None or len(records) != 1 or len(records[0]) < 8:
        raise VepBundleError("VEP sentinel output is missing a valid CSQ record")
    sentinel = contract.sentinel
    record = records[0]
    if record[:5] != [
        sentinel.chromosome,
        str(sentinel.position),
        sentinel.identifier,
        sentinel.reference,
        sentinel.alternate,
    ]:
        raise VepBundleError("VEP sentinel output variant identity changed")
    info = {
        field.split("=", 1)[0]: field.split("=", 1)[1]
        for field in record[7].split(";")
        if "=" in field
    }
    annotations = []
    for raw in info.get("CSQ", "").split(","):
        values = raw.split("|")
        values.extend([""] * (len(format_fields) - len(values)))
        annotations.append(dict(zip(format_fields, values, strict=False)))
    matches = [
        annotation
        for annotation in annotations
        if annotation.get("SYMBOL") == sentinel.gene_symbol
        and annotation.get("Gene") == sentinel.gene_id
        and annotation.get("Feature") == sentinel.transcript_id
        and sentinel.consequence in annotation.get("Consequence", "").split("&")
    ]
    if len(matches) != 1:
        raise VepBundleError("VEP sentinel transcript annotation is missing or ambiguous")
    annotation = matches[0]
    if (
        sentinel.existing_variation
        not in annotation.get("Existing_variation", "").split("&")
        or annotation.get("MANE_SELECT") != sentinel.mane_select
        or not annotation.get("gnomADe_AF")
        or not annotation.get("gnomADg_AF")
        or not annotation.get("CLIN_SIG")
    ):
        raise VepBundleError("VEP sentinel source annotations are incomplete")
    return VepSentinelResult(
        tool_version=contract.tool_version,
        variant=(
            f"{sentinel.chromosome}:{sentinel.position}:"
            f"{sentinel.reference}:{sentinel.alternate}"
        ),
        gene_symbol=annotation["SYMBOL"],
        gene_id=annotation["Gene"],
        transcript_id=annotation["Feature"],
        consequence=sentinel.consequence,
        existing_variation=sentinel.existing_variation,
        mane_select=annotation["MANE_SELECT"],
        gnomad_exomes_af=annotation["gnomADe_AF"],
        gnomad_genomes_af=annotation["gnomADg_AF"],
        clinical_significance=annotation["CLIN_SIG"],
    )


def run_vep_sentinel(
    contract: VepBundleContract, cache_parent: Path, fasta_path: Path
) -> VepSentinelResult:
    executable = "/opt/vep/src/ensembl-vep/vep"
    sentinel = contract.sentinel
    with tempfile.TemporaryDirectory(prefix="genosapiens-vep-") as directory:
        work = Path(directory)
        input_path = work / "sentinel.vcf"
        output_path = work / "annotated.vcf"
        input_path.write_text(
            "##fileformat=VCFv4.2\n"
            f"##contig=<ID={sentinel.chromosome}>\n"
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
            f"{sentinel.chromosome}\t{sentinel.position}\t{sentinel.identifier}\t"
            f"{sentinel.reference}\t{sentinel.alternate}\t.\tPASS\t.\n",
            encoding="ascii",
        )
        try:
            version = subprocess.run(
                [executable, "--help"],
                check=True,
                capture_output=True,
                text=True,
                timeout=120,
            )
            annotation = subprocess.run(
                [
                    executable,
                    "--cache",
                    "--offline",
                    "--dir_cache",
                    str(cache_parent),
                    "--fasta",
                    str(fasta_path),
                    "--species",
                    "homo_sapiens",
                    "--assembly",
                    contract.versions.assembly,
                    "--input_file",
                    str(input_path),
                    "--output_file",
                    str(output_path),
                    "--format",
                    "vcf",
                    "--vcf",
                    "--force_overwrite",
                    "--no_stats",
                    "--symbol",
                    "--mane",
                    "--check_existing",
                    "--af_gnomade",
                    "--af_gnomadg",
                ],
                check=True,
                capture_output=True,
                text=True,
                timeout=600,
            )
            output = output_path.read_text(encoding="utf-8")
        except (OSError, subprocess.SubprocessError, UnicodeDecodeError) as error:
            raise VepBundleError(f"VEP offline sentinel annotation failed: {error}") from error
        return _parse_sentinel_vcf(
            output,
            version.stdout + version.stderr + annotation.stdout + annotation.stderr,
            contract,
        )


def _manifest_path(
    root: Path, resource_id: str, contract_sha256: str | None = None
) -> Path:
    legacy = root / f"manifests/ensembl/vep/{resource_id}.json"
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
    return legacy.with_name(f"{resource_id}.{contract_sha256[:12]}.json")


def _bundle_manifest_path(
    root: Path, contract: VepBundleContract, contract_sha256: str | None = None
) -> Path:
    legacy = root / f"manifests/ensembl/vep/{contract.bundle_id}.bundle.json"
    if contract_sha256 is None or not legacy.exists():
        return legacy
    try:
        document = json.loads(legacy.read_text(encoding="utf-8"))
        if (
            document.get("schema") == "genosapiens.ensembl-vep-bundle-manifest"
            and document.get("contract_sha256") == contract_sha256
        ):
            return legacy
    except (OSError, json.JSONDecodeError, AttributeError):
        pass
    return legacy.with_name(
        f"{contract.bundle_id}.{contract_sha256[:12]}.bundle.json"
    )


def _load_manifest(path: Path) -> PublicResourceManifest:
    try:
        data = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as error:
        raise VepBundleError(f"cannot read VEP resource manifest {path}: {error}") from error
    if not isinstance(data, Mapping):
        raise VepBundleError("VEP resource manifest root must be an object")
    try:
        return PublicResourceManifest.from_mapping(data)
    except ManifestError as error:
        raise VepBundleError(f"invalid VEP resource manifest: {error}") from error


def _verify_resource_manifests(
    root: Path, contract: VepBundleContract, contract_sha256: str
) -> tuple[tuple[Path, str] | None, tuple[Path, str] | None]:
    cache_path = _manifest_path(
        root, contract.cache.source.resource_id, contract_sha256
    )
    fasta_path = _manifest_path(
        root, contract.fasta.source.resource_id, contract_sha256
    )
    cache_source = contract.cache.source
    fasta_source = contract.fasta.source
    expected_cache_paths = {
        f"{contract.cache.output_root}/cache.inventory.json",
        contract.cache.info_path,
    }
    cache_publication: tuple[Path, str] | None = None
    fasta_publication: tuple[Path, str] | None = None
    if cache_path.exists() or cache_path.is_symlink():
        cache_manifest = _load_manifest(cache_path)
        if (
            cache_manifest.resource_id != cache_source.resource_id
            or cache_manifest.release != cache_source.release
            or cache_manifest.source
            != SourceIdentity(cache_source.url, sha256=cache_source.sha256)
            or cache_manifest.license != cache_source.license
            or cache_manifest.assembly != contract.versions.source_assembly
            or cache_manifest.transformation
            != Transformation(CACHE_TRANSFORMATION, contract_sha256)
            or {output.path for output in cache_manifest.outputs}
            != expected_cache_paths
        ):
            raise VepBundleError(
                "VEP cache manifest identity differs from the contract"
            )
        try:
            verify_outputs(cache_manifest, root)
        except ManifestError as error:
            raise VepBundleError(
                f"completed VEP cache output is invalid: {error}"
            ) from error
        cache_publication = (cache_path, _file_identity(cache_path)[0])
    if fasta_path.exists() or fasta_path.is_symlink():
        fasta_manifest = _load_manifest(fasta_path)
        if (
            fasta_manifest.resource_id != fasta_source.resource_id
            or fasta_manifest.release != fasta_source.release
            or fasta_manifest.source
            != SourceIdentity(fasta_source.url, sha256=fasta_source.sha256)
            or fasta_manifest.license != fasta_source.license
            or fasta_manifest.assembly != contract.versions.source_assembly
            or fasta_manifest.transformation
            != Transformation(FASTA_TRANSFORMATION, contract_sha256)
            or {output.path for output in fasta_manifest.outputs}
            != {contract.fasta.output_path, contract.fasta.fai_path}
        ):
            raise VepBundleError(
                "VEP FASTA manifest identity differs from the contract"
            )
        try:
            verify_outputs(fasta_manifest, root)
        except ManifestError as error:
            raise VepBundleError(
                f"completed VEP FASTA output is invalid: {error}"
            ) from error
        fasta_publication = (fasta_path, _file_identity(fasta_path)[0])
    return (
        cache_publication,
        fasta_publication,
    )


def _publish_json(data: Mapping[str, Any], path: Path, root: Path) -> str:
    _require_within_root(root, path, "VEP bundle manifest")
    content = _canonical_json(data)
    digest = hashlib.sha256(content).hexdigest()
    if path.exists() or path.is_symlink():
        if path.is_symlink() or not path.is_file():
            raise VepBundleError(f"VEP bundle manifest is not a regular file: {path}")
        if path.read_bytes() != content:
            raise VepBundleError("completed VEP bundle manifest is immutable")
        return digest
    path.parent.mkdir(parents=True, exist_ok=True)
    partial_root = path.parent / ".partial"
    _require_within_root(root, partial_root, "VEP bundle manifest partial path")
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
                raise VepBundleError("concurrent VEP bundle manifest differs")
        directory_fd = os.open(path.parent, os.O_RDONLY)
        try:
            os.fsync(directory_fd)
        finally:
            os.close(directory_fd)
    finally:
        partial.unlink(missing_ok=True)
    return digest


def _bundle_data(
    contract: VepBundleContract,
    contract_sha256: str,
    publications: tuple[tuple[Path, str], tuple[Path, str]],
    validation: VepBundleValidation,
    sentinel: VepSentinelResult,
    root: Path,
) -> Mapping[str, Any]:
    return {
        "schema": "genosapiens.ensembl-vep-bundle-manifest",
        "classification": "public",
        "bundle_id": contract.bundle_id,
        "contract_sha256": contract_sha256,
        "tool": {
            "name": "ensembl-vep",
            "version": contract.tool_version,
            "upstream_image": contract.tool_image,
            "upstream_image_digest": contract.tool_image_digest,
        },
        "resources": [
            {
                "resource_id": source.resource_id,
                "release": source.release,
                "source_url": source.url,
                "source_sha256": source.sha256,
                "publisher_checksum_url": source.publisher_checksum_url,
                "publisher_bsd_sum": source.publisher_bsd_sum,
                "publisher_blocks": source.publisher_blocks,
                "license": source.license,
                "license_url": source.license_url,
                "manifest_path": str(publication[0].relative_to(root)),
                "manifest_sha256": publication[1],
            }
            for source, publication in zip(
                (contract.cache.source, contract.fasta.source),
                publications,
                strict=True,
            )
        ],
        "versions": asdict(contract.versions),
        "validation": asdict(validation),
        "offline_sentinel": asdict(sentinel),
    }


def _publish_resource_manifests(
    root: Path,
    release: Path,
    contract: VepBundleContract,
    contract_sha256: str,
    image_digest: str,
    retrieved_on: str,
    existing: tuple[tuple[Path, str] | None, tuple[Path, str] | None] = (
        None,
        None,
    ),
) -> tuple[tuple[Path, str], tuple[Path, str]]:
    info_path = root / contract.cache.info_path
    inventory_path = release / "cache.inventory.json"
    inventory_sha, inventory_size = _file_identity(inventory_path)
    info_sha, info_size = _file_identity(info_path)
    cache_source = contract.cache.source
    cache_manifest = PublicResourceManifest(
        resource_id=cache_source.resource_id,
        release=cache_source.release,
        source=SourceIdentity(cache_source.url, sha256=cache_source.sha256),
        license=cache_source.license,
        assembly=contract.versions.source_assembly,
        retrieved_on=retrieved_on,
        transformation=Transformation(CACHE_TRANSFORMATION, contract_sha256),
        image_digest=image_digest,
        outputs=(
            ResourceOutput(
                path=f"{contract.cache.output_root}/cache.inventory.json",
                size_bytes=inventory_size,
                sha256=inventory_sha,
                retention=RetentionPolicy("durable", None),
            ),
            ResourceOutput(
                path=contract.cache.info_path,
                size_bytes=info_size,
                sha256=info_sha,
                retention=RetentionPolicy("durable", None),
            ),
        ),
    )
    fasta_source = contract.fasta.source
    fasta_manifest = PublicResourceManifest(
        resource_id=fasta_source.resource_id,
        release=fasta_source.release,
        source=SourceIdentity(fasta_source.url, sha256=fasta_source.sha256),
        license=fasta_source.license,
        assembly=contract.versions.source_assembly,
        retrieved_on=retrieved_on,
        transformation=Transformation(FASTA_TRANSFORMATION, contract_sha256),
        image_digest=image_digest,
        outputs=(
            ResourceOutput(
                path=contract.fasta.output_path,
                size_bytes=contract.fasta.output_size_bytes,
                sha256=contract.fasta.output_sha256,
                retention=RetentionPolicy("durable", None),
            ),
            ResourceOutput(
                path=contract.fasta.fai_path,
                size_bytes=contract.fasta.fai_size_bytes,
                sha256=contract.fasta.fai_sha256,
                retention=RetentionPolicy("durable", None),
            ),
        ),
    )
    cache_publication = existing[0]
    if cache_publication is None:
        publication = publish_manifest(
            cache_manifest,
            _manifest_path(root, cache_source.resource_id, contract_sha256),
            root,
        )
        cache_publication = (publication.path, publication.manifest_sha256)
    fasta_publication = existing[1]
    if fasta_publication is None:
        publication = publish_manifest(
            fasta_manifest,
            _manifest_path(root, fasta_source.resource_id, contract_sha256),
            root,
        )
        fasta_publication = (publication.path, publication.manifest_sha256)
    return (
        cache_publication,
        fasta_publication,
    )


def bootstrap_vep_bundle(
    *,
    contract_path: str | Path,
    reference_root: str | Path,
    image_digest: str,
    retrieved_on: str | None = None,
    scratch_root: str | Path | None = None,
    reserve_bytes: int = DEFAULT_RESERVE_BYTES,
    downloader: Downloader = download_atomic,
    capacity_checker: CapacityChecker = require_capacity,
    sentinel_runner: SentinelRunner = run_vep_sentinel,
) -> VepBootstrapResult:
    contract_file = Path(contract_path)
    root = Path(reference_root).resolve(strict=True)
    contract = load_vep_bundle_contract(contract_file)
    contract_sha256 = _file_identity(contract_file)[0]
    if not DIGEST.fullmatch(image_digest):
        raise VepBundleError("VEP bootstrap image_digest must be a SHA-256 digest")
    scratch_base = Path(scratch_root) if scratch_root is not None else root / "scratch/public-build"
    if not scratch_base.is_absolute():
        raise VepBundleError("VEP scratch_root must be absolute")
    _require_within_root(root, scratch_base, "VEP scratch_root")
    release = root / contract.cache.output_root
    _require_within_root(root, release, "VEP release")
    bundle_path = _bundle_manifest_path(root, contract, contract_sha256)
    _require_within_root(root, bundle_path, "VEP bundle manifest")
    bundle_existed = bundle_path.exists()
    if bundle_existed and (bundle_path.is_symlink() or not bundle_path.is_file()):
        raise VepBundleError(f"VEP bundle manifest is not a regular file: {bundle_path}")
    publications = _verify_resource_manifests(root, contract, contract_sha256)
    if release.exists():
        validation = _validate_release(release, contract)
        if bundle_existed and any(publication is None for publication in publications):
            raise VepBundleError("completed VEP bundle is missing resource manifests")
        cache_parent = release / str(
            PurePosixPath(contract.cache.cache_root)
            .relative_to(contract.cache.output_root)
            .parent
            .parent
        )
        fasta_path = root / contract.fasta.output_path
        sentinel = sentinel_runner(contract, cache_parent, fasta_path)
        if any(publication is None for publication in publications):
            completed_publications = _publish_resource_manifests(
                root,
                release,
                contract,
                contract_sha256,
                image_digest,
                retrieved_on or date.today().isoformat(),
                existing=publications,
            )
        else:
            completed_publications = (publications[0], publications[1])
        bundle_digest = _publish_json(
            _bundle_data(
                contract,
                contract_sha256,
                completed_publications,
                validation,
                sentinel,
                root,
            ),
            bundle_path,
            root,
        )
        return VepBootstrapResult(
            bundle_id=contract.bundle_id,
            bundle_manifest_path=str(bundle_path),
            bundle_manifest_sha256=bundle_digest,
            cache_manifest_path=str(completed_publications[0][0]),
            cache_manifest_sha256=completed_publications[0][1],
            fasta_manifest_path=str(completed_publications[1][0]),
            fasta_manifest_sha256=completed_publications[1][1],
            validation=validation,
            sentinel=sentinel,
            reused=bundle_existed,
        )
    if any(publication is not None for publication in publications) or bundle_existed:
        raise VepBundleError("completed VEP manifests exist without the release")

    scratch = scratch_base / "vep-discovery"
    _require_within_root(root, scratch, "VEP build scratch")
    scratch.mkdir(parents=True, exist_ok=True, mode=0o700)
    source_paths = {
        "cache": scratch / contract.cache.source.filename,
        "fasta": scratch / contract.fasta.source.filename,
    }
    for path in source_paths.values():
        _require_within_root(root, path, "VEP source scratch")
    missing_sources = [
        source
        for source, path in (
            (contract.cache.source, source_paths["cache"]),
            (contract.fasta.source, source_paths["fasta"]),
        )
        if not path.exists()
    ]
    capacity_checker(
        root,
        CapacityPlan(
            operation="bootstrap Ensembl VEP 116 offline references",
            source_bytes=sum(source.size_bytes for source in missing_sources),
            output_bytes=(
                contract.cache.total_file_bytes
                + contract.fasta.output_size_bytes
                + contract.fasta.fai_size_bytes
                + 10_000_000
            ),
            reserve_bytes=reserve_bytes,
        ),
    )
    for source, destination in (
        (contract.cache.source, source_paths["cache"]),
        (contract.fasta.source, source_paths["fasta"]),
    ):
        downloader(
            DownloadRequest(
                url=source.url,
                destination=destination,
                expected_size_bytes=source.size_bytes,
                expected_sha256=source.sha256,
                attempts=3,
                timeout_seconds=300,
            ),
            root,
        )

    release.parent.mkdir(parents=True, exist_ok=True)
    partial_root = release.parent / ".partial"
    _require_within_root(root, partial_root, "VEP release partial path")
    partial_root.mkdir(mode=0o700, exist_ok=True)
    staging = partial_root / f"{release.name}.{uuid.uuid4().hex}.partial"
    staging.mkdir(mode=0o700)
    try:
        _extract_cache(source_paths["cache"], staging, contract)
        _prepare_fasta(source_paths["fasta"], staging, contract)
        for directory in (path for path in staging.rglob("*") if path.is_dir()):
            os.chmod(directory, 0o755)
        os.chmod(staging, 0o755)
        validation = _validate_release(staging, contract)
        cache_parent = staging / str(
            PurePosixPath(contract.cache.cache_root)
            .relative_to(contract.cache.output_root)
            .parent
            .parent
        )
        fasta_path = staging / str(
            PurePosixPath(contract.fasta.output_path).relative_to(
                contract.cache.output_root
            )
        )
        sentinel = sentinel_runner(contract, cache_parent, fasta_path)
        try:
            os.replace(staging, release)
        except OSError:
            if not release.exists():
                raise
            shutil.rmtree(staging, ignore_errors=True)
            validation = _validate_release(release, contract)
            cache_parent = release / str(
                PurePosixPath(contract.cache.cache_root)
                .relative_to(contract.cache.output_root)
                .parent
                .parent
            )
            sentinel = sentinel_runner(contract, cache_parent, root / contract.fasta.output_path)
        directory_fd = os.open(release.parent, os.O_RDONLY)
        try:
            os.fsync(directory_fd)
        finally:
            os.close(directory_fd)
    except Exception:
        shutil.rmtree(staging, ignore_errors=True)
        raise
    validation = _validate_release(release, contract)
    publications = _publish_resource_manifests(
        root,
        release,
        contract,
        contract_sha256,
        image_digest,
        retrieved_on or date.today().isoformat(),
    )
    for source_path in source_paths.values():
        source_path.unlink(missing_ok=True)
    bundle_digest = _publish_json(
        _bundle_data(
            contract,
            contract_sha256,
            publications,
            validation,
            sentinel,
            root,
        ),
        bundle_path,
        root,
    )
    return VepBootstrapResult(
        bundle_id=contract.bundle_id,
        bundle_manifest_path=str(bundle_path),
        bundle_manifest_sha256=bundle_digest,
        cache_manifest_path=str(publications[0][0]),
        cache_manifest_sha256=publications[0][1],
        fasta_manifest_path=str(publications[1][0]),
        fasta_manifest_sha256=publications[1][1],
        validation=validation,
        sentinel=sentinel,
        reused=False,
    )


def vep_result_as_json(result: VepBootstrapResult) -> str:
    return json.dumps(asdict(result), indent=2, sort_keys=True) + "\n"
