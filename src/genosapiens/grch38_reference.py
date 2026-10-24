"""Reconstruct and validate the masked-v2 GRCh38 plus hs38d1 references."""

from __future__ import annotations

from collections.abc import Callable, Iterator, Mapping
from dataclasses import asdict, dataclass
from datetime import date
import gzip
import hashlib
import json
import os
from pathlib import Path, PurePosixPath
import re
import shutil
import struct
import subprocess
from typing import Any, BinaryIO
from urllib.parse import urlparse
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
from genosapiens.reference_concordance import (
    ReferenceConcordanceError,
    ReferenceConstructionReport,
    ReferenceIdentityContract,
    load_construction_report,
    load_reference_identity,
    verify_reference_construction,
    verify_reference_manifest_contract,
)
from genosapiens.storage import CapacityPlan, require_capacity


CONTIG_ALIAS_SCHEMA = "genosapiens.contig-aliases"
BOOTSTRAP_TRANSFORMATION_SCHEMA = "genosapiens.grch38-bootstrap"
SEQUENCE = re.compile(rb"^[ACGTRYSWKMBDHVN]+$")
BWA_OUTPUT_SUFFIXES = {
    "0123": ".0123",
    "amb": ".amb",
    "ann": ".ann",
    "bwt": ".bwt.2bit.64",
    "pac": ".pac",
}


class Grch38BootstrapError(RuntimeError):
    """Raised when the candidate reference cannot be reproduced exactly."""


@dataclass(frozen=True, slots=True)
class FastaIndexEntry:
    name: str
    length: int
    offset: int
    line_bases: int
    line_width: int

    def as_fai_line(self) -> bytes:
        return (
            f"{self.name}\t{self.length}\t{self.offset}\t"
            f"{self.line_bases}\t{self.line_width}\n"
        ).encode("ascii")


def load_fasta_index(path: str | Path) -> Mapping[str, FastaIndexEntry]:
    index_path = Path(path)
    entries: dict[str, FastaIndexEntry] = {}
    try:
        with index_path.open("rt", encoding="ascii", newline="") as source:
            for line_number, line in enumerate(source, start=1):
                fields = line.rstrip("\r\n").split("\t")
                if len(fields) < 5:
                    raise Grch38BootstrapError(
                        f"invalid FAI line {line_number} in {index_path}"
                    )
                try:
                    entry = FastaIndexEntry(
                        fields[0], *(int(value) for value in fields[1:5])
                    )
                except ValueError as error:
                    raise Grch38BootstrapError(
                        f"invalid FAI integer on line {line_number} in {index_path}"
                    ) from error
                if (
                    not entry.name
                    or entry.length <= 0
                    or entry.offset < 0
                    or entry.line_bases <= 0
                    or entry.line_width < entry.line_bases
                    or entry.name in entries
                ):
                    raise Grch38BootstrapError(
                        f"invalid FAI entry on line {line_number} in {index_path}"
                    )
                entries[entry.name] = entry
    except OSError as error:
        raise Grch38BootstrapError(f"cannot read FAI {index_path}: {error}") from error
    if not entries:
        raise Grch38BootstrapError(f"FAI has no sequences: {index_path}")
    return entries


class IndexedFasta:
    """Minimal 1-based reference lookup backed by a standard FAI."""

    def __init__(self, fasta_path: str | Path, fai_path: str | Path) -> None:
        self.fasta_path = Path(fasta_path)
        self.index = load_fasta_index(fai_path)

    def fetch(self, contig: str, position: int, length: int) -> str:
        if (
            not isinstance(position, int)
            or isinstance(position, bool)
            or position <= 0
            or not isinstance(length, int)
            or isinstance(length, bool)
            or length <= 0
        ):
            raise ValueError("reference lookup requires positive position and length")
        try:
            entry = self.index[contig]
        except KeyError as error:
            raise KeyError(f"unknown reference contig: {contig}") from error
        zero_based = position - 1
        if zero_based + length > entry.length:
            raise IndexError(f"reference lookup exceeds contig {contig}")
        chunks: list[bytes] = []
        remaining = length
        try:
            with self.fasta_path.open("rb") as fasta:
                while remaining:
                    line_number, within_line = divmod(zero_based, entry.line_bases)
                    count = min(remaining, entry.line_bases - within_line)
                    byte_offset = (
                        entry.offset + line_number * entry.line_width + within_line
                    )
                    fasta.seek(byte_offset)
                    chunk = fasta.read(count)
                    if len(chunk) != count or not SEQUENCE.fullmatch(chunk):
                        raise Grch38BootstrapError(
                            f"invalid indexed sequence bytes for {contig}"
                        )
                    chunks.append(chunk)
                    zero_based += count
                    remaining -= count
        except OSError as error:
            raise Grch38BootstrapError(
                f"cannot read indexed FASTA {self.fasta_path}: {error}"
            ) from error
        return b"".join(chunks).decode("ascii")


@dataclass(frozen=True, slots=True)
class SequenceIdentity:
    name: str
    length: int
    md5: str


@dataclass(frozen=True, slots=True)
class _FastaRecord:
    name: str
    sequence: bytes
    body: bytes
    offset: int
    line_bases: int
    line_width: int


@dataclass(frozen=True, slots=True)
class _CompositionResult:
    public_identities: tuple[SequenceIdentity, ...]
    tool_identities: tuple[SequenceIdentity, ...]
    alignment_identities: tuple[SequenceIdentity, ...]
    tool_index: tuple[FastaIndexEntry, ...]
    alignment_index: tuple[FastaIndexEntry, ...]
    mask_delta_intervals: int
    masked_bases: int


@dataclass(frozen=True, slots=True)
class Grch38BootstrapResult:
    manifest_path: str
    manifest_sha256: str
    reference_identity_id: str
    sequence_count: int
    reused: bool


Downloader = Callable[[DownloadRequest, str | Path], DownloadResult]
BwaIndexer = Callable[[Path], None]
CapacityChecker = Callable[[str | Path, CapacityPlan], Any]


def _file_checksums(path: Path) -> tuple[str, str, int]:
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
        raise Grch38BootstrapError(f"cannot read {path}: {error}") from error
    return sha256.hexdigest(), md5.hexdigest(), size


def _file_identity(path: Path) -> tuple[str, int]:
    sha256, _md5, size = _file_checksums(path)
    return sha256, size


def _contract_sha256(path: Path) -> str:
    return _file_identity(path)[0]


def _manifest_relative_path(
    contract: ReferenceIdentityContract, contract_sha256: str
) -> str:
    """Resolve the immutable manifest for one exact reference contract."""

    return (
        f"manifests/grch38/{contract.resource_id}."
        f"{contract_sha256[:12]}.json"
    )


def _bundle_parent(contract: ReferenceIdentityContract) -> PurePosixPath:
    parents = {
        PurePosixPath(output.path).parent for output in contract.expected_outputs
    }
    if len(parents) != 1:
        raise Grch38BootstrapError("reference outputs must share one directory")
    parent = next(iter(parents))
    if any(
        not PurePosixPath(artifact.path).is_relative_to(parent)
        for artifact in contract.source_artifacts
    ):
        raise Grch38BootstrapError("retained reference sources must be under the release")
    return parent


def _load_manifest(path: Path) -> PublicResourceManifest:
    try:
        data = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as error:
        raise Grch38BootstrapError(f"cannot read manifest {path}: {error}") from error
    if not isinstance(data, Mapping):
        raise Grch38BootstrapError("reference manifest root must be an object")
    try:
        return PublicResourceManifest.from_mapping(data)
    except ManifestError as error:
        raise Grch38BootstrapError(f"invalid reference manifest: {error}") from error


def _copy_gzip(source_path: Path, destination: Path) -> None:
    try:
        destination.parent.mkdir(parents=True, exist_ok=True)
        with gzip.open(source_path, "rb") as source, destination.open("xb") as output:
            shutil.copyfileobj(source, output, length=1024 * 1024)
            output.flush()
            os.fsync(output.fileno())
    except (OSError, EOFError, gzip.BadGzipFile) as error:
        raise Grch38BootstrapError(f"cannot decompress {source_path}: {error}") from error


def _copy_file(source_path: Path, destination: Path) -> None:
    try:
        destination.parent.mkdir(parents=True, exist_ok=True)
        with source_path.open("rb") as source, destination.open("xb") as output:
            shutil.copyfileobj(source, output, length=1024 * 1024)
            output.flush()
            os.fsync(output.fileno())
    except OSError as error:
        raise Grch38BootstrapError(
            f"cannot copy {source_path} to {destination}: {error}"
        ) from error


def _record_from_parts(
    name: str | None,
    body: bytearray,
    sequence: bytearray,
    offset: int,
    line_bases: int,
    line_width: int,
    terminal_line_bases: int,
) -> _FastaRecord:
    if name is None or not sequence or not body or line_bases <= 0:
        raise Grch38BootstrapError("FASTA contains an empty record")
    if terminal_line_bases <= 0 or terminal_line_bases > line_bases:
        raise Grch38BootstrapError(f"FASTA terminal line is invalid: {name}")
    return _FastaRecord(
        name,
        bytes(sequence),
        bytes(body),
        offset,
        line_bases,
        line_width,
    )


def _fasta_records(path: Path, *, compressed: bool) -> Iterator[_FastaRecord]:
    name: str | None = None
    body = bytearray()
    sequence = bytearray()
    sequence_offset = 0
    line_bases = 0
    line_width = 0
    previous_line_bases = 0
    previous_line_width = 0
    try:
        stream = gzip.open(path, "rb") if compressed else path.open("rb")
        with stream as source:
            while line := source.readline():
                if line.startswith(b">"):
                    if name is not None:
                        yield _record_from_parts(
                            name,
                            body,
                            sequence,
                            sequence_offset,
                            line_bases,
                            line_width,
                            previous_line_bases,
                        )
                    try:
                        name = line[1:].split(maxsplit=1)[0].decode("ascii")
                    except UnicodeDecodeError as error:
                        raise Grch38BootstrapError("FASTA header is not ASCII") from error
                    if not name:
                        raise Grch38BootstrapError("FASTA header has no sequence name")
                    body = bytearray()
                    sequence = bytearray()
                    sequence_offset = source.tell()
                    line_bases = 0
                    line_width = 0
                    previous_line_bases = 0
                    previous_line_width = 0
                    continue
                value = line.rstrip(b"\r\n")
                if name is None or not value or not SEQUENCE.fullmatch(value):
                    raise Grch38BootstrapError(f"invalid FASTA sequence data in {path}")
                if previous_line_bases and (
                    previous_line_bases != line_bases
                    or previous_line_width != line_width
                ):
                    raise Grch38BootstrapError(
                        f"FASTA has a non-terminal short line: {name}"
                    )
                if line_bases == 0:
                    line_bases = len(value)
                    line_width = len(line)
                previous_line_bases = len(value)
                previous_line_width = len(line)
                body.extend(line)
                sequence.extend(value)
            if name is None:
                raise Grch38BootstrapError(f"FASTA has no records: {path}")
            yield _record_from_parts(
                name,
                body,
                sequence,
                sequence_offset,
                line_bases,
                line_width,
                previous_line_bases,
            )
    except (OSError, EOFError, gzip.BadGzipFile) as error:
        raise Grch38BootstrapError(f"cannot read FASTA {path}: {error}") from error


def _write_record(
    output: BinaryIO, record: _FastaRecord, name: str
) -> FastaIndexEntry:
    output.write(f">{name}\n".encode("ascii"))
    offset = output.tell()
    output.write(record.body)
    return FastaIndexEntry(
        name=name,
        length=len(record.sequence),
        offset=offset,
        line_bases=record.line_bases,
        line_width=record.line_width,
    )


def _mask_delta_intervals(masked: bytes, native: bytes) -> tuple[tuple[int, int], ...]:
    if len(masked) != len(native):
        raise Grch38BootstrapError("masked and unmasked sequence lengths differ")
    intervals: list[tuple[int, int]] = []
    open_start: int | None = None
    chunk_size = 1024 * 1024
    for offset in range(0, len(masked), chunk_size):
        masked_chunk = masked[offset : offset + chunk_size]
        native_chunk = native[offset : offset + chunk_size]
        if masked_chunk == native_chunk:
            if open_start is not None:
                intervals.append((open_start, offset))
                open_start = None
            continue
        for relative, (masked_base, native_base) in enumerate(
            zip(masked_chunk, native_chunk, strict=True)
        ):
            position = offset + relative
            if masked_base != native_base:
                if masked_base != ord("N"):
                    raise Grch38BootstrapError(
                        "masked component contains a non-N sequence difference"
                    )
                if open_start is None:
                    open_start = position
            elif open_start is not None:
                intervals.append((open_start, position))
                open_start = None
    if open_start is not None:
        intervals.append((open_start, len(masked)))
    return tuple(intervals)


def _compose_references(
    *,
    masked_gzip: Path,
    public_native_fasta: Path,
    masked_fai: Path,
    native_fai: Path,
    tool_fasta: Path,
    alignment_fasta: Path,
    mask_delta_path: Path,
    contract: ReferenceIdentityContract,
) -> _CompositionResult:
    masked_index = load_fasta_index(masked_fai)
    native_index = load_fasta_index(native_fai)
    masked_records = iter(_fasta_records(masked_gzip, compressed=True))
    native_records = iter(_fasta_records(public_native_fasta, compressed=False))
    public_identities: list[SequenceIdentity] = []
    tool_identities: list[SequenceIdentity] = []
    alignment_identities: list[SequenceIdentity] = []
    tool_entries: list[FastaIndexEntry] = []
    alignment_entries: list[FastaIndexEntry] = []
    seen_tool: set[str] = set()
    seen_alignment: set[str] = set()
    interval_count = 0
    masked_bases = 0
    sentinel_hits: set[tuple[str, int]] = set()
    sentinels = {
        (sentinel.contig, sentinel.position)
        for sentinel in contract.transformation.mask_sentinels
    }

    try:
        with (
            tool_fasta.open("xb") as tool_output,
            alignment_fasta.open("xb") as alignment_output,
            mask_delta_path.open("x", encoding="ascii", newline="") as delta_output,
        ):
            for index in range(contract.transformation.masked_component_records):
                try:
                    masked = next(masked_records)
                    native = next(native_records)
                except StopIteration as error:
                    raise Grch38BootstrapError(
                        "reference source ended before the masked component boundary"
                    ) from error
                if (
                    masked.name != native.name
                    or masked.name not in masked_index
                    or native.name not in native_index
                    or masked_index[masked.name].length != len(masked.sequence)
                    or native_index[native.name].length != len(native.sequence)
                    or len(masked.sequence) != len(native.sequence)
                ):
                    raise Grch38BootstrapError(
                        f"masked component identity mismatch at record {index + 1}"
                    )
                intervals = _mask_delta_intervals(masked.sequence, native.sequence)
                for start, end in intervals:
                    delta_output.write(f"{masked.name}\t{start}\t{end}\n")
                    interval_count += 1
                    masked_bases += end - start
                for contig, position in sentinels:
                    if contig == masked.name and position <= len(masked.sequence):
                        if (
                            masked.sequence[position - 1] != ord("N")
                            or masked.sequence[position - 1]
                            == native.sequence[position - 1]
                        ):
                            raise Grch38BootstrapError(
                                f"configured mask sentinel is not a measured delta: {contig}:{position}"
                            )
                        sentinel_hits.add((contig, position))
                tool_name = masked.name
                if not tool_name.startswith(contract.transformation.source_prefix):
                    raise Grch38BootstrapError(
                        f"masked component lacks chr prefix: {tool_name}"
                    )
                alignment_name = tool_name.removeprefix(
                    contract.transformation.source_prefix
                )
                if tool_name in seen_tool or alignment_name in seen_alignment:
                    raise Grch38BootstrapError("reference aliases are not unique")
                seen_tool.add(tool_name)
                seen_alignment.add(alignment_name)
                sequence_md5 = hashlib.md5(
                    masked.sequence, usedforsecurity=False
                ).hexdigest()
                native_md5 = hashlib.md5(
                    native.sequence, usedforsecurity=False
                ).hexdigest()
                tool_entries.append(_write_record(tool_output, masked, tool_name))
                alignment_entries.append(
                    _write_record(alignment_output, masked, alignment_name)
                )
                public_identities.append(
                    SequenceIdentity(native.name, len(native.sequence), native_md5)
                )
                tool_identities.append(
                    SequenceIdentity(tool_name, len(masked.sequence), sequence_md5)
                )
                alignment_identities.append(
                    SequenceIdentity(alignment_name, len(masked.sequence), sequence_md5)
                )
            try:
                next(masked_records)
            except StopIteration:
                pass
            else:
                raise Grch38BootstrapError(
                    "masked source contains records beyond the 195-record component"
                )
            if tool_identities[-1].name != contract.transformation.masked_component_terminal:
                raise Grch38BootstrapError("masked component does not end at chrEBV")

            for index in range(contract.transformation.decoy_suffix_records):
                try:
                    native = next(native_records)
                except StopIteration as error:
                    raise Grch38BootstrapError(
                        "NCBI source ended before the hs38d1 suffix was complete"
                    ) from error
                if not native.name.endswith(contract.transformation.decoy_name_suffix):
                    raise Grch38BootstrapError(
                        f"NCBI suffix record is not an hs38d1 decoy: {native.name}"
                    )
                if (
                    native.name not in native_index
                    or native_index[native.name].length != len(native.sequence)
                ):
                    raise Grch38BootstrapError(
                        f"NCBI decoy index mismatch: {native.name}"
                    )
                tool_name = native.name
                if not tool_name.startswith(contract.transformation.source_prefix):
                    raise Grch38BootstrapError(
                        f"NCBI decoy lacks chr prefix: {tool_name}"
                    )
                alignment_name = tool_name.removeprefix(
                    contract.transformation.source_prefix
                )
                if tool_name in seen_tool or alignment_name in seen_alignment:
                    raise Grch38BootstrapError("reference aliases are not unique")
                seen_tool.add(tool_name)
                seen_alignment.add(alignment_name)
                sequence_md5 = hashlib.md5(
                    native.sequence, usedforsecurity=False
                ).hexdigest()
                tool_entries.append(_write_record(tool_output, native, tool_name))
                alignment_entries.append(
                    _write_record(alignment_output, native, alignment_name)
                )
                identity = SequenceIdentity(tool_name, len(native.sequence), sequence_md5)
                public_identities.append(identity)
                tool_identities.append(identity)
                alignment_identities.append(
                    SequenceIdentity(alignment_name, len(native.sequence), sequence_md5)
                )
            try:
                next(native_records)
            except StopIteration:
                pass
            else:
                raise Grch38BootstrapError(
                    "NCBI source contains records beyond the expected hs38d1 suffix"
                )
            if tuple(masked_index) != tuple(
                identity.name
                for identity in tool_identities[
                    : contract.transformation.masked_component_records
                ]
            ):
                raise Grch38BootstrapError(
                    "masked component FASTA and FAI order differ"
                )
            if tuple(native_index) != tuple(
                identity.name for identity in public_identities
            ):
                raise Grch38BootstrapError("NCBI FASTA and FAI order differ")
            if sentinel_hits != sentinels:
                raise Grch38BootstrapError("one or more mask sentinels were not evaluated")
            for output in (tool_output, alignment_output):
                output.flush()
                os.fsync(output.fileno())
            delta_output.flush()
            os.fsync(delta_output.fileno())
    except OSError as error:
        raise Grch38BootstrapError(f"cannot compose reference FASTAs: {error}") from error

    if (
        interval_count == 0
        or masked_bases == 0
        or len(tool_identities) != contract.sequence_count
    ):
        raise Grch38BootstrapError("reference composition did not measure the required mask")
    return _CompositionResult(
        tuple(public_identities),
        tuple(tool_identities),
        tuple(alignment_identities),
        tuple(tool_entries),
        tuple(alignment_entries),
        interval_count,
        masked_bases,
    )


def _write_fai(path: Path, entries: tuple[FastaIndexEntry, ...]) -> None:
    try:
        with path.open("xb") as output:
            for entry in entries:
                output.write(entry.as_fai_line())
            output.flush()
            os.fsync(output.fileno())
    except OSError as error:
        raise Grch38BootstrapError(f"cannot write FASTA index {path}: {error}") from error


def _write_dictionary(path: Path, identities: tuple[SequenceIdentity, ...]) -> None:
    try:
        with path.open("x", encoding="ascii", newline="") as output:
            output.write("@HD\tVN:1.6\tSO:unsorted\n")
            for identity in identities:
                output.write(
                    f"@SQ\tSN:{identity.name}\tLN:{identity.length}\tM5:{identity.md5}"
                    "\tAS:GRCh38-no-alt-plus-hs38d1\n"
                )
            output.flush()
            os.fsync(output.fileno())
    except OSError as error:
        raise Grch38BootstrapError(f"cannot write sequence dictionary {path}: {error}") from error


def _write_m5(path: Path, identities: tuple[SequenceIdentity, ...]) -> None:
    try:
        with path.open("x", encoding="ascii", newline="") as output:
            output.write("sequence\tlength\tm5\n")
            for identity in identities:
                output.write(f"{identity.name}\t{identity.length}\t{identity.md5}\n")
            output.flush()
            os.fsync(output.fileno())
    except OSError as error:
        raise Grch38BootstrapError(f"cannot write M5 table {path}: {error}") from error


def _write_aliases(
    path: Path,
    contract: ReferenceIdentityContract,
    composition: _CompositionResult,
) -> None:
    tool_to_alignment = {
        tool.name: alignment.name
        for tool, alignment in zip(
            composition.tool_identities,
            composition.alignment_identities,
            strict=True,
        )
    }
    document = {
        "schema": CONTIG_ALIAS_SCHEMA,
        "contract_id": contract.contract_id,
        "reference_roles": {
            role.role: {
                "identity_id": role.identity_id,
                "naming": role.naming,
                "sequence_relationship": role.sequence_relationship,
            }
            for role in contract.reference_roles
        },
        "sequence_count": contract.sequence_count,
        "tool_to_alignment": tool_to_alignment,
        "alignment_to_tool": {
            alignment: tool for tool, alignment in tool_to_alignment.items()
        },
        "sequence_lengths": {
            identity.name: identity.length
            for identity in composition.alignment_identities
        },
        "challenge_sequence_md5": {
            identity.name: identity.md5
            for identity in composition.alignment_identities
        },
        "public_native_sequence_md5": {
            identity.name: identity.md5
            for identity in composition.public_identities
        },
        "canonical_submission_contigs": [
            *(str(number) for number in range(1, 23)),
            "X",
            "Y",
            "M",
        ],
    }
    try:
        with path.open("x", encoding="utf-8", newline="") as output:
            json.dump(document, output, indent=2, sort_keys=True)
            output.write("\n")
            output.flush()
            os.fsync(output.fileno())
    except OSError as error:
        raise Grch38BootstrapError(f"cannot write contig aliases {path}: {error}") from error


def _default_bwa_indexer(fasta: Path) -> None:
    try:
        subprocess.run(
            ["bwa-mem2", "index", str(fasta)],
            check=True,
            stdout=subprocess.DEVNULL,
        )
    except (OSError, subprocess.CalledProcessError) as error:
        raise Grch38BootstrapError(f"BWA-MEM2 indexing failed for {fasta}: {error}") from error


def _validate_gzi(path: Path) -> None:
    try:
        data = path.read_bytes()
    except OSError as error:
        raise Grch38BootstrapError(f"cannot read GZI {path}: {error}") from error
    if len(data) < 8:
        raise Grch38BootstrapError("masked component GZI is truncated")
    count = struct.unpack("<Q", data[:8])[0]
    if len(data) != 8 + count * 16:
        raise Grch38BootstrapError("masked component GZI size is inconsistent")
    previous = (-1, -1)
    for offset in range(8, len(data), 16):
        current = struct.unpack("<QQ", data[offset : offset + 16])
        if current[0] <= previous[0] or current[1] <= previous[1]:
            raise Grch38BootstrapError("masked component GZI offsets are unordered")
        previous = current


def _validate_source(path: Path, expected: Any) -> None:
    sha256, md5, size = _file_checksums(path)
    if (
        sha256 != expected.sha256
        or md5 != expected.md5
        or size != expected.size_bytes
    ):
        raise Grch38BootstrapError(f"retained source identity mismatch: {expected.role}")


def _load_dictionary(path: Path) -> Mapping[str, tuple[int, str]]:
    values: dict[str, tuple[int, str]] = {}
    try:
        with path.open("rt", encoding="ascii") as source:
            if source.readline() != "@HD\tVN:1.6\tSO:unsorted\n":
                raise Grch38BootstrapError(f"invalid dictionary header: {path}")
            for line in source:
                fields = dict(
                    value.split(":", 1) for value in line.rstrip("\n").split("\t")[1:]
                )
                values[fields["SN"]] = (int(fields["LN"]), fields["M5"])
    except (KeyError, OSError, ValueError) as error:
        raise Grch38BootstrapError(f"invalid dictionary {path}: {error}") from error
    return values


def _load_m5(path: Path) -> Mapping[str, tuple[int, str]]:
    values: dict[str, tuple[int, str]] = {}
    try:
        with path.open("rt", encoding="ascii") as source:
            if source.readline() != "sequence\tlength\tm5\n":
                raise Grch38BootstrapError(f"invalid M5 table header: {path}")
            for line in source:
                name, length, md5 = line.rstrip("\n").split("\t")
                values[name] = (int(length), md5)
    except (OSError, ValueError) as error:
        raise Grch38BootstrapError(f"invalid M5 table {path}: {error}") from error
    return values


def _inspect_indexed_fasta(
    fasta_path: Path,
    index: Mapping[str, FastaIndexEntry],
) -> tuple[SequenceIdentity, ...]:
    identities: list[SequenceIdentity] = []
    for record in _fasta_records(fasta_path, compressed=False):
        expected_entry = FastaIndexEntry(
            record.name,
            len(record.sequence),
            record.offset,
            record.line_bases,
            record.line_width,
        )
        if index.get(record.name) != expected_entry:
            raise Grch38BootstrapError(
                f"FASTA index lookup metadata differs: {record.name}"
            )
        identities.append(
            SequenceIdentity(
                record.name,
                len(record.sequence),
                hashlib.md5(record.sequence, usedforsecurity=False).hexdigest(),
            )
        )
    if tuple(index) != tuple(identity.name for identity in identities):
        raise Grch38BootstrapError("FASTA and FAI names or order differ")
    return tuple(identities)


def _inspect_source_fasta(
    fasta_path: Path, *, compressed: bool
) -> tuple[SequenceIdentity, ...]:
    return tuple(
        SequenceIdentity(
            record.name,
            len(record.sequence),
            hashlib.md5(record.sequence, usedforsecurity=False).hexdigest(),
        )
        for record in _fasta_records(fasta_path, compressed=compressed)
    )


def _load_bwa_metadata(
    ann_path: Path, amb_path: Path
) -> tuple[int, Mapping[str, int]]:
    try:
        with ann_path.open("rt", encoding="ascii") as source:
            header = source.readline().split()
            if len(header) != 3:
                raise ValueError("invalid ANN header")
            total_length, sequence_count = (int(value) for value in header[:2])
            values: dict[str, int] = {}
            expected_offset = 0
            for _index in range(sequence_count):
                identity = source.readline().rstrip("\n").split(maxsplit=2)
                coordinates = source.readline().split()
                if len(identity) < 2 or len(coordinates) != 3:
                    raise ValueError("truncated ANN sequence metadata")
                offset, length, _ambiguities = (
                    int(value) for value in coordinates
                )
                name = identity[1]
                if offset != expected_offset or length <= 0 or name in values:
                    raise ValueError("invalid ANN sequence metadata")
                values[name] = length
                expected_offset += length
            if source.read().strip() or expected_offset != total_length:
                raise ValueError("inconsistent ANN length or trailing data")
        with amb_path.open("rt", encoding="ascii") as source:
            amb_header = source.readline().split()
            if (
                len(amb_header) != 3
                or int(amb_header[0]) != total_length
                or int(amb_header[1]) != sequence_count
                or int(amb_header[2]) < 0
            ):
                raise ValueError("invalid AMB header")
    except (OSError, UnicodeDecodeError, ValueError) as error:
        raise Grch38BootstrapError(
            f"invalid BWA-MEM2 index metadata: {error}"
        ) from error
    return total_length, values


def _expected_identity_map(
    index: Mapping[str, FastaIndexEntry], m5: Mapping[str, tuple[int, str]]
) -> Mapping[str, tuple[int, str]]:
    if set(index) != set(m5):
        raise Grch38BootstrapError("FAI and M5 sequence names differ")
    if any(index[name].length != m5[name][0] for name in index):
        raise Grch38BootstrapError("FAI and M5 sequence lengths differ")
    return m5


def _validate_bundle_files(
    contract: ReferenceIdentityContract,
    reference_root: Path,
    *,
    contract_sha256: str,
    allow_order_only_difference: bool = False,
) -> int:
    root = reference_root
    for artifact in contract.source_artifacts:
        _validate_source(root / artifact.path, artifact)
    _validate_gzi(root / contract.source_artifact("masked_component_gzi").path)
    if (
        (root / contract.expected_output("public_native_fai").path).read_bytes()
        != (root / contract.source_artifact("unmasked_decoy_source_fai").path).read_bytes()
    ):
        raise Grch38BootstrapError("fresh public-native FAI differs from its source index")

    role_values: dict[str, Mapping[str, tuple[int, str]]] = {}
    indexes: dict[str, Mapping[str, FastaIndexEntry]] = {}
    identities_by_role: dict[str, tuple[SequenceIdentity, ...]] = {}
    for role_name in (
        "public_annotation_benchmark",
        "challenge_tool_view",
        "challenge_alignment",
    ):
        role = contract.reference_role(role_name)
        index = load_fasta_index(root / contract.expected_output(role.fai_output_role).path)
        identities = _inspect_indexed_fasta(
            root / contract.expected_output(role.fasta_output_role).path,
            index,
        )
        m5 = _load_m5(root / contract.expected_output(role.m5_output_role).path)
        expected = _expected_identity_map(index, m5)
        if expected != {
            identity.name: (identity.length, identity.md5)
            for identity in identities
        }:
            raise Grch38BootstrapError(f"FASTA and M5 differ for role {role_name}")
        if _load_dictionary(
            root / contract.expected_output(role.dictionary_output_role).path
        ) != expected:
            raise Grch38BootstrapError(f"dictionary and M5 differ for role {role_name}")
        if len(index) != contract.sequence_count:
            raise Grch38BootstrapError(f"sequence count differs for role {role_name}")
        indexes[role_name] = index
        role_values[role_name] = expected
        identities_by_role[role_name] = identities

    masked_source_identities = _inspect_source_fasta(
        root / contract.source_artifact("masked_component_fasta_gzip").path,
        compressed=True,
    )
    native_source_identities = _inspect_source_fasta(
        root / contract.source_artifact("unmasked_decoy_source_fasta_gzip").path,
        compressed=True,
    )
    masked_records = contract.transformation.masked_component_records
    expected_tool_identities = (
        *masked_source_identities,
        *native_source_identities[masked_records:],
    )
    tool_identity_matches_source = (
        {
            identity.name: (identity.length, identity.md5)
            for identity in identities_by_role["challenge_tool_view"]
        }
        == {
            identity.name: (identity.length, identity.md5)
            for identity in expected_tool_identities
        }
        if allow_order_only_difference
        else identities_by_role["challenge_tool_view"] == expected_tool_identities
    )
    if (
        len(masked_source_identities) != masked_records
        or len(native_source_identities) != contract.sequence_count
        or identities_by_role["public_annotation_benchmark"]
        != native_source_identities
        or not tool_identity_matches_source
    ):
        raise Grch38BootstrapError(
            "installed reference sequences differ from the pinned source objects"
        )
    expected_alignment_identities = tuple(
        SequenceIdentity(
            identity.name.removeprefix(contract.transformation.source_prefix),
            identity.length,
            identity.md5,
        )
        for identity in identities_by_role["challenge_tool_view"]
    )
    if identities_by_role["challenge_alignment"] != expected_alignment_identities:
        raise Grch38BootstrapError(
            "challenge alignment sequences differ from the pinned source objects"
        )

    tool = indexes["challenge_tool_view"]
    alignment = indexes["challenge_alignment"]
    expected_tool_names = tuple(f"chr{name}" for name in alignment)
    if tuple(tool) != expected_tool_names:
        raise Grch38BootstrapError("challenge reference names or order differ")
    for tool_name, alignment_name in zip(tool, alignment, strict=True):
        if (
            tool[tool_name].length != alignment[alignment_name].length
            or role_values["challenge_tool_view"][tool_name]
            != role_values["challenge_alignment"][alignment_name]
        ):
            raise Grch38BootstrapError(
                f"challenge reference views are not sequence-identical: {alignment_name}"
            )
    for name, expected_length in contract.canonical_contig_lengths.items():
        if alignment.get(name, FastaIndexEntry("", 0, 0, 1, 1)).length != expected_length:
            raise Grch38BootstrapError(f"canonical contig length mismatch: {name}")
    public_identities = identities_by_role["public_annotation_benchmark"]
    tool_identities = identities_by_role["challenge_tool_view"]
    if not allow_order_only_difference:
        if tuple(
            identity.name for identity in tool_identities[:masked_records]
        ) != tuple(identity.name for identity in public_identities[:masked_records]):
            raise Grch38BootstrapError(
                "masked and public component names or order differ"
            )
        if any(
            tool_identity.length != public_identity.length
            for tool_identity, public_identity in zip(
                tool_identities[:masked_records],
                public_identities[:masked_records],
                strict=True,
            )
        ):
            raise Grch38BootstrapError("masked and public component lengths differ")
        if tool_identities[masked_records - 1].name != "chrEBV":
            raise Grch38BootstrapError("masked component output does not end at chrEBV")
        if tuple(tool_identities[masked_records:]) != tuple(
            public_identities[masked_records:]
        ) or any(
            not identity.name.endswith(contract.transformation.decoy_name_suffix)
            for identity in tool_identities[masked_records:]
        ):
            raise Grch38BootstrapError(
                "appended decoy identities differ from NCBI source"
            )

    alias_path = root / contract.expected_output("contig_aliases").path
    try:
        aliases = json.loads(alias_path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as error:
        raise Grch38BootstrapError(f"cannot read contig aliases: {error}") from error
    if (
        not isinstance(aliases, Mapping)
        or aliases.get("schema") != CONTIG_ALIAS_SCHEMA
        or aliases.get("contract_id") != contract.contract_id
        or aliases.get("sequence_count") != contract.sequence_count
        or aliases.get("tool_to_alignment")
        != {tool_name: name for tool_name, name in zip(tool, alignment, strict=True)}
    ):
        raise Grch38BootstrapError("contig alias contract is invalid")

    for prefix in ("challenge_tool", "challenge_alignment"):
        for suffix in BWA_OUTPUT_SUFFIXES:
            path = root / contract.expected_output(f"{prefix}_bwa_{suffix}").path
            if not path.is_file() or path.stat().st_size <= 0:
                raise Grch38BootstrapError(f"BWA-MEM2 index is missing: {path}")
        _total_length, bwa_sequences = _load_bwa_metadata(
            root / contract.expected_output(f"{prefix}_bwa_ann").path,
            root / contract.expected_output(f"{prefix}_bwa_amb").path,
        )
        role_name = (
            "challenge_tool_view"
            if prefix == "challenge_tool"
            else "challenge_alignment"
        )
        if bwa_sequences != {
            name: entry.length for name, entry in indexes[role_name].items()
        }:
            raise Grch38BootstrapError(
                f"BWA-MEM2 and FASTA index metadata differ for {role_name}"
            )

    report = load_construction_report(
        root / contract.expected_output("construction_report").path
    )
    expected_output_sha256 = {
        output.role: _file_identity(root / output.path)[0]
        for output in contract.expected_outputs
        if output.role != "construction_report"
    }
    if (
        report.contract_id != contract.contract_id
        or report.contract_sha256 != contract_sha256
        or report.output_sha256 != expected_output_sha256
        or report.source_sha256
        != {artifact.role: artifact.sha256 for artifact in contract.source_artifacts}
        or not report.accepted
    ):
        raise Grch38BootstrapError("construction report does not match installed files")
    return len(alignment)


def validate_grch38_bundle(
    contract: ReferenceIdentityContract,
    reference_root: str | Path,
    *,
    identity_path: str | Path | None = None,
    allow_order_only_difference: bool = False,
) -> int:
    root = Path(reference_root).resolve()
    contract_file = Path(identity_path) if identity_path is not None else None
    if contract_file is None:
        raise Grch38BootstrapError("identity_path is required for immutable validation")
    contract_sha256 = _contract_sha256(contract_file)
    manifest_path = root / _manifest_relative_path(contract, contract_sha256)
    manifest = _load_manifest(manifest_path)
    try:
        verify_reference_manifest_contract(
            contract, manifest, contract_sha256=contract_sha256
        )
        verify_outputs(manifest, root)
        verify_reference_construction(
            contract, manifest, root, contract_sha256=contract_sha256
        )
    except (ManifestError, ReferenceConcordanceError) as error:
        raise Grch38BootstrapError(str(error)) from error
    return _validate_bundle_files(
        contract,
        root,
        contract_sha256=contract_sha256,
        allow_order_only_difference=allow_order_only_difference,
    )


def _resource_outputs(
    root: Path, contract: ReferenceIdentityContract
) -> tuple[ResourceOutput, ...]:
    paths = [
        *(artifact.path for artifact in contract.source_artifacts),
        *(output.path for output in contract.expected_outputs),
    ]
    outputs: list[ResourceOutput] = []
    for relative_path in paths:
        digest, size = _file_identity(root / relative_path)
        outputs.append(
            ResourceOutput(
                path=relative_path,
                size_bytes=size,
                sha256=digest,
                retention=RetentionPolicy("durable", None),
            )
        )
    return tuple(outputs)


def _make_manifest(
    contract: ReferenceIdentityContract,
    identity_path: Path,
    reference_root: Path,
    image_digest: str,
    retrieved_on: str,
) -> PublicResourceManifest:
    source = contract.canonical_source
    return PublicResourceManifest(
        resource_id=contract.resource_id,
        release=contract.release,
        source=SourceIdentity(url=source.url, sha256=source.sha256),
        license="GIAB-and-NCBI-public-data",
        assembly=contract.assembly,
        retrieved_on=retrieved_on,
        transformation=Transformation(
            schema=BOOTSTRAP_TRANSFORMATION_SCHEMA,
            parameters_sha256=_contract_sha256(identity_path),
        ),
        image_digest=image_digest,
        outputs=_resource_outputs(reference_root, contract),
    )


def _download_source(
    contract: ReferenceIdentityContract,
    role: str,
    scratch_root: Path,
    reference_root: Path,
    downloader: Downloader,
) -> Path:
    artifact = contract.source_artifact(role)
    destination = scratch_root / Path(urlparse(artifact.url).path).name
    downloader(
        DownloadRequest(
            url=artifact.url,
            destination=destination,
            expected_size_bytes=artifact.size_bytes,
            expected_sha256=artifact.sha256,
            attempts=3,
            timeout_seconds=120,
        ),
        reference_root,
    )
    _validate_source(destination, artifact)
    return destination


def _write_construction_report(
    path: Path,
    contract: ReferenceIdentityContract,
    contract_sha256: str,
    composition: _CompositionResult,
    staged_outputs: Mapping[str, Path],
) -> None:
    output_sha256 = {
        role: _file_identity(output_path)[0]
        for role, output_path in staged_outputs.items()
        if role != "construction_report"
    }
    report = ReferenceConstructionReport(
        contract_id=contract.contract_id,
        contract_sha256=contract_sha256,
        challenge_alignment_identity_id=contract.reference_role(
            "challenge_alignment"
        ).identity_id,
        challenge_tool_view_identity_id=contract.reference_role(
            "challenge_tool_view"
        ).identity_id,
        public_annotation_benchmark_identity_id=contract.reference_role(
            "public_annotation_benchmark"
        ).identity_id,
        masked_component_records=contract.transformation.masked_component_records,
        decoy_suffix_records=contract.transformation.decoy_suffix_records,
        total_records=contract.sequence_count,
        mask_delta_intervals=composition.mask_delta_intervals,
        masked_bases=composition.masked_bases,
        mask_delta_sha256=output_sha256["mask_delta"],
        source_sha256={
            artifact.role: artifact.sha256 for artifact in contract.source_artifacts
        },
        output_sha256=output_sha256,
        accepted=True,
    )
    try:
        with path.open("x", encoding="utf-8", newline="") as output:
            json.dump(asdict(report), output, indent=2, sort_keys=True)
            output.write("\n")
            output.flush()
            os.fsync(output.fileno())
    except OSError as error:
        raise Grch38BootstrapError(f"cannot write construction report: {error}") from error


def bootstrap_grch38_reference(
    *,
    identity_path: str | Path,
    reference_root: str | Path,
    image_digest: str,
    retrieved_on: str | None = None,
    scratch_root: str | Path | None = None,
    reserve_bytes: int = 100_000_000_000,
    downloader: Downloader = download_atomic,
    bwa_indexer: BwaIndexer = _default_bwa_indexer,
    capacity_checker: CapacityChecker = require_capacity,
) -> Grch38BootstrapResult:
    identity_file = Path(identity_path)
    root = Path(reference_root).resolve(strict=True)
    try:
        contract = load_reference_identity(identity_file)
    except ReferenceConcordanceError as error:
        raise Grch38BootstrapError(str(error)) from error
    contract_sha256 = _contract_sha256(identity_file)
    manifest_path = root / _manifest_relative_path(contract, contract_sha256)
    if manifest_path.exists():
        try:
            sequence_count = validate_grch38_bundle(
                contract, root, identity_path=identity_file
            )
        except Grch38BootstrapError as error:
            raise Grch38BootstrapError(
                f"completed GRCh38 bundle is invalid: {error}"
            ) from error
        return Grch38BootstrapResult(
            str(manifest_path),
            _file_identity(manifest_path)[0],
            contract.identity_id,
            sequence_count,
            True,
        )

    source_bytes = sum(artifact.size_bytes for artifact in contract.source_artifacts)
    capacity_checker(
        root,
        CapacityPlan(
            operation="reconstruct masked-v2 GRCh38 plus hs38d1",
            source_bytes=source_bytes,
            output_bytes=50_000_000_000,
            scratch_bytes=source_bytes,
            reserve_bytes=reserve_bytes,
        ),
    )
    scratch_base = (
        Path(scratch_root).resolve()
        if scratch_root is not None
        else root / "scratch/public-build"
    )
    if not scratch_base.is_relative_to(root):
        raise Grch38BootstrapError("reference scratch must remain under reference_root")
    scratch = scratch_base / "grch38"
    scratch.mkdir(parents=True, exist_ok=True, mode=0o700)
    downloaded = {
        artifact.role: _download_source(
            contract, artifact.role, scratch, root, downloader
        )
        for artifact in contract.source_artifacts
    }

    target_relative = _bundle_parent(contract)
    target = root / target_relative
    if not target.exists():
        target.parent.mkdir(parents=True, exist_ok=True)
        partial_root = target.parent / ".partial"
        partial_root.mkdir(mode=0o700, exist_ok=True)
        staging = partial_root / f"{target.name}.{uuid.uuid4().hex}.partial"
        staging.mkdir(mode=0o700)

        def staged_path(relative: str) -> Path:
            path = staging / PurePosixPath(relative).relative_to(target_relative)
            path.parent.mkdir(parents=True, exist_ok=True)
            return path

        staged_outputs = {
            output.role: staged_path(output.path) for output in contract.expected_outputs
        }
        try:
            for artifact in contract.source_artifacts:
                _copy_file(downloaded[artifact.role], staged_path(artifact.path))
            _copy_gzip(
                downloaded["unmasked_decoy_source_fasta_gzip"],
                staged_outputs["public_native_fasta"],
            )
            _copy_file(
                downloaded["unmasked_decoy_source_fai"],
                staged_outputs["public_native_fai"],
            )
            composition = _compose_references(
                masked_gzip=downloaded["masked_component_fasta_gzip"],
                public_native_fasta=staged_outputs["public_native_fasta"],
                masked_fai=downloaded["masked_component_fai"],
                native_fai=downloaded["unmasked_decoy_source_fai"],
                tool_fasta=staged_outputs["challenge_tool_fasta"],
                alignment_fasta=staged_outputs["challenge_alignment_fasta"],
                mask_delta_path=staged_outputs["mask_delta"],
                contract=contract,
            )
            _write_fai(staged_outputs["challenge_tool_fai"], composition.tool_index)
            _write_fai(
                staged_outputs["challenge_alignment_fai"],
                composition.alignment_index,
            )
            for prefix, identities in (
                ("public_native", composition.public_identities),
                ("challenge_tool", composition.tool_identities),
                ("challenge_alignment", composition.alignment_identities),
            ):
                _write_dictionary(staged_outputs[f"{prefix}_dict"], identities)
                _write_m5(staged_outputs[f"{prefix}_m5"], identities)
            _write_aliases(staged_outputs["contig_aliases"], contract, composition)
            for prefix in ("challenge_tool", "challenge_alignment"):
                bwa_indexer(staged_outputs[f"{prefix}_fasta"])
                for suffix in BWA_OUTPUT_SUFFIXES:
                    expected = staged_outputs[f"{prefix}_bwa_{suffix}"]
                    if not expected.is_file() or expected.stat().st_size <= 0:
                        raise Grch38BootstrapError(
                            f"BWA-MEM2 did not create expected index: {expected.name}"
                        )
            _write_construction_report(
                staged_outputs["construction_report"],
                contract,
                contract_sha256,
                composition,
                staged_outputs,
            )
            for path in staging.rglob("*"):
                if path.is_file():
                    os.chmod(path, 0o444)
            os.replace(staging, target)
            directory_fd = os.open(target.parent, os.O_RDONLY)
            try:
                os.fsync(directory_fd)
            finally:
                os.close(directory_fd)
        except Exception:
            shutil.rmtree(staging, ignore_errors=True)
            raise
    elif not target.is_dir():
        raise Grch38BootstrapError(f"reference release path is not a directory: {target}")

    sequence_count = _validate_bundle_files(
        contract, root, contract_sha256=contract_sha256
    )
    manifest = _make_manifest(
        contract,
        identity_file,
        root,
        image_digest,
        retrieved_on or date.today().isoformat(),
    )
    publication = publish_manifest(manifest, manifest_path, root)
    validate_grch38_bundle(contract, root, identity_path=identity_file)
    for path in downloaded.values():
        path.unlink(missing_ok=True)
    return Grch38BootstrapResult(
        str(publication.path),
        publication.manifest_sha256,
        contract.identity_id,
        sequence_count,
        publication.reused,
    )


def bootstrap_result_as_json(result: Grch38BootstrapResult) -> str:
    return json.dumps(asdict(result), indent=2, sort_keys=True) + "\n"
