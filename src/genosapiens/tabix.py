"""Small dependency-free reader for validating BGZF VCF and TBI queries."""

from __future__ import annotations

from dataclasses import dataclass
import gzip
import io
from pathlib import Path
import struct
from typing import BinaryIO


class TabixError(RuntimeError):
    """Raised when a BGZF stream or tabix index is malformed."""


@dataclass(frozen=True, slots=True)
class TabixReference:
    name: str
    bins: dict[int, tuple[tuple[int, int], ...]]
    linear_offsets: tuple[int, ...]


def _take(data: bytes, offset: int, size: int) -> tuple[bytes, int]:
    end = offset + size
    if end > len(data):
        raise TabixError("truncated tabix index")
    return data[offset:end], end


def _unpack(data: bytes, offset: int, fmt: str) -> tuple[tuple[int, ...], int]:
    size = struct.calcsize(fmt)
    value, end = _take(data, offset, size)
    return struct.unpack(fmt, value), end


class TabixIndex:
    def __init__(
        self,
        references: tuple[TabixReference, ...],
        *,
        format_code: int,
        sequence_column: int,
        begin_column: int,
        end_column: int,
        meta_character: int,
        skipped_lines: int,
    ) -> None:
        self.references = references
        self.by_name = {reference.name: reference for reference in references}
        self.format_code = format_code
        self.sequence_column = sequence_column
        self.begin_column = begin_column
        self.end_column = end_column
        self.meta_character = meta_character
        self.skipped_lines = skipped_lines

    @classmethod
    def load(cls, path: str | Path, *, budget=None) -> "TabixIndex":
        index_path = Path(path)
        try:
            if budget is None:
                compressed = index_path.read_bytes()
            else:
                budget.check()
                with index_path.open("rb") as stream:
                    compressed = stream.read(min(8 * 1024**2, budget.max_row_bytes) + 1)
                budget.charge(len(compressed), "tabix_transport")
                if len(compressed) > 8 * 1024**2:
                    raise TabixError("tabix index exceeds the read limit")
            data = (_bounded_gunzip(compressed, 32 * 1024**2, budget, "tabix_index")
                    if compressed[:2] == b"\x1f\x8b" else compressed)
        except (OSError, gzip.BadGzipFile) as error:
            raise TabixError(f"cannot read tabix index {index_path}: {error}") from error
        if data[:4] != b"TBI\x01":
            raise TabixError("tabix index has invalid magic")
        offset = 4
        header, offset = _unpack(data, offset, "<8i")
        (
            reference_count,
            format_code,
            sequence_column,
            begin_column,
            end_column,
            meta_character,
            skipped_lines,
            names_length,
        ) = header
        if reference_count <= 0 or names_length <= 0:
            raise TabixError("tabix index has invalid header counts")
        names_data, offset = _take(data, offset, names_length)
        try:
            names = tuple(
                name.decode("ascii") for name in names_data.rstrip(b"\x00").split(b"\x00")
            )
        except UnicodeDecodeError as error:
            raise TabixError("tabix reference names are not ASCII") from error
        if len(names) != reference_count or len(set(names)) != reference_count:
            raise TabixError("tabix reference name count is invalid")
        references: list[TabixReference] = []
        for name in names:
            if budget is not None:
                budget.check()
            (bin_count,), offset = _unpack(data, offset, "<i")
            if bin_count < 0:
                raise TabixError("tabix bin count is negative")
            bins: dict[int, tuple[tuple[int, int], ...]] = {}
            for _ in range(bin_count):
                (bin_number, chunk_count), offset = _unpack(data, offset, "<Ii")
                if chunk_count < 0 or bin_number in bins:
                    raise TabixError("tabix bin entry is invalid")
                chunks: list[tuple[int, int]] = []
                for _ in range(chunk_count):
                    (chunk_begin, chunk_end), offset = _unpack(data, offset, "<QQ")
                    # Bin 37450 stores mapped/unmapped record counts in the
                    # chunk fields rather than BGZF virtual offsets.
                    if bin_number != 37450 and chunk_end < chunk_begin:
                        raise TabixError("tabix chunk offsets are invalid")
                    chunks.append((chunk_begin, chunk_end))
                bins[bin_number] = tuple(chunks)
            (linear_count,), offset = _unpack(data, offset, "<i")
            if linear_count < 0:
                raise TabixError("tabix linear index count is negative")
            linear_offsets: tuple[int, ...] = ()
            if linear_count:
                linear_offsets, offset = _unpack(data, offset, f"<{linear_count}Q")
            references.append(TabixReference(name, bins, linear_offsets))
        if len(data) - offset not in {0, 8}:
            raise TabixError("tabix index has unexpected trailing bytes")
        return cls(
            tuple(references),
            format_code=format_code,
            sequence_column=sequence_column,
            begin_column=begin_column,
            end_column=end_column,
            meta_character=meta_character,
            skipped_lines=skipped_lines,
        )

    def require_vcf_layout(self) -> None:
        if (
            self.format_code & 0xFFFF != 2
            or self.sequence_column != 1
            or self.begin_column != 2
            or self.meta_character != ord("#")
            or self.skipped_lines != 0
        ):
            raise TabixError("tabix index does not declare the expected VCF layout")

    def query(
        self,
        bgzf_path: str | Path,
        contig: str,
        start: int,
        end: int,
        *,
        max_records: int = 100,
        budget=None,
        preserve_duplicates: bool = False,
    ) -> tuple[str, ...]:
        if start <= 0 or end < start or max_records <= 0:
            raise ValueError("tabix query coordinates and max_records must be positive")
        try:
            reference = self.by_name[contig]
        except KeyError as error:
            raise KeyError(f"unknown tabix contig: {contig}") from error
        zero_start = start - 1
        zero_end = end
        bins = _region_bins(zero_start, zero_end)
        minimum_offset = 0
        linear_bin = zero_start >> 14
        if linear_bin < len(reference.linear_offsets):
            minimum_offset = reference.linear_offsets[linear_bin]
        chunks = sorted(
            chunk
            for bin_number in bins
            for chunk in reference.bins.get(bin_number, ())
            if chunk[1] > minimum_offset
        )
        merged: list[tuple[int, int]] = []
        for begin, chunk_end in chunks:
            begin = max(begin, minimum_offset)
            if not merged or begin > merged[-1][1]:
                merged.append((begin, chunk_end))
            else:
                merged[-1] = (merged[-1][0], max(merged[-1][1], chunk_end))
        records: list[str] = []
        seen: set[bytes] = set()
        with Path(bgzf_path).open("rb") as source:
            for begin, chunk_end in merged:
                for line in _iter_virtual_lines(source, begin, chunk_end, budget=budget):
                    if not line or line.startswith(b"#") or (not preserve_duplicates and line in seen):
                        continue
                    fields = line.split(b"\t", 3)
                    if len(fields) < 3:
                        raise TabixError("indexed VCF record is malformed")
                    try:
                        record_contig = fields[0].decode("ascii")
                        position = int(fields[1])
                    except (UnicodeDecodeError, ValueError) as error:
                        raise TabixError("indexed VCF coordinates are malformed") from error
                    if record_contig == contig and start <= position <= end:
                        # Merged virtual-offset chunks visit each occurrence once.
                        # Source reconciliation needs identical physical rows too.
                        if not preserve_duplicates:
                            seen.add(line)
                        records.append(line.decode("utf-8"))
                        if len(records) >= max_records:
                            return tuple(records)
        return tuple(records)


def _region_bins(begin: int, end: int) -> tuple[int, ...]:
    if begin < 0 or end <= begin:
        raise ValueError("region must be non-empty and zero-based")
    end -= 1
    bins = [0]
    for start_bin, shift in ((1, 26), (9, 23), (73, 20), (585, 17), (4681, 14)):
        bins.extend(range(start_bin + (begin >> shift), start_bin + (end >> shift) + 1))
    return tuple(bins)


def _read_bgzf_block(source: BinaryIO, compressed_offset: int, *, budget=None) -> tuple[bytes, int]:
    if budget is not None:
        budget.check()
    source.seek(compressed_offset)
    fixed = source.read(12)
    if len(fixed) != 12 or fixed[:3] != b"\x1f\x8b\x08" or not fixed[3] & 4:
        raise TabixError(f"invalid BGZF block at compressed offset {compressed_offset}")
    extra_length = struct.unpack("<H", fixed[10:12])[0]
    extra = source.read(extra_length)
    if len(extra) != extra_length:
        raise TabixError("truncated BGZF extra data")
    block_size: int | None = None
    offset = 0
    while offset < len(extra):
        if offset + 4 > len(extra):
            raise TabixError("malformed BGZF extra field")
        subfield_length = struct.unpack("<H", extra[offset + 2 : offset + 4])[0]
        subfield_end = offset + 4 + subfield_length
        if subfield_end > len(extra):
            raise TabixError("truncated BGZF extra subfield")
        if extra[offset : offset + 2] == b"BC" and subfield_length == 2:
            block_size = struct.unpack("<H", extra[offset + 4 : subfield_end])[0] + 1
        offset = subfield_end
    if block_size is None or block_size < 12 + extra_length + 8:
        raise TabixError("BGZF block lacks a valid BC size field")
    remaining = block_size - 12 - extra_length
    suffix = source.read(remaining)
    if len(suffix) != remaining:
        raise TabixError("truncated BGZF block")
    block = fixed + extra + suffix
    try:
        if budget is not None:
            budget.charge(len(block), "tabix_transport")
        payload = _bounded_gunzip(block, 65536, budget, "tabix_blocks")
    except gzip.BadGzipFile as error:
        raise TabixError("BGZF block checksum or compression is invalid") from error
    return payload, compressed_offset + block_size


def _iter_virtual_lines(
    source: BinaryIO, begin_virtual: int, end_virtual: int, *, budget=None
):
    begin_block, begin_within = divmod(begin_virtual, 1 << 16)
    end_block, end_within = divmod(end_virtual, 1 << 16)
    block_offset = begin_block
    within = begin_within
    pending = b""
    while block_offset < end_block or (
        block_offset == end_block and within < end_within
    ):
        payload, next_block = _read_bgzf_block(source, block_offset, budget=budget)
        stop = end_within if block_offset == end_block else len(payload)
        if within > len(payload) or stop > len(payload) or stop < within:
            raise TabixError("virtual offset exceeds BGZF block")
        pending += payload[within:stop]
        if budget is not None and len(pending) > budget.max_row_bytes:
            raise TabixError("tabix record exceeds the read limit")
        while b"\n" in pending:
            line, pending = pending.split(b"\n", 1)
            yield line
        if next_block <= block_offset:
            raise TabixError("BGZF block offsets do not advance")
        block_offset = next_block
        within = 0


def _bounded_gunzip(content, maximum, budget, representation):
    if budget is None:
        return gzip.decompress(content)
    with gzip.GzipFile(fileobj=io.BytesIO(content)) as stream:
        chunks, size = [], 0
        while True:
            budget.check()
            chunk = stream.read(min(65536, budget.max_row_bytes, maximum - size + 1))
            if not chunk:
                return b"".join(chunks)
            size += len(chunk)
            if size > maximum:
                raise TabixError("tabix decompression exceeds the read limit")
            budget.charge(len(chunk), representation)
            chunks.append(chunk)
