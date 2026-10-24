"""Lossless packed tape for already evaluated complete pair outcomes.

Each position is the corresponding ordinal in CompletePairDomain. Integers name
the owning evaluator's authenticated outcome catalog; this codec neither makes
a scientific decision nor substitutes inputs for unevaluated outcomes. Every
position is present, including pairs that construct no hypothesis.
"""
from __future__ import annotations

from array import array
import base64
from collections.abc import Mapping
from dataclasses import dataclass
import hashlib
import sys
import zlib

from genosapiens.wgs_evidence_store import EvidenceReader, EvidenceReference, EvidenceStoreError
from genosapiens.wgs_pair_domain import MAX_COMPLETE_PAIRS, PairDomainSummary

SCHEMA = "genosapiens.wgs-pair-outcome-tape"
VERSION = 1
CHUNK_PAIRS = 32768
MAX_DECODED_BYTES = CHUNK_PAIRS * 4
MAX_COMPRESSED_BYTES = MAX_DECODED_BYTES + 1024
MAX_CATALOG_ENTRIES = 1_000_000
CODEC = "zlib-u32-little-endian-v1"


def _require(condition, message="invalid complete pair outcome tape", *, capacity=False):
    if not condition:
        raise EvidenceStoreError(message, capacity=capacity)


def _digest(value):
    return type(value) is str and len(value) == 64 and all(c in "0123456789abcdef" for c in value)


@dataclass(frozen=True, slots=True, repr=False)
class PairOutcomeTape:
    domain_sha256: str
    catalog_sha256: str
    catalog_count: int
    pair_count: int
    chunks_name: str
    chunks_sha256: str
    stream_sha256: str
    schema: str = SCHEMA
    version: int = VERSION

    def __post_init__(self):
        _require(self.schema == SCHEMA and type(self.version) is int and self.version == VERSION
            and all(_digest(value) for value in (self.domain_sha256, self.catalog_sha256,
                self.chunks_sha256, self.stream_sha256))
            and type(self.catalog_count) is int and 0 <= self.catalog_count <= MAX_CATALOG_ENTRIES
            and type(self.pair_count) is int and 0 <= self.pair_count <= MAX_COMPLETE_PAIRS
            and bool(self.catalog_count) == bool(self.pair_count)
            and type(self.chunks_name) is str and 0 < len(self.chunks_name) <= 512)


def write_pair_outcomes(writer, name, outcomes, *, domain, catalog_sha256, catalog_count):
    """Exhaust all evaluated outcomes before returning a publishable descriptor."""
    _require(type(domain) is PairDomainSummary and _digest(catalog_sha256)
        and type(catalog_count) is int and 0 <= catalog_count <= MAX_CATALOG_ENTRIES
        and bool(catalog_count) == bool(domain.pair_count))
    stream = hashlib.sha256()
    total = 0
    used = bytearray(catalog_count)
    def chunks():
        nonlocal total
        values = array("I")
        _require(values.itemsize == 4)
        def record():
            first = total - len(values) + 1
            if sys.byteorder != "little":
                values.byteswap()
            content = values.tobytes()
            stream.update(content)
            compressed = zlib.compress(content, level=6)
            _require(len(compressed) <= MAX_COMPRESSED_BYTES, capacity=True)
            return {"first_pair": first, "pair_count": len(values), "codec": CODEC,
                "decoded_sha256": hashlib.sha256(content).hexdigest(),
                "payload": base64.b64encode(compressed).decode("ascii")}
        for outcome in outcomes:
            _require(type(outcome) is int and 0 <= outcome < catalog_count)
            total += 1
            _require(total <= domain.pair_count, "extra pair outcomes")
            used[outcome] = 1
            values.append(outcome)
            if len(values) == CHUNK_PAIRS:
                yield record()
                values = array("I")
        if values:
            yield record()
        _require(total == domain.pair_count, "complete pair outcome coverage is missing")
        _require(all(used), "unreferenced pair outcome catalog entries")
    chunks_name = name + ".chunks"
    digest = writer.add_sequence(chunks_name, "pair-outcome-chunks-v1", chunks())
    tape = PairOutcomeTape(domain.membership_sha256, catalog_sha256, catalog_count,
        total, chunks_name, digest, stream.hexdigest())
    writer.add(name, "pair-outcome-tape-v1", tape)
    return tape


def _decode(chunk, *, first, remaining, catalog_count):
    _require(isinstance(chunk, Mapping) and set(chunk) == {
        "first_pair", "pair_count", "codec", "decoded_sha256", "payload"})
    count = min(CHUNK_PAIRS, remaining)
    _require(type(chunk["first_pair"]) is int and chunk["first_pair"] == first
        and type(chunk["pair_count"]) is int and chunk["pair_count"] == count and count > 0
        and chunk["codec"] == CODEC and _digest(chunk["decoded_sha256"])
        and type(chunk["payload"]) is str
        and len(chunk["payload"]) <= 4 * ((MAX_COMPRESSED_BYTES + 2) // 3))
    try:
        encoded = base64.b64decode(chunk["payload"], validate=True)
        _require(base64.b64encode(encoded).decode("ascii") == chunk["payload"])
        decoder = zlib.decompressobj()
        content = decoder.decompress(encoded, count * 4 + 1)
        # Never flush an unbounded decoder, and reject concatenated/trailing
        # streams even if the first stream's expansion and digest look valid.
        _require(len(content) == count * 4 and decoder.eof
            and not decoder.unconsumed_tail and not decoder.unused_data
            and hashlib.sha256(content).hexdigest() == chunk["decoded_sha256"])
        values = array("I")
        _require(values.itemsize == 4)
        values.frombytes(content)
        if sys.byteorder != "little":
            values.byteswap()
        _require(all(value < catalog_count for value in values))
        return content, values
    except (ValueError, TypeError, zlib.error):
        raise EvidenceStoreError("invalid compressed pair outcome chunk") from None


class VerifiedPairOutcomes:
    """Authenticate all chunks and catalog coverage before yielding an outcome."""
    def __init__(self, reader, reference, *, domain, catalog_sha256, catalog_count):
        _require(type(reader) is EvidenceReader and type(reference) is EvidenceReference
            and reference.kind == "pair-outcome-tape-v1" and type(domain) is PairDomainSummary)
        self.reader, self.reference = reader, reference
        try:
            self.tape = PairOutcomeTape(**reader.reconstruct(reference, max_bytes=8192))
        except (TypeError, ValueError):
            raise EvidenceStoreError() from None
        tape = self.tape
        _require((tape.domain_sha256, tape.catalog_sha256, tape.catalog_count, tape.pair_count)
            == (domain.membership_sha256, catalog_sha256, catalog_count, domain.pair_count),
            "foreign pair domain or outcome catalog")
        self.chunks = reader.reference(tape.chunks_name, kind="pair-outcome-chunks-v1")
        _require(self.chunks.sha256 == tape.chunks_sha256)
        self.counts = array("Q", [0]) * catalog_count
        total, digest = 0, hashlib.sha256()
        for content, values in self._chunks():
            total += len(values)
            digest.update(content)
            for value in values:
                self.counts[value] += 1
        _require(total == tape.pair_count and digest.hexdigest() == tape.stream_sha256
            and all(self.counts), "incomplete pair outcome closure")

    def _chunks(self):
        total = 0
        for chunk in self.reader.iter_array(self.chunks, max_item_bytes=2 * MAX_COMPRESSED_BYTES):
            content, values = _decode(chunk, first=total + 1,
                remaining=self.tape.pair_count - total, catalog_count=self.tape.catalog_count)
            total += len(values)
            yield content, values

    def pages(self):
        total, digest = 0, hashlib.sha256()
        for content, values in self._chunks():
            total += len(values)
            digest.update(content)
            yield values
        _require(total == self.tape.pair_count and digest.hexdigest() == self.tape.stream_sha256)

    def outcomes(self):
        for values in self.pages():
            yield from values
