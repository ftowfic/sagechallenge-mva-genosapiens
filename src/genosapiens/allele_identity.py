"""Immutable class-native allele identities for the controlled evidence model.

Coordinates are 1-based and intervals are closed. These records validate structure,
not reference concordance, normalization, caller validity, or disease mechanism.
Regulatory candidates use SequenceAllele with separately sourced downstream links.
Serialization is explicit: its output remains controlled, despite redacted reprs.
"""

from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass, field, fields
from enum import StrEnum
import math
import re
from typing import TypeAlias, TypeVar, cast


class AlleleIdentityError(ValueError):
    """Invalid identity or evidence, without controlled values in diagnostics."""


class ReferenceRole(StrEnum):
    CHALLENGE_ALIGNMENT = "challenge_alignment"
    CHALLENGE_TOOL_VIEW = "challenge_tool_view"
    PUBLIC_ANNOTATION_BENCHMARK = "public_annotation_benchmark"


class AlleleKind(StrEnum):
    SEQUENCE = "sequence"
    STRUCTURAL = "structural"
    MITOCHONDRIAL = "mitochondrial"
    MOSAIC = "mosaic"
    REPEAT = "repeat"


class StructuralVariantType(StrEnum):
    DEL = "DEL"
    DUP = "DUP"
    INV = "INV"
    INS = "INS"
    BND = "BND"
    CNV = "CNV"
    CPX = "CPX"
    TRA = "TRA"


class Orientation(StrEnum):
    FORWARD = "+"
    REVERSE = "-"


class RepeatSizeUnit(StrEnum):
    REPEAT_UNITS = "repeat-units"
    BASE_PAIRS = "base-pairs"


def _text(value: object) -> bool:
    return (
        type(value) is str
        and bool(value)
        and all(33 <= ord(character) <= 126 for character in value)
    )


def _contig(value: object) -> bool:
    return _text(value) and not any(character in value for character in ",[]")


def _integer(value: object, minimum: int = 0) -> bool:
    return type(value) is int and value >= minimum


def _number(value: object) -> bool:
    return type(value) is int or type(value) is float and math.isfinite(value)


def _sequence(value: object) -> bool:
    return type(value) is str and re.fullmatch(r"[ACGTRYSWKMBDHVN]+", value) is not None


def _reference(role: object, identity: object) -> None:
    if (
        type(role) is not ReferenceRole
        or type(identity) is not str
        or re.fullmatch(r"[a-z0-9][a-z0-9._-]*", identity) is None
    ):
        raise AlleleIdentityError("reference binding is invalid")


def _genotype(source_genotype: object, ploidy: object) -> None:
    if ploidy is not None and not _integer(ploidy, 1):
        raise AlleleIdentityError("ploidy is invalid")
    if source_genotype is None:
        return
    if (
        type(source_genotype) is not str
        or re.fullmatch(r"(?:[0-9]+|\.)(?:[/|](?:[0-9]+|\.))*", source_genotype)
        is None
    ):
        raise AlleleIdentityError("source genotype is invalid")
    # A lone missing GT does not establish ploidy. Never fill or rewrite GT.
    if (
        source_genotype != "."
        and ploidy is not None
        and len(re.split(r"[/|]", source_genotype)) != ploidy
    ):
        raise AlleleIdentityError("source genotype and ploidy disagree")


@dataclass(frozen=True, slots=True, repr=False)
class GenomicInterval:
    """A 1-based closed interval in the owning allele's reference binding."""

    contig: str
    start: int
    end: int

    def __post_init__(self) -> None:
        if (
            not _contig(self.contig)
            or not _integer(self.start, 1)
            or not _integer(self.end, 1)
            or self.end < self.start
        ):
            raise AlleleIdentityError("genomic interval is invalid")


@dataclass(frozen=True, slots=True, repr=False)
class Breakpoint:
    """Source breakpoint anchor; uncertainty bounds are signed position offsets.

    Orientation describes the source's reference strand; null means unknown.
    Breakpoint tuple order is source order, not inferred adjacency or phase.
    """

    contig: str
    position: int
    orientation: Orientation | None = None
    confidence_interval: tuple[int, int] | None = None

    def __post_init__(self) -> None:
        if (
            not _contig(self.contig)
            or not _integer(self.position, 1)
            or self.orientation is not None and type(self.orientation) is not Orientation
        ):
            raise AlleleIdentityError("breakpoint is invalid")
        bounds = self.confidence_interval
        if bounds is not None and (
            type(bounds) is not tuple
            or len(bounds) != 2
            or any(type(value) is not int for value in bounds)
            or not bounds[0] <= 0 <= bounds[1]
            or self.position + bounds[0] < 1
        ):
            raise AlleleIdentityError("breakpoint confidence interval is invalid")


@dataclass(frozen=True, slots=True, repr=False)
class ReadSupport:
    """Source DP and full AD, including unknown entries and other alternates.

    allele_number is the source ALT index (1-based), not a rewritten GT index.
    AD and DP are independent source measurements.
    No fraction is calculated or imputed from these counts.
    """

    depth: int | None = None
    allelic_depths: tuple[int | None, ...] | None = None
    allele_number: int | None = None

    def __post_init__(self) -> None:
        if self.depth is not None and not _integer(self.depth):
            raise AlleleIdentityError("read depth is invalid")
        if self.allele_number is not None and not _integer(self.allele_number, 1):
            raise AlleleIdentityError("source allele number is invalid")
        counts = self.allelic_depths
        if counts is not None:
            if (
                type(counts) is not tuple
                or not counts
                or any(value is not None and not _integer(value) for value in counts)
                or self.allele_number is not None and self.allele_number >= len(counts)
            ):
                raise AlleleIdentityError("allelic read depths are invalid")


@dataclass(frozen=True, slots=True, repr=False)
class CopyNumberEvidence:
    """Measured dosage, not an inferred deletion, inheritance, or mechanism."""

    copy_number: float | None = None
    baseline_copy_number: float | None = None
    source_genotype: str | None = None
    ploidy: int | None = None

    def __post_init__(self) -> None:
        for value in (self.copy_number, self.baseline_copy_number):
            if value is not None and (not _number(value) or value < 0):
                raise AlleleIdentityError("copy number is invalid")
        _genotype(self.source_genotype, self.ploidy)


@dataclass(frozen=True, slots=True, repr=False)
class SequenceAllele:
    """Explicit sequence identity; neither normalization nor gene assignment.

    SNVs, MNVs, indels, and regulatory sequence changes share this identity.
    Symbolic ALT, breakends, reference-only calls, and implicit proxies do not.
    """

    reference_role: ReferenceRole
    reference_identity_id: str
    contig: str
    position: int
    reference: str
    alternate: str
    kind: AlleleKind = field(default=AlleleKind.SEQUENCE, init=False)

    def __post_init__(self) -> None:
        _reference(self.reference_role, self.reference_identity_id)
        if (
            not _contig(self.contig)
            or not _integer(self.position, 1)
            or not _sequence(self.reference)
            or not _sequence(self.alternate)
            or self.reference == self.alternate
        ):
            raise AlleleIdentityError("sequence allele is invalid")


@dataclass(frozen=True, slots=True, repr=False)
class StructuralAllele:
    """Native SV/CNV geometry and dosage; exact_sequence is only nested evidence.

    A single breakend and unknown orientation/dosage are preserved as such.
    Intervals and breakpoint anchors have distinct meanings and are not silently
    interconverted. An exact_sequence does not establish SV equivalence.
    """

    reference_role: ReferenceRole
    reference_identity_id: str
    sv_type: StructuralVariantType
    interval: GenomicInterval | None = None
    breakpoints: tuple[Breakpoint, ...] = ()
    dosage: CopyNumberEvidence | None = None
    exact_sequence: SequenceAllele | None = None
    kind: AlleleKind = field(default=AlleleKind.STRUCTURAL, init=False)

    def __post_init__(self) -> None:
        _reference(self.reference_role, self.reference_identity_id)
        if (
            type(self.sv_type) is not StructuralVariantType
            or self.interval is not None and type(self.interval) is not GenomicInterval
            or type(self.breakpoints) is not tuple
            or any(type(value) is not Breakpoint for value in self.breakpoints)
            or self.interval is None and not self.breakpoints
            or self.dosage is not None and type(self.dosage) is not CopyNumberEvidence
        ):
            raise AlleleIdentityError("structural allele is invalid")
        _exact_sequence(self.exact_sequence, self.reference_role, self.reference_identity_id)
        if self.exact_sequence is not None:
            contigs = {breakpoint.contig for breakpoint in self.breakpoints}
            if self.interval is not None:
                contigs.add(self.interval.contig)
            if self.exact_sequence.contig not in contigs:
                raise AlleleIdentityError("structural and nested sequence contigs disagree")


def _fraction_evidence(
    sequence: object, fraction: object, support: object,
    source_genotype: object, ploidy: object,
) -> None:
    if (
        type(sequence) is not SequenceAllele
        or fraction is not None and (not _number(fraction) or not 0 <= fraction <= 1)
        or type(support) is not ReadSupport
    ):
        raise AlleleIdentityError("sequence fraction evidence is invalid")
    _genotype(source_genotype, ploidy)


@dataclass(frozen=True, slots=True, repr=False)
class MitochondrialAllele:
    """Caller-designated mitochondrial allele with nullable measured heteroplasmy.

    The caller/reference binding, not a contig-name heuristic or a diploid GT,
    establishes class context. Measured fractions need not equal raw AD/DP.
    """

    sequence: SequenceAllele
    heteroplasmy: float | None = None
    support: ReadSupport = field(default_factory=ReadSupport)
    source_genotype: str | None = None
    ploidy: int | None = None
    kind: AlleleKind = field(default=AlleleKind.MITOCHONDRIAL, init=False)

    def __post_init__(self) -> None:
        _fraction_evidence(self.sequence, self.heteroplasmy, self.support, self.source_genotype, self.ploidy)


@dataclass(frozen=True, slots=True, repr=False)
class MosaicAllele:
    """Mosaic candidate evidence; does not assert confirmed mosaicism or rewrite GT."""

    sequence: SequenceAllele
    allele_fraction: float | None = None
    support: ReadSupport = field(default_factory=ReadSupport)
    source_genotype: str | None = None
    ploidy: int | None = None
    kind: AlleleKind = field(default=AlleleKind.MOSAIC, init=False)

    def __post_init__(self) -> None:
        _fraction_evidence(self.sequence, self.allele_fraction, self.support, self.source_genotype, self.ploidy)


@dataclass(frozen=True, slots=True, repr=False)
class RepeatReadSupport:
    """Caller-reported read categories; no assumption that categories are disjoint."""

    spanning_reads: int | None = None
    flanking_reads: int | None = None
    in_repeat_reads: int | None = None

    def __post_init__(self) -> None:
        if any(
            value is not None and not _integer(value)
            for value in (self.spanning_reads, self.flanking_reads, self.in_repeat_reads)
        ):
            raise AlleleIdentityError("repeat read support is invalid")


@dataclass(frozen=True, slots=True, repr=False)
class RepeatSize:
    """One source-ordered size measurement and optional bounds in an explicit unit.

    Null size can retain a no-call or bounds-only estimate. Fractional repeat
    units are allowed; base-pair measurements must be integral. No diagnostic
    expansion threshold, repeat count conversion, or allele sorting is applied.
    """

    size: float | None
    unit: RepeatSizeUnit
    lower_bound: float | None = None
    upper_bound: float | None = None
    support: RepeatReadSupport = field(default_factory=RepeatReadSupport)

    def __post_init__(self) -> None:
        if type(self.unit) is not RepeatSizeUnit or type(self.support) is not RepeatReadSupport:
            raise AlleleIdentityError("repeat size evidence is invalid")
        for value in (self.size, self.lower_bound, self.upper_bound):
            if value is not None and (
                not _number(value)
                or value < 0
                or self.unit is RepeatSizeUnit.BASE_PAIRS and value != int(value)
            ):
                raise AlleleIdentityError("repeat size is invalid")
        if (
            self.lower_bound is not None and self.upper_bound is not None
            and self.lower_bound > self.upper_bound
            or self.size is not None and self.lower_bound is not None and self.size < self.lower_bound
            or self.size is not None and self.upper_bound is not None and self.size > self.upper_bound
        ):
            raise AlleleIdentityError("repeat size bounds are inconsistent")


def _exact_sequence(sequence: object, role: ReferenceRole, identity: str) -> None:
    if sequence is not None and (
        type(sequence) is not SequenceAllele
        or sequence.reference_role is not role
        or sequence.reference_identity_id != identity
    ):
        raise AlleleIdentityError("nested sequence reference binding is invalid")


@dataclass(frozen=True, slots=True, repr=False)
class RepeatAllele:
    """Native repeat payload, not proof of an individual allele or copy count.

    Source-bound interpretation distinguishes a retained locus call from any
    supported individual view. Measurement cardinality is independent of GT
    ploidy; partial measurements never require inventing or rewriting a genotype.
    """

    reference_role: ReferenceRole
    reference_identity_id: str
    locus_id: str
    interval: GenomicInterval
    repeat_unit: str
    allele_sizes: tuple[RepeatSize, ...] | None
    source_genotype: str | None = None
    ploidy: int | None = None
    exact_sequence: SequenceAllele | None = None
    kind: AlleleKind = field(default=AlleleKind.REPEAT, init=False)

    def __post_init__(self) -> None:
        _reference(self.reference_role, self.reference_identity_id)
        if (
            not _text(self.locus_id)
            or type(self.interval) is not GenomicInterval
            or not _sequence(self.repeat_unit)
            or self.allele_sizes is not None and (
                type(self.allele_sizes) is not tuple
                or not self.allele_sizes
                or any(type(value) is not RepeatSize for value in self.allele_sizes)
            )
        ):
            raise AlleleIdentityError("repeat allele is invalid")
        _genotype(self.source_genotype, self.ploidy)
        _exact_sequence(self.exact_sequence, self.reference_role, self.reference_identity_id)
        if self.exact_sequence is not None and self.exact_sequence.contig != self.interval.contig:
            raise AlleleIdentityError("repeat and nested sequence contigs disagree")


AlleleIdentity: TypeAlias = (
    SequenceAllele | StructuralAllele | MitochondrialAllele | MosaicAllele | RepeatAllele
)
_Record = TypeVar("_Record")

_IDENTITIES = {
    AlleleKind.SEQUENCE: SequenceAllele,
    AlleleKind.STRUCTURAL: StructuralAllele,
    AlleleKind.MITOCHONDRIAL: MitochondrialAllele,
    AlleleKind.MOSAIC: MosaicAllele,
    AlleleKind.REPEAT: RepeatAllele,
}
_RECORDS = (
    *_IDENTITIES.values(), GenomicInterval, Breakpoint, ReadSupport,
    CopyNumberEvidence, RepeatReadSupport, RepeatSize,
)
_NESTED = {
    (StructuralAllele, "interval"): GenomicInterval,
    (StructuralAllele, "dosage"): CopyNumberEvidence,
    (StructuralAllele, "exact_sequence"): SequenceAllele,
    (MitochondrialAllele, "sequence"): SequenceAllele,
    (MitochondrialAllele, "support"): ReadSupport,
    (MosaicAllele, "sequence"): SequenceAllele,
    (MosaicAllele, "support"): ReadSupport,
    (RepeatSize, "support"): RepeatReadSupport,
    (RepeatAllele, "interval"): GenomicInterval,
    (RepeatAllele, "exact_sequence"): SequenceAllele,
}
_ARRAYS = {
    (StructuralAllele, "breakpoints"): Breakpoint,
    (RepeatAllele, "allele_sizes"): RepeatSize,
    (Breakpoint, "confidence_interval"): None,
    (ReadSupport, "allelic_depths"): None,
}
_ENUMS = {
    "reference_role": ReferenceRole,
    "sv_type": StructuralVariantType,
    "orientation": Orientation,
    "unit": RepeatSizeUnit,
}


def _encode(value: object) -> object:
    if type(value) in _RECORDS:
        return {item.name: _encode(getattr(value, item.name)) for item in fields(value)}
    if isinstance(value, StrEnum):
        return value.value
    if type(value) is tuple:
        return [_encode(item) for item in value]
    return value


def identity_to_dict(identity: AlleleIdentity) -> dict[str, object]:
    """Return a detached JSON-ready mapping, retaining every native field and null.

    Only the five concrete identity types are accepted, never subclasses or a
    duck-typed sequence view of another class. No submission conversion occurs.
    """
    if type(identity) not in _IDENTITIES.values():
        raise AlleleIdentityError("typed allele identity is required")
    return cast(dict[str, object], _encode(identity))


def _record_from_dict(cls: type[_Record], data: object) -> _Record:
    if not isinstance(data, Mapping) or any(type(key) is not str for key in data):
        raise AlleleIdentityError("identity record must be an object with string keys")
    if set(data) != {item.name for item in fields(cls)}:
        raise AlleleIdentityError("identity record has missing or unexpected fields")
    arguments = {}
    for item in fields(cls):
        value = data[item.name]
        if not item.init:
            if type(value) is not str or value != item.default.value:
                raise AlleleIdentityError("allele kind does not match native identity")
            continue
        key = (cls, item.name)
        if value is not None:
            if key in _NESTED:
                value = _record_from_dict(_NESTED[key], value)
            elif key in _ARRAYS:
                if type(value) not in (list, tuple):
                    raise AlleleIdentityError("identity collection must be an array")
                child = _ARRAYS[key]
                value = tuple(_record_from_dict(child, entry) if child else entry for entry in value)
            elif item.name in _ENUMS:
                enum = _ENUMS[item.name]
                if type(value) is not str or value not in {member.value for member in enum}:
                    raise AlleleIdentityError("identity enum value is invalid")
                value = enum(value)
        arguments[item.name] = value
    return cls(**arguments)


def identity_from_dict(data: Mapping[str, object]) -> AlleleIdentity:
    """Decode a kind-tagged mapping without inferring fields, defaults, or classes.

    All serialized fields, even nullable/defaulted ones, must be present. Nested
    sequence evidence must itself carry kind=sequence. Irrelevant fields fail
    closed rather than being ignored or preserved in an untyped extension bag.
    """
    if not isinstance(data, Mapping):
        raise AlleleIdentityError("allele identity must be an object")
    kind = data.get("kind")
    if type(kind) is not str or kind not in _IDENTITIES:
        raise AlleleIdentityError("allele identity kind is invalid")
    return _record_from_dict(_IDENTITIES[kind], data)
