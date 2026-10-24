"""Lossless, reference-aware normalization of controlled small variants."""

from __future__ import annotations

from collections.abc import Callable, Iterable, Mapping
from contextlib import contextmanager
from dataclasses import asdict, dataclass
from enum import StrEnum
import gzip
import hashlib
import json
import os
from pathlib import Path
import re
from typing import BinaryIO, Iterator
import zlib

from genosapiens.governance import ArtifactClassification
from genosapiens.llm_guard import LLMPolicy
from genosapiens.reference_concordance import (
    PublicValidationConcordanceReport,
    ReferenceConcordanceError,
    ReferenceConcordanceReport,
    concordance_report_as_json,
    load_reference_identity,
    open_installed_reference,
)


VARIANT_NORMALIZATION_SCHEMA = "genosapiens.variant-normalization"
SOURCE_RECORD_NORMALIZATION_SCHEMA = "genosapiens.source-record-normalization"
SHA256 = re.compile(r"^[0-9a-f]{64}$")
IDENTIFIER = re.compile(r"^[A-Za-z0-9][A-Za-z0-9._:-]*$")
REFERENCE_ALLELE = re.compile(r"^[ACGTRYSWKMBDHVN]+$")
ALTERNATE_ALLELE = re.compile(r"^[ACGTN]+$")
SOURCE_RECORD_ID = re.compile(r"^source-record-[0-9a-f]{24}$")
SOURCE_ALLELE_ID = re.compile(r"^source-allele-[0-9a-f]{24}$")
NORMALIZED_ALLELE_ID = re.compile(r"^allele-[0-9a-f]{24}$")

ReferenceFetcher = Callable[[str, int, int], str]


class VariantNormalizationError(RuntimeError):
    """Raised without controlled variant content when normalization fails."""


class SmallVariantKind(StrEnum):
    SNV = "snv"
    MNV = "mnv"
    INSERTION = "insertion"
    DELETION = "deletion"
    DELINS = "delins"


class QuarantineReason(StrEnum):
    REFERENCE_MISMATCH = "reference-mismatch"
    REFERENCE_UNRESOLVABLE = "reference-unresolvable"
    UNSUPPORTED_ALTERNATE = "unsupported-alternate"
    NON_VARIANT_ALTERNATE = "non-variant-alternate"


class _ReferenceUnresolvable(RuntimeError):
    pass


def _matches(pattern: re.Pattern[str], value: object) -> bool:
    return isinstance(value, str) and pattern.fullmatch(value) is not None


def _valid_contig(value: str) -> bool:
    return (
        isinstance(value, str)
        and bool(value)
        and "," not in value
        and all(not character.isspace() and ord(character) >= 33 for character in value)
    )


def _canonical_bytes(value: object) -> bytes:
    return (
        json.dumps(value, sort_keys=True, separators=(",", ":"), ensure_ascii=True)
        + "\n"
    ).encode()


@dataclass(frozen=True, slots=True)
class SourceSmallVariantRecord:
    record_number: int
    contig: str
    position: int
    reference: str
    alternates: tuple[str, ...]

    def __post_init__(self) -> None:
        if (
            not isinstance(self.record_number, int)
            or isinstance(self.record_number, bool)
            or self.record_number <= 0
            or not _valid_contig(self.contig)
            or not isinstance(self.position, int)
            or isinstance(self.position, bool)
            or self.position <= 0
            or not _matches(REFERENCE_ALLELE, self.reference)
            or not isinstance(self.alternates, tuple)
            or not self.alternates
            or any(
                not isinstance(alternate, str)
                or not alternate
                or "," in alternate
                or any(character.isspace() for character in alternate)
                for alternate in self.alternates
            )
        ):
            raise VariantNormalizationError("source small-variant record is invalid")


def source_record_id(input_sha256: str, record: SourceSmallVariantRecord) -> str:
    semantic = {"input_sha256": input_sha256, "record": asdict(record)}
    return f"source-record-{hashlib.sha256(_canonical_bytes(semantic)).hexdigest()[:24]}"


def source_allele_id(source_record_id: str, allele_number: int) -> str:
    digest = hashlib.sha256(
        f"{source_record_id}\0{allele_number}".encode("ascii")
    ).hexdigest()[:24]
    return f"source-allele-{digest}"


@dataclass(frozen=True, slots=True)
class SourceAllele:
    source_allele_id: str
    source_record_id: str
    input_sha256: str
    record: SourceSmallVariantRecord
    allele_number: int
    alternate: str

    def __post_init__(self) -> None:
        if (
            not _matches(SOURCE_ALLELE_ID, self.source_allele_id)
            or not _matches(SOURCE_RECORD_ID, self.source_record_id)
            or not _matches(SHA256, self.input_sha256)
            or not isinstance(self.allele_number, int)
            or isinstance(self.allele_number, bool)
            or self.allele_number <= 0
            or self.allele_number > len(self.record.alternates)
            or self.alternate != self.record.alternates[self.allele_number - 1]
            or self.source_record_id
            != source_record_id(self.input_sha256, self.record)
            or self.source_allele_id
            != source_allele_id(self.source_record_id, self.allele_number)
        ):
            raise VariantNormalizationError("source allele provenance is invalid")


def _normalized_allele_id(
    reference_identity_id: str,
    reference_role: str,
    contig: str,
    position: int,
    reference: str,
    alternate: str,
) -> str:
    semantic = {
        "reference_identity_id": reference_identity_id,
        "reference_role": reference_role,
        "contig": contig,
        "position": position,
        "reference": reference,
        "alternate": alternate,
    }
    return f"allele-{hashlib.sha256(_canonical_bytes(semantic)).hexdigest()[:24]}"


def _variant_kind(reference: str, alternate: str) -> SmallVariantKind:
    if len(reference) == len(alternate):
        return SmallVariantKind.SNV if len(reference) == 1 else SmallVariantKind.MNV
    if len(reference) < len(alternate) and (
        alternate.startswith(reference) or alternate.endswith(reference)
    ):
        return SmallVariantKind.INSERTION
    if len(reference) > len(alternate) and (
        reference.startswith(alternate) or reference.endswith(alternate)
    ):
        return SmallVariantKind.DELETION
    return SmallVariantKind.DELINS


@dataclass(frozen=True, slots=True)
class NormalizedAllele:
    normalized_allele_id: str
    reference_identity_id: str
    reference_role: str
    contig: str
    position: int
    reference: str
    alternate: str
    kind: SmallVariantKind

    def __post_init__(self) -> None:
        if (
            not _matches(NORMALIZED_ALLELE_ID, self.normalized_allele_id)
            or not _matches(IDENTIFIER, self.reference_identity_id)
            or not _matches(IDENTIFIER, self.reference_role)
            or not _valid_contig(self.contig)
            or not isinstance(self.position, int)
            or isinstance(self.position, bool)
            or self.position <= 0
            or not _matches(REFERENCE_ALLELE, self.reference)
            or not _matches(REFERENCE_ALLELE, self.alternate)
            or self.reference == self.alternate
            or self.kind is not _variant_kind(self.reference, self.alternate)
            or self.normalized_allele_id
            != _normalized_allele_id(
                self.reference_identity_id,
                self.reference_role,
                self.contig,
                self.position,
                self.reference,
                self.alternate,
            )
        ):
            raise VariantNormalizationError("normalized allele is invalid")


@dataclass(frozen=True, slots=True)
class SourceToNormalizedMapping:
    source: SourceAllele
    normalized_allele_id: str
    reference_context_start: int
    reference_context: str

    def __post_init__(self) -> None:
        if (
            not _matches(NORMALIZED_ALLELE_ID, self.normalized_allele_id)
            or not _matches(ALTERNATE_ALLELE, self.source.alternate)
            or self.source.alternate == self.source.record.reference
            or not isinstance(self.reference_context_start, int)
            or isinstance(self.reference_context_start, bool)
            or self.reference_context_start <= 0
            or not _matches(REFERENCE_ALLELE, self.reference_context)
        ):
            raise VariantNormalizationError("source-to-normalized mapping is invalid")


@dataclass(frozen=True, slots=True)
class QuarantinedSourceAllele:
    source: SourceAllele
    reason: QuarantineReason

    def __post_init__(self) -> None:
        if not isinstance(self.reason, QuarantineReason):
            raise VariantNormalizationError("source allele quarantine is invalid")


@dataclass(frozen=True, slots=True, repr=False)
class SourceAlleleNormalization:
    """One physical ALT result without repeating its complete source record."""

    allele_number: int
    source_allele_id: str
    alternate: str
    normalized_allele: NormalizedAllele | None
    reference_context_start: int | None
    reference_context: str | None
    quarantine_reason: QuarantineReason | None

    def __post_init__(self) -> None:
        normalized = self.normalized_allele is not None
        if (
            not isinstance(self.allele_number, int)
            or isinstance(self.allele_number, bool)
            or self.allele_number <= 0
            or not _matches(SOURCE_ALLELE_ID, self.source_allele_id)
            or not isinstance(self.alternate, str)
            or not self.alternate
            or "," in self.alternate
            or any(character.isspace() for character in self.alternate)
            or normalized == (self.quarantine_reason is not None)
            or normalized != (self.reference_context_start is not None)
            or normalized != (self.reference_context is not None)
            or normalized and (
                not isinstance(self.reference_context_start, int)
                or isinstance(self.reference_context_start, bool)
                or self.reference_context_start <= 0
                or not _matches(REFERENCE_ALLELE, self.reference_context)
            )
            or not normalized and not isinstance(
                self.quarantine_reason, QuarantineReason
            )
        ):
            raise VariantNormalizationError("source allele normalization is invalid")


@dataclass(frozen=True, slots=True, repr=False)
class SourceRecordNormalization:
    """Bounded normalization facts for one original source record."""

    input_sha256: str
    reference_identity_id: str
    reference_role: str
    source_record_id: str
    record: SourceSmallVariantRecord
    alleles: tuple[SourceAlleleNormalization, ...]
    classification: ArtifactClassification = ArtifactClassification.CONTROLLED
    llm_policy: LLMPolicy = LLMPolicy.NO_LLM
    schema: str = SOURCE_RECORD_NORMALIZATION_SCHEMA

    def __post_init__(self) -> None:
        if (
            self.schema != SOURCE_RECORD_NORMALIZATION_SCHEMA
            or not _matches(SHA256, self.input_sha256)
            or not _matches(IDENTIFIER, self.reference_identity_id)
            or not _matches(IDENTIFIER, self.reference_role)
            or not _matches(SOURCE_RECORD_ID, self.source_record_id)
            or self.source_record_id != source_record_id(self.input_sha256, self.record)
            or not isinstance(self.alleles, tuple)
            or len(self.alleles) != len(self.record.alternates)
            or self.classification is not ArtifactClassification.CONTROLLED
            or self.llm_policy is not LLMPolicy.NO_LLM
        ):
            raise VariantNormalizationError("source record normalization is invalid")
        for allele_number, (alternate, outcome) in enumerate(
            zip(self.record.alternates, self.alleles, strict=True), start=1
        ):
            expected_source = self.source_allele(allele_number)
            if (
                outcome.allele_number != allele_number
                or outcome.source_allele_id != expected_source.source_allele_id
                or outcome.alternate != alternate
            ):
                raise VariantNormalizationError("source record normalization is invalid")
            if outcome.normalized_allele is not None:
                normalized = outcome.normalized_allele
                mapping = SourceToNormalizedMapping(
                    expected_source,
                    normalized.normalized_allele_id,
                    outcome.reference_context_start,
                    outcome.reference_context,
                )
                if (
                    normalized.reference_identity_id != self.reference_identity_id
                    or normalized.reference_role != self.reference_role
                    or not _mapping_is_semantically_consistent(mapping, normalized)
                ):
                    raise VariantNormalizationError(
                        "source record normalization is invalid"
                    )
            elif (
                outcome.quarantine_reason is QuarantineReason.UNSUPPORTED_ALTERNATE
                and _matches(ALTERNATE_ALLELE, outcome.alternate)
            ) or (
                outcome.quarantine_reason is QuarantineReason.NON_VARIANT_ALTERNATE
                and outcome.alternate != self.record.reference
            ):
                raise VariantNormalizationError("source record normalization is invalid")

    def source_allele(self, allele_number: int) -> SourceAllele:
        if (
            not isinstance(allele_number, int)
            or isinstance(allele_number, bool)
            or not 1 <= allele_number <= len(self.record.alternates)
        ):
            raise VariantNormalizationError("source allele ordinal is invalid")
        return SourceAllele(
            source_allele_id=source_allele_id(self.source_record_id, allele_number),
            source_record_id=self.source_record_id,
            input_sha256=self.input_sha256,
            record=self.record,
            allele_number=allele_number,
            alternate=self.record.alternates[allele_number - 1],
        )

    @property
    def source_to_normalized(self) -> tuple[SourceToNormalizedMapping, ...]:
        return tuple(
            SourceToNormalizedMapping(
                self.source_allele(outcome.allele_number),
                outcome.normalized_allele.normalized_allele_id,
                outcome.reference_context_start,
                outcome.reference_context,
            )
            for outcome in self.alleles
            if outcome.normalized_allele is not None
        )

    @property
    def quarantined_source_alleles(self) -> tuple[QuarantinedSourceAllele, ...]:
        return tuple(
            QuarantinedSourceAllele(
                self.source_allele(outcome.allele_number),
                outcome.quarantine_reason,
            )
            for outcome in self.alleles
            if outcome.quarantine_reason is not None
        )


def source_record_normalization_as_json(result: SourceRecordNormalization) -> str:
    if not isinstance(result, SourceRecordNormalization):
        raise VariantNormalizationError("source record normalization is invalid")
    try:
        return json.dumps(
            asdict(result), sort_keys=True, separators=(",", ":"),
            ensure_ascii=True, allow_nan=False,
        ) + "\n"
    except (TypeError, ValueError, UnicodeError, RecursionError, OverflowError):
        raise VariantNormalizationError(
            "cannot encode source record normalization"
        ) from None


@dataclass(frozen=True, slots=True)
class SmallVariantNormalizationResult:
    input_sha256: str
    reference_identity_id: str
    reference_role: str
    reference_manifest_sha256: str
    reference_concordance_sha256: str
    total_source_records: int
    total_source_alleles: int
    normalized_alleles: tuple[NormalizedAllele, ...]
    source_to_normalized: tuple[SourceToNormalizedMapping, ...]
    quarantined_source_alleles: tuple[QuarantinedSourceAllele, ...]
    reference_mismatch_record_ids: tuple[str, ...]
    classification: ArtifactClassification
    llm_policy: LLMPolicy
    schema: str = VARIANT_NORMALIZATION_SCHEMA

    def __post_init__(self) -> None:
        normalized_ids = tuple(
            allele.normalized_allele_id for allele in self.normalized_alleles
        )
        normalized_id_set = set(normalized_ids)
        mapped_ids = tuple(
            mapping.normalized_allele_id for mapping in self.source_to_normalized
        )
        sources = tuple(
            mapping.source for mapping in self.source_to_normalized
        ) + tuple(item.source for item in self.quarantined_source_alleles)
        source_allele_ids = tuple(source.source_allele_id for source in sources)
        records_by_id: dict[str, SourceSmallVariantRecord] = {}
        for source in sources:
            records_by_id.setdefault(source.source_record_id, source.record)
        expected_source_allele_ids = {
            source_allele_id(source_record_id, allele_number)
            for source_record_id, record in records_by_id.items()
            for allele_number in range(1, len(record.alternates) + 1)
        }
        normalized_by_id = {
            allele.normalized_allele_id: allele for allele in self.normalized_alleles
        }
        expected_normalized_order = tuple(dict.fromkeys(mapped_ids))
        expected_mismatch_records = tuple(
            dict.fromkeys(
                item.source.source_record_id
                for item in self.quarantined_source_alleles
                if item.reason is QuarantineReason.REFERENCE_MISMATCH
            )
        )
        if (
            self.schema != VARIANT_NORMALIZATION_SCHEMA
            or not _matches(SHA256, self.input_sha256)
            or not _matches(IDENTIFIER, self.reference_identity_id)
            or not _matches(IDENTIFIER, self.reference_role)
            or not _matches(SHA256, self.reference_manifest_sha256)
            or not _matches(SHA256, self.reference_concordance_sha256)
            or not isinstance(self.total_source_records, int)
            or isinstance(self.total_source_records, bool)
            or self.total_source_records <= 0
            or not isinstance(self.total_source_alleles, int)
            or isinstance(self.total_source_alleles, bool)
            or self.total_source_alleles <= 0
            or self.total_source_records != len(records_by_id)
            or self.total_source_alleles != len(sources)
            or len(source_allele_ids) != len(set(source_allele_ids))
            or set(source_allele_ids) != expected_source_allele_ids
            or len({record.record_number for record in records_by_id.values()})
            != len(records_by_id)
            or len(normalized_ids) != len(set(normalized_ids))
            or expected_normalized_order != normalized_ids
            or any(identifier not in normalized_id_set for identifier in mapped_ids)
            or any(
                not _mapping_is_semantically_consistent(
                    mapping,
                    normalized_by_id[mapping.normalized_allele_id],
                )
                for mapping in self.source_to_normalized
            )
            or self.reference_mismatch_record_ids != expected_mismatch_records
            or any(source.input_sha256 != self.input_sha256 for source in sources)
            or any(
                allele.reference_identity_id != self.reference_identity_id
                or allele.reference_role != self.reference_role
                for allele in self.normalized_alleles
            )
            or any(
                item.reason is QuarantineReason.UNSUPPORTED_ALTERNATE
                and _matches(ALTERNATE_ALLELE, item.source.alternate)
                or item.reason is QuarantineReason.NON_VARIANT_ALTERNATE
                and item.source.alternate != item.source.record.reference
                for item in self.quarantined_source_alleles
            )
            or any(
                source.source_record_id in expected_mismatch_records
                and not any(
                    item.source.source_allele_id == source.source_allele_id
                    and item.reason is QuarantineReason.REFERENCE_MISMATCH
                    for item in self.quarantined_source_alleles
                )
                for source in sources
            )
            or self.classification is not ArtifactClassification.CONTROLLED
            or self.llm_policy is not LLMPolicy.NO_LLM
        ):
            raise VariantNormalizationError("small-variant normalization result is invalid")


def _fetch_reference(
    fetch_reference: ReferenceFetcher,
    contig: str,
    position: int,
    length: int,
) -> str:
    try:
        sequence = fetch_reference(contig, position, length)
    except (IndexError, KeyError, OSError, TypeError, ValueError) as error:
        raise _ReferenceUnresolvable from error
    if (
        not isinstance(sequence, str)
        or len(sequence) != length
        or not _matches(REFERENCE_ALLELE, sequence)
    ):
        raise _ReferenceUnresolvable
    return sequence


def _normalize_sequence_allele(
    source: SourceAllele,
    *,
    reference_identity_id: str,
    reference_role: str,
    fetch_reference: ReferenceFetcher,
) -> NormalizedAllele:
    position = source.record.position
    reference = source.record.reference
    alternate = source.alternate

    # Right trimming plus upstream extension rotates indels to their leftmost
    # equivalent representation while retaining a valid VCF anchor base.
    while reference[-1] == alternate[-1]:
        if len(reference) > 1 and len(alternate) > 1:
            reference = reference[:-1]
            alternate = alternate[:-1]
        elif position > 1:
            preceding = _fetch_reference(
                fetch_reference, source.record.contig, position - 1, 1
            )
            reference = preceding + reference
            alternate = preceding + alternate
            position -= 1
        else:
            break

    while (
        len(reference) > 1
        and len(alternate) > 1
        and reference[0] == alternate[0]
    ):
        reference = reference[1:]
        alternate = alternate[1:]
        position += 1

    if (
        _fetch_reference(
            fetch_reference, source.record.contig, position, len(reference)
        )
        != reference
    ):
        raise VariantNormalizationError(
            "reference changed during small-variant normalization"
        )
    normalized_id = _normalized_allele_id(
        reference_identity_id,
        reference_role,
        source.record.contig,
        position,
        reference,
        alternate,
    )
    return NormalizedAllele(
        normalized_allele_id=normalized_id,
        reference_identity_id=reference_identity_id,
        reference_role=reference_role,
        contig=source.record.contig,
        position=position,
        reference=reference,
        alternate=alternate,
        kind=_variant_kind(reference, alternate),
    )


def _mapping_is_semantically_consistent(
    mapping: SourceToNormalizedMapping, normalized: NormalizedAllele
) -> bool:
    source = mapping.source
    if source.record.contig != normalized.contig:
        return False
    context_start = mapping.reference_context_start
    context_end = context_start + len(mapping.reference_context) - 1
    source_end = source.record.position + len(source.record.reference) - 1
    normalized_end = normalized.position + len(normalized.reference) - 1
    if (
        source.record.position < context_start
        or source_end > context_end
        or normalized.position < context_start
        or normalized_end > context_end
    ):
        return False

    def apply(position: int, reference: str, alternate: str) -> str | None:
        offset = position - context_start
        if mapping.reference_context[offset : offset + len(reference)] != reference:
            return None
        return (
            mapping.reference_context[:offset]
            + alternate
            + mapping.reference_context[offset + len(reference) :]
        )

    source_haplotype = apply(
        source.record.position,
        source.record.reference,
        source.alternate,
    )
    normalized_haplotype = apply(
        normalized.position,
        normalized.reference,
        normalized.alternate,
    )
    return source_haplotype is not None and source_haplotype == normalized_haplotype


def _validate_normalization_binding(
    *, input_sha256: str, reference_identity_id: str,
    reference_role: str, fetch_reference: ReferenceFetcher,
) -> None:
    if (
        not _matches(SHA256, input_sha256)
        or not _matches(IDENTIFIER, reference_identity_id)
        or not _matches(IDENTIFIER, reference_role)
        or not callable(fetch_reference)
    ):
        raise VariantNormalizationError("small-variant normalization binding is invalid")


def normalize_source_record(
    record: SourceSmallVariantRecord,
    *,
    input_sha256: str,
    reference_identity_id: str,
    reference_role: str,
    fetch_reference: ReferenceFetcher,
) -> SourceRecordNormalization:
    """Normalize every physical ALT in one source record with bounded state."""

    _validate_normalization_binding(
        input_sha256=input_sha256,
        reference_identity_id=reference_identity_id,
        reference_role=reference_role,
        fetch_reference=fetch_reference,
    )
    if not isinstance(record, SourceSmallVariantRecord):
        raise VariantNormalizationError("source small-variant record is invalid")
    record_id = source_record_id(input_sha256, record)
    try:
        observed_reference = _fetch_reference(
            fetch_reference, record.contig, record.position, len(record.reference)
        )
    except _ReferenceUnresolvable:
        record_quarantine = QuarantineReason.REFERENCE_UNRESOLVABLE
    else:
        record_quarantine = (
            QuarantineReason.REFERENCE_MISMATCH
            if observed_reference != record.reference
            else None
        )

    outcomes = []
    for allele_number, alternate in enumerate(record.alternates, start=1):
        allele_id = source_allele_id(record_id, allele_number)
        source = SourceAllele(
            source_allele_id=allele_id,
            source_record_id=record_id,
            input_sha256=input_sha256,
            record=record,
            allele_number=allele_number,
            alternate=alternate,
        )
        reason = record_quarantine
        normalized = None
        context_start = None
        reference_context = None
        if reason is None and not _matches(ALTERNATE_ALLELE, alternate):
            reason = QuarantineReason.UNSUPPORTED_ALTERNATE
        elif reason is None and alternate == record.reference:
            reason = QuarantineReason.NON_VARIANT_ALTERNATE
        elif reason is None:
            try:
                normalized = _normalize_sequence_allele(
                    source,
                    reference_identity_id=reference_identity_id,
                    reference_role=reference_role,
                    fetch_reference=fetch_reference,
                )
                context_start = min(record.position, normalized.position)
                context_end = max(
                    record.position + len(record.reference) - 1,
                    normalized.position + len(normalized.reference) - 1,
                )
                reference_context = _fetch_reference(
                    fetch_reference,
                    record.contig,
                    context_start,
                    context_end - context_start + 1,
                )
            except _ReferenceUnresolvable:
                reason = QuarantineReason.REFERENCE_UNRESOLVABLE
                normalized = None
                context_start = None
                reference_context = None
        outcomes.append(SourceAlleleNormalization(
            allele_number=allele_number,
            source_allele_id=allele_id,
            alternate=alternate,
            normalized_allele=normalized,
            reference_context_start=context_start,
            reference_context=reference_context,
            quarantine_reason=reason,
        ))
    return SourceRecordNormalization(
        input_sha256=input_sha256,
        reference_identity_id=reference_identity_id,
        reference_role=reference_role,
        source_record_id=record_id,
        record=record,
        alleles=tuple(outcomes),
    )


def iter_source_record_normalizations(
    records: Iterable[SourceSmallVariantRecord],
    *,
    input_sha256: str,
    reference_identity_id: str,
    reference_role: str,
    fetch_reference: ReferenceFetcher,
) -> Iterator[SourceRecordNormalization]:
    """Yield ordered per-record facts without retaining prior source records."""

    _validate_normalization_binding(
        input_sha256=input_sha256,
        reference_identity_id=reference_identity_id,
        reference_role=reference_role,
        fetch_reference=fetch_reference,
    )
    previous_record_number = 0
    found = False
    for record in records:
        if (
            not isinstance(record, SourceSmallVariantRecord)
            or record.record_number <= previous_record_number
        ):
            raise VariantNormalizationError("source record ordering is invalid")
        found = True
        previous_record_number = record.record_number
        yield normalize_source_record(
            record,
            input_sha256=input_sha256,
            reference_identity_id=reference_identity_id,
            reference_role=reference_role,
            fetch_reference=fetch_reference,
        )
    if not found:
        raise VariantNormalizationError("source record ordering is invalid")


def _normalize_small_variant_records(
    records: Iterable[SourceSmallVariantRecord],
    *,
    input_sha256: str,
    reference_identity_id: str,
    reference_role: str,
    reference_manifest_sha256: str,
    reference_concordance_sha256: str,
    fetch_reference: ReferenceFetcher,
) -> SmallVariantNormalizationResult:
    """Split and normalize source records without discarding source provenance."""

    if (
        not _matches(SHA256, input_sha256)
        or not _matches(IDENTIFIER, reference_identity_id)
        or not _matches(IDENTIFIER, reference_role)
        or not _matches(SHA256, reference_manifest_sha256)
        or not _matches(SHA256, reference_concordance_sha256)
        or not callable(fetch_reference)
    ):
        raise VariantNormalizationError("small-variant normalization binding is invalid")
    normalized_by_id: dict[str, NormalizedAllele] = {}
    mappings: list[SourceToNormalizedMapping] = []
    quarantined: list[QuarantinedSourceAllele] = []
    total_source_records = 0
    total_source_alleles = 0
    for result in iter_source_record_normalizations(
        records,
        input_sha256=input_sha256,
        reference_identity_id=reference_identity_id,
        reference_role=reference_role,
        fetch_reference=fetch_reference,
    ):
        total_source_records += 1
        total_source_alleles += len(result.alleles)
        for outcome in result.alleles:
            if outcome.normalized_allele is not None:
                normalized_by_id.setdefault(
                    outcome.normalized_allele.normalized_allele_id,
                    outcome.normalized_allele,
                )
        mappings.extend(result.source_to_normalized)
        quarantined.extend(result.quarantined_source_alleles)

    mismatch_records = tuple(
        dict.fromkeys(
            item.source.source_record_id
            for item in quarantined
            if item.reason is QuarantineReason.REFERENCE_MISMATCH
        )
    )
    return SmallVariantNormalizationResult(
        input_sha256=input_sha256,
        reference_identity_id=reference_identity_id,
        reference_role=reference_role,
        reference_manifest_sha256=reference_manifest_sha256,
        reference_concordance_sha256=reference_concordance_sha256,
        total_source_records=total_source_records,
        total_source_alleles=total_source_alleles,
        normalized_alleles=tuple(normalized_by_id.values()),
        source_to_normalized=tuple(mappings),
        quarantined_source_alleles=tuple(quarantined),
        reference_mismatch_record_ids=mismatch_records,
        classification=ArtifactClassification.CONTROLLED,
        llm_policy=LLMPolicy.NO_LLM,
    )


def _stat_identity(value: os.stat_result) -> tuple[int, int, int, int, int]:
    return (value.st_dev, value.st_ino, value.st_size, value.st_mtime_ns, value.st_ctime_ns)


def _file_sha256(path: Path) -> tuple[str, tuple[int, int, int, int, int]]:
    digest = hashlib.sha256()
    try:
        initial = _stat_identity(path.stat())
        with path.open("rb") as source:
            while chunk := source.read(1024 * 1024):
                digest.update(chunk)
            final_descriptor = _stat_identity(os.fstat(source.fileno()))
        final_path = _stat_identity(path.stat())
    except OSError as error:
        raise VariantNormalizationError("cannot read normalization input") from error
    if initial != final_descriptor or initial != final_path:
        raise VariantNormalizationError("normalization input changed while being read")
    return digest.hexdigest(), initial


@contextmanager
def _open_vcf_records(path: Path) -> Iterator[BinaryIO]:
    try:
        with path.open("rb") as probe:
            compressed = probe.read(2) == b"\x1f\x8b"
        source = gzip.open(path, "rb") if compressed else path.open("rb")
        with source:
            yield source
    except (EOFError, gzip.BadGzipFile, OSError, zlib.error) as error:
        raise VariantNormalizationError("cannot parse normalization input") from error


def _read_vcf_small_variant_records(path: Path) -> tuple[SourceSmallVariantRecord, ...]:
    records: list[SourceSmallVariantRecord] = []
    with _open_vcf_records(path) as source:
        for raw_line in source:
            if raw_line.startswith(b"#"):
                continue
            try:
                fields = raw_line.rstrip(b"\r\n").decode("ascii").split("\t")
                if len(fields) < 8:
                    raise ValueError
                position = int(fields[1])
                alternates = tuple(fields[4].split(","))
                records.append(
                    SourceSmallVariantRecord(
                        record_number=len(records) + 1,
                        contig=fields[0],
                        position=position,
                        reference=fields[3],
                        alternates=alternates,
                    )
                )
            except (UnicodeDecodeError, ValueError, VariantNormalizationError) as error:
                raise VariantNormalizationError(
                    "normalization input contains an invalid VCF record"
                ) from error
    return tuple(records)


def _normalize_bound_vcf_small_variants(
    vcf_path: str | Path,
    *,
    concordance_report: (
        ReferenceConcordanceReport | PublicValidationConcordanceReport
    ),
    fetch_reference: ReferenceFetcher,
    accepted_reference_roles: tuple[str, ...] = ("challenge_alignment",),
) -> SmallVariantNormalizationResult:
    """Normalize the exact VCF and reference bound by a concordance report."""

    if (
        not isinstance(
            concordance_report,
            (ReferenceConcordanceReport, PublicValidationConcordanceReport),
        )
        or concordance_report.reference_role not in accepted_reference_roles
        or not concordance_report.header_valid
        or not concordance_report.reference_basename_match
        or not concordance_report.contig_names_match
        or not concordance_report.contig_lengths_match
        or not concordance_report.contig_order_match
        or concordance_report.unresolved_records != 0
    ):
        raise VariantNormalizationError("reference concordance binding is invalid")

    path = Path(vcf_path)
    input_sha256, input_identity = _file_sha256(path)
    if input_sha256 != concordance_report.validated_input_sha256:
        raise VariantNormalizationError("normalization input binding is stale")
    records = _read_vcf_small_variant_records(path)
    try:
        final_identity = _stat_identity(path.stat())
    except OSError as error:
        raise VariantNormalizationError("cannot verify normalization input") from error
    if input_identity != final_identity:
        raise VariantNormalizationError("normalization input changed while being read")
    if len(records) != concordance_report.total_records:
        raise VariantNormalizationError("normalization record accounting is inconsistent")

    report_sha256 = hashlib.sha256(
        concordance_report_as_json(concordance_report).encode("ascii")
    ).hexdigest()
    result = _normalize_small_variant_records(
        records,
        input_sha256=input_sha256,
        reference_identity_id=concordance_report.reference_identity_id,
        reference_role=concordance_report.reference_role,
        reference_manifest_sha256=concordance_report.reference_manifest_sha256,
        reference_concordance_sha256=report_sha256,
        fetch_reference=fetch_reference,
    )
    if (
        len(result.reference_mismatch_record_ids)
        != concordance_report.mismatching_records
        or any(
            item.reason is QuarantineReason.REFERENCE_UNRESOLVABLE
            for item in result.quarantined_source_alleles
        )
    ):
        raise VariantNormalizationError(
            "normalization and reference concordance results disagree"
        )
    return result


def normalize_public_validation_vcf_small_variants(
    vcf_path: str | Path,
    *,
    concordance_report: PublicValidationConcordanceReport,
    fetch_reference: ReferenceFetcher,
    source_to_reference_contigs: Mapping[str, str],
) -> SmallVariantNormalizationResult:
    """Normalize a public fixture without granting challenge-reference status."""
    aliases = dict(source_to_reference_contigs)
    alias_sha256 = hashlib.sha256(
        (
            json.dumps(aliases, indent=2, sort_keys=True, ensure_ascii=True) + "\n"
        ).encode("ascii")
    ).hexdigest()
    if (
        type(concordance_report) is not PublicValidationConcordanceReport
        or not aliases
        or any(
            type(source) is not str
            or not source
            or type(target) is not str
            or not target
            for source, target in aliases.items()
        )
        or alias_sha256 != concordance_report.alias_map_sha256
    ):
        raise VariantNormalizationError(
            "public validation concordance binding is invalid"
        )

    def mapped_fetch(contig: str, position: int, length: int) -> str:
        return fetch_reference(aliases[contig], position, length)

    return _normalize_bound_vcf_small_variants(
        vcf_path,
        concordance_report=concordance_report,
        fetch_reference=mapped_fetch,
        accepted_reference_roles=("public_annotation_benchmark",),
    )


def normalize_vcf_small_variants(
    vcf_path: str | Path,
    *,
    concordance_report: ReferenceConcordanceReport,
    identity_path: str | Path,
    reference_manifest_path: str | Path,
    reference_root: str | Path,
) -> SmallVariantNormalizationResult:
    """Normalize a VCF using its checksum-validated installed reference bundle."""

    if not isinstance(concordance_report, ReferenceConcordanceReport):
        raise VariantNormalizationError("reference concordance binding is invalid")
    try:
        contract = load_reference_identity(identity_path)
        role = contract.reference_role("challenge_alignment")
        with open_installed_reference(
            contract,
            identity_path=identity_path,
            reference_manifest_path=reference_manifest_path,
            reference_root=reference_root,
        ) as (manifest_sha256, fetch_reference, _reference_contigs):
            if (
                manifest_sha256 != concordance_report.reference_manifest_sha256
                or role.identity_id != concordance_report.reference_identity_id
                or role.role != concordance_report.reference_role
            ):
                raise VariantNormalizationError(
                    "installed reference differs from concordance binding"
                )
            return _normalize_bound_vcf_small_variants(
                vcf_path,
                concordance_report=concordance_report,
                fetch_reference=fetch_reference,
            )
    except ReferenceConcordanceError as error:
        raise VariantNormalizationError(
            "installed reference binding is invalid"
        ) from error


def variant_normalization_as_json(result: SmallVariantNormalizationResult) -> str:
    return json.dumps(asdict(result), indent=2, sort_keys=True) + "\n"
