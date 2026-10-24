"""Deterministic genotype and call-quality summaries for controlled VCFs."""

from __future__ import annotations

from collections import Counter
from contextlib import contextmanager
from dataclasses import asdict, dataclass
from enum import StrEnum
from decimal import Decimal, InvalidOperation
import gzip
import hashlib
import json
import math
import os
from pathlib import Path
import re
import sqlite3
from typing import BinaryIO, Iterator
import zlib

from genosapiens.vcf_context import prepare_header, FormatLayoutCache
from genosapiens.vcf_contracts import (
    DEFAULT_MAX_HEADER_BYTES, DEFAULT_MAX_RECORD_BYTES, POLICY_SHA256,
    PidPhaseContract, VcfReadLimits, require_current_policy, bounded_line,
    FieldState, PhaseState, FormatFieldEvidence, field_state, expected_cardinality,
)
from genosapiens.governance import ArtifactClassification
from genosapiens.llm_guard import LLMPolicy
from genosapiens.reference_concordance import (
    PublicValidationConcordanceReport,
    ReferenceConcordanceReport,
    concordance_report_as_json,
)
from genosapiens.variant_normalization import (
    ALTERNATE_ALLELE,
    REFERENCE_ALLELE,
    SHA256,
    SourceSmallVariantRecord,
    source_record_id,
)


VCF_QC_SCHEMA = "genosapiens.vcf-qc"
FORMAT_HEADER = re.compile(r"^##FORMAT=<(.+)>$")
CONTIG_HEADER = re.compile(r"^##contig=<(.+)>$")
GENOTYPE = re.compile(r"^(?:\.|[0-9]+)(?:[/|](?:\.|[0-9]+))*$")
FORMAT_CONTRACTS = {
    "GT": ("1", "String"),
    "GQ": ("1", "Integer"),
    "DP": ("1", "Integer"),
    "AD": ("R", "Integer"),
    "PGT": ("1", "String"),
    "PID": ("1", "String"),
    "PS": ("1", "Integer"),
}
CLASS_INFO_FIELDS = frozenset({
    "SVTYPE", "END", "REPID", "RU", "START", "MOTIF", "PERIOD",
    "INEXACT_ALLELE", "BPDIFFS",
})
CLASS_FORMAT_FIELDS = frozenset({"CN", "REPCN", "GB"})
MINIMUM_SUPPORTED_GQ = 20
MINIMUM_SUPPORTED_DP = 10
MINIMUM_HETEROZYGOUS_ALLELE_BALANCE = 0.2
MAXIMUM_HETEROZYGOUS_ALLELE_BALANCE = 0.8
MINIMUM_HOMOZYGOUS_OR_HAPLOID_ALT_BALANCE = 0.8
MINIMUM_POLYPLOID_DOSAGE_BALANCE_RATIO = 0.4
MAXIMUM_POLYPLOID_DOSAGE_BALANCE_RATIO = 1.6


class VcfQcError(RuntimeError):
    """Raised without controlled record content when VCF QC fails."""


class PhaseStatus(StrEnum):
    PHASED = "phased"
    UNPHASED = "unphased"
    PARTIALLY_PHASED = "partially-phased"
    NOT_APPLICABLE = "not-applicable"
    UNKNOWN = "unknown"


class PhaseEvidence(StrEnum):
    GT = "gt"
    PGT_PID = "pgt-pid"
    NONE = "none"


class FilterStatus(StrEnum):
    PASS = "pass"
    FILTERED = "filtered"
    NOT_APPLIED = "not-applied"


class VariantClass(StrEnum):
    SNV = "snv"
    MNV = "mnv"
    INSERTION = "insertion"
    DELETION = "deletion"
    DELINS = "delins"
    SYMBOLIC = "symbolic"
    BREAKEND = "breakend"
    SPANNING_DELETION = "spanning-deletion"
    NON_VARIANT = "non-variant"
    OTHER = "other"


def _valid_text(value: object) -> bool:
    return (
        isinstance(value, str)
        and bool(value)
        and "\t" not in value
        and "\r" not in value
        and "\n" not in value
    )


def _valid_optional_count(value: object) -> bool:
    return value is None or (
        isinstance(value, int) and not isinstance(value, bool) and value >= 0
    )


def _valid_count(value: object) -> bool:
    return isinstance(value, int) and not isinstance(value, bool) and value >= 0


def _valid_class_fields(
    fields: tuple[tuple[str, str | None], ...], allowed: frozenset[str]
) -> bool:
    return (
        isinstance(fields, tuple)
        and all(
            isinstance(item, tuple)
            and len(item) == 2
            and item[0] in allowed
            and (item[1] is None or _valid_text(item[1]))
            for item in fields
        )
        and tuple(key for key, _value in fields)
        == tuple(sorted({key for key, _value in fields}))
    )


def _structured_attributes(value: str) -> dict[str, str] | None:
    attributes: dict[str, str] = {}
    parts: list[str] = []
    current: list[str] = []
    quoted = False
    escaped = False
    for character in value:
        if escaped:
            current.append(character)
            escaped = False
        elif character == "\\" and quoted:
            current.append(character)
            escaped = True
        elif character == '"':
            current.append(character)
            quoted = not quoted
        elif character == "," and not quoted:
            parts.append("".join(current))
            current = []
        else:
            current.append(character)
    if quoted or escaped:
        return None
    parts.append("".join(current))
    for part in parts:
        if "=" not in part:
            return None
        key, item = part.split("=", 1)
        if not key or key in attributes:
            return None
        if len(item) >= 2 and item[0] == item[-1] == '"':
            item = item[1:-1]
        attributes[key] = item
    return attributes


def _parse_contig_header(line: str) -> tuple[str, int]:
    match = CONTIG_HEADER.fullmatch(line)
    attributes = _structured_attributes(match.group(1)) if match is not None else None
    if attributes is None:
        raise VcfQcError("VCF QC contains an invalid contig declaration")
    contig = attributes.get("ID")
    length_text = attributes.get("length")
    try:
        length = int(length_text) if length_text is not None else 0
    except ValueError as error:
        raise VcfQcError("VCF QC contains an invalid contig declaration") from error
    if (
        not _valid_text(contig)
        or any(character.isspace() for character in contig)
        or length <= 0
    ):
        raise VcfQcError("VCF QC contains an invalid contig declaration")
    return contig, length


def _parse_format_header(line: str) -> str:
    match = FORMAT_HEADER.fullmatch(line)
    attributes = _structured_attributes(match.group(1)) if match is not None else None
    if attributes is None:
        raise VcfQcError("VCF QC contains an invalid FORMAT declaration")
    identifier = attributes.get("ID")
    number = attributes.get("Number")
    field_type = attributes.get("Type")
    if (
        not _valid_text(identifier)
        or any(character.isspace() for character in identifier)
        or not _valid_text(number)
        or not _valid_text(field_type)
        or identifier in FORMAT_CONTRACTS
        and (number, field_type) != FORMAT_CONTRACTS[identifier]
    ):
        raise VcfQcError("VCF QC contains an incompatible FORMAT declaration")
    return identifier


def _alternate_allele_count(alternates: tuple[str, ...]) -> int:
    if len(set(alternates)) != len(alternates):
        raise VcfQcError("VCF QC contains duplicate alternate alleles")
    if "." in alternates:
        if alternates == (".",):
            return 0
        raise VcfQcError("VCF QC combines a missing ALT with alternate alleles")
    return len(alternates)


def _gt_tokens(genotype: str | None) -> tuple[str, ...]:
    if genotype is None:
        return ()
    return tuple(re.split(r"[/|]", genotype))


def _ploidy(genotype: str | None) -> int | None:
    if genotype is None:
        return None
    return len(_gt_tokens(genotype))


def _phase(
    genotype: str | None,
    pgt: str | None,
    pid: str | None,
    pgt_pid_declared: bool,
    pid_contract_sha256: str | None = None,
) -> tuple[PhaseStatus, PhaseEvidence]:
    tokens = _gt_tokens(genotype)
    if genotype is None or not tokens or all(token == "." for token in tokens):
        return PhaseStatus.UNKNOWN, PhaseEvidence.NONE
    if "/" in genotype and "|" in genotype:
        return PhaseStatus.PARTIALLY_PHASED, PhaseEvidence.GT
    if "|" in genotype and any(token == "." for token in tokens):
        return PhaseStatus.PARTIALLY_PHASED, PhaseEvidence.GT
    if "|" in genotype:
        return PhaseStatus.PHASED, PhaseEvidence.GT
    if (
        pgt_pid_declared
        and pid_contract_sha256 is not None
        and pgt is not None
        and pid is not None
        and "|" in pgt
        and "/" not in pgt
        and GENOTYPE.fullmatch(pgt) is not None
        and all(token != "." for token in tokens)
        and all(token != "." for token in _gt_tokens(pgt))
        and sorted(tokens) == sorted(_gt_tokens(pgt))
    ):
        return PhaseStatus.PHASED, PhaseEvidence.PGT_PID
    if "/" in genotype:
        return PhaseStatus.UNPHASED, PhaseEvidence.GT
    return PhaseStatus.NOT_APPLICABLE, PhaseEvidence.GT


def _callability(genotype: str | None, depth: int | None) -> bool | None:
    tokens = _gt_tokens(genotype)
    if genotype is None or not tokens:
        return None
    if any(token == "." for token in tokens):
        return False
    if depth is None:
        return None
    return depth > 0


def _expected_allele_balance(
    allelic_depths: tuple[int | None, ...] | None,
    depth: int | None,
) -> tuple[float | None, ...] | None:
    if allelic_depths is None or depth is None or len(allelic_depths) < 2:
        return None
    if any(value is not None and not _valid_count(value) for value in allelic_depths):
        return None
    if depth == 0:
        return None
    return tuple(
        None if value is None else value / depth for value in allelic_depths[1:]
    )


@dataclass(frozen=True, slots=True)
class SampleCallQc:
    sample_name: str
    genotype: str | None
    ploidy: int | None
    phase_status: PhaseStatus
    phase_evidence: PhaseEvidence
    pgt_pid_declared: bool
    pgt: str | None
    pid: str | None
    phase_set: str | None
    genotype_quality: int | None
    depth: int | None
    allelic_depths: tuple[int | None, ...] | None
    alternate_allele_balance: tuple[float | None, ...] | None
    allele_balance_source: str | None
    callable: bool | None
    class_format_fields: tuple[tuple[str, str | None], ...]
    evidence_policy_sha256: str = POLICY_SHA256
    format_fields: tuple[FormatFieldEvidence, ...] = ()
    pid_contract_sha256: str | None = None

    @property
    def phase_state(self):
        return phase_state_for(self)[0]

    @property
    def phase_reasons(self):
        return phase_state_for(self)[1]

    def format_state(self, key):
        return next((field.state for field in self.format_fields if field.key == key), FieldState.ABSENT)

    @property
    def ad_dp_state(self):
        if self.allelic_depths is None or self.depth is None:
            return "unavailable"
        if sum(value for value in self.allelic_depths if value is not None) > self.depth:
            return "unresolved"
        if None in self.allelic_depths:
            return "partial"
        return "available"

    def __post_init__(self) -> None:
        require_current_policy(self.evidence_policy_sha256)
        if (type(self.format_fields) is not tuple
                or any(type(field) is not FormatFieldEvidence for field in self.format_fields)
                or len({field.key for field in self.format_fields}) != len(self.format_fields)):
            raise VcfQcError("invalid sample FORMAT evidence")
        expected_balance = _expected_allele_balance(self.allelic_depths, self.depth)
        expected_phase = _phase(
            self.genotype,
            self.pgt,
            self.pid,
            self.pgt_pid_declared, self.pid_contract_sha256,
        )
        if (
            not _valid_text(self.sample_name)
            or (
                self.genotype is not None
                and GENOTYPE.fullmatch(self.genotype) is None
            )
            or self.ploidy != _ploidy(self.genotype)
            or (self.phase_status, self.phase_evidence) != expected_phase
            or self.pid_contract_sha256 is not None and (
                type(self.pid_contract_sha256) is not str
                or SHA256.fullmatch(self.pid_contract_sha256) is None)
            or not isinstance(self.pgt_pid_declared, bool)
            or any(
                value is not None and not _valid_text(value)
                for value in (self.pgt, self.pid, self.phase_set)
            )
            or not _valid_optional_count(self.genotype_quality)
            or not _valid_optional_count(self.depth)
            or (
                self.allelic_depths is not None
                and (
                    not self.allelic_depths
                    or any(
                        value is not None and not _valid_count(value)
                        for value in self.allelic_depths
                    )
                )
            )
            or (
                self.alternate_allele_balance is not None
                and (
                    not self.alternate_allele_balance
                    or any(
                        value is not None
                        and (
                            not isinstance(value, float)
                            or not math.isfinite(value)
                            or value < 0.0
                        )
                        for value in self.alternate_allele_balance
                    )
                )
            )
            or self.alternate_allele_balance != expected_balance
            or (self.alternate_allele_balance is None)
            != (self.allele_balance_source is None)
            or self.allele_balance_source not in {None, "format-ad-dp"}
            or self.callable != _callability(self.genotype, self.depth)
            or not _valid_class_fields(self.class_format_fields, CLASS_FORMAT_FIELDS)
        ):
            raise VcfQcError("sample call QC record is invalid")


def phase_state_for(call):
    """Preserve uncertainty independently from usable standard GT evidence."""
    gt, pgt = _gt_tokens(call.genotype), _gt_tokens(call.pgt)
    if call.pgt is not None and gt and pgt and "." not in gt and "." not in pgt:
        if sorted(int(t) for t in gt) != sorted(int(t) for t in pgt):
            return PhaseState.CONFLICTING, ("pgt-gt-copy-conflict",)
    reasons = []
    if call.pgt is not None or call.pid is not None:
        if call.pgt is None or call.pid is None:
            reasons.append("pgt-pid-partner-missing")
        if not gt or "." in gt or not pgt or "." in pgt:
            reasons.append("pgt-cannot-complete-missing-genotype")
        if call.pgt is not None and ("|" not in call.pgt or "/" in call.pgt):
            reasons.append("pgt-orientation-incomplete")
        if reasons:
            return PhaseState.INCOMPLETE, tuple(reasons)
        if not call.pgt_pid_declared or call.pid_contract_sha256 is None:
            return PhaseState.UNRESOLVED, ("pgt-pid-producer-contract-unavailable",)
    if call.phase_set is not None and (not call.genotype or "/" in call.genotype or "|" not in call.genotype):
        return PhaseState.INCOMPLETE, ("ps-without-phased-gt",)
    if call.phase_status is PhaseStatus.PARTIALLY_PHASED:
        return PhaseState.INCOMPLETE, ("gt-phase-incomplete",)
    if call.phase_status is PhaseStatus.PHASED:
        return PhaseState.USABLE, ()
    return PhaseState.ABSENT, ()


def phase_reasons_for(call, allele_number):
    """Scope uncertainty to alleles present or unresolved in GT or PGT."""
    tokens = (*_gt_tokens(call.genotype), *_gt_tokens(call.pgt))
    if not tokens or "." in tokens or str(allele_number) in tokens:
        return call.phase_reasons
    return ()


def _called_alternate_without_read_support(
    call: SampleCallQc, allele_number: int
) -> bool:
    if call.genotype is None or call.allelic_depths is None:
        return False
    return (
        str(allele_number) in _gt_tokens(call.genotype)
        and allele_number < len(call.allelic_depths)
        and call.allelic_depths[allele_number] == 0
    )


def _heterozygous_allele_balance_out_of_range(
    call: SampleCallQc, allele_number: int
) -> bool:
    tokens = _gt_tokens(call.genotype)
    if (
        len(tokens) != 2
        or "." in tokens
        or len(set(tokens)) != 2
        or str(allele_number) not in tokens
        or call.alternate_allele_balance is None
        or allele_number > len(call.alternate_allele_balance)
    ):
        return False
    balance = call.alternate_allele_balance[allele_number - 1]
    return balance is not None and not (
        MINIMUM_HETEROZYGOUS_ALLELE_BALANCE
        <= balance
        <= MAXIMUM_HETEROZYGOUS_ALLELE_BALANCE
    )


def _homozygous_or_haploid_alt_balance_below_minimum(
    call: SampleCallQc, allele_number: int,
) -> bool:
    tokens = _gt_tokens(call.genotype)
    if (
        not tokens
        or "." in tokens
        or len(set(tokens)) != 1
        or tokens[0] != str(allele_number)
        or call.alternate_allele_balance is None
        or allele_number > len(call.alternate_allele_balance)
    ):
        return False
    balance = call.alternate_allele_balance[allele_number - 1]
    return (
        balance is not None
        and balance < MINIMUM_HOMOZYGOUS_OR_HAPLOID_ALT_BALANCE
    )


def _mixed_polyploid_allele_balance_out_of_range(
    call: SampleCallQc, allele_number: int
) -> bool:
    tokens = _gt_tokens(call.genotype)
    if (
        len(tokens) <= 2
        or "." in tokens
        or len(set(tokens)) == 1
        or str(allele_number) not in tokens
        or call.alternate_allele_balance is None
        or allele_number > len(call.alternate_allele_balance)
    ):
        return False
    allele_counts = Counter(token for token in tokens if token != "0")
    copies = allele_counts[str(allele_number)]
    balance = call.alternate_allele_balance[allele_number - 1]
    expected_balance = copies / len(tokens)
    return balance is not None and not (
        expected_balance * MINIMUM_POLYPLOID_DOSAGE_BALANCE_RATIO
        <= balance
        <= min(
            1.0,
            expected_balance * MAXIMUM_POLYPLOID_DOSAGE_BALANCE_RATIO,
        )
    )


def _technical_scopes(
    filter_status: FilterStatus,
    sample_calls: tuple[SampleCallQc, ...],
    alternate_slots: int,
) -> tuple[
    tuple[str, ...],
    tuple[tuple[str, ...], ...],
    tuple[tuple[str, ...], ...],
    tuple[tuple[tuple[str, ...], ...], ...],
]:
    record_flags = []
    if filter_status is FilterStatus.FILTERED:
        record_flags.append("source-filtered")
    elif filter_status is FilterStatus.NOT_APPLIED:
        record_flags.append("source-unfiltered")
    sample_flags = []
    sample_alternate_flags = []
    for call in sample_calls:
        flags = []
        if call.ad_dp_state == "unresolved":
            flags.append("ad-dp-relationship-unresolved")
        if call.callable is False:
            flags.append("sample-genotype-uncallable")
        if (
            call.genotype_quality is not None
            and call.genotype_quality < MINIMUM_SUPPORTED_GQ
        ):
            flags.append("genotype-quality-below-20")
        if call.depth is not None and call.depth < MINIMUM_SUPPORTED_DP:
            flags.append("depth-below-10")
        sample_flags.append(tuple(flags))

        alternate_flags = []
        for allele_number in range(1, alternate_slots + 1):
            flags = list(phase_reasons_for(call, allele_number))
            if _heterozygous_allele_balance_out_of_range(call, allele_number):
                flags.append("heterozygous-allele-balance-outside-0.2-0.8")
            if _homozygous_or_haploid_alt_balance_below_minimum(
                call, allele_number
            ):
                flags.append("homozygous-or-haploid-alt-balance-below-0.8")
            if _mixed_polyploid_allele_balance_out_of_range(call, allele_number):
                flags.append("polyploid-allele-balance-outside-0.4x-1.6x-dosage")
            if _called_alternate_without_read_support(call, allele_number):
                flags.append("called-alternate-without-read-support")
            alternate_flags.append(tuple(flags))
        sample_alternate_flags.append(tuple(alternate_flags))
    return (
        tuple(record_flags),
        tuple(() for _alternate in range(alternate_slots)),
        tuple(sample_flags),
        tuple(sample_alternate_flags),
    )


def _technical_scope_required(scopes: tuple[object, ...]) -> bool:
    record_flags, alternate_flags, sample_flags, sample_alternate_flags = scopes
    return bool(
        record_flags
        or any(alternate_flags)
        or any(sample_flags)
        or any(flags for sample in sample_alternate_flags for flags in sample)
    )


@dataclass(frozen=True, slots=True)
class VcfRecordQc:
    record_number: int
    source_record_id: str
    contig: str
    position: int
    reference: str
    alternates: tuple[str, ...]
    variant_classes: tuple[VariantClass, ...]
    filter_value: str
    filter_status: FilterStatus
    filters: tuple[str, ...]
    callset_allele_number: int | None
    callset_allele_frequencies: tuple[float | None, ...] | None
    sample_calls: tuple[SampleCallQc, ...]
    technical_penalty_required: bool
    read_validation_required: bool
    technical_flags: tuple[str, ...]
    alternate_technical_flags: tuple[tuple[str, ...], ...]
    sample_technical_flags: tuple[tuple[str, ...], ...]
    sample_alternate_technical_flags: tuple[tuple[tuple[str, ...], ...], ...]
    class_info_fields: tuple[tuple[str, str | None], ...]
    quality: float | None = None
    raw_line_sha256: str | None = None
    raw_line: str | None = None
    evidence_policy_sha256: str = POLICY_SHA256

    raw_quality: str | None = None

    def __post_init__(self) -> None:
        if self.raw_line is not None:
            raw_fields = self.raw_line.rstrip("\r\n").split("\t")
            if len(raw_fields) < 8 or self.raw_quality != raw_fields[5]:
                raise VcfQcError("VCF raw quality token is inconsistent")
            if self.quality != _parse_quality(self.raw_quality):
                raise VcfQcError("VCF quality conversion is inconsistent")
            if self.filter_value != raw_fields[6]:
                raise VcfQcError("VCF raw FILTER token is inconsistent")
        require_current_policy(self.evidence_policy_sha256)
        alternate_count = _alternate_allele_count(self.alternates)
        expected_scopes = _technical_scopes(
            self.filter_status, self.sample_calls, len(self.alternates)
        )
        if (
            not isinstance(self.record_number, int)
            or isinstance(self.record_number, bool)
            or self.record_number <= 0
            or not _valid_text(self.source_record_id)
            or not _valid_text(self.contig)
            or not isinstance(self.position, int)
            or isinstance(self.position, bool)
            or self.position <= 0
            or REFERENCE_ALLELE.fullmatch(self.reference) is None
            or not self.alternates
            or len(self.alternates) != len(self.variant_classes)
            or not _valid_text(self.filter_value)
            or (self.filter_status is FilterStatus.FILTERED) != bool(self.filters)
            or self.filter_status is FilterStatus.PASS
            and self.filter_value != "PASS"
            or self.filter_status is FilterStatus.NOT_APPLIED
            and self.filter_value != "."
            or not _valid_optional_count(self.callset_allele_number)
            or (
                self.callset_allele_frequencies is not None
                and (
                    not self.callset_allele_frequencies
                    or len(self.callset_allele_frequencies) != alternate_count
                    or any(
                        value is not None
                        and (
                            not isinstance(value, float)
                            or not math.isfinite(value)
                            or value < 0.0
                            or value > 1.0
                        )
                        for value in self.callset_allele_frequencies
                    )
                )
            )
            or len({call.sample_name for call in self.sample_calls})
            != len(self.sample_calls)
            or any(
                call.allelic_depths is not None
                and len(call.allelic_depths) != alternate_count + 1
                for call in self.sample_calls
            )
            or not isinstance(self.technical_penalty_required, bool)
            or not isinstance(self.read_validation_required, bool)
            or (
                self.technical_flags,
                self.alternate_technical_flags,
                self.sample_technical_flags,
                self.sample_alternate_technical_flags,
            ) != expected_scopes
            or self.technical_penalty_required
            != _technical_scope_required(expected_scopes)
            or self.read_validation_required
            != _technical_scope_required(expected_scopes)
            or not _valid_class_fields(self.class_info_fields, CLASS_INFO_FIELDS)
            or self.quality is not None
            and (not isinstance(self.quality, float) or not math.isfinite(self.quality) or self.quality < 0)
            or self.raw_line_sha256 is not None
            and SHA256.fullmatch(self.raw_line_sha256) is None
            or self.raw_line is not None and (
                type(self.raw_line) is not str or not self.raw_line.isascii()
                or hashlib.sha256(self.raw_line.encode("ascii")).hexdigest() != self.raw_line_sha256)
        ):
            raise VcfQcError("VCF record QC is invalid")


def technical_flags_for(
    record: VcfRecordQc,
    sample_call: SampleCallQc,
    allele_number: int | None,
) -> tuple[str, ...]:
    """Return only flags applicable to one sample and optional source ALT."""

    try:
        sample_index = record.sample_calls.index(sample_call)
    except ValueError:
        raise VcfQcError("sample call is not bound to the VCF record") from None
    if allele_number is not None and not 1 <= allele_number <= len(record.alternates):
        raise VcfQcError("source allele number is not bound to the VCF record")
    alternate_indices = (
        range(len(record.alternates))
        if allele_number is None
        else (allele_number - 1,)
    )
    return tuple(dict.fromkeys((
        *record.technical_flags,
        *(flag for index in alternate_indices for flag in record.alternate_technical_flags[index]),
        *record.sample_technical_flags[sample_index],
        *(flag for index in alternate_indices
          for flag in record.sample_alternate_technical_flags[sample_index][index]),
    )))


def quality_acceptance_for(record, sample_call, allele_number):
    """Describe exact occurrence quality without discharging technical review."""
    flags = technical_flags_for(record, sample_call, allele_number)
    if record.filter_value != "PASS":
        state, reasons = "failed", ("record-filter-not-PASS", *flags)
    elif sample_call.ad_dp_state == "unresolved":
        state, reasons = "unresolved", flags
    elif flags:
        state, reasons = "unresolved" if sample_call.phase_reasons else "failed", flags
    else:
        missing = []
        if sample_call.callable is None:
            missing.append("callability-unresolved")
        if sample_call.genotype_quality is None:
            missing.append("genotype-quality-unavailable")
        if sample_call.depth is None:
            missing.append("sample-depth-unavailable")
        if sample_call.alternate_allele_balance is None or (
                allele_number is not None and sample_call.alternate_allele_balance[allele_number - 1] is None):
            missing.append("allele-balance-unavailable")
        state, reasons = ("unresolved", tuple(missing)) if missing else ("accepted", ())
    return {"policy_sha256": POLICY_SHA256, "source_record_id": record.source_record_id,
            "sample_name": sample_call.sample_name, "allele_number": allele_number,
            "filter_value": record.filter_value, "filter_status": record.filter_status.value,
            "state": state, "reasons": tuple(dict.fromkeys(reasons))}


@dataclass(frozen=True, slots=True)
class DeclaredContig:
    contig: str
    length: int

    def __post_init__(self) -> None:
        if not _valid_text(self.contig) or not _valid_count(self.length) or self.length <= 0:
            raise VcfQcError("declared VCF contig is invalid")


@dataclass(frozen=True, slots=True)
class VcfQcStreamHeader:
    input_sha256: str
    input_size_bytes: int
    reference_concordance_sha256: str
    header_sha256: str
    header_bytes: int
    sample_names: tuple[str, ...]
    format_definitions: tuple[str, ...]
    phase_capability_fields: tuple[str, ...]
    declared_contigs: tuple[DeclaredContig, ...]
    class_field_definitions: tuple[tuple[str, str, str, str], ...]
    raw_lines: tuple[str, ...]
    evidence_policy_sha256: str = POLICY_SHA256
    read_limits: VcfReadLimits = VcfReadLimits()
    pid_phase_contract: PidPhaseContract | None = None

    def __post_init__(self) -> None:
        if (type(self.read_limits) is not VcfReadLimits
                or self.pid_phase_contract is not None and (
                    type(self.pid_phase_contract) is not PidPhaseContract
                    or self.pid_phase_contract.input_sha256 != self.input_sha256)):
            raise VcfQcError("VCF header policy binding is invalid")
        require_current_policy(self.evidence_policy_sha256)
        declared = {value.contig: value.length for value in self.declared_contigs}
        if (
            SHA256.fullmatch(self.input_sha256) is None
            or not isinstance(self.input_size_bytes, int)
            or isinstance(self.input_size_bytes, bool)
            or self.input_size_bytes <= 0
            or SHA256.fullmatch(self.reference_concordance_sha256) is None
            or SHA256.fullmatch(self.header_sha256) is None
            or not isinstance(self.header_bytes, int)
            or isinstance(self.header_bytes, bool)
            or self.header_bytes <= 0
            or not self.sample_names
            or len(set(self.sample_names)) != len(self.sample_names)
            or self.format_definitions != tuple(sorted(set(self.format_definitions)))
            or self.phase_capability_fields
            != tuple(sorted(set(self.format_definitions).intersection({"PGT", "PID", "PS"})))
            or not self.declared_contigs
            or len(declared) != len(self.declared_contigs)
            or self.class_field_definitions != tuple(sorted(set(self.class_field_definitions)))
            or len({item[:2] for item in self.class_field_definitions})
            != len(self.class_field_definitions)
            or not self.raw_lines
            or any(type(line) is not str or not line.isascii() for line in self.raw_lines)
            or len("".join(self.raw_lines)) != self.header_bytes
            or hashlib.sha256("".join(self.raw_lines).encode("ascii")).hexdigest() != self.header_sha256
        ):
            raise VcfQcError("streaming VCF QC header is invalid")


@dataclass(frozen=True, slots=True)
class ContigQcSummary:
    contig: str
    declared_length: int
    records: int
    alternate_alleles: int
    sample_calls: int
    callable_calls: int
    uncallable_calls: int
    unknown_callability_calls: int
    minimum_position: int | None
    maximum_position: int | None


@dataclass(frozen=True, slots=True)
class SampleQcSummary:
    sample_name: str
    observed_ploidies: tuple[int, ...]
    records: int
    callable_records: int
    uncallable_records: int
    unknown_callability_records: int
    phased_records: int
    unphased_records: int
    partially_phased_records: int
    phase_not_applicable_records: int
    unknown_phase_records: int


@dataclass(frozen=True, slots=True)
class NamedCount:
    name: str
    count: int


@dataclass(frozen=True, slots=True)
class VcfQcReport:
    input_sha256: str
    reference_concordance_sha256: str
    sample_names: tuple[str, ...]
    sample_count: int
    phase_capability_fields: tuple[str, ...]
    declared_contigs: tuple[DeclaredContig, ...]
    total_records: int
    total_alternate_alleles: int
    records: tuple[VcfRecordQc, ...]
    contig_summaries: tuple[ContigQcSummary, ...]
    sample_summaries: tuple[SampleQcSummary, ...]
    variant_class_counts: tuple[NamedCount, ...]
    filter_counts: tuple[NamedCount, ...]
    class_field_definitions: tuple[tuple[str, str, str, str], ...]
    classification: ArtifactClassification
    llm_policy: LLMPolicy
    schema: str = VCF_QC_SCHEMA
    evidence_policy_sha256: str = POLICY_SHA256

    pid_phase_contract: PidPhaseContract | None = None

    def __post_init__(self) -> None:
        if self.pid_phase_contract is not None and (
                type(self.pid_phase_contract) is not PidPhaseContract
                or self.pid_phase_contract.input_sha256 != self.input_sha256):
            raise VcfQcError("VCF report phase contract is invalid")
        expected_pid = None if self.pid_phase_contract is None else self.pid_phase_contract.sha256
        if any(call.pid_contract_sha256 != expected_pid for record in self.records for call in record.sample_calls):
            raise VcfQcError("VCF sample phase contract is invalid")
        require_current_policy(self.evidence_policy_sha256)
        expected_record_numbers = tuple(range(1, len(self.records) + 1))
        actual_record_numbers = tuple(record.record_number for record in self.records)
        declared_contigs = {item.contig: item.length for item in self.declared_contigs}
        expected_contigs = _contig_summaries(self.records, declared_contigs)
        expected_samples = _sample_summaries(self.records, self.sample_names)
        expected_classes = _named_counts(
            variant_class.value
            for record in self.records
            for variant_class in record.variant_classes
        )
        expected_filters = _named_counts(record.filter_value for record in self.records)
        if (
            self.schema != VCF_QC_SCHEMA
            or SHA256.fullmatch(self.input_sha256) is None
            or SHA256.fullmatch(self.reference_concordance_sha256) is None
            or not self.sample_names
            or len(set(self.sample_names)) != len(self.sample_names)
            or self.sample_count != len(self.sample_names)
            or self.phase_capability_fields
            != tuple(sorted(set(self.phase_capability_fields)))
            or any(
                field not in {"PGT", "PID", "PS"}
                for field in self.phase_capability_fields
            )
            or not self.declared_contigs
            or len(declared_contigs) != len(self.declared_contigs)
            or not _valid_count(self.total_records)
            or self.total_records <= 0
            or not _valid_count(self.total_alternate_alleles)
            or self.total_records != len(self.records)
            or self.total_alternate_alleles
            != sum(_alternate_allele_count(record.alternates) for record in self.records)
            or actual_record_numbers != expected_record_numbers
            or len({record.source_record_id for record in self.records})
            != len(self.records)
            or any(
                record.source_record_id
                != source_record_id(
                    self.input_sha256,
                    SourceSmallVariantRecord(
                        record.record_number,
                        record.contig,
                        record.position,
                        record.reference,
                        record.alternates,
                    ),
                )
                for record in self.records
            )
            or any(
                tuple(call.sample_name for call in record.sample_calls)
                != self.sample_names
                for record in self.records
            )
            or any(record.contig not in declared_contigs for record in self.records)
            or any(
                record.position > declared_contigs[record.contig]
                for record in self.records
                if record.contig in declared_contigs
            )
            or any(
                call.pgt_pid_declared
                != {"PGT", "PID"}.issubset(self.phase_capability_fields)
                for record in self.records
                for call in record.sample_calls
            )
            or self.contig_summaries != expected_contigs
            or self.sample_summaries != expected_samples
            or self.variant_class_counts != expected_classes
            or self.filter_counts != expected_filters
            or any(
                len(item) != 4
                or item[0] not in {"INFO", "FORMAT"}
                or item[1] not in (
                    CLASS_INFO_FIELDS if item[0] == "INFO" else CLASS_FORMAT_FIELDS
                )
                or not _valid_text(item[2])
                or not _valid_text(item[3])
                for item in self.class_field_definitions
            )
            or self.class_field_definitions != tuple(sorted(set(self.class_field_definitions)))
            or len({item[:2] for item in self.class_field_definitions})
            != len(self.class_field_definitions)
            or self.classification is not ArtifactClassification.CONTROLLED
            or self.llm_policy is not LLMPolicy.NO_LLM
        ):
            raise VcfQcError("VCF QC report is invalid")


def _named_counts(values: Iterator[str]) -> tuple[NamedCount, ...]:
    counts = Counter(values)
    return tuple(NamedCount(name, counts[name]) for name in sorted(counts))


def _contig_summaries(
    records: tuple[VcfRecordQc, ...], declared_contigs: dict[str, int]
) -> tuple[ContigQcSummary, ...]:
    contigs: dict[str, list[VcfRecordQc]] = {}
    for record in records:
        contigs.setdefault(record.contig, []).append(record)
    summaries = []
    for contig, declared_length in declared_contigs.items():
        values = contigs.get(contig, [])
        calls = tuple(call for record in values for call in record.sample_calls)
        summaries.append(
            ContigQcSummary(
                contig=contig,
                declared_length=declared_length,
                records=len(values),
                alternate_alleles=sum(
                    _alternate_allele_count(record.alternates) for record in values
                ),
                sample_calls=len(calls),
                callable_calls=sum(call.callable is True for call in calls),
                uncallable_calls=sum(call.callable is False for call in calls),
                unknown_callability_calls=sum(call.callable is None for call in calls),
                minimum_position=min(
                    (record.position for record in values), default=None
                ),
                maximum_position=max(
                    (record.position for record in values), default=None
                ),
            )
        )
    return tuple(summaries)


def _sample_summaries(
    records: tuple[VcfRecordQc, ...], sample_names: tuple[str, ...]
) -> tuple[SampleQcSummary, ...]:
    summaries = []
    for sample_index, sample_name in enumerate(sample_names):
        calls = tuple(record.sample_calls[sample_index] for record in records)
        summaries.append(
            SampleQcSummary(
                sample_name=sample_name,
                observed_ploidies=tuple(
                    sorted({call.ploidy for call in calls if call.ploidy is not None})
                ),
                records=len(calls),
                callable_records=sum(call.callable is True for call in calls),
                uncallable_records=sum(call.callable is False for call in calls),
                unknown_callability_records=sum(call.callable is None for call in calls),
                phased_records=sum(
                    call.phase_status is PhaseStatus.PHASED for call in calls
                ),
                unphased_records=sum(
                    call.phase_status is PhaseStatus.UNPHASED for call in calls
                ),
                partially_phased_records=sum(
                    call.phase_status is PhaseStatus.PARTIALLY_PHASED for call in calls
                ),
                phase_not_applicable_records=sum(
                    call.phase_status is PhaseStatus.NOT_APPLICABLE for call in calls
                ),
                unknown_phase_records=sum(
                    call.phase_status is PhaseStatus.UNKNOWN for call in calls
                ),
            )
        )
    return tuple(summaries)


@dataclass(frozen=True, slots=True, repr=False)
class VcfQcIncrementalSummary:
    input_sha256: str
    reference_concordance_sha256: str
    header_sha256: str
    sample_names: tuple[str, ...]
    phase_capability_fields: tuple[str, ...]
    declared_contigs: tuple[DeclaredContig, ...]
    total_records: int
    total_physical_alternate_slots: int
    total_alternate_alleles: int
    contig_summaries: tuple[ContigQcSummary, ...]
    sample_summaries: tuple[SampleQcSummary, ...]
    variant_class_counts: tuple[NamedCount, ...]
    filter_counts: tuple[NamedCount, ...]
    class_field_definitions: tuple[tuple[str, str, str, str], ...]
    classification: ArtifactClassification = ArtifactClassification.CONTROLLED
    llm_policy: LLMPolicy = LLMPolicy.NO_LLM
    evidence_policy_sha256: str = POLICY_SHA256

    def __post_init__(self) -> None:
        require_current_policy(self.evidence_policy_sha256)
        declared = {item.contig: item.length for item in self.declared_contigs}
        contig_names = tuple(item.contig for item in self.contig_summaries)
        sample_names = tuple(item.sample_name for item in self.sample_summaries)

        def valid_named_counts(values: tuple[NamedCount, ...]) -> bool:
            return (
                values == tuple(sorted(values, key=lambda item: item.name))
                and len({item.name for item in values}) == len(values)
                and all(_valid_text(item.name) and _valid_count(item.count) and item.count > 0
                        for item in values)
            )

        if (
            SHA256.fullmatch(self.input_sha256) is None
            or SHA256.fullmatch(self.reference_concordance_sha256) is None
            or SHA256.fullmatch(self.header_sha256) is None
            or not self.sample_names
            or len(set(self.sample_names)) != len(self.sample_names)
            or self.phase_capability_fields
            != tuple(sorted(set(self.phase_capability_fields)))
            or any(field not in {"PGT", "PID", "PS"}
                   for field in self.phase_capability_fields)
            or not self.declared_contigs
            or len(declared) != len(self.declared_contigs)
            or contig_names != tuple(declared)
            or sample_names != self.sample_names
            or not _valid_count(self.total_records)
            or self.total_records <= 0
            or not _valid_count(self.total_physical_alternate_slots)
            or self.total_physical_alternate_slots < self.total_records
            or not _valid_count(self.total_alternate_alleles)
            or not 0 <= self.total_alternate_alleles <= self.total_physical_alternate_slots
            or sum(item.records for item in self.contig_summaries) != self.total_records
            or sum(item.alternate_alleles for item in self.contig_summaries)
            != self.total_alternate_alleles
            or any(
                item.declared_length != declared[item.contig]
                or any(not _valid_count(value) for value in (
                    item.records, item.alternate_alleles, item.sample_calls,
                    item.callable_calls, item.uncallable_calls,
                    item.unknown_callability_calls,
                ))
                or item.sample_calls != item.records * len(self.sample_names)
                or item.callable_calls + item.uncallable_calls
                + item.unknown_callability_calls != item.sample_calls
                or (item.records == 0) != (item.minimum_position is None)
                or (item.records == 0) != (item.maximum_position is None)
                or item.records > 0 and (
                    not isinstance(item.minimum_position, int)
                    or isinstance(item.minimum_position, bool)
                    or not isinstance(item.maximum_position, int)
                    or isinstance(item.maximum_position, bool)
                    or not 0 < item.minimum_position <= item.maximum_position
                    <= item.declared_length
                )
                for item in self.contig_summaries
            )
            or any(
                item.records != self.total_records
                or item.observed_ploidies
                != tuple(sorted(set(item.observed_ploidies)))
                or any(not isinstance(value, int) or isinstance(value, bool) or value <= 0
                       for value in item.observed_ploidies)
                or any(not _valid_count(value) for value in (
                    item.records, item.callable_records, item.uncallable_records,
                    item.unknown_callability_records, item.phased_records,
                    item.unphased_records, item.partially_phased_records,
                    item.phase_not_applicable_records, item.unknown_phase_records,
                ))
                or item.callable_records + item.uncallable_records
                + item.unknown_callability_records != item.records
                or item.phased_records + item.unphased_records
                + item.partially_phased_records + item.phase_not_applicable_records
                + item.unknown_phase_records != item.records
                for item in self.sample_summaries
            )
            or not valid_named_counts(self.variant_class_counts)
            or sum(item.count for item in self.variant_class_counts)
            != self.total_physical_alternate_slots
            or any(item.name not in {value.value for value in VariantClass}
                   for item in self.variant_class_counts)
            or not valid_named_counts(self.filter_counts)
            or sum(item.count for item in self.filter_counts) != self.total_records
            or any(
                len(item) != 4
                or item[0] not in {"INFO", "FORMAT"}
                or item[1] not in (
                    CLASS_INFO_FIELDS if item[0] == "INFO" else CLASS_FORMAT_FIELDS
                )
                or not _valid_text(item[2])
                or not _valid_text(item[3])
                for item in self.class_field_definitions
            )
            or self.class_field_definitions
            != tuple(sorted(set(self.class_field_definitions)))
            or len({item[:2] for item in self.class_field_definitions})
            != len(self.class_field_definitions)
            or self.classification is not ArtifactClassification.CONTROLLED
            or self.llm_policy is not LLMPolicy.NO_LLM
        ):
            raise VcfQcError("incremental VCF QC summary is invalid")


def _sqlite_rows(cursor, batch_size=1024):
    while True:
        batch = cursor.fetchmany(batch_size)
        if not batch:
            return
        yield from batch


class QcSummaryAccumulator:
    """Disk-backed QC reducers; no source record object is retained."""

    __slots__ = (
        "header", "_header_context", "expected_records", "scratch_path", "_connection", "_records",
        "_physical", "_biological", "_finalized", "_closed", "_owns_connection",
        "_pending_seen", "_pending_contigs", "_pending_samples",
        "_pending_ploidies", "_pending_named", "_flushed_records",
        "_flushed_physical", "_flushed_biological",
    )

    def __init__(self, header: VcfQcStreamHeader, *, expected_records: int,
                 scratch_path: str | Path | None = None,
                 connection: sqlite3.Connection | None = None,
                 resume: bool = False) -> None:
        if (
            not isinstance(header, VcfQcStreamHeader)
            or not isinstance(expected_records, int) or isinstance(expected_records, bool)
            or expected_records <= 0
            or (connection is None) == (scratch_path is None)
            or type(resume) is not bool or (resume and connection is None)
        ):
            raise VcfQcError("incremental VCF QC input is invalid")
        self.header = header
        self._header_context = prepare_header(header)
        self.expected_records = expected_records
        self.scratch_path = Path(scratch_path) if scratch_path is not None else None
        self._owns_connection = connection is None
        descriptor = -1
        try:
            if connection is None:
                descriptor = os.open(
                    self.scratch_path, os.O_WRONLY | os.O_CREAT | os.O_EXCL | os.O_NOFOLLOW,
                    0o600,
                )
                os.fchmod(descriptor, 0o600)
                os.close(descriptor)
                descriptor = -1
                connection = sqlite3.connect(self.scratch_path)
                connection.execute("PRAGMA foreign_keys=ON")
                connection.execute("PRAGMA cache_size=-8192")
                connection.execute("PRAGMA temp_store=FILE")
                connection.execute("PRAGMA mmap_size=0")
                connection.execute("PRAGMA trusted_schema=OFF")
            self._connection = connection
            if not resume:
                self._connection.executescript("""
                CREATE TABLE qc_seen(record_number INTEGER PRIMARY KEY, source_record_id TEXT UNIQUE);
                CREATE TABLE qc_contigs(
                    ordinal INTEGER UNIQUE, contig TEXT PRIMARY KEY, declared_length INTEGER,
                    records INTEGER, alternate_alleles INTEGER, sample_calls INTEGER,
                    callable_calls INTEGER, uncallable_calls INTEGER,
                    unknown_callability_calls INTEGER, minimum_position INTEGER,
                    maximum_position INTEGER
                );
                CREATE TABLE qc_samples(
                    ordinal INTEGER UNIQUE, sample_name TEXT PRIMARY KEY, records INTEGER,
                    callable_records INTEGER, uncallable_records INTEGER,
                    unknown_callability_records INTEGER, phased_records INTEGER,
                    unphased_records INTEGER, partially_phased_records INTEGER,
                    phase_not_applicable_records INTEGER, unknown_phase_records INTEGER
                );
                CREATE TABLE qc_ploidies(
                    sample_name TEXT, ploidy INTEGER,
                    PRIMARY KEY(sample_name,ploidy),
                    FOREIGN KEY(sample_name) REFERENCES qc_samples(sample_name)
                );
                CREATE TABLE qc_named_counts(
                    category TEXT, name TEXT, count INTEGER,
                    PRIMARY KEY(category,name)
                );
                CREATE TABLE qc_totals(
                    singleton INTEGER PRIMARY KEY CHECK(singleton=1),
                    records INTEGER NOT NULL, physical INTEGER NOT NULL,
                    biological INTEGER NOT NULL
                );
                CREATE TABLE qc_meta(
                    singleton INTEGER PRIMARY KEY CHECK(singleton=1),
                    input_sha256 TEXT NOT NULL,
                    concordance_sha256 TEXT NOT NULL,
                    header_sha256 TEXT NOT NULL,
                    expected_records INTEGER NOT NULL
                );
            """)
                self._connection.executemany(
                    "INSERT INTO qc_contigs VALUES (?,?,?,?,?,?,?,?,?,?,?)",
                    ((ordinal, item.contig, item.length, 0, 0, 0, 0, 0, 0, None, None)
                     for ordinal, item in enumerate(header.declared_contigs, 1)),
                )
                self._connection.executemany(
                    "INSERT INTO qc_samples VALUES (?,?,?,?,?,?,?,?,?,?,?)",
                    ((ordinal, sample, 0, 0, 0, 0, 0, 0, 0, 0, 0)
                     for ordinal, sample in enumerate(header.sample_names, 1)),
                )
                self._connection.execute("INSERT INTO qc_totals VALUES (1,0,0,0)")
                self._connection.execute(
                    "INSERT INTO qc_meta VALUES (1,?,?,?,?)",
                    (header.input_sha256, header.reference_concordance_sha256,
                     header.header_sha256, expected_records),
                )
                if self._owns_connection:
                    self._connection.commit()
                totals = (0, 0, 0)
            else:
                totals = self._connection.execute(
                    "SELECT records,physical,biological FROM qc_totals WHERE singleton=1"
                ).fetchone()
                if (totals is None or len(totals) != 3
                        or not all(type(value) is int and value >= 0 for value in totals)
                        or totals[0] > expected_records
                        or self._connection.execute(
                            "SELECT input_sha256,concordance_sha256,header_sha256,expected_records "
                            "FROM qc_meta WHERE singleton=1"
                        ).fetchone() != (
                            header.input_sha256, header.reference_concordance_sha256,
                            header.header_sha256, expected_records,
                        )
                        or self._connection.execute(
                            "SELECT count(*) FROM qc_seen"
                        ).fetchone() != (totals[0],)
                        or self._connection.execute(
                            "SELECT contig,declared_length FROM qc_contigs ORDER BY ordinal"
                        ).fetchall() != [
                            (item.contig, item.length) for item in header.declared_contigs
                        ]
                        or self._connection.execute(
                            "SELECT sample_name FROM qc_samples ORDER BY ordinal"
                        ).fetchall() != [(name,) for name in header.sample_names]):
                    raise VcfQcError("incremental VCF QC state is inconsistent")
        except (OSError, sqlite3.Error):
            if descriptor >= 0:
                os.close(descriptor)
            raise VcfQcError("cannot create incremental VCF QC state") from None
        self._records, self._physical, self._biological = totals
        self._finalized = False
        self._closed = False
        self._pending_seen = []
        self._pending_contigs = {}
        self._pending_samples = {}
        self._pending_ploidies = set()
        self._pending_named = Counter()
        self._flushed_records, self._flushed_physical, self._flushed_biological = totals

    def flush(self) -> None:
        """Insert bounded source uniqueness facts in one coordinator batch."""
        if self._closed or self._finalized:
            raise VcfQcError("incremental VCF QC state is unavailable")
        if self._pending_seen:
            try:
                self._connection.executemany(
                    "INSERT INTO qc_seen VALUES (?,?)", self._pending_seen,
                )
                self._connection.executemany("""
                    UPDATE qc_contigs SET records=records+?,
                        alternate_alleles=alternate_alleles+?,
                        sample_calls=sample_calls+?, callable_calls=callable_calls+?,
                        uncallable_calls=uncallable_calls+?,
                        unknown_callability_calls=unknown_callability_calls+?,
                        minimum_position=CASE WHEN minimum_position IS NULL OR minimum_position>?
                                              THEN ? ELSE minimum_position END,
                        maximum_position=CASE WHEN maximum_position IS NULL OR maximum_position<?
                                              THEN ? ELSE maximum_position END
                    WHERE contig=?
                """, ((values[0], values[1], values[2], values[3], values[4], values[5],
                       values[6], values[6], values[7], values[7], contig)
                      for contig, values in self._pending_contigs.items()))
                self._connection.executemany("""
                    UPDATE qc_samples SET records=records+?,
                        callable_records=callable_records+?,
                        uncallable_records=uncallable_records+?,
                        unknown_callability_records=unknown_callability_records+?,
                        phased_records=phased_records+?, unphased_records=unphased_records+?,
                        partially_phased_records=partially_phased_records+?,
                        phase_not_applicable_records=phase_not_applicable_records+?,
                        unknown_phase_records=unknown_phase_records+?
                    WHERE sample_name=?
                """, ((*values, sample) for sample, values in self._pending_samples.items()))
                self._connection.executemany(
                    "INSERT OR IGNORE INTO qc_ploidies VALUES (?,?)",
                    self._pending_ploidies,
                )
                self._connection.executemany("""
                    INSERT INTO qc_named_counts VALUES (?,?,?)
                    ON CONFLICT(category,name) DO UPDATE SET count=count+excluded.count
                """, ((category, name, count)
                      for (category, name), count in self._pending_named.items()))
                self._connection.execute(
                    "UPDATE qc_totals SET records=records+?, physical=physical+?, "
                    "biological=biological+? WHERE singleton=1",
                    (self._records - self._flushed_records,
                     self._physical - self._flushed_physical,
                     self._biological - self._flushed_biological),
                )
            except sqlite3.Error:
                raise VcfQcError("cannot update incremental VCF QC state") from None
            self._pending_seen.clear()
            self._pending_contigs.clear()
            self._pending_samples.clear()
            self._pending_ploidies.clear()
            self._pending_named.clear()
            self._flushed_records = self._records
            self._flushed_physical = self._physical
            self._flushed_biological = self._biological

    def add(self, record: VcfRecordQc) -> None:
        if self._closed or self._finalized or not isinstance(record, VcfRecordQc):
            raise VcfQcError("incremental VCF QC state is unavailable")
        expected_id = source_record_id(
            self.header.input_sha256,
            SourceSmallVariantRecord(
                record.record_number, record.contig, record.position,
                record.reference, record.alternates,
            ),
        )
        declared = self._header_context.contig_lengths
        if (
            record.record_number != self._records + 1
            or record.contig not in declared
            or record.position > declared.get(record.contig, 0)
            or record.source_record_id != expected_id
        ):
            raise VcfQcError("incremental VCF QC source order is invalid")
        if tuple(call.sample_name for call in record.sample_calls) != self.header.sample_names:
            raise VcfQcError("incremental VCF QC sample order is invalid")
        biological = _alternate_allele_count(record.alternates)
        physical = len(record.alternates)
        calls = record.sample_calls
        contig = self._pending_contigs.setdefault(record.contig, [0, 0, 0, 0, 0, 0,
                                                                 record.position,
                                                                 record.position])
        contig[0] += 1
        contig[1] += biological
        contig[2] += len(calls)
        contig[3] += sum(call.callable is True for call in calls)
        contig[4] += sum(call.callable is False for call in calls)
        contig[5] += sum(call.callable is None for call in calls)
        contig[6] = min(contig[6], record.position)
        contig[7] = max(contig[7], record.position)
        for call in calls:
            sample = self._pending_samples.setdefault(call.sample_name, [0] * 9)
            sample[0] += 1
            sample[1] += int(call.callable is True)
            sample[2] += int(call.callable is False)
            sample[3] += int(call.callable is None)
            sample[4] += int(call.phase_status is PhaseStatus.PHASED)
            sample[5] += int(call.phase_status is PhaseStatus.UNPHASED)
            sample[6] += int(call.phase_status is PhaseStatus.PARTIALLY_PHASED)
            sample[7] += int(call.phase_status is PhaseStatus.NOT_APPLICABLE)
            sample[8] += int(call.phase_status is PhaseStatus.UNKNOWN)
            if call.ploidy is not None:
                self._pending_ploidies.add((call.sample_name, call.ploidy))
        for category, name in (
            *(("variant", item.value) for item in record.variant_classes),
            ("filter", record.filter_value),
        ):
            self._pending_named[(category, name)] += 1
        self._records += 1
        self._physical += physical
        self._biological += biological
        self._pending_seen.append((record.record_number, record.source_record_id))
        if len(self._pending_seen) >= 1024:
            self.flush()
        if self._owns_connection and self._records % 1024 == 0:
            self._connection.commit()

    def finalize(self) -> VcfQcIncrementalSummary:
        if self._closed or self._finalized or self._records != self.expected_records:
            raise VcfQcError("incremental VCF QC record accounting is inconsistent")
        self.flush()
        if self._owns_connection:
            self._connection.commit()
        contigs = tuple(ContigQcSummary(*row[1:]) for row in _sqlite_rows(
            self._connection.execute("SELECT * FROM qc_contigs ORDER BY ordinal")
        ))
        samples = []
        for row in _sqlite_rows(self._connection.execute("SELECT * FROM qc_samples ORDER BY ordinal")):
            ploidies = tuple(value[0] for value in _sqlite_rows(self._connection.execute(
                "SELECT ploidy FROM qc_ploidies WHERE sample_name=? ORDER BY ploidy", (row[1],)
            )))
            samples.append(SampleQcSummary(row[1], ploidies, *row[2:]))
        named = {
            category: tuple(NamedCount(name, count) for _category, name, count in _sqlite_rows(
                self._connection.execute(
                    "SELECT category,name,count FROM qc_named_counts WHERE category=? ORDER BY name",
                    (category,),
                )
            ))
            for category in ("variant", "filter")
        }
        self._finalized = True
        return VcfQcIncrementalSummary(
            self.header.input_sha256, self.header.reference_concordance_sha256,
            self.header.header_sha256, self.header.sample_names,
            self.header.phase_capability_fields, self.header.declared_contigs,
            self._records, self._physical, self._biological, contigs, tuple(samples),
            named["variant"], named["filter"], self.header.class_field_definitions,
        )

    def close(self) -> None:
        if not self._closed:
            if self._owns_connection:
                self._connection.close()
            self._closed = True

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        incomplete = exc_type is None and not self._finalized
        self.close()
        if incomplete:
            raise VcfQcError("incremental VCF QC state was not finalized")
        return False


def _variant_class(reference: str, alternate: str) -> VariantClass:
    if alternate == ".":
        return VariantClass.NON_VARIANT
    if alternate == "*":
        return VariantClass.SPANNING_DELETION
    if "[" in alternate or "]" in alternate:
        return VariantClass.BREAKEND
    if alternate.startswith("<") and alternate.endswith(">"):
        return VariantClass.SYMBOLIC
    if ALTERNATE_ALLELE.fullmatch(alternate) is None:
        return VariantClass.OTHER
    if alternate == reference:
        return VariantClass.NON_VARIANT
    trimmed_reference = reference
    trimmed_alternate = alternate
    while trimmed_reference and trimmed_alternate and trimmed_reference[-1] == trimmed_alternate[-1]:
        trimmed_reference = trimmed_reference[:-1]
        trimmed_alternate = trimmed_alternate[:-1]
    while trimmed_reference and trimmed_alternate and trimmed_reference[0] == trimmed_alternate[0]:
        trimmed_reference = trimmed_reference[1:]
        trimmed_alternate = trimmed_alternate[1:]
    if not trimmed_reference:
        return VariantClass.INSERTION
    if not trimmed_alternate:
        return VariantClass.DELETION
    if len(trimmed_reference) == len(trimmed_alternate) == 1:
        return VariantClass.SNV
    if len(trimmed_reference) == len(trimmed_alternate):
        return VariantClass.MNV
    return VariantClass.DELINS


def _parse_quality(token):
    if token == ".":
        return None
    try:
        value = float(token)
    except (TypeError, ValueError):
        raise VcfQcError("VCF QC quality is invalid") from None
    if not math.isfinite(value) or value < 0:
        raise VcfQcError("VCF QC quality is invalid")
    return value


def compare_quality_rendering(record, rendered_token, *, adapter_version):
    """Compare one rendered QUAL token without replacing any source evidence."""
    if not adapter_version or record.raw_quality is None:
        raise VcfQcError("numeric rendering comparison lacks source or adapter identity")
    rendered = _parse_quality(rendered_token)
    try:
        same_decimal = (rendered_token == record.raw_quality if rendered is None or record.quality is None
                        else Decimal(rendered_token) == Decimal(record.raw_quality))
    except InvalidOperation:
        raise VcfQcError("numeric rendering token is invalid") from None
    return {"adapter_version": adapter_version, "source_record_id": record.source_record_id,
            "source_line_sha256": record.raw_line_sha256, "raw_qual": record.raw_quality,
            "typed_qual": record.quality, "typed_conversion": "source-decimal-to-binary64",
            "rendered_qual": rendered_token, "decimal_values_equal": same_decimal,
            "inspected_records": 1, "whole_file_equivalence": "not-established"}


def _parse_optional_integer(value: str | None) -> int | None:
    if value in {None, "."}:
        return None
    if not value.isdigit():
        raise VcfQcError("VCF QC contains an invalid integer field")
    return int(value)


def _parse_allelic_depths(value: str | None) -> tuple[int | None, ...] | None:
    if value in {None, "."}:
        return None
    fields = value.split(",")
    if not fields or any(field != "." and not field.isdigit() for field in fields):
        raise VcfQcError("VCF QC contains invalid allelic depths")
    return tuple(None if field == "." else int(field) for field in fields)


def _parse_info(value: str) -> dict[str, str | None]:
    if value == ".":
        return {}
    fields: dict[str, str | None] = {}
    for item in value.split(";"):
        if not item:
            raise VcfQcError("VCF QC contains an invalid INFO field")
        key, separator, field_value = item.partition("=")
        if not key or key in fields:
            raise VcfQcError("VCF QC contains duplicate or empty INFO fields")
        fields[key] = field_value if separator else None
    return fields


def _parse_callset_frequencies(value: str | None) -> tuple[float | None, ...] | None:
    if value is None or value == ".":
        return None
    frequencies = []
    for field in value.split(","):
        if field == ".":
            frequencies.append(None)
            continue
        try:
            frequency = float(field)
        except ValueError as error:
            raise VcfQcError("VCF QC contains an invalid callset frequency") from error
        if not math.isfinite(frequency) or frequency < 0.0 or frequency > 1.0:
            raise VcfQcError("VCF QC contains an invalid callset frequency")
        frequencies.append(frequency)
    return tuple(frequencies)


def _filter(value: str) -> tuple[FilterStatus, tuple[str, ...]]:
    if value == "PASS":
        return FilterStatus.PASS, ()
    if value == ".":
        return FilterStatus.NOT_APPLIED, ()
    filters = tuple(value.split(";"))
    if any(not _valid_text(item) for item in filters) or len(set(filters)) != len(filters):
        raise VcfQcError("VCF QC contains invalid FILTER values")
    return FilterStatus.FILTERED, filters


def _sample_call(
    sample_name: str,
    format_keys: tuple[str, ...],
    sample_value: str,
    alternate_count: int,
    pgt_pid_declared: bool,
    pid_phase_contract: PidPhaseContract | None = None,
    positions=None,
) -> SampleCallQc:
    if (len(set(format_keys)) != len(format_keys)
            or any(not _valid_text(key) for key in format_keys)
            or "GT" in format_keys and format_keys[0] != "GT"):
        raise VcfQcError("VCF QC FORMAT keys are invalid")
    sample_fields = sample_value.split(":")
    if not format_keys and sample_value == ".":
        sample_fields = []
    if len(sample_fields) > len(format_keys):
        raise VcfQcError("VCF QC FORMAT and sample fields are inconsistent")
    format_fields = tuple(FormatFieldEvidence(key,
        sample_fields[index] if index < len(sample_fields) else None,
        field_state(sample_fields[index] if index < len(sample_fields) else None),
        key in FORMAT_CONTRACTS) for index, key in enumerate(format_keys))
    positions = positions if positions is not None else {key: index for index, key in enumerate(format_keys)}
    consumed = FORMAT_CONTRACTS.keys() | CLASS_FORMAT_FIELDS
    values = {key: None if format_fields[index].raw_value in (None, ".") else format_fields[index].raw_value
              for key, index in positions.items() if key in consumed}
    genotype = values.get("GT")
    if genotype is not None:
        if GENOTYPE.fullmatch(genotype) is None:
            raise VcfQcError("VCF QC contains an invalid genotype")
        if any(
            token != "." and int(token) > alternate_count
            for token in _gt_tokens(genotype)
        ):
            raise VcfQcError("VCF QC genotype allele number is out of range")
    pgt = values.get("PGT")
    if pgt is not None and GENOTYPE.fullmatch(pgt) is None:
        raise VcfQcError("VCF QC contains invalid physical-phase data")
    phase_set = values.get("PS")
    if phase_set is not None and not phase_set.isdigit():
        raise VcfQcError("VCF QC contains invalid phase-set data")
    genotype_quality = _parse_optional_integer(values.get("GQ"))
    depth = _parse_optional_integer(values.get("DP"))
    allelic_depths = _parse_allelic_depths(values.get("AD"))
    if allelic_depths is not None and len(allelic_depths) != expected_cardinality("R", alternate_count):
        raise VcfQcError("VCF QC allelic depths do not match the allele count")
    allele_balance = _expected_allele_balance(allelic_depths, depth)
    phase_status, phase_evidence = _phase(
        genotype,
        pgt,
        values.get("PID"),
        pgt_pid_declared,
        None if pid_phase_contract is None else pid_phase_contract.sha256,
    )
    return SampleCallQc(
        sample_name=sample_name,
        format_fields=format_fields,
        pid_contract_sha256=None if pid_phase_contract is None else pid_phase_contract.sha256,
        genotype=genotype,
        ploidy=_ploidy(genotype),
        phase_status=phase_status,
        phase_evidence=phase_evidence,
        pgt_pid_declared=pgt_pid_declared,
        pgt=pgt,
        pid=values.get("PID"),
        phase_set=phase_set,
        genotype_quality=genotype_quality,
        depth=depth,
        allelic_depths=allelic_depths,
        alternate_allele_balance=allele_balance,
        allele_balance_source="format-ad-dp" if allele_balance is not None else None,
        callable=_callability(genotype, depth),
        class_format_fields=tuple(
            (key, values[key]) for key in sorted(CLASS_FORMAT_FIELDS.intersection(values))
        ),
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
        raise VcfQcError("cannot read VCF QC input") from error
    if initial != final_descriptor or initial != final_path:
        raise VcfQcError("VCF QC input changed while being read")
    return digest.hexdigest(), initial


@contextmanager
def _open_vcf(path: Path) -> Iterator[BinaryIO]:
    try:
        with path.open("rb") as raw:
            compressed = raw.read(2) == b"\x1f\x8b"
            raw.seek(0)
            if compressed:
                with gzip.GzipFile(fileobj=raw, mode="rb") as source:
                    yield source
            else:
                yield raw
    except (EOFError, gzip.BadGzipFile, OSError, zlib.error) as error:
        raise VcfQcError("cannot parse VCF QC input") from error


def _validate_concordance(concordance_report) -> None:
    if (
        not isinstance(
            concordance_report,
            (ReferenceConcordanceReport, PublicValidationConcordanceReport),
        )
        or not concordance_report.header_valid
        or not concordance_report.reference_basename_match
        or not concordance_report.contig_names_match
        or not concordance_report.contig_lengths_match
        or not concordance_report.contig_order_match
        or concordance_report.unresolved_records != 0
    ):
        raise VcfQcError("VCF QC concordance binding is invalid")


def _bounded_line(source: BinaryIO, limit: int, *, header: bool) -> bytes:
    try:
        return bounded_line(source, limit, header=header, error_type=VcfQcError)
    except (EOFError, gzip.BadGzipFile, OSError, zlib.error) as error:
        raise VcfQcError("cannot parse VCF QC input") from error


def _ascii_line(raw_line: bytes) -> str:
    try:
        return raw_line.rstrip(b"\r\n").decode("ascii")
    except UnicodeDecodeError as error:
        raise VcfQcError("VCF QC input is not ASCII") from error


def _sample_columns(line: str) -> tuple[str, ...]:
    columns = line.split("\t")
    if columns[:9] != [
        "#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT",
    ]:
        raise VcfQcError("VCF QC column header is invalid")
    samples = tuple(columns[9:])
    if (not samples or len(set(samples)) != len(samples)
            or any(not name or len(name) > 1024
                   or any(char.isspace() or ord(char) < 33 or ord(char) == 127
                          for char in name) for name in samples)):
        raise VcfQcError("VCF QC sample header is invalid")
    return samples


def read_vcf_samples(vcf_path: str | Path, *, max_header_bytes: int = DEFAULT_MAX_HEADER_BYTES
                     ) -> tuple[str, ...]:
    """Read only a bounded header, sharing the source parser's syntax rules.

    Names are private operator data. This performs no body or reference check.
    """
    if type(max_header_bytes) is not int or max_header_bytes <= 0:
        raise VcfQcError("VCF QC streaming limits are invalid")
    path = Path(vcf_path)
    before = _stat_identity(path.stat())
    formats, contigs = set(), set()
    used = 0
    with _open_vcf(path) as source:
        while True:
            raw = _bounded_line(source, max_header_bytes - used, header=True)
            used += len(raw)
            line = _ascii_line(raw)
            if used == len(raw) and not re.fullmatch(r"##fileformat=VCFv4\.[0-9]+", line):
                raise VcfQcError("VCF QC header ordering is invalid")
            if line.startswith("#CHROM\t"):
                samples = _sample_columns(line)
                break
            if not line.startswith("##"):
                raise VcfQcError("VCF QC header ordering is invalid")
            if line.startswith("##FORMAT="):
                name = _parse_format_header(line)
                if name in formats:
                    raise VcfQcError("VCF QC FORMAT declaration is invalid")
                formats.add(name)
            if line.startswith("##contig="):
                name, _ = _parse_contig_header(line)
                if name in contigs:
                    raise VcfQcError("VCF QC contig declaration is duplicated")
                contigs.add(name)
        if _stat_identity(os.fstat(source.fileno())) != before:
            raise VcfQcError("VCF QC input changed while being read")
    if _stat_identity(path.stat()) != before:
        raise VcfQcError("VCF QC input changed while being read")
    return samples


def _parse_stream_header(
    source: BinaryIO, *, input_sha256: str, input_size_bytes: int,
    concordance_report, max_header_bytes: int,
    max_record_bytes: int = DEFAULT_MAX_RECORD_BYTES,
    pid_phase_contract: PidPhaseContract | None = None,
) -> VcfQcStreamHeader:
    sample_names: tuple[str, ...] | None = None
    format_definitions: set[str] = set()
    declared_contigs: dict[str, int] = {}
    class_field_definitions: dict[tuple[str, str], tuple[str, str, str, str]] = {}
    header_digest = hashlib.sha256()
    header_bytes = 0
    raw_lines = []
    while sample_names is None:
        raw_line = _bounded_line(
            source, max_header_bytes - header_bytes, header=True,
        )
        if not raw_line:
            raise VcfQcError("VCF QC header ordering is invalid")
        header_digest.update(raw_line)
        header_bytes += len(raw_line)
        line = _ascii_line(raw_line)
        raw_lines.append(raw_line.decode("ascii"))
        if line.startswith(("##INFO=<", "##FORMAT=<")):
            namespace = line[2:].split("=", 1)[0]
            attributes = _structured_attributes(line.split("=<", 1)[1][:-1])
            allowed = CLASS_INFO_FIELDS if namespace == "INFO" else CLASS_FORMAT_FIELDS
            if attributes is not None and attributes.get("ID") in allowed:
                identifier = attributes["ID"]
                number, field_type = attributes.get("Number"), attributes.get("Type")
                key = (namespace, identifier)
                if (
                    not line.endswith(">")
                    or not _valid_text(number)
                    or not _valid_text(field_type)
                    or key in class_field_definitions
                ):
                    raise VcfQcError("VCF QC class field declaration is invalid")
                class_field_definitions[key] = (namespace, identifier, number, field_type)
        if line.startswith("##FORMAT="):
            identifier = _parse_format_header(line)
            if identifier in format_definitions:
                raise VcfQcError("VCF QC FORMAT declaration is invalid")
            format_definitions.add(identifier)
            continue
        if line.startswith("##contig="):
            contig, length = _parse_contig_header(line)
            if contig in declared_contigs:
                raise VcfQcError("VCF QC contig declaration is duplicated")
            declared_contigs[contig] = length
            continue
        if line.startswith("##"):
            continue
        if line.startswith("#CHROM\t"):
            sample_names = _sample_columns(line)
            continue
        raise VcfQcError("VCF QC header ordering is invalid")
    if len(declared_contigs) != concordance_report.declared_contigs:
        raise VcfQcError("VCF QC record accounting is inconsistent")
    return VcfQcStreamHeader(
        input_sha256=input_sha256,
        input_size_bytes=input_size_bytes,
        reference_concordance_sha256=hashlib.sha256(
            concordance_report_as_json(concordance_report).encode("ascii")
        ).hexdigest(),
        header_sha256=header_digest.hexdigest(),
        header_bytes=header_bytes,
        sample_names=sample_names,
        format_definitions=tuple(sorted(format_definitions)),
        phase_capability_fields=tuple(
            sorted(format_definitions.intersection({"PGT", "PID", "PS"}))
        ),
        declared_contigs=tuple(
            DeclaredContig(contig, length) for contig, length in declared_contigs.items()
        ),
        class_field_definitions=tuple(sorted(class_field_definitions.values())),
        raw_lines=tuple(raw_lines),
        read_limits=VcfReadLimits(max_header_bytes, max_record_bytes, max_record_bytes),
        pid_phase_contract=pid_phase_contract,
    )


def _parse_stream_record(
    line: str, *, record_number: int, header: VcfQcStreamHeader,
    raw_line_sha256: str, raw_line: str,
    prepared_header=None, format_cache=None,
) -> VcfRecordQc:
    prepared = prepare_header(header) if prepared_header is None else prepared_header
    if prepared.header is not header:
        raise VcfQcError("prepared VCF header belongs to another read owner")
    fields = line.split("\t")
    if len(fields) != 9 + len(header.sample_names):
        raise VcfQcError("VCF QC record column count is invalid")
    try:
        position = int(fields[1])
        source_record = SourceSmallVariantRecord(
            record_number, fields[0], position, fields[3], tuple(fields[4].split(",")),
        )
    except (ValueError, RuntimeError) as error:
        raise VcfQcError("VCF QC source record is invalid") from error
    declared = prepared.contig_lengths
    if (
        source_record.contig not in declared
        or source_record.position > declared[source_record.contig]
    ):
        raise VcfQcError("VCF QC source record is outside the declared contigs")
    format_keys = () if fields[8] == "." else tuple(fields[8].split(":"))
    format_definitions = prepared.format_definitions
    cache = FormatLayoutCache(max_entries=0) if format_cache is None else format_cache
    try:
        positions = cache.positions(format_keys, allowed=format_definitions)
    except ValueError:
        raise VcfQcError("VCF QC FORMAT keys are invalid") from None
    info = _parse_info(fields[7])
    filter_status, filters = _filter(fields[6])
    callset_frequencies = _parse_callset_frequencies(info.get("AF"))
    biological_alternates = _alternate_allele_count(source_record.alternates)
    sample_calls = tuple(
        _sample_call(
            sample_name, format_keys, sample_value, biological_alternates,
            prepared.pgt_pid_declared, header.pid_phase_contract, positions,
        )
        for sample_name, sample_value in zip(
            header.sample_names, fields[9:], strict=True,
        )
    )
    technical_scopes = _technical_scopes(
        filter_status, sample_calls, len(source_record.alternates),
    )
    technical_required = _technical_scope_required(technical_scopes)
    quality = _parse_quality(fields[5])
    return VcfRecordQc(
        record_number=source_record.record_number,
        source_record_id=source_record_id(header.input_sha256, source_record),
        contig=source_record.contig,
        position=source_record.position,
        reference=source_record.reference,
        alternates=source_record.alternates,
        variant_classes=tuple(
            _variant_class(source_record.reference, alternate)
            for alternate in source_record.alternates
        ),
        filter_value=fields[6],
        filter_status=filter_status,
        filters=filters,
        callset_allele_number=_parse_optional_integer(info.get("AN")),
        callset_allele_frequencies=callset_frequencies,
        sample_calls=sample_calls,
        technical_penalty_required=technical_required,
        read_validation_required=technical_required,
        technical_flags=technical_scopes[0],
        alternate_technical_flags=technical_scopes[1],
        sample_technical_flags=technical_scopes[2],
        sample_alternate_technical_flags=technical_scopes[3],
        class_info_fields=tuple(
            (key, info[key]) for key in sorted(CLASS_INFO_FIELDS.intersection(info))
        ),
        quality=quality,
        raw_quality=fields[5],
        raw_line_sha256=raw_line_sha256,
        raw_line=raw_line,
    )


class VcfQcRecordIterator(Iterator[VcfRecordQc]):
    """One bounded sequential VCF pass; exhaust before accepting source closure."""

    def __init__(
        self, *, path: Path, context, source: BinaryIO, header: VcfQcStreamHeader,
        input_identity: tuple[int, int, int, int, int], expected_records: int,
        max_record_bytes: int,
    ) -> None:
        self.path = path
        self._context = context
        self._source = source
        self.header = header
        self._input_identity = input_identity
        self._expected_records = expected_records
        self._max_record_bytes = max_record_bytes
        self.prepared_header = prepare_header(header)
        self.format_cache = FormatLayoutCache()
        self._record_number = 0
        self.complete = False
        self.closed = False

    def __iter__(self):
        return self

    def skip_committed_prefix(self, records: int) -> None:
        """Skip only raw bounded lines already committed by the source stage."""
        if (type(records) is not int or not 0 <= records <= self._expected_records
                or self._record_number != 0 or self.complete or self.closed):
            raise VcfQcError("VCF QC resume cursor is invalid")
        try:
            for _ in range(records):
                if not _bounded_line(self._source, self._max_record_bytes, header=False):
                    raise VcfQcError("VCF QC resume cursor exceeds source")
                self._record_number += 1
        except Exception:
            self.close()
            raise

    def __next__(self) -> VcfRecordQc:
        if self.complete or self.closed:
            raise StopIteration
        try:
            raw_line = _bounded_line(self._source, self._max_record_bytes, header=False)
            if not raw_line:
                self._finish()
                raise StopIteration
            line = _ascii_line(raw_line)
            if line.startswith("#"):
                raise VcfQcError("VCF QC header ordering is invalid")
            self._record_number += 1
            if self._record_number > self._expected_records:
                raise VcfQcError("VCF QC record accounting is inconsistent")
            return _parse_stream_record(
                line, record_number=self._record_number, header=self.header,
                raw_line_sha256=hashlib.sha256(raw_line).hexdigest(),
                raw_line=raw_line.decode("ascii"),
                prepared_header=self.prepared_header, format_cache=self.format_cache,
            )
        except StopIteration:
            raise
        except Exception:
            self.close()
            raise

    def _finish(self) -> None:
        try:
            descriptor_identity = _stat_identity(os.fstat(self._source.fileno()))
            path_identity = _stat_identity(self.path.stat())
        except OSError as error:
            self.close()
            raise VcfQcError("cannot verify VCF QC input") from error
        if (
            self._record_number != self._expected_records
            or descriptor_identity != self._input_identity
            or path_identity != self._input_identity
        ):
            self.close()
            if self._record_number != self._expected_records:
                raise VcfQcError("VCF QC record accounting is inconsistent")
            raise VcfQcError("VCF QC input changed while being read")
        self.complete = True
        self.close()

    def close(self) -> None:
        if not self.closed:
            self.closed = True
            self._context.__exit__(None, None, None)

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        incomplete = exc_type is None and not self.complete
        self.close()
        if incomplete:
            raise VcfQcError("VCF QC stream was not exhausted")
        return False


def iter_vcf_qc_records(
    vcf_path: str | Path, *,
    concordance_report: ReferenceConcordanceReport | PublicValidationConcordanceReport,
    max_record_bytes: int = DEFAULT_MAX_RECORD_BYTES,
    max_header_bytes: int = DEFAULT_MAX_HEADER_BYTES,
    resume_input_identity: tuple[int, int, int, int, int] | None = None,
    skip_records: int = 0,
    pid_phase_contract: PidPhaseContract | None = None,
) -> VcfQcRecordIterator:
    """Open a bounded streaming parser that shares the legacy per-record kernel."""
    _validate_concordance(concordance_report)
    if (
        not isinstance(max_record_bytes, int) or isinstance(max_record_bytes, bool)
        or max_record_bytes <= 0
        or not isinstance(max_header_bytes, int) or isinstance(max_header_bytes, bool)
        or max_header_bytes <= 0 or max_header_bytes < max_record_bytes
        or type(skip_records) is not int or not 0 <= skip_records <= concordance_report.total_records
        or resume_input_identity is not None and (
            type(resume_input_identity) is not tuple
            or len(resume_input_identity) != 5
            or any(type(value) is not int or value < 0 for value in resume_input_identity)
        )
        or skip_records and resume_input_identity is None
    ):
        raise VcfQcError("VCF QC streaming limits are invalid")
    path = Path(vcf_path)
    input_sha256, input_identity = _file_sha256(path)
    if resume_input_identity is not None and input_identity != resume_input_identity:
        raise VcfQcError("VCF QC resume input identity changed")
    if input_sha256 != concordance_report.validated_input_sha256:
        raise VcfQcError("VCF QC input binding is stale")
    context = _open_vcf(path)
    source = context.__enter__()
    try:
        if _stat_identity(os.fstat(source.fileno())) != input_identity:
            raise VcfQcError("VCF QC input changed while being read")
        header = _parse_stream_header(
            source, input_sha256=input_sha256, input_size_bytes=input_identity[2],
            concordance_report=concordance_report, max_header_bytes=max_header_bytes,
            max_record_bytes=max_record_bytes, pid_phase_contract=pid_phase_contract,
        )
        stream = VcfQcRecordIterator(
            path=path, context=context, source=source, header=header,
            input_identity=input_identity, expected_records=concordance_report.total_records,
            max_record_bytes=max_record_bytes,
        )
        stream.skip_committed_prefix(skip_records)
        return stream
    except Exception:
        context.__exit__(None, None, None)
        raise


def summarize_vcf_genotype_qc(
    vcf_path: str | Path,
    *,
    concordance_report: ReferenceConcordanceReport | PublicValidationConcordanceReport,
    pid_phase_contract: PidPhaseContract | None = None,
) -> VcfQcReport:
    """Legacy fixture API built from the shared bounded streaming record kernel."""
    stream = iter_vcf_qc_records(vcf_path, concordance_report=concordance_report,
                                 pid_phase_contract=pid_phase_contract)
    with stream:
        record_tuple = tuple(stream)
    header = stream.header
    declared_contigs = {item.contig: item.length for item in header.declared_contigs}
    return VcfQcReport(
        pid_phase_contract=header.pid_phase_contract,
        input_sha256=header.input_sha256,
        reference_concordance_sha256=header.reference_concordance_sha256,
        sample_names=header.sample_names,
        sample_count=len(header.sample_names),
        phase_capability_fields=header.phase_capability_fields,
        declared_contigs=header.declared_contigs,
        total_records=len(record_tuple),
        total_alternate_alleles=sum(
            _alternate_allele_count(record.alternates) for record in record_tuple
        ),
        records=record_tuple,
        contig_summaries=_contig_summaries(record_tuple, declared_contigs),
        sample_summaries=_sample_summaries(record_tuple, header.sample_names),
        variant_class_counts=_named_counts(
            variant_class.value
            for record in record_tuple
            for variant_class in record.variant_classes
        ),
        filter_counts=_named_counts(record.filter_value for record in record_tuple),
        class_field_definitions=header.class_field_definitions,
        classification=ArtifactClassification.CONTROLLED,
        llm_policy=LLMPolicy.NO_LLM,
    )


def vcf_qc_as_json(report: VcfQcReport) -> str:
    return json.dumps(asdict(report), indent=2, sort_keys=True) + "\n"
