"""Inventory supplied call representations; never infer genome-wide negatives."""

from __future__ import annotations

from collections import Counter
from dataclasses import asdict, dataclass, field
from enum import StrEnum
import hashlib
import json
import math
import os
from pathlib import Path
import re
import sqlite3

from genosapiens.governance import ArtifactClassification
from genosapiens.llm_guard import LLMPolicy
from genosapiens.variant_normalization import (
    SHA256,
    SOURCE_RECORD_ID,
    SourceSmallVariantRecord,
    source_allele_id,
    source_record_id,
)
from genosapiens.vcf_qc import (
    SampleCallQc,
    VariantClass,
    VcfQcReport,
    VcfRecordQc,
    VcfQcStreamHeader,
    technical_flags_for,
    vcf_qc_as_json,
)


INVENTORY_SCHEMA = "genosapiens.variant-class-inventory"
AUTOSOMES = frozenset(str(number) for number in range(1, 23))
SEQUENCE_CLASSES = frozenset({
    VariantClass.SNV, VariantClass.MNV, VariantClass.INSERTION,
    VariantClass.DELETION, VariantClass.DELINS,
})
SV_TYPES = frozenset({"DEL", "DUP", "INV", "INS", "BND", "CNV", "CPX", "TRA"})
REPEAT_SYMBOL = re.compile(r"<(?:STR(?:[0-9]+)?|VNTR)(?::[A-Za-z0-9:._-]+)?>")
CN_SYMBOL = re.compile(r"<CN[0-9]+>")
BREAKEND = re.compile(
    r"(?:[ACGTN]+([\[\]])[^\s:\[\]]+:[1-9][0-9]*\1|"
    r"([\[\]])[^\s:\[\]]+:[1-9][0-9]*\2[ACGTN]+|"
    r"\.[ACGTN]+|[ACGTN]+\.)"
)


class VariantInventoryError(RuntimeError):
    """Invalid inventory input, without controlled content in the message."""


class WgsClass(StrEnum):
    SMALL_VARIANT = "small-variant"
    STRUCTURAL = "structural-variant"
    COPY_NUMBER = "copy-number"
    MITOCHONDRIAL = "mitochondrial"
    MOSAIC = "mosaic-candidate"
    REPEAT = "repeat-allele"
    REGULATORY = "regulatory"


class AssessmentStatus(StrEnum):
    ASSESSED = "assessed"
    ABSENT = "absent-in-assessed-scope"
    NOT_ASSESSED = "not-assessed"
    UNSUPPORTED = "unsupported"


class ScreenOutcome(StrEnum):
    CANDIDATE = "candidate-requires-read-validation"
    NO_CANDIDATE = "no-candidate"
    EXCLUDED = "excluded"
    UNEVALUABLE = "unevaluable"


@dataclass(frozen=True, slots=True)
class AlleleBalanceScreenPolicy:
    policy_id: str = field(default="genosapiens.allele-balance-screen", init=False)
    minimum_depth: int = field(default=20, init=False)
    minimum_alternate_depth: int = field(default=3, init=False)
    minimum_balance: float = field(default=0.05, init=False)
    maximum_balance_exclusive: float = field(default=0.20, init=False)
    scope: str = field(default="eligible-supplied-autosomal-biallelic-calls", init=False)
    validation_status: str = field(default="development-not-clinically-validated", init=False)


@dataclass(frozen=True, slots=True, repr=False)
class ScreenCall:
    source_record_id: str
    source_allele_id: str | None
    sample_call: SampleCallQc
    filter_value: str
    technical_flags: tuple[str, ...]
    outcome: ScreenOutcome
    reason: str
    allele_balance: float | None
    read_validation_required: bool


@dataclass(frozen=True, slots=True, repr=False)
class AlleleBalanceScreen:
    enabled: bool
    policy: AlleleBalanceScreenPolicy
    policy_sha256: str
    total_sample_calls: int
    evaluated_calls: int
    excluded_calls: int
    unevaluable_calls: int
    candidate_calls: int
    calls: tuple[ScreenCall, ...]

    def __post_init__(self) -> None:
        counts = Counter(call.outcome for call in self.calls)
        if (
            type(self.enabled) is not bool
            or self.policy_sha256 != _sha256(asdict(self.policy))
            or self.evaluated_calls != counts[ScreenOutcome.CANDIDATE] + counts[ScreenOutcome.NO_CANDIDATE]
            or self.excluded_calls != counts[ScreenOutcome.EXCLUDED]
            or self.unevaluable_calls != counts[ScreenOutcome.UNEVALUABLE]
            or self.candidate_calls != counts[ScreenOutcome.CANDIDATE]
            or (self.enabled and len(self.calls) != self.total_sample_calls)
            or (not self.enabled and self.calls)
            or any(
                call.read_validation_required != (call.outcome is ScreenOutcome.CANDIDATE)
                for call in self.calls
            )
        ):
            raise VariantInventoryError("allele-balance screen accounting is invalid")


@dataclass(frozen=True, slots=True, repr=False)
class ClassInventoryEntry:
    variant_class: WgsClass
    status: AssessmentStatus
    scope: str
    source_record_ids: tuple[str, ...]
    unsupported_record_ids: tuple[str, ...]
    reason: str


@dataclass(frozen=True, slots=True, repr=False)
class VariantClassInventory:
    input_sha256: str
    reference_concordance_sha256: str
    vcf_qc_sha256: str
    total_records: int
    classes: tuple[ClassInventoryEntry, ...]
    unclassified_record_ids: tuple[str, ...]
    non_variant_record_ids: tuple[str, ...]
    allele_balance_screen: AlleleBalanceScreen
    schema: str = field(default=INVENTORY_SCHEMA, init=False)
    classification: ArtifactClassification = field(default=ArtifactClassification.CONTROLLED, init=False)
    llm_policy: LLMPolicy = field(default=LLMPolicy.NO_LLM, init=False)
    genome_wide_assessment: bool = field(default=False, init=False)


def _sha256(value: object) -> str:
    return hashlib.sha256(
        json.dumps(value, sort_keys=True, separators=(",", ":"), ensure_ascii=True).encode("ascii")
    ).hexdigest()


def _record_classes(
    record: VcfRecordQc,
    definitions: dict[tuple[str, str], tuple[str, str]],
    contig_length: int,
) -> tuple[set[WgsClass], set[WgsClass], bool]:
    represented: set[WgsClass] = set()
    unsupported: set[WgsClass] = set()
    unclassified = False
    info = dict(record.class_info_fields)
    svtype = info.get("SVTYPE")
    longtr_hint = "MOTIF" in info and "START" in info and "INEXACT_ALLELE" in info
    repeat_hint = info.get("RU") not in {None, "."} or svtype in {"STR", "VNTR"} or any(
        dict(call.class_format_fields).get("REPCN") not in {None, "."}
        for call in record.sample_calls
    )
    repeat_hint = repeat_hint or longtr_hint
    if longtr_hint:
        valid = (definitions.get(("INFO", "START")) == ("1", "Integer")
                 and definitions.get(("INFO", "END")) == ("1", "Integer")
                 and definitions.get(("INFO", "MOTIF")) == (".", "String")
                 and isinstance(info["START"], str) and info["START"].isdigit()
                 and isinstance(info.get("END"), str) and info["END"].isdigit()
                 and record.position <= int(info["START"]) <= int(info["END"]) <= contig_length
                 and isinstance(info["MOTIF"], str) and re.fullmatch(r"[ACGT]+(?:,[ACGT]+)*", info["MOTIF"]) is not None)
        (represented if valid else unsupported).add(WgsClass.REPEAT)
    cn_hint = any(
        dict(call.class_format_fields).get("CN") not in {None, "."}
        for call in record.sample_calls
    )
    if "SVTYPE" in info and svtype is None:
        unsupported.add(WgsClass.STRUCTURAL)
    if "RU" in info and info["RU"] is None:
        unsupported.add(WgsClass.REPEAT)
    if svtype not in {None, "."}:
        if definitions.get(("INFO", "SVTYPE")) != ("1", "String"):
            unsupported.add(WgsClass.STRUCTURAL)
        elif svtype in SV_TYPES:
            represented.add(WgsClass.STRUCTURAL)
            if svtype in {"DEL", "DUP", "CNV"}:
                represented.add(WgsClass.COPY_NUMBER)
        elif svtype not in {"STR", "VNTR"}:
            unsupported.add(WgsClass.STRUCTURAL)

    for alternate, kind in zip(record.alternates, record.variant_classes, strict=True):
        if alternate in {".", "*", "<NON_REF>", "<*>"} or alternate == record.reference:
            continue
        if REPEAT_SYMBOL.fullmatch(alternate):
            represented.add(WgsClass.REPEAT)
        elif CN_SYMBOL.fullmatch(alternate):
            represented.update({WgsClass.STRUCTURAL, WgsClass.COPY_NUMBER})
        elif alternate.startswith("<") and alternate.endswith(">"):
            base = alternate[1:-1].split(":", 1)[0]
            if base in SV_TYPES:
                represented.add(WgsClass.STRUCTURAL)
                if base in {"DEL", "DUP", "CNV"}:
                    represented.add(WgsClass.COPY_NUMBER)
            else:
                unsupported.add(WgsClass.STRUCTURAL)
                unclassified = True
        elif BREAKEND.fullmatch(alternate):
            represented.add(WgsClass.STRUCTURAL)
        elif kind in SEQUENCE_CLASSES:
            if not repeat_hint and not cn_hint and svtype in {None, "."}:
                represented.add(WgsClass.SMALL_VARIANT)
        else:
            unclassified = True
            if kind is VariantClass.BREAKEND:
                unsupported.add(WgsClass.STRUCTURAL)

    # Numeric CN alone does not establish dosage change, but represents an
    # interval copy-number assessment when both field and interval are usable.
    end = info.get("END")
    valid_interval = (
        definitions.get(("INFO", "END")) == ("1", "Integer")
        and isinstance(end, str) and end.isdigit()
        and record.position <= int(end) <= contig_length
    )
    if cn_hint:
        valid_cn = False
        for call in record.sample_calls:
            value = dict(call.class_format_fields).get("CN")
            if value in {None, "."}:
                continue
            try:
                number = float(value)
            except ValueError:
                continue
            contract = definitions.get(("FORMAT", "CN"))
            if (
                contract in {("1", "Integer"), ("1", "Float")}
                and math.isfinite(number) and number >= 0
                and (contract[1] != "Integer" or value.isdigit())
            ):
                valid_cn = True
        (represented if valid_cn and valid_interval else unsupported).add(WgsClass.COPY_NUMBER)

    if repeat_hint and WgsClass.REPEAT not in represented:
        repeat_unit = info.get("RU")
        valid_repeat = (
            definitions.get(("INFO", "RU")) == ("1", "String")
            and isinstance(repeat_unit, str)
            and re.fullmatch(r"[ACGT]+", repeat_unit) is not None
            and definitions.get(("FORMAT", "REPCN")) == ("1", "String")
            and any(
                isinstance(value := dict(call.class_format_fields).get("REPCN"), str)
                and re.fullmatch(r"[0-9]+(?:/[0-9]+)*", value) is not None
                for call in record.sample_calls
            )
        )
        (represented if valid_repeat else unsupported).add(WgsClass.REPEAT)

    if record.contig in {"M", "chrM"} and represented:
        represented.add(WgsClass.MITOCHONDRIAL)
    return represented, unsupported, unclassified


def _screen_call(
    record: VcfRecordQc,
    call: SampleCallQc,
    represented: set[WgsClass],
    unsupported: set[WgsClass],
    policy: AlleleBalanceScreenPolicy,
) -> ScreenCall:
    contig = record.contig.removeprefix("chr")
    genotype = re.split(r"[/|]", call.genotype) if call.genotype is not None else ()
    balance = None
    outcome = ScreenOutcome.EXCLUDED
    if contig not in AUTOSOMES:
        reason = "not-autosomal"
    elif len(record.alternates) != 1:
        reason = "not-biallelic"
    elif represented != {WgsClass.SMALL_VARIANT} or unsupported:
        reason = "not-small-sequence-variant"
    elif not genotype or "." in genotype:
        reason = "missing-genotype"
    elif len(genotype) != 2 or sorted(genotype) not in (["0", "0"], ["0", "1"]):
        reason = "unsupported-genotype-or-ploidy"
    elif call.depth is None or call.allelic_depths is None or any(
        depth is None for depth in call.allelic_depths
    ):
        outcome, reason = ScreenOutcome.UNEVALUABLE, "missing-ad-or-dp"
    elif call.depth < policy.minimum_depth:
        outcome, reason = ScreenOutcome.UNEVALUABLE, "depth-below-screen-minimum"
    else:
        alternate_depth = call.allelic_depths[1]
        balance = alternate_depth / call.depth
        # Integer comparisons preserve the approved inclusive/exclusive boundaries.
        candidate = (
            alternate_depth >= policy.minimum_alternate_depth
            and alternate_depth * 20 >= call.depth
            and alternate_depth * 5 < call.depth
        )
        outcome = ScreenOutcome.CANDIDATE if candidate else ScreenOutcome.NO_CANDIDATE
        reason = "low-allele-balance-candidate" if candidate else "candidate-thresholds-not-met"
    return ScreenCall(
        source_record_id=record.source_record_id,
        source_allele_id=(
            source_allele_id(record.source_record_id, 1)
            if len(record.alternates) == 1 and record.variant_classes[0] in SEQUENCE_CLASSES
            else None
        ),
        sample_call=call,
        filter_value=record.filter_value,
        technical_flags=technical_flags_for(
            record, call, 1 if len(record.alternates) == 1 else None
        ),
        outcome=outcome,
        reason=reason,
        allele_balance=balance,
        read_validation_required=outcome is ScreenOutcome.CANDIDATE,
    )


@dataclass(frozen=True, slots=True, repr=False)
class IncrementalClassEntry:
    variant_class: WgsClass
    status: AssessmentStatus
    scope: str
    represented_records: int
    unsupported_records: int
    reason: str

    def __post_init__(self) -> None:
        represented = self.represented_records
        unsupported = self.unsupported_records
        expected_status = (
            AssessmentStatus.ASSESSED if represented > 0
            else AssessmentStatus.UNSUPPORTED if unsupported > 0
            else None
        )
        if (
            not isinstance(self.variant_class, WgsClass)
            or not isinstance(self.status, AssessmentStatus)
            or not _summary_text(self.scope)
            or not _summary_count(represented)
            or not _summary_count(unsupported)
            or not _summary_text(self.reason)
            or expected_status is not None and self.status is not expected_status
            or expected_status is None
            and self.status not in {AssessmentStatus.NOT_ASSESSED, AssessmentStatus.ABSENT}
        ):
            raise VariantInventoryError("incremental class entry is invalid")


@dataclass(frozen=True, slots=True, repr=False)
class IncrementalScreenSummary:
    enabled: bool
    policy: AlleleBalanceScreenPolicy
    policy_sha256: str
    total_sample_calls: int
    evaluated_calls: int
    excluded_calls: int
    unevaluable_calls: int
    candidate_calls: int
    outcome_counts: tuple[tuple[str, int], ...]
    reason_counts: tuple[tuple[str, int], ...]

    def __post_init__(self) -> None:
        outcomes = _summary_counts(self.outcome_counts)
        reasons = _summary_counts(self.reason_counts)
        if outcomes is None or reasons is None:
            raise VariantInventoryError("incremental allele-balance summary is invalid")
        candidate = outcomes.get(ScreenOutcome.CANDIDATE.value, 0)
        no_candidate = outcomes.get(ScreenOutcome.NO_CANDIDATE.value, 0)
        if (
            type(self.enabled) is not bool
            or not isinstance(self.policy, AlleleBalanceScreenPolicy)
            or SHA256.fullmatch(self.policy_sha256) is None
            or self.policy_sha256 != _sha256(asdict(self.policy))
            or any(not _summary_count(value) for value in (
                self.total_sample_calls, self.evaluated_calls, self.excluded_calls,
                self.unevaluable_calls, self.candidate_calls,
            ))
            or any(name not in {value.value for value in ScreenOutcome} for name in outcomes)
            or self.evaluated_calls != candidate + no_candidate
            or self.excluded_calls != outcomes.get(ScreenOutcome.EXCLUDED.value, 0)
            or self.unevaluable_calls != outcomes.get(ScreenOutcome.UNEVALUABLE.value, 0)
            or self.candidate_calls != candidate
            or self.enabled and (
                sum(outcomes.values()) != self.total_sample_calls
                or sum(reasons.values()) != self.total_sample_calls
            )
            or not self.enabled and (
                outcomes or reasons or self.evaluated_calls or self.excluded_calls
                or self.unevaluable_calls or self.candidate_calls
            )
        ):
            raise VariantInventoryError("incremental allele-balance summary is invalid")


@dataclass(frozen=True, slots=True, repr=False)
class VariantRecordInventoryFacts:
    record_number: int
    source_record_id: str
    represented_classes: tuple[WgsClass, ...]
    unsupported_classes: tuple[WgsClass, ...]
    unclassified: bool
    non_variant: bool
    screen_calls: tuple[ScreenCall, ...]

    def __post_init__(self) -> None:
        if (
            not isinstance(self.record_number, int)
            or isinstance(self.record_number, bool)
            or self.record_number <= 0
            or SOURCE_RECORD_ID.fullmatch(self.source_record_id) is None
            or self.represented_classes
            != tuple(sorted(set(self.represented_classes), key=lambda value: value.value))
            or self.unsupported_classes
            != tuple(sorted(set(self.unsupported_classes), key=lambda value: value.value))
            or type(self.unclassified) is not bool
            or type(self.non_variant) is not bool
            or self.non_variant != (
                not self.represented_classes
                and not self.unsupported_classes
                and not self.unclassified
            )
            or any(
                not isinstance(call, ScreenCall)
                or call.source_record_id != self.source_record_id
                for call in self.screen_calls
            )
        ):
            raise VariantInventoryError("incremental record inventory facts are invalid")


@dataclass(frozen=True, slots=True, repr=False)
class IncrementalVariantClassSummary:
    input_sha256: str
    reference_concordance_sha256: str
    qc_header_sha256: str
    total_records: int
    classes: tuple[IncrementalClassEntry, ...]
    unclassified_records: int
    non_variant_records: int
    allele_balance_screen: IncrementalScreenSummary
    classification: ArtifactClassification = ArtifactClassification.CONTROLLED
    llm_policy: LLMPolicy = LLMPolicy.NO_LLM
    genome_wide_assessment: bool = False

    def __post_init__(self) -> None:
        expected_classes = tuple(WgsClass)
        actual_classes = tuple(item.variant_class for item in self.classes)
        expected_entries = []
        for item in self.classes:
            if item.represented_records:
                status = AssessmentStatus.ASSESSED
                reason = (
                    "screen-candidates-require-read-validation"
                    if item.variant_class is WgsClass.MOSAIC
                    else "represented-in-supplied-records"
                )
            elif item.unsupported_records:
                status = AssessmentStatus.UNSUPPORTED
                reason = "class-fields-or-representation-unsupported"
            elif (
                item.variant_class is WgsClass.MOSAIC
                and self.allele_balance_screen.evaluated_calls
            ):
                status = AssessmentStatus.ABSENT
                reason = "no-candidate-among-eligible-supplied-calls"
            else:
                status = AssessmentStatus.NOT_ASSESSED
                reason = "no-supported-class-evidence"
            scope = (
                self.allele_balance_screen.policy.scope
                if item.variant_class is WgsClass.MOSAIC
                else "supplied-record-representations-only"
            )
            expected_entries.append((status, scope, reason))
        if (
            SHA256.fullmatch(self.input_sha256) is None
            or SHA256.fullmatch(self.reference_concordance_sha256) is None
            or SHA256.fullmatch(self.qc_header_sha256) is None
            or not isinstance(self.total_records, int)
            or isinstance(self.total_records, bool)
            or self.total_records <= 0
            or actual_classes != expected_classes
            or any(
                (item.status, item.scope, item.reason) != expected
                or item.represented_records > self.total_records
                or item.unsupported_records > self.total_records
                for item, expected in zip(self.classes, expected_entries, strict=True)
            )
            or not _summary_count(self.unclassified_records)
            or self.unclassified_records > self.total_records
            or not _summary_count(self.non_variant_records)
            or self.non_variant_records > self.total_records
            or self.allele_balance_screen.total_sample_calls < self.total_records
            or self.classification is not ArtifactClassification.CONTROLLED
            or self.llm_policy is not LLMPolicy.NO_LLM
            or self.genome_wide_assessment is not False
        ):
            raise VariantInventoryError("incremental variant class summary is invalid")


def _summary_text(value: object) -> bool:
    return (
        isinstance(value, str)
        and bool(value)
        and "\t" not in value
        and "\r" not in value
        and "\n" not in value
    )


def _summary_count(value: object) -> bool:
    return isinstance(value, int) and not isinstance(value, bool) and value >= 0


def _summary_counts(
    values: tuple[tuple[str, int], ...],
) -> dict[str, int] | None:
    if not isinstance(values, tuple) or any(
        not isinstance(item, tuple) or len(item) != 2 for item in values
    ):
        return None
    if any(
        not _summary_text(item[0])
        or not _summary_count(item[1])
        or item[1] <= 0
        for item in values
    ):
        return None
    if values != tuple(sorted(values)) or len({item[0] for item in values}) != len(values):
        return None
    return dict(values)


def _bounded_rows(cursor):
    while True:
        batch = cursor.fetchmany(1024)
        if not batch:
            return
        yield from batch


def classify_record_facts(
    record: VcfRecordQc, *, input_sha256: str, sample_names: tuple[str, ...],
    definitions: dict, lengths: dict[str, int],
    policy: AlleleBalanceScreenPolicy, allele_balance_screen: bool,
) -> VariantRecordInventoryFacts:
    """Pure per-record class projection for workers and the ordered reducer."""
    if not isinstance(record, VcfRecordQc) or type(allele_balance_screen) is not bool:
        raise VariantInventoryError("incremental variant inventory input is invalid")
    expected_id = source_record_id(
        input_sha256,
        SourceSmallVariantRecord(
            record.record_number, record.contig, record.position,
            record.reference, record.alternates,
        ),
    )
    if (
        record.contig not in lengths
        or record.position > lengths.get(record.contig, 0)
        or record.source_record_id != expected_id
        or tuple(call.sample_name for call in record.sample_calls) != sample_names
    ):
        raise VariantInventoryError("incremental variant inventory source order is invalid")
    represented, unsupported, unclassified = _record_classes(
        record, definitions, lengths[record.contig],
    )
    non_variant = not represented and not unsupported and not unclassified
    calls = tuple(
        _screen_call(record, call, represented, unsupported, policy)
        for call in record.sample_calls
    ) if allele_balance_screen else ()
    return VariantRecordInventoryFacts(
        record.record_number, record.source_record_id,
        tuple(sorted(represented, key=lambda value: value.value)),
        tuple(sorted(unsupported, key=lambda value: value.value)),
        unclassified, non_variant, calls,
    )


class IncrementalVariantInventoryAccumulator:
    """Count class/screen facts on disk while returning bounded per-record facts."""

    __slots__ = (
        "header", "expected_records", "allele_balance_screen", "scratch_path",
        "_connection", "_definitions", "_lengths", "_policy", "_records",
        "_finalized", "_closed", "_owns_connection", "_pending_seen",
        "_pending_classes", "_pending_scalars", "_pending_screen",
        "_pending_mosaic",
    )

    def __init__(
        self, header: VcfQcStreamHeader, *, expected_records: int,
        allele_balance_screen: bool, scratch_path: str | Path | None = None,
        connection: sqlite3.Connection | None = None,
        resume: bool = False,
    ) -> None:
        if (
            not isinstance(header, VcfQcStreamHeader)
            or not isinstance(expected_records, int) or isinstance(expected_records, bool)
            or expected_records <= 0 or type(allele_balance_screen) is not bool
            or (connection is None) == (scratch_path is None)
            or type(resume) is not bool or (resume and connection is None)
        ):
            raise VariantInventoryError("incremental variant inventory input is invalid")
        self.header = header
        self.expected_records = expected_records
        self.allele_balance_screen = allele_balance_screen
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
                connection.execute("PRAGMA cache_size=-8192")
                connection.execute("PRAGMA temp_store=FILE")
                connection.execute("PRAGMA mmap_size=0")
                connection.execute("PRAGMA trusted_schema=OFF")
            self._connection = connection
            if not resume:
                self._connection.executescript("""
                CREATE TABLE class_seen(record_number INTEGER PRIMARY KEY, source_record_id TEXT UNIQUE);
                CREATE TABLE class_classes(
                    kind TEXT PRIMARY KEY, represented INTEGER, unsupported INTEGER
                );
                CREATE TABLE class_mosaic_candidates(source_record_id TEXT PRIMARY KEY);
                CREATE TABLE class_scalars(name TEXT PRIMARY KEY, count INTEGER);
                CREATE TABLE class_screen_counts(
                    category TEXT, name TEXT, count INTEGER, PRIMARY KEY(category,name)
                );
                CREATE TABLE class_meta(
                    singleton INTEGER PRIMARY KEY CHECK(singleton=1),
                    input_sha256 TEXT NOT NULL,
                    concordance_sha256 TEXT NOT NULL,
                    header_sha256 TEXT NOT NULL,
                    expected_records INTEGER NOT NULL,
                    allele_balance_screen INTEGER NOT NULL CHECK(allele_balance_screen IN (0,1))
                );
            """)
                self._connection.executemany(
                    "INSERT INTO class_classes VALUES (?,?,?)",
                    ((kind.value, 0, 0) for kind in WgsClass),
                )
                self._connection.executemany(
                    "INSERT INTO class_scalars VALUES (?,0)",
                    ((name,) for name in ("unclassified", "non_variant")),
                )
                self._connection.execute(
                    "INSERT INTO class_meta VALUES (1,?,?,?,?,?)",
                    (header.input_sha256, header.reference_concordance_sha256,
                     header.header_sha256, expected_records, int(allele_balance_screen)),
                )
                if self._owns_connection:
                    self._connection.commit()
                records = 0
            else:
                records = self._connection.execute(
                    "SELECT count(*) FROM class_seen"
                ).fetchone()[0]
                if (type(records) is not int or records < 0 or records > expected_records
                        or self._connection.execute(
                            "SELECT input_sha256,concordance_sha256,header_sha256,"
                            "expected_records,allele_balance_screen FROM class_meta WHERE singleton=1"
                        ).fetchone() != (
                            header.input_sha256, header.reference_concordance_sha256,
                            header.header_sha256, expected_records, int(allele_balance_screen),
                        )
                        or self._connection.execute(
                            "SELECT kind FROM class_classes ORDER BY kind"
                        ).fetchall() != sorted((kind.value,) for kind in WgsClass)):
                    raise VariantInventoryError("incremental variant inventory state is inconsistent")
        except (OSError, sqlite3.Error):
            if descriptor >= 0:
                os.close(descriptor)
            raise VariantInventoryError("cannot create incremental variant inventory") from None
        self._definitions = {
            (namespace, key): (number, field_type)
            for namespace, key, number, field_type in header.class_field_definitions
        }
        self._lengths = {item.contig: item.length for item in header.declared_contigs}
        self._policy = AlleleBalanceScreenPolicy()
        self._records = records
        self._finalized = False
        self._closed = False
        self._pending_seen = []
        self._pending_classes = {}
        self._pending_scalars = Counter()
        self._pending_screen = Counter()
        self._pending_mosaic = set()

    def flush(self) -> None:
        """Insert bounded source uniqueness facts in one coordinator batch."""
        if self._closed or self._finalized:
            raise VariantInventoryError("incremental variant inventory is unavailable")
        if self._pending_seen:
            try:
                self._connection.executemany(
                    "INSERT INTO class_seen VALUES (?,?)", self._pending_seen,
                )
                self._connection.executemany(
                    "UPDATE class_classes SET represented=represented+?, "
                    "unsupported=unsupported+? WHERE kind=?",
                    ((values[0], values[1], kind)
                     for kind, values in self._pending_classes.items()),
                )
                self._connection.executemany(
                    "UPDATE class_scalars SET count=count+? WHERE name=?",
                    ((count, name) for name, count in self._pending_scalars.items()),
                )
                self._connection.executemany("""
                    INSERT INTO class_screen_counts VALUES (?,?,?)
                    ON CONFLICT(category,name) DO UPDATE SET count=count+excluded.count
                """, ((category, name, count)
                      for (category, name), count in self._pending_screen.items()))
                self._connection.executemany(
                    "INSERT OR IGNORE INTO class_mosaic_candidates VALUES (?)",
                    ((source_record_id,) for source_record_id in self._pending_mosaic),
                )
            except sqlite3.Error:
                raise VariantInventoryError("cannot update incremental variant inventory") from None
            self._pending_seen.clear()
            self._pending_classes.clear()
            self._pending_scalars.clear()
            self._pending_screen.clear()
            self._pending_mosaic.clear()

    def facts(self, record: VcfRecordQc) -> VariantRecordInventoryFacts:
        """Compute one record's class facts without mutating reducer state."""
        if self._closed or self._finalized or not isinstance(record, VcfRecordQc):
            raise VariantInventoryError("incremental variant inventory is unavailable")
        if record.record_number != self._records + 1:
            raise VariantInventoryError("incremental variant inventory source order is invalid")
        return classify_record_facts(
            record, input_sha256=self.header.input_sha256,
            sample_names=self.header.sample_names,
            definitions=self._definitions, lengths=self._lengths,
            policy=self._policy,
            allele_balance_screen=self.allele_balance_screen,
        )

    def add_facts(self, facts: VariantRecordInventoryFacts) -> None:
        if (self._closed or self._finalized
                or not isinstance(facts, VariantRecordInventoryFacts)
                or facts.record_number != self._records + 1):
            raise VariantInventoryError("incremental variant inventory source order is invalid")
        for kind in facts.represented_classes:
            self._pending_classes.setdefault(kind.value, [0, 0])[0] += 1
        for kind in facts.unsupported_classes:
            self._pending_classes.setdefault(kind.value, [0, 0])[1] += 1
        if facts.unclassified:
            self._pending_scalars["unclassified"] += 1
        if facts.non_variant:
            self._pending_scalars["non_variant"] += 1
        for call in facts.screen_calls:
            self._pending_screen[("outcome", call.outcome.value)] += 1
            self._pending_screen[("reason", call.reason)] += 1
            if call.outcome is ScreenOutcome.CANDIDATE:
                self._pending_mosaic.add(facts.source_record_id)
        self._records += 1
        self._pending_seen.append((facts.record_number, facts.source_record_id))
        if len(self._pending_seen) >= 1024:
            self.flush()
        if self._owns_connection and self._records % 1024 == 0:
            self._connection.commit()

    def add_projection(self, projection: dict) -> None:
        """Reduce worker-produced class facts without rerunning classification."""
        try:
            if (self._closed or self._finalized or type(projection) is not dict
                    or set(projection) != {
                        "record_number", "source_record_id", "represented_classes",
                        "unsupported_classes", "unclassified", "non_variant",
                        "screen_calls",
                    }
                    or projection["record_number"] != self._records + 1
                    or SOURCE_RECORD_ID.fullmatch(projection["source_record_id"]) is None
                    or type(projection["unclassified"]) is not bool
                    or type(projection["non_variant"]) is not bool):
                raise ValueError
            represented = tuple(WgsClass(value) for value in projection["represented_classes"])
            unsupported = tuple(WgsClass(value) for value in projection["unsupported_classes"])
            if (represented != tuple(sorted(set(represented), key=lambda value: value.value))
                    or unsupported != tuple(sorted(set(unsupported), key=lambda value: value.value))
                    or projection["non_variant"] != (
                        not represented and not unsupported and not projection["unclassified"]
                    )):
                raise ValueError
            calls = projection["screen_calls"]
            if (type(calls) is not list
                    or len(calls) != (len(self.header.sample_names)
                                     if self.allele_balance_screen else 0)):
                raise ValueError
            compact_calls = []
            for sample, call in zip(self.header.sample_names, calls):
                if (type(call) is not dict
                        or call.get("source_record_id") != projection["source_record_id"]
                        or type(call.get("sample_call")) is not dict
                        or call["sample_call"].get("sample_name") != sample
                        or type(call.get("reason")) is not str):
                    raise ValueError
                compact_calls.append((ScreenOutcome(call["outcome"]), call["reason"]))
        except (KeyError, TypeError, ValueError):
            raise VariantInventoryError("source class projection is invalid") from None
        for kind in represented:
            self._pending_classes.setdefault(kind.value, [0, 0])[0] += 1
        for kind in unsupported:
            self._pending_classes.setdefault(kind.value, [0, 0])[1] += 1
        if projection["unclassified"]:
            self._pending_scalars["unclassified"] += 1
        if projection["non_variant"]:
            self._pending_scalars["non_variant"] += 1
        for outcome, reason in compact_calls:
            self._pending_screen[("outcome", outcome.value)] += 1
            self._pending_screen[("reason", reason)] += 1
            if outcome is ScreenOutcome.CANDIDATE:
                self._pending_mosaic.add(projection["source_record_id"])
        self._records += 1
        self._pending_seen.append((projection["record_number"], projection["source_record_id"]))
        if len(self._pending_seen) >= 1024:
            self.flush()

    def add(self, record: VcfRecordQc) -> VariantRecordInventoryFacts:
        facts = self.facts(record)
        self.add_facts(facts)
        return facts

    def finalize(self) -> IncrementalVariantClassSummary:
        if self._closed or self._finalized or self._records != self.expected_records:
            raise VariantInventoryError("incremental variant inventory accounting is invalid")
        self.flush()
        if self._owns_connection:
            self._connection.commit()
        counts = {
            WgsClass(kind): (represented, unsupported)
            for kind, represented, unsupported in _bounded_rows(
                self._connection.execute("SELECT kind,represented,unsupported FROM class_classes")
            )
        }
        mosaic_count = self._connection.execute(
            "SELECT count(*) FROM class_mosaic_candidates"
        ).fetchone()[0]
        counts[WgsClass.MOSAIC] = (mosaic_count, counts[WgsClass.MOSAIC][1])
        outcome_counts = tuple((name, count) for name, count in _bounded_rows(
            self._connection.execute(
                "SELECT name,count FROM class_screen_counts WHERE category='outcome' ORDER BY name"
            )
        ))
        reason_counts = tuple((name, count) for name, count in _bounded_rows(
            self._connection.execute(
                "SELECT name,count FROM class_screen_counts WHERE category='reason' ORDER BY name"
            )
        ))
        outcomes = dict(outcome_counts)
        candidate = outcomes.get(ScreenOutcome.CANDIDATE.value, 0)
        no_candidate = outcomes.get(ScreenOutcome.NO_CANDIDATE.value, 0)
        evaluated = candidate + no_candidate
        entries = []
        for kind in WgsClass:
            represented, unsupported = counts[kind]
            if represented:
                status, reason = AssessmentStatus.ASSESSED, "represented-in-supplied-records"
            elif unsupported:
                status, reason = (
                    AssessmentStatus.UNSUPPORTED,
                    "class-fields-or-representation-unsupported",
                )
            elif kind is WgsClass.MOSAIC and evaluated:
                status, reason = (
                    AssessmentStatus.ABSENT,
                    "no-candidate-among-eligible-supplied-calls",
                )
            else:
                status, reason = AssessmentStatus.NOT_ASSESSED, "no-supported-class-evidence"
            if kind is WgsClass.MOSAIC and represented:
                reason = "screen-candidates-require-read-validation"
            entries.append(IncrementalClassEntry(
                kind, status,
                self._policy.scope if kind is WgsClass.MOSAIC
                else "supplied-record-representations-only",
                represented, unsupported, reason,
            ))
        scalars = dict(_bounded_rows(self._connection.execute(
            "SELECT name,count FROM class_scalars"
        )))
        screen = IncrementalScreenSummary(
            self.allele_balance_screen, self._policy, _sha256(asdict(self._policy)),
            len(self.header.sample_names) * self._records, evaluated,
            outcomes.get(ScreenOutcome.EXCLUDED.value, 0),
            outcomes.get(ScreenOutcome.UNEVALUABLE.value, 0), candidate,
            outcome_counts, reason_counts,
        )
        self._finalized = True
        return IncrementalVariantClassSummary(
            self.header.input_sha256, self.header.reference_concordance_sha256,
            self.header.header_sha256, self._records, tuple(entries),
            scalars["unclassified"], scalars["non_variant"], screen,
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
            raise VariantInventoryError("incremental variant inventory was not finalized")
        return False


def build_variant_class_inventory(
    qc_report: VcfQcReport, *, allele_balance_screen: bool = False
) -> VariantClassInventory:
    """Classify supplied representations, with an optional candidate-only screen.

    Classes overlap (e.g. a symbolic deletion is also copy-number evidence).
    An assessed representation is not a validated patient genotype or a claim
    of comprehensive detection. Empty classes stay not-assessed, not negative.
    """
    if not isinstance(qc_report, VcfQcReport) or type(allele_balance_screen) is not bool:
        raise VariantInventoryError("variant inventory input is invalid")
    represented_ids: dict[WgsClass, list[str]] = {kind: [] for kind in WgsClass}
    unsupported_ids: dict[WgsClass, list[str]] = {kind: [] for kind in WgsClass}
    unclassified_ids, non_variant_ids, screen_calls = [], [], []
    definitions = {(ns, key): (number, field_type) for ns, key, number, field_type in qc_report.class_field_definitions}
    lengths = {item.contig: item.length for item in qc_report.declared_contigs}
    policy = AlleleBalanceScreenPolicy()
    for record in qc_report.records:
        represented, unsupported, unclassified = _record_classes(record, definitions, lengths[record.contig])
        for kind in represented:
            represented_ids[kind].append(record.source_record_id)
        for kind in unsupported:
            unsupported_ids[kind].append(record.source_record_id)
        if unclassified:
            unclassified_ids.append(record.source_record_id)
        if not represented and not unsupported and not unclassified:
            non_variant_ids.append(record.source_record_id)
        if allele_balance_screen:
            screen_calls.extend(
                _screen_call(record, call, represented, unsupported, policy)
                for call in record.sample_calls
            )
    counts = Counter(call.outcome for call in screen_calls)
    screen = AlleleBalanceScreen(
        enabled=allele_balance_screen,
        policy=policy,
        policy_sha256=_sha256(asdict(policy)),
        total_sample_calls=qc_report.sample_count * qc_report.total_records,
        evaluated_calls=counts[ScreenOutcome.CANDIDATE] + counts[ScreenOutcome.NO_CANDIDATE],
        excluded_calls=counts[ScreenOutcome.EXCLUDED],
        unevaluable_calls=counts[ScreenOutcome.UNEVALUABLE],
        candidate_calls=counts[ScreenOutcome.CANDIDATE],
        calls=tuple(screen_calls),
    )
    represented_ids[WgsClass.MOSAIC] = list(dict.fromkeys(
        call.source_record_id for call in screen_calls if call.outcome is ScreenOutcome.CANDIDATE
    ))
    entries = []
    for kind in WgsClass:
        records, unsupported = represented_ids[kind], unsupported_ids[kind]
        if records:
            status, reason = AssessmentStatus.ASSESSED, "represented-in-supplied-records"
        elif unsupported:
            status, reason = AssessmentStatus.UNSUPPORTED, "class-fields-or-representation-unsupported"
        elif kind is WgsClass.MOSAIC and screen.evaluated_calls:
            status, reason = AssessmentStatus.ABSENT, "no-candidate-among-eligible-supplied-calls"
        else:
            status, reason = AssessmentStatus.NOT_ASSESSED, "no-supported-class-evidence"
        if kind is WgsClass.MOSAIC and records:
            reason = "screen-candidates-require-read-validation"
        entries.append(ClassInventoryEntry(
            variant_class=kind, status=status,
            scope=policy.scope if kind is WgsClass.MOSAIC else "supplied-record-representations-only",
            source_record_ids=tuple(records), unsupported_record_ids=tuple(unsupported), reason=reason,
        ))
    return VariantClassInventory(
        input_sha256=qc_report.input_sha256,
        reference_concordance_sha256=qc_report.reference_concordance_sha256,
        vcf_qc_sha256=hashlib.sha256(vcf_qc_as_json(qc_report).encode("ascii")).hexdigest(),
        total_records=qc_report.total_records, classes=tuple(entries),
        unclassified_record_ids=tuple(unclassified_ids), non_variant_record_ids=tuple(non_variant_ids),
        allele_balance_screen=screen,
    )


def variant_inventory_as_json(inventory: VariantClassInventory) -> str:
    return json.dumps(asdict(inventory), indent=2, sort_keys=True, ensure_ascii=True) + "\n"
