"""Bounded immutable storage for WGS records; no scientific eligibility policy.

Canonical shards own scientific data; stage SQLite state owns committed cursors
and reducers, while derived query indexes may be disposable. All explicit
serialization is controlled/no-LLM, including the state and journal. These
helpers do not issue preparation, release, or controlled-execution authority.
"""

from __future__ import annotations

from collections.abc import Iterable, Iterator
from contextlib import contextmanager
from dataclasses import asdict, dataclass, field, fields
import fcntl
import gzip
import hashlib
import json
import os
from pathlib import Path
import re
import sqlite3
import stat
import tempfile
from typing import ClassVar
import uuid
import zlib
from genosapiens.vcf_contracts import POLICY_SHA256 as VCF_POLICY_SHA256


SCHEMA = "genosapiens.wgs-ledger"
VERSION = 1
MAX_ROW_BYTES = 4 * 1024 * 1024
MAX_SHARD_BYTES = 64 * 1024 * 1024
MAX_SHARD_ROWS = 1024
MAX_MANIFEST_BYTES = 8 * 1024 * 1024
MAX_SHARDS = 16384
MAX_JOURNAL_EVENTS = MAX_SHARDS * 8 + 1
_SHA = re.compile(r"[0-9a-f]{64}")
_NAME = re.compile(r"[a-z0-9][a-z0-9._-]{0,127}")
_SCHEMA = re.compile(r"genosapiens\.[a-z0-9][a-z0-9.-]{0,95}")
_SOURCE_RECORD = re.compile(r"source-record-[0-9a-f]{24}")
_SOURCE_ALLELE = re.compile(r"source-allele-[0-9a-f]{24}")
_NORMALIZED_ALLELE = re.compile(r"allele-[0-9a-f]{24}")
_IDENTIFIER = re.compile(r"[A-Za-z0-9][A-Za-z0-9._:/-]{0,127}")

SOURCE_ROW_SCHEMA = "genosapiens.wgs-source-row"
SOURCE_BINDING_SCHEMA = "genosapiens.wgs-source-binding"
ARTIFACT_SCHEMA = "genosapiens.wgs-artifact"
PARENT_IDENTITY_SCHEMA = "genosapiens.wgs-parent-identity"
COMPRESSION_SCHEMA = "genosapiens.wgs-compression"
SHARD_SCHEMA = "genosapiens.wgs-shard"
SOURCE_SUMMARY_SCHEMA = "genosapiens.wgs-source-summary"
SOURCE_MANIFEST_SCHEMA = "genosapiens.wgs-source-manifest"
FILTER_DECISION_SCHEMA = "genosapiens.wgs-filter-decision"
FILTER_SUMMARY_SCHEMA = "genosapiens.wgs-filter-summary"
FILTER_MANIFEST_SCHEMA = "genosapiens.wgs-filter-manifest"
WORKING_SET_ROW_SCHEMA = "genosapiens.wgs-working-set-row"
WORKING_SET_SUMMARY_SCHEMA = "genosapiens.wgs-working-set-summary"
WORKING_SET_MANIFEST_SCHEMA = "genosapiens.wgs-working-set-manifest"
JOURNAL_ENTRY_SCHEMA = "genosapiens.wgs-journal-entry"
STAGE_RUN_BINDING_SCHEMA = "genosapiens.wgs-stage-run-binding"
PREPARED_CASE_SCHEMA = "genosapiens.wgs-prepared-case"
INDEX_TABLE_BINDING_SCHEMA = "genosapiens.wgs-index-table-binding"
INDEX_CLOSURE_SCHEMA = "genosapiens.wgs-index-closure"
NORMALIZATION_INDEX_ROW_SCHEMA = "genosapiens.wgs-index-normalization-row"
DECISION_INDEX_ROW_SCHEMA = "genosapiens.wgs-index-decision-row"
WORKING_IDENTITY_INDEX_ROW_SCHEMA = "genosapiens.wgs-index-working-identity-row"
RESPONSE_INDEX_ROW_SCHEMA = "genosapiens.wgs-index-response-row"
CONTEXT_INDEX_ROW_SCHEMA = "genosapiens.wgs-index-context-row"
HEADER_SCHEMA = "genosapiens.wgs-header"
DIAGNOSTIC_RECEIPT_SCHEMA = "genosapiens.wgs-preparation-diagnostics"
INVENTORY_SCHEMA = "genosapiens.wgs-preparation-inventory"

PREPARATION_DISPOSITIONS = ("normalized", "native_class", "quarantined", "non_variant")
EARLY_FILTER_DISPOSITIONS = ("retain", "retain_unknown", "exclude_common")
SCOPE_DISPOSITIONS = ('in_scope', 'outside_scope', 'unresolved_scope', 'not_applicable')
PREPARATION_ANALYSIS_OUTCOMES = (
    "pending", "excluded_common_by_policy", "quarantined", "non_variant",
    "excluded_outside_scope_by_policy", "deferred_unresolved_scope",
)
PREPARED_SUMMARY_ROLES = frozenset(
    {"qc_summary", "normalization_summary", "class_summary", "filter_summary"}
)
PREPARED_RESOURCE_ROLES = frozenset(
    {"population_manifest", "clinical_manifest", "exception_manifest", "precision_manifest",
     "early_scope_policy", "early_scope_manifest"}
)


class LedgerError(ValueError):
    """Static error messages only; row values and paths must not be disclosed."""


def _require(condition):
    if not condition:
        raise LedgerError("invalid WGS ledger contract")


def _sha(value):
    return hashlib.sha256(value).hexdigest()


def _digest(value):
    return type(value) is str and _SHA.fullmatch(value) is not None


def _count(value, minimum=0):
    return type(value) is int and minimum <= value <= 2**63 - 1


def _name(value):
    return type(value) is str and _NAME.fullmatch(value) is not None


def _schema(value):
    return type(value) is str and _SCHEMA.fullmatch(value) is not None


def _identifier(value):
    return type(value) is str and _IDENTIFIER.fullmatch(value) is not None


def _strings(values, *, maximum=256, pattern=_IDENTIFIER, allow_empty=True):
    return (type(values) is tuple and (allow_empty or bool(values)) and len(values) <= maximum
            and len(values) == len(set(values))
            and all(type(value) is str and pattern.fullmatch(value) is not None for value in values))


def _pairs(items):
    result = {}
    for key, value in items:
        _require(key not in result)
        result[key] = value
    return result


def _constant(_value):
    raise LedgerError("nonfinite WGS ledger value")


def canonical(value) -> bytes:
    try:
        return (json.dumps(value, sort_keys=True, separators=(",", ":"),
                           ensure_ascii=True, allow_nan=False) + "\n").encode("ascii")
    except (ValueError, TypeError, UnicodeError, RecursionError, OverflowError):
        raise LedgerError("cannot encode WGS ledger value") from None


def _json(content: bytes, limit=MAX_MANIFEST_BYTES):
    try:
        _require(type(content) is bytes and len(content) <= limit)
        value = json.loads(content, object_pairs_hook=_pairs, parse_constant=_constant)
        _require(canonical(value) == content)
        return value
    except (ValueError, TypeError, UnicodeError, RecursionError, OverflowError):
        raise LedgerError("invalid canonical WGS ledger bytes") from None


class _Controlled:
    # Dataclass field traversal and equality ignore this per-instance owner cache.
    __slots__ = ("_identity", "_publication")
    classification: ClassVar[str] = "controlled"
    llm_policy: ClassVar[str] = "no-llm"
    cache_identity: ClassVar[bool] = False
    schema: ClassVar[str]
    version: ClassVar[int] = 1
    discriminator: ClassVar[str]

    def __repr__(self):
        return "<controlled WGS ledger; no-llm>"

    def as_bytes(self):
        if self.cache_identity:
            cached = getattr(self, "_identity", None)
            if cached is not None:
                return cached[0]
        content = canonical({"schema": self.schema, "version": self.version,
                             "kind": self.discriminator, "classification": "controlled",
                             "llm_policy": "no-llm", "data": asdict(self)})
        _require(len(content) <= MAX_MANIFEST_BYTES)
        if self.cache_identity:
            object.__setattr__(self, "_identity", (content, _sha(content)))
        return content

    @property
    def sha256(self):
        if self.cache_identity:
            cached = getattr(self, "_identity", None)
            if cached is not None:
                return cached[1]
            self.as_bytes()
            return self._identity[1]
        return _sha(self.as_bytes())


@dataclass(frozen=True, slots=True, repr=False)
class LedgerRow(_Controlled):
    schema: ClassVar[str] = SOURCE_ROW_SCHEMA
    discriminator: ClassVar[str] = "source_row"
    record_number: int
    source_record_id: str
    source_allele_ids: tuple[str, ...]
    biological_alternate_count: int
    # Canonical JSON text freezes nested values without materializing a file.
    payload_json: str

    def __post_init__(self):
        _require(_count(self.record_number, 1))
        _require(type(self.source_record_id) is str and _SOURCE_RECORD.fullmatch(self.source_record_id))
        _require(type(self.source_allele_ids) is tuple and 0 < len(self.source_allele_ids) <= 65536)
        _require(len(set(self.source_allele_ids)) == len(self.source_allele_ids))
        _require(all(type(v) is str and _SOURCE_ALLELE.fullmatch(v) for v in self.source_allele_ids))
        _require(_count(self.biological_alternate_count)
                 and self.biological_alternate_count <= len(self.source_allele_ids))
        _require(type(self.payload_json) is str and len(self.payload_json) <= MAX_ROW_BYTES)
        try:
            payload = _json(self.payload_json.encode("ascii"), MAX_ROW_BYTES)
        except UnicodeError:
            raise LedgerError("invalid WGS ledger payload") from None
        _require(type(payload) is dict)

    def as_bytes(self):
        content = super(LedgerRow, self).as_bytes()
        _require(len(content) <= MAX_ROW_BYTES)
        return content


@dataclass(frozen=True, slots=True, repr=False)
class EncodedLedgerRow:
    """One bounded encoding with the source ordinal and accounting facts."""

    row: LedgerRow
    content: bytes
    record_number: int
    physical_alternate_slots: int
    biological_alternates: int


def encode_ledger_row(row: LedgerRow) -> EncodedLedgerRow:
    _require(type(row) is LedgerRow)
    content = row.as_bytes()
    return EncodedLedgerRow(
        row, content, row.record_number, len(row.source_allele_ids),
        row.biological_alternate_count,
    )


@dataclass(frozen=True, slots=True, repr=False)
class SourceBinding(_Controlled):
    schema: ClassVar[str] = SOURCE_BINDING_SCHEMA
    discriminator: ClassVar[str] = "source_binding"
    input_sha256: str
    input_size_bytes: int
    case_manifest_sha256: str
    selected_sample_sha256: str
    concordance_sha256: str
    reference_manifest_sha256: str
    executable_stack_sha256: str

    vcf_evidence_policy_sha256: str = VCF_POLICY_SHA256

    def __post_init__(self):
        _require(self.vcf_evidence_policy_sha256 == VCF_POLICY_SHA256)
        _require(all(_digest(v) for v in (
            self.input_sha256, self.case_manifest_sha256, self.selected_sample_sha256,
            self.concordance_sha256, self.reference_manifest_sha256,
            self.executable_stack_sha256,
        )))
        _require(_count(self.input_size_bytes, 1))


@dataclass(frozen=True, slots=True, repr=False)
class Artifact(_Controlled):
    schema: ClassVar[str] = ARTIFACT_SCHEMA
    discriminator: ClassVar[str] = "artifact"
    name: str
    sha256: str = field()
    size_bytes: int

    def __post_init__(self):
        _require(_name(self.name) and _digest(self.sha256) and _count(self.size_bytes, 1))


@dataclass(frozen=True, slots=True, repr=False)
class ParentIdentity(_Controlled):
    schema: ClassVar[str] = PARENT_IDENTITY_SCHEMA
    discriminator: ClassVar[str] = "parent_identity"
    role: str
    parent_schema: str
    parent_version: int
    sha256: str = field()

    def __post_init__(self):
        _require(_name(self.role) and _schema(self.parent_schema)
                 and _count(self.parent_version, 1) and self.parent_version <= 2**31 - 1
                 and _digest(self.sha256))


def _parents(values, *, required=frozenset()):
    _require(type(values) is tuple and 0 < len(values) <= 64)
    _require(all(type(value) is ParentIdentity for value in values))
    roles = tuple(value.role for value in values)
    _require(len(roles) == len(set(roles)) and required.issubset(roles))


@dataclass(frozen=True, slots=True, repr=False)
class CompressionMetadata(_Controlled):
    schema: ClassVar[str] = COMPRESSION_SCHEMA
    discriminator: ClassVar[str] = "compression"
    algorithm: str
    level: int
    mtime: int
    filename: str
    text_encoding: str
    final_newline: bool

    def __post_init__(self):
        _require(self.algorithm == "gzip" and type(self.level) is int and self.level == 6
                 and type(self.mtime) is int and self.mtime == 0 and self.filename == ""
                 and self.text_encoding == "ascii" and self.final_newline is True)


GZIP_JSONL = CompressionMetadata("gzip", 6, 0, "", "ascii", True)


@dataclass(frozen=True, slots=True, repr=False)
class Shard(_Controlled):
    schema: ClassVar[str] = SHARD_SCHEMA
    discriminator: ClassVar[str] = "shard"
    role: str
    row_schema: str
    row_version: int
    ordinal: int
    artifact: Artifact
    compression: CompressionMetadata
    uncompressed_sha256: str
    uncompressed_bytes: int
    first_record: int
    last_record: int
    rows: int
    physical_alternate_slots: int
    biological_alternates: int
    previous_shard_sha256: str | None

    def __post_init__(self):
        _require(self.role in {"source", "filter", "working_set"}
                 and _schema(self.row_schema) and _count(self.row_version, 1)
                 and self.row_version <= 2**31 - 1)
        _require(_count(self.ordinal, 1) and self.ordinal <= MAX_SHARDS)
        _require(type(self.artifact) is Artifact and self.compression == GZIP_JSONL
                 and _digest(self.uncompressed_sha256))
        _require(_count(self.uncompressed_bytes, 1) and self.uncompressed_bytes <= MAX_SHARD_BYTES)
        _require(_count(self.first_record, 1) and _count(self.last_record, 1))
        _require(_count(self.rows, 1) and self.rows <= MAX_SHARD_ROWS)
        _require(self.last_record - self.first_record + 1 == self.rows)
        _require(_count(self.physical_alternate_slots)
                 and _count(self.biological_alternates)
                 and self.physical_alternate_slots >= self.rows
                 and self.biological_alternates <= self.physical_alternate_slots)
        _require(self.previous_shard_sha256 is None if self.ordinal == 1
                 else _digest(self.previous_shard_sha256))
        prefix = {"source": "shard", "filter": "decisions", "working_set": "retained"}[self.role]
        _require(self.artifact.name == f"{prefix}-{self.ordinal:06d}.jsonl.gz")
        _require(self.artifact.size_bytes <= MAX_SHARD_BYTES + 65536)


@dataclass(frozen=True, slots=True, repr=False)
class LedgerSummary(_Controlled):
    schema: ClassVar[str] = SOURCE_SUMMARY_SCHEMA
    discriminator: ClassVar[str] = "source_summary"
    total_records: int
    total_physical_alternate_slots: int
    total_biological_alternates: int
    # Individual stage summaries are separately bounded and authenticated.
    reports: tuple[Artifact, ...]

    def __post_init__(self):
        _require(_count(self.total_records, 1)
                 and _count(self.total_physical_alternate_slots, self.total_records)
                 and _count(self.total_biological_alternates)
                 and self.total_biological_alternates <= self.total_physical_alternate_slots)
        _artifacts(self.reports)


def _artifacts(values):
    _require(type(values) is tuple and len(values) <= 64)
    _require(all(type(v) is Artifact for v in values))
    _require(len({v.name for v in values}) == len(values))


@dataclass(frozen=True, slots=True, repr=False)
class LedgerManifest(_Controlled):
    schema: ClassVar[str] = SOURCE_MANIFEST_SCHEMA
    discriminator: ClassVar[str] = "source_manifest"
    cache_identity: ClassVar[bool] = True
    source: SourceBinding
    header: Artifact
    parents: tuple[ParentIdentity, ...]
    summary: LedgerSummary
    shards: tuple[Shard, ...]

    def __post_init__(self):
        _require(type(self.source) is SourceBinding and type(self.header) is Artifact)
        _parents(self.parents, required={
            "original_input", "case_manifest", "selected_sample", "complete_concordance",
            "reference_manifest", "execution_limits", "executable_stack",
        })
        parents = {value.role: value.sha256 for value in self.parents}
        schemas = {value.role: value.parent_schema for value in self.parents}
        _require(parents["original_input"] == self.source.input_sha256
                 and parents["case_manifest"] == self.source.case_manifest_sha256
                 and parents["selected_sample"] == self.source.selected_sample_sha256
                 and parents["complete_concordance"] == self.source.concordance_sha256
                 and parents["reference_manifest"] == self.source.reference_manifest_sha256
                 and parents["executable_stack"] == self.source.executable_stack_sha256
                 and schemas["original_input"] == "genosapiens.original-input"
                 and schemas["selected_sample"] == "genosapiens.selected-sample-identity")
        _require(type(self.summary) is LedgerSummary and type(self.shards) is tuple)
        _require(0 < len(self.shards) <= MAX_SHARDS)
        next_record, previous, physical, biological = 1, None, 0, 0
        for ordinal, shard in enumerate(self.shards, 1):
            _require(type(shard) is Shard and shard.role == "source"
                     and shard.row_schema == SOURCE_ROW_SCHEMA and shard.row_version == 1
                     and shard.ordinal == ordinal)
            _require(shard.first_record == next_record and shard.previous_shard_sha256 == previous)
            next_record, previous = shard.last_record + 1, shard.sha256
            physical += shard.physical_alternate_slots
            biological += shard.biological_alternates
        _require(next_record - 1 == self.summary.total_records and
                 physical == self.summary.total_physical_alternate_slots and
                 biological == self.summary.total_biological_alternates)
        names = (self.header.name, *(v.name for v in self.summary.reports),
                 *(v.artifact.name for v in self.shards))
        _require(len(names) == len(set(names)))


@dataclass(frozen=True, slots=True, repr=False)
class NamedCount(_Controlled):
    schema: ClassVar[str] = "genosapiens.wgs-named-count"
    discriminator: ClassVar[str] = "named_count"
    name: str
    count: int

    def __post_init__(self):
        _require(_name(self.name) and _count(self.count))


def _named_counts(values, *, exact_names=None):
    _require(type(values) is tuple and len(values) <= 256
             and all(type(value) is NamedCount for value in values))
    names = tuple(value.name for value in values)
    _require(len(names) == len(set(names)))
    if exact_names is not None:
        _require(names == tuple(exact_names))


@dataclass(frozen=True, slots=True, repr=False)
class FilterDecision(_Controlled):
    version: ClassVar[int] = 2
    schema: ClassVar[str] = FILTER_DECISION_SCHEMA
    discriminator: ClassVar[str] = "filter_decision"
    record_number: int
    allele_number: int
    source_allele_id: str
    normalized_allele_id: str | None
    preparation_disposition: str
    early_filter_disposition: str | None
    bypass_reason: str | None
    reason_codes: tuple[str, ...]
    parent_source_manifest_sha256: str
    policy_sha256: str
    evidence_json: str
    scope_disposition: str = field(kw_only=True)
    scope_policy_sha256: str = field(kw_only=True)
    scope_manifest_sha256: str = field(kw_only=True)

    def __post_init__(self):
        _require(_count(self.record_number, 1) and _count(self.allele_number, 1)
                 and type(self.source_allele_id) is str
                 and _SOURCE_ALLELE.fullmatch(self.source_allele_id)
                 and self.preparation_disposition in PREPARATION_DISPOSITIONS
                 and _strings(self.reason_codes, allow_empty=True)
                 and _digest(self.parent_source_manifest_sha256) and _digest(self.policy_sha256)
                 and _digest(self.scope_policy_sha256) and _digest(self.scope_manifest_sha256)
                 and self.scope_disposition in SCOPE_DISPOSITIONS
                 and (self.scope_disposition == 'not_applicable') ==
                     (self.preparation_disposition in ('non_variant','quarantined')))
        if self.preparation_disposition == "normalized":
            _require(type(self.normalized_allele_id) is str
                     and _NORMALIZED_ALLELE.fullmatch(self.normalized_allele_id)
                     and (self.early_filter_disposition in EARLY_FILTER_DISPOSITIONS
                          if self.scope_disposition == 'in_scope' else self.early_filter_disposition is None)
                     and self.bypass_reason is None)
        else:
            _require(self.normalized_allele_id is None and self.early_filter_disposition is None
                     and _name(self.bypass_reason))
        if self.early_filter_disposition in {"retain_unknown", "exclude_common"}:
            _require(bool(self.reason_codes))
        _require(type(self.evidence_json) is str and len(self.evidence_json) <= MAX_ROW_BYTES)
        try:
            evidence = _json(self.evidence_json.encode("ascii"), MAX_ROW_BYTES)
        except UnicodeError:
            raise LedgerError("invalid WGS ledger payload") from None
        _require(type(evidence) is dict and len(self.as_bytes()) <= MAX_ROW_BYTES)


@dataclass(frozen=True, slots=True, repr=False)
class FilterSummary(_Controlled):
    version: ClassVar[int] = 2
    schema: ClassVar[str] = FILTER_SUMMARY_SCHEMA
    discriminator: ClassVar[str] = "filter_summary"
    total_source_occurrences: int
    unique_normalized_alleles: int
    preparation_counts: tuple[NamedCount, ...]
    early_filter_counts: tuple[NamedCount, ...]
    reason_counts: tuple[NamedCount, ...]
    scope_counts: tuple[NamedCount, ...] = field(kw_only=True)
    in_scope_normalized_occurrences: int = field(kw_only=True)
    population_request_keys: int = field(kw_only=True)
    clinical_query_keys: int = field(kw_only=True)

    def __post_init__(self):
        _require(_count(self.total_source_occurrences, 1)
                 and _count(self.unique_normalized_alleles))
        _named_counts(self.preparation_counts, exact_names=PREPARATION_DISPOSITIONS)
        _named_counts(self.early_filter_counts, exact_names=EARLY_FILTER_DISPOSITIONS)
        _named_counts(self.reason_counts)
        _named_counts(self.scope_counts, exact_names=SCOPE_DISPOSITIONS)
        preparation = {value.name: value.count for value in self.preparation_counts}
        _require(sum(preparation.values()) == self.total_source_occurrences)
        scope = {value.name:value.count for value in self.scope_counts}
        _require(sum(scope.values()) == self.total_source_occurrences
            and scope['not_applicable'] == preparation['quarantined'] + preparation['non_variant']
            and _count(self.in_scope_normalized_occurrences)
            and self.in_scope_normalized_occurrences <= min(preparation['normalized'],scope['in_scope'])
            and sum(value.count for value in self.early_filter_counts) == self.in_scope_normalized_occurrences
            and _count(self.population_request_keys) and _count(self.clinical_query_keys)
            and self.population_request_keys <= self.in_scope_normalized_occurrences
            and self.clinical_query_keys <= self.in_scope_normalized_occurrences)


@dataclass(frozen=True, slots=True, repr=False)
class FilterManifest(_Controlled):
    version: ClassVar[int] = 2
    schema: ClassVar[str] = FILTER_MANIFEST_SCHEMA
    discriminator: ClassVar[str] = "filter_manifest"
    cache_identity: ClassVar[bool] = True
    parent_source_manifest_sha256: str
    policy_sha256: str
    parents: tuple[ParentIdentity, ...]
    summary: FilterSummary
    shards: tuple[Shard, ...]
    scope_policy_sha256: str = field(kw_only=True)
    scope_manifest_sha256: str = field(kw_only=True)

    def __post_init__(self):
        _require(_digest(self.parent_source_manifest_sha256) and _digest(self.policy_sha256)
            and _digest(self.scope_policy_sha256) and _digest(self.scope_manifest_sha256))
        _parents(self.parents, required={
            "source_manifest", "prefilter_policy", "population_manifest",
            "clinical_manifest", "exception_manifest", "precision_manifest",
            "early_scope_policy", "early_scope_manifest",
        })
        parents = {value.role: value for value in self.parents}
        _require(parents["source_manifest"].sha256 == self.parent_source_manifest_sha256
                 and parents["source_manifest"].parent_schema == SOURCE_MANIFEST_SCHEMA
                 and parents["source_manifest"].parent_version == 1
                 and parents["prefilter_policy"].sha256 == self.policy_sha256
                 and parents["prefilter_policy"].parent_schema
                 == "genosapiens.common-allele-prefilter-policy"
                 and parents["prefilter_policy"].parent_version == 1)
        for role, digest, schema in (
            ('early_scope_policy', self.scope_policy_sha256, 'genosapiens.track1-early-scope-policy'),
            ('early_scope_manifest', self.scope_manifest_sha256, 'genosapiens.track1-early-scope-manifest')):
            _require(parents[role].sha256 == digest and parents[role].parent_schema == schema
                and parents[role].parent_version == 1)
        _require(type(self.summary) is FilterSummary and type(self.shards) is tuple
                 and 0 < len(self.shards) <= MAX_SHARDS)
        _validate_shard_chain(self.shards, "filter", FILTER_DECISION_SCHEMA)
        _require(sum(shard.rows for shard in self.shards) == self.summary.total_source_occurrences)


@dataclass(frozen=True, slots=True, repr=False)
class WorkingSetRow(_Controlled):
    schema: ClassVar[str] = WORKING_SET_ROW_SCHEMA
    discriminator: ClassVar[str] = "working_set_row"
    ordinal: int
    identity_id: str
    route: str
    source_allele_ids: tuple[str, ...]
    parent_filter_manifest_sha256: str
    payload_json: str

    def __post_init__(self):
        _require(_count(self.ordinal, 1) and self.route in {"sequence", "native"}
                 and _strings(self.source_allele_ids, maximum=65536,
                              pattern=_SOURCE_ALLELE, allow_empty=False)
                 and _digest(self.parent_filter_manifest_sha256))
        _require((self.route == "sequence" and type(self.identity_id) is str
                  and _NORMALIZED_ALLELE.fullmatch(self.identity_id) is not None)
                 or (self.route == "native" and _identifier(self.identity_id)))
        _require(type(self.payload_json) is str and len(self.payload_json) <= MAX_ROW_BYTES)
        try:
            payload = _json(self.payload_json.encode("ascii"), MAX_ROW_BYTES)
        except UnicodeError:
            raise LedgerError("invalid WGS ledger payload") from None
        _require(type(payload) is dict and len(self.as_bytes()) <= MAX_ROW_BYTES)


@dataclass(frozen=True, slots=True, repr=False)
class WorkingSetSummary(_Controlled):
    version: ClassVar[int] = 2
    schema: ClassVar[str] = WORKING_SET_SUMMARY_SCHEMA
    discriminator: ClassVar[str] = "working_set_summary"
    total_source_occurrences: int
    retained_occurrences: int
    retained_unique_identities: int
    excluded_common_occurrences: int
    bypassed_occurrences: int
    analysis_outcome_counts: tuple[NamedCount, ...]
    excluded_scope_occurrences: int = field(kw_only=True)
    unresolved_scope_occurrences: int = field(kw_only=True)

    def __post_init__(self):
        _require(_count(self.total_source_occurrences, 1)
                 and _count(self.retained_occurrences)
                 and _count(self.retained_unique_identities)
                 and _count(self.excluded_common_occurrences)
                 and _count(self.bypassed_occurrences)
                 and _count(self.excluded_scope_occurrences)
                 and _count(self.unresolved_scope_occurrences)
                 and self.retained_unique_identities <= self.retained_occurrences
                 and self.retained_occurrences + self.excluded_common_occurrences
                 + self.excluded_scope_occurrences + self.unresolved_scope_occurrences
                 + self.bypassed_occurrences == self.total_source_occurrences)
        _named_counts(
            self.analysis_outcome_counts,
            exact_names=PREPARATION_ANALYSIS_OUTCOMES,
        )
        outcomes = {value.name: value.count for value in self.analysis_outcome_counts}
        _require(
            outcomes["pending"] == self.retained_occurrences
            and outcomes["excluded_common_by_policy"]
            == self.excluded_common_occurrences
            and outcomes['excluded_outside_scope_by_policy'] == self.excluded_scope_occurrences
            and outcomes['deferred_unresolved_scope'] == self.unresolved_scope_occurrences
            and outcomes["quarantined"] + outcomes["non_variant"]
            == self.bypassed_occurrences
            and sum(outcomes.values()) == self.total_source_occurrences
        )


@dataclass(frozen=True, slots=True, repr=False)
class WorkingSetManifest(_Controlled):
    version: ClassVar[int] = 2
    schema: ClassVar[str] = WORKING_SET_MANIFEST_SCHEMA
    discriminator: ClassVar[str] = "working_set_manifest"
    cache_identity: ClassVar[bool] = True
    parent_filter_manifest_sha256: str
    parents: tuple[ParentIdentity, ...]
    summary: WorkingSetSummary
    shards: tuple[Shard, ...]

    def __post_init__(self):
        _require(_digest(self.parent_filter_manifest_sha256))
        _parents(self.parents, required={"source_manifest", "filter_manifest"})
        parents = {value.role: value for value in self.parents}
        _require(parents["source_manifest"].parent_schema == SOURCE_MANIFEST_SCHEMA
                 and parents["source_manifest"].parent_version == 1
                 and parents["filter_manifest"].parent_schema == FILTER_MANIFEST_SCHEMA
                 and parents["filter_manifest"].parent_version == 2
                 and parents["filter_manifest"].sha256 == self.parent_filter_manifest_sha256)
        _require(type(self.summary) is WorkingSetSummary and type(self.shards) is tuple)
        if self.shards:
            _require(len(self.shards) <= MAX_SHARDS)
            _validate_shard_chain(self.shards, "working_set", WORKING_SET_ROW_SCHEMA)
            _require(sum(shard.physical_alternate_slots for shard in self.shards)
                     == self.summary.retained_occurrences)
        else:
            _require(self.summary.retained_occurrences == 0
                     and self.summary.retained_unique_identities == 0)


def _validate_shard_chain(shards, role, row_schema):
    next_record, previous = 1, None
    for ordinal, shard in enumerate(shards, 1):
        _require(type(shard) is Shard and shard.role == role
                 and shard.row_schema == row_schema
                 and shard.row_version == (2 if role == 'filter' else 1)
                 and shard.ordinal == ordinal and shard.first_record == next_record
                 and shard.previous_shard_sha256 == previous)
        next_record, previous = shard.last_record + 1, shard.sha256


@dataclass(frozen=True, slots=True, repr=False)
class IndexTableBinding(_Controlled):
    schema: ClassVar[str] = INDEX_TABLE_BINDING_SCHEMA
    discriminator: ClassVar[str] = "index_table_binding"
    parent_sha256: str | None
    rows: int
    logical_sha256: str

    def __post_init__(self):
        _require(_count(self.rows) and _digest(self.logical_sha256))
        _require(_digest(self.parent_sha256) if self.rows else
                 self.parent_sha256 is None or _digest(self.parent_sha256))
        if self.rows == 0:
            _require(self.logical_sha256 == hashlib.sha256(b"").hexdigest())


@dataclass(frozen=True, slots=True, repr=False)
class IndexClosure(_Controlled):
    schema: ClassVar[str] = INDEX_CLOSURE_SCHEMA
    discriminator: ClassVar[str] = "index_closure"
    source_manifest_sha256: str
    normalization: IndexTableBinding
    decisions: IndexTableBinding
    working_identities: IndexTableBinding
    responses: IndexTableBinding
    context_memberships: IndexTableBinding

    def __post_init__(self):
        _require(_digest(self.source_manifest_sha256))
        _require(all(type(value) is IndexTableBinding for value in (
            self.normalization, self.decisions, self.working_identities,
            self.responses, self.context_memberships,
        )))


@dataclass(frozen=True, slots=True, repr=False)
class NormalizationIndexRow(_Controlled):
    schema: ClassVar[str] = NORMALIZATION_INDEX_ROW_SCHEMA
    discriminator: ClassVar[str] = "normalization_index_row"
    source_allele_id: str
    normalized_identity_id: str | None
    preparation_disposition: str

    def __post_init__(self):
        _require(type(self.source_allele_id) is str
                 and _SOURCE_ALLELE.fullmatch(self.source_allele_id)
                 and self.preparation_disposition in PREPARATION_DISPOSITIONS)
        if self.preparation_disposition == "normalized":
            _require(type(self.normalized_identity_id) is str
                     and _NORMALIZED_ALLELE.fullmatch(self.normalized_identity_id))
        else:
            _require(self.normalized_identity_id is None)


@dataclass(frozen=True, slots=True, repr=False)
class DecisionIndexRow(_Controlled):
    version: ClassVar[int] = 2
    schema: ClassVar[str] = DECISION_INDEX_ROW_SCHEMA
    discriminator: ClassVar[str] = "decision_index_row"
    source_allele_id: str
    normalized_identity_id: str | None
    preparation_disposition: str
    early_filter_disposition: str | None
    scope_disposition: str = field(kw_only=True)

    def __post_init__(self):
        _require(type(self.source_allele_id) is str
                 and _SOURCE_ALLELE.fullmatch(self.source_allele_id)
                 and self.preparation_disposition in PREPARATION_DISPOSITIONS
                 and self.scope_disposition in SCOPE_DISPOSITIONS
                 and (self.scope_disposition == 'not_applicable') ==
                     (self.preparation_disposition in ('non_variant','quarantined')))
        if self.preparation_disposition == "normalized":
            _require(type(self.normalized_identity_id) is str
                     and _NORMALIZED_ALLELE.fullmatch(self.normalized_identity_id)
                     and (self.early_filter_disposition in EARLY_FILTER_DISPOSITIONS
                          if self.scope_disposition == 'in_scope' else self.early_filter_disposition is None))
        else:
            _require(self.normalized_identity_id is None
                     and self.early_filter_disposition is None)


@dataclass(frozen=True, slots=True, repr=False)
class WorkingIdentityIndexRow(_Controlled):
    schema: ClassVar[str] = WORKING_IDENTITY_INDEX_ROW_SCHEMA
    discriminator: ClassVar[str] = "working_identity_index_row"
    identity_id: str
    route: str

    def __post_init__(self):
        _require(self.route in {"sequence", "native"})
        _require((self.route == "sequence" and type(self.identity_id) is str
                  and _NORMALIZED_ALLELE.fullmatch(self.identity_id))
                 or (self.route == "native" and _identifier(self.identity_id)))


@dataclass(frozen=True, slots=True, repr=False)
class ResponseIndexRow(_Controlled):
    schema: ClassVar[str] = RESPONSE_INDEX_ROW_SCHEMA
    discriminator: ClassVar[str] = "response_index_row"
    identity_id: str
    response_ordinal: int
    response_id: str
    response_sha256: str

    def __post_init__(self):
        _require(_identifier(self.identity_id) and _count(self.response_ordinal, 1)
                 and _identifier(self.response_id) and _digest(self.response_sha256))


@dataclass(frozen=True, slots=True, repr=False)
class ContextMembershipIndexRow(_Controlled):
    schema: ClassVar[str] = CONTEXT_INDEX_ROW_SCHEMA
    discriminator: ClassVar[str] = "context_membership_index_row"
    context_id: str
    member_ordinal: int
    identity_id: str

    def __post_init__(self):
        _require(_identifier(self.context_id) and _count(self.member_ordinal, 1)
                 and _identifier(self.identity_id))


@dataclass(frozen=True, slots=True, repr=False)
class JournalEntry(_Controlled):
    schema: ClassVar[str] = JOURNAL_ENTRY_SCHEMA
    discriminator: ClassVar[str] = "journal_entry"
    stage: str
    ordinal: int
    event_type: str
    attempt_id: str | None
    attempt_number: int | None
    process_id: str
    partition_id: str | None
    scientific: bool
    stage_parent_sha256: str
    parents: tuple[ParentIdentity, ...]
    artifacts: tuple[Artifact, ...]
    failure_code: str | None
    previous_entry_sha256: str | None

    def __post_init__(self):
        _require(_name(self.stage) and _count(self.ordinal, 1)
                 and self.ordinal <= MAX_JOURNAL_EVENTS)
        _require(self.event_type in {
            "attempt_started", "partition_completed", "attempt_failed",
            "attempt_interrupted", "stage_completed",
        })
        _require(type(self.process_id) is str
                 and re.fullmatch(r"process-[0-9a-f]{32}", self.process_id)
                 and type(self.scientific) is bool and _digest(self.stage_parent_sha256))
        _parents(self.parents, required={"stage_parent", "partition_plan", "execution_limits",
                                         "executable_stack", "runtime_image"})
        parent_by_role = {value.role: value for value in self.parents}
        _require(parent_by_role["stage_parent"].sha256 == self.stage_parent_sha256)
        _artifacts(self.artifacts)
        _require(self.previous_entry_sha256 is None if self.ordinal == 1
                 else _digest(self.previous_entry_sha256))
        if self.event_type == "stage_completed":
            _require(self.attempt_id is None and self.attempt_number is None
                     and self.partition_id is None and self.scientific is False
                     and len(self.artifacts) == 1 and self.failure_code is None)
        else:
            _require(type(self.attempt_id) is str
                     and re.fullmatch(r"attempt-[0-9a-f]{32}", self.attempt_id)
                     and _count(self.attempt_number, 1) and _name(self.partition_id))
            if self.event_type == "attempt_started":
                _require(not self.artifacts and self.failure_code is None)
            elif self.event_type == "partition_completed":
                _require(len(self.artifacts) == 1 and self.failure_code is None)
            else:
                _require(_name(self.failure_code))


@dataclass(frozen=True, slots=True, repr=False)
class StageRunBinding(_Controlled):
    schema: ClassVar[str] = STAGE_RUN_BINDING_SCHEMA
    discriminator: ClassVar[str] = "stage_run_binding"
    cache_identity: ClassVar[bool] = True
    stage: str
    stage_parent_sha256: str
    executable_stack_sha256: str
    runtime_image_digest: str
    execution_config_sha256: str
    run_options_sha256: str
    input_sha256: str
    reference_sha256: str
    policy_sha256: str
    resource_sha256s: tuple[str, ...]
    parents: tuple[ParentIdentity, ...]

    vcf_evidence_policy_sha256: str = VCF_POLICY_SHA256

    def __post_init__(self):
        _require(self.vcf_evidence_policy_sha256 == VCF_POLICY_SHA256)
        _require(_name(self.stage) and len(self.stage) <= 80)
        _require(all(_digest(value) for value in (
            self.stage_parent_sha256, self.executable_stack_sha256,
            self.execution_config_sha256, self.run_options_sha256,
            self.input_sha256, self.reference_sha256, self.policy_sha256,
        )))
        _require(type(self.runtime_image_digest) is str and
                 re.fullmatch(r"sha256:[0-9a-f]{64}", self.runtime_image_digest))
        _require(type(self.resource_sha256s) is tuple
                 and 0 < len(self.resource_sha256s) <= 64
                 and all(_digest(value) for value in self.resource_sha256s))
        _parents(self.parents, required={
            "stage_parent", "partition_plan", "execution_limits",
            "executable_stack", "runtime_image",
        })
        parents = {item.role: item.sha256 for item in self.parents}
        _require(parents["stage_parent"] == self.stage_parent_sha256
                 and parents["execution_limits"] == self.execution_config_sha256
                 and parents["executable_stack"] == self.executable_stack_sha256
                 and parents["runtime_image"] == self.runtime_image_digest[7:])


@dataclass(frozen=True, slots=True, repr=False)
class StageCommittedPartition:
    """One bounded committed cursor and its publication-backed outputs."""

    completion: JournalEntry
    cursor: tuple[int, ...]
    outputs: tuple[Artifact, ...]
    counts: dict[str, int]


@dataclass(frozen=True, slots=True, repr=False)
class PreparedCaseManifest(_Controlled):
    schema: ClassVar[str] = PREPARED_CASE_SCHEMA
    version: ClassVar[int] = 4
    discriminator: ClassVar[str] = "prepared_case_manifest"
    cache_identity: ClassVar[bool] = True
    status: str
    case_manifest_sha256: str
    original_input_sha256: str
    original_input_size_bytes: int
    selected_sample_sha256: str
    concordance_sha256: str
    header: ParentIdentity
    source_manifest: ParentIdentity
    summaries: tuple[ParentIdentity, ...]
    prefilter_policy_sha256: str
    prefilter_manifest: ParentIdentity
    resources: tuple[ParentIdentity, ...]
    journal_heads: tuple[ParentIdentity, ...]
    working_sets: tuple[ParentIdentity, ...]
    hpo_input_sha256: str
    ontology_manifest_sha256: str
    diagnostic_receipt: ParentIdentity
    execution_limits_sha256: str
    inventory: ParentIdentity
    base_release_sha256: str | None
    executable_stack_sha256: str
    runtime_image_digest: str
    source_accounting: ParentIdentity
    source_audit: ParentIdentity
    preparation_run_options: ParentIdentity
    prepared_index: ParentIdentity
    completed_source_input: ParentIdentity
    completed_source_input_path: str
    source_origin: str
    public_validation_authority_sha256: str | None = None

    def __post_init__(self):
        _require(self.status == "passed" and _count(self.original_input_size_bytes, 1))
        _require(all(_digest(v) for v in (
            self.case_manifest_sha256, self.original_input_sha256, self.selected_sample_sha256,
            self.concordance_sha256, self.prefilter_policy_sha256, self.hpo_input_sha256,
            self.ontology_manifest_sha256, self.execution_limits_sha256,
            self.executable_stack_sha256,
        )))
        _require((self.base_release_sha256 is None) != (self.public_validation_authority_sha256 is None))
        _require(_digest(self.base_release_sha256 or self.public_validation_authority_sha256))
        _require(type(self.runtime_image_digest) is str and
                 re.fullmatch(r"sha256:[0-9a-f]{64}", self.runtime_image_digest))
        _require(type(self.header) is ParentIdentity and self.header.role == "header"
                 and self.header.parent_schema == HEADER_SCHEMA and self.header.parent_version == 1
                 and type(self.source_manifest) is ParentIdentity
                 and self.source_manifest.role == "source_manifest"
                 and self.source_manifest.parent_schema == SOURCE_MANIFEST_SCHEMA
                 and self.source_manifest.parent_version == 1
                 and type(self.prefilter_manifest) is ParentIdentity
                 and self.prefilter_manifest.role == "prefilter_manifest"
                 and self.prefilter_manifest.parent_schema == FILTER_MANIFEST_SCHEMA
                 and self.prefilter_manifest.parent_version == 2
                 and type(self.diagnostic_receipt) is ParentIdentity
                 and self.diagnostic_receipt.role == "diagnostic_receipt"
                 and self.diagnostic_receipt.parent_schema == DIAGNOSTIC_RECEIPT_SCHEMA
                 and self.diagnostic_receipt.parent_version == 4
                 and type(self.inventory) is ParentIdentity and self.inventory.role == "inventory"
                 and self.inventory.parent_schema == INVENTORY_SCHEMA
                 and self.inventory.parent_version == 4)
        _require(self.source_origin in {"produced_this_run", "completed_input"}
                 and type(self.completed_source_input_path) is str
                 and Path(self.completed_source_input_path).is_absolute()
                 and ".." not in Path(self.completed_source_input_path).parts)
        for value, role, schema, version in (
            (self.source_accounting, "source_accounting", "genosapiens.wgs-source-accounting", 1),
            (self.source_audit, "source_audit", "genosapiens.wgs-source-audit", 1),
            (self.preparation_run_options, "preparation_run_options", "genosapiens.wgs-preparation-run-options", 1),
            (self.prepared_index, "prepared_index", "genosapiens.wgs-prepared-index", 3),
            (self.completed_source_input, "completed_source_input", "genosapiens.wgs-completed-source-input", 1),
        ):
            _require(type(value) is ParentIdentity and value.role == role
                     and value.parent_schema == schema and value.parent_version == version)
        for values in (self.summaries, self.resources, self.journal_heads, self.working_sets):
            _parents(values)
        _require({value.role for value in self.summaries} == PREPARED_SUMMARY_ROLES
                 and {value.role for value in self.resources} == PREPARED_RESOURCE_ROLES
                 and {value.role for value in self.working_sets} == {"working_set_manifest"}
                 and self.working_sets[0].parent_schema == WORKING_SET_MANIFEST_SCHEMA
                 and self.working_sets[0].parent_version == 2
                 and bool(self.journal_heads))


_MODELS = (
    LedgerRow, SourceBinding, Artifact, ParentIdentity, CompressionMetadata, Shard,
    LedgerSummary, LedgerManifest, NamedCount, FilterDecision, FilterSummary,
    FilterManifest, WorkingSetRow, WorkingSetSummary, WorkingSetManifest,
    IndexTableBinding, IndexClosure, NormalizationIndexRow, DecisionIndexRow,
    WorkingIdentityIndexRow, ResponseIndexRow, ContextMembershipIndexRow,
    JournalEntry, StageRunBinding, PreparedCaseManifest,
)
_CHILDREN = {
    Shard: {"artifact": Artifact, "compression": CompressionMetadata},
    LedgerSummary: {"reports": (Artifact,)},
    LedgerManifest: {"source": SourceBinding, "header": Artifact,
                     "parents": (ParentIdentity,), "summary": LedgerSummary, "shards": (Shard,)},
    FilterSummary: {"preparation_counts": (NamedCount,), "early_filter_counts": (NamedCount,),
                    "reason_counts": (NamedCount,), "scope_counts": (NamedCount,)},
    FilterManifest: {"parents": (ParentIdentity,), "summary": FilterSummary, "shards": (Shard,)},
    WorkingSetManifest: {"parents": (ParentIdentity,), "summary": WorkingSetSummary,
                         "shards": (Shard,)},
    WorkingSetSummary: {"analysis_outcome_counts": (NamedCount,)},
    IndexClosure: {
        "normalization": IndexTableBinding, "decisions": IndexTableBinding,
        "working_identities": IndexTableBinding, "responses": IndexTableBinding,
        "context_memberships": IndexTableBinding,
    },
    JournalEntry: {"parents": (ParentIdentity,), "artifacts": (Artifact,)},
    StageRunBinding: {"parents": (ParentIdentity,)},
    PreparedCaseManifest: {
        "header": ParentIdentity, "source_manifest": ParentIdentity,
        "summaries": (ParentIdentity,), "prefilter_manifest": ParentIdentity,
        "resources": (ParentIdentity,), "journal_heads": (ParentIdentity,),
        "working_sets": (ParentIdentity,), "diagnostic_receipt": ParentIdentity,
        "inventory": ParentIdentity, "source_accounting": ParentIdentity,
        "source_audit": ParentIdentity, "preparation_run_options": ParentIdentity,
        "prepared_index": ParentIdentity, "completed_source_input": ParentIdentity,
    },
}


def _decode(data, cls):
    _require(type(data) is dict and set(data) == {f.name for f in fields(cls)})
    values = dict(data)
    for name, child in _CHILDREN.get(cls, {}).items():
        if type(child) is tuple:
            _require(type(values[name]) is list)
            values[name] = tuple(_decode(v, child[0]) for v in values[name])
        else:
            values[name] = _decode(values[name], child)
    tuple_fields = {
        LedgerRow: ("source_allele_ids",), FilterDecision: ("reason_codes",),
        WorkingSetRow: ("source_allele_ids",),
        StageRunBinding: ("resource_sha256s",),
    }
    for name in tuple_fields.get(cls, ()):
        _require(type(values[name]) is list)
        values[name] = tuple(values[name])
    return cls(**values)


def load_model(content: bytes, cls):
    try:
        _require(cls in _MODELS)
        value = _json(content, MAX_ROW_BYTES if cls is LedgerRow else MAX_MANIFEST_BYTES)
        _require(type(value) is dict and set(value) ==
                 {"schema", "version", "kind", "classification", "llm_policy", "data"})
        _require(value["schema"] == cls.schema and type(value["version"]) is int
                 and value["version"] == cls.version)
        _require(value["kind"] == cls.discriminator and value["classification"] == "controlled"
                 and value["llm_policy"] == "no-llm")
        return _decode(value["data"], cls)
    except (ValueError, TypeError, KeyError, AttributeError, RecursionError, OverflowError):
        raise LedgerError("invalid WGS ledger envelope") from None


@contextmanager
def _directory(path):
    """Walk with descriptors: no ancestor symlinks or path lookup races."""
    fd = -1
    try:
        path = Path(path)
        _require(path.is_absolute() and ".." not in path.parts and len(path.parts) > 1)
        flags = os.O_RDONLY | os.O_DIRECTORY | os.O_NOFOLLOW
        fd = os.open(path.anchor, flags)
        for part in path.parts[1:]:
            child = os.open(part, flags, dir_fd=fd)
            os.close(fd)
            fd = child
        metadata = os.fstat(fd)
        _require(stat.S_IMODE(metadata.st_mode) == 0o700 and metadata.st_uid == os.getuid())
    except (OSError, ValueError):
        if fd >= 0: os.close(fd)
        raise LedgerError("WGS ledger directory is unavailable or unsafe") from None
    try:
        yield fd
    finally:
        if fd >= 0:
            os.close(fd)


@contextmanager
def _read(fd, name):
    _require(_name(name))
    descriptor = -1
    try:
        descriptor = os.open(name, os.O_RDONLY | os.O_NOFOLLOW | os.O_NONBLOCK, dir_fd=fd)
        metadata = os.fstat(descriptor)
        _require(stat.S_ISREG(metadata.st_mode) and stat.S_IMODE(metadata.st_mode) == 0o600
                 and metadata.st_uid == os.getuid() and metadata.st_nlink == 1)
        with os.fdopen(descriptor, "rb") as stream:
            descriptor = -1
            yield stream, metadata
            final = os.fstat(stream.fileno())
            _require((metadata.st_ino, metadata.st_size, metadata.st_mtime_ns, metadata.st_ctime_ns)
                     == (final.st_ino, final.st_size, final.st_mtime_ns, final.st_ctime_ns))
    except OSError:
        raise LedgerError("WGS ledger artifact is unavailable or unsafe") from None
    finally:
        if descriptor >= 0:
            os.close(descriptor)


@contextmanager
def _atomic(fd, name):
    _require(_name(name))
    temporary = ".pending-" + uuid.uuid4().hex
    destination_linked = False
    try:
        descriptor = os.open(temporary, os.O_WRONLY | os.O_CREAT | os.O_EXCL | os.O_NOFOLLOW,
                             0o600, dir_fd=fd)
        with os.fdopen(descriptor, "wb") as stream:
            os.fchmod(stream.fileno(), 0o600)
            yield stream
            stream.flush()
            os.fsync(stream.fileno())
        # link is atomic and refuses an existing destination, including symlinks.
        os.link(temporary, name, src_dir_fd=fd, dst_dir_fd=fd, follow_symlinks=False)
        destination_linked = True
        os.unlink(temporary, dir_fd=fd)
        os.fsync(fd)
    except OSError:
        if destination_linked:
            try:
                os.unlink(name, dir_fd=fd)
                os.fsync(fd)
            except OSError:
                pass
        raise LedgerError("cannot publish immutable WGS ledger artifact") from None
    finally:
        try:
            os.unlink(temporary, dir_fd=fd)
        except FileNotFoundError:
            pass


def _authenticate(fd, artifact):
    _require(type(artifact) is Artifact)
    digest = hashlib.sha256()
    with _read(fd, artifact.name) as (stream, metadata):
        _require(metadata.st_size == artifact.size_bytes)
        while block := stream.read(1024 * 1024):
            digest.update(block)
    _require(digest.hexdigest() == artifact.sha256)


class _HashingWriter:
    """Observe the exact compressed bytes below gzip while they are written."""

    def __init__(self, stream):
        self.stream = stream
        self.digest = hashlib.sha256()
        self.size = 0

    def write(self, content):
        written = self.stream.write(content)
        _require(written == len(content))
        self.digest.update(content)
        self.size += written
        _require(self.size <= MAX_SHARD_BYTES + 65536)
        return written

    def flush(self):
        self.stream.flush()


class _HashingReader:
    def __init__(self, stream, maximum_bytes):
        self.stream = stream
        self.maximum_bytes = maximum_bytes
        self.digest = hashlib.sha256()
        self.size = 0

    def read(self, amount=-1):
        content = self.stream.read(amount)
        self.size += len(content)
        _require(self.size <= self.maximum_bytes)
        self.digest.update(content)
        return content


def _remember_published_shard(fd, shard):
    metadata = os.stat(shard.artifact.name, dir_fd=fd, follow_symlinks=False)
    _require(stat.S_ISREG(metadata.st_mode) and stat.S_IMODE(metadata.st_mode) == 0o600
             and metadata.st_uid == os.getuid() and metadata.st_nlink == 1
             and metadata.st_size == shard.artifact.size_bytes)
    directory = os.fstat(fd)
    proof = (
        directory.st_dev, directory.st_ino, metadata.st_dev, metadata.st_ino,
        metadata.st_size, metadata.st_mtime_ns, metadata.st_ctime_ns,
    )
    object.__setattr__(shard, "_publication", proof)
    object.__setattr__(shard.artifact, "_publication", proof)


def _remember_published_artifact(fd, artifact):
    metadata = os.stat(artifact.name, dir_fd=fd, follow_symlinks=False)
    _require(stat.S_ISREG(metadata.st_mode) and stat.S_IMODE(metadata.st_mode) == 0o600
             and metadata.st_uid == os.getuid() and metadata.st_nlink == 1
             and metadata.st_size == artifact.size_bytes)
    directory = os.fstat(fd)
    object.__setattr__(artifact, "_publication", (
        directory.st_dev, directory.st_ino, metadata.st_dev, metadata.st_ino,
        metadata.st_size, metadata.st_mtime_ns, metadata.st_ctime_ns,
    ))


def _require_published_shard(fd, shard):
    _require(type(shard) is Shard)
    proof = getattr(shard, "_publication", None)
    _require(type(proof) is tuple and len(proof) == 7)
    metadata = os.stat(shard.artifact.name, dir_fd=fd, follow_symlinks=False)
    directory = os.fstat(fd)
    _require(stat.S_ISREG(metadata.st_mode) and stat.S_IMODE(metadata.st_mode) == 0o600
             and metadata.st_uid == os.getuid() and metadata.st_nlink == 1
             and metadata.st_size == shard.artifact.size_bytes
             and proof == (directory.st_dev, directory.st_ino, metadata.st_dev,
                           metadata.st_ino, metadata.st_size, metadata.st_mtime_ns,
                           metadata.st_ctime_ns))


def _require_published_artifact(fd, artifact):
    _require(type(artifact) is Artifact)
    proof = getattr(artifact, "_publication", None)
    _require(type(proof) is tuple and len(proof) == 7)
    metadata = os.stat(artifact.name, dir_fd=fd, follow_symlinks=False)
    directory = os.fstat(fd)
    _require(stat.S_ISREG(metadata.st_mode) and stat.S_IMODE(metadata.st_mode) == 0o600
             and metadata.st_uid == os.getuid() and metadata.st_nlink == 1
             and metadata.st_size == artifact.size_bytes
             and proof == (directory.st_dev, directory.st_ino, metadata.st_dev,
                           metadata.st_ino, metadata.st_size, metadata.st_mtime_ns,
                           metadata.st_ctime_ns))


def require_published_shards(root, shards):
    """Check writer-supplied descriptors without rereading shard contents."""
    _require(type(shards) is tuple and len(shards) <= MAX_SHARDS)
    with _directory(root) as fd:
        for shard in shards:
            _require_published_shard(fd, shard)


def attest_existing_artifact(root, artifact):
    """Adopt a small immutable artifact after checking its canonical bytes."""
    _require(type(artifact) is Artifact)
    read_artifact(root, artifact)
    with _directory(root) as fd:
        _remember_published_artifact(fd, artifact)
    return artifact


def publish_bytes(root, name, content):
    _require(type(content) is bytes and 0 < len(content) <= MAX_MANIFEST_BYTES)
    artifact = Artifact(name, _sha(content), len(content))
    with _directory(root) as fd, _atomic(fd, name) as stream:
        stream.write(content)
    with _directory(root) as fd:
        _remember_published_artifact(fd, artifact)
    return artifact


def publish_stream(root, name, chunks, *, maximum_bytes=MAX_SHARD_BYTES):
    """Atomically publish a bounded artifact without materializing all bytes."""

    _require(
        _name(name)
        and _count(maximum_bytes, 1)
        and maximum_bytes <= MAX_SHARD_BYTES
    )
    digest = hashlib.sha256()
    size = 0
    with _directory(root) as fd, _atomic(fd, name) as stream:
        for chunk in chunks:
            _require(type(chunk) is bytes and bool(chunk))
            size += len(chunk)
            _require(size <= maximum_bytes)
            digest.update(chunk)
            stream.write(chunk)
        _require(size > 0)
    artifact = Artifact(name, digest.hexdigest(), size)
    with _directory(root) as fd:
        _remember_published_artifact(fd, artifact)
    return artifact


def verify_artifact(root, artifact, *, maximum_bytes=MAX_SHARD_BYTES):
    """Authenticate a bounded artifact by size, file safety and complete digest."""

    _require(
        type(artifact) is Artifact
        and _count(maximum_bytes, 1)
        and maximum_bytes <= MAX_SHARD_BYTES
        and artifact.size_bytes <= maximum_bytes
    )
    with _directory(root) as fd:
        _authenticate(fd, artifact)


def iter_artifact_lines(
    root, artifact, *, maximum_bytes=MAX_SHARD_BYTES, maximum_line_bytes=MAX_ROW_BYTES,
):
    """Authenticate one bounded read, then yield from a private temporary spool."""

    _require(
        _count(maximum_line_bytes, 1)
        and maximum_line_bytes <= MAX_ROW_BYTES
    )
    _require(type(artifact) is Artifact and artifact.size_bytes <= maximum_bytes)
    with tempfile.TemporaryFile(mode="w+b") as spool:
        with _directory(root) as fd, _read(fd, artifact.name) as (stream, metadata):
            _require(metadata.st_size == artifact.size_bytes)
            total = 0
            digest = hashlib.sha256()
            while content := stream.readline(maximum_line_bytes + 1):
                total += len(content)
                _require(len(content) <= maximum_line_bytes
                         and content.endswith(b"\n")
                         and total <= artifact.size_bytes)
                digest.update(content)
                _require(spool.write(content) == len(content))
            _require(total == artifact.size_bytes and digest.hexdigest() == artifact.sha256)
        spool.seek(0)
        while content := spool.readline(maximum_line_bytes + 1):
            yield content


def read_artifact(root, artifact, *, maximum_bytes=MAX_MANIFEST_BYTES):
    _require(_count(maximum_bytes, 1) and maximum_bytes <= MAX_MANIFEST_BYTES)
    _require(type(artifact) is Artifact and artifact.size_bytes <= maximum_bytes)
    with _directory(root) as fd, _read(fd, artifact.name) as (stream, metadata):
        _require(metadata.st_size == artifact.size_bytes)
        content = stream.read(maximum_bytes + 1)
    _require(len(content) == artifact.size_bytes and _sha(content) == artifact.sha256)
    return content


def publish_shard(root, rows: Iterable[EncodedLedgerRow], *, ordinal, first_record, previous=None):
    """Publish one bounded partition. The caller must rotate at declared limits."""
    _require(_count(ordinal, 1) and ordinal <= MAX_SHARDS and _count(first_record, 1))
    _require(previous is None if ordinal == 1 else type(previous) is Shard)
    if previous is not None:
        _require(previous.ordinal + 1 == ordinal and previous.last_record + 1 == first_record)
    else:
        _require(first_record == 1)
    name = f"shard-{ordinal:06d}.jsonl.gz"
    digest = hashlib.sha256()
    size, count, physical, biological = 0, 0, 0, 0
    with _directory(root) as fd:
        if previous is not None:
            _require_published_shard(fd, previous)
        with _atomic(fd, name) as stream:
            writer = _HashingWriter(stream)
            with gzip.GzipFile(filename="", mode="wb", fileobj=writer, mtime=0, compresslevel=6) as compressed:
                for encoded in rows:
                    _require(type(encoded) is EncodedLedgerRow
                             and type(encoded.row) is LedgerRow
                             and encoded.record_number == first_record + count
                             and encoded.record_number == encoded.row.record_number
                             and encoded.physical_alternate_slots == len(encoded.row.source_allele_ids)
                             and encoded.biological_alternates == encoded.row.biological_alternate_count
                             and type(encoded.content) is bytes
                             and 0 < len(encoded.content) <= MAX_ROW_BYTES)
                    content = encoded.content
                    size += len(content)
                    count += 1
                    _require(count <= MAX_SHARD_ROWS and size <= MAX_SHARD_BYTES)
                    digest.update(content)
                    compressed.write(content)
                    physical += encoded.physical_alternate_slots
                    biological += encoded.biological_alternates
                _require(count > 0)
        shard = Shard(
            role="source", row_schema=SOURCE_ROW_SCHEMA, row_version=1, ordinal=ordinal,
            artifact=Artifact(name, writer.digest.hexdigest(), writer.size), compression=GZIP_JSONL,
            uncompressed_sha256=digest.hexdigest(), uncompressed_bytes=size,
            first_record=first_record, last_record=first_record + count - 1, rows=count,
            physical_alternate_slots=physical, biological_alternates=biological,
            previous_shard_sha256=None if previous is None else previous.sha256,
        )
        _remember_published_shard(fd, shard)
        _remember_published_artifact(fd, shard.artifact)
        return shard


def _gzip_lines(stream, maximum_bytes):
    """Yield one exact gzip member without allowing unbounded inflation."""
    _require(_count(maximum_bytes, 1) and maximum_bytes <= MAX_SHARD_BYTES)
    block = stream.read(10)
    _require(block == b"\x1f\x8b\x08\x00\x00\x00\x00\x00\x00\xff")
    inflater = zlib.decompressobj(16 + zlib.MAX_WBITS)
    pending = bytearray()
    total = 0
    ended = False
    while block:
        _require(not ended)
        compressed = block
        while True:
            output = inflater.decompress(compressed, 64 * 1024)
            compressed = inflater.unconsumed_tail
            if output:
                total += len(output)
                _require(total <= maximum_bytes)
                pending.extend(output)
                while True:
                    newline = pending.find(b"\n")
                    if newline < 0:
                        break
                    _require(newline + 1 <= MAX_ROW_BYTES)
                    content = bytes(pending[:newline + 1])
                    del pending[:newline + 1]
                    yield content
                _require(len(pending) <= MAX_ROW_BYTES)
            if inflater.eof:
                _require(not inflater.unused_data and not compressed)
                ended = True
                break
            if compressed:
                continue
            if len(output) == 64 * 1024:
                compressed = b""
                continue
            break
        if ended:
            _require(stream.read(1) == b"")
            break
        block = stream.read(64 * 1024)
    _require(ended and not pending)


def _scan_shard(root, shard, *, on_content=None):
    _require(type(shard) is Shard)
    with _directory(root) as fd, _read(fd, shard.artifact.name) as (stream, metadata):
        _require(metadata.st_size == shard.artifact.size_bytes)
        compressed = _HashingReader(stream, shard.artifact.size_bytes)
        digest = hashlib.sha256()
        count, size, physical, biological = 0, 0, 0, 0
        for content in _gzip_lines(compressed, shard.uncompressed_bytes):
            size += len(content)
            _require(count < shard.rows)
            row = load_model(content, LedgerRow)
            _require(row.record_number == shard.first_record + count)
            digest.update(content)
            count += 1
            physical += len(row.source_allele_ids)
            biological += row.biological_alternate_count
            if on_content is not None:
                _require(on_content(content) == len(content))
        _require(count == shard.rows
                 and physical == shard.physical_alternate_slots
                 and biological == shard.biological_alternates
                 and size == shard.uncompressed_bytes
                 and digest.hexdigest() == shard.uncompressed_sha256
                 and compressed.size == shard.artifact.size_bytes
                 and compressed.digest.hexdigest() == shard.artifact.sha256)


def verify_shard(root, shard):
    """Authenticate an entire bounded shard before any dependent output is accepted."""
    try:
        _scan_shard(root, shard)
    except (OSError, EOFError, zlib.error):
        raise LedgerError("cannot authenticate WGS ledger shard") from None


def iter_shard(root, shard) -> Iterator[LedgerRow]:
    """Close one bounded shard before exposing rows from its temporary spool."""
    try:
        with tempfile.TemporaryFile(mode="w+b") as spool:
            _scan_shard(root, shard, on_content=spool.write)
            spool.seek(0)
            for _ordinal in range(shard.rows):
                content = spool.readline(MAX_ROW_BYTES + 1)
                _require(content and len(content) <= MAX_ROW_BYTES)
                yield load_model(content, LedgerRow)
            _require(spool.read(1) == b"")
    except (OSError, EOFError, zlib.error):
        raise LedgerError("cannot authenticate WGS ledger shard") from None


def _stage_row_counts(row, role):
    if role == "filter":
        _require(type(row) is FilterDecision)
        return 1, int(row.preparation_disposition == "normalized")
    _require(role == "working_set" and type(row) is WorkingSetRow)
    return len(row.source_allele_ids), len(row.source_allele_ids)


def _scan_stage_shard(root, shard, *, cls, role, on_content=None):
    _require(type(shard) is Shard and shard.role == role and shard.row_schema == cls.schema)
    with _directory(root) as fd, _read(fd, shard.artifact.name) as (stream, metadata):
        _require(metadata.st_size == shard.artifact.size_bytes)
        compressed = _HashingReader(stream, shard.artifact.size_bytes)
        digest = hashlib.sha256()
        count, size, physical, biological = 0, 0, 0, 0
        for content in _gzip_lines(compressed, shard.uncompressed_bytes):
            size += len(content)
            _require(count < shard.rows)
            row = load_model(content, cls)
            if role == "working_set":
                _require(row.ordinal == shard.first_record + count)
            digest.update(content)
            count += 1
            row_physical, row_biological = _stage_row_counts(row, role)
            physical += row_physical
            biological += row_biological
            if on_content is not None:
                _require(on_content(content) == len(content))
        _require(count == shard.rows
                 and physical == shard.physical_alternate_slots
                 and biological == shard.biological_alternates
                 and size == shard.uncompressed_bytes
                 and digest.hexdigest() == shard.uncompressed_sha256
                 and compressed.size == shard.artifact.size_bytes
                 and compressed.digest.hexdigest() == shard.artifact.sha256)


def _publish_stage_shard(
    root, rows, *, cls, role, ordinal, first_record, previous=None,
):
    _require(role in {"filter", "working_set"}
             and _count(ordinal, 1) and ordinal <= MAX_SHARDS
             and _count(first_record, 1))
    _require(previous is None if ordinal == 1 else type(previous) is Shard)
    if previous is not None:
        _require(previous.role == role and previous.ordinal + 1 == ordinal
                 and previous.last_record + 1 == first_record)
    else:
        _require(first_record == 1)
    prefix = "decisions" if role == "filter" else "retained"
    name = f"{prefix}-{ordinal:06d}.jsonl.gz"
    digest = hashlib.sha256()
    size, count, physical, biological = 0, 0, 0, 0
    with _directory(root) as fd:
        if previous is not None:
            _require_published_shard(fd, previous)
        with _atomic(fd, name) as stream:
            writer = _HashingWriter(stream)
            with gzip.GzipFile(
                filename="", mode="wb", fileobj=writer, mtime=0, compresslevel=6,
            ) as compressed:
                for row in rows:
                    _require(type(row) is cls)
                    if role == "working_set":
                        _require(row.ordinal == first_record + count)
                    content = row.as_bytes()
                    count += 1
                    size += len(content)
                    _require(count <= MAX_SHARD_ROWS and size <= MAX_SHARD_BYTES)
                    digest.update(content)
                    compressed.write(content)
                    row_physical, row_biological = _stage_row_counts(row, role)
                    physical += row_physical
                    biological += row_biological
                _require(count > 0)
        shard = Shard(
            role=role, row_schema=cls.schema, row_version=cls.version, ordinal=ordinal,
            artifact=Artifact(name, writer.digest.hexdigest(), writer.size), compression=GZIP_JSONL,
            uncompressed_sha256=digest.hexdigest(), uncompressed_bytes=size,
            first_record=first_record, last_record=first_record + count - 1, rows=count,
            physical_alternate_slots=physical, biological_alternates=biological,
            previous_shard_sha256=None if previous is None else previous.sha256,
        )
        _remember_published_shard(fd, shard)
        return shard


def publish_filter_shard(root, rows, *, ordinal, first_record, previous=None):
    return _publish_stage_shard(
        root, rows, cls=FilterDecision, role="filter", ordinal=ordinal,
        first_record=first_record, previous=previous,
    )


def publish_working_set_shard(root, rows, *, ordinal, first_record, previous=None):
    return _publish_stage_shard(
        root, rows, cls=WorkingSetRow, role="working_set", ordinal=ordinal,
        first_record=first_record, previous=previous,
    )


def iter_filter_shard(root, shard):
    yield from _iter_stage_shard(root, shard, cls=FilterDecision, role="filter")


def iter_working_set_shard(root, shard):
    yield from _iter_stage_shard(root, shard, cls=WorkingSetRow, role="working_set")


def _iter_stage_shard(root, shard, *, cls, role):
    try:
        with tempfile.TemporaryFile(mode="w+b") as spool:
            _scan_stage_shard(root, shard, cls=cls, role=role, on_content=spool.write)
            spool.seek(0)
            for _ordinal in range(shard.rows):
                content = spool.readline(MAX_ROW_BYTES + 1)
                _require(content and len(content) <= MAX_ROW_BYTES)
                yield load_model(content, cls)
            _require(spool.read(1) == b"")
    except (OSError, EOFError, zlib.error):
        raise LedgerError("cannot authenticate WGS stage shard") from None


class _StageMutationLock:
    __slots__ = ("stage", "directory_device", "directory_inode", "lock_inode", "fd", "active")

    def __init__(self, stage, directory_device, directory_inode, lock_inode, fd):
        self.stage = stage
        self.directory_device = directory_device
        self.directory_inode = directory_inode
        self.lock_inode = lock_inode
        self.fd = fd
        self.active = True

    def __repr__(self):
        return "<active private WGS stage lock>" if self.active else "<released WGS stage lock>"


@contextmanager
def stage_lock(root, stage):
    """Hold one fail-fast advisory lock for all mutations of a stage root."""
    _require(_name(stage) and len(stage) <= 80)
    name = f"{stage}.lock"
    with _directory(root) as directory_fd:
        descriptor = -1
        locked = False
        lock = None
        try:
            descriptor = os.open(name, os.O_RDWR | os.O_CREAT | os.O_NOFOLLOW, 0o600,
                                 dir_fd=directory_fd)
            metadata = os.fstat(descriptor)
            _require(stat.S_ISREG(metadata.st_mode) and metadata.st_uid == os.getuid()
                     and stat.S_IMODE(metadata.st_mode) == 0o600 and metadata.st_nlink == 1)
            try:
                fcntl.flock(descriptor, fcntl.LOCK_EX | fcntl.LOCK_NB)
            except BlockingIOError:
                raise LedgerError("WGS stage is already locked") from None
            locked = True
            directory = os.fstat(directory_fd)
            token = canonical({
                "schema": "genosapiens.wgs-stage-lock", "version": 1,
                "stage": stage, "process_id": "process-" + uuid.uuid4().hex,
            })
            os.ftruncate(descriptor, 0)
            os.lseek(descriptor, 0, os.SEEK_SET)
            view = memoryview(token)
            while view:
                written = os.write(descriptor, view)
                _require(written > 0)
                view = view[written:]
            os.fsync(descriptor)
            os.fsync(directory_fd)
            lock = _StageMutationLock(
                stage, directory.st_dev, directory.st_ino, metadata.st_ino, descriptor,
            )
        except (OSError, LedgerError) as error:
            if descriptor >= 0:
                os.close(descriptor)
            if isinstance(error, LedgerError):
                raise
            raise LedgerError("cannot acquire private WGS stage lock") from None
        try:
            yield lock
        finally:
            if lock is not None:
                lock.active = False
            if locked:
                try:
                    fcntl.flock(descriptor, fcntl.LOCK_UN)
                except OSError:
                    pass
            if descriptor >= 0:
                os.close(descriptor)


def _require_lock(root, stage, lock):
    _require(type(lock) is _StageMutationLock and lock.active and lock.stage == stage
             and lock.fd >= 0)
    try:
        with _directory(root) as directory_fd:
            directory = os.fstat(directory_fd)
            held = os.fstat(lock.fd)
            named = os.stat(f"{stage}.lock", dir_fd=directory_fd, follow_symlinks=False)
            _require((directory.st_dev, directory.st_ino)
                     == (lock.directory_device, lock.directory_inode)
                     and stat.S_ISREG(held.st_mode) and stat.S_ISREG(named.st_mode)
                     and held.st_ino == named.st_ino == lock.lock_inode
                     and held.st_nlink == named.st_nlink == 1
                     and held.st_uid == named.st_uid == os.getuid()
                     and stat.S_IMODE(held.st_mode) == stat.S_IMODE(named.st_mode) == 0o600)
    except OSError:
        raise LedgerError("invalid private WGS stage lock") from None


_STAGE_STATE_DDL = """
PRAGMA application_id=1196972365;
PRAGMA user_version=1;
CREATE TABLE run_binding (
    singleton INTEGER PRIMARY KEY CHECK(singleton=1),
    canonical_bytes BLOB NOT NULL
);
CREATE TABLE events (
    ordinal INTEGER PRIMARY KEY CHECK(ordinal > 0),
    event_bytes BLOB NOT NULL,
    event_sha256 TEXT NOT NULL,
    event_type TEXT NOT NULL,
    partition_id TEXT,
    attempt_id TEXT
);
CREATE TABLE commits (
    ordinal INTEGER PRIMARY KEY REFERENCES events(ordinal),
    partition_id TEXT NOT NULL UNIQUE,
    cursor_bytes BLOB NOT NULL,
    outputs_bytes BLOB NOT NULL,
    counts_bytes BLOB NOT NULL
);
CREATE TABLE committed_outputs (
    name TEXT PRIMARY KEY,
    commit_ordinal INTEGER NOT NULL REFERENCES commits(ordinal),
    metadata_bytes BLOB NOT NULL
);
"""
_STAGE_STATE_TABLES = {"run_binding", "events", "commits", "committed_outputs"}


def _replace_private_pointer(root, name, content):
    _require(_name(name) and type(content) is bytes and len(content) <= 4096)
    pending = ".pending-" + uuid.uuid4().hex
    with _directory(root) as fd:
        try:
            try:
                metadata = os.stat(name, dir_fd=fd, follow_symlinks=False)
            except FileNotFoundError:
                metadata = None
            if metadata is not None:
                _require(stat.S_ISREG(metadata.st_mode)
                         and stat.S_IMODE(metadata.st_mode) == 0o600
                         and metadata.st_uid == os.getuid() and metadata.st_nlink == 1)
            descriptor = os.open(
                pending, os.O_WRONLY | os.O_CREAT | os.O_EXCL | os.O_NOFOLLOW,
                0o600, dir_fd=fd,
            )
            with os.fdopen(descriptor, "wb") as stream:
                stream.write(content)
                stream.flush()
                os.fsync(stream.fileno())
            os.replace(pending, name, src_dir_fd=fd, dst_dir_fd=fd)
            os.fsync(fd)
        except OSError:
            raise LedgerError("cannot publish WGS stage head") from None
        finally:
            try:
                os.unlink(pending, dir_fd=fd)
            except FileNotFoundError:
                pass


class StageCommitState:
    """One locked stage database coordinates journal events and output cursors."""

    def __init__(self, stage_root, state_path, binding, *, lock, create=False):
        _require(type(binding) is StageRunBinding and type(create) is bool)
        self.stage_root = Path(stage_root)
        self.state_path = Path(state_path)
        self.binding = binding
        self.lock = lock
        self.connection = None
        _require_lock(self.stage_root, binding.stage, lock)
        _require(self.state_path.is_absolute() and _name(self.state_path.name)
                 and self.state_path.suffix == ".sqlite")
        with _directory(self.state_path.parent):
            pass
        exists = self.state_path.exists() or self.state_path.is_symlink()
        _require((create and not exists) or (not create and exists))
        try:
            if create:
                descriptor = os.open(
                    self.state_path, os.O_WRONLY | os.O_CREAT | os.O_EXCL | os.O_NOFOLLOW,
                    0o600,
                )
                os.close(descriptor)
            metadata = self.state_path.lstat()
            _require(not self.state_path.is_symlink()
                     and stat.S_ISREG(metadata.st_mode)
                     and stat.S_IMODE(metadata.st_mode) == 0o600
                     and metadata.st_uid == os.getuid() and metadata.st_nlink == 1)
            connection = sqlite3.connect(self.state_path, isolation_level=None)
            self.connection = connection
            connection.enable_load_extension(False)
            connection.execute("PRAGMA journal_mode=DELETE")
            connection.execute("PRAGMA synchronous=FULL")
            connection.execute("PRAGMA foreign_keys=ON")
            connection.execute(f"PRAGMA cache_size=-{SQLITE_CACHE_KIB}")
            connection.execute("PRAGMA temp_store=FILE")
            connection.execute("PRAGMA mmap_size=0")
            connection.execute("PRAGMA trusted_schema=OFF")
            if create:
                connection.executescript(_STAGE_STATE_DDL)
                connection.execute(f"PRAGMA user_version={2 if binding.stage == 'prefilter' else 1}")
                connection.execute(
                    "INSERT INTO run_binding VALUES (1,?)", (binding.as_bytes(),),
                )
            else:
                self._check_schema()
                _require(connection.execute(
                    "SELECT canonical_bytes FROM run_binding WHERE singleton=1"
                ).fetchall() == [(binding.as_bytes(),)])
            final = self.state_path.lstat()
            _require((metadata.st_dev, metadata.st_ino) == (final.st_dev, final.st_ino))
            self._recover_head()
        except (OSError, sqlite3.Error, LedgerError):
            self.close()
            raise LedgerError("WGS stage checkpoint is inconsistent") from None

    def __enter__(self):
        return self

    def __exit__(self, _kind, _value, _traceback):
        self.close()

    def close(self):
        if self.connection is not None:
            self.connection.close()
            self.connection = None

    def _check_schema(self):
        connection = self.connection
        _require(connection.execute("PRAGMA application_id").fetchone() == (1196972365,)
                 and connection.execute("PRAGMA user_version").fetchone() == (2 if self.binding.stage == 'prefilter' else 1,)
                 and connection.execute("PRAGMA journal_mode").fetchone() == ("delete",)
                 and connection.execute("PRAGMA synchronous").fetchone() == (2,)
                 and connection.execute("PRAGMA foreign_keys").fetchone() == (1,))
        tables = {
            row[0] for row in connection.execute(
                "SELECT name FROM sqlite_master WHERE type='table' AND name NOT LIKE 'sqlite_%'"
            )
        }
        _require(_STAGE_STATE_TABLES <= tables)

    def _latest(self):
        return self.connection.execute(
            "SELECT ordinal,event_bytes,event_sha256 FROM events ORDER BY ordinal DESC LIMIT 1"
        ).fetchone()

    def _head_pointer(self):
        name = f"{self.binding.stage}-head.json"
        with _directory(self.stage_root) as fd:
            try:
                os.stat(name, dir_fd=fd, follow_symlinks=False)
            except FileNotFoundError:
                return None
            try:
                with _read(fd, name) as (stream, metadata):
                    _require(metadata.st_size <= 4096)
                    content = stream.read(4097)
            except LedgerError:
                raise
        value = _json(content, 4096)
        _require(type(value) is dict and set(value) == {
            "schema", "version", "stage", "ordinal", "sha256",
        } and value["schema"] == "genosapiens.wgs-stage-journal-head"
                 and value["version"] == 1 and value["stage"] == self.binding.stage
                 and _count(value["ordinal"], 1) and _digest(value["sha256"]))
        return value["ordinal"], value["sha256"]

    def _publish_event_and_head(self, ordinal, content, digest):
        entry = load_model(content, JournalEntry)
        _require(entry.ordinal == ordinal and entry.sha256 == digest
                 and entry.stage == self.binding.stage
                 and entry.stage_parent_sha256 == self.binding.stage_parent_sha256
                 and entry.parents == self.binding.parents)
        name = f"{entry.stage}-{ordinal:06d}.journal.json"
        expected = Artifact(name, _sha(content), len(content))
        path = self.stage_root / name
        if path.exists() or path.is_symlink():
            _require(read_artifact(self.stage_root, expected) == content)
        else:
            publish_bytes(self.stage_root, name, content)
        _replace_private_pointer(self.stage_root, f"{entry.stage}-head.json", canonical({
            "schema": "genosapiens.wgs-stage-journal-head", "version": 1,
            "stage": entry.stage, "ordinal": ordinal, "sha256": digest,
        }))

    def _recover_head(self):
        latest = self._latest()
        pointer = self._head_pointer()
        if latest is None:
            _require(pointer is None)
            _require(not (self.stage_root / f"{self.binding.stage}-000001.journal.json").exists())
            return
        ordinal, content, digest = latest
        _require(_count(ordinal, 1) and ordinal <= MAX_JOURNAL_EVENTS
                 and type(content) is bytes and len(content) <= MAX_MANIFEST_BYTES
                 and _digest(digest) and _sha(content) == digest)
        entry = load_model(content, JournalEntry)
        _require(entry.ordinal == ordinal and entry.sha256 == digest
                 and entry.stage == self.binding.stage
                 and entry.stage_parent_sha256 == self.binding.stage_parent_sha256
                 and entry.parents == self.binding.parents)
        _require(not (self.stage_root / f"{self.binding.stage}-{ordinal + 1:06d}.journal.json").exists())
        if pointer == (ordinal, digest):
            expected = Artifact(
                f"{self.binding.stage}-{ordinal:06d}.journal.json",
                digest, len(content),
            )
            _require(read_artifact(self.stage_root, expected) == content)
        else:
            if ordinal == 1:
                _require(pointer is None)
            else:
                previous = self.connection.execute(
                    "SELECT event_sha256 FROM events WHERE ordinal=?", (ordinal - 1,)
                ).fetchone()
                _require(previous is not None and pointer == (ordinal - 1, previous[0]))
            self._publish_event_and_head(ordinal, content, digest)
        self._check_latest_outputs()

    def _check_latest_outputs(self):
        from genosapiens.wgs_migration import matches_migrated_identity
        latest = self.connection.execute(
            "SELECT ordinal FROM commits ORDER BY ordinal DESC LIMIT 1"
        ).fetchone()
        if latest is None:
            return
        for name, metadata_bytes in self.connection.execute(
            "SELECT name,metadata_bytes FROM committed_outputs WHERE commit_ordinal=?",
            latest,
        ):
            metadata = _json(metadata_bytes, 4096)
            _require(type(metadata) is dict and set(metadata) == {
                "name", "sha256", "size_bytes", "publication",
            } and metadata["name"] == name and _digest(metadata["sha256"])
                     and _count(metadata["size_bytes"], 1))
            with _directory(self.stage_root) as fd:
                actual = os.stat(name, dir_fd=fd, follow_symlinks=False)
                directory = os.fstat(fd)
            _require(stat.S_ISREG(actual.st_mode) and stat.S_IMODE(actual.st_mode) == 0o600
                     and actual.st_uid == os.getuid() and actual.st_nlink == 1
                     and (metadata["publication"] == [
                         directory.st_dev, directory.st_ino, actual.st_dev,
                         actual.st_ino, actual.st_size, actual.st_mtime_ns,
                         actual.st_ctime_ns,
                     ] or matches_migrated_identity(self.stage_root / name,
                         metadata["publication"], kind="publication",
                         sha256=metadata["sha256"], size_bytes=metadata["size_bytes"])))

    def committed(self, partition_id):
        """Read one committed partition without scanning previous partitions."""
        _require_lock(self.stage_root, self.binding.stage, self.lock)
        _require(_name(partition_id))
        self._recover_head()
        row = self.connection.execute(
            "SELECT commits.ordinal,cursor_bytes,outputs_bytes,counts_bytes,event_bytes "
            "FROM commits JOIN events ON events.ordinal=commits.ordinal "
            "WHERE commits.partition_id=?", (partition_id,),
        ).fetchone()
        if row is None:
            return None
        ordinal, cursor_bytes, outputs_bytes, counts_bytes, event_bytes = row
        completion = load_model(event_bytes, JournalEntry)
        _require(completion.ordinal == ordinal and completion.event_type == "partition_completed"
                 and completion.partition_id == partition_id
                 and completion.stage == self.binding.stage
                 and completion.parents == self.binding.parents)
        cursor = _json(cursor_bytes, 4096)
        names = _json(outputs_bytes, 4096)
        counts = _json(counts_bytes, 4096)
        _require(type(cursor) is list and 0 < len(cursor) <= 4
                 and all(_count(value) for value in cursor)
                 and type(names) is list and 0 < len(names) <= 64
                 and len(set(names)) == len(names) and all(_name(name) for name in names)
                 and type(counts) is dict and all(
                     _name(key) and _count(value) for key, value in counts.items()
                 ))
        outputs = []
        with _directory(self.stage_root) as fd:
            for name in names:
                metadata_row = self.connection.execute(
                    "SELECT metadata_bytes FROM committed_outputs "
                    "WHERE name=? AND commit_ordinal=?", (name, ordinal),
                ).fetchone()
                _require(metadata_row is not None)
                metadata = _json(metadata_row[0], 4096)
                _require(type(metadata) is dict and set(metadata) == {
                    "name", "sha256", "size_bytes", "publication",
                } and metadata["name"] == name)
                artifact = Artifact(name, metadata["sha256"], metadata["size_bytes"])
                proof = metadata["publication"]
                _require(type(proof) is list and len(proof) == 7
                         and all(_count(value) for value in proof))
                from genosapiens.wgs_migration import matches_migrated_identity
                _require(matches_migrated_identity(self.stage_root / name, proof,
                    kind="publication", sha256=artifact.sha256, size_bytes=artifact.size_bytes))
                _remember_published_artifact(fd, artifact)
                _require_published_artifact(fd, artifact)
                outputs.append(artifact)
        _require(completion.artifacts == (outputs[-1],))
        return StageCommittedPartition(completion, tuple(cursor), tuple(outputs), counts)

    def open_attempt(self, partition_id):
        """Return the durable incomplete attempt, if one exists."""
        _require_lock(self.stage_root, self.binding.stage, self.lock)
        _require(_name(partition_id))
        row = self.connection.execute(
            "SELECT event_bytes FROM events WHERE partition_id=? "
            "AND event_type='attempt_started' AND attempt_id NOT IN "
            "(SELECT attempt_id FROM events WHERE event_type IN "
            "('partition_completed','attempt_failed','attempt_interrupted')) "
            "ORDER BY ordinal DESC LIMIT 1", (partition_id,),
        ).fetchone()
        return None if row is None else load_model(row[0], JournalEntry)

    def _save_event(self, entry, *, commit=None, reducer=None,
                    prepared_transaction=False):
        self._recover_head()
        latest = self._latest()
        _require(entry.ordinal == (1 if latest is None else latest[0] + 1)
                 and entry.previous_entry_sha256 == (
                     None if latest is None else latest[2]
                 ))
        connection = self.connection
        try:
            _require(type(prepared_transaction) is bool
                     and connection.in_transaction is prepared_transaction)
            if not prepared_transaction:
                connection.execute("BEGIN IMMEDIATE")
            if reducer is not None:
                reducer(connection)
            content = entry.as_bytes()
            connection.execute(
                "INSERT INTO events VALUES (?,?,?,?,?,?)",
                (entry.ordinal, content, entry.sha256, entry.event_type,
                 entry.partition_id, entry.attempt_id),
            )
            if commit is not None:
                cursor, outputs, counts = commit
                output_rows = []
                for artifact in outputs:
                    proof = getattr(artifact, "_publication", None)
                    metadata = canonical({
                        "name": artifact.name, "sha256": artifact.sha256,
                        "size_bytes": artifact.size_bytes,
                        "publication": proof,
                    })
                    output_rows.append((artifact.name, entry.ordinal, metadata))
                connection.execute(
                    "INSERT INTO commits VALUES (?,?,?,?,?)",
                    (entry.ordinal, entry.partition_id, canonical(cursor),
                     canonical([item.name for item in outputs]), canonical(counts)),
                )
                connection.executemany(
                    "INSERT INTO committed_outputs VALUES (?,?,?)", output_rows,
                )
            connection.execute("COMMIT")
        except Exception:
            connection.execute("ROLLBACK") if connection.in_transaction else None
            raise LedgerError("cannot commit WGS stage state") from None
        self._publish_event_and_head(entry.ordinal, entry.as_bytes(), entry.sha256)
        return entry

    def begin(self, *, partition_id, scientific, max_attempts,
              attempt_id=None, process_id=None):
        _require_lock(self.stage_root, self.binding.stage, self.lock)
        _require(_name(partition_id) and type(scientific) is bool
                 and _count(max_attempts, 1))
        starts = self.connection.execute(
            "SELECT attempt_id FROM events WHERE event_type='attempt_started' "
            "AND partition_id=? ORDER BY ordinal", (partition_id,),
        ).fetchall()
        _require(len(starts) < max_attempts)
        for (prior_attempt,) in starts:
            terminal = self.connection.execute(
                "SELECT 1 FROM events WHERE attempt_id=? AND event_type IN "
                "('partition_completed','attempt_failed','attempt_interrupted') LIMIT 1",
                (prior_attempt,),
            ).fetchone()
            _require(terminal is not None)
        _require(self.connection.execute(
            "SELECT 1 FROM commits WHERE partition_id=?", (partition_id,),
        ).fetchone() is None)
        latest = self._latest()
        _require(latest is None or load_model(latest[1], JournalEntry).event_type
                 != "stage_completed")
        entry = JournalEntry(
            self.binding.stage, 1 if latest is None else latest[0] + 1,
            "attempt_started", attempt_id or "attempt-" + uuid.uuid4().hex,
            len(starts) + 1, process_id or "process-" + uuid.uuid4().hex,
            partition_id, scientific, self.binding.stage_parent_sha256,
            self.binding.parents, (), None, None if latest is None else latest[2],
        )
        return self._save_event(entry)

    def commit_partition(self, start, *, cursor, outputs, counts, reducer=None,
                         prepared_transaction=False):
        _require_lock(self.stage_root, self.binding.stage, self.lock)
        _require(type(start) is JournalEntry and start.event_type == "attempt_started"
                 and type(cursor) is tuple and 0 < len(cursor) <= 4
                 and all(_count(value) for value in cursor)
                 and type(outputs) is tuple and bool(outputs)
                 and len(outputs) <= 64 and all(type(value) is Artifact for value in outputs)
                 and type(counts) is dict and all(
                     _name(key) and _count(value) for key, value in counts.items()
                 ))
        _require(self.connection.execute(
            "SELECT event_bytes FROM events WHERE attempt_id=? AND event_type='attempt_started'",
            (start.attempt_id,),
        ).fetchone() == (start.as_bytes(),))
        _require(self.connection.execute(
            "SELECT 1 FROM events WHERE attempt_id=? AND event_type!='attempt_started'",
            (start.attempt_id,),
        ).fetchone() is None)
        prior = self.connection.execute(
            "SELECT cursor_bytes FROM commits ORDER BY ordinal DESC LIMIT 1"
        ).fetchone()
        _require(prior is None or tuple(_json(prior[0])) < cursor)
        with _directory(self.stage_root) as fd:
            for artifact in outputs:
                _require_published_artifact(fd, artifact)
        latest = self._latest()
        entry = JournalEntry(
            start.stage, latest[0] + 1, "partition_completed",
            start.attempt_id, start.attempt_number, start.process_id,
            start.partition_id, start.scientific, start.stage_parent_sha256,
            start.parents, (outputs[-1],), None, latest[2],
        )
        return self._save_event(
            entry, commit=(cursor, outputs, counts), reducer=reducer,
            prepared_transaction=prepared_transaction,
        )

    def fail(self, start, *, failure_code, interrupted=False):
        _require_lock(self.stage_root, self.binding.stage, self.lock)
        _require(type(start) is JournalEntry and start.event_type == "attempt_started"
                 and _name(failure_code) and type(interrupted) is bool)
        _require(self.connection.execute(
            "SELECT event_bytes FROM events WHERE attempt_id=? AND event_type='attempt_started'",
            (start.attempt_id,),
        ).fetchone() == (start.as_bytes(),))
        _require(self.connection.execute(
            "SELECT 1 FROM events WHERE attempt_id=? AND event_type!='attempt_started'",
            (start.attempt_id,),
        ).fetchone() is None)
        latest = self._latest()
        entry = JournalEntry(
            start.stage, latest[0] + 1,
            "attempt_interrupted" if interrupted else "attempt_failed",
            start.attempt_id, start.attempt_number, start.process_id,
            start.partition_id, start.scientific, start.stage_parent_sha256,
            start.parents, (), failure_code, latest[2],
        )
        return self._save_event(entry)

    def finish(self, *, expected_partitions, output, process_id, closure_check):
        """Seal a stage only after every named partition and reducer closes."""
        _require_lock(self.stage_root, self.binding.stage, self.lock)
        _require(type(expected_partitions) is tuple
                 and len(expected_partitions) <= MAX_SHARDS
                 and len(set(expected_partitions)) == len(expected_partitions)
                 and all(_name(value) for value in expected_partitions)
                 and type(output) is Artifact
                 and type(process_id) is str
                 and re.fullmatch(r"process-[0-9a-f]{32}", process_id)
                 and callable(closure_check))
        self._recover_head()
        latest = self._latest()
        _require(latest is None or load_model(latest[1], JournalEntry).event_type
                 != "stage_completed")
        _require(latest is not None or not expected_partitions)
        committed = {
            item[0] for item in self.connection.execute(
                "SELECT partition_id FROM commits"
            )
        }
        _require(committed == set(expected_partitions))
        _require(self.connection.execute(
            "SELECT 1 FROM events AS starts WHERE starts.event_type='attempt_started' "
            "AND NOT EXISTS (SELECT 1 FROM events AS ends WHERE "
            "ends.attempt_id=starts.attempt_id AND ends.event_type IN "
            "('partition_completed','attempt_failed','attempt_interrupted')) LIMIT 1"
        ).fetchone() is None)
        _require(closure_check(self.connection) is True)
        with _directory(self.stage_root) as fd:
            _require_published_artifact(fd, output)
        entry = JournalEntry(
            self.binding.stage, 1 if latest is None else latest[0] + 1,
            "stage_completed", None, None,
            process_id, None, False, self.binding.stage_parent_sha256,
            self.binding.parents, (output,), None,
            None if latest is None else latest[2],
        )
        return self._save_event(entry)

    def archive_orphans(self, *, names, attempt_id):
        """Move only declared uncommitted output names to private history."""
        _require_lock(self.stage_root, self.binding.stage, self.lock)
        _require(type(names) is tuple and len(names) <= 64 and _name(attempt_id)
                 and len(names) == len(set(names)) and all(_name(name) for name in names))
        present = []
        with _directory(self.stage_root) as fd:
            for name in names:
                _require(self.connection.execute(
                    "SELECT 1 FROM committed_outputs WHERE name=?", (name,),
                ).fetchone() is None)
                try:
                    metadata = os.stat(name, dir_fd=fd, follow_symlinks=False)
                except FileNotFoundError:
                    continue
                _require(stat.S_ISREG(metadata.st_mode)
                         and stat.S_IMODE(metadata.st_mode) == 0o600
                         and metadata.st_uid == os.getuid() and metadata.st_nlink == 1)
                present.append(name)
            history_name = f"{self.binding.stage}-orphan-{attempt_id}"
            if not present:
                return ()
            try:
                try:
                    os.mkdir(history_name, 0o700, dir_fd=fd)
                except FileExistsError:
                    history_metadata = os.stat(
                        history_name, dir_fd=fd, follow_symlinks=False,
                    )
                    _require(stat.S_ISDIR(history_metadata.st_mode)
                             and stat.S_IMODE(history_metadata.st_mode) == 0o700
                             and history_metadata.st_uid == os.getuid())
                history_fd = os.open(
                    history_name, os.O_RDONLY | os.O_DIRECTORY | os.O_NOFOLLOW,
                    dir_fd=fd,
                )
                try:
                    for name in present:
                        try:
                            os.stat(name, dir_fd=history_fd, follow_symlinks=False)
                        except FileNotFoundError:
                            pass
                        else:
                            raise LedgerError("conflicting orphan WGS output")
                        os.rename(name, name, src_dir_fd=fd, dst_dir_fd=history_fd)
                    os.fsync(history_fd)
                    os.fsync(fd)
                finally:
                    os.close(history_fd)
            except OSError:
                raise LedgerError("cannot preserve orphan WGS output") from None
        return tuple(present)


def _journal_state(entries, *, stage, parent_sha256):
    previous = None
    fixed_parents = None
    attempts = {}
    starts_by_partition = {}
    completed_partitions = set()
    stage_complete = False
    for ordinal, entry in enumerate(entries, 1):
        _require(type(entry) is JournalEntry and entry.stage == stage
                 and entry.stage_parent_sha256 == parent_sha256 and entry.ordinal == ordinal
                 and entry.previous_entry_sha256 == previous and not stage_complete)
        if fixed_parents is None:
            fixed_parents = entry.parents
        else:
            _require(entry.parents == fixed_parents)
        if entry.event_type == "attempt_started":
            _require(entry.attempt_id not in attempts
                     and entry.partition_id not in completed_partitions
                     and not any(start.partition_id == entry.partition_id and terminal is None
                                 for start, terminal in attempts.values()))
            starts = starts_by_partition.setdefault(entry.partition_id, [])
            _require(entry.attempt_number == len(starts) + 1)
            starts.append(entry)
            attempts[entry.attempt_id] = (entry, None)
        elif entry.event_type in {
            "partition_completed", "attempt_failed", "attempt_interrupted",
        }:
            _require(entry.attempt_id in attempts)
            start, terminal = attempts[entry.attempt_id]
            _require(terminal is None and entry.attempt_number == start.attempt_number
                     and entry.partition_id == start.partition_id
                     and entry.process_id == start.process_id
                     and entry.scientific == start.scientific)
            attempts[entry.attempt_id] = (start, entry)
            if entry.event_type == "partition_completed":
                _require(entry.partition_id not in completed_partitions)
                completed_partitions.add(entry.partition_id)
        else:
            _require(not any(terminal is None for _start, terminal in attempts.values()))
            stage_complete = True
        previous = entry.sha256
    return {
        "attempts": attempts, "starts_by_partition": starts_by_partition,
        "completed_partitions": frozenset(completed_partitions),
        "stage_complete": stage_complete, "parents": fixed_parents,
    }


def append_journal(root, entry, *, lock):
    """Append one authenticated event while holding the stage mutation lock."""
    _require(type(entry) is JournalEntry)
    _require_lock(root, entry.stage, lock)
    prefix = read_journal(root, stage=entry.stage, parent_sha256=entry.stage_parent_sha256)
    _journal_state((*prefix, entry), stage=entry.stage,
                   parent_sha256=entry.stage_parent_sha256)
    with _directory(root) as fd:
        for artifact in entry.artifacts:
            _authenticate(fd, artifact)
    return publish_bytes(root, f"{entry.stage}-{entry.ordinal:06d}.journal.json",
                         entry.as_bytes())


def read_journal_head(root, *, stage, expected_head_sha256, parent_sha256=None):
    """Read externally bound terminal metadata without replaying completed work.

    Producers and recovery retain read_journal's complete validation. Completed
    consumers use this head and validate referenced evidence when they use it.
    """
    _require(_name(stage) and len(stage)<=80 and _digest(expected_head_sha256)
             and (parent_sha256 is None or _digest(parent_sha256)))
    with _directory(root) as fd:
        # append_journal publishes numbered events; StageCommitState also has a
        # head pointer. Both current producers have externally bound head hashes.
        entries=os.listdir(fd)
        names=[name for name in entries if re.fullmatch(re.escape(stage)+r'-[0-9]{6}\.journal\.json',name)]
        _require(0<len(names)<=MAX_JOURNAL_EVENTS)
        name=max(names);ordinal=int(name[-19:-13])
        if f'{stage}-head.json' in entries:
            with _read(fd,f'{stage}-head.json') as (stream,metadata):
                _require(metadata.st_size<=4096)
                head=_json(stream.read(4097),4096)
            _require(head==dict(schema='genosapiens.wgs-stage-journal-head',version=1,
                stage=stage,ordinal=ordinal,sha256=expected_head_sha256))
        with _read(fd,name) as (stream,metadata):
            _require(metadata.st_size<=MAX_MANIFEST_BYTES)
            entry=load_model(stream.read(MAX_MANIFEST_BYTES+1),JournalEntry)
    _require(entry.sha256==expected_head_sha256 and entry.ordinal==ordinal
        and entry.stage==stage and (parent_sha256 is None or entry.stage_parent_sha256==parent_sha256))
    return entry


def read_journal(root, *, stage, parent_sha256, expected_head_sha256=None):
    _require(_name(stage) and len(stage) <= 80 and _digest(parent_sha256))
    _require(expected_head_sha256 is None or _digest(expected_head_sha256))
    prefix = []
    with _directory(root) as fd:
        pattern = re.compile(re.escape(stage) + r"-([0-9]{6})\.journal\.json")
        names = sorted(name for name in os.listdir(fd) if pattern.fullmatch(name))
        _require(len(names) <= MAX_JOURNAL_EVENTS)
        for ordinal, name in enumerate(names, 1):
            _require(name == f"{stage}-{ordinal:06d}.journal.json")
            with _read(fd, name) as (stream, metadata):
                _require(metadata.st_size <= MAX_MANIFEST_BYTES)
                entry = load_model(stream.read(MAX_MANIFEST_BYTES + 1), JournalEntry)
            for artifact in entry.artifacts:
                _authenticate(fd, artifact)
            prefix.append(entry)
    result = tuple(prefix)
    _journal_state(result, stage=stage, parent_sha256=parent_sha256)
    if expected_head_sha256 is not None:
        _require(bool(result) and result[-1].sha256 == expected_head_sha256)
    return result


def begin_attempt(root, *, lock, stage, partition_id, parent_sha256, parents,
                  scientific, max_attempts, attempt_id=None, process_id=None):
    """Durably count an attempt before its child or scientific kernel is invoked."""
    _require_lock(root, stage, lock)
    _require(_name(partition_id) and type(scientific) is bool and _count(max_attempts, 1))
    _parents(parents, required={"stage_parent", "partition_plan", "execution_limits",
                                "executable_stack", "runtime_image"})
    prefix = read_journal(root, stage=stage, parent_sha256=parent_sha256)
    state = _journal_state(prefix, stage=stage, parent_sha256=parent_sha256)
    _require(not state["stage_complete"] and partition_id not in state["completed_partitions"])
    starts = state["starts_by_partition"].get(partition_id, [])
    _require(len(starts) < max_attempts
             and not any(start.partition_id == partition_id and terminal is None
                         for start, terminal in state["attempts"].values()))
    attempt_id = attempt_id or "attempt-" + uuid.uuid4().hex
    process_id = process_id or "process-" + uuid.uuid4().hex
    entry = JournalEntry(
        stage, len(prefix) + 1, "attempt_started", attempt_id, len(starts) + 1,
        process_id, partition_id, scientific, parent_sha256, parents, (), None,
        prefix[-1].sha256 if prefix else None,
    )
    append_journal(root, entry, lock=lock)
    return entry


def _terminal_attempt(root, *, lock, start, event_type, artifacts, failure_code):
    _require(type(start) is JournalEntry and start.event_type == "attempt_started")
    _require_lock(root, start.stage, lock)
    prefix = read_journal(root, stage=start.stage,
                          parent_sha256=start.stage_parent_sha256)
    entry = JournalEntry(
        start.stage, len(prefix) + 1, event_type, start.attempt_id, start.attempt_number,
        start.process_id, start.partition_id, start.scientific, start.stage_parent_sha256,
        start.parents, artifacts, failure_code, prefix[-1].sha256 if prefix else None,
    )
    append_journal(root, entry, lock=lock)
    return entry


def complete_attempt(root, *, lock, start, output):
    _require(type(output) is Artifact)
    return _terminal_attempt(root, lock=lock, start=start, event_type="partition_completed",
                             artifacts=(output,), failure_code=None)


def fail_attempt(root, *, lock, start, failure_code, interrupted=False, evidence=()):
    _require(type(interrupted) is bool and type(evidence) is tuple)
    return _terminal_attempt(
        root, lock=lock, start=start,
        event_type="attempt_interrupted" if interrupted else "attempt_failed",
        artifacts=evidence, failure_code=failure_code,
    )


def complete_stage(root, *, lock, stage, parent_sha256, parents, output,
                   process_id=None):
    _require_lock(root, stage, lock)
    _require(type(output) is Artifact)
    prefix = read_journal(root, stage=stage, parent_sha256=parent_sha256)
    entry = JournalEntry(
        stage, len(prefix) + 1, "stage_completed", None, None,
        process_id or "process-" + uuid.uuid4().hex, None, False, parent_sha256,
        parents, (output,), None, prefix[-1].sha256 if prefix else None,
    )
    append_journal(root, entry, lock=lock)
    return entry


def journal_attempts_used(entries, partition_id):
    _require(type(entries) is tuple and _name(partition_id))
    return sum(entry.event_type == "attempt_started" and entry.partition_id == partition_id
               for entry in entries)


def verify_ledger(root, manifest):
    """Authenticate every manifest artifact and shard without accepting a prefix."""
    _require(type(manifest) is LedgerManifest)
    with _directory(root) as fd:
        for artifact in (manifest.header, *manifest.summary.reports):
            _authenticate(fd, artifact)
    for shard in manifest.shards:
        verify_shard(root, shard)


def publish_manifest(root, manifest, *, name="manifest.json"):
    _require(type(manifest) is LedgerManifest and name == "manifest.json")
    with _directory(root) as fd:
        for artifact in (manifest.header, *manifest.summary.reports):
            _authenticate(fd, artifact)
        for shard in manifest.shards:
            _require_published_shard(fd, shard)
    return publish_bytes(root, name, manifest.as_bytes())


def load_manifest(root, artifact):
    _require(type(artifact) is Artifact and artifact.name == "manifest.json")
    return load_model(read_artifact(root, artifact), LedgerManifest)


def iter_ledger(root, manifest):
    _require(type(manifest) is LedgerManifest)
    with _directory(root) as fd:
        for artifact in (manifest.header, *manifest.summary.reports):
            _authenticate(fd, artifact)
    try:
        for shard in manifest.shards:
            yield from iter_shard(root, shard)
    except (OSError, EOFError, zlib.error):
        raise LedgerError("cannot authenticate WGS ledger shard") from None


SQLITE_CACHE_KIB = 64 * 1024
SQLITE_CURSOR_BATCH_ROWS = 1024
_INDEX_DDL = """
PRAGMA application_id=1196972364;
PRAGMA user_version=2;
CREATE TABLE binding (
    singleton INTEGER PRIMARY KEY NOT NULL CHECK(singleton=1),
    closure BLOB NOT NULL
);
CREATE TABLE records (
    record_number INTEGER PRIMARY KEY NOT NULL,
    source_record_id TEXT UNIQUE NOT NULL,
    canonical_row BLOB NOT NULL
);
CREATE TABLE alleles (
    source_allele_id TEXT PRIMARY KEY NOT NULL,
    record_number INTEGER NOT NULL REFERENCES records(record_number),
    allele_number INTEGER NOT NULL
);
CREATE UNIQUE INDEX allele_record_slot ON alleles(record_number, allele_number);
CREATE TABLE normalization_map (
    source_allele_id TEXT PRIMARY KEY NOT NULL REFERENCES alleles(source_allele_id),
    normalized_identity_id TEXT,
    preparation_disposition TEXT NOT NULL
);
CREATE INDEX normalization_identity ON normalization_map(normalized_identity_id, source_allele_id);
CREATE TABLE decisions (
    source_allele_id TEXT PRIMARY KEY NOT NULL REFERENCES normalization_map(source_allele_id),
    normalized_identity_id TEXT,
    preparation_disposition TEXT NOT NULL,
    early_filter_disposition TEXT,
    scope_disposition TEXT NOT NULL
);
CREATE TABLE working_identities (
    identity_id TEXT PRIMARY KEY NOT NULL,
    route TEXT NOT NULL
);
CREATE TABLE responses (
    response_id TEXT PRIMARY KEY NOT NULL,
    identity_id TEXT NOT NULL REFERENCES working_identities(identity_id),
    response_ordinal INTEGER NOT NULL,
    response_sha256 TEXT NOT NULL
);
CREATE UNIQUE INDEX response_identity_ordinal ON responses(identity_id, response_ordinal);
CREATE TABLE context_memberships (
    context_id TEXT NOT NULL,
    member_ordinal INTEGER NOT NULL,
    identity_id TEXT NOT NULL REFERENCES working_identities(identity_id),
    PRIMARY KEY(context_id, member_ordinal)
);
CREATE UNIQUE INDEX context_identity ON context_memberships(context_id, identity_id);
"""
_INDEX_TABLE_COLUMNS = {
    "binding": (("singleton", "INTEGER", 1, 1), ("closure", "BLOB", 1, 0)),
    "records": (("record_number", "INTEGER", 1, 1), ("source_record_id", "TEXT", 1, 0),
                ("canonical_row", "BLOB", 1, 0)),
    "alleles": (("source_allele_id", "TEXT", 1, 1), ("record_number", "INTEGER", 1, 0),
                ("allele_number", "INTEGER", 1, 0)),
    "normalization_map": (("source_allele_id", "TEXT", 1, 1),
                          ("normalized_identity_id", "TEXT", 0, 0),
                          ("preparation_disposition", "TEXT", 1, 0)),
    "decisions": (("source_allele_id", "TEXT", 1, 1),
                  ("normalized_identity_id", "TEXT", 0, 0),
                  ("preparation_disposition", "TEXT", 1, 0),
                  ("early_filter_disposition", "TEXT", 0, 0),
                  ("scope_disposition", "TEXT", 1, 0)),
    "working_identities": (("identity_id", "TEXT", 1, 1), ("route", "TEXT", 1, 0)),
    "responses": (("response_id", "TEXT", 1, 1), ("identity_id", "TEXT", 1, 0),
                  ("response_ordinal", "INTEGER", 1, 0),
                  ("response_sha256", "TEXT", 1, 0)),
    "context_memberships": (("context_id", "TEXT", 1, 1),
                            ("member_ordinal", "INTEGER", 1, 2),
                            ("identity_id", "TEXT", 1, 0)),
}
_INDEX_OBJECTS = {
    *(('table', name) for name in _INDEX_TABLE_COLUMNS),
    ("index", "allele_record_slot"), ("index", "normalization_identity"),
    ("index", "response_identity_ordinal"), ("index", "context_identity"),
}
_INDEX_FOREIGN_KEYS = {
    "alleles": {("records", "record_number", "record_number")},
    "normalization_map": {("alleles", "source_allele_id", "source_allele_id")},
    "decisions": {("normalization_map", "source_allele_id", "source_allele_id")},
    "responses": {("working_identities", "identity_id", "identity_id")},
    "context_memberships": {("working_identities", "identity_id", "identity_id")},
}


def _cursor_rows(cursor):
    while True:
        batch = cursor.fetchmany(SQLITE_CURSOR_BATCH_ROWS)
        _require(len(batch) <= SQLITE_CURSOR_BATCH_ROWS)
        if not batch:
            return
        yield from batch


def _configure_index(connection, *, writable):
    connection.enable_load_extension(False)
    connection.execute("PRAGMA foreign_keys=ON")
    connection.execute(f"PRAGMA cache_size=-{SQLITE_CACHE_KIB}")
    connection.execute("PRAGMA temp_store=FILE")
    connection.execute("PRAGMA mmap_size=0")
    connection.execute("PRAGMA trusted_schema=OFF")
    if writable:
        connection.execute("PRAGMA journal_mode=DELETE")
    else:
        connection.execute("PRAGMA query_only=ON")


def _insert_projection(connection, rows, cls, *, sql, values, key, parent_sha256):
    digest = hashlib.sha256()
    count = 0
    previous = None
    for row in rows:
        _require(type(row) is cls)
        row_key = key(row)
        _require(previous is None or row_key > previous)
        previous = row_key
        connection.execute(sql, values(row))
        digest.update(row.as_bytes())
        count += 1
        _require(count <= 2**63 - 1)
    return IndexTableBinding(parent_sha256, count, digest.hexdigest())


def _verify_index_schema(connection):
    objects = set(_cursor_rows(connection.execute(
        "SELECT type,name FROM sqlite_master WHERE name NOT LIKE 'sqlite_%' ORDER BY type,name"
    )))
    _require(objects == _INDEX_OBJECTS)
    for table, expected in _INDEX_TABLE_COLUMNS.items():
        actual = tuple((row[1], row[2], row[3], row[5]) for row in _cursor_rows(
            connection.execute(f"PRAGMA table_info({table})")
        ))
        _require(actual == expected)
        foreign = {(row[2], row[3], row[4]) for row in _cursor_rows(
            connection.execute(f"PRAGMA foreign_key_list({table})")
        )}
        _require(foreign == _INDEX_FOREIGN_KEYS.get(table, set()))
    _require(connection.execute("PRAGMA application_id").fetchone() == (1196972364,)
             and connection.execute("PRAGMA user_version").fetchone() == (2,))


def _verify_logical_table(cursor, cls, binding, factory):
    digest = hashlib.sha256()
    count = 0
    for values in _cursor_rows(cursor):
        row = factory(values)
        _require(type(row) is cls)
        digest.update(row.as_bytes())
        count += 1
    _require(count == binding.rows and digest.hexdigest() == binding.logical_sha256)


def build_index(
    root, manifest, *, name="ledger.sqlite", normalization_rows=(), decision_rows=(),
    working_identity_rows=(), response_rows=(), context_membership_rows=(),
    normalization_parent_sha256=None, decision_parent_sha256=None,
    working_set_parent_sha256=None, response_parent_sha256=None,
    context_parent_sha256=None, _replace=False, _lock=None, _stage=None,
):
    """Build and reconcile a disposable index from canonical streamed rows."""
    _require(type(manifest) is LedgerManifest and _name(name))
    if _replace:
        _require(type(_stage) is str)
        _require_lock(root, _stage, _lock)
    temporary = ".pending-index-" + uuid.uuid4().hex
    connection = None
    closure = None
    with _directory(root) as fd:
        try:
            descriptor = os.open(temporary, os.O_WRONLY | os.O_CREAT | os.O_EXCL | os.O_NOFOLLOW,
                                 0o600, dir_fd=fd)
            os.fchmod(descriptor, 0o600)
            os.close(descriptor)
            connection = sqlite3.connect(f"/proc/self/fd/{fd}/{temporary}")
            _configure_index(connection, writable=True)
            connection.executescript(_INDEX_DDL)
            for row in iter_ledger(root, manifest):
                connection.execute("INSERT INTO records VALUES (?,?,?)",
                                   (row.record_number, row.source_record_id, row.as_bytes()))
                connection.executemany(
                    "INSERT INTO alleles VALUES (?,?,?)",
                    ((value, row.record_number, number)
                     for number, value in enumerate(row.source_allele_ids, 1)),
                )
            normalization = _insert_projection(
                connection, normalization_rows, NormalizationIndexRow,
                sql="INSERT INTO normalization_map VALUES (?,?,?)",
                values=lambda row: (row.source_allele_id, row.normalized_identity_id,
                                    row.preparation_disposition),
                key=lambda row: (row.source_allele_id,),
                parent_sha256=normalization_parent_sha256,
            )
            decisions = _insert_projection(
                connection, decision_rows, DecisionIndexRow,
                sql="INSERT INTO decisions VALUES (?,?,?,?,?)",
                values=lambda row: (row.source_allele_id, row.normalized_identity_id,
                                    row.preparation_disposition, row.early_filter_disposition, row.scope_disposition),
                key=lambda row: (row.source_allele_id,), parent_sha256=decision_parent_sha256,
            )
            working = _insert_projection(
                connection, working_identity_rows, WorkingIdentityIndexRow,
                sql="INSERT INTO working_identities VALUES (?,?)",
                values=lambda row: (row.identity_id, row.route),
                key=lambda row: (row.identity_id,), parent_sha256=working_set_parent_sha256,
            )
            responses = _insert_projection(
                connection, response_rows, ResponseIndexRow,
                sql="INSERT INTO responses VALUES (?,?,?,?)",
                values=lambda row: (row.response_id, row.identity_id, row.response_ordinal,
                                    row.response_sha256),
                key=lambda row: (row.identity_id, row.response_ordinal, row.response_id),
                parent_sha256=response_parent_sha256,
            )
            contexts = _insert_projection(
                connection, context_membership_rows, ContextMembershipIndexRow,
                sql="INSERT INTO context_memberships VALUES (?,?,?)",
                values=lambda row: (row.context_id, row.member_ordinal, row.identity_id),
                key=lambda row: (row.context_id, row.member_ordinal, row.identity_id),
                parent_sha256=context_parent_sha256,
            )
            closure = IndexClosure(manifest.sha256, normalization, decisions, working,
                                   responses, contexts)
            connection.execute("INSERT INTO binding VALUES (1,?)", (closure.as_bytes(),))
            connection.commit()
            connection.close()
            connection = None
            verify_index(root, manifest, closure, name=temporary, _temporary=True)
            descriptor = os.open(temporary, os.O_RDONLY | os.O_NOFOLLOW, dir_fd=fd)
            try:
                os.fsync(descriptor)
            finally:
                os.close(descriptor)
            if _replace:
                try:
                    existing = os.open(name, os.O_RDONLY | os.O_NOFOLLOW | os.O_NONBLOCK,
                                       dir_fd=fd)
                except FileNotFoundError:
                    existing = -1
                if existing >= 0:
                    try:
                        metadata = os.fstat(existing)
                        _require(stat.S_ISREG(metadata.st_mode) and metadata.st_uid == os.getuid()
                                 and stat.S_IMODE(metadata.st_mode) == 0o600
                                 and metadata.st_nlink == 1)
                    finally:
                        os.close(existing)
                os.rename(temporary, name, src_dir_fd=fd, dst_dir_fd=fd)
            else:
                os.link(temporary, name, src_dir_fd=fd, dst_dir_fd=fd, follow_symlinks=False)
                os.unlink(temporary, dir_fd=fd)
            os.fsync(fd)
            return closure
        except (OSError, sqlite3.Error):
            raise LedgerError("cannot build WGS ledger index") from None
        finally:
            if connection is not None:
                connection.close()
            for suffix in ("", "-journal", "-wal", "-shm"):
                try:
                    os.unlink(temporary + suffix, dir_fd=fd)
                except FileNotFoundError:
                    pass


def rebuild_index(root, manifest, *, lock, stage, **kwargs):
    return build_index(root, manifest, _replace=True, _lock=lock, _stage=stage, **kwargs)


def verify_index(root, manifest, closure, *, name="ledger.sqlite", _temporary=False):
    """Reconcile every logical row; SQLite page bytes are never authority."""
    _require(type(manifest) is LedgerManifest and type(closure) is IndexClosure
             and closure.source_manifest_sha256 == manifest.sha256)
    _require(_name(name) or _temporary and re.fullmatch(r"\.pending-index-[0-9a-f]{32}", name))
    connection = None
    with _directory(root) as fd:
        descriptor = -1
        try:
            descriptor = os.open(name, os.O_RDONLY | os.O_NOFOLLOW | os.O_NONBLOCK, dir_fd=fd)
            metadata = os.fstat(descriptor)
            _require(stat.S_ISREG(metadata.st_mode) and metadata.st_uid == os.getuid()
                     and stat.S_IMODE(metadata.st_mode) == 0o600 and metadata.st_nlink == 1)
            connection = sqlite3.connect(f"file:/proc/self/fd/{descriptor}?mode=ro&immutable=1", uri=True)
            _configure_index(connection, writable=False)
            integrity = connection.execute("PRAGMA integrity_check")
            _require(integrity.fetchone() == ("ok",) and integrity.fetchone() is None)
            _require(connection.execute("PRAGMA foreign_key_check").fetchone() is None)
            _verify_index_schema(connection)
            binding = connection.execute("SELECT singleton,closure FROM binding ORDER BY singleton")
            _require(binding.fetchone() == (1, closure.as_bytes()) and binding.fetchone() is None)
            records = iter(_cursor_rows(connection.execute(
                "SELECT record_number,source_record_id,canonical_row FROM records ORDER BY record_number"
            )))
            alleles = iter(_cursor_rows(connection.execute(
                "SELECT source_allele_id,record_number,allele_number FROM alleles "
                "ORDER BY record_number,allele_number"
            )))
            for row in iter_ledger(root, manifest):
                _require(next(records, None) == (row.record_number, row.source_record_id, row.as_bytes()))
                for number, value in enumerate(row.source_allele_ids, 1):
                    _require(next(alleles, None) == (value, row.record_number, number))
            _require(next(records, None) is None and next(alleles, None) is None)
            _verify_logical_table(
                connection.execute(
                    "SELECT source_allele_id,normalized_identity_id,preparation_disposition "
                    "FROM normalization_map ORDER BY source_allele_id"
                ), NormalizationIndexRow, closure.normalization,
                lambda values: NormalizationIndexRow(*values),
            )
            _verify_logical_table(
                connection.execute(
                    "SELECT source_allele_id,normalized_identity_id,preparation_disposition,"
                    "early_filter_disposition,scope_disposition FROM decisions ORDER BY source_allele_id"
                ), DecisionIndexRow, closure.decisions,
                lambda values: DecisionIndexRow(*values[:4], scope_disposition=values[4]),
            )
            _verify_logical_table(
                connection.execute("SELECT identity_id,route FROM working_identities ORDER BY identity_id"),
                WorkingIdentityIndexRow, closure.working_identities,
                lambda values: WorkingIdentityIndexRow(*values),
            )
            _verify_logical_table(
                connection.execute(
                    "SELECT identity_id,response_ordinal,response_id,response_sha256 FROM responses "
                    "ORDER BY identity_id,response_ordinal,response_id"
                ), ResponseIndexRow, closure.responses,
                lambda values: ResponseIndexRow(*values),
            )
            _verify_logical_table(
                connection.execute(
                    "SELECT context_id,member_ordinal,identity_id FROM context_memberships "
                    "ORDER BY context_id,member_ordinal,identity_id"
                ), ContextMembershipIndexRow, closure.context_memberships,
                lambda values: ContextMembershipIndexRow(*values),
            )
            if closure.normalization.parent_sha256 is not None:
                _require(connection.execute(
                    "SELECT 1 FROM alleles a LEFT JOIN normalization_map n USING(source_allele_id) "
                    "WHERE n.source_allele_id IS NULL LIMIT 1"
                ).fetchone() is None)
            if closure.decisions.parent_sha256 is not None:
                _require(connection.execute(
                    "SELECT 1 FROM normalization_map n LEFT JOIN decisions d USING(source_allele_id) "
                    "WHERE d.source_allele_id IS NULL OR d.normalized_identity_id IS NOT "
                    "n.normalized_identity_id OR d.preparation_disposition IS NOT "
                    "n.preparation_disposition LIMIT 1"
                ).fetchone() is None)
            _require(connection.execute(
                "SELECT 1 FROM working_identities w LEFT JOIN normalization_map n "
                "ON n.normalized_identity_id=w.identity_id WHERE w.route='sequence' "
                "AND n.source_allele_id IS NULL LIMIT 1"
            ).fetchone() is None)
            if closure.responses.parent_sha256 is not None:
                _require(connection.execute(
                    "SELECT 1 FROM working_identities w LEFT JOIN responses r "
                    "ON r.identity_id=w.identity_id WHERE w.route='sequence' "
                    "AND r.response_id IS NULL LIMIT 1"
                ).fetchone() is None)
            final = os.fstat(descriptor)
            _require((metadata.st_dev, metadata.st_ino, metadata.st_size,
                      metadata.st_mtime_ns, metadata.st_ctime_ns) ==
                     (final.st_dev, final.st_ino, final.st_size,
                      final.st_mtime_ns, final.st_ctime_ns))
        except (OSError, sqlite3.Error):
            raise LedgerError("invalid WGS ledger index") from None
        finally:
            if connection is not None:
                connection.close()
            if descriptor >= 0:
                os.close(descriptor)
