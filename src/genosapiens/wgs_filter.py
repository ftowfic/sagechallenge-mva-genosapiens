"""Immutable early-filter decisions and retained working-set closure."""

from __future__ import annotations

from dataclasses import asdict, dataclass
from collections import deque
from contextlib import nullcontext
import hashlib
import json
import os
from pathlib import Path
import sqlite3
import stat
from typing import Iterable, Iterator

from genosapiens.common_allele_filter import PrefilterDecision, PrefilterOccurrence
from genosapiens.wgs_ledger import (
    EARLY_FILTER_DISPOSITIONS,
    FILTER_DECISION_SCHEMA,
    FILTER_MANIFEST_SCHEMA,
    MAX_MANIFEST_BYTES,
    MAX_ROW_BYTES,
    MAX_SHARD_BYTES,
    MAX_SHARD_ROWS,
    PREPARATION_DISPOSITIONS,
    SCOPE_DISPOSITIONS,
    SOURCE_MANIFEST_SCHEMA,
    WORKING_SET_MANIFEST_SCHEMA,
    Artifact,
    FilterDecision,
    FilterManifest,
    FilterSummary,
    LedgerManifest,
    NamedCount,
    ParentIdentity,
    Shard,
    StageCommitState,
    StageRunBinding,
    WorkingSetManifest,
    WorkingSetRow,
    WorkingSetSummary,
    attest_existing_artifact,
    canonical,
    iter_filter_shard,
    iter_ledger,
    iter_shard,
    iter_working_set_shard,
    load_model,
    publish_bytes,
    publish_filter_shard,
    publish_working_set_shard,
    read_artifact,
    require_published_shards,
    stage_lock,
)
from genosapiens.wgs_preparation import SourceAccountingClosure
from genosapiens.wgs_failures import WgsFailure, failure_boundary, static_failure
from genosapiens.wgs_runtime import RuntimePhaseRecorder, runtime_phase
from genosapiens.wgs_source_routing import (
    SourceRoutingContext,
    WgsSourceRoutingError,
    classify_source_occurrence,
    decode_source_occurrence,
)


FILTER_EVIDENCE_SCHEMA = "genosapiens.common-allele-prefilter-evidence"
FILTER_SUMMARY_ARTIFACT = "summary.json"
MANIFEST_ARTIFACT = "manifest.json"
class WgsFilterError(RuntimeError):
    """Static filter publication/replay failure."""


class _ClassifiedFilterError(WgsFailure, WgsFilterError):
    def __init__(self, message, *, origin, operation):
        WgsFailure.__init__(
            self,
            origin,
            "malformed_evidence",
            operation=operation,
            exception_kind="project_failure",
        )
        RuntimeError.__init__(self, message)


def _digest(value: object) -> bool:
    return (
        isinstance(value, str)
        and len(value) == 64
        and all(character in "0123456789abcdef" for character in value)
    )


def _count(value: object, minimum: int = 0) -> bool:
    return isinstance(value, int) and not isinstance(value, bool) and value >= minimum


def _artifact_for(name: str, content: bytes) -> Artifact:
    return Artifact(name, hashlib.sha256(content).hexdigest(), len(content))


def _validate_source_routing(
    source_manifest: LedgerManifest,
    routing: SourceRoutingContext,
) -> None:
    if (
        not isinstance(routing, SourceRoutingContext)
        or routing.contract.parent_source_manifest_sha256 != source_manifest.sha256
        or routing.contract.reference_manifest_sha256
        != source_manifest.source.reference_manifest_sha256
    ):
        raise WgsFilterError("source routing context is invalid")


@dataclass(frozen=True, slots=True, repr=False)
class FilterPublication:
    manifest: FilterManifest
    manifest_artifact: Artifact
    summary_artifact: Artifact


@dataclass(frozen=True, slots=True, repr=False)
class WorkingSetPublication:
    manifest: WorkingSetManifest
    manifest_artifact: Artifact
    summary_artifact: Artifact


def ledger_filter_decision(
    occurrence: PrefilterOccurrence,
    decision: PrefilterDecision | None,
    *,
    parent_source_manifest_sha256: str,
    scope_decision,
    policy_sha256=None,
    resource_sha256=None,
    population_request_keys=0,
    clinical_query_keys=0,
) -> FilterDecision:
    """Convert a pure policy result without changing its disposition semantics."""

    from genosapiens.track1_early_scope import EarlyScopeDecision
    if (type(scope_decision) is not EarlyScopeDecision
            or scope_decision.source_allele_id != occurrence.source_allele_id
            or scope_decision.status not in SCOPE_DISPOSITIONS):
        raise WgsFilterError('filter scope decision binding is invalid')
    deferred = scope_decision.status in ('outside_scope', 'unresolved_scope')
    if (
        not isinstance(occurrence, PrefilterOccurrence)
        or deferred != (decision is None)
        or decision is not None and (not isinstance(decision, PrefilterDecision)
        or occurrence.source_allele_id != decision.source_allele_id
        or occurrence.record_number != decision.record_number
        or occurrence.allele_number != decision.allele_number
        or occurrence.normalized_allele_id != decision.normalized_allele_id)
        or not _digest(parent_source_manifest_sha256)
    ):
        raise WgsFilterError("prefilter decision binding is invalid")
    normalized = occurrence.preparation_disposition == "normalized"
    if decision is not None:
        policy_sha256 = decision.policy_sha256
        resource_sha256 = tuple(getattr(decision, key) for key in (
            'population_manifest_sha256','clinical_manifest_sha256',
            'exception_manifest_sha256','precision_manifest_sha256'))
    if (not _digest(policy_sha256) or type(resource_sha256) is not tuple
            or len(resource_sha256) != 4 or any(not _digest(v) for v in resource_sha256)
            or type(population_request_keys) is not int or population_request_keys not in (0,1)
            or type(clinical_query_keys) is not int or clinical_query_keys not in (0,1)
            or (deferred or not normalized) and (population_request_keys or clinical_query_keys)):
        raise WgsFilterError('filter resource or request binding is invalid')
    evidence = canonical({
        "schema": FILTER_EVIDENCE_SCHEMA,
        "version": 2,
        "classification": "controlled",
        "llm_policy": "no-llm",
        "occurrence": asdict(occurrence),
        "scope_decision": asdict(scope_decision),
        "decision": None if decision is None else asdict(decision),
        'resource_sha256': resource_sha256,
        'population_lookup': 'performed' if population_request_keys else 'not_performed',
        'clinical_lookup': 'performed' if clinical_query_keys else 'not_performed',
        'request_counts': dict(population_request_keys=population_request_keys, clinical_query_keys=clinical_query_keys),
    }).decode("ascii")
    return FilterDecision(
        occurrence.record_number,
        occurrence.allele_number,
        occurrence.source_allele_id,
        occurrence.normalized_allele_id if normalized else None,
        occurrence.preparation_disposition,
        decision.disposition if normalized and not deferred else None,
        None if normalized else ('preparation_'+occurrence.preparation_disposition if deferred else decision.primary_reason),
        scope_decision.reason_codes if deferred else decision.reason_codes,
        parent_source_manifest_sha256,
        policy_sha256,
        evidence,
        scope_disposition=scope_decision.status,
        scope_policy_sha256=scope_decision.scope_policy_sha256,
        scope_manifest_sha256=scope_decision.scope_manifest_sha256,
    )


class SourceDecisionCursor:
    """One source-owned decode cursor with a bounded independent decision oracle."""

    def __init__(self, source_root: Path, manifest: LedgerManifest,
                 routing: SourceRoutingContext, *, max_pending: int):
        if (not isinstance(source_root, Path) or not source_root.is_absolute()
                or type(manifest) is not LedgerManifest
                or type(routing) is not SourceRoutingContext
                or routing.contract.parent_source_manifest_sha256 != manifest.sha256
                or not _count(max_pending, 1)):
            raise WgsFilterError("source decision cursor binding is invalid")
        self.source_root = source_root
        self.manifest = manifest
        self.routing = routing
        self.max_pending = max_pending
        self._pending = deque()
        self._started = False
        self._completed = False
        self._resume_after = (0, 0)

    def resume_after(self, record_number: int, allele_number: int) -> None:
        """Position the decision oracle at the committed source occurrence."""
        if (self._started or self._completed or self._pending
                or not _count(record_number, 1) or not _count(allele_number, 1)
                or record_number > self.manifest.summary.total_records):
            raise WgsFilterError("source decision resume cursor is invalid")
        self._resume_after = (record_number, allele_number)

    def __iter__(self):
        if self._started:
            raise WgsFilterError("source decision cursor was reused")
        self._started = True
        if self._resume_after == (0, 0):
            source_rows = iter_ledger(self.source_root, self.manifest)
        else:
            source_rows = (
                row for shard in self.manifest.shards
                if shard.last_record >= self._resume_after[0]
                for row in iter_shard(self.source_root, shard)
            )
        for row in source_rows:
            for allele_number, source_id in enumerate(row.source_allele_ids, 1):
                if (row.record_number, allele_number) <= self._resume_after:
                    continue
                with failure_boundary(
                    "source_decode", "source_payload_decode",
                    known_errors=(WgsSourceRoutingError,),
                ):
                    facts = decode_source_occurrence(
                        row.payload_json,
                        source_allele_id=source_id,
                        source_manifest=self.manifest,
                        contract=self.routing.contract,
                        expected_record_number=row.record_number,
                        expected_source_record_id=row.source_record_id,
                        expected_source_allele_ids=row.source_allele_ids,
                    )
                with failure_boundary(
                    "source_classification", "source_disposition_classify",
                    known_errors=(WgsSourceRoutingError,),
                ):
                    route = classify_source_occurrence(
                        facts,
                        contract=self.routing.contract,
                        renderer=self.routing.renderer,
                    )
                if len(self._pending) >= self.max_pending:
                    raise WgsFilterError("source decision cursor exceeded its batch bound")
                self._pending.append((
                    row.record_number, facts.source.allele_number, source_id,
                    route.preparation_disposition, route.normalized_allele_id,
                ))
                yield facts, route
        self._completed = True

    def pop_expected(self):
        if not self._pending:
            raise WgsFilterError("filter decision has no source-owned occurrence")
        return self._pending.popleft()

    def require_complete(self):
        if not self._started or not self._completed or self._pending:
            raise WgsFilterError("source-owned decision coverage is incomplete")


def _replay_source_tuples(
    source_root: Path, manifest: LedgerManifest, routing: SourceRoutingContext,
):
    """Explicit low-level replay uses the same sole routing cursor as production."""
    cursor = SourceDecisionCursor(
        source_root, manifest, routing, max_pending=1,
    )
    for _facts, _route in cursor:
        yield cursor.pop_expected()
    cursor.require_complete()


@dataclass(frozen=True, slots=True)
class CompletedFilterBatch:
    decisions: tuple[FilterDecision, ...]
    population_request_keys: int
    clinical_query_keys: int

    def __post_init__(self):
        if (type(self.decisions) is not tuple or not self.decisions
                or any(type(row) is not FilterDecision for row in self.decisions)):
            raise WgsFilterError("completed filter batch is invalid")
        for key in ("population_request_keys", "clinical_query_keys"):
            total = getattr(self, key)
            if (type(total) is not int or not 0 <= total <= len(self.decisions)
                    or total != sum(json.loads(row.evidence_json)["request_counts"][key]
                                    for row in self.decisions)):
                raise WgsFilterError("completed filter batch counts differ")


@dataclass(frozen=True, slots=True)
class SourceBoundDecisionStream:
    source_cursor: SourceDecisionCursor
    decisions: Iterable[FilterDecision | CompletedFilterBatch]
    completed_batches: bool = False

    def __iter__(self):
        for item in self.decisions:
            if self.completed_batches:
                yield from item.decisions
            else:
                yield item


def _rows(cursor, size: int = 1024):
    while True:
        batch = cursor.fetchmany(size)
        if not batch:
            return
        yield from batch


class _FilterState:
    __slots__ = (
        "path", "expected", "source_manifest_sha256", "policy_sha256",
        "resource_sha256", "_connection", "_owns_connection", "_count",
        "scope_sha256",
        "reference_sha256",
        "_previous", "_sealed", "_closed",
    )

    def __init__(
        self, path: Path, *, expected: int, source_manifest_sha256: str,
        policy_sha256: str, resource_sha256: tuple[str, str, str, str],
        scope_sha256: tuple[str, str],
        reference_sha256: str,
        connection: sqlite3.Connection | None = None,
    ) -> None:
        if (
            not _count(expected, 1)
            or not _digest(source_manifest_sha256)
            or not _digest(policy_sha256)
            or any(not _digest(value) for value in resource_sha256)
            or type(scope_sha256) is not tuple or len(scope_sha256) != 2
            or any(not _digest(value) for value in scope_sha256)
            or not _digest(reference_sha256)
        ):
            raise WgsFilterError("filter state binding is invalid")
        self.path = path
        self.expected = expected
        self.source_manifest_sha256 = source_manifest_sha256
        self.policy_sha256 = policy_sha256
        self.resource_sha256 = resource_sha256
        self.scope_sha256 = scope_sha256
        self.reference_sha256 = reference_sha256
        self._owns_connection = connection is None
        descriptor = -1
        try:
            if connection is None:
                descriptor = os.open(
                    path, os.O_WRONLY | os.O_CREAT | os.O_EXCL | os.O_NOFOLLOW, 0o600,
                )
                os.fchmod(descriptor, 0o600)
                os.close(descriptor)
                descriptor = -1
                connection = sqlite3.connect(path)
                connection.enable_load_extension(False)
                connection.execute("PRAGMA foreign_keys=ON")
                connection.execute("PRAGMA cache_size=-8192")
                connection.execute("PRAGMA temp_store=FILE")
                connection.execute("PRAGMA mmap_size=0")
                connection.execute("PRAGMA trusted_schema=OFF")
            self._connection = connection
            self._connection.executescript("""
                CREATE TABLE decisions(
                    source_allele_id TEXT PRIMARY KEY,
                    record_number INTEGER NOT NULL,
                    allele_number INTEGER NOT NULL,
                    normalized_allele_id TEXT,
                    preparation_disposition TEXT NOT NULL,
                    early_filter_disposition TEXT,
                    bypass_reason TEXT,
                    scope_disposition TEXT NOT NULL,
                    population_request_keys INTEGER NOT NULL,
                    clinical_query_keys INTEGER NOT NULL,
                    UNIQUE(record_number,allele_number)
                );
                CREATE INDEX decision_normalized ON decisions(
                    normalized_allele_id,record_number,allele_number
                );
                CREATE TABLE normalized_aggregates(
                    normalized_allele_id TEXT PRIMARY KEY,
                    total_occurrences INTEGER NOT NULL,
                    excluded_occurrences INTEGER NOT NULL,
                    advancing_occurrences INTEGER NOT NULL,
                    first_record_number INTEGER NOT NULL,
                    first_allele_number INTEGER NOT NULL
                );
                CREATE TABLE reasons(
                    source_allele_id TEXT NOT NULL REFERENCES decisions(source_allele_id),
                    reason TEXT NOT NULL,
                    PRIMARY KEY(source_allele_id,reason)
                );
                CREATE TABLE filter_shards(
                    ordinal INTEGER PRIMARY KEY,
                    canonical_bytes BLOB NOT NULL
                );
                CREATE TABLE filter_meta(
                    singleton INTEGER PRIMARY KEY CHECK(singleton=1),
                    source_manifest_sha256 TEXT NOT NULL,
                    policy_sha256 TEXT NOT NULL,
                    population_sha256 TEXT NOT NULL,
                    clinical_sha256 TEXT NOT NULL,
                    exception_sha256 TEXT NOT NULL,
                    precision_sha256 TEXT NOT NULL,
                    expected_occurrences INTEGER NOT NULL,
                    committed_occurrences INTEGER NOT NULL,
                    sealed INTEGER NOT NULL CHECK(sealed IN (0,1)),
                    filter_manifest_sha256 TEXT,
                    working_manifest_sha256 TEXT,
                    working_summary_sha256 TEXT,
                    scope_policy_sha256 TEXT NOT NULL,
                    scope_manifest_sha256 TEXT NOT NULL,
                    reference_sha256 TEXT NOT NULL
                );
            """)
            self._connection.execute("PRAGMA user_version=2")
            self._connection.execute(
                "INSERT INTO filter_meta VALUES (1,?,?,?,?,?,?,?,0,0,NULL,NULL,NULL,?,?,?)",
                (source_manifest_sha256, policy_sha256, *resource_sha256, expected, *scope_sha256, reference_sha256),
            )
            self._connection.commit()
        except (OSError, sqlite3.Error):
            if descriptor >= 0:
                os.close(descriptor)
            raise WgsFilterError("cannot create filter replay state") from None
        self._count = 0
        self._previous = (0, 0)
        self._sealed = False
        self._closed = False

    @failure_boundary(
        "source_accounting", "filter_row_validate",
        known_errors=(WgsFilterError,),
    )
    def validate(self, decision: FilterDecision, expected_source, *, previous=None) -> None:
        prior = self._previous if previous is None else previous
        if (
            self._closed or self._sealed or not isinstance(decision, FilterDecision)
            or expected_source != (
                decision.record_number, decision.allele_number,
                decision.source_allele_id, decision.preparation_disposition,
                decision.normalized_allele_id,
            )
            or (decision.record_number, decision.allele_number) <= prior
            or decision.parent_source_manifest_sha256 != self.source_manifest_sha256
            or decision.policy_sha256 != self.policy_sha256
            or (decision.scope_policy_sha256, decision.scope_manifest_sha256) != self.scope_sha256
        ):
                raise _ClassifiedFilterError(
                    "filter decision source closure is invalid",
                    origin="source_accounting",
                    operation="filter_source_tuple_compare",
                )
        try:
            evidence = json.loads(decision.evidence_json)
            if (
                type(evidence) is not dict
                or set(evidence) != {
                    "schema", "version", "classification", "llm_policy",
                    "occurrence", "decision", "scope_decision", "resource_sha256",
                    "population_lookup", "clinical_lookup", "request_counts",
                }
                or evidence["schema"] != FILTER_EVIDENCE_SCHEMA
                or evidence["version"] != 2
                or evidence["classification"] != "controlled"
                or evidence["llm_policy"] != "no-llm"
                or evidence["occurrence"]["source_allele_id"]
                != decision.source_allele_id
                or tuple(evidence['resource_sha256']) != self.resource_sha256
            ):
                raise ValueError
            from genosapiens.track1_early_scope import EarlyScopeDecision
            scope_value = dict(evidence['scope_decision'])
            for key in ('reason_codes','matched_routes','source_span','normalized_span','matched_region','matched_exact_key'):
                if scope_value[key] is not None:
                    scope_value[key] = tuple(scope_value[key])
            scope = EarlyScopeDecision(**scope_value)
            if (scope.source_allele_id != decision.source_allele_id
                    or scope.status != decision.scope_disposition
                    or (scope.scope_policy_sha256,scope.scope_manifest_sha256) != self.scope_sha256
                    or scope.reference_manifest_sha256 != self.reference_sha256
                    or decision.preparation_disposition == 'normalized'
                    and scope.normalized_allele_id != decision.normalized_allele_id):
                raise ValueError
            counts = evidence['request_counts']
            if (type(counts) is not dict or set(counts) != {'population_request_keys','clinical_query_keys'}
                    or any(type(v) is not int or v not in (0,1) for v in counts.values())):
                raise ValueError
            for key, lookup in (('population_request_keys','population_lookup'),('clinical_query_keys','clinical_lookup')):
                if evidence[lookup] != ('performed' if counts[key] else 'not_performed'):
                    raise ValueError
            if scope.status in ('outside_scope','unresolved_scope'):
                if evidence['decision'] is not None or any(counts.values()):
                    raise ValueError
            else:
                pure = evidence['decision']
                if (pure['source_allele_id'] != decision.source_allele_id
                        or pure['policy_sha256'] != self.policy_sha256
                        or tuple(pure[key] for key in ('population_manifest_sha256','clinical_manifest_sha256',
                            'exception_manifest_sha256','precision_manifest_sha256')) != self.resource_sha256
                        or decision.preparation_disposition == 'normalized'
                        and pure['disposition'] != decision.early_filter_disposition):
                    raise ValueError
        except (KeyError, TypeError, ValueError, json.JSONDecodeError):
            raise _ClassifiedFilterError(
                "filter decision evidence is invalid",
                origin="source_accounting",
                operation="filter_row_validate",
            ) from None

    @failure_boundary(
        "source_accounting", "filter_row_validate",
        known_errors=(WgsFilterError,),
    )
    def add(self, decision: FilterDecision, expected_source) -> None:
        self.validate(decision, expected_source)
        try:
            self._connection.execute(
                "INSERT INTO decisions VALUES (?,?,?,?,?,?,?,?,?,?)",
                (
                    decision.source_allele_id, decision.record_number,
                    decision.allele_number, decision.normalized_allele_id,
                    decision.preparation_disposition,
                    decision.early_filter_disposition, decision.bypass_reason,
                    decision.scope_disposition,
                    json.loads(decision.evidence_json)['request_counts']['population_request_keys'],
                    json.loads(decision.evidence_json)['request_counts']['clinical_query_keys'],
                ),
            )
            reasons = decision.reason_codes or (decision.bypass_reason,)
            self._connection.executemany(
                "INSERT INTO reasons VALUES (?,?)",
                ((decision.source_allele_id, reason) for reason in reasons),
            )
        except sqlite3.Error:
            raise _ClassifiedFilterError(
                "filter decision evidence is invalid",
                origin="source_accounting",
                operation="filter_row_validate",
            ) from None
        self._count += 1
        self._previous = (decision.record_number, decision.allele_number)
        if self._owns_connection and self._count % 1024 == 0:
            self._connection.commit()

    def seal(self) -> None:
        if self._closed or self._sealed or self._count != self.expected:
            raise WgsFilterError("filter decision accounting is invalid")
        self._connection.execute("""
            INSERT INTO normalized_aggregates
            WITH ranked AS (
                SELECT normalized_allele_id,record_number,allele_number,
                       early_filter_disposition,scope_disposition,
                       row_number() OVER (
                           PARTITION BY normalized_allele_id
                           ORDER BY record_number,allele_number
                       ) AS source_ordinal
                FROM decisions WHERE preparation_disposition='normalized'
            )
            SELECT normalized_allele_id,count(*),
                   sum(CASE WHEN early_filter_disposition='exclude_common' THEN 1 ELSE 0 END),
                   sum(CASE WHEN scope_disposition='in_scope' AND early_filter_disposition IN ('retain','retain_unknown') THEN 1 ELSE 0 END),
                   max(CASE WHEN source_ordinal=1 THEN record_number END),
                   max(CASE WHEN source_ordinal=1 THEN allele_number END)
            FROM ranked GROUP BY normalized_allele_id
        """)
        self._connection.execute(
            "UPDATE filter_meta SET committed_occurrences=?,sealed=1 WHERE singleton=1",
            (self._count,),
        )
        self._connection.commit()
        self._sealed = True

    def record_shard(self, shard: Shard) -> None:
        if self._closed or self._sealed or type(shard) is not Shard:
            raise WgsFilterError("filter shard state is invalid")
        self._connection.execute(
            "INSERT INTO filter_shards VALUES (?,?)",
            (shard.ordinal, shard.as_bytes()),
        )
        self._connection.execute(
            "UPDATE filter_meta SET committed_occurrences=? WHERE singleton=1",
            (self._count,),
        )

    def bind_manifest(self, manifest_sha256: str) -> None:
        if self._closed or not self._sealed or not _digest(manifest_sha256):
            raise WgsFilterError("filter state manifest binding is invalid")
        updated = self._connection.execute(
            "UPDATE filter_meta SET filter_manifest_sha256=? WHERE singleton=1 "
            "AND filter_manifest_sha256 IS NULL", (manifest_sha256,),
        )
        if updated.rowcount != 1:
            raise WgsFilterError("filter state manifest binding is invalid")
        self._connection.commit()

    @classmethod
    def open_active(
        cls, path: Path, *, expected: int, source_manifest_sha256: str,
        policy_sha256: str, resource_sha256: tuple[str, str, str, str],
        scope_sha256: tuple[str, str],
        reference_sha256: str,
        stage_state: StageCommitState,
    ):
        """Restore only committed shard descriptors and the durable source cursor."""
        connection = stage_state.connection
        try:
            if connection.execute('SELECT scope_policy_sha256,scope_manifest_sha256,reference_sha256 FROM filter_meta WHERE singleton=1').fetchone() != (*scope_sha256,reference_sha256):
                raise ValueError
            meta = connection.execute(
                "SELECT source_manifest_sha256,policy_sha256,population_sha256,"
                "clinical_sha256,exception_sha256,precision_sha256,"
                "expected_occurrences,committed_occurrences,sealed,"
                "filter_manifest_sha256 FROM filter_meta WHERE singleton=1"
            ).fetchone()
            if (meta is None or meta[:7] != (
                    source_manifest_sha256, policy_sha256, *resource_sha256, expected)
                    or not _count(meta[7]) or meta[7] > expected
                    or meta[8] not in (0, 1)
                    or meta[8] == 1 and meta[7] != expected
                    or meta[8] == 0 and meta[9] is not None
                    or meta[8] == 1 and not _digest(meta[9])):
                raise ValueError
            rows = connection.execute(
                "SELECT ordinal,canonical_bytes FROM filter_shards ORDER BY ordinal"
            ).fetchall()
            partition_ids = [row[0] for row in connection.execute(
                "SELECT partition_id FROM commits ORDER BY ordinal")]
            shards = []
            total = 0
            previous_cursor = (0, 0)
            position = 0
            for partition_id in partition_ids:
                completed = stage_state.committed(partition_id)
                if completed is None or not completed.outputs:
                    raise ValueError
                group_rows = rows[position:position + len(completed.outputs)]
                group_total = 0
                for (ordinal, content), output in zip(group_rows, completed.outputs, strict=True):
                    shard = load_model(content, Shard)
                    if (ordinal != len(shards) + 1 or shard.ordinal != ordinal
                            or output != shard.artifact):
                        raise ValueError
                    object.__setattr__(shard, "_publication", getattr(output, "_publication"))
                    group_total += shard.rows
                    shards.append(shard)
                if (completed.counts.get("decisions") != group_total
                        or completed.cursor[2] != len(shards)
                        or completed.cursor[:2] <= previous_cursor
                        or partition_id != f"decision-shard-{position + 1:06d}"):
                    raise ValueError
                total += group_total
                previous_cursor = completed.cursor[:2]
                position += len(completed.outputs)
            if position != len(rows):
                raise ValueError
            if total != meta[7]:
                raise ValueError
            result = object.__new__(cls)
            result.path = path
            result.expected = expected
            result.source_manifest_sha256 = source_manifest_sha256
            result.policy_sha256 = policy_sha256
            result.resource_sha256 = resource_sha256
            result.scope_sha256 = scope_sha256
            result.reference_sha256 = reference_sha256
            result._connection = connection
            result._owns_connection = False
            result._count = total
            result._previous = previous_cursor
            result._sealed = bool(meta[8])
            result._closed = False
            return result, shards, meta[9]
        except (IndexError, TypeError, ValueError, sqlite3.Error):
            raise WgsFilterError("incomplete filter state is unavailable or stale") from None

    @classmethod
    def open_sealed(
        cls, path: Path, *, expected: int, source_manifest_sha256: str,
        policy_sha256: str, resource_sha256: tuple[str, str, str, str],
        scope_sha256: tuple[str, str],
        reference_sha256: str,
        filter_manifest_sha256: str,
    ):
        """Reuse private completed state without replaying canonical decision shards."""
        connection = None
        try:
            metadata = path.lstat()
            if (not path.is_absolute() or path.is_symlink()
                    or not stat.S_ISREG(metadata.st_mode)
                    or stat.S_IMODE(metadata.st_mode) != 0o600
                    or metadata.st_uid != os.getuid()
                    or not _digest(filter_manifest_sha256)):
                raise ValueError
            connection = sqlite3.connect(path.as_uri() + "?mode=ro&immutable=1", uri=True)
            connection.execute("PRAGMA query_only=ON")
            connection.execute("PRAGMA trusted_schema=OFF")
            if (connection.execute("PRAGMA user_version").fetchone() != (2,)
                    or connection.execute('SELECT scope_policy_sha256,scope_manifest_sha256,reference_sha256 FROM filter_meta WHERE singleton=1').fetchone() != (*scope_sha256,reference_sha256)
                    or {row[0] for row in connection.execute(
                        "SELECT name FROM sqlite_master WHERE type='table' "
                        "AND name IN ('decisions','reasons','filter_shards',"
                        "'filter_meta','normalized_aggregates')"
                    )} != {
                        "decisions", "reasons", "filter_shards", "filter_meta",
                        "normalized_aggregates",
                    }
                    or connection.execute(
                        "SELECT singleton,source_manifest_sha256,policy_sha256,"
                        "population_sha256,clinical_sha256,exception_sha256,precision_sha256,"
                        "expected_occurrences,committed_occurrences,sealed,"
                        "filter_manifest_sha256 FROM filter_meta WHERE singleton=1"
                    ).fetchone()
                    != (1, source_manifest_sha256, policy_sha256, *resource_sha256,
                        expected, expected, 1, filter_manifest_sha256)):
                raise ValueError
            result = object.__new__(cls)
            result.path = path
            result.expected = expected
            result.source_manifest_sha256 = source_manifest_sha256
            result.policy_sha256 = policy_sha256
            result.resource_sha256 = resource_sha256
            result.scope_sha256 = scope_sha256
            result.reference_sha256 = reference_sha256
            result._connection = connection
            result._owns_connection = True
            result._count = expected
            result._previous = (0, 0)
            result._sealed = True
            result._closed = False
            return result
        except (OSError, sqlite3.Error, ValueError):
            if connection is not None:
                connection.close()
            raise WgsFilterError("completed filter state is unavailable or stale") from None

    def working_binding(self) -> tuple[str, str]:
        if self._closed or not self._sealed:
            raise WgsFilterError("completed filter state is unavailable")
        row = self._connection.execute(
            "SELECT working_manifest_sha256,working_summary_sha256 "
            "FROM filter_meta WHERE singleton=1"
        ).fetchone()
        if row is None or not all(_digest(value) for value in row):
            raise WgsFilterError("completed working state is unavailable")
        return row

    @classmethod
    def bind_working(cls, path: Path, *, filter_manifest_sha256: str,
                     working_manifest_sha256: str, working_summary_sha256: str) -> None:
        if any(not _digest(value) for value in (
                filter_manifest_sha256, working_manifest_sha256,
                working_summary_sha256)):
            raise WgsFilterError("working state binding is invalid")
        connection = None
        try:
            metadata = path.lstat()
            if (not path.is_absolute() or path.is_symlink()
                    or not stat.S_ISREG(metadata.st_mode)
                    or stat.S_IMODE(metadata.st_mode) != 0o600
                    or metadata.st_uid != os.getuid()):
                raise ValueError
            connection = sqlite3.connect(path)
            connection.execute("PRAGMA journal_mode=DELETE")
            connection.execute("PRAGMA synchronous=FULL")
            updated = connection.execute(
                "UPDATE filter_meta SET working_manifest_sha256=?,"
                "working_summary_sha256=? WHERE singleton=1 AND sealed=1 "
                "AND filter_manifest_sha256=? AND working_manifest_sha256 IS NULL",
                (working_manifest_sha256, working_summary_sha256,
                 filter_manifest_sha256),
            )
            if updated.rowcount != 1:
                raise ValueError
            connection.commit()
        except (OSError, sqlite3.Error, ValueError):
            raise WgsFilterError("working state binding is unavailable") from None
        finally:
            if connection is not None:
                connection.close()

    def summary(self) -> FilterSummary:
        if self._closed or not self._sealed:
            raise WgsFilterError("filter replay state is unavailable")
        preparation = dict(_rows(self._connection.execute(
            "SELECT preparation_disposition,count(*) FROM decisions "
            "GROUP BY preparation_disposition"
        )))
        early = dict(_rows(self._connection.execute(
            "SELECT early_filter_disposition,count(*) FROM decisions "
            "WHERE early_filter_disposition IS NOT NULL GROUP BY early_filter_disposition"
        )))
        reasons = tuple(
            NamedCount(reason, count)
            for reason, count in _rows(self._connection.execute(
                "SELECT reason,count(*) FROM reasons GROUP BY reason ORDER BY reason"
            ))
        )
        unique_normalized = self._connection.execute(
            "SELECT count(DISTINCT normalized_allele_id) FROM decisions "
            "WHERE normalized_allele_id IS NOT NULL"
        ).fetchone()[0]
        scope = dict(self._connection.execute('SELECT scope_disposition,count(*) FROM decisions GROUP BY scope_disposition'))
        counts = self._connection.execute('SELECT coalesce(sum(population_request_keys),0),coalesce(sum(clinical_query_keys),0) FROM decisions').fetchone()
        return FilterSummary(
            self._count,
            unique_normalized,
            tuple(NamedCount(name, preparation.get(name, 0))
                  for name in PREPARATION_DISPOSITIONS),
            tuple(NamedCount(name, early.get(name, 0))
                  for name in EARLY_FILTER_DISPOSITIONS),
            reasons,
            scope_counts=tuple(NamedCount(name,scope.get(name,0)) for name in SCOPE_DISPOSITIONS),
            in_scope_normalized_occurrences=sum(early.values()),
            population_request_keys=counts[0], clinical_query_keys=counts[1],
        )

    def expected_working_summary(self) -> WorkingSetSummary:
        if self._closed or not self._sealed:
            raise WgsFilterError("filter replay state is unavailable")
        retained_normalized, retained_identities = self._connection.execute("""
            SELECT coalesce(sum(CASE WHEN advancing_occurrences>0
                                     THEN total_occurrences ELSE 0 END),0),
                   coalesce(sum(advancing_occurrences>0),0)
            FROM normalized_aggregates
        """).fetchone()
        native = self._connection.execute(
            "SELECT count(*) FROM decisions WHERE preparation_disposition='native_class' AND scope_disposition='in_scope'"
        ).fetchone()[0]
        quarantined = self._connection.execute(
            "SELECT count(*) FROM decisions WHERE preparation_disposition='quarantined'"
        ).fetchone()[0]
        non_variant = self._connection.execute(
            "SELECT count(*) FROM decisions WHERE preparation_disposition='non_variant'"
        ).fetchone()[0]
        bypassed = quarantined + non_variant
        identities = retained_identities + native
        excluded = outside = unresolved = 0
        for scope, early, count in self._connection.execute('''
            SELECT d.scope_disposition,d.early_filter_disposition,count(*) FROM decisions d
            LEFT JOIN normalized_aggregates a USING(normalized_allele_id)
            WHERE (d.preparation_disposition='normalized' AND a.advancing_occurrences=0)
               OR (d.preparation_disposition='native_class' AND d.scope_disposition!='in_scope')
            GROUP BY d.scope_disposition,d.early_filter_disposition'''):
            if scope == 'outside_scope': outside += count
            elif scope == 'unresolved_scope': unresolved += count
            elif scope == 'in_scope' and early == 'exclude_common': excluded += count
            else: raise WgsFilterError('working scope accounting is invalid')
        return WorkingSetSummary(
            self._count, retained_normalized + native, identities, excluded, bypassed,
            tuple(NamedCount(name, count) for name, count in (
                ("pending", retained_normalized + native),
                ("excluded_common_by_policy", excluded),
                ("quarantined", quarantined),
                ("non_variant", non_variant),
                ('excluded_outside_scope_by_policy', outside),
                ('deferred_unresolved_scope', unresolved),
            )),
            excluded_scope_occurrences=outside, unresolved_scope_occurrences=unresolved,
        )

    def iter_expected_working_rows(self, parent_filter_sha256: str):
        if self._closed or not self._sealed or not _digest(parent_filter_sha256):
            raise WgsFilterError("filter replay state is unavailable")
        candidates = self._connection.execute("""
            SELECT normalized_allele_id AS identity_id,'sequence' AS route,
                   first_record_number AS first_record,
                   first_allele_number AS first_allele
            FROM normalized_aggregates
            WHERE advancing_occurrences>0
            UNION ALL
            SELECT source_allele_id,'native',record_number,allele_number
            FROM decisions WHERE preparation_disposition='native_class' AND scope_disposition='in_scope'
            ORDER BY first_record,first_allele,identity_id
        """)
        ordinal = 0
        for identity_id, route, first_record, first_allele in _rows(candidates):
            ordinal += 1
            if route == "sequence":
                source_ids = tuple(row[0] for row in _rows(self._connection.execute("""
                    SELECT source_allele_id FROM decisions
                    WHERE normalized_allele_id=? ORDER BY record_number,allele_number
                """, (identity_id,))))
            else:
                source_ids = (identity_id,)
            if len(source_ids) > 65536:
                raise WgsFilterError("working identity occurrence capacity exceeded")
            yield WorkingSetRow(
                ordinal,
                identity_id,
                route,
                source_ids,
                parent_filter_sha256,
                canonical({
                    "schema": "genosapiens.wgs-working-set-payload",
                    "version": 1,
                    "classification": "controlled",
                    "llm_policy": "no-llm",
                    "first_record_number": first_record,
                    "first_allele_number": first_allele,
                    "source_occurrences": len(source_ids),
                    "closure": "any_non_excluded_occurrence_retains_identity",
                    "analysis_outcome": "pending",
                }).decode("ascii"),
            )

    def close(self) -> None:
        if not self._closed:
            if self._owns_connection:
                self._connection.close()
            self._closed = True

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        incomplete = exc_type is None and not self._sealed
        self.close()
        if incomplete:
            raise WgsFilterError("filter replay state was not sealed")
        return False


def _validate_source_accounting(
    manifest: LedgerManifest, closure: SourceAccountingClosure,
) -> None:
    if (
        not isinstance(manifest, LedgerManifest)
        or not isinstance(closure, SourceAccountingClosure)
        or closure.source_manifest_sha256 != manifest.sha256
        or closure.input_sha256 != manifest.source.input_sha256
        or closure.total_records != manifest.summary.total_records
        or closure.total_physical_alternate_slots
        != manifest.summary.total_physical_alternate_slots
        or closure.total_biological_alternates
        != manifest.summary.total_biological_alternates
    ):
        raise WgsFilterError("complete source accounting is required")


def _filter_parents(
    source_manifest: LedgerManifest, *, policy_sha256: str,
    population_manifest_sha256: str, clinical_manifest_sha256: str,
    exception_manifest_sha256: str, precision_manifest_sha256: str,
    scope_policy_sha256: str, scope_manifest_sha256: str,
    completed_source_input: ParentIdentity | None = None,
) -> tuple[ParentIdentity, ...]:
    if completed_source_input is not None and (
            completed_source_input.role != "completed_source_input"
            or completed_source_input.parent_schema != "genosapiens.wgs-completed-source-input"
            or completed_source_input.parent_version != 1):
        raise WgsFilterError("filter source input parent is invalid")
    return (
        ParentIdentity(
            "source_manifest", SOURCE_MANIFEST_SCHEMA, 1, source_manifest.sha256,
        ),
        ParentIdentity(
            "prefilter_policy", "genosapiens.common-allele-prefilter-policy", 1,
            policy_sha256,
        ),
        ParentIdentity(
            "population_manifest", "genosapiens.population-manifest", 1,
            population_manifest_sha256,
        ),
        ParentIdentity(
            "clinical_manifest", "genosapiens.clinical-manifest", 1,
            clinical_manifest_sha256,
        ),
        ParentIdentity(
            "exception_manifest", "genosapiens.exception-manifest", 1,
            exception_manifest_sha256,
        ),
        ParentIdentity(
            "precision_manifest", "genosapiens.precision-manifest", 1,
            precision_manifest_sha256,
        ),
        ParentIdentity('early_scope_policy', 'genosapiens.track1-early-scope-policy', 1, scope_policy_sha256),
        ParentIdentity('early_scope_manifest', 'genosapiens.track1-early-scope-manifest', 1, scope_manifest_sha256),
    ) + (() if completed_source_input is None else (completed_source_input,))


def _verified_filter_resources(
    source_manifest: LedgerManifest,
    filter_manifest: FilterManifest,
    *,
    policy_sha256: str,
    population_manifest_sha256: str,
    clinical_manifest_sha256: str,
    exception_manifest_sha256: str,
    precision_manifest_sha256: str,
    scope_policy_sha256: str, scope_manifest_sha256: str,
) -> tuple[str, str, str, str]:
    resources = (
        population_manifest_sha256, clinical_manifest_sha256,
        exception_manifest_sha256, precision_manifest_sha256,
    )
    if (
        not _digest(policy_sha256)
        or any(not _digest(value) for value in resources)
        or filter_manifest.parent_source_manifest_sha256 != source_manifest.sha256
        or filter_manifest.policy_sha256 != policy_sha256
        or filter_manifest.scope_policy_sha256 != scope_policy_sha256
        or filter_manifest.scope_manifest_sha256 != scope_manifest_sha256
        or filter_manifest.parents != _filter_parents(
            source_manifest,
            policy_sha256=policy_sha256,
            scope_policy_sha256=filter_manifest.scope_policy_sha256,
            scope_manifest_sha256=filter_manifest.scope_manifest_sha256,
            population_manifest_sha256=population_manifest_sha256,
            clinical_manifest_sha256=clinical_manifest_sha256,
            exception_manifest_sha256=exception_manifest_sha256,
            precision_manifest_sha256=precision_manifest_sha256,
            completed_source_input=next((parent for parent in filter_manifest.parents
                if parent.role == "completed_source_input"), None),
        )
    ):
        raise WgsFilterError("filter manifest parent closure is invalid")
    return resources


def _publish_filter_manifest(root: Path, manifest: FilterManifest):
    require_published_shards(root, manifest.shards)
    summary_content = manifest.summary.as_bytes()
    summary_artifact = publish_bytes(root, FILTER_SUMMARY_ARTIFACT, summary_content)
    manifest_content = manifest.as_bytes()
    manifest_artifact = publish_bytes(root, MANIFEST_ARTIFACT, manifest_content)
    return manifest_artifact, summary_artifact


@failure_boundary(
    "filter_publication", "filter_artifact_publish",
    known_errors=(WgsFilterError,),
)
def publish_filter_decisions(
    *,
    source_root: str | Path,
    source_manifest: LedgerManifest,
    source_accounting: SourceAccountingClosure,
    source_routing_context: SourceRoutingContext,
    source_cursor: SourceDecisionCursor | None = None,
    decisions: Iterable[FilterDecision],
    output_root: str | Path,
    policy_sha256: str,
    population_manifest_sha256: str,
    clinical_manifest_sha256: str,
    exception_manifest_sha256: str,
    precision_manifest_sha256: str,
    scope_policy_sha256: str,
    scope_manifest_sha256: str,
    retained_state_path: str | Path | None = None,
    stage_binding: StageRunBinding | None = None,
    resume: bool = False,
    max_attempts: int = 1,
    technical_classification_resume: bool = False,
    max_shard_rows: int = MAX_SHARD_ROWS,
    max_shard_bytes: int = MAX_SHARD_BYTES,
    phase_recorder: RuntimePhaseRecorder | None = None,
) -> FilterPublication:
    """Publish one exact decision for every physical source occurrence."""

    scope_sha256 = (scope_policy_sha256, scope_manifest_sha256)

    _validate_source_accounting(source_manifest, source_accounting)
    _validate_source_routing(source_manifest, source_routing_context)
    source_root = Path(source_root)
    root = Path(output_root)
    resources = (
        population_manifest_sha256, clinical_manifest_sha256,
        exception_manifest_sha256, precision_manifest_sha256,
    )
    if (
        not source_root.is_absolute() or not root.is_absolute()
        or not _digest(policy_sha256)
        or any(not _digest(value) for value in resources)
        or not _count(max_shard_rows, 1) or max_shard_rows > MAX_SHARD_ROWS
        or not _count(max_shard_bytes, 1) or max_shard_bytes > MAX_SHARD_BYTES
        or type(resume) is not bool or not _count(max_attempts, 1)
        or type(technical_classification_resume) is not bool
        or technical_classification_resume and (not resume or stage_binding is None)
        or phase_recorder is not None
        and type(phase_recorder) is not RuntimePhaseRecorder
    ):
        raise WgsFilterError("filter publication binding is invalid")
    state_path = (Path(retained_state_path) if retained_state_path is not None
                  else root / "filter-state.sqlite")
    if retained_state_path is not None and (
        not state_path.is_absolute() or state_path.name != "filter.sqlite"
        or state_path.parent.name != "state"
    ):
        raise WgsFilterError("retained filter state path is invalid")
    if stage_binding is not None and (
        type(stage_binding) is not StageRunBinding
        or stage_binding.stage != "prefilter"
        or stage_binding.policy_sha256 != policy_sha256
        or stage_binding.input_sha256 != source_manifest.source.input_sha256
        or stage_binding.reference_sha256 != source_manifest.source.reference_manifest_sha256
        or retained_state_path is None
    ):
        raise WgsFilterError("filter stage binding is invalid")
    source_bound = type(decisions) is SourceBoundDecisionStream
    if source_bound and (
        decisions.source_cursor is not source_cursor
        or source_cursor.source_root != source_root
        or source_cursor.manifest != source_manifest
        or source_cursor.routing != source_routing_context
    ):
        raise WgsFilterError("filter decision cursor binding is invalid")
    shards: list[Shard] = []
    batch: list[FilterDecision] = []
    expected_batch = []
    batch_bytes = 0
    previous = None
    first_row = 1
    decision_iterator = None
    primary_failure = None
    state_exists = state_path.exists() or state_path.is_symlink()
    if state_exists and (not resume or stage_binding is None):
        raise WgsFilterError("filter state already exists without exact resume binding")
    try:
        with stage_lock(root, "prefilter") as lock, (
            StageCommitState(root, state_path, stage_binding, lock=lock,
                             create=not state_exists)
            if stage_binding is not None else nullcontext(None)
        ) as stage_state:
            bound_manifest = None
            if state_exists:
                state, restored_shards, bound_manifest = _FilterState.open_active(
                    state_path,
                    expected=source_accounting.total_physical_alternate_slots,
                    source_manifest_sha256=source_manifest.sha256,
                    policy_sha256=policy_sha256,
                    resource_sha256=resources,
                    scope_sha256=scope_sha256,
                    reference_sha256=source_manifest.source.reference_manifest_sha256,
                    stage_state=stage_state,
                )
                shards.extend(restored_shards)
                previous = shards[-1] if shards else None
                first_row = state._count + 1
            else:
                state = _FilterState(
                    state_path,
                    expected=source_accounting.total_physical_alternate_slots,
                    source_manifest_sha256=source_manifest.sha256,
                    policy_sha256=policy_sha256,
                    resource_sha256=resources,
                    scope_sha256=scope_sha256,
                    reference_sha256=source_manifest.source.reference_manifest_sha256,
                    connection=stage_state.connection if stage_state is not None else None,
                )
            with state:
                partition_ids = ([] if stage_state is None else [row[0] for row in
                    stage_state.connection.execute("SELECT partition_id FROM commits ORDER BY ordinal")])
                current_start = None
                if stage_state is not None and not state._sealed:
                    next_partition = f"decision-shard-{len(shards) + 1:06d}"
                    incomplete = stage_state.open_attempt(next_partition)
                    if incomplete is not None:
                        stage_state.fail(
                            incomplete, failure_code="prefilter_interrupted",
                            interrupted=True,
                        )
                        stage_state.archive_orphans(
                            names=tuple(f"decisions-{ordinal:06d}.jsonl.gz"
                                for ordinal in range(len(shards) + 1, len(shards) + 65)
                                if (root / f"decisions-{ordinal:06d}.jsonl.gz").exists()),
                            attempt_id=incomplete.attempt_id,
                        )

                def flush_batch():
                    nonlocal previous, first_row, batch, expected_batch, batch_bytes, current_start
                    if not batch:
                        return
                    published = []
                    chunk = []
                    size = 0
                    def publish_chunk():
                        # The journal accepts at most 64 outputs per atomic batch.
                        # Stop before an orphan can exceed the recovery window.
                        if stage_state is not None and len(published) >= 64:
                            raise WgsFilterError("completed filter batch exceeds journal output bound")
                        return publish_filter_shard(root, tuple(chunk),
                            ordinal=len(shards) + len(published) + 1,
                            first_record=first_row + sum(row.rows for row in published),
                            previous=published[-1] if published else previous)
                    for decision in batch:
                        encoded_size = len(decision.as_bytes())
                        if chunk and (len(chunk) >= max_shard_rows or size + encoded_size > max_shard_bytes):
                            published.append(publish_chunk())
                            chunk, size = [], 0
                        chunk.append(decision)
                        size += encoded_size
                    if chunk:
                        published.append(publish_chunk())
                    if stage_state is not None:
                        def commit_rows(_connection):
                            for decision, expected in zip(batch, expected_batch, strict=True):
                                state.add(decision, expected)
                            for shard in published:
                                state.record_shard(shard)

                        counts = {"decisions": len(batch)}
                        for key in ("population_request_keys", "clinical_query_keys"):
                            counts[key] = sum(json.loads(row.evidence_json)["request_counts"][key]
                                              for row in batch)
                        with runtime_phase(phase_recorder, "filter_index_commit"):
                            stage_state.commit_partition(current_start,
                                cursor=(batch[-1].record_number, batch[-1].allele_number,
                                        published[-1].ordinal),
                                outputs=tuple(shard.artifact for shard in published),
                                counts=counts, reducer=commit_rows)
                        partition_ids.append(current_start.partition_id)
                    previous = published[-1]
                    shards.extend(published)
                    first_row += len(batch)
                    batch, expected_batch, batch_bytes = [], [], 0
                    if stage_state is not None:
                        current_start = None

                if state._sealed:
                    if bound_manifest is None:
                        raise WgsFilterError("sealed filter state lacks its manifest binding")
                    with (root / MANIFEST_ARTIFACT).open("rb") as stream:
                        manifest_content = stream.read(MAX_MANIFEST_BYTES + 1)
                    manifest = load_model(manifest_content, FilterManifest)
                    if (manifest.sha256 != bound_manifest
                            or manifest.shards != tuple(shards)
                            or manifest.parent_source_manifest_sha256 != source_manifest.sha256
                            or manifest.policy_sha256 != policy_sha256):
                        raise WgsFilterError("sealed filter manifest changed")
                    _verify_filter_artifacts(root, manifest)
                    manifest_artifact = attest_existing_artifact(
                        root, _artifact_for(MANIFEST_ARTIFACT, manifest_content),
                    )
                    summary_artifact = _artifact_for(
                        FILTER_SUMMARY_ARTIFACT, manifest.summary.as_bytes(),
                    )
                    latest = stage_state.connection.execute(
                        "SELECT event_type FROM events ORDER BY ordinal DESC LIMIT 1"
                    ).fetchone()
                    if latest != ("stage_completed",):
                        stage_state.finish(
                            expected_partitions=tuple(partition_ids),
                            output=manifest_artifact,
                            process_id="process-" + os.urandom(16).hex(),
                            closure_check=lambda connection: connection.execute(
                                "SELECT committed_occurrences,sealed,filter_manifest_sha256 "
                                "FROM filter_meta WHERE singleton=1"
                            ).fetchone() == (state.expected, 1, manifest.sha256),
                        )
                    return FilterPublication(
                        manifest, manifest_artifact, summary_artifact,
                    )

                if state_exists and state._previous != (0, 0):
                    if not source_bound:
                        raise WgsFilterError("filter resume requires source-owned decision cursor")
                    source_cursor.resume_after(*state._previous)
                completed_batches = source_bound and decisions.completed_batches
                decision_iterator = iter(decisions.decisions if completed_batches else decisions)

                def matched_decisions():
                    if source_bound:
                        for item in decision_iterator:
                            if completed_batches:
                                if type(item) is not CompletedFilterBatch:
                                    raise WgsFilterError("completed batch boundary is missing")
                                for index, decision in enumerate(item.decisions):
                                    yield source_cursor.pop_expected(), decision, index == len(item.decisions) - 1
                            else:
                                yield source_cursor.pop_expected(), item, False
                        source_cursor.require_complete()
                    else:
                        for expected in _replay_source_tuples(
                            source_root, source_manifest, source_routing_context,
                        ):
                            try:
                                decision = next(decision_iterator)
                            except StopIteration:
                                raise WgsFilterError(
                                    "filter decision coverage is incomplete"
                                ) from None
                            yield expected, decision, False
                        try:
                            next(decision_iterator)
                        except StopIteration:
                            pass
                        else:
                            raise WgsFilterError("filter decision coverage has extra rows")

                for expected_source, decision, batch_end in matched_decisions():
                    if stage_state is None:
                        state.add(decision, expected_source)
                    else:
                        prior = (
                            (batch[-1].record_number, batch[-1].allele_number)
                            if batch else state._previous
                        )
                        state.validate(decision, expected_source, previous=prior)
                    size = len(decision.as_bytes())
                    if size > max_shard_bytes:
                        raise WgsFilterError("filter decision exceeds its bounded size")
                    if not completed_batches and batch and (
                        len(batch) >= max_shard_rows
                        or batch_bytes + size > max_shard_bytes
                    ):
                        flush_batch()
                    if stage_state is not None and current_start is None:
                        partition_id = f"decision-shard-{len(shards) + 1:06d}"
                        partition_attempts = max_attempts
                        if technical_classification_resume:
                            starts = stage_state.connection.execute(
                                "SELECT count(*) FROM events WHERE partition_id=? "
                                "AND event_type='attempt_started'", (partition_id,),
                            ).fetchone()[0]
                            partition_attempts = max(max_attempts, starts + 1)
                        current_start = stage_state.begin(
                            partition_id=partition_id,
                            scientific=False, max_attempts=partition_attempts,
                        )
                    batch.append(decision)
                    expected_batch.append(expected_source)
                    batch_bytes += size
                    if batch_end:
                        flush_batch()
                close = getattr(decision_iterator, "close", None)
                decision_iterator = None
                if close is not None:
                    try:
                        close()
                    except (Exception, KeyboardInterrupt) as error:
                        raise static_failure(
                            error,
                            "resource_session",
                            operation="resource_session_close",
                        ) from None
                if batch:
                    flush_batch()
                state.seal()
                summary = state.summary()
                parents = _filter_parents(
                    source_manifest,
                    policy_sha256=policy_sha256,
                    scope_policy_sha256=scope_policy_sha256,
                    scope_manifest_sha256=scope_manifest_sha256,
                    population_manifest_sha256=population_manifest_sha256,
                    clinical_manifest_sha256=clinical_manifest_sha256,
                    exception_manifest_sha256=exception_manifest_sha256,
                    precision_manifest_sha256=precision_manifest_sha256,
                    completed_source_input=None if stage_binding is None else next(
                        parent for parent in stage_binding.parents if parent.role == "completed_source_input"),
                )
                manifest = FilterManifest(
                    source_manifest.sha256, policy_sha256, parents, summary, tuple(shards),
                    scope_policy_sha256=scope_policy_sha256,
                    scope_manifest_sha256=scope_manifest_sha256,
                )
                manifest_artifact, summary_artifact = _publish_filter_manifest(root, manifest)
                if retained_state_path is not None:
                    state.bind_manifest(manifest.sha256)
                if stage_state is not None:
                    stage_state.finish(
                        expected_partitions=tuple(partition_ids),
                        output=manifest_artifact,
                        process_id="process-" + os.urandom(16).hex(),
                        closure_check=lambda connection: connection.execute(
                            "SELECT committed_occurrences,sealed,filter_manifest_sha256 "
                            "FROM filter_meta WHERE singleton=1"
                        ).fetchone() == (
                            source_accounting.total_physical_alternate_slots,
                            1, manifest.sha256,
                        ),
                    )
                return FilterPublication(manifest, manifest_artifact, summary_artifact)
    except (Exception, KeyboardInterrupt) as error:
        primary_failure = static_failure(
            error,
            "filter_publication",
            operation="filter_artifact_publish",
        )
        raise primary_failure from None
    finally:
        close = getattr(decision_iterator, "close", None)
        if close is not None:
            try:
                close()
            except (Exception, KeyboardInterrupt) as error:
                cleanup = static_failure(
                    error,
                    "resource_session",
                    operation="resource_session_close",
                )
                if primary_failure is not None:
                    primary_failure.add_secondary(cleanup)
                else:
                    raise cleanup from None
        if retained_state_path is None:
            try:
                state_path.unlink()
            except (FileNotFoundError, OSError):
                pass


def _iter_filter_rows(root: Path, manifest: FilterManifest):
    for shard in manifest.shards:
        yield from iter_filter_shard(root, shard)


def _verify_filter_artifacts(root: Path, manifest: FilterManifest) -> None:
    manifest_content = manifest.as_bytes()
    if (
        read_artifact(root, _artifact_for(MANIFEST_ARTIFACT, manifest_content))
        != manifest_content
        or read_artifact(
            root,
            _artifact_for(FILTER_SUMMARY_ARTIFACT, manifest.summary.as_bytes()),
        ) != manifest.summary.as_bytes()
    ):
        raise WgsFilterError("filter manifest artifacts are invalid")


def _check_filter_read_boundary(
    *, source_manifest, source_accounting, source_routing_context, filter_root,
    filter_manifest, policy_sha256, population_manifest_sha256,
    clinical_manifest_sha256, exception_manifest_sha256, precision_manifest_sha256,
    scope_policy_sha256: str, scope_manifest_sha256: str,
):
    _validate_source_accounting(source_manifest, source_accounting)
    _validate_source_routing(source_manifest, source_routing_context)
    root = Path(filter_root)
    if not root.is_absolute():
        raise WgsFilterError("filter verification binding is invalid")
    resources = _verified_filter_resources(
        source_manifest, filter_manifest,
        policy_sha256=policy_sha256,
        population_manifest_sha256=population_manifest_sha256,
        clinical_manifest_sha256=clinical_manifest_sha256,
        exception_manifest_sha256=exception_manifest_sha256,
        precision_manifest_sha256=precision_manifest_sha256,
        scope_policy_sha256=scope_policy_sha256, scope_manifest_sha256=scope_manifest_sha256,
    )
    _verify_filter_artifacts(root, filter_manifest)
    return resources


@failure_boundary(
    "source_accounting", "filter_row_validate",
    known_errors=(WgsFilterError,),
)
def verify_filter_decisions(
    *,
    source_root: str | Path,
    source_manifest: LedgerManifest,
    source_accounting: SourceAccountingClosure,
    source_routing_context: SourceRoutingContext,
    filter_root: str | Path,
    filter_manifest: FilterManifest,
    policy_sha256: str,
    population_manifest_sha256: str,
    clinical_manifest_sha256: str,
    exception_manifest_sha256: str,
    precision_manifest_sha256: str,
    scope_policy_sha256: str, scope_manifest_sha256: str,
) -> FilterSummary:
    """Replay every source occurrence and recompute the embedded filter summary."""
    scope_sha256 = (filter_manifest.scope_policy_sha256, filter_manifest.scope_manifest_sha256)

    _validate_source_accounting(source_manifest, source_accounting)
    _validate_source_routing(source_manifest, source_routing_context)
    source_root, root = Path(source_root), Path(filter_root)
    if not source_root.is_absolute() or not root.is_absolute():
        raise WgsFilterError("filter verification binding is invalid")
    resources = _verified_filter_resources(
        source_manifest,
        filter_manifest,
        policy_sha256=policy_sha256,
        population_manifest_sha256=population_manifest_sha256,
        clinical_manifest_sha256=clinical_manifest_sha256,
        exception_manifest_sha256=exception_manifest_sha256,
        precision_manifest_sha256=precision_manifest_sha256,
        scope_policy_sha256=scope_policy_sha256, scope_manifest_sha256=scope_manifest_sha256,
    )
    _verify_filter_artifacts(root, filter_manifest)
    state_path = root / "verify-filter-state.sqlite"
    try:
        with stage_lock(root, "prefilter-verify"):
            with _FilterState(
                state_path,
                expected=source_accounting.total_physical_alternate_slots,
                source_manifest_sha256=source_manifest.sha256,
                policy_sha256=filter_manifest.policy_sha256,
                resource_sha256=resources,
                    scope_sha256=scope_sha256,
                    reference_sha256=source_manifest.source.reference_manifest_sha256,
            ) as state:
                rows = _iter_filter_rows(root, filter_manifest)
                for expected_source in _replay_source_tuples(
                    source_root, source_manifest, source_routing_context,
                ):
                    try:
                        state.add(next(rows), expected_source)
                    except StopIteration:
                        raise WgsFilterError("filter decision coverage is incomplete") from None
                try:
                    next(rows)
                except StopIteration:
                    pass
                else:
                    raise WgsFilterError("filter decision coverage has extra rows")
                state.seal()
                summary = state.summary()
                if summary != filter_manifest.summary:
                    raise WgsFilterError("filter summary replay is invalid")
                return summary
    finally:
        try:
            state_path.unlink()
        except (FileNotFoundError, OSError):
            pass


def _working_parents(
    source_manifest: LedgerManifest, filter_manifest: FilterManifest,
) -> tuple[ParentIdentity, ...]:
    return (
        ParentIdentity(
            "source_manifest", SOURCE_MANIFEST_SCHEMA, 1, source_manifest.sha256,
        ),
        ParentIdentity(
            "filter_manifest", FILTER_MANIFEST_SCHEMA, 2, filter_manifest.sha256,
        ),
    )


def _replay_filter_rows(
    state: _FilterState,
    *,
    source_root: Path,
    source_manifest: LedgerManifest,
    source_routing_context: SourceRoutingContext,
    filter_root: Path,
    filter_manifest: FilterManifest,
) -> None:
    rows = _iter_filter_rows(filter_root, filter_manifest)
    for expected_source in _replay_source_tuples(
        source_root, source_manifest, source_routing_context,
    ):
        try:
            state.add(next(rows), expected_source)
        except StopIteration:
            raise WgsFilterError("filter decision coverage is incomplete") from None
    try:
        next(rows)
    except StopIteration:
        return
    raise WgsFilterError("filter decision coverage has extra rows")


@failure_boundary(
    "working_set_publication", "working_set_publish",
    known_errors=(WgsFilterError,),
)
def publish_working_set(
    *,
    source_root: str | Path,
    source_manifest: LedgerManifest,
    source_accounting: SourceAccountingClosure,
    source_routing_context: SourceRoutingContext,
    filter_root: str | Path,
    filter_manifest: FilterManifest,
    output_root: str | Path,
    policy_sha256: str,
    population_manifest_sha256: str,
    clinical_manifest_sha256: str,
    exception_manifest_sha256: str,
    precision_manifest_sha256: str,
    scope_policy_sha256: str, scope_manifest_sha256: str,
    retained_state_path: str | Path | None = None,
    max_shard_rows: int = MAX_SHARD_ROWS,
    max_shard_bytes: int = MAX_SHARD_BYTES,
) -> WorkingSetPublication:
    scope_sha256 = (filter_manifest.scope_policy_sha256, filter_manifest.scope_manifest_sha256)
    _check_filter_read_boundary(
        source_manifest=source_manifest,
        source_accounting=source_accounting,
        source_routing_context=source_routing_context,
        filter_root=filter_root,
        filter_manifest=filter_manifest,
        policy_sha256=policy_sha256,
        population_manifest_sha256=population_manifest_sha256,
        clinical_manifest_sha256=clinical_manifest_sha256,
        exception_manifest_sha256=exception_manifest_sha256,
        precision_manifest_sha256=precision_manifest_sha256,
        scope_policy_sha256=scope_policy_sha256, scope_manifest_sha256=scope_manifest_sha256,
    )
    filter_root, root = Path(filter_root), Path(output_root)
    if (
        not root.is_absolute()
        or not _count(max_shard_rows, 1) or max_shard_rows > MAX_SHARD_ROWS
        or not _count(max_shard_bytes, 1) or max_shard_bytes > MAX_SHARD_BYTES
    ):
        raise WgsFilterError("working-set publication binding is invalid")
    resources = (
        population_manifest_sha256, clinical_manifest_sha256,
        exception_manifest_sha256, precision_manifest_sha256,
    )
    state_path = (Path(retained_state_path) if retained_state_path is not None
                  else root / "working-state.sqlite")
    if retained_state_path is not None and (
        not state_path.is_absolute() or state_path.name != "filter.sqlite"
        or state_path.parent.name != "state"
    ):
        raise WgsFilterError("retained filter state path is invalid")
    shards: list[Shard] = []
    try:
        with stage_lock(root, "working-set"):
            state_context = (
                _FilterState.open_sealed(
                    state_path,
                    expected=source_accounting.total_physical_alternate_slots,
                    source_manifest_sha256=source_manifest.sha256,
                    policy_sha256=filter_manifest.policy_sha256,
                    resource_sha256=resources,
                    scope_sha256=scope_sha256,
                    reference_sha256=source_manifest.source.reference_manifest_sha256,
                    filter_manifest_sha256=filter_manifest.sha256,
                ) if retained_state_path is not None else _FilterState(
                    state_path,
                    expected=source_accounting.total_physical_alternate_slots,
                    source_manifest_sha256=source_manifest.sha256,
                    policy_sha256=filter_manifest.policy_sha256,
                    resource_sha256=resources,
                    scope_sha256=scope_sha256,
                    reference_sha256=source_manifest.source.reference_manifest_sha256,
                )
            )
            with state_context as state:
                if retained_state_path is None:
                    _replay_filter_rows(
                        state,
                        source_root=Path(source_root),
                        source_manifest=source_manifest,
                        source_routing_context=source_routing_context,
                        filter_root=filter_root,
                        filter_manifest=filter_manifest,
                    )
                    state.seal()
                summary = state.expected_working_summary()
                previous = None
                batch: list[WorkingSetRow] = []
                batch_bytes = 0
                first_row = 1
                for row in state.iter_expected_working_rows(filter_manifest.sha256):
                    size = len(row.as_bytes())
                    if size > max_shard_bytes:
                        raise WgsFilterError("working-set row exceeds its bounded size")
                    if batch and (
                        len(batch) >= max_shard_rows
                        or batch_bytes + size > max_shard_bytes
                    ):
                        previous = publish_working_set_shard(
                            root, tuple(batch), ordinal=len(shards) + 1,
                            first_record=first_row, previous=previous,
                        )
                        shards.append(previous)
                        first_row += len(batch)
                        batch, batch_bytes = [], 0
                    batch.append(row)
                    batch_bytes += size
                if batch:
                    previous = publish_working_set_shard(
                        root, tuple(batch), ordinal=len(shards) + 1,
                        first_record=first_row, previous=previous,
                    )
                    shards.append(previous)
            manifest = WorkingSetManifest(
                filter_manifest.sha256,
                _working_parents(source_manifest, filter_manifest),
                summary,
                tuple(shards),
            )
            require_published_shards(root, manifest.shards)
            summary_artifact = publish_bytes(
                root, FILTER_SUMMARY_ARTIFACT, summary.as_bytes(),
            )
            manifest_artifact = publish_bytes(
                root, MANIFEST_ARTIFACT, manifest.as_bytes(),
            )
            if retained_state_path is not None:
                _FilterState.bind_working(
                    state_path, filter_manifest_sha256=filter_manifest.sha256,
                    working_manifest_sha256=manifest.sha256,
                    working_summary_sha256=summary.sha256,
                )
            return WorkingSetPublication(manifest, manifest_artifact, summary_artifact)
    finally:
        if retained_state_path is None:
            try:
                state_path.unlink()
            except (FileNotFoundError, OSError):
                pass


@failure_boundary(
    "working_set_publication", "working_set_validate",
    known_errors=(WgsFilterError,),
)
def verify_working_set(
    *,
    source_root: str | Path,
    source_manifest: LedgerManifest,
    source_accounting: SourceAccountingClosure,
    source_routing_context: SourceRoutingContext,
    filter_root: str | Path,
    filter_manifest: FilterManifest,
    working_root: str | Path,
    working_manifest: WorkingSetManifest,
    policy_sha256: str,
    population_manifest_sha256: str,
    clinical_manifest_sha256: str,
    exception_manifest_sha256: str,
    precision_manifest_sha256: str,
    scope_policy_sha256: str, scope_manifest_sha256: str,
    retained_state_path: str | Path | None = None,
) -> WorkingSetSummary:
    scope_sha256 = (filter_manifest.scope_policy_sha256, filter_manifest.scope_manifest_sha256)
    _check_filter_read_boundary(
        source_manifest=source_manifest,
        source_accounting=source_accounting,
        source_routing_context=source_routing_context,
        filter_root=filter_root,
        filter_manifest=filter_manifest,
        policy_sha256=policy_sha256,
        population_manifest_sha256=population_manifest_sha256,
        clinical_manifest_sha256=clinical_manifest_sha256,
        exception_manifest_sha256=exception_manifest_sha256,
        precision_manifest_sha256=precision_manifest_sha256,
        scope_policy_sha256=scope_policy_sha256, scope_manifest_sha256=scope_manifest_sha256,
    )
    filter_root, root = Path(filter_root), Path(working_root)
    if (
        not root.is_absolute()
        or working_manifest.parent_filter_manifest_sha256 != filter_manifest.sha256
        or working_manifest.parents
        != _working_parents(source_manifest, filter_manifest)
    ):
        raise WgsFilterError("working-set parent closure is invalid")
    manifest_content = working_manifest.as_bytes()
    if (
        read_artifact(root, _artifact_for(MANIFEST_ARTIFACT, manifest_content))
        != manifest_content
        or read_artifact(
            root,
            _artifact_for(FILTER_SUMMARY_ARTIFACT, working_manifest.summary.as_bytes()),
        ) != working_manifest.summary.as_bytes()
    ):
        raise WgsFilterError("working-set manifest artifacts are invalid")
    resources = (
        population_manifest_sha256, clinical_manifest_sha256,
        exception_manifest_sha256, precision_manifest_sha256,
    )
    if retained_state_path is not None:
        retained = Path(retained_state_path)
        if (
            not retained.is_absolute()
            or retained.name != "filter.sqlite"
            or retained.parent.name != "state"
            or retained.parent.parent != Path(filter_root).parent
        ):
            raise WgsFilterError("retained working state path is invalid")
        with _FilterState.open_sealed(
            retained,
            expected=source_accounting.total_physical_alternate_slots,
            source_manifest_sha256=source_manifest.sha256,
            policy_sha256=filter_manifest.policy_sha256,
            resource_sha256=resources,
                    scope_sha256=scope_sha256,
                    reference_sha256=source_manifest.source.reference_manifest_sha256,
            filter_manifest_sha256=filter_manifest.sha256,
        ) as state:
            if state.working_binding() != (
                    working_manifest.sha256, working_manifest.summary.sha256):
                raise WgsFilterError("retained working state binding is invalid")
        for shard in working_manifest.shards:
            try:
                metadata = (root / shard.artifact.name).lstat()
            except OSError:
                raise WgsFilterError("retained working shard is unavailable") from None
            if (not stat.S_ISREG(metadata.st_mode)
                    or stat.S_IMODE(metadata.st_mode) != 0o600
                    or metadata.st_uid != os.getuid() or metadata.st_nlink != 1
                    or metadata.st_size != shard.artifact.size_bytes):
                raise WgsFilterError("retained working shard is unsafe")
        return working_manifest.summary
    state_path = root / "verify-working-state.sqlite"
    try:
        with stage_lock(root, "working-set-verify"):
            with _FilterState(
                state_path,
                expected=source_accounting.total_physical_alternate_slots,
                source_manifest_sha256=source_manifest.sha256,
                policy_sha256=filter_manifest.policy_sha256,
                resource_sha256=resources,
                    scope_sha256=scope_sha256,
                    reference_sha256=source_manifest.source.reference_manifest_sha256,
            ) as state:
                _replay_filter_rows(
                    state,
                    source_root=Path(source_root),
                    source_manifest=source_manifest,
                    source_routing_context=source_routing_context,
                    filter_root=filter_root,
                    filter_manifest=filter_manifest,
                )
                state.seal()
                if state.expected_working_summary() != working_manifest.summary:
                    raise WgsFilterError("working-set summary replay is invalid")
                actual = (
                    row
                    for shard in working_manifest.shards
                    for row in iter_working_set_shard(root, shard)
                )
                expected = state.iter_expected_working_rows(filter_manifest.sha256)
                for expected_row in expected:
                    try:
                        if next(actual) != expected_row:
                            raise WgsFilterError("working-set replay is invalid")
                    except StopIteration:
                        raise WgsFilterError("working-set replay is incomplete") from None
                try:
                    next(actual)
                except StopIteration:
                    pass
                else:
                    raise WgsFilterError("working-set replay has extra rows")
                return working_manifest.summary
    finally:
        try:
            state_path.unlink()
        except (FileNotFoundError, OSError):
            pass
