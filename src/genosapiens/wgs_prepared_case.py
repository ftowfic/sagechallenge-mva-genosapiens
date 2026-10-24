"""Immutable prepared-case authority and ledger-bound early snapshot."""

from __future__ import annotations

from dataclasses import asdict, dataclass, replace
import hashlib
import json
import os
from pathlib import Path
import re
import sqlite3
import stat
from typing import Literal

from genosapiens.bounded_wgs_config import (
    BoundedWgsExecutionConfig,
    bounded_wgs_execution_config_as_json,
)
from genosapiens.track1_installation import load_run_versions as load_base_release
from genosapiens.wgs_ledger import (
    DIAGNOSTIC_RECEIPT_SCHEMA,
    FILTER_MANIFEST_SCHEMA,
    HEADER_SCHEMA,
    INVENTORY_SCHEMA,
    JOURNAL_ENTRY_SCHEMA,
    SOURCE_MANIFEST_SCHEMA,
    StageRunBinding,
    WORKING_SET_MANIFEST_SCHEMA,
    Artifact,
    FilterManifest,
    LedgerManifest,
    ParentIdentity,
    PreparedCaseManifest,
    WorkingSetManifest,
    canonical,
    load_model,
    publish_bytes,
    read_artifact,
    read_journal,
)
from genosapiens.wgs_preparation import (
    COMPLETED_SOURCE_SCHEMA, CompletedSourceArtifact, ResolvedCompletedSource,
    PreparationRunOptions, SourceAccountingClosure, SourceAuditResult,
    load_preparation_run_options, load_source_audit, load_completed_source_input,
    resolve_completed_source,
)


EARLY_SNAPSHOT_SCHEMA = "genosapiens.wgs-ledger-bound-early-snapshot"
PREPARED_INDEX_SCHEMA = "genosapiens.wgs-prepared-index"
PREPARED_INDEX_NAME = "prepared-index.json"
VERSION = 1
HANDOFF_VERSION = 4
PREPARED_INDEX_VERSION = 3
PREPARED_MANIFEST_NAME = "prepared-case-manifest.json"
DIAGNOSTIC_NAME = "preparation-diagnostics.json"
INVENTORY_NAME = "inventory-preparation.json"
EARLY_SNAPSHOT_NAME = "ledger-bound-early-snapshot.json"
_DIGEST = re.compile(r"^[0-9a-f]{64}$")
_IMAGE = re.compile(r"^sha256:[0-9a-f]{64}$")
_NAME = re.compile(r"^[a-z0-9][a-z0-9._-]{0,127}$")
_ERRORS = (
    ValueError, TypeError, RuntimeError, AttributeError, KeyError, IndexError,
    OSError, UnicodeError, OverflowError, RecursionError,
)


class WgsPreparedCaseError(RuntimeError):
    """Static prepared-case failure without controlled values or paths."""


def _digest(value: object) -> bool:
    return type(value) is str and _DIGEST.fullmatch(value) is not None


def _bytes(value: object) -> bytes:
    return canonical(value)


@dataclass(frozen=True, slots=True, repr=False)
class PreparedStateDescriptor:
    """Measured private state and its completed writer commitment, not a page hash."""

    stage: str
    relative_path: str
    schema: str
    user_version: int
    size_bytes: int
    run_binding_sha256: str | None
    committed_cursor: tuple[int, ...]
    committed_counts: tuple[tuple[str, int], ...]
    completed_stage_binding_sha256: str
    storage_owner: Literal["completed_source", "preparation"]

    def __post_init__(self) -> None:
        expected = {
            "source": ("source.sqlite", "genosapiens.wgs-stage-checkpoint", 1),
            "source_audit": ("source-audit.sqlite", "genosapiens.wgs-stage-checkpoint", 1),
            "prefilter": ("filter.sqlite", "genosapiens.wgs-stage-checkpoint", 2),
        }
        if (
            self.stage not in expected
            or self.storage_owner != ("preparation" if self.stage == "prefilter" else "completed_source")
            or (self.relative_path, self.schema, self.user_version) != expected[self.stage]
            or type(self.size_bytes) is not int or self.size_bytes <= 0
            or not _digest(self.completed_stage_binding_sha256)
            or not _digest(self.run_binding_sha256)
            or type(self.committed_cursor) is not tuple
            or any(type(value) is not int or value < 0 for value in self.committed_cursor)
            or not self.committed_cursor
            or type(self.committed_counts) is not tuple
            or any(type(name) is not str or _NAME.fullmatch(name) is None
                   or type(count) is not int or count < 0
                   for name, count in self.committed_counts)
            or self.committed_counts != tuple(sorted(self.committed_counts))
            or len({name for name, _ in self.committed_counts}) != len(self.committed_counts)
        ):
            raise WgsPreparedCaseError("prepared state descriptor is invalid")


@dataclass(frozen=True, slots=True, repr=False)
class PreparedIndex:
    """Retained computation descriptor with explicit stage and product parents."""

    source_manifest_sha256: str
    filter_manifest_sha256: str
    working_manifest_sha256: str
    source_accounting_sha256: str
    source_audit_sha256: str
    preparation_run_options_sha256: str
    states: tuple[PreparedStateDescriptor, ...]
    completed_stage_bindings: tuple[tuple[str, str], ...]
    completed_source_input_sha256: str
    completed_source_state_root: str
    preparation_state_root: str
    schema: Literal["genosapiens.wgs-prepared-index"] = PREPARED_INDEX_SCHEMA
    version: Literal[3] = PREPARED_INDEX_VERSION
    classification: Literal["controlled"] = "controlled"
    llm_policy: Literal["no-llm"] = "no-llm"

    def __post_init__(self) -> None:
        stage_names = tuple(value.stage for value in self.states)
        stages = tuple(name for name, _ in self.completed_stage_bindings)
        if (
            any(not _digest(value) for value in (
                self.source_manifest_sha256, self.filter_manifest_sha256,
                self.working_manifest_sha256, self.source_accounting_sha256,
                self.source_audit_sha256, self.preparation_run_options_sha256,
                self.completed_source_input_sha256,
            ))
            or type(self.states) is not tuple
            or any(type(value) is not PreparedStateDescriptor for value in self.states)
            or stage_names != tuple(sorted(set(stage_names)))
            or not {"source", "prefilter"}.issubset(stage_names)
            or type(self.completed_stage_bindings) is not tuple
            or stages != tuple(sorted(set(stages)))
            or set(stages) != {"source", "source_audit", "prefilter"}
            or any(not _digest(digest) for _, digest in self.completed_stage_bindings)
            or not set(stage_names).issubset(stages)
            or self.schema != PREPARED_INDEX_SCHEMA or self.version != PREPARED_INDEX_VERSION
            or any(type(path) is not str or not Path(path).is_absolute()
                   or ".." in Path(path).parts for path in (
                       self.completed_source_state_root, self.preparation_state_root))
            or self.classification != "controlled" or self.llm_policy != "no-llm"
        ):
            raise WgsPreparedCaseError("prepared index is invalid")

    def as_bytes(self) -> bytes:
        return _bytes({
            "classification": self.classification,
            "data": asdict(self),
            "kind": "prepared_index",
            "llm_policy": self.llm_policy,
            "schema": self.schema,
            "version": self.version,
        })

    @property
    def sha256(self) -> str:
        return hashlib.sha256(self.as_bytes()).hexdigest()

    def state_path(self, descriptor: PreparedStateDescriptor) -> Path:
        if descriptor not in self.states:
            raise WgsPreparedCaseError("prepared state is undeclared")
        root = (self.completed_source_state_root if descriptor.storage_owner == "completed_source"
                else self.preparation_state_root)
        return Path(root) / descriptor.relative_path


def _unique_pairs(items):
    value = {}
    for key, child in items:
        if key in value:
            raise ValueError
        value[key] = child
    return value


def _state_descriptor(
    state_root: Path, stage: str, stage_parent_sha256: str,
    completed_stage_sha256: str,
    *, source_manifest: LedgerManifest, filter_manifest: FilterManifest,
    working_manifest: WorkingSetManifest,
    preparation_run_options_sha256: str,
    completed_source: ResolvedCompletedSource,
    executable_stack_sha256: str,
    execution_config_sha256: str,
    runtime_image_digest: str,
) -> PreparedStateDescriptor:
    if stage in {"source", "source_audit"}:
        descriptor = completed_source.descriptor
        state = descriptor.source_state if stage == "source" else descriptor.source_audit_state
        journal = next(item for item in completed_source.journals if item.stage == stage)
        if (state is None or state_root != Path(descriptor.source_state_root)
                or journal.stage_parent_sha256 != stage_parent_sha256
                or journal.head_sha256 != completed_stage_sha256):
            raise WgsPreparedCaseError("completed source state ownership is invalid")
        counts = state.committed_counts if stage == "source" else (
            ("records", dict(state.committed_counts)["records"]),)
        return PreparedStateDescriptor(stage, state.relative_path, state.schema,
            state.user_version, state.size_bytes, state.run_binding_sha256,
            state.committed_cursor, counts, state.completed_event_sha256, "completed_source")
    name = {"source": "source.sqlite", "source_audit": "source-audit.sqlite",
            "prefilter": "filter.sqlite"}[stage]
    path = state_root / name
    connection = None
    try:
        metadata = path.lstat()
        if (not state_root.is_absolute() or path.is_symlink()
                or not stat.S_ISREG(metadata.st_mode)
                or stat.S_IMODE(metadata.st_mode) != 0o600
                or metadata.st_uid != os.getuid() or metadata.st_nlink != 1
                or metadata.st_size <= 0):
            raise ValueError
        connection = sqlite3.connect(path.as_uri() + "?mode=ro&immutable=1", uri=True)
        connection.execute("PRAGMA query_only=ON")
        connection.execute("PRAGMA trusted_schema=OFF")
        version = connection.execute("PRAGMA user_version").fetchone()[0]
        if version != (2 if stage == "prefilter" else 1):
            raise ValueError
        binding_row = connection.execute(
            "SELECT canonical_bytes FROM run_binding WHERE singleton=1"
        ).fetchone()
        if binding_row is None:
            raise ValueError
        binding = load_model(binding_row[0], StageRunBinding)
        source_parents = {parent.role: parent.sha256 for parent in source_manifest.parents}
        if (binding.stage != stage
                or binding.stage_parent_sha256 != stage_parent_sha256
                or binding.input_sha256 != source_manifest.source.input_sha256
                or binding.run_options_sha256 != preparation_run_options_sha256
                or binding.executable_stack_sha256 != executable_stack_sha256
                or binding.execution_config_sha256 != execution_config_sha256
                or binding.runtime_image_digest != runtime_image_digest
                or {parent.role: parent.sha256 for parent in binding.parents}.get("completed_source_input")
                != completed_source.descriptor.sha256
                or binding.reference_sha256
                != source_manifest.source.reference_manifest_sha256):
            raise ValueError
        latest = connection.execute(
            "SELECT event_sha256,event_type FROM events ORDER BY ordinal DESC LIMIT 1"
        ).fetchone()
        if latest is None or latest[1] != "stage_completed" or not _digest(latest[0]):
            raise ValueError
        cursor = ()
        counts: dict[str, int] = {}
        ordinal = 0
        for commit_ordinal, cursor_bytes, counts_bytes in connection.execute(
            "SELECT ordinal,cursor_bytes,counts_bytes FROM commits ORDER BY ordinal"
        ):
            decoded_cursor = json.loads(cursor_bytes, object_pairs_hook=_unique_pairs)
            decoded_counts = json.loads(counts_bytes, object_pairs_hook=_unique_pairs)
            if (type(commit_ordinal) is not int or commit_ordinal <= ordinal
                    or type(decoded_cursor) is not list
                    or len(decoded_cursor) != (3 if stage == "prefilter" else 2)
                    or any(type(value) is not int or value < 1 for value in decoded_cursor)
                    or tuple(decoded_cursor) <= cursor
                    or type(decoded_counts) is not dict
                    or any(_NAME.fullmatch(key) is None or type(value) is not int
                           or value < 0 for key, value in decoded_counts.items())):
                raise ValueError
            ordinal = commit_ordinal
            cursor = tuple(decoded_cursor)
            for key, value in decoded_counts.items():
                counts[key] = counts.get(key, 0) + value
        if not cursor:
            raise ValueError
        if stage == "prefilter":
            seal = connection.execute(
                "SELECT sealed,committed_occurrences,filter_manifest_sha256,"
                "working_manifest_sha256 FROM filter_meta WHERE singleton=1"
            ).fetchone()
            if (seal is None or seal[0] != 1 or seal[1] != counts.get("decisions")
                    or seal[2] != filter_manifest.sha256
                    or seal[3] != working_manifest.sha256):
                raise ValueError
        return PreparedStateDescriptor(stage, name,
            "genosapiens.wgs-stage-checkpoint", version, metadata.st_size,
            hashlib.sha256(binding_row[0]).hexdigest(), cursor,
            tuple(sorted(counts.items())), latest[0], "preparation")
    except _ERRORS + (sqlite3.Error,):
        raise WgsPreparedCaseError("completed prepared state is unavailable or stale") from None
    finally:
        if connection is not None:
            connection.close()


def build_prepared_index(
    *, state_root: Path, source_manifest: LedgerManifest,
    filter_manifest: FilterManifest, working_manifest: WorkingSetManifest,
    source_accounting_sha256: str, source_audit_sha256: str,
    preparation_run_options_sha256: str,
    journals: tuple["CompletedJournalBinding", ...], audit_requested: bool,
    completed_source: ResolvedCompletedSource,
    executable_stack_sha256: str,
    execution_config_sha256: str,
    runtime_image_digest: str,
) -> PreparedIndex:
    """Record completed writer metadata without hashing retained SQLite pages."""
    try:
        heads = {value.stage: value.head_sha256 for value in journals}
        stage_parents = {value.stage: value.stage_parent_sha256 for value in journals}
        if (type(audit_requested) is not bool or not state_root.is_absolute()
                or set(heads) != {"source", "source_audit", "prefilter"}
                or len(heads) != len(journals)):
            raise ValueError
        names = ("prefilter", "source") + (("source_audit",) if audit_requested else ())
        states = tuple(sorted((
            _state_descriptor(state_root if name == "prefilter"
                              else Path(completed_source.descriptor.source_state_root),
                name, stage_parents[name], heads[name],
                source_manifest=source_manifest, filter_manifest=filter_manifest,
                working_manifest=working_manifest,
                preparation_run_options_sha256=preparation_run_options_sha256,
                completed_source=completed_source,
                executable_stack_sha256=executable_stack_sha256,
                execution_config_sha256=execution_config_sha256,
                runtime_image_digest=runtime_image_digest)
                for name in names
        ), key=lambda value: value.stage))
        by_stage = {value.stage: value for value in states}
        if (dict(by_stage["source"].committed_counts) != {
                "records": source_manifest.summary.total_records,
                "physical_alternates": source_manifest.summary.total_physical_alternate_slots,
                "biological_alternates": source_manifest.summary.total_biological_alternates,
            } or by_stage["source"].committed_cursor != (
                source_manifest.summary.total_records, len(source_manifest.shards))
            or dict(by_stage["prefilter"].committed_counts) != {
                "decisions": filter_manifest.summary.total_source_occurrences,
                "population_request_keys": filter_manifest.summary.population_request_keys,
                "clinical_query_keys": filter_manifest.summary.clinical_query_keys,
            } or by_stage["prefilter"].committed_cursor[-1] != len(filter_manifest.shards)
            or filter_manifest.parent_source_manifest_sha256 != source_manifest.sha256
            or working_manifest.parent_filter_manifest_sha256 != filter_manifest.sha256):
            raise ValueError
        if audit_requested and (
                by_stage["source_audit"].committed_cursor != (
                    source_manifest.summary.total_records, len(source_manifest.shards))
                or dict(by_stage["source_audit"].committed_counts) != {
                    "records": source_manifest.summary.total_records,
                }):
            raise ValueError
        return PreparedIndex(source_manifest.sha256, filter_manifest.sha256,
            working_manifest.sha256, source_accounting_sha256,
            source_audit_sha256, preparation_run_options_sha256, states,
            tuple(sorted(heads.items())), completed_source.descriptor.sha256,
            completed_source.descriptor.source_state_root, str(state_root))
    except WgsPreparedCaseError:
        raise
    except _ERRORS:
        raise WgsPreparedCaseError("prepared index binding is invalid") from None


def load_prepared_index(content: bytes) -> PreparedIndex:
    try:
        value = json.loads(content, object_pairs_hook=_unique_pairs)
        if (type(value) is not dict or set(value) != {
                "classification", "data", "kind", "llm_policy", "schema", "version"}
                or value["schema"] != PREPARED_INDEX_SCHEMA
                or value["version"] != PREPARED_INDEX_VERSION or value["kind"] != "prepared_index"
                or value["classification"] != "controlled"
                or value["llm_policy"] != "no-llm"):
            raise ValueError
        data = value["data"]
        if type(data) is not dict:
            raise ValueError
        data["states"] = tuple(PreparedStateDescriptor(
            **{**row,
               "committed_cursor": tuple(row["committed_cursor"]),
               "committed_counts": tuple(tuple(count) for count in row["committed_counts"])}
        ) for row in data["states"])
        data["completed_stage_bindings"] = tuple(
            tuple(item) for item in data["completed_stage_bindings"])
        result = PreparedIndex(**data)
        if result.as_bytes() != content:
            raise ValueError
        return result
    except _ERRORS:
        raise WgsPreparedCaseError("prepared index envelope is invalid") from None


def publish_prepared_index(state_root: Path, index: PreparedIndex) -> Artifact:
    """Publish an immutable state descriptor after completed state is inspected."""
    try:
        if (type(index) is not PreparedIndex or not state_root.is_absolute()
                or state_root.name != "state"):
            raise ValueError
        metadata = state_root.lstat()
        if (state_root.is_symlink() or not stat.S_ISDIR(metadata.st_mode)
                or stat.S_IMODE(metadata.st_mode) != 0o700
                or metadata.st_uid != os.getuid()):
            raise ValueError
        return _publish_or_reuse(state_root, PREPARED_INDEX_NAME, index.as_bytes())[0]
    except WgsPreparedCaseError:
        raise
    except _ERRORS:
        raise WgsPreparedCaseError("prepared index publication failed") from None


@dataclass(frozen=True, slots=True, repr=False)
class PreparationDiagnosticReceipt:
    case_manifest_sha256: str
    selected_sample_sha256: str
    concordance_sha256: str
    source_manifest_sha256: str
    source_accounting_sha256: str
    prefilter_manifest_sha256: str
    working_set_manifest_sha256: str
    hpo_input_sha256: str
    phenotype_ingestion_sha256: str
    ontology_manifest_sha256: str
    diagnostic_policy_sha256: str
    source_status: Literal["complete"]
    concordance_status: Literal["complete"]
    qc_status: Literal["passed"]
    class_and_technical_status: Literal["complete_original_source"]
    prefilter_status: Literal["complete"]
    hpo_status: Literal["validated"]
    blockers: tuple[str, ...]
    source_audit_sha256: str
    preparation_run_options_sha256: str
    prepared_index_sha256: str
    audit_requested: bool
    audit_status: Literal["not_requested", "passed"]
    completed_source_input_sha256: str
    schema: Literal["genosapiens.wgs-preparation-diagnostics"] = (
        DIAGNOSTIC_RECEIPT_SCHEMA
    )
    version: Literal[4] = HANDOFF_VERSION
    status: Literal["passed"] = "passed"
    classification: Literal["controlled"] = "controlled"
    llm_policy: Literal["no-llm"] = "no-llm"

    def __post_init__(self) -> None:
        if (
            any(not _digest(value) for value in (
                self.case_manifest_sha256,
                self.selected_sample_sha256,
                self.concordance_sha256,
                self.source_manifest_sha256,
                self.source_accounting_sha256,
                self.source_audit_sha256,
                self.preparation_run_options_sha256,
                self.prepared_index_sha256,
                self.completed_source_input_sha256,
                self.prefilter_manifest_sha256,
                self.working_set_manifest_sha256,
                self.hpo_input_sha256,
                self.phenotype_ingestion_sha256,
                self.ontology_manifest_sha256,
                self.diagnostic_policy_sha256,
            ))
            or self.blockers
            or self.schema != DIAGNOSTIC_RECEIPT_SCHEMA
            or self.version != HANDOFF_VERSION
            or type(self.audit_requested) is not bool
            or self.audit_status != ("passed" if self.audit_requested else "not_requested")
            or self.status != "passed"
            or self.classification != "controlled"
            or self.llm_policy != "no-llm"
            or self.source_status != "complete"
            or self.concordance_status != "complete"
            or self.qc_status != "passed"
            or self.class_and_technical_status != "complete_original_source"
            or self.prefilter_status != "complete"
            or self.hpo_status != "validated"
        ):
            raise WgsPreparedCaseError("preparation diagnostics are not passing")

    def as_bytes(self) -> bytes:
        return _bytes({
            "classification": "controlled",
            "data": asdict(self),
            "kind": "preparation_diagnostics",
            "llm_policy": "no-llm",
            "schema": self.schema,
            "version": self.version,
        })

    @property
    def sha256(self) -> str:
        return hashlib.sha256(self.as_bytes()).hexdigest()


@dataclass(frozen=True, slots=True, repr=False)
class PreparationInventoryEntry:
    role: str
    parent_schema: str
    sha256: str | None
    size_bytes: int | None
    kind: Literal["canonical_artifact", "identity_binding", "trusted_stage_state"]
    relative_path: str | None = None
    stage: str | None = None
    user_version: int | None = None
    commit_binding_sha256: str | None = None
    run_binding_sha256: str | None = None
    storage_owner: Literal["completed_source", "preparation"] = "preparation"

    def __post_init__(self) -> None:
        if (
            type(self.role) is not str
            or self.storage_owner not in {"completed_source", "preparation"}
            or _NAME.fullmatch(self.role) is None
            or type(self.parent_schema) is not str
            or not self.parent_schema.startswith("genosapiens.")
            or self.kind not in {"canonical_artifact", "identity_binding", "trusted_stage_state"}
            or (self.kind == "canonical_artifact" and (
                not _digest(self.sha256) or type(self.size_bytes) is not int
                or self.size_bytes <= 0))
            or (self.kind == "identity_binding" and (
                not _digest(self.sha256) or self.size_bytes is not None))
            or (self.kind == "trusted_stage_state" and (
                self.sha256 is not None or type(self.size_bytes) is not int
                or self.size_bytes <= 0 or self.stage not in {"source", "source_audit", "prefilter"}
                or self.relative_path != "state/" + {
                    "source": "source.sqlite", "source_audit": "source-audit.sqlite",
                    "prefilter": "filter.sqlite",
                }[self.stage]
                or type(self.user_version) is not int or self.user_version != (2 if self.stage == "prefilter" else 1)
                or not _digest(self.commit_binding_sha256)
                or (self.run_binding_sha256 is not None and not _digest(self.run_binding_sha256))))
            or (self.kind != "trusted_stage_state" and any(value is not None for value in (
                self.relative_path, self.stage, self.user_version,
                self.commit_binding_sha256, self.run_binding_sha256)))
        ):
            raise WgsPreparedCaseError("preparation inventory entry is invalid")


@dataclass(frozen=True, slots=True, repr=False)
class PreparationInventory:
    case_manifest_sha256: str
    entries: tuple[PreparationInventoryEntry, ...]
    journal_stages: tuple[str, ...]
    completed_source_input_sha256: str
    completed_source_input_path: str
    completed_source_owner_root: str
    preparation_owner_root: str
    schema: Literal["genosapiens.wgs-preparation-inventory"] = INVENTORY_SCHEMA
    version: Literal[4] = HANDOFF_VERSION
    status: Literal["complete"] = "complete"
    classification: Literal["controlled"] = "controlled"
    llm_policy: Literal["no-llm"] = "no-llm"

    def __post_init__(self) -> None:
        roles = tuple(value.role for value in self.entries)
        required = {
            "case_manifest", "original_input", "header", "source_manifest",
            "qc_summary", "normalization_summary", "class_summary",
            "prefilter_policy", "prefilter_manifest", "filter_summary",
            "working_set_manifest", "hpo_input", "ontology_manifest",
            "phenotype_ingestion",
            "population_manifest", "clinical_manifest", "exception_manifest",
            "precision_manifest", "early_scope_policy", "early_scope_manifest",
            "diagnostic_receipt", "execution_limits", "execution_authority",
            "executable_stack", "runtime_image",
            "source_accounting", "source_audit", "preparation_run_options",
            "prepared_index", "state_source", "state_prefilter",
            "completed_source_input",
        }
        if (
            not _digest(self.case_manifest_sha256)
            or not _digest(self.completed_source_input_sha256)
            or any(type(path) is not str or not Path(path).is_absolute()
                   or ".." in Path(path).parts for path in (
                       self.completed_source_input_path, self.completed_source_owner_root,
                       self.preparation_owner_root))
            or self.schema != INVENTORY_SCHEMA
            or self.version != HANDOFF_VERSION
            or self.status != "complete"
            or self.classification != "controlled"
            or self.llm_policy != "no-llm"
            or self.entries
            != tuple(sorted(self.entries, key=lambda value: value.role))
            or len(roles) != len(set(roles))
            or not required.issubset(roles)
            or self.journal_stages
            != tuple(sorted(set(self.journal_stages)))
            or not self.journal_stages
            or any(_NAME.fullmatch(value) is None for value in self.journal_stages)
            or {f"journal_{value}" for value in self.journal_stages}
            - set(roles)
        ):
            raise WgsPreparedCaseError("preparation inventory is incomplete")

    def as_bytes(self) -> bytes:
        return _bytes({
            "classification": "controlled",
            "data": asdict(self),
            "kind": "preparation_inventory",
            "llm_policy": "no-llm",
            "schema": self.schema,
            "version": self.version,
        })

    @property
    def sha256(self) -> str:
        return hashlib.sha256(self.as_bytes()).hexdigest()


@dataclass(frozen=True, slots=True, repr=False)
class CompletedJournalBinding:
    stage: str
    root: Path
    stage_parent_sha256: str
    head_sha256: str

    def __post_init__(self) -> None:
        if (
            type(self.stage) is not str
            or _NAME.fullmatch(self.stage) is None
            or not self.root.is_absolute()
            or not _digest(self.stage_parent_sha256)
            or not _digest(self.head_sha256)
        ):
            raise WgsPreparedCaseError("completed journal binding is invalid")

    def __repr__(self) -> str:
        return "<controlled completed journal; no-llm>"


@dataclass(frozen=True, slots=True, repr=False)
class LedgerBoundEarlySnapshot:
    prepared_case_sha256: str
    case_manifest_sha256: str
    source_manifest_sha256: str
    prefilter_manifest_sha256: str
    working_set_manifest_sha256: str
    diagnostic_receipt_sha256: str
    inventory_sha256: str
    hpo_input_sha256: str
    ontology_manifest_sha256: str
    execution_limits_sha256: str
    base_release_sha256: str | None
    executable_stack_sha256: str
    runtime_image_digest: str
    source_occurrences: int
    retained_occurrences: int
    excluded_common_occurrences: int
    excluded_scope_occurrences: int
    unresolved_scope_occurrences: int
    source_accounting_sha256: str
    source_audit_sha256: str
    preparation_run_options_sha256: str
    prepared_index_sha256: str
    completed_source_input_sha256: str
    public_validation_authority_sha256: str | None = None
    schema: Literal["genosapiens.wgs-ledger-bound-early-snapshot"] = (
        EARLY_SNAPSHOT_SCHEMA
    )
    version: Literal[4] = HANDOFF_VERSION
    preparation_status: Literal["passed"] = "passed"
    analysis_authority: Literal["prepared-case-manifest"] = "prepared-case-manifest"
    evidence_storage: Literal["authenticated-shards-no-embedded-document"] = (
        "authenticated-shards-no-embedded-document"
    )
    competitive_eligibility: Literal["not-established"] = "not-established"
    classification: Literal["controlled"] = "controlled"
    llm_policy: Literal["no-llm"] = "no-llm"

    def __post_init__(self) -> None:
        if (
            any(not _digest(value) for value in (
                self.prepared_case_sha256,
                self.case_manifest_sha256,
                self.source_manifest_sha256,
                self.prefilter_manifest_sha256,
                self.working_set_manifest_sha256,
                self.diagnostic_receipt_sha256,
                self.inventory_sha256,
                self.hpo_input_sha256,
                self.ontology_manifest_sha256,
                self.execution_limits_sha256,
                self.executable_stack_sha256,
                self.source_accounting_sha256,
                self.source_audit_sha256,
                self.preparation_run_options_sha256,
                self.prepared_index_sha256,
                self.completed_source_input_sha256,
            ))
            or (self.base_release_sha256 is None) == (self.public_validation_authority_sha256 is None)
            or not _digest(self.base_release_sha256 or self.public_validation_authority_sha256)
            or _IMAGE.fullmatch(self.runtime_image_digest) is None
            or any(type(value) is not int or value < 0 for value in (
                self.source_occurrences,
                self.retained_occurrences,
                self.excluded_common_occurrences,
                self.excluded_scope_occurrences, self.unresolved_scope_occurrences,
            ))
            or (self.retained_occurrences + self.excluded_common_occurrences
                + self.excluded_scope_occurrences + self.unresolved_scope_occurrences)
            > self.source_occurrences
            or self.schema != EARLY_SNAPSHOT_SCHEMA
            or self.version != HANDOFF_VERSION
            or self.preparation_status != "passed"
            or self.analysis_authority != "prepared-case-manifest"
            or self.evidence_storage
            != "authenticated-shards-no-embedded-document"
            or self.competitive_eligibility != "not-established"
            or self.classification != "controlled"
            or self.llm_policy != "no-llm"
        ):
            raise WgsPreparedCaseError("ledger-bound early snapshot is invalid")

    def as_bytes(self) -> bytes:
        return _bytes({
            "classification": "controlled",
            "data": asdict(self),
            "kind": "ledger_bound_early_snapshot",
            "llm_policy": "no-llm",
            "schema": self.schema,
            "version": self.version,
        })

    @property
    def sha256(self) -> str:
        return hashlib.sha256(self.as_bytes()).hexdigest()


@dataclass(frozen=True, slots=True, repr=False)
class PreparedCasePublication:
    root: Path
    manifest: PreparedCaseManifest
    diagnostic: PreparationDiagnosticReceipt
    inventory: PreparationInventory
    early_snapshot: LedgerBoundEarlySnapshot
    prepared_index: PreparedIndex
    artifacts: tuple[Artifact, ...]
    reused: bool
    completed_source: ResolvedCompletedSource

    def __post_init__(self) -> None:
        if (
            not self.root.is_absolute()
            or type(self.completed_source) is not ResolvedCompletedSource
            or self.completed_source.descriptor.sha256 != self.manifest.completed_source_input.sha256
            or tuple(value.name for value in self.artifacts) != (
                DIAGNOSTIC_NAME,
                INVENTORY_NAME,
                PREPARED_MANIFEST_NAME,
                EARLY_SNAPSHOT_NAME,
            )
            or self.early_snapshot.prepared_case_sha256 != self.manifest.sha256
            or self.early_snapshot.diagnostic_receipt_sha256
            != self.diagnostic.sha256
            or self.early_snapshot.inventory_sha256 != self.inventory.sha256
            or self.manifest.prepared_index.sha256 != self.prepared_index.sha256
            or self.early_snapshot.prepared_index_sha256 != self.prepared_index.sha256
        ):
            raise WgsPreparedCaseError("prepared-case publication is invalid")

    def __repr__(self) -> str:
        return "<controlled prepared case; no-llm>"


def _source_accounting_sha256(value: SourceAccountingClosure) -> str:
    return hashlib.sha256(_bytes(asdict(value))).hexdigest()


def _verify_source_accounting(
    source_manifest: LedgerManifest,
    source_accounting: SourceAccountingClosure,
) -> None:
    reports = {value.name: value for value in source_manifest.summary.reports}
    if (
        set(reports) != {"qc.json", "normalization.json", "classes.json"}
        or source_accounting.source_manifest_sha256 != source_manifest.sha256
        or source_accounting.input_sha256 != source_manifest.source.input_sha256
        or source_accounting.header_sha256 != source_manifest.header.sha256
        or source_accounting.qc_summary_sha256 != reports["qc.json"].sha256
        or source_accounting.normalization_summary_sha256
        != reports["normalization.json"].sha256
        or source_accounting.class_summary_sha256 != reports["classes.json"].sha256
        or source_accounting.total_records != source_manifest.summary.total_records
        or source_accounting.total_physical_alternate_slots
        != source_manifest.summary.total_physical_alternate_slots
        or source_accounting.total_biological_alternates
        != source_manifest.summary.total_biological_alternates
    ):
        raise WgsPreparedCaseError("source accounting is stale")


def _journal_parents(
    journals: tuple[CompletedJournalBinding, ...],
    *,
    execution_limits_sha256: str,
    execution_limits_version: int,
    executable_stack_sha256: str,
    runtime_image_digest: str,
    completed_source: ResolvedCompletedSource,
) -> tuple[ParentIdentity, ...]:
    if (
        journals != tuple(sorted(journals, key=lambda value: value.stage))
        or len({value.stage for value in journals}) != len(journals)
        or not journals
    ):
        raise WgsPreparedCaseError("prepared journal set is invalid")
    values = []
    runtime_sha256 = runtime_image_digest.removeprefix("sha256:")
    for journal in journals:
        if journal.stage in {"source", "source_audit"}:
            if journal not in completed_source.journals:
                raise WgsPreparedCaseError("prepared source journal ownership is stale")
            values.append(ParentIdentity(f"{journal.stage}_journal", JOURNAL_ENTRY_SCHEMA, 1, journal.head_sha256))
            continue
        from genosapiens.wgs_ledger import read_journal_head
        terminal = read_journal_head(
            journal.root,
            stage=journal.stage,
            parent_sha256=journal.stage_parent_sha256,
            expected_head_sha256=journal.head_sha256,
        )
        if terminal.event_type != "stage_completed":
            raise WgsPreparedCaseError("prepared journal is not complete")
        parents = {value.role: value.sha256 for value in terminal.parents}
        parent_models = {value.role: value for value in terminal.parents}
        if (
            parents.get("execution_limits") != execution_limits_sha256
            or parents.get("completed_source_input") != completed_source.descriptor.sha256
            or parents.get("executable_stack") != executable_stack_sha256
            or parents.get("runtime_image") != runtime_sha256
            or parent_models["execution_limits"].parent_schema
            != "genosapiens.bounded-wgs-execution"
            or parent_models["execution_limits"].parent_version
            != execution_limits_version
            or parent_models["executable_stack"].parent_schema
            != "genosapiens.executable-stack"
            or parent_models["executable_stack"].parent_version != VERSION
            or parent_models["runtime_image"].parent_schema
            != "genosapiens.runtime-image"
            or parent_models["runtime_image"].parent_version != VERSION
        ):
            raise WgsPreparedCaseError("prepared journal bindings are stale")
        values.append(ParentIdentity(
            f"{journal.stage}_journal",
            JOURNAL_ENTRY_SCHEMA,
            1,
            journal.head_sha256,
        ))
    return tuple(values)


def _resource_parents(filter_manifest: FilterManifest) -> tuple[ParentIdentity, ...]:
    by_role = {value.role: value for value in filter_manifest.parents}
    return tuple(
        ParentIdentity(
            role,
            by_role[role].parent_schema,
            by_role[role].parent_version,
            by_role[role].sha256,
        )
        for role in (
            "population_manifest", "clinical_manifest",
            "exception_manifest", "precision_manifest", "early_scope_policy", "early_scope_manifest",
        )
    )


def _summary_parents(
    source: LedgerManifest,
    filter_manifest: FilterManifest,
) -> tuple[ParentIdentity, ...]:
    reports = {value.name: value for value in source.summary.reports}
    return (
        ParentIdentity("qc_summary", "genosapiens.wgs-qc-summary", 1,
                       reports["qc.json"].sha256),
        ParentIdentity("normalization_summary", "genosapiens.wgs-normalization-summary", 1,
                       reports["normalization.json"].sha256),
        ParentIdentity("class_summary", "genosapiens.wgs-class-summary", 1,
                       reports["classes.json"].sha256),
        ParentIdentity("filter_summary", "genosapiens.wgs-filter-summary", 2,
                       filter_manifest.summary.sha256),
    )


def _entry(role: str, schema: str, sha256: str, size: int | None = None) -> PreparationInventoryEntry:
    return PreparationInventoryEntry(
        role, schema, sha256, size,
        "canonical_artifact" if size is not None else "identity_binding",
    )


def _state_entry(state: PreparedStateDescriptor) -> PreparationInventoryEntry:
    return PreparationInventoryEntry(
        "state_" + state.stage, state.schema, None, state.size_bytes,
        "trusted_stage_state", "state/" + state.relative_path, state.stage,
        state.user_version, state.completed_stage_binding_sha256,
        state.run_binding_sha256,
        state.storage_owner,
    )


def _build_inventory(
    manifest: PreparedCaseManifest,
    source: LedgerManifest,
    filter_manifest: FilterManifest,
    working: WorkingSetManifest,
    diagnostic: PreparationDiagnosticReceipt,
    index: PreparedIndex,
    source_accounting: SourceAccountingClosure,
    audit: SourceAuditResult,
    options: PreparationRunOptions,
    completed_source: ResolvedCompletedSource,
    preparation_root: Path,
) -> PreparationInventory:
    summaries = {value.role: value for value in manifest.summaries}
    entries = [
        _entry("completed_source_input", COMPLETED_SOURCE_SCHEMA,
               completed_source.descriptor.sha256, len(completed_source.descriptor.as_bytes())),
        _entry("case_manifest", "genosapiens.case-manifest", manifest.case_manifest_sha256),
        _entry("original_input", "genosapiens.original-input", manifest.original_input_sha256,
               manifest.original_input_size_bytes),
        _entry("header", HEADER_SCHEMA, manifest.header.sha256, source.header.size_bytes),
        _entry("source_manifest", SOURCE_MANIFEST_SCHEMA, source.sha256,
               len(source.as_bytes())),
        _entry("qc_summary", summaries["qc_summary"].parent_schema,
               summaries["qc_summary"].sha256, source.summary.reports[0].size_bytes),
        _entry("normalization_summary", summaries["normalization_summary"].parent_schema,
               summaries["normalization_summary"].sha256,
               source.summary.reports[1].size_bytes),
        _entry("class_summary", summaries["class_summary"].parent_schema,
               summaries["class_summary"].sha256, source.summary.reports[2].size_bytes),
        _entry("prefilter_policy", "genosapiens.common-allele-prefilter-policy",
               manifest.prefilter_policy_sha256),
        _entry("prefilter_manifest", FILTER_MANIFEST_SCHEMA, filter_manifest.sha256,
               len(filter_manifest.as_bytes())),
        _entry("filter_summary", summaries["filter_summary"].parent_schema,
               summaries["filter_summary"].sha256,
               len(filter_manifest.summary.as_bytes())),
        _entry("working_set_manifest", WORKING_SET_MANIFEST_SCHEMA, working.sha256,
               len(working.as_bytes())),
        _entry("hpo_input", "genosapiens.hpo-input", manifest.hpo_input_sha256),
        _entry(
            "phenotype_ingestion", "genosapiens.phenotype-ingestion",
            diagnostic.phenotype_ingestion_sha256,
        ),
        _entry("ontology_manifest", "genosapiens.ontology-manifest",
               manifest.ontology_manifest_sha256),
        _entry("diagnostic_receipt", DIAGNOSTIC_RECEIPT_SCHEMA, diagnostic.sha256,
               len(diagnostic.as_bytes())),
        _entry("source_accounting", "genosapiens.wgs-source-accounting",
               manifest.source_accounting.sha256,
               len(source_accounting.as_bytes())),
        _entry("source_audit", "genosapiens.wgs-source-audit",
               manifest.source_audit.sha256, len(audit.as_bytes())),
        _entry("preparation_run_options", "genosapiens.wgs-preparation-run-options",
               manifest.preparation_run_options.sha256, len(options.as_bytes())),
        _entry("prepared_index", PREPARED_INDEX_SCHEMA, index.sha256,
               len(index.as_bytes())),
        _entry("execution_limits", "genosapiens.bounded-wgs-execution",
               manifest.execution_limits_sha256),
        _entry("execution_authority", "genosapiens.track1-base-release" if manifest.base_release_sha256
               else "genosapiens.wgs-public-validation-authority",
               manifest.base_release_sha256 or manifest.public_validation_authority_sha256),
        _entry("executable_stack", "genosapiens.executable-stack",
               manifest.executable_stack_sha256),
        _entry("runtime_image", "genosapiens.runtime-image",
               manifest.runtime_image_digest.removeprefix("sha256:")),
    ]
    entries.extend(
        _entry(value.role, value.parent_schema, value.sha256)
        for value in manifest.resources
    )
    entries.extend(
        _entry(f"journal_{value.role.removesuffix('_journal')}",
               JOURNAL_ENTRY_SCHEMA, value.sha256)
        for value in manifest.journal_heads
    )
    entries.extend(_state_entry(state) for state in index.states)
    source_roles = {"header", "source_manifest", "qc_summary", "normalization_summary",
                    "class_summary", "source_accounting", "source_audit",
                    "journal_source", "journal_source_audit", "state_source", "state_source_audit"}
    entries = [replace(entry, storage_owner="completed_source")
               if entry.role in source_roles else entry for entry in entries]
    return PreparationInventory(
        manifest.case_manifest_sha256,
        tuple(sorted(entries, key=lambda value: value.role)),
        tuple(sorted(value.role.removesuffix("_journal")
                     for value in manifest.journal_heads)),
        completed_source.descriptor.sha256, manifest.completed_source_input_path,
        str(completed_source.root.parent), str(preparation_root),
    )


def _publish_or_reuse(root: Path, name: str, content: bytes) -> tuple[Artifact, bool]:
    expected = Artifact(name, hashlib.sha256(content).hexdigest(), len(content))
    path = root / name
    if path.exists():
        if read_artifact(root, expected) != content:
            raise WgsPreparedCaseError("prepared-case artifact is immutable")
        return expected, True
    return publish_bytes(root, name, content), False


def _validate_authority(
    *,
    publication_root: Path,
    original_input_path: Path,
    source_root: Path,
    source_manifest: LedgerManifest,
    source_accounting: SourceAccountingClosure,
    filter_root: Path,
    filter_manifest: FilterManifest,
    working_root: Path,
    working_manifest: WorkingSetManifest,
    policy_sha256: str,
    population_manifest_sha256: str,
    clinical_manifest_sha256: str,
    exception_manifest_sha256: str,
    precision_manifest_sha256: str,
    scope_policy_sha256: str,
    scope_manifest_sha256: str,
    case_manifest_sha256: str,
    selected_sample_sha256: str,
    hpo_input_sha256: str,
    phenotype_ingestion_sha256: str,
    ontology_manifest_sha256: str,
    diagnostic: PreparationDiagnosticReceipt,
    execution_config: BoundedWgsExecutionConfig,
    base_release_bytes: bytes | None,
    expected_base_release_sha256: str | None,
    executable_stack_sha256: str,
    runtime_image_digest: str,
    journals: tuple[CompletedJournalBinding, ...],
    completed_source_input: CompletedSourceArtifact,
    source_origin: str,
    hpo_input_path: Path | None = None,
    public_qualification_authority_bytes: bytes | None = None,
    expected_public_qualification_authority_sha256: str | None = None,
    public_validation_authority_bytes: bytes | None = None,
    expected_public_validation_authority_sha256: str | None = None,
) -> tuple[PreparedCaseManifest, PreparationInventory, PreparedIndex, ResolvedCompletedSource]:
    for value in (
        policy_sha256,
        population_manifest_sha256,
        clinical_manifest_sha256,
        exception_manifest_sha256,
        precision_manifest_sha256,
        scope_policy_sha256, scope_manifest_sha256,
        case_manifest_sha256,
        selected_sample_sha256,
        hpo_input_sha256,
        phenotype_ingestion_sha256,
        ontology_manifest_sha256,
        executable_stack_sha256,
    ):
        if not _digest(value):
            raise WgsPreparedCaseError("prepared-case binding is invalid")
    if _IMAGE.fullmatch(runtime_image_digest) is None:
        raise WgsPreparedCaseError("prepared-case runtime image is invalid")
    _verify_source_accounting(source_manifest, source_accounting)
    preparation_root = filter_root.parent
    if working_root.parent != preparation_root:
        raise WgsPreparedCaseError("prepared-case roots are inconsistent")
    index = load_prepared_index(_read_local(
        preparation_root / "state", PREPARED_INDEX_NAME,
        diagnostic.prepared_index_sha256,
    ))
    options = load_preparation_run_options(_read_local(
        preparation_root, "run-options.json",
        index.preparation_run_options_sha256,
    ))
    if (type(completed_source_input) is not CompletedSourceArtifact
            or completed_source_input.schema != COMPLETED_SOURCE_SCHEMA
            or completed_source_input.sha256 != index.completed_source_input_sha256):
        raise WgsPreparedCaseError("prepared completed source input is invalid")
    completed_source = resolve_completed_source(load_completed_source_input(completed_source_input.read()),
        expected_source=source_manifest.source, audit_requested=options.verify_source_ledger)
    if (source_root != completed_source.root or completed_source.manifest != source_manifest
            or completed_source.accounting != source_accounting
            or source_origin not in {"produced_this_run", "completed_input"}
            or source_origin == "produced_this_run" and source_root.parent != preparation_root):
        raise WgsPreparedCaseError("prepared completed source ownership is invalid")
    producer = completed_source.descriptor.producer
    if source_origin == "produced_this_run" and (
            producer.source.executable_stack_sha256 != executable_stack_sha256
            or producer.runtime_image_digest != runtime_image_digest
            or producer.execution_config_sha256 != hashlib.sha256(
                bounded_wgs_execution_config_as_json(execution_config)).hexdigest()):
        raise WgsPreparedCaseError("new source producer differs from preparation")
    audit = completed_source.audit
    execution_limits_sha256 = hashlib.sha256(
        bounded_wgs_execution_config_as_json(execution_config)
    ).hexdigest()
    if (index != build_prepared_index(
            state_root=preparation_root / "state",
            source_manifest=source_manifest,
            filter_manifest=filter_manifest,
            working_manifest=working_manifest,
            source_accounting_sha256=source_accounting.sha256,
            source_audit_sha256=audit.sha256,
            preparation_run_options_sha256=options.sha256,
            journals=journals,
            audit_requested=options.verify_source_ledger,
            completed_source=completed_source,
            executable_stack_sha256=executable_stack_sha256,
            execution_config_sha256=execution_limits_sha256,
            runtime_image_digest=runtime_image_digest,
        ) or audit.status != ("passed" if options.verify_source_ledger else "not_requested")
            or audit.requested != options.verify_source_ledger
            or audit.accounting_sha256 != source_accounting.sha256
            or audit.source_manifest_sha256 != source_manifest.sha256
            or audit.run_options_sha256 != options.sha256
            or audit.executable_stack_sha256 != producer.source.executable_stack_sha256
            or audit.runtime_image_digest != producer.runtime_image_digest
            or audit.reference_manifest_sha256
            != source_manifest.source.reference_manifest_sha256
            or audit.input_sha256 != source_manifest.source.input_sha256):
        raise WgsPreparedCaseError("prepared run products are inconsistent")
    working_summary = working_manifest.summary
    execution_limits_sha256 = hashlib.sha256(
        bounded_wgs_execution_config_as_json(execution_config)
    ).hexdigest()
    if public_validation_authority_bytes is not None:
        from genosapiens.wgs_validation_authority import load_public_validation_authority
        if (base_release_bytes is not None or expected_base_release_sha256 is not None
                or public_qualification_authority_bytes is not None
                or expected_public_qualification_authority_sha256 is not None):
            raise WgsPreparedCaseError("mixed preparation execution authorities")
        release = load_public_validation_authority(public_validation_authority_bytes,
            expected_sha256=expected_public_validation_authority_sha256)
        if (source_origin == "completed_input" and release.completed_source_input != completed_source_input
                or source_origin == "produced_this_run" and release.completed_source_input is not None):
            raise WgsPreparedCaseError("prepared completed input differs from validation authority")
        release.require_execution(original_input_path=original_input_path,
            hpo_input_sha256=hpo_input_sha256, output_path=publication_root,
            executable_stack_sha256=executable_stack_sha256,
            runtime_image_digest=runtime_image_digest, verify_input_bytes=False,
            preparation_run_options_sha256=options.sha256)
        concordance_schema = "genosapiens.public-validation-reference-concordance"
        release_authorized = True
    else:
        if expected_public_validation_authority_sha256 is not None:
            raise WgsPreparedCaseError("public validation authority is incomplete")
        release = load_base_release(
            base_release_bytes,
            expected_sha256=expected_base_release_sha256,
        )
        if (
            (public_qualification_authority_bytes is None)
            != (expected_public_qualification_authority_sha256 is None)
        ):
            raise WgsPreparedCaseError(
                "public qualification authority binding is incomplete"
            )
        concordance_schema = "genosapiens.reference-concordance"
        if public_qualification_authority_bytes is None:
            release_authorized = True  # Version and resource closure is checked below.
        else:
            from genosapiens.wgs_public_qualification import (
                load_public_qualification_authority,
                require_public_qualification_execution,
                require_public_qualification_release,
            )

            authority = load_public_qualification_authority(
                public_qualification_authority_bytes,
                expected_sha256=expected_public_qualification_authority_sha256,
            )
            if (source_origin == "completed_input" and authority.completed_source_input != completed_source_input
                    or source_origin == "produced_this_run" and authority.completed_source_input is not None):
                raise WgsPreparedCaseError("prepared completed input differs from qualification authority")
            require_public_qualification_execution(
                authority,
                original_input_path=original_input_path,
                hpo_input_path=hpo_input_path,
                output_path=publication_root,
                input_sha256=source_manifest.source.input_sha256,
                input_size_bytes=source_manifest.source.input_size_bytes,
                hpo_input_sha256=hpo_input_sha256,
                execution_config=execution_config,
                executable_stack_sha256=executable_stack_sha256,
                runtime_image_digest=runtime_image_digest,
                policy_sha256=policy_sha256,
                population_manifest_sha256=population_manifest_sha256,
                clinical_manifest_sha256=clinical_manifest_sha256,
                exception_manifest_sha256=exception_manifest_sha256,
                precision_manifest_sha256=precision_manifest_sha256,
                preparation_run_options_sha256=options.sha256,
            )
            require_public_qualification_release(authority, release)
            if authority.invocation_id in {"complete-hg002-deepvariant-wgs", "focused-fresh-subsample"}:
                concordance_schema = "genosapiens.public-validation-reference-concordance"
            release_authorized = True
    bindings = release.policy_bindings
    filter_parents = {value.role: value for value in filter_manifest.parents}
    filter_resources = {role: value.sha256 for role, value in filter_parents.items()}
    source_parent_models = {value.role: value for value in source_manifest.parents}
    if (
        not release_authorized
        or filter_parents.get("completed_source_input") != ParentIdentity(
            "completed_source_input", COMPLETED_SOURCE_SCHEMA, 1, completed_source.descriptor.sha256)
        or bindings.common_allele_policy_sha256 != policy_sha256
        or bindings.common_allele_exception_manifest_sha256
        != exception_manifest_sha256
        or bindings.gnomad_projection_precision_sha256
        != precision_manifest_sha256
        or bindings.bounded_wgs_execution_sha256 != execution_limits_sha256
        or bindings.common_allele_public_index_manifest_sha256
        not in release.resource_manifest_sha256s
        or any(value not in release.resource_manifest_sha256s for value in (
            population_manifest_sha256,
            clinical_manifest_sha256,
            exception_manifest_sha256,
            precision_manifest_sha256,
            ontology_manifest_sha256,
        ))
        or release.reference_manifest_sha256
        != source_manifest.source.reference_manifest_sha256
        or release.executable_stack_sha256 != executable_stack_sha256
        or release.runtime_image_digest != runtime_image_digest
        or source_manifest.source.case_manifest_sha256 != case_manifest_sha256
        or source_manifest.source.selected_sample_sha256 != selected_sample_sha256
        or source_manifest.source.concordance_sha256
        != diagnostic.concordance_sha256
        or source_manifest.source.executable_stack_sha256 != producer.source.executable_stack_sha256
        or source_parent_models["case_manifest"].parent_schema
        != "genosapiens.case-manifest"
        or source_parent_models["case_manifest"].parent_version != VERSION
        or source_parent_models["complete_concordance"].parent_schema
        != concordance_schema
        or source_parent_models["complete_concordance"].parent_version != VERSION
        or source_parent_models["reference_manifest"].parent_schema
        != "genosapiens.public-resource-manifest"
        or source_parent_models["reference_manifest"].parent_version != VERSION
        or filter_resources.get("population_manifest")
        != population_manifest_sha256
        or filter_resources.get("clinical_manifest") != clinical_manifest_sha256
        or filter_resources.get("exception_manifest") != exception_manifest_sha256
        or filter_resources.get("early_scope_policy") != scope_policy_sha256
        or filter_resources.get("early_scope_manifest") != scope_manifest_sha256
        or filter_parents["early_scope_policy"].parent_schema != "genosapiens.track1-early-scope-policy"
        or filter_parents["early_scope_manifest"].parent_schema != "genosapiens.track1-early-scope-manifest"
        or filter_resources.get("precision_manifest") != precision_manifest_sha256
        or filter_parents["population_manifest"].parent_schema
        != "genosapiens.population-manifest"
        or filter_parents["clinical_manifest"].parent_schema
        != "genosapiens.clinical-manifest"
        or filter_parents["exception_manifest"].parent_schema
        != "genosapiens.exception-manifest"
        or filter_parents["precision_manifest"].parent_schema
        != "genosapiens.precision-manifest"
        or any(
            filter_parents[role].parent_version != VERSION
            for role in (
                "population_manifest", "clinical_manifest",
                "exception_manifest", "precision_manifest", "early_scope_policy", "early_scope_manifest",
            )
        )
    ):
        raise WgsPreparedCaseError("prepared-case parent closure is invalid")
    source_parents = {value.role: value.sha256 for value in source_manifest.parents}
    if (
        source_parents.get("execution_limits") != producer.execution_config_sha256
        or source_parents.get("executable_stack") != producer.source.executable_stack_sha256
        or source_parent_models["execution_limits"].parent_schema
        != execution_config.schema
        or source_parent_models["execution_limits"].parent_version
        != producer.execution_config_version
        or source_parent_models["executable_stack"].parent_schema
        != "genosapiens.executable-stack"
        or source_parent_models["executable_stack"].parent_version != VERSION
    ):
        raise WgsPreparedCaseError("prepared-case execution binding is stale")
    journal_heads = _journal_parents(
        journals,
        execution_limits_sha256=execution_limits_sha256,
        execution_limits_version=execution_config.version,
        executable_stack_sha256=executable_stack_sha256,
        runtime_image_digest=runtime_image_digest,
        completed_source=completed_source,
    )
    expected_diagnostic = PreparationDiagnosticReceipt(
        case_manifest_sha256,
        selected_sample_sha256,
        source_manifest.source.concordance_sha256,
        source_manifest.sha256,
        _source_accounting_sha256(source_accounting),
        filter_manifest.sha256,
        working_manifest.sha256,
        hpo_input_sha256,
        phenotype_ingestion_sha256,
        ontology_manifest_sha256,
        bindings.diagnostic_policy_sha256,
        "complete",
        "complete",
        "passed",
        "complete_original_source",
        "complete",
        "validated",
        (),
        audit.sha256,
        options.sha256,
        index.sha256,
        options.verify_source_ledger,
        audit.status,
        completed_source.descriptor.sha256,
    )
    if diagnostic != expected_diagnostic:
        raise WgsPreparedCaseError("prepared diagnostics are stale")
    summaries = _summary_parents(source_manifest, filter_manifest)
    resources = _resource_parents(filter_manifest)
    manifest_without_inventory = dict(
        status="passed",
        case_manifest_sha256=case_manifest_sha256,
        original_input_sha256=source_manifest.source.input_sha256,
        original_input_size_bytes=source_manifest.source.input_size_bytes,
        selected_sample_sha256=selected_sample_sha256,
        concordance_sha256=source_manifest.source.concordance_sha256,
        header=ParentIdentity("header", HEADER_SCHEMA, 1, source_manifest.header.sha256),
        source_manifest=ParentIdentity(
            "source_manifest", SOURCE_MANIFEST_SCHEMA, 1, source_manifest.sha256,
        ),
        summaries=summaries,
        prefilter_policy_sha256=policy_sha256,
        prefilter_manifest=ParentIdentity(
            "prefilter_manifest", FILTER_MANIFEST_SCHEMA, 2,
            filter_manifest.sha256,
        ),
        resources=resources,
        journal_heads=journal_heads,
        working_sets=(ParentIdentity(
            "working_set_manifest", WORKING_SET_MANIFEST_SCHEMA, 2,
            working_manifest.sha256,
        ),),
        hpo_input_sha256=hpo_input_sha256,
        ontology_manifest_sha256=ontology_manifest_sha256,
        diagnostic_receipt=ParentIdentity(
            "diagnostic_receipt", DIAGNOSTIC_RECEIPT_SCHEMA, HANDOFF_VERSION,
            diagnostic.sha256,
        ),
        source_accounting=ParentIdentity(
            "source_accounting", "genosapiens.wgs-source-accounting", 1,
            source_accounting.sha256,
        ),
        source_audit=ParentIdentity(
            "source_audit", "genosapiens.wgs-source-audit", 1, audit.sha256,
        ),
        preparation_run_options=ParentIdentity(
            "preparation_run_options", "genosapiens.wgs-preparation-run-options", 1,
            options.sha256,
        ),
        prepared_index=ParentIdentity(
            "prepared_index", PREPARED_INDEX_SCHEMA, PREPARED_INDEX_VERSION, index.sha256,
        ),
        completed_source_input=ParentIdentity("completed_source_input", COMPLETED_SOURCE_SCHEMA,
                                              1, completed_source.descriptor.sha256),
        completed_source_input_path=completed_source_input.path,
        source_origin=source_origin,
        execution_limits_sha256=execution_limits_sha256,
        base_release_sha256=expected_base_release_sha256,
        public_validation_authority_sha256=expected_public_validation_authority_sha256,
        executable_stack_sha256=executable_stack_sha256,
        runtime_image_digest=runtime_image_digest,
    )
    placeholder = ParentIdentity("inventory", INVENTORY_SCHEMA, HANDOFF_VERSION, "0" * 64)
    provisional = PreparedCaseManifest(
        **manifest_without_inventory,
        inventory=placeholder,
    )
    inventory = _build_inventory(
        provisional, source_manifest, filter_manifest, working_manifest, diagnostic,
        index, source_accounting, audit, options, completed_source, preparation_root,
    )
    manifest = PreparedCaseManifest(
        **manifest_without_inventory,
        inventory=ParentIdentity(
            "inventory", INVENTORY_SCHEMA, HANDOFF_VERSION, inventory.sha256,
        ),
    )
    if (
        working_summary.total_source_occurrences
        != source_accounting.total_physical_alternate_slots
    ):
        raise WgsPreparedCaseError("prepared occurrence accounting is incomplete")
    return manifest, inventory, index, completed_source


def publish_prepared_case(
    *,
    output_root: str | Path,
    original_input_path: str | Path,
    source_root: str | Path,
    source_manifest: LedgerManifest,
    source_accounting: SourceAccountingClosure,
    filter_root: str | Path,
    filter_manifest: FilterManifest,
    working_root: str | Path,
    working_manifest: WorkingSetManifest,
    policy_sha256: str,
    population_manifest_sha256: str,
    clinical_manifest_sha256: str,
    exception_manifest_sha256: str,
    precision_manifest_sha256: str,
    scope_policy_sha256: str,
    scope_manifest_sha256: str,
    case_manifest_sha256: str,
    selected_sample_sha256: str,
    hpo_input_sha256: str,
    phenotype_ingestion_sha256: str,
    ontology_manifest_sha256: str,
    diagnostic: PreparationDiagnosticReceipt,
    execution_config: BoundedWgsExecutionConfig,
    base_release_bytes: bytes | None,
    expected_base_release_sha256: str | None,
    executable_stack_sha256: str,
    runtime_image_digest: str,
    journals: tuple[CompletedJournalBinding, ...],
    completed_source_input: CompletedSourceArtifact,
    source_origin: str,
    hpo_input_path: str | Path | None = None,
    public_qualification_authority_bytes: bytes | None = None,
    expected_public_qualification_authority_sha256: str | None = None,
    public_validation_authority_bytes: bytes | None = None,
    expected_public_validation_authority_sha256: str | None = None,
) -> PreparedCasePublication:
    """Publish passing preparation authority after full closure replay."""

    try:
        root = Path(output_root)
        if not root.is_absolute():
            raise ValueError
        if not root.exists():
            root.mkdir(mode=0o700, parents=True)
        manifest, inventory, index, completed_source = _validate_authority(
            publication_root=root,
            original_input_path=Path(original_input_path),
            source_root=Path(source_root),
            source_manifest=source_manifest,
            source_accounting=source_accounting,
            filter_root=Path(filter_root),
            filter_manifest=filter_manifest,
            working_root=Path(working_root),
            working_manifest=working_manifest,
            policy_sha256=policy_sha256,
            population_manifest_sha256=population_manifest_sha256,
            clinical_manifest_sha256=clinical_manifest_sha256,
            exception_manifest_sha256=exception_manifest_sha256,
            precision_manifest_sha256=precision_manifest_sha256,
            scope_policy_sha256=scope_policy_sha256,
            scope_manifest_sha256=scope_manifest_sha256,
            case_manifest_sha256=case_manifest_sha256,
            selected_sample_sha256=selected_sample_sha256,
            hpo_input_sha256=hpo_input_sha256,
            phenotype_ingestion_sha256=phenotype_ingestion_sha256,
            ontology_manifest_sha256=ontology_manifest_sha256,
            diagnostic=diagnostic,
            execution_config=execution_config,
            base_release_bytes=base_release_bytes,
            expected_base_release_sha256=expected_base_release_sha256,
            executable_stack_sha256=executable_stack_sha256,
            runtime_image_digest=runtime_image_digest,
            journals=journals,
            completed_source_input=completed_source_input,
            source_origin=source_origin,
            hpo_input_path=(
                None if hpo_input_path is None else Path(hpo_input_path)
            ),
            public_qualification_authority_bytes=(
                public_qualification_authority_bytes
            ),
            expected_public_qualification_authority_sha256=(
                expected_public_qualification_authority_sha256
            ),
            public_validation_authority_bytes=public_validation_authority_bytes,
            expected_public_validation_authority_sha256=expected_public_validation_authority_sha256,
        )
        working = working_manifest.summary
        snapshot = LedgerBoundEarlySnapshot(
            manifest.sha256,
            manifest.case_manifest_sha256,
            manifest.source_manifest.sha256,
            manifest.prefilter_manifest.sha256,
            manifest.working_sets[0].sha256,
            diagnostic.sha256,
            inventory.sha256,
            manifest.hpo_input_sha256,
            manifest.ontology_manifest_sha256,
            manifest.execution_limits_sha256,
            manifest.base_release_sha256,
            manifest.executable_stack_sha256,
            manifest.runtime_image_digest,
            working.total_source_occurrences,
            working.retained_occurrences,
            working.excluded_common_occurrences,
            working.excluded_scope_occurrences, working.unresolved_scope_occurrences,
            manifest.source_accounting.sha256,
            manifest.source_audit.sha256,
            manifest.preparation_run_options.sha256,
            manifest.prepared_index.sha256,
            manifest.completed_source_input.sha256,
            public_validation_authority_sha256=manifest.public_validation_authority_sha256,
        )
        contents = (
            (DIAGNOSTIC_NAME, diagnostic.as_bytes()),
            (INVENTORY_NAME, inventory.as_bytes()),
            (PREPARED_MANIFEST_NAME, manifest.as_bytes()),
            (EARLY_SNAPSHOT_NAME, snapshot.as_bytes()),
        )
        published = tuple(_publish_or_reuse(root, name, content)
                          for name, content in contents)
        return PreparedCasePublication(
            root,
            manifest,
            diagnostic,
            inventory,
            snapshot,
            index,
            tuple(value[0] for value in published),
            all(value[1] for value in published),
            completed_source,
        )
    except WgsPreparedCaseError:
        raise
    except _ERRORS:
        raise WgsPreparedCaseError("prepared-case publication failed") from None


def _read_local(root: Path, name: str, expected_sha256: str) -> bytes:
    try:
        path = root / name
        metadata = path.lstat()
        if metadata.st_size <= 0 or metadata.st_size > 8 * 1024 * 1024:
            raise ValueError
        artifact = Artifact(name, expected_sha256, metadata.st_size)
        return read_artifact(root, artifact)
    except _ERRORS:
        raise WgsPreparedCaseError("prepared-case artifact readback failed") from None


def _decode_envelope(content: bytes, cls, *, schema: str, kind: str):
    try:
        def unique_pairs(items):
            result = {}
            for key, child in items:
                if key in result:
                    raise ValueError
                result[key] = child
            return result

        value = json.loads(content, object_pairs_hook=unique_pairs)
        if (
            type(value) is not dict
            or set(value) != {
                "classification", "data", "kind", "llm_policy", "schema", "version",
            }
            or value["classification"] != "controlled"
            or value["llm_policy"] != "no-llm"
            or value["schema"] != schema
            or value["version"] != HANDOFF_VERSION
            or value["kind"] != kind
            or type(value["data"]) is not dict
        ):
            raise ValueError
        if cls is PreparationDiagnosticReceipt:
            data = dict(value["data"])
            data["blockers"] = tuple(data["blockers"])
        elif cls is PreparationInventory:
            data = dict(value["data"])
            data["entries"] = tuple(
                PreparationInventoryEntry(**item) for item in data["entries"]
            )
            data["journal_stages"] = tuple(data["journal_stages"])
        else:
            data = value["data"]
        result = cls(**data)
        if result.as_bytes() != content:
            raise ValueError
        return result
    except _ERRORS:
        raise WgsPreparedCaseError("prepared-case envelope is invalid") from None


def load_prepared_case(
    *,
    publication_root: str | Path,
    expected_manifest_sha256: str,
    expected_early_snapshot_sha256: str,
    original_input_path: str | Path,
    source_root: str | Path,
    source_manifest: LedgerManifest,
    source_accounting: SourceAccountingClosure,
    filter_root: str | Path,
    filter_manifest: FilterManifest,
    working_root: str | Path,
    working_manifest: WorkingSetManifest,
    policy_sha256: str,
    population_manifest_sha256: str,
    clinical_manifest_sha256: str,
    exception_manifest_sha256: str,
    precision_manifest_sha256: str,
    scope_policy_sha256: str,
    scope_manifest_sha256: str,
    case_manifest_sha256: str,
    selected_sample_sha256: str,
    hpo_input_sha256: str,
    phenotype_ingestion_sha256: str,
    ontology_manifest_sha256: str,
    execution_config: BoundedWgsExecutionConfig,
    base_release_bytes: bytes | None,
    expected_base_release_sha256: str | None,
    executable_stack_sha256: str,
    runtime_image_digest: str,
    journals: tuple[CompletedJournalBinding, ...],
    hpo_input_path: str | Path | None = None,
    public_qualification_authority_bytes: bytes | None = None,
    expected_public_qualification_authority_sha256: str | None = None,
    public_validation_authority_bytes: bytes | None = None,
    expected_public_validation_authority_sha256: str | None = None,
) -> PreparedCasePublication:
    """Load preparation authority without reparsing the original VCF."""

    try:
        root = Path(publication_root)
        if (
            not root.is_absolute()
            or not _digest(expected_manifest_sha256)
            or not _digest(expected_early_snapshot_sha256)
        ):
            raise ValueError
        manifest_content = _read_local(
            root, PREPARED_MANIFEST_NAME, expected_manifest_sha256,
        )
        manifest = load_model(manifest_content, PreparedCaseManifest)
        diagnostic_content = _read_local(
            root, DIAGNOSTIC_NAME, manifest.diagnostic_receipt.sha256,
        )
        diagnostic = _decode_envelope(
            diagnostic_content,
            PreparationDiagnosticReceipt,
            schema=DIAGNOSTIC_RECEIPT_SCHEMA,
            kind="preparation_diagnostics",
        )
        inventory_content = _read_local(
            root, INVENTORY_NAME, manifest.inventory.sha256,
        )
        inventory = _decode_envelope(
            inventory_content,
            PreparationInventory,
            schema=INVENTORY_SCHEMA,
            kind="preparation_inventory",
        )
        snapshot_content = _read_local(
            root, EARLY_SNAPSHOT_NAME, expected_early_snapshot_sha256,
        )
        snapshot = _decode_envelope(
            snapshot_content,
            LedgerBoundEarlySnapshot,
            schema=EARLY_SNAPSHOT_SCHEMA,
            kind="ledger_bound_early_snapshot",
        )
        descriptor_path = Path(manifest.completed_source_input_path)
        descriptor_content = _read_local(descriptor_path.parent, descriptor_path.name,
                                         manifest.completed_source_input.sha256)
        completed_source_input = CompletedSourceArtifact(str(descriptor_path), COMPLETED_SOURCE_SCHEMA,
            1, manifest.completed_source_input.sha256, len(descriptor_content))
        expected_manifest, expected_inventory, index, completed_source = _validate_authority(
            publication_root=root,
            original_input_path=Path(original_input_path),
            source_root=Path(source_root),
            source_manifest=source_manifest,
            source_accounting=source_accounting,
            filter_root=Path(filter_root),
            filter_manifest=filter_manifest,
            working_root=Path(working_root),
            working_manifest=working_manifest,
            policy_sha256=policy_sha256,
            population_manifest_sha256=population_manifest_sha256,
            clinical_manifest_sha256=clinical_manifest_sha256,
            exception_manifest_sha256=exception_manifest_sha256,
            precision_manifest_sha256=precision_manifest_sha256,
            scope_policy_sha256=scope_policy_sha256,
            scope_manifest_sha256=scope_manifest_sha256,
            case_manifest_sha256=case_manifest_sha256,
            selected_sample_sha256=selected_sample_sha256,
            hpo_input_sha256=hpo_input_sha256,
            phenotype_ingestion_sha256=phenotype_ingestion_sha256,
            ontology_manifest_sha256=ontology_manifest_sha256,
            diagnostic=diagnostic,
            execution_config=execution_config,
            base_release_bytes=base_release_bytes,
            expected_base_release_sha256=expected_base_release_sha256,
            executable_stack_sha256=executable_stack_sha256,
            runtime_image_digest=runtime_image_digest,
            journals=journals,
            completed_source_input=completed_source_input,
            source_origin=manifest.source_origin,
            hpo_input_path=(
                None if hpo_input_path is None else Path(hpo_input_path)
            ),
            public_qualification_authority_bytes=(
                public_qualification_authority_bytes
            ),
            expected_public_qualification_authority_sha256=(
                expected_public_qualification_authority_sha256
            ),
            public_validation_authority_bytes=public_validation_authority_bytes,
            expected_public_validation_authority_sha256=expected_public_validation_authority_sha256,
        )
        expected_snapshot = LedgerBoundEarlySnapshot(
            expected_manifest.sha256,
            expected_manifest.case_manifest_sha256,
            expected_manifest.source_manifest.sha256,
            expected_manifest.prefilter_manifest.sha256,
            expected_manifest.working_sets[0].sha256,
            diagnostic.sha256,
            expected_inventory.sha256,
            expected_manifest.hpo_input_sha256,
            expected_manifest.ontology_manifest_sha256,
            expected_manifest.execution_limits_sha256,
            expected_manifest.base_release_sha256,
            expected_manifest.executable_stack_sha256,
            expected_manifest.runtime_image_digest,
            working_manifest.summary.total_source_occurrences,
            working_manifest.summary.retained_occurrences,
            working_manifest.summary.excluded_common_occurrences,
            working_manifest.summary.excluded_scope_occurrences,
            working_manifest.summary.unresolved_scope_occurrences,
            expected_manifest.source_accounting.sha256,
            expected_manifest.source_audit.sha256,
            expected_manifest.preparation_run_options.sha256,
            expected_manifest.prepared_index.sha256,
            expected_manifest.completed_source_input.sha256,
            public_validation_authority_sha256=expected_manifest.public_validation_authority_sha256,
        )
        if (
            manifest != expected_manifest
            or inventory != expected_inventory
            or snapshot != expected_snapshot
        ):
            raise ValueError
        artifacts = (
            Artifact(DIAGNOSTIC_NAME, diagnostic.sha256, len(diagnostic_content)),
            Artifact(INVENTORY_NAME, inventory.sha256, len(inventory_content)),
            Artifact(PREPARED_MANIFEST_NAME, manifest.sha256, len(manifest_content)),
            Artifact(EARLY_SNAPSHOT_NAME, snapshot.sha256, len(snapshot_content)),
        )
        return PreparedCasePublication(
            root, manifest, diagnostic, inventory, snapshot, index, artifacts, True, completed_source,
        )
    except WgsPreparedCaseError:
        raise
    except _ERRORS:
        raise WgsPreparedCaseError("prepared-case readback failed") from None
