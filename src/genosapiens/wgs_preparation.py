"""Bounded construction of the immutable whole-genome source ledger."""

from __future__ import annotations

from dataclasses import asdict, dataclass, fields, replace
from contextlib import ExitStack, closing, contextmanager, nullcontext
from concurrent.futures import ProcessPoolExecutor
from collections import deque
import hashlib
import json
import math
import multiprocessing
import os
from pathlib import Path
import pickle
import re
import sqlite3
import stat
import tempfile
from typing import Callable, Iterable, Literal, Protocol, TYPE_CHECKING

from genosapiens.bounded_wgs_config import (
    BoundedWgsExecutionConfig,
    bounded_wgs_execution_config_as_json,
)
from genosapiens.wgs_failures import SOURCE_OPERATIONS, publish_failure_evidence, static_failure
from genosapiens.governance import ArtifactClassification
from genosapiens.llm_guard import LLMPolicy
from genosapiens.ontology_bundle import ParsedOntology
from genosapiens.phenotype_ingestion import (
    ingest_phenotypes,
    phenotype_ingestion_as_json,
)
from genosapiens.reference_concordance import (
    PublicValidationConcordanceReport,
    ReferenceConcordanceReport,
    concordance_report_as_json,
)
from genosapiens.variant_inventory import (
    IncrementalVariantClassSummary,
    IncrementalVariantInventoryAccumulator,
    VariantRecordInventoryFacts,
)
from genosapiens.variant_normalization import (
    QuarantineReason,
    SourceRecordNormalization,
    SourceSmallVariantRecord,
    normalize_source_record,
    source_allele_id,
)
from genosapiens.vcf_qc import (
    QcSummaryAccumulator,
    VcfQcIncrementalSummary,
    VcfQcError,
    VcfRecordQc,
    _alternate_allele_count,
    iter_vcf_qc_records,
)
from genosapiens.wgs_ledger import (
    HEADER_SCHEMA,
    MAX_MANIFEST_BYTES,
    MAX_ROW_BYTES,
    MAX_SHARD_BYTES,
    MAX_SHARD_ROWS,
    Artifact,
    EncodedLedgerRow,
    FilterManifest,
    JournalEntry,
    LedgerManifest,
    LedgerRow,
    LedgerSummary,
    ParentIdentity,
    SourceBinding,
    Shard,
    StageCommitState,
    StageRunBinding,
    WorkingSetManifest,
    attest_existing_artifact,
    canonical,
    encode_ledger_row,
    iter_ledger,
    iter_shard,
    load_model,
    publish_bytes,
    publish_manifest,
    publish_shard,
    read_artifact,
    read_journal,
    stage_lock,
)
from genosapiens.wgs_runtime import RuntimePhaseRecorder, runtime_phase

if TYPE_CHECKING:
    from genosapiens.wgs_ledger import FilterDecision
    from genosapiens.wgs_prepared_case import PreparedCasePublication
    from genosapiens.wgs_filter import SourceDecisionCursor


NORMALIZATION_SUMMARY_SCHEMA = "genosapiens.wgs-normalization-summary"
SOURCE_PAYLOAD_SCHEMA = "genosapiens.wgs-source-payload"
from genosapiens.vcf_contracts import SOURCE_PAYLOAD_VERSION


class WgsPreparationError(RuntimeError):
    """Static bounded-preparation failure without controlled row content."""


class SourceWorkerBatchError(WgsPreparationError):
    """Static failure with a bounded source interval for local diagnostics."""

    def __init__(self, first_record: int, last_record: int):
        super().__init__("source worker batch failed")
        self.first_record = first_record
        self.last_record = last_record


def _source_call(operation, function, *args, **kwargs):
    """Retain a static failure location without changing low-level exceptions."""
    if operation not in SOURCE_OPERATIONS:
        raise ValueError("source operation is invalid")
    try:
        return function(*args, **kwargs)
    except (Exception, KeyboardInterrupt) as error:
        previous = getattr(error, "_wgs_source_operation", None)
        if type(previous) is not str or previous not in SOURCE_OPERATIONS:
            error._wgs_source_operation = operation
        raise


@dataclass(frozen=True, slots=True, repr=False)
class SourceWorkerReference:
    """Small authenticated file identity for independent spawned FASTA handles."""

    fasta_path: Path
    fai_path: Path
    fasta_identity: tuple[int, int, int, int, int]
    fai_sha256: str
    aliases: tuple[tuple[str, str], ...] = ()

    def __post_init__(self):
        if (
            not self.fasta_path.is_absolute() or not self.fai_path.is_absolute()
            or len(self.fasta_identity) != 5
            or not all(type(value) is int and value >= 0 for value in self.fasta_identity)
            or not _digest(self.fai_sha256)
            or len(self.aliases) != len(set(source for source, _ in self.aliases))
            or not all(type(source) is str and type(target) is str
                       and source and target for source, target in self.aliases)
        ):
            raise WgsPreparationError("source worker reference binding is invalid")


def source_worker_reference(
    fasta_path: str | Path, fai_path: str | Path, *, aliases=(),
) -> SourceWorkerReference:
    """Bind worker handles to already authenticated reference files."""
    fasta = Path(fasta_path)
    fai = Path(fai_path)
    try:
        fasta_stat = fasta.lstat()
        fai_stat = fai.lstat()
        if (any(parent.is_symlink() for parent in (fasta, *fasta.parents, fai, *fai.parents))
                or not stat.S_ISREG(fasta_stat.st_mode)
                or not stat.S_ISREG(fai_stat.st_mode)
                or fai_stat.st_size > 8 * 1024 * 1024):
            raise ValueError
        content = fai.read_bytes()
        if len(content) != fai_stat.st_size:
            raise ValueError
    except (OSError, ValueError):
        raise WgsPreparationError("source worker reference is unavailable") from None
    return SourceWorkerReference(
        fasta, fai,
        (fasta_stat.st_dev, fasta_stat.st_ino, fasta_stat.st_size,
         fasta_stat.st_mtime_ns, fasta_stat.st_ctime_ns),
        hashlib.sha256(content).hexdigest(), tuple(aliases),
    )


_SOURCE_WORKER_FETCH = None
_SOURCE_WORKER_STREAM = None
_SOURCE_WORKER_ERROR = False
_SOURCE_WORKER_CONTEXT = None
_SOURCE_WORKER_DEFINITIONS = None
_SOURCE_WORKER_LENGTHS = None


@dataclass(frozen=True, slots=True, repr=False)
class _SourceWorkerContext:
    input_sha256: str
    reference_identity_id: str
    reference_role: str
    selected_sample_name: str
    sample_names: tuple[str, ...]
    definitions: tuple[tuple[tuple[str, str], tuple[object, object]], ...]
    lengths: tuple[tuple[str, int], ...]
    allele_balance_screen: bool
    max_record_bytes: int


def _initialize_source_worker(reference: SourceWorkerReference,
                              context: _SourceWorkerContext):
    global _SOURCE_WORKER_FETCH, _SOURCE_WORKER_STREAM, _SOURCE_WORKER_ERROR
    global _SOURCE_WORKER_CONTEXT, _SOURCE_WORKER_DEFINITIONS, _SOURCE_WORKER_LENGTHS
    from genosapiens.reference_concordance import _OpenIndexedFasta, _parse_fai_bytes

    _SOURCE_WORKER_ERROR = False
    _SOURCE_WORKER_CONTEXT = context
    from types import MappingProxyType
    _SOURCE_WORKER_DEFINITIONS = MappingProxyType(dict(context.definitions))
    _SOURCE_WORKER_LENGTHS = MappingProxyType(dict(context.lengths))
    try:
        fasta_stat = reference.fasta_path.lstat()
        fai_stat = reference.fai_path.lstat()
        if (reference.fasta_path.is_symlink() or reference.fai_path.is_symlink()
                or not stat.S_ISREG(fasta_stat.st_mode)
                or not stat.S_ISREG(fai_stat.st_mode)
                or (fasta_stat.st_dev, fasta_stat.st_ino, fasta_stat.st_size,
                    fasta_stat.st_mtime_ns, fasta_stat.st_ctime_ns)
                != reference.fasta_identity):
            raise ValueError
        fai_content = reference.fai_path.read_bytes()
        if (len(fai_content) != fai_stat.st_size
                or len(fai_content) > 8 * 1024 * 1024
                or hashlib.sha256(fai_content).hexdigest() != reference.fai_sha256):
            raise ValueError
        descriptor = os.open(reference.fasta_path, os.O_RDONLY | os.O_NOFOLLOW)
        _SOURCE_WORKER_STREAM = os.fdopen(descriptor, "rb")
        indexed = _OpenIndexedFasta(_SOURCE_WORKER_STREAM, _parse_fai_bytes(fai_content))
        aliases = dict(reference.aliases)

        def fetch(contig, position, length):
            return indexed.fetch(aliases.get(contig, contig), position, length)

        _SOURCE_WORKER_FETCH = fetch
    except Exception:
        if _SOURCE_WORKER_STREAM is not None:
            _SOURCE_WORKER_STREAM.close()
            _SOURCE_WORKER_STREAM = None
        _SOURCE_WORKER_FETCH = None
        _SOURCE_WORKER_ERROR = True


@dataclass(frozen=True, slots=True, repr=False)
class _SourceJobBlock:
    artifact: Artifact
    rows: int


@dataclass(frozen=True, slots=True, repr=False)
class _SourceJobResult:
    first_record: int
    last_record: int
    job_name: str
    blocks: tuple[_SourceJobBlock, ...]


def _run_source_job(
    first_record: int, records: tuple[VcfRecordQc, ...], staging_root: Path,
) -> _SourceJobResult:
    """Normalize one bounded batch and leave canonical rows in private staging."""
    from genosapiens.variant_inventory import (
        AlleleBalanceScreenPolicy, classify_record_facts,
    )

    context = _SOURCE_WORKER_CONTEXT
    if (_SOURCE_WORKER_ERROR or _SOURCE_WORKER_FETCH is None
            or type(context) is not _SourceWorkerContext or not records
            or len(records) > 256 or records[0].record_number != first_record):
        raise WgsPreparationError("source worker batch is invalid")
    job_name = f"source-job-{first_record:09d}"
    job_root = staging_root / job_name
    block_stream = None
    try:
        os.mkdir(job_root, 0o700)
        blocks = []
        block_name = None
        block_digest = hashlib.sha256()
        block_size = block_rows = 0

        def close_block():
            nonlocal block_stream, block_name, block_digest, block_size, block_rows
            if block_stream is None:
                return
            block_stream.flush()
            os.fsync(block_stream.fileno())
            block_stream.close()
            blocks.append(_SourceJobBlock(
                Artifact(block_name, block_digest.hexdigest(), block_size), block_rows,
            ))
            block_stream = None

        bounded_reference = _BoundedReferenceFetcher(
            _SOURCE_WORKER_FETCH, context.max_record_bytes, first_record=first_record,
        )
        policy = AlleleBalanceScreenPolicy()
        definitions = _SOURCE_WORKER_DEFINITIONS
        lengths = _SOURCE_WORKER_LENGTHS
        for offset, record in enumerate(records):
            if record.record_number != first_record + offset:
                raise WgsPreparationError("source worker order is invalid")
            bounded_reference.begin_record(record.record_number)
            normalization = normalize_source_record(
                _source_record(record), input_sha256=context.input_sha256,
                reference_identity_id=context.reference_identity_id,
                reference_role=context.reference_role, fetch_reference=bounded_reference,
            )
            facts = classify_record_facts(
                record, input_sha256=context.input_sha256,
                sample_names=context.sample_names,
                definitions=definitions, lengths=lengths, policy=policy,
                allele_balance_screen=context.allele_balance_screen,
            )
            row = LedgerRow(
                record.record_number, record.source_record_id,
                tuple(item.source_allele_id for item in normalization.alleles),
                _alternate_allele_count(record.alternates),
                _source_payload(record, normalization, facts,
                                selected_sample_name=context.selected_sample_name),
            )
            content = encode_ledger_row(row).content
            if len(content) > MAX_ROW_BYTES:
                raise WgsPreparationError("source worker row exceeds its bounded size")
            if block_stream is None or block_size + len(content) > MAX_SHARD_BYTES:
                close_block()
                block_name = f"rows-{len(blocks) + 1:06d}.jsonl"
                descriptor = os.open(
                    job_root / block_name,
                    os.O_WRONLY | os.O_CREAT | os.O_EXCL | os.O_NOFOLLOW, 0o600,
                )
                block_stream = os.fdopen(descriptor, "wb")
                block_digest = hashlib.sha256()
                block_size = block_rows = 0
            block_stream.write(content)
            block_digest.update(content)
            block_size += len(content)
            block_rows += 1
        close_block()
        directory_fd = os.open(job_root, os.O_RDONLY | os.O_DIRECTORY | os.O_NOFOLLOW)
        try:
            os.fsync(directory_fd)
        finally:
            os.close(directory_fd)
        return _SourceJobResult(first_record, records[-1].record_number,
                                job_name, tuple(blocks))
    except (OSError, ValueError, TypeError, KeyError, IndexError):
        raise WgsPreparationError("source worker batch failed") from None
    finally:
        if block_stream is not None and not block_stream.closed:
            block_stream.close()


def _read_source_job(staging_root: Path, result: _SourceJobResult):
    """Authenticate each completed staging block before yielding its rows."""
    if not isinstance(result, _SourceJobResult):
        raise WgsPreparationError("source worker descriptor is invalid")
    job_root = staging_root / result.job_name
    try:
        directory = job_root.lstat()
        if (job_root.is_symlink() or not stat.S_ISDIR(directory.st_mode)
                or stat.S_IMODE(directory.st_mode) != 0o700
                or directory.st_uid != os.getuid()):
            raise ValueError
        record_number = result.first_record
        for block in result.blocks:
            if (block.rows <= 0 or block.rows > 256
                    or block.artifact.size_bytes > MAX_SHARD_BYTES):
                raise ValueError
            descriptor = os.open(job_root / block.artifact.name,
                                 os.O_RDONLY | os.O_NOFOLLOW)
            with os.fdopen(descriptor, "rb") as source, tempfile.TemporaryFile(mode="w+b") as spool:
                metadata = os.fstat(source.fileno())
                if (not stat.S_ISREG(metadata.st_mode)
                        or stat.S_IMODE(metadata.st_mode) != 0o600
                        or metadata.st_uid != os.getuid() or metadata.st_nlink != 1
                        or metadata.st_size != block.artifact.size_bytes):
                    raise ValueError
                digest = hashlib.sha256()
                size = 0
                while content := source.read(1024 * 1024):
                    size += len(content)
                    if size > block.artifact.size_bytes:
                        raise ValueError
                    digest.update(content)
                    spool.write(content)
                if (size != block.artifact.size_bytes
                        or digest.hexdigest() != block.artifact.sha256):
                    raise ValueError
                spool.seek(0)
                rows = 0
                while content := spool.readline(MAX_ROW_BYTES + 1):
                    if len(content) > MAX_ROW_BYTES or not content.endswith(b"\n"):
                        raise ValueError
                    row = load_model(content, LedgerRow)
                    if row.record_number != record_number:
                        raise ValueError
                    yield EncodedLedgerRow(
                        row, content, row.record_number, len(row.source_allele_ids),
                        row.biological_alternate_count,
                    )
                    rows += 1
                    record_number += 1
                if rows != block.rows:
                    raise ValueError
        if record_number != result.last_record + 1:
            raise ValueError
    except (OSError, ValueError, TypeError):
        raise WgsPreparationError("source worker staging is inconsistent") from None


def _iter_parallel_source_rows(
    stream, *, root: Path, reference: SourceWorkerReference,
    workers: int, max_inflight: int, batch_rows: int, batch_bytes: int,
    input_sha256: str, reference_identity_id: str, reference_role: str,
    selected_sample_name: str, sample_names: tuple[str, ...],
    definitions: dict, lengths: dict,
    allele_balance_screen: bool, max_record_bytes: int,
):
    """Bounded spawned work with original-order drain and disk-backed results."""
    if workers not in (1, 4) or max_inflight != 4 or batch_rows != 256 or batch_bytes != 64 * 1024 * 1024:
        raise WgsPreparationError("source worker limits are invalid")
    staging_root = root / ".staging"
    _private_directory(staging_root)
    # Interrupted workers can leave partial blocks, including directories in
    # the old flat layout. Never reuse those unauthenticated results or collide
    # with their deterministic job names. Committed progress lives in SQLite
    # and immutable shards; retain unfinished staging for local inspection.
    staging_root = Path(tempfile.mkdtemp(prefix="source-run-", dir=staging_root))
    pending = deque()
    current = []
    current_bytes = 0
    context = _SourceWorkerContext(
        input_sha256, reference_identity_id, reference_role,
        selected_sample_name, sample_names,
        tuple(sorted(definitions.items())), tuple(sorted(lengths.items())),
        allele_balance_screen, max_record_bytes,
    )
    executor = ProcessPoolExecutor(
        max_workers=workers, mp_context=multiprocessing.get_context("spawn"),
        initializer=_initialize_source_worker, initargs=(reference, context),
    )

    def submit():
        nonlocal current, current_bytes
        if not current:
            return
        records = tuple(current)
        if len(pickle.dumps(records, protocol=5)) > batch_bytes:
            raise WgsPreparationError("source worker batch exceeds transport bound")
        first = records[0].record_number
        future = executor.submit(
            _run_source_job, first, records, staging_root,
        )
        pending.append((records, future))
        current = []
        current_bytes = 0

    def drain():
        records, future = pending.popleft()
        try:
            result = future.result()
        except Exception:
            error = SourceWorkerBatchError(
                records[0].record_number, records[-1].record_number,
            )
            error._wgs_source_operation = "source_worker_execute"
            raise error from None
        if (not isinstance(result, _SourceJobResult)
                or result.first_record != records[0].record_number
                or result.last_record != records[-1].record_number
                or result.job_name != f"source-job-{result.first_record:09d}"
                or sum(block.rows for block in result.blocks) != len(records)):
            raise WgsPreparationError("source worker descriptor is inconsistent")
        staged = _read_source_job(staging_root, result)
        for record in records:
            try:
                encoded = next(staged)
            except StopIteration:
                raise WgsPreparationError("source worker row coverage is invalid") from None
            if (encoded.record_number != record.record_number
                    or encoded.row.source_record_id != record.source_record_id
                    or encoded.row.source_allele_ids != tuple(
                        source_allele_id(record.source_record_id, number)
                        for number in range(1, len(record.alternates) + 1)
                    )
                    or encoded.biological_alternates != _alternate_allele_count(record.alternates)):
                raise WgsPreparationError("source worker source tuple is invalid")
            yield record, encoded
        try:
            next(staged)
        except StopIteration:
            pass
        else:
            raise WgsPreparationError("source worker row coverage is invalid")
        job_root = staging_root / result.job_name
        for block in result.blocks:
            (job_root / block.artifact.name).unlink()
        job_root.rmdir()

    failed = True
    try:
        for record in stream:
            try:
                transport = len(pickle.dumps(record, protocol=5)) + 64
            except (pickle.PicklingError, TypeError, ValueError):
                raise WgsPreparationError("source worker transport is invalid") from None
            if transport > batch_bytes:
                raise WgsPreparationError("source worker record exceeds transport bound")
            if current and (len(current) >= batch_rows or current_bytes + transport > batch_bytes):
                if len(pending) >= max_inflight:
                    yield from drain()
                submit()
            current.append(record)
            current_bytes += transport
        if current:
            if len(pending) >= max_inflight:
                yield from drain()
            submit()
        while pending:
            yield from drain()
        failed = False
    finally:
        if failed:
            for _records, future in pending:
                try:
                    future.cancel()
                except Exception:
                    pass
            processes = tuple((getattr(executor, "_processes", None) or {}).values())
            try:
                executor.shutdown(wait=False, cancel_futures=True)
            except Exception:
                pass
            for process in processes:
                try:
                    if process.is_alive():
                        process.terminate()
                        process.join(timeout=2)
                        if process.is_alive():
                            process.kill()
                            process.join(timeout=2)
                except Exception:
                    pass
        else:
            executor.shutdown(wait=True)
            staging_root.rmdir()


def _read_bounded_manifest(path: Path) -> bytes:
    descriptor = -1
    try:
        descriptor = os.open(path, os.O_RDONLY | os.O_NOFOLLOW)
        metadata = os.fstat(descriptor)
        if (
            not stat.S_ISREG(metadata.st_mode)
            or metadata.st_uid != os.geteuid()
            or stat.S_IMODE(metadata.st_mode) != 0o600
            or metadata.st_size <= 0
            or metadata.st_size > MAX_MANIFEST_BYTES
        ):
            raise WgsPreparationError("completed manifest is not a private bounded file")
        with os.fdopen(descriptor, "rb", closefd=True) as stream:
            descriptor = -1
            content = stream.read(MAX_MANIFEST_BYTES + 1)
        if len(content) != metadata.st_size:
            raise WgsPreparationError("completed manifest changed while it was read")
        return content
    except OSError:
        raise WgsPreparationError("completed manifest cannot be read safely") from None
    finally:
        if descriptor >= 0:
            os.close(descriptor)


def _digest(value: object) -> bool:
    return (
        isinstance(value, str)
        and len(value) == 64
        and all(character in "0123456789abcdef" for character in value)
    )


def _count(value: object, minimum: int = 0) -> bool:
    return isinstance(value, int) and not isinstance(value, bool) and value >= minimum


@dataclass(frozen=True, slots=True, repr=False)
class NormalizationReasonCount:
    reason: str
    count: int

    def __post_init__(self) -> None:
        if (
            not isinstance(self.reason, str)
            or not self.reason
            or any(character in self.reason for character in "\t\r\n")
            or not _count(self.count, 1)
        ):
            raise WgsPreparationError("normalization summary is invalid")


@dataclass(frozen=True, slots=True, repr=False)
class WgsNormalizationSummary:
    input_sha256: str
    reference_concordance_sha256: str
    total_source_records: int
    total_physical_alternate_slots: int
    normalized_occurrences: int
    unique_normalized_alleles: int
    quarantined_occurrences: int
    reference_mismatch_records: int
    quarantine_reason_counts: tuple[NormalizationReasonCount, ...]
    classification: ArtifactClassification = ArtifactClassification.CONTROLLED
    llm_policy: LLMPolicy = LLMPolicy.NO_LLM
    schema: str = NORMALIZATION_SUMMARY_SCHEMA

    def __post_init__(self) -> None:
        reasons = tuple(item.reason for item in self.quarantine_reason_counts)
        if (
            self.schema != NORMALIZATION_SUMMARY_SCHEMA
            or not _digest(self.input_sha256)
            or not _digest(self.reference_concordance_sha256)
            or not _count(self.total_source_records, 1)
            or not _count(self.total_physical_alternate_slots, self.total_source_records)
            or not _count(self.normalized_occurrences)
            or not _count(self.unique_normalized_alleles)
            or self.unique_normalized_alleles > self.normalized_occurrences
            or not _count(self.quarantined_occurrences)
            or self.normalized_occurrences + self.quarantined_occurrences
            != self.total_physical_alternate_slots
            or not _count(self.reference_mismatch_records)
            or self.reference_mismatch_records > self.total_source_records
            or reasons != tuple(sorted(set(reasons)))
            or sum(item.count for item in self.quarantine_reason_counts)
            != self.quarantined_occurrences
            or any(item.reason not in {reason.value for reason in QuarantineReason}
                   for item in self.quarantine_reason_counts)
            or self.classification is not ArtifactClassification.CONTROLLED
            or self.llm_policy is not LLMPolicy.NO_LLM
        ):
            raise WgsPreparationError("normalization summary is invalid")


@dataclass(frozen=True, slots=True, repr=False)
class WgsSourceBuildResult:
    manifest: LedgerManifest
    manifest_artifact: Artifact
    qc_summary: VcfQcIncrementalSummary
    normalization_summary: WgsNormalizationSummary
    class_summary: IncrementalVariantClassSummary


@dataclass(frozen=True, slots=True, repr=False)
class PreparationRunOptions:
    """The explicit, immutable replay choice for one preparation run."""

    verify_source_ledger: bool = False
    schema: Literal["genosapiens.wgs-preparation-run-options"] = (
        "genosapiens.wgs-preparation-run-options"
    )
    version: Literal[1] = 1
    classification: Literal["controlled"] = "controlled"
    llm_policy: Literal["no-llm"] = "no-llm"

    def __post_init__(self) -> None:
        if (type(self.verify_source_ledger) is not bool
                or self.schema != "genosapiens.wgs-preparation-run-options"
                or type(self.version) is not int or self.version != 1
                or self.classification != "controlled"
                or self.llm_policy != "no-llm"):
            raise WgsPreparationError("preparation run options are invalid")

    def as_bytes(self) -> bytes:
        return canonical(asdict(self))

    @property
    def sha256(self) -> str:
        return hashlib.sha256(self.as_bytes()).hexdigest()


def _load_source_contract(content: bytes, cls):
    """Reject obsolete, incomplete, duplicate-key and noncanonical envelopes."""
    if type(content) is not bytes or len(content) > 64 * 1024:
        raise WgsPreparationError("source contract is invalid")

    def unique(pairs):
        value = {}
        for key, item in pairs:
            if key in value:
                raise ValueError
            value[key] = item
        return value

    try:
        value = json.loads(content, object_pairs_hook=unique,
                           parse_constant=lambda _: (_ for _ in ()).throw(ValueError()))
        if type(value) is not dict or set(value) != {item.name for item in fields(cls)}:
            raise ValueError
        result = cls(**value)
        if result.as_bytes() != content:
            raise ValueError
        return result
    except (ValueError, TypeError, KeyError, OverflowError):
        raise WgsPreparationError("source contract is invalid") from None


def load_preparation_run_options(content: bytes) -> PreparationRunOptions:
    return _load_source_contract(content, PreparationRunOptions)


@dataclass(frozen=True, slots=True, repr=False)
class SourceAccountingClosure:
    """Mandatory producer-side counts, without an independent replay claim."""

    source_manifest_sha256: str
    input_sha256: str
    header_sha256: str
    qc_summary_sha256: str
    normalization_summary_sha256: str
    class_summary_sha256: str
    total_records: int
    total_physical_alternate_slots: int
    total_biological_alternates: int
    accounting_status: Literal["passed"] = "passed"
    schema: Literal["genosapiens.wgs-source-accounting"] = "genosapiens.wgs-source-accounting"
    version: Literal[1] = 1
    classification: Literal["controlled"] = "controlled"
    llm_policy: Literal["no-llm"] = "no-llm"

    def __post_init__(self) -> None:
        if (any(not _digest(value) for value in (
                self.source_manifest_sha256, self.input_sha256, self.header_sha256,
                self.qc_summary_sha256, self.normalization_summary_sha256,
                self.class_summary_sha256))
                or not _count(self.total_records, 1)
                or not _count(self.total_physical_alternate_slots, self.total_records)
                or not _count(self.total_biological_alternates)
                or self.total_biological_alternates > self.total_physical_alternate_slots
                or self.accounting_status != "passed"
                or self.schema != "genosapiens.wgs-source-accounting"
                or type(self.version) is not int or self.version != 1
                or self.classification != "controlled"
                or self.llm_policy != "no-llm"):
            raise WgsPreparationError("source accounting closure is invalid")

    def as_bytes(self) -> bytes:
        return canonical(asdict(self))

    @property
    def sha256(self) -> str:
        return hashlib.sha256(self.as_bytes()).hexdigest()


def load_source_accounting(content: bytes) -> SourceAccountingClosure:
    return _load_source_contract(content, SourceAccountingClosure)


@dataclass(frozen=True, slots=True, repr=False)
class SourceAuditResult:
    requested: bool
    status: Literal["not_requested", "passed", "failed"]
    accounting_sha256: str
    source_manifest_sha256: str
    run_options_sha256: str
    executable_stack_sha256: str
    runtime_image_digest: str
    reference_manifest_sha256: str
    input_sha256: str
    failure_category: str | None = None
    schema: Literal["genosapiens.wgs-source-audit"] = "genosapiens.wgs-source-audit"
    version: Literal[1] = 1
    classification: Literal["controlled"] = "controlled"
    llm_policy: Literal["no-llm"] = "no-llm"

    def __post_init__(self) -> None:
        if (type(self.requested) is not bool
                or self.status not in {"not_requested", "passed", "failed"}
                or (self.status == "not_requested") != (not self.requested)
                or (self.status == "failed") != (self.failure_category is not None)
                or self.failure_category is not None
                and self.failure_category not in {"semantic_mismatch", "invalid_source", "unexpected"}
                or any(not _digest(value) for value in (
                    self.accounting_sha256, self.source_manifest_sha256,
                    self.run_options_sha256, self.executable_stack_sha256,
                    self.reference_manifest_sha256, self.input_sha256))
                or not isinstance(self.runtime_image_digest, str)
                or not self.runtime_image_digest.startswith("sha256:")
                or not _digest(self.runtime_image_digest.removeprefix("sha256:"))
                or self.schema != "genosapiens.wgs-source-audit"
                or type(self.version) is not int or self.version != 1
                or self.classification != "controlled"
                or self.llm_policy != "no-llm"):
            raise WgsPreparationError("source audit result is invalid")

    def as_bytes(self) -> bytes:
        return canonical(asdict(self))

    @property
    def sha256(self) -> str:
        return hashlib.sha256(self.as_bytes()).hexdigest()


def load_source_audit(content: bytes) -> SourceAuditResult:
    return _load_source_contract(content, SourceAuditResult)


COMPLETED_SOURCE_SCHEMA = "genosapiens.wgs-completed-source-input"


def _completed_path(value: str) -> Path:
    if (type(value) is not str or not value or "\x00" in value
            or not Path(value).is_absolute() or str(Path(value)) != value
            or ".." in Path(value).parts):
        raise WgsPreparationError("completed source path is invalid")
    return Path(value)


@dataclass(frozen=True, slots=True, repr=False)
class CompletedSourceArtifact:
    path: str
    schema: str
    version: int
    sha256: str
    size_bytes: int

    def __post_init__(self):
        _completed_path(self.path)
        if (type(self.schema) is not str or not self.schema.startswith("genosapiens.")
                or type(self.version) is not int or self.version != 1
                or not _digest(self.sha256) or not _count(self.size_bytes, 1)
                or self.size_bytes > MAX_MANIFEST_BYTES):
            raise WgsPreparationError("completed source artifact is invalid")

    def read(self) -> bytes:
        path = _completed_path(self.path)
        return read_artifact(path.parent, Artifact(path.name, self.sha256, self.size_bytes))


@dataclass(frozen=True, slots=True, repr=False)
class CompletedSourceState:
    """Trusted completed SQLite metadata, without a database-page digest."""

    stage: Literal["source", "source_audit"]
    relative_path: str
    schema: str
    user_version: int
    size_bytes: int
    file_identity: tuple[int, int, int, int, int]
    run_binding_sha256: str
    committed_ordinal: int
    committed_cursor: tuple[int, int]
    committed_counts: tuple[tuple[str, int], ...]
    completed_event_ordinal: int
    completed_event_sha256: str
    journal_head_sha256: str

    def __post_init__(self):
        names = {"source": "source.sqlite", "source_audit": "source-audit.sqlite"}
        if (self.stage not in names or self.relative_path != names[self.stage]
                or self.schema != "genosapiens.wgs-stage-checkpoint"
                or type(self.user_version) is not int or self.user_version != 1
                or not _count(self.size_bytes, 1)
                or type(self.file_identity) is not tuple or len(self.file_identity) != 5
                or not all(_count(value) for value in self.file_identity)
                or self.file_identity[2] != self.size_bytes
                or not all(_digest(value) for value in (
                    self.run_binding_sha256, self.completed_event_sha256,
                    self.journal_head_sha256))
                or not _count(self.committed_ordinal, 1)
                or not _count(self.completed_event_ordinal, 1)
                or self.completed_event_ordinal != self.committed_ordinal + 1
                or type(self.committed_cursor) is not tuple or len(self.committed_cursor) != 2
                or not all(_count(value, 1) for value in self.committed_cursor)
                or type(self.committed_counts) is not tuple
                or self.committed_counts != tuple(sorted(self.committed_counts))
                or len(dict(self.committed_counts)) != len(self.committed_counts)
                or set(dict(self.committed_counts)) != {
                    "records", "physical_alternates", "biological_alternates"}
                or any(not _count(value) for _, value in self.committed_counts)):
            raise WgsPreparationError("completed source state is invalid")


@dataclass(frozen=True, slots=True, repr=False)
class CompletedSourceProducer:
    candidate_sha256: str | None
    source: SourceBinding
    runtime_image_digest: str
    execution_config_sha256: str
    execution_config_version: int
    run_options_sha256: str
    source_stage_parent_sha256: str
    source_audit_stage_parent_sha256: str
    source_journal_head_sha256: str
    source_audit_journal_head_sha256: str

    def __post_init__(self):
        if (type(self.source) is not SourceBinding
                or self.candidate_sha256 is not None and not _digest(self.candidate_sha256)
                or type(self.runtime_image_digest) is not str
                or not re.fullmatch(r"sha256:[0-9a-f]{64}", self.runtime_image_digest)
                or not _count(self.execution_config_version, 1)
                or not all(_digest(value) for value in (
                    self.execution_config_sha256, self.run_options_sha256,
                    self.source_stage_parent_sha256, self.source_audit_stage_parent_sha256,
                    self.source_journal_head_sha256, self.source_audit_journal_head_sha256))):
            raise WgsPreparationError("completed source producer is invalid")


@dataclass(frozen=True, slots=True, repr=False)
class CompletedSourceInput:
    source_root: str
    source_state_root: str
    source_journal_root: str
    source_audit_journal_root: str
    source_manifest: CompletedSourceArtifact
    source_accounting: CompletedSourceArtifact
    source_audit: CompletedSourceArtifact
    source_run_options: CompletedSourceArtifact
    source_stage_output: CompletedSourceArtifact
    source_audit_stage_output: CompletedSourceArtifact
    producer: CompletedSourceProducer
    source_state: CompletedSourceState
    source_audit_state: CompletedSourceState | None
    total_records: int
    total_physical_alternate_slots: int
    total_biological_alternates: int
    audit_requested: bool
    audit_status: Literal["passed", "not_requested"]
    source_runtime_checkpoint_path: str | None = None
    schema: Literal["genosapiens.wgs-completed-source-input"] = COMPLETED_SOURCE_SCHEMA
    version: Literal[1] = 1
    classification: Literal["controlled"] = "controlled"
    llm_policy: Literal["no-llm"] = "no-llm"

    def __post_init__(self):
        source = _completed_path(self.source_root)
        state = _completed_path(self.source_state_root)
        source_journal = _completed_path(self.source_journal_root)
        audit_journal = _completed_path(self.source_audit_journal_root)
        if (self.source_runtime_checkpoint_path is not None
                and _completed_path(self.source_runtime_checkpoint_path).name
                != ".runtime-phases-preparation.json"):
            raise WgsPreparationError("completed source runtime checkpoint path is invalid")
        if (source.name != "source" or state != source.parent / "state"
                or source_journal != source.parent / "journals/source"
                or audit_journal != source.parent / "journals/source_audit"
                or type(self.producer) is not CompletedSourceProducer
                or type(self.source_state) is not CompletedSourceState
                or self.source_state.stage != "source"
                or type(self.audit_requested) is not bool
                or self.audit_status != ("passed" if self.audit_requested else "not_requested")
                or self.audit_requested != (self.source_audit_state is not None)
                or self.source_audit_state is not None and (
                    type(self.source_audit_state) is not CompletedSourceState
                    or self.source_audit_state.stage != "source_audit")
                or not _count(self.total_records, 1)
                or not _count(self.total_physical_alternate_slots, self.total_records)
                or not _count(self.total_biological_alternates)
                or self.total_biological_alternates > self.total_physical_alternate_slots
                or self.schema != COMPLETED_SOURCE_SCHEMA
                or type(self.version) is not int or self.version != 1
                or self.classification != "controlled" or self.llm_policy != "no-llm"):
            raise WgsPreparationError("completed source input is invalid")
        for artifact, path, schema in (
            (self.source_manifest, source / "manifest.json", LedgerManifest.schema),
            (self.source_accounting, source / "accounting.json", "genosapiens.wgs-source-accounting"),
            (self.source_audit, source / "audit.json", "genosapiens.wgs-source-audit"),
            (self.source_run_options, source.parent / "run-options.json", "genosapiens.wgs-preparation-run-options"),
            (self.source_stage_output, source_journal / "stage-output.json", "genosapiens.wgs-preparation-stage-output"),
            (self.source_audit_stage_output, audit_journal / "stage-output.json", "genosapiens.wgs-preparation-stage-output"),
        ):
            if (type(artifact) is not CompletedSourceArtifact
                    or artifact.path != str(path) or artifact.schema != schema):
                raise WgsPreparationError("completed source artifact ownership is invalid")

    def as_bytes(self) -> bytes:
        return canonical(asdict(self))

    @property
    def sha256(self) -> str:
        return hashlib.sha256(self.as_bytes()).hexdigest()


@dataclass(frozen=True, slots=True, repr=False)
class ResolvedCompletedSource:
    descriptor: CompletedSourceInput
    manifest: LedgerManifest
    accounting: SourceAccountingClosure
    audit: SourceAuditResult
    run_options: PreparationRunOptions
    journals: tuple

    @property
    def root(self) -> Path:
        return Path(self.descriptor.source_root)


def load_completed_source_input(content: bytes) -> CompletedSourceInput:
    """Decode the sole contract without opening any declared dependency."""
    try:
        from genosapiens.wgs_ledger import _json
        value = _json(content, 64 * 1024)
        if type(value) is not dict or set(value) != {item.name for item in fields(CompletedSourceInput)}:
            raise ValueError
        for name in ("source_manifest", "source_accounting", "source_audit",
                     "source_run_options", "source_stage_output", "source_audit_stage_output"):
            value[name] = CompletedSourceArtifact(**value[name])
        producer = value["producer"]
        producer["source"] = SourceBinding(**producer["source"])
        value["producer"] = CompletedSourceProducer(**producer)
        for name in ("source_state", "source_audit_state"):
            item = value[name]
            if item is not None:
                item["file_identity"] = tuple(item["file_identity"])
                item["committed_cursor"] = tuple(item["committed_cursor"])
                item["committed_counts"] = tuple(tuple(pair) for pair in item["committed_counts"])
                value[name] = CompletedSourceState(**item)
        result = CompletedSourceInput(**value)
        if result.as_bytes() != content:
            raise ValueError
        return result
    except (ValueError, TypeError, KeyError, AttributeError, OverflowError):
        raise WgsPreparationError("completed source input is invalid") from None


def publish_completed_source_input(root: Path, value: CompletedSourceInput) -> Artifact:
    """The descriptor identity never includes its destination or its consumer."""
    if type(value) is not CompletedSourceInput:
        raise WgsPreparationError("completed source input is invalid")
    return _publish_or_reuse_small(root, "completed-source.json", value.as_bytes())


def completed_source_readonly_paths(binding: CompletedSourceArtifact) -> tuple[str, ...]:
    """Exact mount/retention dependencies, with no access to predecessor filter state."""
    value = load_completed_source_input(binding.read())
    return tuple(sorted({binding.path, value.source_root, value.source_state_root,
        value.source_journal_root, value.source_audit_journal_root, value.source_run_options.path,
        *((value.source_runtime_checkpoint_path,) if value.source_runtime_checkpoint_path is not None else ())}))


def load_preparation_source_binding(workspace: Path, *, expected_sha256: str) -> CompletedSourceArtifact:
    from genosapiens.wgs_ledger import _json
    content = _read_bounded_manifest(workspace / "source-input-binding.json")
    value = CompletedSourceArtifact(**_json(content))
    if (value.schema != COMPLETED_SOURCE_SCHEMA or value.sha256 != expected_sha256
            or canonical(asdict(value)) != content):
        raise WgsPreparationError("preparation source binding changed")
    value.read()
    return value


def require_completed_source_readonly(binding: CompletedSourceArtifact) -> None:
    from genosapiens.vep_annotation import _readonly
    try:
        for path in completed_source_readonly_paths(binding):
            _readonly(Path(path))
    except (OSError, ValueError, RuntimeError):
        raise WgsPreparationError("completed source dependencies require read-only mounts") from None


def _completed_source_state(
    source_root: Path, state_root: Path, stage: str, *, manifest: LedgerManifest,
    audit: SourceAuditResult, stage_parent: str, outer_head: str,
) -> tuple[CompletedSourceState, StageRunBinding]:
    """Read indexed completion metadata only, with no recovery or page scan."""
    from genosapiens.wgs_ledger import _directory, _json

    name = {"source": "source.sqlite", "source_audit": "source-audit.sqlite"}[stage]
    path = state_root / name
    connection = None
    try:
        with _directory(state_root):
            pass
        metadata = path.lstat()
        identity = (metadata.st_dev, metadata.st_ino, metadata.st_size,
                    metadata.st_mtime_ns, metadata.st_ctime_ns)
        if (not stat.S_ISREG(metadata.st_mode) or metadata.st_uid != os.getuid()
                or stat.S_IMODE(metadata.st_mode) != 0o600 or metadata.st_nlink != 1):
            raise ValueError
        connection = sqlite3.connect(path.as_uri() + "?mode=ro&immutable=1", uri=True)
        connection.execute("PRAGMA query_only=ON")
        connection.execute("PRAGMA trusted_schema=OFF")
        connection.execute("PRAGMA cache_size=-65536")
        connection.execute("PRAGMA mmap_size=0")
        if (connection.execute("PRAGMA application_id").fetchone() != (1196972365,)
                or connection.execute("PRAGMA user_version").fetchone() != (1,)):
            raise ValueError
        row = connection.execute("SELECT canonical_bytes FROM run_binding WHERE singleton=1").fetchone()
        binding = load_model(row[0], StageRunBinding)
        parents = {parent.role: parent for parent in manifest.parents}
        if (binding.stage != stage or binding.stage_parent_sha256 != stage_parent
                or binding.input_sha256 != manifest.source.input_sha256
                or binding.reference_sha256 != manifest.source.reference_manifest_sha256
                or binding.executable_stack_sha256 != manifest.source.executable_stack_sha256
                or binding.execution_config_sha256 != parents["execution_limits"].sha256
                or binding.runtime_image_digest != audit.runtime_image_digest
                or binding.run_options_sha256 != audit.run_options_sha256):
            raise ValueError
        latest = connection.execute(
            "SELECT ordinal,event_bytes,event_sha256 FROM events ORDER BY ordinal DESC LIMIT 1"
        ).fetchone()
        event = load_model(latest[1], JournalEntry)
        if (event.event_type != "stage_completed" or event.ordinal != latest[0]
                or event.sha256 != latest[2] or event.stage != stage
                or event.stage_parent_sha256 != stage_parent or event.parents != binding.parents
                or len(event.artifacts) != 1):
            raise ValueError
        final_artifact = (Artifact("manifest.json", manifest.sha256, len(manifest.as_bytes()))
                          if stage == "source" else Artifact(
                              "audit-closure.json", audit.accounting_sha256,
                              len(read_artifact(source_root, event.artifacts[0]))))
        if event.artifacts != (final_artifact,):
            raise ValueError
        head = _json(_read_bounded_manifest(source_root / f"{stage}-head.json"), 4096)
        if head != {"schema": "genosapiens.wgs-stage-journal-head", "version": 1,
                    "stage": stage, "ordinal": event.ordinal, "sha256": event.sha256}:
            raise ValueError
        if read_artifact(source_root, Artifact(
                f"{stage}-{event.ordinal:06d}.journal.json", event.sha256, len(latest[1]))) != latest[1]:
            raise ValueError
        commit = connection.execute(
            "SELECT ordinal,cursor_bytes,counts_bytes FROM commits ORDER BY ordinal DESC LIMIT 1"
        ).fetchone()
        cursor = tuple(_json(commit[1], 4096))
        counts = connection.execute(
            "SELECT records,physical,biological FROM qc_totals WHERE singleton=1"
        ).fetchone()
        expected_counts = (manifest.summary.total_records,
                           manifest.summary.total_physical_alternate_slots,
                           manifest.summary.total_biological_alternates)
        if (counts != expected_counts or cursor != (counts[0], len(manifest.shards))
                or latest[0] != commit[0] + 1):
            raise ValueError
        final_counts = _json(commit[2], 4096)
        last = manifest.shards[-1]
        if final_counts != ({"records": last.rows,
                             "physical_alternates": last.physical_alternate_slots,
                             "biological_alternates": last.biological_alternates}
                            if stage == "source" else {"records": last.rows}):
            raise ValueError
        if stage == "source":
            saved = connection.execute("SELECT ordinal,descriptor FROM source_shards ORDER BY ordinal DESC LIMIT 1").fetchone()
            if saved != (last.ordinal, last.as_bytes()):
                raise ValueError
        after = path.lstat()
        if identity != (after.st_dev, after.st_ino, after.st_size, after.st_mtime_ns, after.st_ctime_ns):
            raise ValueError
        return CompletedSourceState(stage, name, "genosapiens.wgs-stage-checkpoint", 1,
            metadata.st_size, identity, binding.sha256, commit[0], cursor,
            tuple(sorted(zip(("records", "physical_alternates", "biological_alternates"), counts))),
            event.ordinal, event.sha256, outer_head), binding
    except (OSError, sqlite3.Error, ValueError, TypeError, IndexError, KeyError):
        raise WgsPreparationError("completed source state is unavailable or stale") from None
    finally:
        if connection is not None:
            connection.close()


def _describe_completed_source(
    source_root: Path, *, candidate_sha256: str | None,
    journals: tuple,
    runtime_checkpoint_path: Path | None = None,
) -> CompletedSourceInput:
    """Capture immutable metadata after the two source stages close independently."""
    def artifact(path, schema):
        content = _read_bounded_manifest(path)
        value = CompletedSourceArtifact(str(path), schema, 1, hashlib.sha256(content).hexdigest(), len(content))
        if value.read() != content:
            raise WgsPreparationError("completed source artifact changed")
        return value, content

    by_stage = {journal.stage: journal for journal in journals}
    if len(by_stage) != len(journals) or set(by_stage) != {"source", "source_audit"}:
        raise WgsPreparationError("completed source journals are invalid")
    manifest_artifact, content = artifact(source_root / "manifest.json", LedgerManifest.schema)
    manifest = load_model(content, LedgerManifest)
    accounting_artifact, _ = artifact(source_root / "accounting.json", "genosapiens.wgs-source-accounting")
    audit_artifact, content = artifact(source_root / "audit.json", "genosapiens.wgs-source-audit")
    audit = load_source_audit(content)
    options_artifact, _ = artifact(source_root.parent / "run-options.json", "genosapiens.wgs-preparation-run-options")
    source_output, _ = artifact(by_stage["source"].root / "stage-output.json", "genosapiens.wgs-preparation-stage-output")
    audit_output, _ = artifact(by_stage["source_audit"].root / "stage-output.json", "genosapiens.wgs-preparation-stage-output")
    states = {}
    for stage in ("source", "source_audit") if audit.requested else ("source",):
        states[stage], _ = _completed_source_state(source_root, source_root.parent / "state", stage,
            manifest=manifest, audit=audit, stage_parent=by_stage[stage].stage_parent_sha256,
            outer_head=by_stage[stage].head_sha256)
    config = next(parent for parent in manifest.parents if parent.role == "execution_limits")
    producer = CompletedSourceProducer(candidate_sha256, manifest.source, audit.runtime_image_digest,
        config.sha256, config.parent_version, audit.run_options_sha256,
        by_stage["source"].stage_parent_sha256, by_stage["source_audit"].stage_parent_sha256,
        by_stage["source"].head_sha256, by_stage["source_audit"].head_sha256)
    value = CompletedSourceInput(str(source_root), str(source_root.parent / "state"),
        str(by_stage["source"].root), str(by_stage["source_audit"].root),
        manifest_artifact, accounting_artifact, audit_artifact, options_artifact,
        source_output, audit_output, producer, states["source"], states.get("source_audit"),
        manifest.summary.total_records, manifest.summary.total_physical_alternate_slots,
        manifest.summary.total_biological_alternates, audit.requested, audit.status,
        source_runtime_checkpoint_path=None if runtime_checkpoint_path is None else str(runtime_checkpoint_path))
    resolve_completed_source(value, expected_source=manifest.source, audit_requested=audit.requested)
    return value


def resolve_completed_source(
    value: CompletedSourceInput, *, expected_source: SourceBinding, audit_requested: bool,
) -> ResolvedCompletedSource:
    """Validate original completion read-only; never rebind, replay or rebuild it."""
    from genosapiens.wgs_ledger import _directory, _json
    from genosapiens.wgs_prepared_case import CompletedJournalBinding

    try:
        if (type(value) is not CompletedSourceInput or type(expected_source) is not SourceBinding
                or type(audit_requested) is not bool or value.audit_requested != audit_requested):
            raise ValueError
        for root in (value.source_root, value.source_state_root,
                     value.source_journal_root, value.source_audit_journal_root):
            with _directory(root):
                pass
        source = value.producer.source
        if any(getattr(source, field.name) != getattr(expected_source, field.name)
               for field in fields(SourceBinding) if field.name != "executable_stack_sha256"):
            raise ValueError
        manifest = load_model(value.source_manifest.read(), LedgerManifest)
        accounting = load_source_accounting(value.source_accounting.read())
        audit = load_source_audit(value.source_audit.read())
        options = load_preparation_run_options(value.source_run_options.read())
        _require_source_accounting(manifest, accounting)
        producer = value.producer
        if (source != manifest.source
                or (value.total_records, value.total_physical_alternate_slots, value.total_biological_alternates)
                != (accounting.total_records, accounting.total_physical_alternate_slots, accounting.total_biological_alternates)
                or audit != SourceAuditResult(audit_requested, value.audit_status,
                    accounting.sha256, manifest.sha256, options.sha256,
                    source.executable_stack_sha256, producer.runtime_image_digest,
                    source.reference_manifest_sha256, source.input_sha256)
                or options.verify_source_ledger != audit_requested
                or producer.run_options_sha256 != options.sha256):
            raise ValueError
        config = next(parent for parent in manifest.parents if parent.role == "execution_limits")
        if (config.sha256, config.parent_version) != (producer.execution_config_sha256, producer.execution_config_version):
            raise ValueError
        # Header and summary bytes are bounded source metadata. Shards remain for
        # their necessary filter reads; the original VCF is never opened here.
        read_artifact(value.source_root, manifest.header)
        for report in manifest.summary.reports:
            read_artifact(value.source_root, report)
        journals = []
        for stage, root, parent, head, output, result_sha in (
            ("source", value.source_journal_root, producer.source_stage_parent_sha256,
             producer.source_journal_head_sha256, value.source_stage_output, manifest.sha256),
            ("source_audit", value.source_audit_journal_root, producer.source_audit_stage_parent_sha256,
             producer.source_audit_journal_head_sha256, value.source_audit_stage_output, audit.sha256),
        ):
            expected_parent = _preparation_stage_parent(stage,
                case_manifest_sha256=source.case_manifest_sha256,
                selected_sample_sha256=source.selected_sample_sha256,
                preceding_sha256=source.concordance_sha256 if stage == "source" else accounting.sha256,
                execution_limits_sha256=producer.execution_config_sha256,
                executable_stack_sha256=source.executable_stack_sha256,
                runtime_image_digest=producer.runtime_image_digest,
                run_options_sha256=options.sha256)
            if expected_parent != parent:
                raise ValueError
            from genosapiens.wgs_ledger import read_journal_head
            final = read_journal_head(root, stage=stage, parent_sha256=parent, expected_head_sha256=head)
            if final.event_type != "stage_completed":
                raise ValueError
            if final.artifacts != (Artifact("stage-output.json", output.sha256, output.size_bytes),):
                raise ValueError
            parents = {item.role: item for item in final.parents}
            for role, sha, schema, version in (
                ("execution_limits", config.sha256, config.parent_schema, config.parent_version),
                ("executable_stack", source.executable_stack_sha256, "genosapiens.executable-stack", 1),
                ("runtime_image", producer.runtime_image_digest[7:], "genosapiens.runtime-image", 1),
            ):
                if parents[role] != ParentIdentity(role, schema, version, sha):
                    raise ValueError
            output_content = _json(output.read())
            if (set(output_content) != {"schema", "version", "stage", "partition_output_sha256", "output_sha256", "status"}
                    or output_content["schema"] != output.schema or output_content["version"] != 1
                    or output_content["stage"] != stage or output_content["status"] != "complete"
                    or output_content["output_sha256"] != result_sha):
                raise ValueError
            if not _digest(output_content['partition_output_sha256']):
                raise ValueError
            journals.append(CompletedJournalBinding(stage, Path(root), parent, head))
            saved_state = value.source_state if stage == "source" else value.source_audit_state
            if saved_state is not None:
                current, _ = _completed_source_state(Path(value.source_root), Path(value.source_state_root), stage,
                    manifest=manifest, audit=audit, stage_parent=parent, outer_head=head)
                if current != saved_state:
                    from genosapiens.wgs_migration import matches_migrated_identity
                    if (replace(current, file_identity=saved_state.file_identity) != saved_state
                            or not matches_migrated_identity(
                                Path(value.source_state_root) / saved_state.relative_path,
                                saved_state.file_identity)):
                        raise ValueError
        return ResolvedCompletedSource(value, manifest, accounting, audit, options, tuple(journals))
    except (OSError, ValueError, TypeError, KeyError, IndexError, AttributeError, StopIteration):
        raise WgsPreparationError("completed source input is unavailable or incompatible") from None


def source_accounting_from_build(result: WgsSourceBuildResult) -> SourceAccountingClosure:
    """Seal counters supplied by the completed source writer and reducers."""
    if (type(result) is not WgsSourceBuildResult
            or not isinstance(result.qc_summary, VcfQcIncrementalSummary)
            or not isinstance(result.normalization_summary, WgsNormalizationSummary)
            or not isinstance(result.class_summary, IncrementalVariantClassSummary)):
        raise WgsPreparationError("completed source accounting is unavailable")
    manifest = result.manifest
    reports = manifest.summary.reports
    if (tuple(item.name for item in reports)
            != ("qc.json", "normalization.json", "classes.json")
            or tuple(hashlib.sha256(canonical(asdict(item))).hexdigest() for item in (
                result.qc_summary, result.normalization_summary, result.class_summary,
            )) != tuple(item.sha256 for item in reports)
            or result.qc_summary.total_records != manifest.summary.total_records
            or result.qc_summary.total_physical_alternate_slots
            != manifest.summary.total_physical_alternate_slots
            or result.qc_summary.total_alternate_alleles
            != manifest.summary.total_biological_alternates
            or result.class_summary.total_records != manifest.summary.total_records
            or result.normalization_summary.total_source_records
            != manifest.summary.total_records
            or result.normalization_summary.total_physical_alternate_slots
            != manifest.summary.total_physical_alternate_slots):
        raise WgsPreparationError("completed source accounting is inconsistent")
    return SourceAccountingClosure(
        manifest.sha256, manifest.source.input_sha256,
        manifest.header.sha256, *(item.sha256 for item in reports),
        manifest.summary.total_records,
        manifest.summary.total_physical_alternate_slots,
        manifest.summary.total_biological_alternates,
    )


def _require_source_accounting(
    manifest: LedgerManifest, accounting: SourceAccountingClosure,
) -> None:
    reports = manifest.summary.reports
    if (type(accounting) is not SourceAccountingClosure
            or tuple(item.name for item in reports)
            != ("qc.json", "normalization.json", "classes.json")
            or (accounting.source_manifest_sha256, accounting.input_sha256,
                accounting.header_sha256)
            != (manifest.sha256, manifest.source.input_sha256, manifest.header.sha256)
            or (accounting.qc_summary_sha256,
                accounting.normalization_summary_sha256,
                accounting.class_summary_sha256)
            != tuple(item.sha256 for item in reports)
            or (accounting.total_records,
                accounting.total_physical_alternate_slots,
                accounting.total_biological_alternates)
            != (manifest.summary.total_records,
                manifest.summary.total_physical_alternate_slots,
                manifest.summary.total_biological_alternates)):
        raise WgsPreparationError("source accounting is inconsistent")


class _BoundedReferenceFetcher:
    __slots__ = (
        "_fetch", "_max_request_bytes", "_max_record_bytes", "_record_number",
        "_requested_bytes", "_requests",
    )

    def __init__(self, fetch: Callable[[str, int, int], str], maximum: int,
                 *, first_record: int = 1) -> None:
        if not callable(fetch) or not _count(maximum, 1) or not _count(first_record, 1):
            raise WgsPreparationError("reference fetch limits are invalid")
        self._fetch = fetch
        self._max_request_bytes = maximum
        self._max_record_bytes = maximum
        self._record_number = first_record - 1
        self._requested_bytes = 0
        self._requests = 0

    def begin_record(self, record_number: int) -> None:
        if record_number != self._record_number + 1:
            raise WgsPreparationError("reference fetch source order is invalid")
        self._record_number = record_number
        self._requested_bytes = 0
        self._requests = 0

    def __call__(self, contig: str, position: int, length: int) -> str:
        if (
            self._record_number <= 0
            or not _count(length, 1)
            or length > self._max_request_bytes
            or self._requested_bytes + length > self._max_record_bytes
            or self._requests >= self._max_record_bytes
        ):
            raise WgsPreparationError("reference fetch exceeds its bounded size")
        self._requested_bytes += length
        self._requests += 1
        return self._fetch(contig, position, length)


def _rows(cursor, batch_size: int):
    while True:
        batch = cursor.fetchmany(batch_size)
        if not batch:
            return
        yield from batch


@contextmanager
def _source_accounting_connection(path: Path, *, cache_bytes: int):
    """One private retained database for source QC, class and normalization state."""
    _private_directory(path.parent)
    descriptor = -1
    connection = None
    try:
        descriptor = os.open(path, os.O_WRONLY | os.O_CREAT | os.O_EXCL | os.O_NOFOLLOW,
                             0o600)
        os.close(descriptor)
        descriptor = -1
        connection = sqlite3.connect(path)
        connection.execute("PRAGMA journal_mode=DELETE")
        connection.execute("PRAGMA synchronous=FULL")
        connection.execute("PRAGMA foreign_keys=ON")
        connection.execute(f"PRAGMA cache_size=-{max(1, cache_bytes // 1024)}")
        connection.execute("PRAGMA temp_store=FILE")
        connection.execute("PRAGMA mmap_size=0")
        connection.execute("PRAGMA trusted_schema=OFF")
        yield connection
    except (OSError, sqlite3.Error):
        raise WgsPreparationError("source accounting state is unavailable") from None
    finally:
        if descriptor >= 0:
            os.close(descriptor)
        if connection is not None:
            connection.close()


@contextmanager
def _source_stage_state(root, state_path, *, cache_bytes, lock, binding, resume=False):
    if binding is None:
        with _source_accounting_connection(state_path, cache_bytes=cache_bytes) as connection:
            yield connection, None
    else:
        _private_directory(state_path.parent)
        with _source_call("source_checkpoint_open", StageCommitState,
                         root, state_path, binding, lock=lock, create=not resume) as state:
            yield state.connection, state


class _NormalizationSummaryAccumulator:
    __slots__ = (
        "input_sha256", "concordance_sha256", "expected_records", "batch_size",
        "path", "_connection", "_records", "_physical", "_closed", "_finalized",
        "_owns_connection", "_pending_records", "_pending_occurrences",
        "_pending_alleles",
    )

    def __init__(
        self, path: Path | None = None, *, input_sha256: str, concordance_sha256: str,
        expected_records: int, cache_bytes: int, batch_size: int,
        connection: sqlite3.Connection | None = None,
        resume: bool = False,
    ) -> None:
        self.input_sha256 = input_sha256
        self.concordance_sha256 = concordance_sha256
        self.expected_records = expected_records
        self.batch_size = batch_size
        self.path = path
        if (path is None) == (connection is None) or type(resume) is not bool or (resume and connection is None):
            raise WgsPreparationError("normalization summary state binding is invalid")
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
                connection.execute("PRAGMA foreign_keys=ON")
                connection.execute(f"PRAGMA cache_size=-{max(1, cache_bytes // 1024)}")
                connection.execute("PRAGMA temp_store=FILE")
                connection.execute("PRAGMA mmap_size=0")
                connection.execute("PRAGMA trusted_schema=OFF")
            self._connection = connection
            if not resume:
                self._connection.executescript("""
                CREATE TABLE normalization_records(
                    record_number INTEGER PRIMARY KEY,
                    source_record_id TEXT UNIQUE NOT NULL,
                    reference_mismatch INTEGER NOT NULL CHECK(reference_mismatch IN (0,1))
                );
                CREATE TABLE normalization_occurrences(
                    source_allele_id TEXT PRIMARY KEY,
                    record_number INTEGER NOT NULL REFERENCES normalization_records(record_number),
                    normalized_allele_id TEXT,
                    quarantine_reason TEXT,
                    shard_ordinal INTEGER,
                    shard_row INTEGER,
                    CHECK((shard_ordinal IS NULL) = (shard_row IS NULL)),
                    CHECK((normalized_allele_id IS NULL) != (quarantine_reason IS NULL))
                );
                CREATE TABLE normalization_alleles(normalized_allele_id TEXT PRIMARY KEY);
                CREATE TABLE normalization_meta(
                    singleton INTEGER PRIMARY KEY CHECK(singleton=1),
                    input_sha256 TEXT NOT NULL,
                    concordance_sha256 TEXT NOT NULL,
                    expected_records INTEGER NOT NULL
                );
            """)
                self._connection.execute(
                    "INSERT INTO normalization_meta VALUES (1,?,?,?)",
                    (input_sha256, concordance_sha256, expected_records),
                )
                if self._owns_connection:
                    self._connection.commit()
                records = physical = 0
            else:
                records = self._connection.execute(
                    "SELECT count(*) FROM normalization_records"
                ).fetchone()[0]
                physical = self._connection.execute(
                    "SELECT count(*) FROM normalization_occurrences"
                ).fetchone()[0]
                if (type(records) is not int or type(physical) is not int
                        or records < 0 or records > expected_records or physical < records
                        or self._connection.execute(
                            "SELECT input_sha256,concordance_sha256,expected_records "
                            "FROM normalization_meta WHERE singleton=1"
                        ).fetchone() != (
                            input_sha256, concordance_sha256, expected_records,
                        )):
                    raise WgsPreparationError("normalization summary state is inconsistent")
        except (OSError, sqlite3.Error):
            if descriptor >= 0:
                os.close(descriptor)
            raise WgsPreparationError("cannot create normalization summary state") from None
        self._records = records
        self._physical = physical
        self._closed = False
        self._finalized = False

        self._pending_records = []
        self._pending_occurrences = []
        self._pending_alleles = []

    def flush(self) -> None:
        """Apply bounded per-shard normalization facts with batched uniqueness."""
        if self._closed or self._finalized:
            raise WgsPreparationError("normalization summary state is unavailable")
        if not self._pending_records:
            return
        try:
            self._connection.executemany(
                "INSERT INTO normalization_records VALUES (?,?,?)",
                self._pending_records,
            )
            self._connection.executemany(
                "INSERT INTO normalization_occurrences VALUES (?,?,?,?,?,?)",
                self._pending_occurrences,
            )
            self._connection.executemany(
                "INSERT OR IGNORE INTO normalization_alleles VALUES (?)",
                self._pending_alleles,
            )
        except sqlite3.Error:
            raise WgsPreparationError("cannot update normalization summary state") from None
        self._pending_records.clear()
        self._pending_occurrences.clear()
        self._pending_alleles.clear()

    def add(self, result: SourceRecordNormalization, *,
            shard_ordinal: int | None = None, shard_row: int | None = None) -> None:
        if (
            self._closed
            or self._finalized
            or not isinstance(result, SourceRecordNormalization)
            or result.input_sha256 != self.input_sha256
            or result.record.record_number != self._records + 1
            or (shard_ordinal is None) != (shard_row is None)
            or (shard_ordinal is not None and (
                not _count(shard_ordinal, 1) or not _count(shard_row, 1)
            ))
        ):
            raise WgsPreparationError("normalization source order is invalid")
        mismatch = any(
            item.quarantine_reason is QuarantineReason.REFERENCE_MISMATCH
            for item in result.alleles
        )
        self._pending_records.append(
            (result.record.record_number, result.source_record_id, int(mismatch))
        )
        for item in result.alleles:
            normalized_id = (
                item.normalized_allele.normalized_allele_id
                if item.normalized_allele is not None else None
            )
            reason = item.quarantine_reason.value if item.quarantine_reason else None
            self._pending_occurrences.append((
                item.source_allele_id, result.record.record_number, normalized_id,
                reason, shard_ordinal, shard_row,
            ))
            if normalized_id is not None:
                self._pending_alleles.append((normalized_id,))
        self._records += 1
        self._physical += len(result.alleles)
        if len(self._pending_records) >= self.batch_size or len(self._pending_occurrences) >= 8192:
            self.flush()
            if self._owns_connection:
                self._connection.commit()

    def add_projection(self, *, record_number: int, source_record_id: str,
                       physical_alternates: int, projection: dict,
                       shard_ordinal: int, shard_row: int) -> None:
        """Reduce typed compact fields from a closed worker row once."""
        try:
            if (self._closed or self._finalized
                    or record_number != self._records + 1
                    or not _count(shard_ordinal, 1) or not _count(shard_row, 1)
                    or type(projection) is not dict
                    or projection.get("reference_identity_id") is None
                    or projection.get("reference_role") is None
                    or type(projection.get("alleles")) is not list
                    or len(projection["alleles"]) != physical_alternates):
                raise ValueError
            pending = []
            normalized = []
            mismatch = False
            for allele_number, item in enumerate(projection["alleles"], 1):
                expected = source_allele_id(source_record_id, allele_number)
                if (type(item) is not dict
                        or item.get("source_allele_id") != expected
                        or item.get("allele_number") != allele_number):
                    raise ValueError
                allele = item.get("normalized_allele")
                reason = item.get("quarantine_reason")
                if (allele is None) == (reason is None):
                    raise ValueError
                if allele is not None:
                    if type(allele) is not dict:
                        raise ValueError
                    normalized_id = allele.get("normalized_allele_id")
                    if (type(normalized_id) is not str
                            or re.fullmatch(r"allele-[0-9a-f]{24}", normalized_id) is None):
                        raise ValueError
                    normalized.append((normalized_id,))
                else:
                    normalized_id = None
                    reason = QuarantineReason(reason).value
                    mismatch |= reason == QuarantineReason.REFERENCE_MISMATCH.value
                pending.append((expected, record_number, normalized_id, reason,
                                shard_ordinal, shard_row))
        except (ValueError, TypeError, KeyError):
            raise WgsPreparationError("source normalization projection is invalid") from None
        self._pending_records.append((record_number, source_record_id, int(mismatch)))
        self._pending_occurrences.extend(pending)
        self._pending_alleles.extend(normalized)
        self._records += 1
        self._physical += physical_alternates
        if len(self._pending_records) >= self.batch_size or len(self._pending_occurrences) >= 8192:
            self.flush()

    def finalize(self) -> WgsNormalizationSummary:
        if self._closed or self._finalized or self._records != self.expected_records:
            raise WgsPreparationError("normalization summary accounting is invalid")
        self.flush()
        if self._owns_connection:
            self._connection.commit()
        grouped = tuple(
            NormalizationReasonCount(reason, count)
            for reason, count in _rows(
                self._connection.execute(
                    "SELECT quarantine_reason,count(*) FROM normalization_occurrences "
                    "WHERE quarantine_reason IS NOT NULL GROUP BY quarantine_reason "
                    "ORDER BY quarantine_reason"
                ),
                self.batch_size,
            )
        )
        normalized_occurrences = self._connection.execute(
            "SELECT count(*) FROM normalization_occurrences WHERE normalized_allele_id IS NOT NULL"
        ).fetchone()[0]
        unique_normalized = self._connection.execute(
            "SELECT count(*) FROM normalization_alleles"
        ).fetchone()[0]
        mismatch_records = self._connection.execute(
            "SELECT count(*) FROM normalization_records WHERE reference_mismatch=1"
        ).fetchone()[0]
        self._finalized = True
        return WgsNormalizationSummary(
            self.input_sha256,
            self.concordance_sha256,
            self._records,
            self._physical,
            normalized_occurrences,
            unique_normalized,
            self._physical - normalized_occurrences,
            mismatch_records,
            grouped,
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
            raise WgsPreparationError("normalization summary was not finalized")
        return False


def _source_record(record: VcfRecordQc) -> SourceSmallVariantRecord:
    return SourceSmallVariantRecord(
        record.record_number, record.contig, record.position,
        record.reference, record.alternates,
    )


def _source_payload(
    record: VcfRecordQc,
    normalization: SourceRecordNormalization,
    inventory: VariantRecordInventoryFacts,
    *, selected_sample_name: str,
) -> str:
    if record.raw_line_sha256 is None or record.raw_line is None:
        raise WgsPreparationError("source record lacks raw-line authentication")
    selected_index = tuple(call.sample_name for call in record.sample_calls).index(
        selected_sample_name
    )
    value = {
        "schema": SOURCE_PAYLOAD_SCHEMA,
        "version": SOURCE_PAYLOAD_VERSION,
        "classification": "controlled",
        "llm_policy": "no-llm",
        "selected_sample_name": selected_sample_name,
        "selected_sample_index": selected_index,
        "raw_line_sha256": record.raw_line_sha256,
        "qc": asdict(record),
        "normalization": {
            "reference_identity_id": normalization.reference_identity_id,
            "reference_role": normalization.reference_role,
            "alleles": [asdict(item) for item in normalization.alleles],
        },
        "inventory": asdict(inventory),
    }
    return canonical(value).decode("ascii")


def _publish_report(root: Path, name: str, value: object) -> Artifact:
    return _publish_or_reuse_small(root, name, canonical(asdict(value)))


def _publish_or_reuse_small(root: Path, name: str, content: bytes) -> Artifact:
    artifact = Artifact(name, hashlib.sha256(content).hexdigest(), len(content))
    path = root / name
    if path.exists() or path.is_symlink():
        return attest_existing_artifact(root, artifact)
    return publish_bytes(root, name, content)


_SOURCE_REDUCER_TABLES = (
    "qc_seen", "qc_contigs", "qc_samples", "qc_ploidies", "qc_named_counts",
    "qc_totals", "qc_meta", "class_seen", "class_classes",
    "class_mosaic_candidates", "class_scalars", "class_screen_counts", "class_meta",
    "normalization_records", "normalization_occurrences", "normalization_alleles",
    "normalization_meta",
)


def _source_reducer_positions(connection):
    """Read singleton counters and indexed tails, never scan source records."""
    positions = {"qc_totals": connection.execute(
        "SELECT records FROM qc_totals WHERE singleton=1"
    ).fetchone()[0]}
    for table in ("qc_seen", "class_seen", "normalization_records"):
        row = connection.execute(
            f"SELECT record_number FROM {table} ORDER BY record_number DESC LIMIT 1"
        ).fetchone()
        positions[table] = row[0] if row is not None else 0
    return positions


def _recover_source_reducers(root, connection, *, shards, binding, header,
                             expected_records, allele_balance_screen, config):
    """Rebuild legacy non-atomic reducers from retained, authenticated shards.

    The original database and all its tables are retained. A separate recovery
    database commits each replayed shard, so an interrupted recovery resumes.
    Installation renames the old reducer tables and copies the repaired tables
    in one transaction. Source journals, shard cursors and files stay intact.
    """
    from genosapiens.evidence_model import _typed
    import sys
    import time

    committed_records = shards[-1].last_record if shards else 0
    positions = _source_reducer_positions(connection)
    if all(value == committed_records for value in positions.values()):
        return
    if (any(type(value) is not int or value < committed_records
            or value > expected_records for value in positions.values())
            or connection.in_transaction
            or connection.execute("SELECT 1 FROM sqlite_master WHERE type='table' "
                                  "AND name='source_legacy_qc_seen'").fetchone()):
        raise WgsPreparationError("source reducer recovery binding is invalid")
    # A reducer can flush before the current shard commits, but cannot advance
    # beyond that one bounded shard. Refuse unrelated or unbounded damage.
    if max(positions.values()) - committed_records > config.source_shards.max_uncompressed_rows:
        raise WgsPreparationError("source reducer recovery exceeds the pending shard")
    print("Recovering interrupted source accounting from completed shards; "
          "original checkpoint tables will be retained.", file=sys.stderr, flush=True)
    recovery_root = root / ".recovery"
    _private_directory(recovery_root)
    recovery_path = recovery_root / "source-reducers.sqlite"
    recovery_binding = canonical({
        "schema": "genosapiens.source-reducer-recovery", "version": 1,
        "stage_binding_sha256": hashlib.sha256(binding.as_bytes()).hexdigest(),
        "header_sha256": header.header_sha256, "expected_records": expected_records,
        "allele_balance_screen": allele_balance_screen,
        "committed_records": committed_records, "committed_shards": len(shards),
        "last_shard_sha256": shards[-1].sha256 if shards else None,
    })
    exists = recovery_path.exists() or recovery_path.is_symlink()
    if exists:
        metadata = recovery_path.lstat()
        if (not stat.S_ISREG(metadata.st_mode) or stat.S_IMODE(metadata.st_mode) != 0o600
                or metadata.st_uid != os.getuid() or metadata.st_nlink != 1):
            raise WgsPreparationError("source reducer recovery file is unsafe")
    else:
        descriptor = os.open(recovery_path, os.O_WRONLY | os.O_CREAT | os.O_EXCL | os.O_NOFOLLOW, 0o600)
        os.close(descriptor)
    recovered = sqlite3.connect(recovery_path)
    try:
        recovered.execute("PRAGMA foreign_keys=ON")
        recovered.execute("PRAGMA synchronous=FULL")
        recovered.execute("PRAGMA journal_mode=DELETE")
        recovered.execute("PRAGMA trusted_schema=OFF")
        recovered.execute("PRAGMA temp_store=FILE")
        recovered.execute(f"PRAGMA cache_size=-{max(1, config.sqlite.cache_bytes // 1024)}")
        ready = recovered.execute("SELECT 1 FROM sqlite_master WHERE type='table' "
                                  "AND name='recovery_cursor'").fetchone()
        if ready:
            row = recovered.execute("SELECT binding,ordinal,records FROM recovery_cursor WHERE singleton=1").fetchone()
            if (row is None or row[0] != recovery_binding or type(row[1]) is not int
                    or not 0 <= row[1] <= len(shards)
                    or row[2] != (shards[row[1] - 1].last_record if row[1] else 0)
                    or any(value != row[2] for value in _source_reducer_positions(recovered).values())):
                raise WgsPreparationError("source reducer recovery checkpoint is inconsistent")
            recovered_shards = row[1]
        else:
            # Initialization contains no committed recovery work. Preserve an
            # interrupted initializer before creating a fresh recovery database.
            if exists:
                recovered.close()
                os.rename(recovery_path, recovery_root / ("incomplete-" + os.urandom(16).hex() + ".sqlite"))
                return _recover_source_reducers(root, connection, shards=shards, binding=binding,
                    header=header, expected_records=expected_records,
                    allele_balance_screen=allele_balance_screen, config=config)
            recovered_shards = 0
        with (
            closing(QcSummaryAccumulator(header, expected_records=expected_records,
                                 connection=recovered, resume=bool(ready))) as qc,
            closing(IncrementalVariantInventoryAccumulator(header, expected_records=expected_records,
                allele_balance_screen=allele_balance_screen, connection=recovered, resume=bool(ready))) as classes,
            closing(_NormalizationSummaryAccumulator(input_sha256=header.input_sha256,
                concordance_sha256=header.reference_concordance_sha256, expected_records=expected_records,
                cache_bytes=config.sqlite.cache_bytes, batch_size=config.sqlite.cursor_batch_rows,
                connection=recovered, resume=bool(ready))) as normalization,
        ):
            if not ready:
                recovered.execute("CREATE TABLE recovery_cursor(singleton INTEGER PRIMARY KEY CHECK(singleton=1), "
                                  "binding BLOB NOT NULL,ordinal INTEGER NOT NULL,records INTEGER NOT NULL)")
                recovered.execute("INSERT INTO recovery_cursor VALUES (1,?,0,0)", (recovery_binding,))
                recovered.commit()
            next_notice = time.monotonic() + 30
            for shard in shards[recovered_shards:]:
                recovered.execute("BEGIN IMMEDIATE")
                for shard_row, row in enumerate(iter_shard(root, shard), 1):
                    payload = json.loads(row.payload_json)
                    record = _typed(payload["qc"], VcfRecordQc, decode=True)
                    if (payload.get("schema") != SOURCE_PAYLOAD_SCHEMA or payload.get("version") != SOURCE_PAYLOAD_VERSION
                            or record.record_number != row.record_number
                            or record.source_record_id != row.source_record_id
                            or row.source_allele_ids != tuple(source_allele_id(row.source_record_id, number)
                                for number in range(1, len(record.alternates) + 1))):
                        raise WgsPreparationError("source recovery row binding is invalid")
                    qc.add(record)
                    classes.add_projection(payload["inventory"])
                    normalization.add_projection(record_number=row.record_number,
                        source_record_id=row.source_record_id, physical_alternates=len(record.alternates),
                        projection=payload["normalization"], shard_ordinal=shard.ordinal, shard_row=shard_row)
                qc.flush()
                classes.flush()
                normalization.flush()
                recovered.execute("UPDATE recovery_cursor SET ordinal=?,records=? WHERE singleton=1",
                                  (shard.ordinal, shard.last_record))
                recovered.commit()
                if time.monotonic() >= next_notice or shard.ordinal == len(shards):
                    print(f"Source accounting recovery: {shard.ordinal}/{len(shards)} "
                          "completed shards.", file=sys.stderr, flush=True)
                    next_notice = time.monotonic() + 30
        if any(value != committed_records for value in _source_reducer_positions(recovered).values()):
            raise WgsPreparationError("source reducer recovery did not close")
        connection.execute("BEGIN IMMEDIATE")
        try:
            for table in _SOURCE_REDUCER_TABLES:
                connection.execute(f"ALTER TABLE {table} RENAME TO source_legacy_{table}")
            for table in _SOURCE_REDUCER_TABLES:
                statement = recovered.execute("SELECT sql FROM sqlite_master WHERE type='table' AND name=?",
                                               (table,)).fetchone()[0]
                connection.execute(statement)
                columns = len(recovered.execute(f"PRAGMA table_info({table})").fetchall())
                parameters = ",".join("?" for _ in range(columns))
                cursor = recovered.execute(f"SELECT * FROM {table}")
                while batch := cursor.fetchmany(config.sqlite.cursor_batch_rows):
                    connection.executemany(f"INSERT INTO {table} VALUES ({parameters})", batch)
            connection.execute("CREATE TABLE source_reducer_recovery(binding BLOB NOT NULL)")
            connection.execute("INSERT INTO source_reducer_recovery VALUES (?)", (recovery_binding,))
            connection.commit()
        except BaseException:
            connection.rollback()
            raise
    finally:
        recovered.close()


def build_wgs_source_ledger(
    vcf_path: str | Path,
    *,
    concordance_report: ReferenceConcordanceReport | PublicValidationConcordanceReport,
    fetch_reference: Callable[[str, int, int], str],
    output_root: str | Path,
    selected_sample_name: str,
    selected_sample_sha256: str,
    case_manifest_sha256: str,
    executable_stack_sha256: str,
    execution_config: BoundedWgsExecutionConfig,
    allele_balance_screen: bool = False,
    source_state_path: str | Path | None = None,
    stage_binding: StageRunBinding | None = None,
    resume: bool = False,
    worker_reference: SourceWorkerReference | None = None,
    worker_count: int | None = None,
    operator_case_binding=None,
) -> WgsSourceBuildResult:
    """Stream the complete source into immutable shards before any early filter."""

    if (
        not isinstance(execution_config, BoundedWgsExecutionConfig)
        or type(allele_balance_screen) is not bool
        or type(resume) is not bool
        or resume and stage_binding is None
        or not isinstance(selected_sample_name, str)
        or not selected_sample_name
        or any(not _digest(value) for value in (
            selected_sample_sha256, case_manifest_sha256, executable_stack_sha256,
        ))
        or execution_config.source_shards.max_uncompressed_rows > MAX_SHARD_ROWS
        or execution_config.source_shards.max_uncompressed_bytes > MAX_SHARD_BYTES
        or execution_config.input.max_record_bytes > MAX_ROW_BYTES
        or (stage_binding is not None and (
            type(stage_binding) is not StageRunBinding
            or stage_binding.stage != "source"
            or stage_binding.executable_stack_sha256 != executable_stack_sha256
        ))
        or (worker_reference is not None and type(worker_reference) is not SourceWorkerReference)
        or (worker_count is not None and worker_count not in (1, 4))
        or (worker_count is not None and worker_reference is None)
    ):
        raise WgsPreparationError("source ledger build binding is invalid")
    root = Path(output_root)
    if not root.is_absolute():
        raise WgsPreparationError("source ledger root is invalid")
    execution_limits_sha256 = hashlib.sha256(
        bounded_wgs_execution_config_as_json(execution_config)
    ).hexdigest()
    state_path = Path(source_state_path) if source_state_path is not None else root / "state/source.sqlite"
    if not state_path.is_absolute() or state_path.name != "source.sqlite":
        raise WgsPreparationError("source accounting state path is invalid")
    with ExitStack() as contexts:
        source_lock = contexts.enter_context(stage_lock(root, "source"))
        source_connection, stage_state = contexts.enter_context(_source_stage_state(
            root, state_path, cache_bytes=execution_config.sqlite.cache_bytes,
            lock=source_lock, binding=stage_binding, resume=resume,
        ))
        if stage_state is not None and resume:
            tail = source_connection.execute(
                "SELECT event_bytes FROM events ORDER BY ordinal DESC LIMIT 1"
            ).fetchone()
            if tail is not None:
                last_event = load_model(tail[0], JournalEntry)
                if last_event.event_type == "stage_completed":
                    manifest_content = _read_bounded_manifest(root / "manifest.json")
                    manifest = load_model(manifest_content, LedgerManifest)
                    manifest_artifact = Artifact(
                        "manifest.json", hashlib.sha256(manifest_content).hexdigest(),
                        len(manifest_content),
                    )
                    if last_event.artifacts != (manifest_artifact,):
                        raise WgsPreparationError("completed source checkpoint is inconsistent")
                    return WgsSourceBuildResult(
                        manifest, manifest_artifact, None, None, None,
                    )
        shards: list[Shard] = []
        resume_identity = None
        if operator_case_binding is not None:
            from genosapiens.wgs_operator_invocation import OperatorCaseBinding
            if type(operator_case_binding) is not OperatorCaseBinding:
                raise WgsPreparationError("source authentication binding is invalid")
            resume_identity = operator_case_binding.authenticated_source_identity(vcf_path, concordance_report)
        if resume:
            try:
                identity_row = source_connection.execute(
                    "SELECT identity_bytes FROM source_input WHERE singleton=1"
                ).fetchone()
                if identity_row is None:
                    raise WgsPreparationError("source resume input identity is missing")
                identity = json.loads(identity_row[0])
                if (type(identity) is not list or len(identity) != 5
                        or any(type(value) is not int or value < 0 for value in identity)):
                    raise WgsPreparationError("source resume input identity is invalid")
                resume_identity = tuple(identity)
                rows = source_connection.execute(
                    "SELECT ordinal,descriptor FROM source_shards ORDER BY ordinal"
                )
                for ordinal, descriptor in rows:
                    shard = load_model(descriptor, Shard)
                    if (shard.ordinal != ordinal or ordinal != len(shards) + 1
                            or shard.first_record != (1 if not shards else shards[-1].last_record + 1)
                            or shard.previous_shard_sha256 != (
                                None if not shards else shards[-1].sha256
                            )):
                        raise WgsPreparationError("source committed shard cursor is invalid")
                    shards.append(shard)
                if shards:
                    latest = stage_state.committed(f"source-shard-{len(shards):06d}")
                    if (latest is None or latest.cursor != (
                            shards[-1].last_record, shards[-1].ordinal)
                            or latest.outputs != (shards[-1].artifact,)):
                        raise WgsPreparationError("source committed shard cursor is inconsistent")
                    proof = getattr(latest.outputs[0], "_publication", None)
                    object.__setattr__(shards[-1], "_publication", proof)
                    object.__setattr__(shards[-1].artifact, "_publication", proof)
            except (sqlite3.Error, ValueError, TypeError, KeyError):
                raise WgsPreparationError("source resume state is inconsistent") from None
        stream = contexts.enter_context(_source_call("source_input_open", iter_vcf_qc_records,
            vcf_path,
            concordance_report=concordance_report,
            max_record_bytes=execution_config.input.max_record_bytes,
            max_header_bytes=execution_config.input.max_header_bytes,
            resume_input_identity=resume_identity,
            skip_records=shards[-1].last_record if shards else 0,
        ))
        header = stream.header
        if selected_sample_name not in header.sample_names:
            raise WgsPreparationError("selected sample is absent from the source header")
        if stage_binding is not None and (
            stage_binding.execution_config_sha256 != execution_limits_sha256
            or stage_binding.input_sha256 != header.input_sha256
            or stage_binding.reference_sha256 != concordance_report.reference_manifest_sha256
        ):
            raise WgsPreparationError("source stage binding is invalid")
        if stage_state is not None and not resume:
            source_connection.execute(
                "CREATE TABLE source_input(singleton INTEGER PRIMARY KEY CHECK(singleton=1),"
                "identity_bytes BLOB NOT NULL)"
            )
            source_connection.execute(
                "CREATE TABLE source_shards(ordinal INTEGER PRIMARY KEY,descriptor BLOB NOT NULL)"
            )
            source_connection.execute(
                "INSERT INTO source_input VALUES (1,?)",
                (canonical(stream._input_identity),),
            )
        if stage_state is not None and resume:
            _source_call("source_reducer_recovery", _recover_source_reducers,
                root, source_connection, shards=shards, binding=stage_binding,
                header=header, expected_records=concordance_report.total_records,
                allele_balance_screen=allele_balance_screen, config=execution_config)
        batch: list[EncodedLedgerRow] = []
        batch_bytes = 0
        previous: Shard | None = shards[-1] if shards else None
        first_record = previous.last_record + 1 if previous else 1
        bounded_reference = _BoundedReferenceFetcher(
            fetch_reference, execution_config.input.max_record_bytes,
            first_record=first_record,
        )
        with nullcontext():
            with (
                _source_call("source_qc_initialize", QcSummaryAccumulator,
                    header,
                    expected_records=concordance_report.total_records,
                    connection=source_connection,
                    resume=resume,
                ) as qc_accumulator,
                _source_call("source_class_initialize", IncrementalVariantInventoryAccumulator,
                    header,
                    expected_records=concordance_report.total_records,
                    allele_balance_screen=allele_balance_screen,
                    connection=source_connection,
                    resume=resume,
                ) as class_accumulator,
                _source_call("source_normalization_initialize", _NormalizationSummaryAccumulator,
                    input_sha256=header.input_sha256,
                    concordance_sha256=header.reference_concordance_sha256,
                    expected_records=concordance_report.total_records,
                    cache_bytes=execution_config.sqlite.cache_bytes,
                    batch_size=execution_config.sqlite.cursor_batch_rows,
                    connection=source_connection,
                    resume=resume,
                ) as normalization_accumulator,
            ):
                source_connection.commit()
                next_partition = f"source-shard-{len(shards) + 1:06d}"
                if stage_state is not None and resume:
                    incomplete = stage_state.open_attempt(next_partition)
                    if incomplete is not None:
                        stage_state.archive_orphans(
                            names=(f"shard-{len(shards) + 1:06d}.jsonl.gz",),
                            attempt_id=incomplete.attempt_id,
                        )
                        stage_state.fail(incomplete, failure_code="process_interrupted",
                                         interrupted=True)
                current_start = _begin_source_checkpoint_attempt(
                    stage_state, next_partition, resume=resume,
                ) if stage_state is not None and first_record <= concordance_report.total_records else None
                header_artifact = _publish_or_reuse_small(root, "header.json", canonical({
                    "schema": HEADER_SCHEMA,
                    "version": 1,
                    "classification": "controlled",
                    "llm_policy": "no-llm",
                    "selected_sample_name": selected_sample_name,
                    "header": asdict(header),
                }))
                if current_start is not None:
                    # Reducers may flush before a shard reaches its row/byte
                    # limit. Their writes must share the shard's transaction.
                    source_connection.execute("BEGIN IMMEDIATE")

                def commit_source_shard(shard: Shard) -> None:
                    def apply_reducers(_connection):
                        qc_accumulator.flush()
                        class_accumulator.flush()
                        normalization_accumulator.flush()
                        _connection.execute(
                            "INSERT INTO source_shards VALUES (?,?)",
                            (shard.ordinal, shard.as_bytes()),
                        ) if stage_state is not None else None

                    if stage_state is None:
                        apply_reducers(source_connection)
                        source_connection.commit()
                    else:
                        _source_call("source_shard_commit", stage_state.commit_partition,
                            current_start,
                            cursor=(shard.last_record, shard.ordinal),
                            outputs=(shard.artifact,),
                            counts={
                                "records": shard.rows,
                                "physical_alternates": shard.physical_alternate_slots,
                                "biological_alternates": shard.biological_alternates,
                            },
                            reducer=apply_reducers,
                            prepared_transaction=True,
                        )

                def append_source_row(qc_record, encoded, apply_facts):
                    nonlocal previous, first_record, batch, batch_bytes, current_start
                    if qc_record.quality is not None and (
                        not math.isfinite(qc_record.quality)
                        or qc_record.quality < 0.0
                    ):
                        raise WgsPreparationError("source quality is invalid")
                    if (encoded.record_number != qc_record.record_number
                            or encoded.row.source_record_id != qc_record.source_record_id
                            or encoded.row.source_allele_ids != tuple(
                                source_allele_id(qc_record.source_record_id, number)
                                for number in range(1, len(qc_record.alternates) + 1)
                            )
                            or encoded.biological_alternates
                            != _alternate_allele_count(qc_record.alternates)):
                        raise WgsPreparationError("source worker source tuple is invalid")
                    row_bytes = len(encoded.content)
                    if row_bytes > execution_config.source_shards.max_uncompressed_bytes:
                        raise WgsPreparationError("source ledger row exceeds its bounded size")
                    if batch and (
                        len(batch) >= execution_config.source_shards.max_uncompressed_rows
                        or batch_bytes + row_bytes
                        > execution_config.source_shards.max_uncompressed_bytes
                    ):
                        previous = _source_call("source_shard_publish", publish_shard,
                            root, tuple(batch), ordinal=len(shards) + 1,
                            first_record=first_record, previous=previous,
                        )
                        shards.append(previous)
                        commit_source_shard(previous)
                        if stage_state is not None:
                            current_start = _begin_source_checkpoint_attempt(stage_state,
                                f"source-shard-{len(shards) + 1:06d}", resume=resume)
                            source_connection.execute("BEGIN IMMEDIATE")
                        first_record = qc_record.record_number
                        batch = []
                        batch_bytes = 0
                    qc_accumulator.add(qc_record)
                    apply_facts(len(shards) + 1, len(batch) + 1)
                    batch.append(encoded)
                    batch_bytes += row_bytes

                if worker_reference is None:
                    for qc_record in stream:
                        bounded_reference.begin_record(qc_record.record_number)
                        normalization = normalize_source_record(
                            _source_record(qc_record),
                            input_sha256=header.input_sha256,
                            reference_identity_id=concordance_report.reference_identity_id,
                            reference_role=concordance_report.reference_role,
                            fetch_reference=bounded_reference,
                        )
                        inventory = class_accumulator.facts(qc_record)
                        row = LedgerRow(
                            record_number=qc_record.record_number,
                            source_record_id=qc_record.source_record_id,
                            source_allele_ids=tuple(
                                item.source_allele_id for item in normalization.alleles
                            ),
                            biological_alternate_count=_alternate_allele_count(
                                qc_record.alternates
                            ),
                            payload_json=_source_payload(
                                qc_record, normalization, inventory,
                                selected_sample_name=selected_sample_name,
                            ),
                        )
                        encoded = encode_ledger_row(row)

                        def apply_serial(shard_ordinal, shard_row):
                            class_accumulator.add_facts(inventory)
                            normalization_accumulator.add(
                                normalization, shard_ordinal=shard_ordinal,
                                shard_row=shard_row,
                            )

                        _source_call("source_row_reduce", append_source_row, qc_record, encoded, apply_serial)
                else:
                    jobs = _iter_parallel_source_rows(
                        stream, root=root, reference=worker_reference,
                        workers=worker_count or execution_config.parallelism.source_workers,
                        max_inflight=execution_config.parallelism.max_inflight_cpu_batches,
                        batch_rows=execution_config.parallelism.source_batch_rows,
                        batch_bytes=execution_config.parallelism.max_source_batch_bytes,
                        input_sha256=header.input_sha256,
                        reference_identity_id=concordance_report.reference_identity_id,
                        reference_role=concordance_report.reference_role,
                        selected_sample_name=selected_sample_name,
                        sample_names=header.sample_names,
                        definitions=class_accumulator._definitions,
                        lengths=class_accumulator._lengths,
                        allele_balance_screen=allele_balance_screen,
                        max_record_bytes=execution_config.input.max_record_bytes,
                    )
                    with closing(jobs):
                        for qc_record, encoded in jobs:
                            payload = json.loads(encoded.row.payload_json)
                            if (payload.get("schema") != SOURCE_PAYLOAD_SCHEMA
                                    or payload.get("version") != SOURCE_PAYLOAD_VERSION
                                    or payload.get("selected_sample_name") != selected_sample_name
                                    or payload.get("raw_line_sha256") != qc_record.raw_line_sha256
                                    or type(payload.get("qc")) is not dict
                                    or payload["qc"].get("source_record_id")
                                    != qc_record.source_record_id
                                    or type(payload.get("inventory")) is not dict
                                    or type(payload.get("normalization")) is not dict
                                    or payload["normalization"].get("reference_identity_id")
                                    != concordance_report.reference_identity_id
                                    or payload["normalization"].get("reference_role")
                                    != concordance_report.reference_role):
                                raise WgsPreparationError("source worker payload is invalid")

                            def apply_parallel(shard_ordinal, shard_row):
                                class_accumulator.add_projection(payload["inventory"])
                                normalization_accumulator.add_projection(
                                    record_number=qc_record.record_number,
                                    source_record_id=qc_record.source_record_id,
                                    physical_alternates=len(qc_record.alternates),
                                    projection=payload["normalization"],
                                    shard_ordinal=shard_ordinal, shard_row=shard_row,
                                )

                            _source_call("source_row_reduce", append_source_row, qc_record, encoded, apply_parallel)
                if batch:
                    previous = _source_call("source_shard_publish", publish_shard,
                        root, tuple(batch), ordinal=len(shards) + 1,
                        first_record=first_record, previous=previous,
                    )
                    shards.append(previous)
                    commit_source_shard(previous)
                qc_summary = qc_accumulator.finalize()
                class_summary = class_accumulator.finalize()
                normalization_summary = normalization_accumulator.finalize()
                source_connection.commit()

            reports = (
                _publish_report(root, "qc.json", qc_summary),
                _publish_report(root, "normalization.json", normalization_summary),
                _publish_report(root, "classes.json", class_summary),
            )
            source = SourceBinding(
                header.input_sha256,
                header.input_size_bytes,
                case_manifest_sha256,
                selected_sample_sha256,
                header.reference_concordance_sha256,
                concordance_report.reference_manifest_sha256,
                executable_stack_sha256,
            )
            parents = (
                ParentIdentity(
                    "original_input", "genosapiens.original-input", 1,
                    header.input_sha256,
                ),
                ParentIdentity(
                    "case_manifest", "genosapiens.case-manifest", 1,
                    case_manifest_sha256,
                ),
                ParentIdentity(
                    "selected_sample", "genosapiens.selected-sample-identity", 1,
                    selected_sample_sha256,
                ),
                ParentIdentity(
                    "complete_concordance", concordance_report.schema, 1,
                    header.reference_concordance_sha256,
                ),
                ParentIdentity(
                    "reference_manifest", "genosapiens.public-resource-manifest", 1,
                    concordance_report.reference_manifest_sha256,
                ),
                ParentIdentity(
                    "execution_limits", execution_config.schema,
                    execution_config.version, execution_limits_sha256,
                ),
                ParentIdentity(
                    "executable_stack", "genosapiens.executable-stack", 1,
                    executable_stack_sha256,
                ),
            )
            manifest = LedgerManifest(
                source,
                header_artifact,
                parents,
                LedgerSummary(
                    qc_summary.total_records,
                    qc_summary.total_physical_alternate_slots,
                    qc_summary.total_alternate_alleles,
                    reports,
                ),
                tuple(shards),
            )
            if (root / "manifest.json").exists() or (root / "manifest.json").is_symlink():
                manifest_artifact = attest_existing_artifact(
                    root, Artifact("manifest.json", manifest.sha256, len(manifest.as_bytes())),
                )
            else:
                manifest_artifact = (
                    _publish_or_reuse_small(root, "manifest.json", manifest.as_bytes())
                    if stage_state is not None else publish_manifest(root, manifest)
                )
            if stage_state is not None:
                stage_state.finish(
                    expected_partitions=tuple(
                        f"source-shard-{shard.ordinal:06d}" for shard in shards
                    ),
                    output=manifest_artifact,
                    process_id="process-" + os.urandom(16).hex(),
                    closure_check=lambda connection: (
                        connection.execute(
                            "SELECT records,physical,biological FROM qc_totals "
                            "WHERE singleton=1"
                        ).fetchone() == (
                            manifest.summary.total_records,
                            manifest.summary.total_physical_alternate_slots,
                            manifest.summary.total_biological_alternates,
                        )
                        and connection.execute(
                            "SELECT count(*) FROM class_seen"
                        ).fetchone() == (manifest.summary.total_records,)
                        and connection.execute(
                            "SELECT count(*) FROM normalization_records"
                        ).fetchone() == (manifest.summary.total_records,)
                        and connection.execute(
                            "SELECT count(*) FROM normalization_occurrences"
                        ).fetchone() == (
                            manifest.summary.total_physical_alternate_slots,
                        )
                    ),
                )
            return WgsSourceBuildResult(
                manifest, manifest_artifact, qc_summary,
                normalization_summary, class_summary,
            )


def verify_source_ledger(
    vcf_path: str | Path,
    *,
    concordance_report: ReferenceConcordanceReport | PublicValidationConcordanceReport,
    fetch_reference: Callable[[str, int, int], str],
    ledger_root: str | Path,
    manifest: LedgerManifest,
    selected_sample_name: str,
    selected_sample_sha256: str,
    case_manifest_sha256: str,
    executable_stack_sha256: str,
    execution_config: BoundedWgsExecutionConfig,
    allele_balance_screen: bool = False,
    audit_state_path: str | Path | None = None,
    stage_binding: StageRunBinding | None = None,
    resume: bool = False,
) -> SourceAccountingClosure:
    """Independently replay source bytes and canonical rows to close preparation."""

    if (
        not isinstance(manifest, LedgerManifest)
        or not isinstance(execution_config, BoundedWgsExecutionConfig)
        or type(allele_balance_screen) is not bool
        or not isinstance(selected_sample_name, str)
        or not selected_sample_name
        or any(not _digest(value) for value in (
            selected_sample_sha256, case_manifest_sha256, executable_stack_sha256,
        ))
        or execution_config.source_shards.max_uncompressed_rows > MAX_SHARD_ROWS
        or execution_config.source_shards.max_uncompressed_bytes > MAX_SHARD_BYTES
        or execution_config.input.max_record_bytes > MAX_ROW_BYTES
        or type(resume) is not bool
        or resume and stage_binding is None
        or stage_binding is not None and (
            type(stage_binding) is not StageRunBinding
            or stage_binding.stage != "source_audit"
            or audit_state_path is None
        )
        or audit_state_path is not None and (
            not Path(audit_state_path).is_absolute()
            or Path(audit_state_path).name != "source-audit.sqlite")
    ):
        raise WgsPreparationError("source ledger verification binding is invalid")
    root = Path(ledger_root)
    if not root.is_absolute():
        raise WgsPreparationError("source ledger root is invalid")
    if stage_binding is not None:
        return _verify_source_ledger_checkpointed(
            vcf_path, concordance_report=concordance_report,
            fetch_reference=fetch_reference, root=root, manifest=manifest,
            selected_sample_name=selected_sample_name,
            selected_sample_sha256=selected_sample_sha256,
            case_manifest_sha256=case_manifest_sha256,
            executable_stack_sha256=executable_stack_sha256,
            execution_config=execution_config,
            allele_balance_screen=allele_balance_screen,
            audit_state_path=Path(audit_state_path),
            stage_binding=stage_binding, resume=resume,
        )
    manifest_bytes = manifest.as_bytes()
    manifest_artifact = Artifact(
        "manifest.json", hashlib.sha256(manifest_bytes).hexdigest(), len(manifest_bytes),
    )
    if read_artifact(root, manifest_artifact) != manifest_bytes:
        raise WgsPreparationError("source ledger manifest closure is invalid")
    stream = iter_vcf_qc_records(
        vcf_path,
        concordance_report=concordance_report,
        max_record_bytes=execution_config.input.max_record_bytes,
        max_header_bytes=execution_config.input.max_header_bytes,
    )
    header = stream.header
    if selected_sample_name not in header.sample_names:
        stream.close()
        raise WgsPreparationError("selected sample is absent from the source header")
    execution_limits_sha256 = hashlib.sha256(
        bounded_wgs_execution_config_as_json(execution_config)
    ).hexdigest()
    expected_source = SourceBinding(
        header.input_sha256,
        header.input_size_bytes,
        case_manifest_sha256,
        selected_sample_sha256,
        header.reference_concordance_sha256,
        concordance_report.reference_manifest_sha256,
        executable_stack_sha256,
    )
    expected_parents = (
        ParentIdentity(
            "original_input", "genosapiens.original-input", 1,
            header.input_sha256,
        ),
        ParentIdentity(
            "case_manifest", "genosapiens.case-manifest", 1,
            case_manifest_sha256,
        ),
        ParentIdentity(
            "selected_sample", "genosapiens.selected-sample-identity", 1,
            selected_sample_sha256,
        ),
        ParentIdentity(
            "complete_concordance", concordance_report.schema, 1,
            header.reference_concordance_sha256,
        ),
        ParentIdentity(
            "reference_manifest", "genosapiens.public-resource-manifest", 1,
            concordance_report.reference_manifest_sha256,
        ),
        ParentIdentity(
            "execution_limits", execution_config.schema,
            execution_config.version, execution_limits_sha256,
        ),
        ParentIdentity(
            "executable_stack", "genosapiens.executable-stack", 1,
            executable_stack_sha256,
        ),
    )
    expected_header = canonical({
        "schema": HEADER_SCHEMA,
        "version": 1,
        "classification": "controlled",
        "llm_policy": "no-llm",
        "selected_sample_name": selected_sample_name,
        "header": asdict(header),
    })
    if (
        manifest.source != expected_source
        or manifest.parents != expected_parents
        or manifest.header.name != "header.json"
        or read_artifact(root, manifest.header) != expected_header
        or tuple(report.name for report in manifest.summary.reports)
        != ("qc.json", "normalization.json", "classes.json")
    ):
        stream.close()
        raise WgsPreparationError("source ledger parent closure is invalid")

    scratch_paths = tuple(
        root / name for name in (
            "verify-qc-state.sqlite",
            "verify-class-state.sqlite",
            "verify-normalization-state.sqlite",
        )
    )
    bounded_reference = _BoundedReferenceFetcher(
        fetch_reference, execution_config.input.max_record_bytes,
    )
    try:
        state_context = (
            _source_accounting_connection(Path(audit_state_path),
                cache_bytes=execution_config.sqlite.cache_bytes)
            if audit_state_path is not None else nullcontext(None)
        )
        with stage_lock(root, "source-verify"), stream, state_context as audit_connection:
            with (
                QcSummaryAccumulator(
                    header,
                    expected_records=concordance_report.total_records,
                    **({"connection": audit_connection} if audit_connection is not None
                       else {"scratch_path": scratch_paths[0]}),
                ) as qc_accumulator,
                IncrementalVariantInventoryAccumulator(
                    header,
                    expected_records=concordance_report.total_records,
                    allele_balance_screen=allele_balance_screen,
                    **({"connection": audit_connection} if audit_connection is not None
                       else {"scratch_path": scratch_paths[1]}),
                ) as class_accumulator,
                _NormalizationSummaryAccumulator(
                    None if audit_connection is not None else scratch_paths[2],
                    input_sha256=header.input_sha256,
                    concordance_sha256=header.reference_concordance_sha256,
                    expected_records=concordance_report.total_records,
                    cache_bytes=execution_config.sqlite.cache_bytes,
                    batch_size=execution_config.sqlite.cursor_batch_rows,
                    connection=audit_connection,
                ) as normalization_accumulator,
            ):
                rows = iter_ledger(root, manifest)
                for qc_record in stream:
                    try:
                        actual_row = next(rows)
                    except StopIteration:
                        raise WgsPreparationError(
                            "source ledger record coverage is invalid"
                        ) from None
                    bounded_reference.begin_record(qc_record.record_number)
                    normalization = normalize_source_record(
                        _source_record(qc_record),
                        input_sha256=header.input_sha256,
                        reference_identity_id=concordance_report.reference_identity_id,
                        reference_role=concordance_report.reference_role,
                        fetch_reference=bounded_reference,
                    )
                    qc_accumulator.add(qc_record)
                    inventory = class_accumulator.add(qc_record)
                    normalization_accumulator.add(normalization)
                    expected_row = LedgerRow(
                        record_number=qc_record.record_number,
                        source_record_id=qc_record.source_record_id,
                        source_allele_ids=tuple(
                            item.source_allele_id for item in normalization.alleles
                        ),
                        biological_alternate_count=_alternate_allele_count(
                            qc_record.alternates
                        ),
                        payload_json=_source_payload(
                            qc_record,
                            normalization,
                            inventory,
                            selected_sample_name=selected_sample_name,
                        ),
                    )
                    if actual_row != expected_row:
                        raise WgsPreparationError(
                            "source ledger semantic replay is invalid"
                        )
                try:
                    next(rows)
                except StopIteration:
                    pass
                else:
                    raise WgsPreparationError(
                        "source ledger record coverage is invalid"
                    )
                qc_summary = qc_accumulator.finalize()
                class_summary = class_accumulator.finalize()
                normalization_summary = normalization_accumulator.finalize()
                if audit_connection is not None:
                    audit_connection.commit()

            report_values = (qc_summary, normalization_summary, class_summary)
            for artifact, value in zip(
                manifest.summary.reports, report_values, strict=True,
            ):
                if read_artifact(root, artifact) != canonical(asdict(value)):
                    raise WgsPreparationError(
                        "source ledger summary replay is invalid"
                    )
            if manifest.summary != LedgerSummary(
                qc_summary.total_records,
                qc_summary.total_physical_alternate_slots,
                qc_summary.total_alternate_alleles,
                manifest.summary.reports,
            ):
                raise WgsPreparationError("source ledger summary closure is invalid")
            return SourceAccountingClosure(
                manifest.sha256,
                header.input_sha256,
                manifest.header.sha256,
                manifest.summary.reports[0].sha256,
                manifest.summary.reports[1].sha256,
                manifest.summary.reports[2].sha256,
                qc_summary.total_records,
                qc_summary.total_physical_alternate_slots,
                qc_summary.total_alternate_alleles,
            )
    finally:
        if not stream.closed:
            stream.close()
        for path in scratch_paths:
            try:
                path.unlink()
            except FileNotFoundError:
                pass
            except OSError:
                pass


def _verify_source_ledger_checkpointed(
    vcf_path, *, concordance_report, fetch_reference, root, manifest,
    selected_sample_name, selected_sample_sha256, case_manifest_sha256,
    executable_stack_sha256, execution_config, allele_balance_screen,
    audit_state_path, stage_binding, resume,
) -> SourceAccountingClosure:
    """Audit one source shard per durable cursor without replaying older shards."""
    limits_sha256 = hashlib.sha256(
        bounded_wgs_execution_config_as_json(execution_config)
    ).hexdigest()
    if (stage_binding.execution_config_sha256 != limits_sha256
            or stage_binding.input_sha256 != concordance_report.input_sha256
            or stage_binding.reference_sha256
            != concordance_report.reference_manifest_sha256
            or stage_binding.executable_stack_sha256 != executable_stack_sha256
            or manifest.source.input_sha256 != concordance_report.input_sha256):
        raise WgsPreparationError("source audit stage binding is invalid")
    manifest_bytes = manifest.as_bytes()
    if read_artifact(root, Artifact(
            "manifest.json", manifest.sha256, len(manifest_bytes))) != manifest_bytes:
        raise WgsPreparationError("source audit manifest changed")
    _private_directory(audit_state_path.parent)
    with ExitStack() as contexts:
        lock = contexts.enter_context(stage_lock(root, "source_audit"))
        stage_state = contexts.enter_context(StageCommitState(
            root, audit_state_path, stage_binding, lock=lock, create=not resume,
        ))
        connection = stage_state.connection
        if resume:
            tail = connection.execute(
                "SELECT event_bytes FROM events ORDER BY ordinal DESC LIMIT 1"
            ).fetchone()
            if tail is not None:
                completion = load_model(tail[0], JournalEntry)
                if completion.event_type == "stage_completed":
                    if len(completion.artifacts) != 1:
                        raise WgsPreparationError("source audit completion is invalid")
                    closure = load_source_accounting(
                        read_artifact(root, completion.artifacts[0]))
                    _require_source_accounting(manifest, closure)
                    return closure
            try:
                saved = connection.execute(
                    "SELECT identity_bytes FROM audit_input WHERE singleton=1"
                ).fetchone()
                if saved is None:
                    raise ValueError
                identity = json.loads(saved[0])
                if (type(identity) is not list or len(identity) != 5
                        or any(type(value) is not int or value < 0 for value in identity)):
                    raise ValueError
                resume_identity = tuple(identity)
                committed = connection.execute(
                    "SELECT partition_id,cursor_bytes FROM commits ORDER BY ordinal"
                ).fetchall()
                if len(committed) > len(manifest.shards):
                    raise ValueError
                for ordinal, (partition_id, cursor_bytes) in enumerate(committed, 1):
                    shard = manifest.shards[ordinal - 1]
                    if (partition_id != f"audit-shard-{ordinal:06d}"
                            or json.loads(cursor_bytes) != [shard.last_record, ordinal]):
                        raise ValueError
                if committed:
                    latest = stage_state.committed(
                        f"audit-shard-{len(committed):06d}")
                    if latest is None or latest.cursor != (
                            manifest.shards[len(committed) - 1].last_record,
                            len(committed)):
                        raise ValueError
                committed_shards = len(committed)
            except (sqlite3.Error, ValueError, TypeError, IndexError, KeyError):
                raise WgsPreparationError("source audit checkpoint is inconsistent") from None
        else:
            resume_identity = None
            committed_shards = 0
        skip_records = (
            manifest.shards[committed_shards - 1].last_record
            if committed_shards else 0
        )
        stream = contexts.enter_context(iter_vcf_qc_records(
            vcf_path, concordance_report=concordance_report,
            max_record_bytes=execution_config.input.max_record_bytes,
            max_header_bytes=execution_config.input.max_header_bytes,
            resume_input_identity=resume_identity, skip_records=skip_records,
        ))
        header = stream.header
        expected_source = SourceBinding(
            header.input_sha256, header.input_size_bytes, case_manifest_sha256,
            selected_sample_sha256, header.reference_concordance_sha256,
            concordance_report.reference_manifest_sha256, executable_stack_sha256,
        )
        expected_parents = (
            ParentIdentity("original_input", "genosapiens.original-input", 1,
                           header.input_sha256),
            ParentIdentity("case_manifest", "genosapiens.case-manifest", 1,
                           case_manifest_sha256),
            ParentIdentity("selected_sample", "genosapiens.selected-sample-identity", 1,
                           selected_sample_sha256),
            ParentIdentity("complete_concordance", concordance_report.schema, 1,
                           header.reference_concordance_sha256),
            ParentIdentity("reference_manifest", "genosapiens.public-resource-manifest", 1,
                           concordance_report.reference_manifest_sha256),
            ParentIdentity("execution_limits", execution_config.schema,
                           execution_config.version, limits_sha256),
            ParentIdentity("executable_stack", "genosapiens.executable-stack", 1,
                           executable_stack_sha256),
        )
        expected_header = canonical({
            "schema": HEADER_SCHEMA, "version": 1, "classification": "controlled",
            "llm_policy": "no-llm", "selected_sample_name": selected_sample_name,
            "header": asdict(header),
        })
        if (selected_sample_name not in header.sample_names
                or manifest.source != expected_source
                or manifest.parents != expected_parents
                or read_artifact(root, manifest.header) != expected_header
                or tuple(report.name for report in manifest.summary.reports)
                != ("qc.json", "normalization.json", "classes.json")):
            raise WgsPreparationError("source audit parent closure is invalid")
        if not resume:
            connection.execute(
                "CREATE TABLE audit_input(singleton INTEGER PRIMARY KEY CHECK(singleton=1),"
                "identity_bytes BLOB NOT NULL)"
            )
            connection.execute(
                "INSERT INTO audit_input VALUES (1,?)",
                (canonical(stream._input_identity),),
            )
        bounded_reference = _BoundedReferenceFetcher(
            fetch_reference, execution_config.input.max_record_bytes,
            first_record=skip_records + 1,
        )
        with (
            QcSummaryAccumulator(
                header, expected_records=concordance_report.total_records,
                connection=connection, resume=resume,
            ) as qc_accumulator,
            IncrementalVariantInventoryAccumulator(
                header, expected_records=concordance_report.total_records,
                allele_balance_screen=allele_balance_screen,
                connection=connection, resume=resume,
            ) as class_accumulator,
            _NormalizationSummaryAccumulator(
                input_sha256=header.input_sha256,
                concordance_sha256=header.reference_concordance_sha256,
                expected_records=concordance_report.total_records,
                cache_bytes=execution_config.sqlite.cache_bytes,
                batch_size=execution_config.sqlite.cursor_batch_rows,
                connection=connection, resume=resume,
            ) as normalization_accumulator,
        ):
            connection.commit()
            for shard in manifest.shards[committed_shards:]:
                partition_id = f"audit-shard-{shard.ordinal:06d}"
                if resume:
                    incomplete = stage_state.open_attempt(partition_id)
                    if incomplete is not None:
                        stage_state.archive_orphans(
                            names=(f"audit-shard-{shard.ordinal:06d}.json",),
                            attempt_id=incomplete.attempt_id,
                        )
                        stage_state.fail(
                            incomplete, failure_code="process_interrupted",
                            interrupted=True,
                        )
                start = _begin_source_checkpoint_attempt(stage_state, partition_id, resume=resume)
                connection.execute("BEGIN IMMEDIATE")
                for actual_row in iter_shard(root, shard):
                    try:
                        qc_record = next(stream)
                    except StopIteration:
                        raise WgsPreparationError(
                            "source audit record coverage is invalid") from None
                    bounded_reference.begin_record(qc_record.record_number)
                    normalization = normalize_source_record(
                        _source_record(qc_record),
                        input_sha256=header.input_sha256,
                        reference_identity_id=concordance_report.reference_identity_id,
                        reference_role=concordance_report.reference_role,
                        fetch_reference=bounded_reference,
                    )
                    qc_accumulator.add(qc_record)
                    inventory = class_accumulator.add(qc_record)
                    normalization_accumulator.add(normalization)
                    expected_row = LedgerRow(
                        qc_record.record_number, qc_record.source_record_id,
                        tuple(item.source_allele_id for item in normalization.alleles),
                        _alternate_allele_count(qc_record.alternates),
                        _source_payload(qc_record, normalization, inventory,
                                        selected_sample_name=selected_sample_name),
                    )
                    if actual_row != expected_row:
                        raise WgsPreparationError("source ledger semantic replay is invalid")
                marker = _publish_or_reuse_small(
                    root, f"audit-shard-{shard.ordinal:06d}.json", canonical({
                        "schema": "genosapiens.wgs-source-audit-shard", "version": 1,
                        "source_manifest_sha256": manifest.sha256,
                        "shard_sha256": shard.sha256,
                        "last_record": shard.last_record,
                    }),
                )
                stage_state.commit_partition(
                    start, cursor=(shard.last_record, shard.ordinal),
                    outputs=(marker,), counts={"records": shard.rows},
                    reducer=lambda _connection: (
                        qc_accumulator.flush(), class_accumulator.flush(),
                        normalization_accumulator.flush(),
                    ),
                    prepared_transaction=True,
                )
            try:
                next(stream)
            except StopIteration:
                pass
            else:
                raise WgsPreparationError("source audit record coverage is invalid")
            qc_summary = qc_accumulator.finalize()
            class_summary = class_accumulator.finalize()
            normalization_summary = normalization_accumulator.finalize()
            connection.commit()
        for artifact, value in zip(manifest.summary.reports,
                                   (qc_summary, normalization_summary, class_summary),
                                   strict=True):
            if read_artifact(root, artifact) != canonical(asdict(value)):
                raise WgsPreparationError("source ledger summary replay is invalid")
        if manifest.summary != LedgerSummary(
                qc_summary.total_records, qc_summary.total_physical_alternate_slots,
                qc_summary.total_alternate_alleles, manifest.summary.reports):
            raise WgsPreparationError("source ledger summary closure is invalid")
        closure = SourceAccountingClosure(
            manifest.sha256, header.input_sha256, manifest.header.sha256,
            *(item.sha256 for item in manifest.summary.reports),
            qc_summary.total_records, qc_summary.total_physical_alternate_slots,
            qc_summary.total_alternate_alleles,
        )
        closure_artifact = _publish_or_reuse_small(
            root, "audit-closure.json", closure.as_bytes(),
        )
        stage_state.finish(
            expected_partitions=tuple(
                f"audit-shard-{shard.ordinal:06d}" for shard in manifest.shards
            ), output=closure_artifact,
            process_id="process-" + os.urandom(16).hex(),
            closure_check=lambda checked: checked.execute(
                "SELECT records,physical,biological FROM qc_totals WHERE singleton=1"
            ).fetchone() == (
                closure.total_records, closure.total_physical_alternate_slots,
                closure.total_biological_alternates,
            ),
        )
        return closure


PREPARATION_RECEIPT_SCHEMA = "genosapiens.wgs-bounded-preparation-result"
PREPARATION_RECEIPT_NAME = "controlled-preparation-result.json"
PREPARATION_PLAN_SCHEMA = "genosapiens.wgs-preparation-partition-plan"
PREPARATION_STAGE_PARENT_SCHEMA = "genosapiens.wgs-preparation-stage-parent"



@dataclass(frozen=True, slots=True, repr=False)
class PreparationFilterContext:
    """The complete early-filter input surface; it deliberately has no HPO."""

    source_root: Path
    source_manifest: LedgerManifest
    source_accounting: SourceAccountingClosure
    policy_sha256: str
    population_manifest_sha256: str
    clinical_manifest_sha256: str
    exception_manifest_sha256: str
    precision_manifest_sha256: str
    scope_policy_sha256: str
    scope_manifest_sha256: str
    completed_source: ResolvedCompletedSource | None = None
    source_cursor: SourceDecisionCursor | None = None

    def __post_init__(self) -> None:
        if (
            not self.source_root.is_absolute()
            or self.completed_source is not None and (
                type(self.completed_source) is not ResolvedCompletedSource
                or self.source_root != self.completed_source.root
                or self.source_manifest != self.completed_source.manifest)
            or self.source_manifest.sha256
            != self.source_accounting.source_manifest_sha256
            or any(not _digest(value) for value in (
                self.policy_sha256,
                self.population_manifest_sha256,
                self.clinical_manifest_sha256,
                self.exception_manifest_sha256,
                self.precision_manifest_sha256,
                self.scope_policy_sha256, self.scope_manifest_sha256,
            ))
        ):
            raise WgsPreparationError("preparation filter context is invalid")


class PreparationDecisionProvider(Protocol):
    def __call__(
        self, context: PreparationFilterContext,
    ) -> Iterable["FilterDecision"]: ...


class SourceRoutingContextProvider(Protocol):
    def __call__(self, context: PreparationFilterContext) -> "SourceRoutingContext": ...


@dataclass(frozen=True, slots=True, repr=False)
class PreparationRunReceipt:
    status: Literal["passed", "failed", "interrupted"]
    terminal_stage: str
    base_release_sha256: str | None
    prepared_case_sha256: str | None
    early_snapshot_sha256: str | None
    completed_journal_heads: tuple[str, ...]
    failure_code: str | None
    preparation_run_options_sha256: str
    audit_requested: bool
    source_accounting_sha256: str | None
    source_audit_sha256: str | None
    prepared_index_sha256: str | None
    completed_source_input_sha256: str | None
    source_origin: Literal["produced_this_run", "completed_input"]
    public_validation_authority_sha256: str | None = None
    schema: Literal["genosapiens.wgs-bounded-preparation-result"] = (
        PREPARATION_RECEIPT_SCHEMA
    )
    version: Literal[3] = 3
    classification: Literal["controlled"] = "controlled"
    llm_policy: Literal["no-llm"] = "no-llm"
    competitive_readiness_claimed: Literal[False] = False

    def __post_init__(self) -> None:
        passed = self.status == "passed"
        if (
            self.schema != PREPARATION_RECEIPT_SCHEMA
            or self.version != 3
            or self.source_origin not in {"produced_this_run", "completed_input"}
            or self.completed_source_input_sha256 is not None
            and not _digest(self.completed_source_input_sha256)
            or passed and self.completed_source_input_sha256 is None
            or self.classification != "controlled"
            or self.llm_policy != "no-llm"
            or self.competitive_readiness_claimed is not False
            or self.status not in {"passed", "failed", "interrupted"}
            or type(self.terminal_stage) is not str
            or not self.terminal_stage
            or (self.base_release_sha256 is None) == (self.public_validation_authority_sha256 is None)
            or not _digest(self.base_release_sha256 or self.public_validation_authority_sha256)
            or any(not _digest(value) for value in self.completed_journal_heads)
            or not _digest(self.preparation_run_options_sha256)
            or type(self.audit_requested) is not bool
            or self.source_accounting_sha256 is not None
            and not _digest(self.source_accounting_sha256)
            or self.source_audit_sha256 is not None
            and not _digest(self.source_audit_sha256)
            or self.prepared_index_sha256 is not None
            and not _digest(self.prepared_index_sha256)
            or self.source_audit_sha256 is not None
            and self.source_accounting_sha256 is None
            or passed and (self.source_accounting_sha256 is None
                           or self.source_audit_sha256 is None
                           or self.prepared_index_sha256 is None)
            or passed != (self.prepared_case_sha256 is not None)
            or passed != (self.early_snapshot_sha256 is not None)
            or passed == (self.failure_code is not None)
            or self.prepared_case_sha256 is not None
            and not _digest(self.prepared_case_sha256)
            or self.early_snapshot_sha256 is not None
            and not _digest(self.early_snapshot_sha256)
            or self.failure_code is not None
            and self.failure_code not in {
                "preparation_failed", "preparation_interrupted",
            }
        ):
            raise WgsPreparationError("preparation result receipt is invalid")

    def as_bytes(self) -> bytes:
        return canonical(asdict(self))


@dataclass(frozen=True, slots=True, repr=False)
class BoundedPreparationResult:
    root: Path
    receipt: PreparationRunReceipt
    receipt_artifact: Artifact
    prepared_case: "PreparedCasePublication | None"

    def __post_init__(self) -> None:
        if (
            not self.root.is_absolute()
            or self.receipt_artifact.name != PREPARATION_RECEIPT_NAME
            or self.receipt_artifact.sha256
            != hashlib.sha256(self.receipt.as_bytes()).hexdigest()
            or (self.receipt.status == "passed")
            != (self.prepared_case is not None)
            or self.prepared_case is not None and (
                self.receipt.prepared_case_sha256 != self.prepared_case.manifest.sha256
                or self.receipt.early_snapshot_sha256
                != self.prepared_case.early_snapshot.sha256
                or self.receipt.prepared_index_sha256
                != self.prepared_case.prepared_index.sha256
                or self.receipt.source_accounting_sha256
                != self.prepared_case.manifest.source_accounting.sha256
                or self.receipt.source_audit_sha256
                != self.prepared_case.manifest.source_audit.sha256
                or self.receipt.preparation_run_options_sha256
                != self.prepared_case.manifest.preparation_run_options.sha256
                or self.receipt.audit_requested
                != self.prepared_case.diagnostic.audit_requested
            )
        ):
            raise WgsPreparationError("bounded preparation result is invalid")


class _PreparationInterrupted(Exception):
    pass


def _private_directory(path: Path) -> None:
    try:
        if not path.exists():
            path.mkdir(mode=0o700, parents=True)
        metadata = path.lstat()
        if (
            not path.is_absolute()
            or path.is_symlink()
            or not stat.S_ISDIR(metadata.st_mode)
            or stat.S_IMODE(metadata.st_mode) != 0o700
            or metadata.st_uid != os.getuid()
        ):
            raise ValueError
    except (OSError, ValueError):
        raise WgsPreparationError("private preparation directory is unsafe") from None


def _authenticate_hpo_input(
    path: Path, *, expected_sha256: str, expected_size_bytes: int,
) -> None:
    try:
        metadata = path.lstat()
        if (
            not path.is_absolute()
            or path.is_symlink()
            or not stat.S_ISREG(metadata.st_mode)
            or metadata.st_size != expected_size_bytes
            or metadata.st_uid != os.getuid()
            or stat.S_IMODE(metadata.st_mode) != 0o600
            or not _digest(expected_sha256)
            or not _count(expected_size_bytes, 1)
        ):
            raise ValueError
        digest = hashlib.sha256()
        with path.open("rb") as stream:
            for block in iter(lambda: stream.read(1024 * 1024), b""):
                digest.update(block)
        if digest.hexdigest() != expected_sha256:
            raise ValueError
    except (OSError, ValueError):
        raise WgsPreparationError("HPO input authentication failed") from None


def _preparation_stage_parent(
    stage: str,
    *,
    case_manifest_sha256: str,
    selected_sample_sha256: str,
    preceding_sha256: str,
    execution_limits_sha256: str,
    executable_stack_sha256: str,
    runtime_image_digest: str,
    run_options_sha256: str,
    completed_source_input_sha256: str | None = None,
) -> str:
    value = {
        "schema": PREPARATION_STAGE_PARENT_SCHEMA,
        "version": 1,
        "stage": stage,
        "case_manifest_sha256": case_manifest_sha256,
        "selected_sample_sha256": selected_sample_sha256,
        "preceding_sha256": preceding_sha256,
        "execution_limits_sha256": execution_limits_sha256,
        "executable_stack_sha256": executable_stack_sha256,
        "runtime_image_digest": runtime_image_digest,
        "run_options_sha256": run_options_sha256,
    }
    if completed_source_input_sha256 is not None:
        if stage != "prefilter" or not _digest(completed_source_input_sha256):
            raise WgsPreparationError("completed source stage parent is invalid")
        value["completed_source_input_sha256"] = completed_source_input_sha256
    return hashlib.sha256(canonical(value)).hexdigest()


def _begin_source_checkpoint_attempt(state, partition_id, *, resume):
    """Allow one technical attempt in this explicit invocation, under its lock."""
    starts = state.connection.execute(
        "SELECT COUNT(*) FROM events WHERE event_type='attempt_started' AND partition_id=?",
        (partition_id,),
    ).fetchone()[0]
    return state.begin(partition_id=partition_id, scientific=False,
                       max_attempts=starts + 1 if resume else 1)


def _execute_preparation_stage(
    *,
    journal_root: Path,
    stage: str,
    partition_id: str,
    stage_parent_sha256: str,
    execution_config: BoundedWgsExecutionConfig,
    execution_limits_sha256: str,
    executable_stack_sha256: str,
    runtime_image_digest: str,
    scientific: bool,
    action: Callable[[], tuple[object, str]],
    resume_action: Callable[[str], tuple[object, str]],
    resume: bool,
    prepare_retry: Callable[[int], None] | None = None,
    input_parents: tuple[ParentIdentity, ...] = (),
    technical_classification_resume: bool = False,
):
    from genosapiens.wgs_ledger import (
        begin_attempt,
        complete_attempt,
        complete_stage,
        fail_attempt,
    )
    from genosapiens.wgs_prepared_case import CompletedJournalBinding

    _private_directory(journal_root)
    plan_content = canonical({
        "schema": PREPARATION_PLAN_SCHEMA,
        "version": 1,
        "stage": stage,
        "partitions": (partition_id,),
        "complete_source_required": True,
        "classification": "controlled",
        "llm_policy": "no-llm",
    })
    plan_expected = Artifact(
        "partition-plan.json", hashlib.sha256(plan_content).hexdigest(),
        len(plan_content),
    )
    plan_path = journal_root / plan_expected.name
    if plan_path.exists():
        if not resume or read_artifact(journal_root, plan_expected) != plan_content:
            raise WgsPreparationError("preparation resume plan is invalid")
        plan = plan_expected
    else:
        plan = publish_bytes(journal_root, plan_expected.name, plan_content)
    parents = (
        ParentIdentity(
            "stage_parent", PREPARATION_STAGE_PARENT_SCHEMA, 1,
            stage_parent_sha256,
        ),
        ParentIdentity(
            "partition_plan", PREPARATION_PLAN_SCHEMA, 1, plan.sha256,
        ),
        ParentIdentity(
            "execution_limits", execution_config.schema,
            execution_config.version, execution_limits_sha256,
        ),
        ParentIdentity(
            "executable_stack", "genosapiens.executable-stack", 1,
            executable_stack_sha256,
        ),
        ParentIdentity(
            "runtime_image", "genosapiens.runtime-image", 1,
            runtime_image_digest.removeprefix("sha256:"),
        ),
    )
    parents += input_parents
    with stage_lock(journal_root, stage) as lock:
        entries = read_journal(
            journal_root, stage=stage, parent_sha256=stage_parent_sha256,
        )
        if entries and not resume:
            raise WgsPreparationError("preparation stage requires explicit resume")
        if entries and entries[0].parents != parents:
            raise WgsPreparationError("preparation resume parents changed")
        terminals = {
            entry.attempt_id for entry in entries
            if entry.event_type in {
                "partition_completed", "attempt_failed", "attempt_interrupted",
            }
        }
        open_starts = [
            entry for entry in entries
            if entry.event_type == "attempt_started"
            and entry.attempt_id not in terminals
        ]
        for open_start in open_starts:
            fail_attempt(
                journal_root,
                lock=lock,
                start=open_start,
                failure_code="process_lost",
                interrupted=True,
            )
        if open_starts:
            entries = read_journal(
                journal_root, stage=stage, parent_sha256=stage_parent_sha256,
            )
        stage_events = [
            entry for entry in entries if entry.event_type == "stage_completed"
        ]
        completed = [
            entry for entry in entries
            if entry.event_type == "partition_completed"
            and entry.partition_id == partition_id
        ]
        if stage_events:
            if (
                not resume
                or len(stage_events) != 1
                or stage_events[0] != entries[-1]
                or len(completed) != 1
            ):
                raise WgsPreparationError("preparation completed journal is invalid")
            stage_output = read_artifact(
                journal_root, stage_events[0].artifacts[0],
            )
            try:
                record = json.loads(stage_output)
                output_sha256 = record["output_sha256"]
            except (UnicodeError, json.JSONDecodeError, KeyError, TypeError):
                raise WgsPreparationError("preparation completed output is invalid") from None
            if (
                record.get("schema") != "genosapiens.wgs-preparation-stage-output"
                or record.get("version") != 1
                or record.get("stage") != stage
                or record.get("status") != "complete"
                or not _digest(output_sha256)
            ):
                raise WgsPreparationError("preparation completed output is invalid")
            result, actual_sha256 = resume_action(output_sha256)
            if actual_sha256 != output_sha256:
                raise WgsPreparationError("preparation resumed output changed")
            return result, CompletedJournalBinding(
                stage, journal_root, stage_parent_sha256, stage_events[0].sha256,
            )
        if completed:
            if not resume or len(completed) != 1:
                raise WgsPreparationError("preparation partition journal is invalid")
            partition_output = read_artifact(
                journal_root, completed[0].artifacts[0],
            )
            try:
                record = json.loads(partition_output)
                output_sha256 = record["output_sha256"]
            except (UnicodeError, json.JSONDecodeError, KeyError, TypeError):
                raise WgsPreparationError("preparation partition output is invalid") from None
            if (
                record.get("schema")
                != "genosapiens.wgs-preparation-partition-output"
                or record.get("version") != 1
                or record.get("stage") != stage
                or record.get("partition_id") != partition_id
                or record.get("status") != "complete"
                or not _digest(output_sha256)
            ):
                raise WgsPreparationError("preparation partition output is invalid")
            result, actual_sha256 = resume_action(output_sha256)
            if actual_sha256 != output_sha256:
                raise WgsPreparationError("preparation resumed output changed")
            stage_content = canonical({
                "schema": "genosapiens.wgs-preparation-stage-output",
                "version": 1,
                "stage": stage,
                "partition_output_sha256": completed[0].artifacts[0].sha256,
                "output_sha256": output_sha256,
                "status": "complete",
            })
            expected_stage_output = Artifact(
                "stage-output.json", hashlib.sha256(stage_content).hexdigest(),
                len(stage_content),
            )
            stage_path = journal_root / expected_stage_output.name
            stage_output_artifact = (
                expected_stage_output
                if stage_path.exists()
                and read_artifact(journal_root, expected_stage_output) == stage_content
                else publish_bytes(journal_root, "stage-output.json", stage_content)
            )
            head = complete_stage(
                journal_root,
                lock=lock,
                stage=stage,
                parent_sha256=stage_parent_sha256,
                parents=parents,
                output=stage_output_artifact,
            )
            return result, CompletedJournalBinding(
                stage, journal_root, stage_parent_sha256, head.sha256,
            )
        starts_used = sum(
            entry.event_type == "attempt_started"
            and entry.partition_id == partition_id
            for entry in entries
        )
        maximum_attempts = (
            execution_config.journal.default_scientific_attempts_per_partition
            if scientific else 1
        )
        if resume and not scientific and stage in {'source', 'source_audit'}:
            maximum_attempts = starts_used + 1
        if (resume and stage == "prefilter" and technical_classification_resume
                and any(entry.event_type == "attempt_failed"
                        and entry.failure_code == "source_classification_malformed_evidence"
                        for entry in entries)):
            maximum_attempts = max(maximum_attempts, starts_used + 1)
        if starts_used >= maximum_attempts:
            raise WgsPreparationError("preparation attempt budget is exhausted")
        if starts_used and prepare_retry is not None:
            prepare_retry(starts_used)
        start = begin_attempt(
            journal_root,
            lock=lock,
            stage=stage,
            partition_id=partition_id,
            parent_sha256=stage_parent_sha256,
            parents=parents,
            scientific=scientific,
            max_attempts=maximum_attempts,
        )
        try:
            result, output_sha256 = action()
            if not _digest(output_sha256):
                raise WgsPreparationError("preparation stage output is invalid")
        except (Exception, KeyboardInterrupt) as error:
            operation = getattr(error, "_wgs_source_operation", None)
            if (stage not in {"source", "source_audit"} or type(operation) is not str
                    or operation not in SOURCE_OPERATIONS):
                operation = "preparation_stage_action"
            failure = static_failure(
                error, stage, operation=operation,
            )
            interrupted_failure = (
                isinstance(error, KeyboardInterrupt)
                or failure.category == "interrupted"
            )
            try:
                head = fail_attempt(
                    journal_root, lock=lock, start=start,
                    failure_code=failure.failure_code,
                    interrupted=interrupted_failure,
                )
                failure.journal_heads = (head.sha256,)
            except Exception as journal_error:
                # The durable start still consumes the attempt if recording fails.
                failure.add_secondary(static_failure(
                    journal_error, stage, operation="journal_publish",
                ))
            if interrupted_failure:
                interrupted = _PreparationInterrupted()
                interrupted.failure = failure
                raise interrupted from None
            raise failure from None
        partition_output = publish_bytes(
            journal_root,
            "partition-output.json",
            canonical({
                "schema": "genosapiens.wgs-preparation-partition-output",
                "version": 1,
                "stage": stage,
                "partition_id": partition_id,
                "output_sha256": output_sha256,
                "status": "complete",
            }),
        )
        complete_attempt(
            journal_root, lock=lock, start=start, output=partition_output,
        )
        stage_output = publish_bytes(
            journal_root,
            "stage-output.json",
            canonical({
                "schema": "genosapiens.wgs-preparation-stage-output",
                "version": 1,
                "stage": stage,
                "partition_output_sha256": partition_output.sha256,
                "output_sha256": output_sha256,
                "status": "complete",
            }),
        )
        head = complete_stage(
            journal_root,
            lock=lock,
            stage=stage,
            parent_sha256=stage_parent_sha256,
            parents=parents,
            output=stage_output,
        )
    return result, CompletedJournalBinding(
        stage, journal_root, stage_parent_sha256, head.sha256,
    )


def _publish_preparation_receipt(
    root: Path,
    *,
    status: Literal["passed", "failed", "interrupted"],
    terminal_stage: str,
    base_release_sha256: str | None,
    journals: tuple[object, ...],
    preparation_run_options_sha256: str,
    audit_requested: bool,
    source_accounting_sha256: str | None,
    source_audit_sha256: str | None,
    completed_source_input_sha256: str | None,
    source_origin: Literal["produced_this_run", "completed_input"],
    prepared_index_sha256: str | None = None,
    prepared_case: "PreparedCasePublication | None" = None,
    failure=None,
    public_validation_authority_sha256: str | None = None,
) -> BoundedPreparationResult:
    receipt = PreparationRunReceipt(
        status,
        terminal_stage,
        base_release_sha256,
        prepared_case.manifest.sha256 if prepared_case is not None else None,
        prepared_case.early_snapshot.sha256 if prepared_case is not None else None,
        tuple(value.head_sha256 for value in journals) + (
            () if failure is None else tuple(failure.journal_heads)
        ),
        None if status == "passed" else (
            "preparation_interrupted" if status == "interrupted"
            else "preparation_failed"
        ),
        preparation_run_options_sha256,
        audit_requested,
        source_accounting_sha256,
        source_audit_sha256,
        prepared_index_sha256,
        completed_source_input_sha256,
        source_origin,
        public_validation_authority_sha256=(None if base_release_sha256 is not None
                                            else public_validation_authority_sha256),
    )
    try:
        content = receipt.as_bytes()
        artifact = Artifact(
            PREPARATION_RECEIPT_NAME, hashlib.sha256(content).hexdigest(),
            len(content),
        )
        if (root / PREPARATION_RECEIPT_NAME).exists():
            if read_artifact(root, artifact) != content:
                raise WgsPreparationError("preparation receipt changed on resume")
        else:
            artifact = publish_bytes(root, PREPARATION_RECEIPT_NAME, content)
    except (Exception, KeyboardInterrupt) as error:
        if failure is not None:
            failure.add_secondary(static_failure(
                error, "final_publication", operation="receipt_publish",
            ))
            raise failure from None
        raise
    if failure is not None:
        publish_failure_evidence(
            root, stage="preparation", receipt_content=receipt.as_bytes(), failure=failure,
        )
    return BoundedPreparationResult(root, receipt, artifact, prepared_case)


def _preserve_retry_directories(
    root: Path,
    *,
    stage: str,
    attempt_number: int,
    directories: tuple[Path, ...],
) -> None:
    """Move an incomplete attempt aside before any allowed retry mutates paths."""

    history_root = root / "resume-history"
    _private_directory(history_root)
    history = history_root / f"{stage}-attempt-{attempt_number:02d}"
    if history.exists() or history.is_symlink():
        raise WgsPreparationError("preparation retry history already exists")
    _private_directory(history)
    records = []
    try:
        for source in directories:
            metadata = source.lstat()
            if (
                source.is_symlink()
                or not stat.S_ISDIR(metadata.st_mode)
                or stat.S_IMODE(metadata.st_mode) != 0o700
                or metadata.st_uid != os.getuid()
            ):
                raise ValueError
            destination = history / source.name
            os.rename(source, destination)
            _private_directory(source)
            for path in sorted(destination.rglob("*")):
                relative = path.relative_to(history).as_posix()
                item = path.lstat()
                if path.is_symlink() or item.st_uid != os.getuid():
                    raise ValueError
                if path.is_dir():
                    if stat.S_IMODE(item.st_mode) != 0o700:
                        raise ValueError
                    records.append({"path": relative, "kind": "directory"})
                elif (
                    stat.S_ISREG(item.st_mode)
                    and stat.S_IMODE(item.st_mode) == 0o600
                    and item.st_nlink == 1
                ):
                    digest = hashlib.sha256()
                    with path.open("rb") as stream:
                        for block in iter(lambda: stream.read(1024 * 1024), b""):
                            digest.update(block)
                    records.append({
                        "path": relative,
                        "kind": "file",
                        "sha256": digest.hexdigest(),
                        "size_bytes": item.st_size,
                    })
                else:
                    raise ValueError
        publish_bytes(history, "attempt-evidence.json", canonical({
            "schema": "genosapiens.wgs-preparation-retry-evidence",
            "version": 1,
            "classification": "controlled",
            "llm_policy": "no-llm",
            "stage": stage,
            "attempt_number": attempt_number,
            "products": records,
        }))
    except (OSError, ValueError):
        raise WgsPreparationError("preparation partial attempt is unsafe") from None


def prepare_bounded_case(
    *,
    output_root: str | Path,
    original_input_path: str | Path,
    hpo_input_path: str | Path,
    hpo_input_sha256: str,
    hpo_input_size_bytes: int,
    hpo_ontology: ParsedOntology,
    hpo_reference_root: str | Path,
    concordance_report: ReferenceConcordanceReport | PublicValidationConcordanceReport,
    fetch_reference: Callable[[str, int, int], str],
    selected_sample_name: str,
    selected_sample_sha256: str,
    case_manifest_sha256: str,
    ontology_manifest_sha256: str,
    executable_stack_sha256: str,
    runtime_image_digest: str,
    execution_config: BoundedWgsExecutionConfig,
    run_options: PreparationRunOptions,
    policy_sha256: str,
    population_manifest_sha256: str,
    clinical_manifest_sha256: str,
    exception_manifest_sha256: str,
    precision_manifest_sha256: str,
    scope_policy_sha256: str,
    scope_manifest_sha256: str,
    decision_provider: PreparationDecisionProvider,
    source_routing_context_provider: SourceRoutingContextProvider,
    base_release_bytes: bytes | None,
    expected_base_release_sha256: str | None,
    allele_balance_screen: bool = True,
    source_worker_reference: SourceWorkerReference | None = None,
    resume: bool = False,
    public_qualification_authority_bytes: bytes | None = None,
    expected_public_qualification_authority_sha256: str | None = None,
    public_validation_authority_bytes: bytes | None = None,
    expected_public_validation_authority_sha256: str | None = None,
    phase_recorder: RuntimePhaseRecorder | None = None,
    completed_source_input: CompletedSourceArtifact | None = None,
    operator_case_binding=None,
) -> BoundedPreparationResult:
    """Run the bounded preparation DAG and publish one immutable static receipt."""

    from genosapiens.wgs_filter import (
        FilterPublication,
        WorkingSetPublication,
        SourceDecisionCursor,
        _FilterState,
        _check_filter_read_boundary,
        _working_parents,
        publish_filter_decisions,
        publish_working_set,
    )
    from genosapiens.wgs_prepared_case import (
        PreparationDiagnosticReceipt,
        _source_accounting_sha256,
        build_prepared_index,
        publish_prepared_index,
        publish_prepared_case,
    )
    from genosapiens.wgs_source_routing import SourceRoutingContext

    root = Path(output_root)
    if (
        not root.is_absolute()
        or type(resume) is not bool
        or type(run_options) is not PreparationRunOptions
        or (source_worker_reference is not None
            and type(source_worker_reference) is not SourceWorkerReference)
        or phase_recorder is not None
        and type(phase_recorder) is not RuntimePhaseRecorder
    ):
        raise WgsPreparationError("bounded preparation root must be absolute")
    _private_directory(root)
    journals = []
    stage = "input_authentication"
    options_sha256 = run_options.sha256
    accounting_sha256 = None
    audit_result_sha256 = None
    prepared_index_sha256 = None
    source_input_sha256 = None
    source_origin = "completed_input" if completed_source_input is not None else "produced_this_run"
    try:
        original = Path(original_input_path)
        hpo = Path(hpo_input_path)
        with runtime_phase(phase_recorder, "input_validation"):
            _authenticate_hpo_input(
                hpo,
                expected_sha256=hpo_input_sha256,
                expected_size_bytes=hpo_input_size_bytes,
            )
        if (
            not original.is_absolute()
            or not callable(fetch_reference)
            or not callable(decision_provider)
            or not callable(source_routing_context_provider)
            or not (
                isinstance(concordance_report, ReferenceConcordanceReport)
                or (public_qualification_authority_bytes is not None or public_validation_authority_bytes is not None)
                and isinstance(
                    concordance_report, PublicValidationConcordanceReport,
                )
            )
            or not concordance_report.accepted
            or any(not _digest(value) for value in (
                selected_sample_sha256,
                case_manifest_sha256,
                ontology_manifest_sha256,
                executable_stack_sha256,
                policy_sha256,
                population_manifest_sha256,
                clinical_manifest_sha256,
                exception_manifest_sha256,
                precision_manifest_sha256,
                scope_policy_sha256, scope_manifest_sha256,
                expected_base_release_sha256 or expected_public_validation_authority_sha256,
            ))
            or type(runtime_image_digest) is not str
            or not runtime_image_digest.startswith("sha256:")
            or not _digest(runtime_image_digest.removeprefix("sha256:"))
        ):
            raise WgsPreparationError("bounded preparation input is invalid")
        validation_authority = None
        qualification = None
        if public_validation_authority_bytes is not None:
            from genosapiens.wgs_validation_authority import load_public_validation_authority
            if (base_release_bytes is not None or expected_base_release_sha256 is not None
                    or public_qualification_authority_bytes is not None or resume):
                raise WgsPreparationError("mixed or resumed public validation authority")
            validation_authority = load_public_validation_authority(public_validation_authority_bytes,
                expected_sha256=expected_public_validation_authority_sha256)
            validation_authority.require_execution(original_input_path=original,
                hpo_input_sha256=hpo_input_sha256, output_path=root,
                executable_stack_sha256=executable_stack_sha256,
                runtime_image_digest=runtime_image_digest,
                preparation_run_options_sha256=options_sha256,
                verify_input_bytes=validation_authority.completed_source_input is None)
        elif expected_public_validation_authority_sha256 is not None:
            raise WgsPreparationError("public validation authority is incomplete")
        if (
            (public_qualification_authority_bytes is None)
            != (expected_public_qualification_authority_sha256 is None)
        ):
            raise WgsPreparationError(
                "public qualification authority binding is incomplete"
            )
        if public_qualification_authority_bytes is not None:
            from genosapiens.wgs_public_qualification import (
                load_public_qualification_authority,
                require_public_qualification_execution,
            )

            qualification = load_public_qualification_authority(
                public_qualification_authority_bytes,
                expected_sha256=expected_public_qualification_authority_sha256,
            )
            require_public_qualification_execution(
                qualification,
                original_input_path=original,
                hpo_input_path=hpo,
                output_path=root,
                input_sha256=qualification.input_sha256,
                input_size_bytes=qualification.input_size_bytes,
                hpo_input_sha256=hpo_input_sha256,
                execution_config=execution_config,
                executable_stack_sha256=executable_stack_sha256,
                runtime_image_digest=runtime_image_digest,
                policy_sha256=policy_sha256,
                population_manifest_sha256=population_manifest_sha256,
                clinical_manifest_sha256=clinical_manifest_sha256,
                exception_manifest_sha256=exception_manifest_sha256,
                precision_manifest_sha256=precision_manifest_sha256,
                preparation_run_options_sha256=options_sha256,
            )
        authority = validation_authority or qualification
        if authority is not None:
            if completed_source_input is not None and completed_source_input != authority.completed_source_input:
                raise WgsPreparationError("completed source input differs from authority")
            completed_source_input = authority.completed_source_input
            source_origin = "completed_input" if completed_source_input is not None else "produced_this_run"
        execution_limits_sha256 = hashlib.sha256(
            bounded_wgs_execution_config_as_json(execution_config)
        ).hexdigest()
        options_path = root / "run-options.json"
        if options_path.exists():
            if not resume or options_path.read_bytes() != run_options.as_bytes():
                raise WgsPreparationError("preparation run options changed")
        else:
            publish_bytes(root, "run-options.json", run_options.as_bytes())
        for child in (
            "prefilter", "working-set", "phenotype", "prepared",
            "journals", "journals/prefilter",
            *(("source", "journals/source", "journals/source_audit")
              if completed_source_input is None else ()),
        ):
            _private_directory(root / child)

        with runtime_phase(phase_recorder, "hpo_initialization"):
            phenotype = ingest_phenotypes(
                hpo,
                expected_input_sha256=hpo_input_sha256,
                hpo_ontology=hpo_ontology,
                ontology_manifest_sha256=ontology_manifest_sha256,
                reference_root=Path(hpo_reference_root),
                expected_input_size_bytes=hpo_input_size_bytes,
            )
        phenotype_content = phenotype_ingestion_as_json(phenotype).encode("ascii")
        phenotype_expected = Artifact(
            "phenotype-ingestion.json",
            hashlib.sha256(phenotype_content).hexdigest(),
            len(phenotype_content),
        )
        phenotype_path = root / "phenotype" / phenotype_expected.name
        if phenotype_path.exists():
            if (
                not resume
                or read_artifact(root / "phenotype", phenotype_expected)
                != phenotype_content
            ):
                raise WgsPreparationError("phenotype resume artifact is invalid")
            phenotype_artifact = phenotype_expected
        else:
            phenotype_artifact = publish_bytes(
                root / "phenotype",
                phenotype_expected.name,
                phenotype_content,
            )

        if completed_source_input is None:
            stage = "source"
            source_parent = _preparation_stage_parent(
                stage,
                case_manifest_sha256=case_manifest_sha256,
                selected_sample_sha256=selected_sample_sha256,
                preceding_sha256=hashlib.sha256(
                    concordance_report_as_json(concordance_report).encode("ascii")
                ).hexdigest(),
                execution_limits_sha256=execution_limits_sha256,
                executable_stack_sha256=executable_stack_sha256,
                runtime_image_digest=runtime_image_digest,
                run_options_sha256=options_sha256,
            )

            def source_action():
                source_plan_sha256 = hashlib.sha256(canonical({
                    "schema": PREPARATION_PLAN_SCHEMA,
                    "version": 1,
                    "stage": "source",
                    "partitions": ("complete-source",),
                    "complete_source_required": True,
                    "classification": "controlled",
                    "llm_policy": "no-llm",
                })).hexdigest()
                source_binding = StageRunBinding(
                    "source", source_parent, executable_stack_sha256,
                    runtime_image_digest, execution_limits_sha256,
                    options_sha256,
                    concordance_report.input_sha256,
                    concordance_report.reference_manifest_sha256,
                    policy_sha256,
                    (population_manifest_sha256, clinical_manifest_sha256,
                     exception_manifest_sha256, precision_manifest_sha256,
                     ontology_manifest_sha256),
                    (
                        ParentIdentity("stage_parent", PREPARATION_STAGE_PARENT_SCHEMA,
                                       1, source_parent),
                        ParentIdentity("partition_plan", PREPARATION_PLAN_SCHEMA,
                                       1, source_plan_sha256),
                        ParentIdentity("execution_limits", execution_config.schema,
                                       execution_config.version, execution_limits_sha256),
                        ParentIdentity("executable_stack", "genosapiens.executable-stack",
                                       1, executable_stack_sha256),
                        ParentIdentity("runtime_image", "genosapiens.runtime-image",
                                       1, runtime_image_digest.removeprefix("sha256:")),
                    ),
                )
                with runtime_phase(phase_recorder, "source_ingestion"):
                    source_result = build_wgs_source_ledger(
                        original,
                        concordance_report=concordance_report,
                        fetch_reference=fetch_reference,
                        output_root=root / "source",
                        selected_sample_name=selected_sample_name,
                        selected_sample_sha256=selected_sample_sha256,
                        case_manifest_sha256=case_manifest_sha256,
                        executable_stack_sha256=executable_stack_sha256,
                        execution_config=execution_config,
                        allele_balance_screen=allele_balance_screen,
                        source_state_path=root / "state/source.sqlite",
                        stage_binding=source_binding,
                        resume=resume,
                        worker_reference=source_worker_reference,
                        operator_case_binding=operator_case_binding,
                    )
                with runtime_phase(phase_recorder, "source_accounting"):
                    if source_result.qc_summary is None:
                        manifest = source_result.manifest
                        reports = manifest.summary.reports
                        if tuple(item.name for item in reports) != (
                                "qc.json", "normalization.json", "classes.json"):
                            raise WgsPreparationError("completed source reports are inconsistent")
                        accounting = SourceAccountingClosure(
                            manifest.sha256, manifest.source.input_sha256,
                            manifest.header.sha256, *(item.sha256 for item in reports),
                            manifest.summary.total_records,
                            manifest.summary.total_physical_alternate_slots,
                            manifest.summary.total_biological_alternates,
                        )
                    else:
                        accounting = source_accounting_from_build(source_result)
                    _publish_or_reuse_small(root / "source", "accounting.json", accounting.as_bytes())
                return (source_result, accounting), source_result.manifest.sha256

            def source_resume_action(output_sha256):
                try:
                    content = _read_bounded_manifest(root / "source/manifest.json")
                    manifest = load_model(content, LedgerManifest)
                except (WgsPreparationError, ValueError, TypeError):
                    raise WgsPreparationError("completed source manifest is invalid") from None
                if manifest.sha256 != output_sha256:
                    raise WgsPreparationError("completed source manifest changed")
                artifact = Artifact(
                    "manifest.json", hashlib.sha256(content).hexdigest(), len(content),
                )
                if read_artifact(root / "source", artifact) != content:
                    raise WgsPreparationError("completed source manifest changed")
                accounting_path = root / "source/accounting.json"
                try:
                    accounting = load_source_accounting(
                        _read_bounded_manifest(accounting_path))
                except (WgsPreparationError, OSError):
                    raise WgsPreparationError("completed source accounting is invalid") from None
                _require_source_accounting(manifest, accounting)
                return (
                    WgsSourceBuildResult(manifest, artifact, None, None, None), accounting,
                ), manifest.sha256

            source_pair, journal = _execute_preparation_stage(
                journal_root=root / "journals/source",
                stage="source",
                partition_id="complete-source",
                stage_parent_sha256=source_parent,
                execution_config=execution_config,
                execution_limits_sha256=execution_limits_sha256,
                executable_stack_sha256=executable_stack_sha256,
                runtime_image_digest=runtime_image_digest,
                scientific=False,
                action=source_action,
                resume_action=source_resume_action,
                resume=resume,
            )
            source_result, source_accounting = source_pair
            accounting_sha256 = source_accounting.sha256
            journals.append(journal)

            stage = "source_audit"
            audit_parent = _preparation_stage_parent(
                stage, case_manifest_sha256=case_manifest_sha256,
                selected_sample_sha256=selected_sample_sha256,
                preceding_sha256=source_accounting.sha256,
                execution_limits_sha256=execution_limits_sha256,
                executable_stack_sha256=executable_stack_sha256,
                runtime_image_digest=runtime_image_digest,
                run_options_sha256=options_sha256,
            )
            audit_plan_sha256 = hashlib.sha256(canonical({
                "schema": PREPARATION_PLAN_SCHEMA, "version": 1,
                "stage": "source_audit",
                "partitions": tuple(
                    f"audit-shard-{shard.ordinal:06d}"
                    for shard in source_result.manifest.shards
                ),
                "source_manifest_sha256": source_result.manifest.sha256,
                "classification": "controlled", "llm_policy": "no-llm",
            })).hexdigest()
            audit_binding = StageRunBinding(
                "source_audit", audit_parent, executable_stack_sha256,
                runtime_image_digest, execution_limits_sha256, options_sha256,
                concordance_report.input_sha256,
                concordance_report.reference_manifest_sha256,
                policy_sha256,
                (source_result.manifest.sha256, population_manifest_sha256,
                 clinical_manifest_sha256, exception_manifest_sha256,
                 precision_manifest_sha256, ontology_manifest_sha256),
                (
                    ParentIdentity("stage_parent", PREPARATION_STAGE_PARENT_SCHEMA,
                                   1, audit_parent),
                    ParentIdentity("partition_plan", PREPARATION_PLAN_SCHEMA,
                                   1, audit_plan_sha256),
                    ParentIdentity("execution_limits", execution_config.schema,
                                   execution_config.version, execution_limits_sha256),
                    ParentIdentity("executable_stack", "genosapiens.executable-stack",
                                   1, executable_stack_sha256),
                    ParentIdentity("runtime_image", "genosapiens.runtime-image",
                                   1, runtime_image_digest.removeprefix("sha256:")),
                ),
            )

            def audit_result(status, failure_category=None):
                return SourceAuditResult(
                    run_options.verify_source_ledger, status, source_accounting.sha256,
                    source_result.manifest.sha256, options_sha256,
                    executable_stack_sha256, runtime_image_digest,
                    concordance_report.reference_manifest_sha256,
                    concordance_report.input_sha256, failure_category,
                )

            def audit_action():
                nonlocal audit_result_sha256
                if run_options.verify_source_ledger:
                    try:
                        with runtime_phase(phase_recorder, "source_audit"):
                            replay = verify_source_ledger(
                                original, concordance_report=concordance_report,
                                fetch_reference=fetch_reference, ledger_root=root / "source",
                                manifest=source_result.manifest,
                                selected_sample_name=selected_sample_name,
                                selected_sample_sha256=selected_sample_sha256,
                                case_manifest_sha256=case_manifest_sha256,
                                executable_stack_sha256=executable_stack_sha256,
                                execution_config=execution_config,
                                allele_balance_screen=allele_balance_screen,
                                audit_state_path=root / "state/source-audit.sqlite",
                                stage_binding=audit_binding,
                                resume=resume,
                            )
                        if replay != source_accounting:
                            raise WgsPreparationError("source audit accounting mismatch")
                    except KeyboardInterrupt:
                        raise
                    except Exception as error:
                        failure = static_failure(
                            error, "source_audit", operation="preparation_stage_action",
                        )
                        if failure.category == "interrupted":
                            raise
                        category = (
                            "invalid_source" if isinstance(error, VcfQcError)
                            else "semantic_mismatch"
                            if isinstance(error, (WgsPreparationError, ValueError))
                            else "unexpected"
                        )
                        failed = audit_result("failed", category)
                        _publish_or_reuse_small(root / "source", "audit.json", failed.as_bytes())
                        audit_result_sha256 = failed.sha256
                        raise WgsPreparationError("source audit failed") from None
                    outcome = audit_result("passed")
                else:
                    outcome = audit_result("not_requested")
                _publish_or_reuse_small(root / "source", "audit.json", outcome.as_bytes())
                audit_result_sha256 = outcome.sha256
                return outcome, outcome.sha256

            def audit_resume_action(output_sha256):
                try:
                    outcome = load_source_audit(_read_bounded_manifest(root / "source/audit.json"))
                except (WgsPreparationError, OSError):
                    raise WgsPreparationError("completed source audit is invalid") from None
                if (outcome.sha256 != output_sha256
                        or outcome != audit_result(
                            "passed" if run_options.verify_source_ledger else "not_requested")):
                    raise WgsPreparationError("completed source audit changed")
                return outcome, outcome.sha256

            source_audit, audit_journal = _execute_preparation_stage(
                journal_root=root / "journals/source_audit",
                stage="source_audit", partition_id="complete-source-audit",
                stage_parent_sha256=audit_parent,
                execution_config=execution_config,
                execution_limits_sha256=execution_limits_sha256,
                executable_stack_sha256=executable_stack_sha256,
                runtime_image_digest=runtime_image_digest,
                scientific=False,
                action=audit_action,
                resume_action=audit_resume_action,
                resume=resume,
            )
            audit_result_sha256 = source_audit.sha256
            journals.append(audit_journal)

            descriptor_path = root / "source-input/completed-source.json"
            if resume and (descriptor_path.exists() or descriptor_path.is_symlink()):
                # This is the producer's immutable identity, already bound into
                # prefilter checkpoints. A new invocation's runtime path must
                # never replace the original producer's runtime checkpoint.
                descriptor = load_completed_source_input(_read_bounded_manifest(descriptor_path))
                if (descriptor.source_root != str(root / "source")
                        or descriptor.producer.source != source_result.manifest.source
                        or descriptor.source_manifest.sha256 != source_result.manifest.sha256
                        or descriptor.source_accounting.sha256 != source_accounting.sha256
                        or descriptor.source_audit.sha256 != source_audit.sha256):
                    raise WgsPreparationError("completed source resume binding changed")
            else:
                descriptor = _describe_completed_source(root / "source",
                    candidate_sha256=qualification.engineering_candidate_sha256 if qualification is not None else None,
                    journals=tuple(journals),
                    runtime_checkpoint_path=None if phase_recorder is None else phase_recorder.path)
            _private_directory(root / "source-input")
            artifact = publish_completed_source_input(root / "source-input", descriptor)
            completed_source_input = CompletedSourceArtifact(
                str(root / "source-input" / artifact.name), COMPLETED_SOURCE_SCHEMA,
                1, artifact.sha256, artifact.size_bytes)
        else:
            if (type(completed_source_input) is not CompletedSourceArtifact
                    or completed_source_input.schema != COMPLETED_SOURCE_SCHEMA):
                raise WgsPreparationError("completed source input binding is invalid")
            descriptor = load_completed_source_input(completed_source_input.read())
        if authority is not None:
            expected_input_size = (qualification.input_size_bytes if qualification is not None else
                next(row.size_bytes for row in authority.inputs if row.path == str(original)))
        else:
            expected_input_size = descriptor.producer.source.input_size_bytes
        expected_source = SourceBinding(
            concordance_report.input_sha256, expected_input_size,
            case_manifest_sha256, selected_sample_sha256,
            hashlib.sha256(concordance_report_as_json(concordance_report).encode("ascii")).hexdigest(),
            concordance_report.reference_manifest_sha256, executable_stack_sha256)
        completed_source = resolve_completed_source(descriptor,
            expected_source=expected_source, audit_requested=run_options.verify_source_ledger)
        source_input_sha256 = descriptor.sha256
        _publish_or_reuse_small(root, "source-input-binding.json", canonical(asdict(completed_source_input)))
        if source_origin == "completed_input":
            source_result = WgsSourceBuildResult(completed_source.manifest,
                Artifact("manifest.json", descriptor.source_manifest.sha256,
                         descriptor.source_manifest.size_bytes), None, None, None)
            source_accounting = completed_source.accounting
            source_audit = completed_source.audit
            journals.extend(completed_source.journals)
            accounting_sha256 = source_accounting.sha256
            audit_result_sha256 = source_audit.sha256
        source_root = completed_source.root

        stage = "prefilter"
        filter_parent = _preparation_stage_parent(
            stage,
            case_manifest_sha256=case_manifest_sha256,
            selected_sample_sha256=selected_sample_sha256,
            preceding_sha256=source_result.manifest.sha256,
            execution_limits_sha256=execution_limits_sha256,
            executable_stack_sha256=executable_stack_sha256,
            runtime_image_digest=runtime_image_digest,
            run_options_sha256=options_sha256,
            completed_source_input_sha256=source_input_sha256,
        )
        context = PreparationFilterContext(
            source_root,
            source_result.manifest,
            source_accounting,
            policy_sha256,
            population_manifest_sha256,
            clinical_manifest_sha256,
            exception_manifest_sha256,
            precision_manifest_sha256,
            scope_policy_sha256, scope_manifest_sha256,
            completed_source=completed_source,
        )
        # Only the authenticated operator continuation of the diagnosed source
        # classification failure may exceed the normal prefilter attempt cap.
        # Retain the original events and all policy/checkpoint parent identities.
        from genosapiens.wgs_operator_invocation import OperatorCaseBinding
        technical_classification_resume = bool(
            resume and type(operator_case_binding) is OperatorCaseBinding
            and operator_case_binding.classification_resume_authorized
            and (root / "journals/prefilter").exists()
            and any(entry.event_type == "attempt_failed"
                    and entry.failure_code == "source_classification_malformed_evidence"
                    for entry in read_journal(root / "journals/prefilter",
                        stage="prefilter", parent_sha256=filter_parent))
        )
        source_routing_context = source_routing_context_provider(context)
        if not isinstance(source_routing_context, SourceRoutingContext):
            raise WgsPreparationError("source routing context provider is invalid")
        context = replace(
            context,
            source_cursor=SourceDecisionCursor(
                source_root, source_result.manifest, source_routing_context,
                max_pending=(
                    2 * execution_config.public_queries.concurrent_batches
                    * execution_config.public_queries.batch_rows
                ),
            ),
        )

        filter_plan_sha256 = hashlib.sha256(canonical({
            "schema": PREPARATION_PLAN_SCHEMA, "version": 1,
            "stage": "prefilter", "partitions": ("complete-source-filter",),
            "complete_source_required": True,
            "classification": "controlled", "llm_policy": "no-llm",
        })).hexdigest()
        filter_binding = StageRunBinding(
            "prefilter", filter_parent, executable_stack_sha256,
            runtime_image_digest, execution_limits_sha256, options_sha256,
            source_result.manifest.source.input_sha256,
            source_result.manifest.source.reference_manifest_sha256,
            policy_sha256,
            (population_manifest_sha256, clinical_manifest_sha256,
             exception_manifest_sha256, precision_manifest_sha256,
             scope_policy_sha256, scope_manifest_sha256),
            (
                ParentIdentity("stage_parent", PREPARATION_STAGE_PARENT_SCHEMA,
                               1, filter_parent),
                ParentIdentity("partition_plan", PREPARATION_PLAN_SCHEMA,
                               1, filter_plan_sha256),
                ParentIdentity("execution_limits", execution_config.schema,
                               execution_config.version, execution_limits_sha256),
                ParentIdentity("executable_stack", "genosapiens.executable-stack",
                               1, executable_stack_sha256),
                ParentIdentity("runtime_image", "genosapiens.runtime-image",
                               1, runtime_image_digest.removeprefix("sha256:")),
                ParentIdentity("completed_source_input", COMPLETED_SOURCE_SCHEMA,
                               1, source_input_sha256),
            ),
        )

        def filter_action():
            _private_directory(root / "state")
            decisions = decision_provider(context)
            with runtime_phase(phase_recorder, "filter_closure"):
                with runtime_phase(phase_recorder, "filter_publication"):
                    filtered = publish_filter_decisions(
                        source_root=source_root,
                        source_manifest=source_result.manifest,
                        source_accounting=source_accounting,
                        source_routing_context=source_routing_context,
                        source_cursor=context.source_cursor,
                        decisions=decisions,
                        output_root=root / "prefilter",
                        policy_sha256=policy_sha256,
                        population_manifest_sha256=population_manifest_sha256,
                        clinical_manifest_sha256=clinical_manifest_sha256,
                        exception_manifest_sha256=exception_manifest_sha256,
                        precision_manifest_sha256=precision_manifest_sha256,
                        scope_policy_sha256=scope_policy_sha256,
                        scope_manifest_sha256=scope_manifest_sha256,
                        retained_state_path=root / "state/filter.sqlite",
                        stage_binding=filter_binding,
                        resume=resume,
                        technical_classification_resume=technical_classification_resume,
                        max_attempts=(
                            execution_config.journal.default_scientific_attempts_per_partition
                        ),
                        phase_recorder=phase_recorder,
                    )
            with runtime_phase(phase_recorder, "working_set_closure"):
                with runtime_phase(phase_recorder, "working_set_publication"):
                    working = publish_working_set(
                        source_root=source_root,
                        source_manifest=source_result.manifest,
                        source_accounting=source_accounting,
                        source_routing_context=source_routing_context,
                        filter_root=root / "prefilter",
                        filter_manifest=filtered.manifest,
                        output_root=root / "working-set",
                        policy_sha256=policy_sha256,
                        population_manifest_sha256=population_manifest_sha256,
                        clinical_manifest_sha256=clinical_manifest_sha256,
                        exception_manifest_sha256=exception_manifest_sha256,
                        precision_manifest_sha256=precision_manifest_sha256,
                        scope_policy_sha256=scope_policy_sha256,
                        scope_manifest_sha256=scope_manifest_sha256,
                        retained_state_path=root / "state/filter.sqlite",
                    )
            return (filtered, working), working.manifest.sha256

        def filter_resume_action(output_sha256):
            try:
                filter_content = _read_bounded_manifest(root / "prefilter/manifest.json")
                working_content = _read_bounded_manifest(root / "working-set/manifest.json")
                filter_manifest = load_model(filter_content, FilterManifest)
                working_manifest = load_model(working_content, WorkingSetManifest)
            except (WgsPreparationError, ValueError, TypeError):
                raise WgsPreparationError("completed prefilter manifests are invalid") from None
            if working_manifest.sha256 != output_sha256:
                raise WgsPreparationError("completed prefilter output changed")
            resources = _check_filter_read_boundary(
                source_manifest=source_result.manifest,
                source_accounting=source_accounting,
                source_routing_context=source_routing_context,
                filter_root=root / "prefilter",
                filter_manifest=filter_manifest,
                policy_sha256=policy_sha256,
                population_manifest_sha256=population_manifest_sha256,
                clinical_manifest_sha256=clinical_manifest_sha256,
                exception_manifest_sha256=exception_manifest_sha256,
                precision_manifest_sha256=precision_manifest_sha256,
                scope_policy_sha256=scope_policy_sha256,
                scope_manifest_sha256=scope_manifest_sha256,
            )
            if (working_manifest.parent_filter_manifest_sha256 != filter_manifest.sha256
                    or working_manifest.parents
                    != _working_parents(source_result.manifest, filter_manifest)
                    or working_manifest.summary.total_source_occurrences
                    != source_accounting.total_physical_alternate_slots):
                raise WgsPreparationError("completed working manifest is inconsistent")
            with _FilterState.open_sealed(
                root / "state/filter.sqlite",
                expected=source_accounting.total_physical_alternate_slots,
                source_manifest_sha256=source_result.manifest.sha256,
                policy_sha256=filter_manifest.policy_sha256,
                resource_sha256=resources,
                scope_sha256=(scope_policy_sha256, scope_manifest_sha256),
                reference_sha256=source_result.manifest.source.reference_manifest_sha256,
                filter_manifest_sha256=filter_manifest.sha256,
            ) as retained_filter:
                if retained_filter.working_binding() != (
                    working_manifest.sha256, working_manifest.summary.sha256,
                ):
                    raise WgsPreparationError("completed working state changed")
            filter_summary = filter_manifest.summary.as_bytes()
            working_summary = working_manifest.summary.as_bytes()
            filtered = FilterPublication(
                filter_manifest,
                Artifact(
                    "manifest.json", hashlib.sha256(filter_content).hexdigest(),
                    len(filter_content),
                ),
                Artifact(
                    "summary.json", hashlib.sha256(filter_summary).hexdigest(),
                    len(filter_summary),
                ),
            )
            working = WorkingSetPublication(
                working_manifest,
                Artifact(
                    "manifest.json", hashlib.sha256(working_content).hexdigest(),
                    len(working_content),
                ),
                Artifact(
                    "summary.json", hashlib.sha256(working_summary).hexdigest(),
                    len(working_summary),
                ),
            )
            return (filtered, working), working_manifest.sha256

        filter_pair, journal = _execute_preparation_stage(
            journal_root=root / "journals/prefilter",
            stage="prefilter",
            partition_id="complete-source-filter",
            stage_parent_sha256=filter_parent,
            execution_config=execution_config,
            execution_limits_sha256=execution_limits_sha256,
            executable_stack_sha256=executable_stack_sha256,
            runtime_image_digest=runtime_image_digest,
            scientific=True,
            technical_classification_resume=technical_classification_resume,
            input_parents=(ParentIdentity("completed_source_input", COMPLETED_SOURCE_SCHEMA,
                                          1, source_input_sha256),),
            action=filter_action,
            resume_action=filter_resume_action,
            resume=resume,
            prepare_retry=lambda attempt: _preserve_retry_directories(
                root,
                stage="prefilter",
                attempt_number=attempt,
                directories=(
                    (root / "working-set",)
                    if (root / "state/filter.sqlite").exists()
                    else (root / "prefilter", root / "working-set")
                ),
            ),
        )
        filtered, working = filter_pair
        journals.append(journal)

        stage = "prepared_case"
        prepared_index = build_prepared_index(
            state_root=root / "state",
            completed_source=completed_source,
            executable_stack_sha256=executable_stack_sha256,
            execution_config_sha256=execution_limits_sha256,
            runtime_image_digest=runtime_image_digest,
            source_manifest=source_result.manifest,
            filter_manifest=filtered.manifest,
            working_manifest=working.manifest,
            source_accounting_sha256=accounting_sha256,
            source_audit_sha256=audit_result_sha256,
            preparation_run_options_sha256=options_sha256,
            journals=tuple(sorted(journals, key=lambda value: value.stage)),
            audit_requested=run_options.verify_source_ledger,
        )
        prepared_index_sha256 = publish_prepared_index(
            root / "state", prepared_index,
        ).sha256
        from genosapiens.track1_installation import load_run_versions as load_base_release
        release = validation_authority or load_base_release(
            base_release_bytes,
            expected_sha256=expected_base_release_sha256,
        )
        bindings = release.policy_bindings
        diagnostic = PreparationDiagnosticReceipt(
            case_manifest_sha256,
            selected_sample_sha256,
            source_result.manifest.source.concordance_sha256,
            source_result.manifest.sha256,
            _source_accounting_sha256(source_accounting),
            filtered.manifest.sha256,
            working.manifest.sha256,
            hpo_input_sha256,
            phenotype_artifact.sha256,
            ontology_manifest_sha256,
            bindings.diagnostic_policy_sha256,
            "complete",
            "complete",
            "passed",
            "complete_original_source",
            "complete",
            "validated",
            (),
            audit_result_sha256,
            options_sha256,
            prepared_index_sha256,
            run_options.verify_source_ledger,
            source_audit.status,
            source_input_sha256,
        )
        with runtime_phase(phase_recorder, "preparation_terminal_publication"):
            prepared = publish_prepared_case(
                output_root=root / "prepared",
                completed_source_input=completed_source_input,
                source_origin=source_origin,
                original_input_path=original,
                source_root=source_root,
                source_manifest=source_result.manifest,
                source_accounting=source_accounting,
                filter_root=root / "prefilter",
                filter_manifest=filtered.manifest,
                working_root=root / "working-set",
                working_manifest=working.manifest,
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
                phenotype_ingestion_sha256=phenotype_artifact.sha256,
                ontology_manifest_sha256=ontology_manifest_sha256,
                diagnostic=diagnostic,
                execution_config=execution_config,
                base_release_bytes=base_release_bytes,
                expected_base_release_sha256=expected_base_release_sha256,
                executable_stack_sha256=executable_stack_sha256,
                runtime_image_digest=runtime_image_digest,
                journals=tuple(sorted(journals, key=lambda value: value.stage)),
                hpo_input_path=hpo,
                public_qualification_authority_bytes=(
                    public_qualification_authority_bytes
                ),
                expected_public_qualification_authority_sha256=(
                    expected_public_qualification_authority_sha256
                ),
                public_validation_authority_bytes=public_validation_authority_bytes,
                expected_public_validation_authority_sha256=expected_public_validation_authority_sha256,
            )
            return _publish_preparation_receipt(
                root,
                status="passed",
                terminal_stage=stage,
                base_release_sha256=expected_base_release_sha256,
                journals=tuple(journals),
                preparation_run_options_sha256=options_sha256,
                audit_requested=run_options.verify_source_ledger,
                source_accounting_sha256=accounting_sha256,
                source_audit_sha256=audit_result_sha256,
                completed_source_input_sha256=source_input_sha256,
                source_origin=source_origin,
                prepared_index_sha256=prepared_index_sha256,
                prepared_case=prepared,
                public_validation_authority_sha256=expected_public_validation_authority_sha256,
            )
    except _PreparationInterrupted as error:
        with runtime_phase(phase_recorder, "preparation_terminal_publication"):
            return _publish_preparation_receipt(
                root,
                status="interrupted",
                terminal_stage=stage,
                base_release_sha256=expected_base_release_sha256,
                journals=tuple(journals),
                preparation_run_options_sha256=options_sha256,
                audit_requested=run_options.verify_source_ledger,
                source_accounting_sha256=accounting_sha256,
                source_audit_sha256=audit_result_sha256,
                completed_source_input_sha256=source_input_sha256,
                source_origin=source_origin,
                prepared_index_sha256=prepared_index_sha256,
                failure=error.failure,
                public_validation_authority_sha256=expected_public_validation_authority_sha256,
            )
    except Exception as error:
        with runtime_phase(phase_recorder, "preparation_terminal_publication"):
            return _publish_preparation_receipt(
                root,
                status="failed",
                terminal_stage=stage,
                base_release_sha256=expected_base_release_sha256,
                journals=tuple(journals),
                preparation_run_options_sha256=options_sha256,
                audit_requested=run_options.verify_source_ledger,
                source_accounting_sha256=accounting_sha256,
                source_audit_sha256=audit_result_sha256,
                completed_source_input_sha256=source_input_sha256,
                source_origin=source_origin,
                prepared_index_sha256=prepared_index_sha256,
                failure=static_failure(
                    error, stage, operation={
                        "input_authentication": "runner_input_authentication",
                        "source": "source_stage",
                        "source_audit": "preparation_stage_action",
                        "prefilter": "preparation_stage_action",
                        "prepared_case": "prepared_case_publish",
                    }[stage],
                ),
                public_validation_authority_sha256=expected_public_validation_authority_sha256,
            )
