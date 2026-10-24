"""Bounded, allowlisted WGS failures; never persist exception detail."""

from __future__ import annotations

from contextlib import contextmanager
from dataclasses import dataclass
import hashlib
import json
import os
from pathlib import Path
import sqlite3
import stat


ORIGINS = frozenset({
    "population_lookup", "population_adaptation", "clinical_index",
    "clinical_replay", "clinical_adaptation", "decision_conversion",
    "resource_authentication", "source_decode", "source_classification",
    "population_request", "clinical_request", "source_accounting",
    "filter_publication", "working_set_publication", "resource_session",
    "source", "source_audit", "prefilter", "input_authentication", "prepared_case",
    "prepared_authority", "bounded_products", "final_publication",
    "annotation_child", "unexpected",
})
CATEGORIES = frozenset({
    "authentication", "malformed_evidence", "missing_response", "response_capacity",
    "timeout", "interrupted", "memory_exhaustion", "storage", "unexpected",
})
FAILURE_CODES = frozenset(
    origin + "_" + category for origin in ORIGINS for category in CATEGORIES
)

# Static code locations only. Values derived from a record, path, sample, or
# exception are never admitted into this vocabulary.
SOURCE_OPERATIONS = frozenset({
    "source_checkpoint_open", "source_input_open", "source_reducer_recovery",
    "source_qc_initialize", "source_class_initialize", "source_normalization_initialize",
    "source_worker_execute", "source_row_reduce", "source_shard_publish", "source_shard_commit",
})
OPERATIONS = frozenset({
    "outermost_unknown",
    "runner_authority_authentication", "runner_input_authentication",
    "runner_resource_authentication", "runner_index_reconciliation",
    "preparation_stage_action", "analysis_stage_action",
    "source_stage", "source_payload_decode", "source_disposition_classify",
    "population_request_build", "clinical_request_build",
    "population_query", "population_response_bind", "population_response_adapt",
    "clinical_index_query", "clinical_response_bind", "clinical_response_adapt",
    "decision_convert", "source_replay_decode", "filter_source_tuple_compare",
    "filter_row_validate", "filter_artifact_publish",
    "working_set_validate", "working_set_publish",
    "resource_session_construct", "resource_session_open",
    "resource_session_work", "resource_session_close",
    "analysis_product_build", "annotation_replay",
    "clinical_association_batch", "candidate_admission_convert",
    "context_assemble", "context_annotate", "context_annotation_inputs",
    "context_population_query", "context_population_bind",
    "context_clinical_query", "context_clinical_bind",
    "context_annotation_finalize", "hypothesis_generate",
    "population_assess", "evidence_adapt", "scientific_manifest_write",
    "prepared_case_validate", "prepared_case_publish",
    "prepared_authority_validate", "bounded_product_validate", "pair_decision_header",
    "annotation_query", "annotation_child_execute", "final_publish",
    "journal_publish", "receipt_publish", "diagnostic_publish",
    "runtime_finalize",
    "context_job_encode", "context_worker_execute", "context_result_encode",
    "context_result_validate", "context_store", "context_merge",
    "context_assessment_publish", "graph_job_encode", "graph_worker_execute",
    "graph_result_validate", "worker_shutdown",
    "context_graph_evaluate", "context_evidence_attach", "context_payload_encode",
    "context_fragment_encode", "context_record_encode",
}) | SOURCE_OPERATIONS
EXCEPTION_KINDS = frozenset({
    "project_failure", "memory_error", "os_error", "value_error", "type_error",
    "key_error", "sqlite_error", "interrupted", "other",
})
MAX_SECONDARY_FAILURES = 4
MAX_FAILURE_EVIDENCE_BYTES = 16 * 1024
MAX_FAILURE_JOURNAL_DIRECTORIES = 64
_EVIDENCE_FIELDS = frozenset({
    "schema", "version", "classification", "llm_policy", "stage",
    "receipt_sha256", "failure_code", "operation", "exception_kind",
    "journal_heads", "secondary_failures", "secondary_failures_truncated",
})
_SECONDARY_FIELDS = frozenset({"failure_code", "operation", "exception_kind"})


def _digest(value: object) -> bool:
    return (
        type(value) is str
        and len(value) == 64
        and all(character in "0123456789abcdef" for character in value)
    )


def _read_bounded_regular(path: Path, maximum_bytes: int) -> bytes:
    descriptor = -1
    try:
        descriptor = os.open(path, os.O_RDONLY | os.O_NOFOLLOW)
        before = os.fstat(descriptor)
        if (
            not stat.S_ISREG(before.st_mode)
            or before.st_size < 1
            or before.st_size > maximum_bytes
        ):
            raise OSError
        with os.fdopen(descriptor, "rb", closefd=True) as stream:
            descriptor = -1
            content = stream.read(maximum_bytes + 1)
            after = os.fstat(stream.fileno())
        if (
            len(content) != before.st_size
            or before.st_dev != after.st_dev
            or before.st_ino != after.st_ino
            or before.st_size != after.st_size
            or before.st_mtime_ns != after.st_mtime_ns
            or before.st_ctime_ns != after.st_ctime_ns
        ):
            raise OSError
        return content
    except OSError:
        raise ValueError("bounded failure input is invalid") from None
    finally:
        if descriptor >= 0:
            os.close(descriptor)


def exception_kind(error: BaseException) -> str:
    """Map explicit built-in exception families before discarding detail."""

    if isinstance(error, WgsFailure):
        return error.exception_kind
    if isinstance(error, MemoryError):
        return "memory_error"
    if isinstance(error, (KeyboardInterrupt, SystemExit, GeneratorExit)):
        return "interrupted"
    if isinstance(error, sqlite3.Error):
        return "sqlite_error"
    if isinstance(error, OSError):
        return "os_error"
    if isinstance(error, ValueError):
        return "value_error"
    if isinstance(error, TypeError):
        return "type_error"
    if isinstance(error, KeyError):
        return "key_error"
    return "other"


@dataclass(frozen=True, slots=True)
class SecondaryFailure:
    failure_code: str
    operation: str
    exception_kind: str

    def __post_init__(self) -> None:
        if (
            self.failure_code not in FAILURE_CODES
            or self.operation not in OPERATIONS
            or self.exception_kind not in EXCEPTION_KINDS
        ):
            raise ValueError("invalid static secondary failure")

    def as_dict(self) -> dict[str, str]:
        return {
            "failure_code": self.failure_code,
            "operation": self.operation,
            "exception_kind": self.exception_kind,
        }


class WgsFailure(RuntimeError):
    def __init__(
        self,
        origin=None,
        category=None,
        *,
        operation=None,
        exception_kind=None,
        journal_heads=(),
        secondary_failures=(),
        secondary_failures_truncated=False,
    ):
        # The no-argument form is reserved for a true outermost boundary whose
        # operation could not be classified. All classified constructors must
        # supply every static field.
        if origin is category is operation is exception_kind is None:
            origin = category = "unexpected"
            operation = "outermost_unknown"
            exception_kind = "other"
        if (
            origin not in ORIGINS
            or category not in CATEGORIES
            or operation not in OPERATIONS
            or exception_kind not in EXCEPTION_KINDS
            or type(journal_heads) not in {tuple, list}
            or len(journal_heads) > 16
            or len(set(journal_heads)) != len(journal_heads)
            or any(not _digest(head) for head in journal_heads)
            or type(secondary_failures) not in {tuple, list}
            or len(secondary_failures) > MAX_SECONDARY_FAILURES
            or any(type(value) is not SecondaryFailure for value in secondary_failures)
            or type(secondary_failures_truncated) is not bool
            or secondary_failures_truncated and len(secondary_failures) != MAX_SECONDARY_FAILURES
        ):
            raise ValueError("invalid static failure")
        self.origin = origin
        self.category = category
        self.failure_code = origin + "_" + category
        self.operation = operation
        self.exception_kind = exception_kind
        self.journal_heads = tuple(journal_heads)
        self.secondary_failures = tuple(secondary_failures)
        self.secondary_failures_truncated = secondary_failures_truncated
        super().__init__(self.failure_code)

    def add_secondary(self, failure: "WgsFailure") -> None:
        """Retain a bounded cleanup/reporting failure without replacing primary."""

        if not isinstance(failure, WgsFailure):
            raise TypeError("secondary failure must be static")
        if len(self.secondary_failures) < MAX_SECONDARY_FAILURES:
            self.secondary_failures += (SecondaryFailure(
                failure.failure_code, failure.operation, failure.exception_kind,
            ),)
        else:
            self.secondary_failures_truncated = True


def _cleanup_details(error):
    values = getattr(error, '_wgs_secondary_failures', ())
    truncated = getattr(error, '_wgs_secondary_failures_truncated', False)
    if (type(values) is not tuple or len(values) > MAX_SECONDARY_FAILURES
            or any(type(v) is not SecondaryFailure for v in values)
            or type(truncated) is not bool
            or truncated and len(values) != MAX_SECONDARY_FAILURES):
        return (), False
    return values, truncated


def static_failure(error, origin=None, *, operation=None):
    if isinstance(error, WgsFailure):
        return error
    if origin is None and operation is None:
        origin = "unexpected"
        operation = "outermost_unknown"
    if origin not in ORIGINS or operation not in OPERATIONS:
        raise ValueError("classified static failure requires allowlisted fields")
    kind = exception_kind(error)
    category = (
        "memory_exhaustion" if kind == "memory_error" else
        "interrupted" if kind == "interrupted" else
        "timeout" if isinstance(error, TimeoutError) else
        "storage" if kind == "os_error" else
        "storage" if isinstance(error, sqlite3.Error) and getattr(error, "sqlite_errorcode", -1) & 255
        in {sqlite3.SQLITE_FULL, sqlite3.SQLITE_IOERR, sqlite3.SQLITE_CANTOPEN, sqlite3.SQLITE_READONLY} else
        "unexpected"
    )
    secondary, truncated = _cleanup_details(error)
    return WgsFailure(origin, category, operation=operation, exception_kind=kind,
        secondary_failures=secondary, secondary_failures_truncated=truncated)


def worker_failure_message(error, *, worker, job, digest, operation):
    """Private job-bound transport contains no exception text or class name."""
    failure = static_failure(error, "bounded_products", operation=operation)
    return ("error", worker, job, digest, failure.category,
            failure.operation, failure.exception_kind)


def checked_worker_failure(message, *, worker, job, digest, operations):
    if (type(message) is not tuple or len(message) != 7
            or message[:4] != ("error", worker, job, digest)
            or type(message[1]) is not int or type(message[2]) is not int
            or any(type(v) is not str for v in message[3:])
            or message[4] not in CATEGORIES or message[5] not in operations
            or message[6] not in EXCEPTION_KINDS):
        raise ValueError("worker failure envelope is invalid")
    return WgsFailure("bounded_products", message[4], operation=message[5],
                      exception_kind=message[6])


def run_cleanup(actions, *, primary=None):
    """Attempt every bounded cleanup and preserve the first failure."""
    failure = primary
    for operation, action in actions:
        try:
            action()
        except BaseException as error:
            secondary = static_failure(error, "bounded_products", operation=operation)
            if failure is None:
                failure = secondary
            elif isinstance(failure, WgsFailure):
                failure.add_secondary(secondary)
            else:
                # Keep a raw interruption/domain error's original type until
                # its owning boundary classifies it, carrying only safe fields.
                values, truncated = _cleanup_details(failure)
                if len(values) < MAX_SECONDARY_FAILURES:
                    values += (SecondaryFailure(secondary.failure_code,
                        secondary.operation, secondary.exception_kind),)
                else:
                    truncated = True
                failure._wgs_secondary_failures = values
                failure._wgs_secondary_failures_truncated = truncated
    if primary is None and failure is not None:
        raise failure from None


@contextmanager
def managed_resource_session(manager):
    """Classify lifecycle failures without relabelling the managed operation."""

    primary = None
    entered = False
    resource = None
    try:
        with failure_boundary("resource_session", "resource_session_open"):
            resource = manager.__enter__()
            entered = True
        try:
            yield resource
        except (Exception, KeyboardInterrupt) as error:
            primary = static_failure(error)
            raise primary from None
    finally:
        if entered:
            try:
                manager.__exit__(None, None, None)
            except (Exception, KeyboardInterrupt) as error:
                cleanup = static_failure(
                    error,
                    "resource_session",
                    operation="resource_session_close",
                )
                if primary is not None:
                    primary.add_secondary(cleanup)
                else:
                    raise cleanup from None


@contextmanager
def failure_boundary(
    origin,
    operation,
    *,
    known_errors=(),
    category="malformed_evidence",
):
    """Label an actual boundary; outer scopes preserve prior classification."""

    if origin not in ORIGINS or operation not in OPERATIONS:
        raise ValueError("invalid static failure boundary")
    try:
        yield
    except WgsFailure:
        raise
    except (Exception, KeyboardInterrupt) as error:
        kind = exception_kind(error)
        known = isinstance(error, known_errors)
        if known and kind == "other":
            kind = "project_failure"
        classified_category = (
            "memory_exhaustion" if kind == "memory_error" else
            "interrupted" if kind == "interrupted" else
            "timeout" if isinstance(error, TimeoutError) else
            "storage" if kind == "os_error" else
            category if known else
            "unexpected"
        )
        secondary, truncated = _cleanup_details(error)
        raise WgsFailure(
            origin,
            classified_category,
            operation=operation,
            exception_kind=kind,
            secondary_failures=secondary,
            secondary_failures_truncated=truncated,
        ) from None


def _secondary(value: object) -> SecondaryFailure:
    if type(value) is not dict or set(value) != _SECONDARY_FIELDS:
        raise ValueError
    return SecondaryFailure(
        value["failure_code"], value["operation"], value["exception_kind"],
    )


def load_failure_evidence(
    content,
    *,
    receipt_content,
    expected_failure_code,
    expected_journal_heads,
):
    """Read only the current diagnostic and authenticate caller-known bindings."""

    from genosapiens.wgs_ledger import canonical

    try:
        if (
            type(content) is not bytes
            or not content
            or len(content) > MAX_FAILURE_EVIDENCE_BYTES
            or type(receipt_content) is not bytes
            or expected_failure_code is not None
            and expected_failure_code not in FAILURE_CODES
            or type(expected_journal_heads) not in {tuple, list}
            or len(expected_journal_heads) > 16
            or len(set(expected_journal_heads)) != len(expected_journal_heads)
            or any(not _digest(value) for value in expected_journal_heads)
        ):
            raise ValueError

        def unique(pairs):
            result = dict(pairs)
            if len(result) != len(pairs):
                raise ValueError
            return result

        value = json.loads(content, object_pairs_hook=unique)
        if type(value) is not dict or set(value) != _EVIDENCE_FIELDS:
            raise ValueError
        secondaries = value["secondary_failures"]
        decoded_secondaries = tuple(_secondary(item) for item in secondaries)
        if (
            value["schema"] != "genosapiens.wgs-failure-evidence"
            or type(value["version"]) is not int or value["version"] != 2
            or value["classification"] != "controlled"
            or value["llm_policy"] != "no-llm"
            or value["stage"] not in {"preparation", "analysis"}
            or value["receipt_sha256"] != hashlib.sha256(receipt_content).hexdigest()
            or value["failure_code"] not in FAILURE_CODES
            or expected_failure_code is not None
            and value["failure_code"] != expected_failure_code
            or value["operation"] not in OPERATIONS
            or value["exception_kind"] not in EXCEPTION_KINDS
            or type(value["journal_heads"]) is not list
            or value["journal_heads"] != list(expected_journal_heads)
            or type(secondaries) is not list
            or len(secondaries) > MAX_SECONDARY_FAILURES
            or type(value["secondary_failures_truncated"]) is not bool
            or value["secondary_failures_truncated"]
            and len(secondaries) != MAX_SECONDARY_FAILURES
            or canonical(value) != content
        ):
            raise ValueError
        if len(decoded_secondaries) != len(secondaries):
            raise ValueError
        return value
    except (ValueError, TypeError, KeyError, UnicodeError, json.JSONDecodeError):
        raise ValueError("invalid static failure evidence") from None


def terminal_failure_binding(workspace, *, stage, receipt_content):
    """Derive a failed journal code/head from receipt-bound journal bytes."""

    from genosapiens.wgs_ledger import (
        JournalEntry,
        MAX_JOURNAL_EVENTS,
        MAX_MANIFEST_BYTES,
        load_model,
        read_journal,
    )

    try:
        root = Path(workspace)
        receipt = json.loads(
            receipt_content,
            object_pairs_hook=lambda pairs: (
                dict(pairs)
                if len(pairs) == len({key for key, _value in pairs})
                else (_ for _ in ()).throw(ValueError())
            ),
        )
        if type(receipt) is not dict or stage not in {"preparation", "analysis"}:
            raise ValueError
        if stage == "analysis":
            return None, ()
        heads = receipt.get("completed_journal_heads")
        if (
            type(heads) is not list
            or len(heads) > 16
            or len(set(heads)) != len(heads)
            or any(not _digest(head) for head in heads)
        ):
            raise ValueError
        if not heads:
            return None, ()
        journals_root = root / "journals"
        found = {}
        directories = sorted(journals_root.iterdir())
        if len(directories) > MAX_FAILURE_JOURNAL_DIRECTORIES:
            raise ValueError
        for directory in directories:
            if directory.is_symlink() or not directory.is_dir():
                raise ValueError
            names = sorted(directory.glob("*.journal.json"))
            if not names:
                continue
            if len(names) > MAX_JOURNAL_EVENTS:
                raise ValueError
            content = _read_bounded_regular(names[0], MAX_MANIFEST_BYTES)
            first = load_model(content, JournalEntry)
            entries = read_journal(
                directory,
                stage=first.stage,
                parent_sha256=first.stage_parent_sha256,
            )
            head = entries[-1]
            if head.sha256 in heads:
                if head.sha256 in found:
                    raise ValueError
                found[head.sha256] = head
        if set(found) != set(heads):
            raise ValueError
        failed_heads = tuple(
            head for head in heads
            if found[head].event_type in {"attempt_failed", "attempt_interrupted"}
        )
        codes = {found[head].failure_code for head in failed_heads}
        if len(codes) > 1 or any(code not in FAILURE_CODES for code in codes):
            raise ValueError
        return (next(iter(codes)) if codes else None), failed_heads
    except (OSError, ValueError, TypeError, KeyError, json.JSONDecodeError):
        raise ValueError("invalid terminal failure binding") from None


def publish_failure_evidence(root, *, stage, receipt_content, failure):
    """Best effort: a diagnostic write must never replace the failed exit."""

    from genosapiens.wgs_ledger import canonical, publish_bytes

    if not isinstance(failure, WgsFailure):
        return None
    content = canonical({
        "schema": "genosapiens.wgs-failure-evidence", "version": 2,
        "classification": "controlled", "llm_policy": "no-llm", "stage": stage,
        "receipt_sha256": hashlib.sha256(receipt_content).hexdigest(),
        "failure_code": failure.failure_code,
        "operation": failure.operation,
        "exception_kind": failure.exception_kind,
        "journal_heads": failure.journal_heads,
        "secondary_failures": tuple(
            value.as_dict() for value in failure.secondary_failures
        ),
        "secondary_failures_truncated": failure.secondary_failures_truncated,
    })
    try:
        load_failure_evidence(
            content,
            receipt_content=receipt_content,
            expected_failure_code=failure.failure_code,
            expected_journal_heads=failure.journal_heads,
        )
        return publish_bytes(root, "failure-" + stage + ".json", content)
    except Exception as error:
        failure.add_secondary(static_failure(
            error,
            "final_publication",
            operation="diagnostic_publish",
        ))
        return None
