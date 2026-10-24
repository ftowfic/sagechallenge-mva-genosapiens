"""Local, append-only disclosure stops; no payload collection or network access."""

from __future__ import annotations

from contextlib import contextmanager
from datetime import datetime, timezone
from enum import StrEnum
import json
import os
from pathlib import Path
import re
import stat
from types import MappingProxyType
from typing import Callable, Iterator, NoReturn, TypeVar
from uuid import uuid4

from genosapiens.governance import GovernanceError


# Official rules pinned by config/challenge-rules/mva-hackathon-2026-track2-v1.json.
NOTIFICATION_REQUIREMENT = (
    "Notify Sage Bionetworks' Privacy and Compliance Office via Sage's Help Center "
    "whenever unauthorized disclosure occurs or is suspected. The pinned rules "
    "do not specify a Help Center URL; obtain it from an official Sage source."
)
_NOTIFICATION = MappingProxyType({
    "notification_office": "sage_privacy_and_compliance_office",
    "notification_channel": "sage_help_center",
    "notification_help_center_url": None,
    "notification_url_state": "not_specified_in_pinned_rules",
    "notification_rule_sha256":
        "dccb3740f6301b8629c704c375739e23391d19dee62513587ee461d6c8637f8f",
})
_UUID4 = r"[0-9a-f]{12}4[0-9a-f]{3}[89ab][0-9a-f]{15}"
_READY = b"genosapiens.incident-environment/v1\n"
_MAX_RECORD_BYTES = 4096
_T = TypeVar("_T")


class EventCategory(StrEnum):
    CANARY_ALERT = "canary_alert"
    REPOSITORY_SCAN = "repository_scan"
    LOG_REVIEW = "log_review"
    OPERATOR_REPORT = "operator_report"
    QUERY_BOUNDARY = "query_boundary"
    OTHER_CONTROL = "other_control"


class ProcessingBlocked(GovernanceError):
    """Safe to display as text, without exception locals or caller arguments."""

    def __init__(self) -> None:
        super().__init__("affected_processing_blocked; " + NOTIFICATION_REQUIREMENT)


class DisclosureAlert(Exception):
    """Signal from a detector to execute_guarded; never attach detector content."""

    def __init__(self, category: EventCategory) -> None:
        # Invalid alert metadata must still stop and latch a guarded operation.
        self.category = category
        super().__init__("suspected_disclosure")


def _validate_id(value: str, prefix: str) -> None:
    if type(value) is not str or re.fullmatch(prefix + _UUID4, value) is None:
        raise ProcessingBlocked() from None


@contextmanager
def _root(incident_root: str | Path) -> Iterator[int]:
    """Anchor operations to a private directory, rejecting every symlink component."""
    descriptor = -1
    try:
        path = Path(incident_root)
        if not path.is_absolute() or ".." in path.parts:
            raise ProcessingBlocked() from None
        flags = os.O_RDONLY | os.O_DIRECTORY | os.O_NOFOLLOW
        descriptor = os.open(path.anchor, flags)
        for component in path.parts[1:]:
            child = os.open(component, flags, dir_fd=descriptor)
            os.close(descriptor)
            descriptor = child
        metadata = os.fstat(descriptor)
        if metadata.st_uid != os.geteuid() or stat.S_IMODE(metadata.st_mode) != 0o700:
            raise ProcessingBlocked() from None
        yield descriptor
    except (OSError, ValueError, TypeError, RecursionError):
        raise ProcessingBlocked() from None
    finally:
        if descriptor >= 0:
            os.close(descriptor)


def _read(descriptor: int, name: str) -> bytes:
    # Unlike a plain JSON open, O_NONBLOCK cannot hang on a substituted FIFO.
    child = os.open(
        name, os.O_RDONLY | os.O_NOFOLLOW | os.O_NONBLOCK, dir_fd=descriptor
    )
    try:
        metadata = os.fstat(child)
        if (
            not stat.S_ISREG(metadata.st_mode)
            or metadata.st_nlink != 1
            or metadata.st_uid != os.geteuid()
            or stat.S_IMODE(metadata.st_mode) != 0o600
            or metadata.st_size > _MAX_RECORD_BYTES
        ):
            raise ProcessingBlocked() from None
        with os.fdopen(child, "rb", closefd=False) as handle:
            content = handle.read(_MAX_RECORD_BYTES + 1)
        if len(content) > _MAX_RECORD_BYTES:
            raise ProcessingBlocked() from None
        return content
    finally:
        os.close(child)


def _write(descriptor: int, name: str, content: bytes) -> None:
    child = os.open(
        name, os.O_WRONLY | os.O_CREAT | os.O_EXCL | os.O_NOFOLLOW,
        0o600, dir_fd=descriptor,
    )
    try:
        with os.fdopen(child, "wb", closefd=False) as handle:
            handle.write(content)
            handle.flush()
            os.fsync(child)
        os.fsync(descriptor)
    finally:
        os.close(child)


def _json(record: dict) -> bytes:
    return (json.dumps(record, sort_keys=True, separators=(",", ":")) + "\n").encode()


def _timestamp() -> str:
    return datetime.now(timezone.utc).isoformat(timespec="microseconds")


def register_environment(incident_root: str | Path) -> str:
    """Provision once in an existing private root; persist the returned ID externally."""
    with _root(incident_root) as descriptor:
        environment_id = "env_" + uuid4().hex
        _write(descriptor, environment_id + ".ready", _READY)
        return environment_id


def assert_processing_allowed(incident_root: str | Path, environment_id: str) -> None:
    """Deny missing, unreadable, unregistered or stopped state; never initialize it."""
    _validate_id(environment_id, "env_")
    with _root(incident_root) as descriptor:
        if _read(descriptor, environment_id + ".ready") != _READY:
            raise ProcessingBlocked() from None
        try:
            os.stat(environment_id + ".halted", dir_fd=descriptor, follow_symlinks=False)
        except FileNotFoundError:
            return
        # Existence alone blocks, including empty, damaged or symlink stop markers.
        raise ProcessingBlocked() from None


def report_suspected_disclosure(
    incident_root: str | Path,
    environment_id: str,
    *,
    category: EventCategory,
) -> NoReturn:
    """Durably latch, save only fixed-schema evidence, then ALWAYS raise ProcessingBlocked.

    Even invalid category metadata leaves a registered environment latched. A failed
    evidence write never removes a stop marker. No reset or notification bypass exists.
    """
    _validate_id(environment_id, "env_")
    with _root(incident_root) as descriptor:
        try:
            if _read(descriptor, environment_id + ".ready") != _READY:
                raise ProcessingBlocked() from None
        except (FileNotFoundError, ProcessingBlocked):
            # Another reporter may unlink readiness after this reporter opened
            # it, so _read rejects an inode whose link count became zero. Only
            # an existing stop marker permits appending evidence in that state.
            # Repeated reports remain possible after readiness has been revoked.
            os.stat(environment_id + ".halted", dir_fd=descriptor, follow_symlinks=False)
        try:
            _write(descriptor, environment_id + ".halted", b"")
        except FileExistsError:
            os.fsync(descriptor)
        try:
            os.unlink(environment_id + ".ready", dir_fd=descriptor)
        except FileNotFoundError:
            pass
        os.fsync(descriptor)
        if type(category) is not EventCategory:
            raise ProcessingBlocked() from None
        incident_id = "inc_" + uuid4().hex
        record = {
            "schema": "genosapiens.disclosure-incident/v1",
            "incident_id": incident_id,
            "environment_id": environment_id,
            "event_category": category.value,
            "detected_at": _timestamp(),
            "processing_state": "halted",
            "notification_state": "required",
            **_NOTIFICATION,
        }
        _write(descriptor, incident_id + ".json", _json(record))
    raise ProcessingBlocked() from None


def record_manual_notification(
    incident_root: str | Path, environment_id: str, incident_id: str
) -> None:
    """Append a participant's sent attestation, not proof of receipt; NEVER unblock."""
    _validate_id(environment_id, "env_")
    _validate_id(incident_id, "inc_")
    with _root(incident_root) as descriptor:
        content = _read(descriptor, incident_id + ".json")
        record = json.loads(content)
        expected_keys = {
            "schema", "incident_id", "environment_id", "event_category",
            "detected_at", "processing_state", "notification_state", *_NOTIFICATION,
        }
        if type(record) is not dict or set(record) != expected_keys:
            raise ProcessingBlocked() from None
        if (
            record["schema"] != "genosapiens.disclosure-incident/v1"
            or record["incident_id"] != incident_id
            or record["environment_id"] != environment_id
            or record["processing_state"] != "halted"
            or record["notification_state"] != "required"
            or type(record["event_category"]) is not str
            or record["event_category"] not in EventCategory
            or any(record[key] != value for key, value in _NOTIFICATION.items())
            or type(record["detected_at"]) is not str
            or re.fullmatch(
                r"\d{4}-\d{2}-\d{2}T\d{2}:\d{2}:\d{2}\.\d{6}\+00:00",
                record["detected_at"], flags=re.ASCII,
            ) is None
            or _json(record) != content
        ):
            raise ProcessingBlocked() from None
        detected_at = datetime.fromisoformat(record["detected_at"])
        sent_at = _timestamp()
        if detected_at > datetime.fromisoformat(sent_at):
            raise ProcessingBlocked() from None
        # Neither a notification nor a malformed marker can clear a stop.
        os.stat(environment_id + ".halted", dir_fd=descriptor, follow_symlinks=False)
        _write(descriptor, incident_id + ".notification.json", _json({
            "schema": "genosapiens.disclosure-notification/v1",
            "incident_id": incident_id,
            "environment_id": environment_id,
            "notification_state": "participant_attested_sent",
            "attested_at": sent_at,
            **_NOTIFICATION,
        }))


def execute_guarded(
    incident_root: str | Path, environment_id: str, operation: Callable[[], _T]
) -> _T:
    """Check before/after one cooperative step and latch a DisclosureAlert immediately."""
    assert_processing_allowed(incident_root, environment_id)
    try:
        result = operation()
    except DisclosureAlert as alert:
        report_suspected_disclosure(incident_root, environment_id, category=alert.category)
    except Exception:
        # Never echo arbitrary callback errors, paths, arguments or query contents.
        raise ProcessingBlocked() from None
    assert_processing_allowed(incident_root, environment_id)
    return result
