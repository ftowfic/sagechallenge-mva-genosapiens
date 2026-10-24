"""Bounded host checks, phase aggregation, and cgroup-v2 runtime evidence."""

from __future__ import annotations

from contextlib import contextmanager
from dataclasses import asdict, dataclass, fields, replace
import hashlib
import json
import os
from pathlib import Path
import signal
import stat
import subprocess
import sys
import threading
import time

from genosapiens.bounded_wgs_config import BoundedWgsExecutionConfig


class WgsRuntimeError(RuntimeError):
    """Static capacity/runtime failure without host paths or case-derived values."""


@dataclass(frozen=True, slots=True)
class FixedHostCapacity:
    host_available_memory_bytes: int
    filesystem_available_bytes: int
    projected_scratch_bytes: int
    required_filesystem_bytes: int


def require_fixed_host_capacity(
    config: BoundedWgsExecutionConfig,
    *,
    host_available_memory_bytes: int,
    allow_low_host_memory: bool = False,
    filesystem_available_bytes: int = 0,
    projected_scratch_bytes: int = 0,
) -> FixedHostCapacity:
    """Require fixed memory and container limits. Legacy disk arguments are inert."""

    if (
        type(config) is not BoundedWgsExecutionConfig
        or type(allow_low_host_memory) is not bool
        or type(host_available_memory_bytes) is not int
        or host_available_memory_bytes < 0
        or config.runtime.container_memory_bytes != 24 * 1024**3
        or config.runtime.container_memory_swap_bytes
        != config.runtime.container_memory_bytes
        or config.runtime.host_available_memory_bytes != 26 * 1024**3
        or config.runtime.container_cpus != 12
        or config.runtime.pids_limit != 512
        or config.runtime.tmpfs_bytes != 1024**3
        or config.runtime.network != "none"
        or config.runtime.non_root_uid_gid != "65532:65532"
        or config.runtime.read_only_root is not True
    ):
        raise WgsRuntimeError("bounded runtime contract is invalid")
    if host_available_memory_bytes < config.runtime.host_available_memory_bytes and not allow_low_host_memory:
        raise WgsRuntimeError("insufficient host memory for bounded execution")
    return FixedHostCapacity(
        host_available_memory_bytes,
        filesystem_available_bytes,
        projected_scratch_bytes,
        0,
    )


RUNTIME_SAMPLE_SCHEMA = "genosapiens.wgs-bounded-runtime-sample"
RUNTIME_EVIDENCE_SCHEMA = "genosapiens.wgs-bounded-runtime-evidence"
RUNTIME_PHASE_SCHEMA = "genosapiens.wgs-bounded-runtime-phases"
MAX_RUNTIME_DOCUMENT_BYTES = 64 * 1024
MAX_CGROUP_FILE_BYTES = 64 * 1024
MAX_PHASE_DEPTH = 8
RUNTIME_STAGE_ENV = "GENOSAPIENS_RUNTIME_STAGE"
RUNTIME_WORKSPACE_ENV = "GENOSAPIENS_RUNTIME_WORKSPACE"
RUNTIME_ATTEMPT_ENV = "GENOSAPIENS_RUNTIME_ATTEMPT_ID"
RUNTIME_CONFIG_ENV = "GENOSAPIENS_RUNTIME_CONFIG_SHA256"
STAGES = frozenset({"preparation", "analysis"})
PHASES = (
    "authority_validation",
    "input_validation",
    "reference_initialization",
    "resource_initialization",
    "hpo_initialization",
    "public_index_reconciliation",
    "source_ingestion",
    "source_accounting",
    "source_audit",
    "prefilter_session_authentication",
    "prefilter_request",
    "prefilter_query",
    "prefilter_worker_query",
    "prefilter_adaptation",
    "filter_publication",
    "filter_index_commit",
    "filter_closure",
    "working_set_publication",
    "working_set_closure",
    "preparation_terminal_publication",
    "prepared_authority_validation",
    "analysis_resource_initialization",
    "annotation_planning",
    "annotation_execution",
    "annotation_response_reconciliation",
    "admission_projection",
    "context_formation",
    "scientific_analysis",
    "sensitivity_analysis",
    "report_publication",
    "final_publication",
    "analysis_terminal_publication",
)
_PHASE_ORDER = {name: ordinal for ordinal, name in enumerate(PHASES)}
_CORE_EVENT_NAMES = ("low", "high", "max", "oom", "oom_kill")
_EVENT_NAMES = (*_CORE_EVENT_NAMES, "oom_group_kill")


def _canonical(value: object) -> bytes:
    return (
        json.dumps(
            value,
            sort_keys=True,
            separators=(",", ":"),
            ensure_ascii=True,
            allow_nan=False,
        )
        + "\n"
    ).encode("ascii")


def _count(value: object, minimum: int = 0) -> bool:
    return (
        type(value) is int
        and not isinstance(value, bool)
        and minimum <= value <= 2**63 - 1
    )


def _optional_count(value: object, minimum: int = 0) -> bool:
    return value is None or _count(value, minimum)


def _digest(value: object) -> bool:
    return (
        type(value) is str
        and len(value) == 64
        and all(character in "0123456789abcdef" for character in value)
    )


@dataclass(frozen=True, slots=True)
class CgroupEvents:
    low: int
    high: int
    max: int
    oom: int
    oom_kill: int
    oom_group_kill: int | None = None

    def __post_init__(self) -> None:
        if (
            any(not _count(getattr(self, name)) for name in _CORE_EVENT_NAMES)
            or not _optional_count(self.oom_group_kill)
        ):
            raise WgsRuntimeError("cgroup event counters are invalid")


@dataclass(frozen=True, slots=True)
class CgroupSnapshot:
    memory_current_bytes: int
    memory_peak_bytes: int
    memory_limit_bytes: int | None
    anon_bytes: int
    file_bytes: int
    shmem_bytes: int
    events: CgroupEvents

    def __post_init__(self) -> None:
        if (
            any(not _count(value) for value in (
                self.memory_current_bytes,
                self.memory_peak_bytes,
                self.anon_bytes,
                self.file_bytes,
                self.shmem_bytes,
            ))
            or not _optional_count(self.memory_limit_bytes, 1)
            or self.memory_peak_bytes < self.memory_current_bytes
            or self.shmem_bytes > self.file_bytes
            or type(self.events) is not CgroupEvents
        ):
            raise WgsRuntimeError("cgroup snapshot is invalid")


def _event_binding(
    start: CgroupEvents | None,
    end: CgroupEvents | None,
) -> tuple[CgroupEvents | None, str]:
    if start is None or end is None:
        return None, "unavailable"
    values = {}
    for name in _EVENT_NAMES:
        before, after = getattr(start, name), getattr(end, name)
        if before is None or after is None:
            values[name] = None
        elif after < before:
            return None, "counter_reset"
        else:
            values[name] = after - before
    return CgroupEvents(**values), "measured"


@dataclass(frozen=True, slots=True)
class RuntimePhaseAggregate:
    phase: str
    entry_count: int
    completed_count: int
    incomplete_count: int
    elapsed_ns: int
    sample_count: int
    sampled_memory_current_max_bytes: int | None
    sampled_anon_max_bytes: int | None
    sampled_file_max_bytes: int | None
    sampled_shmem_max_bytes: int | None
    event_start: CgroupEvents | None
    event_end: CgroupEvents | None
    event_delta: CgroupEvents | None
    event_quality: str

    def __post_init__(self) -> None:
        sampled = self.sample_count > 0
        if (
            self.phase not in _PHASE_ORDER
            or not _count(self.entry_count, 1)
            or not _count(self.completed_count)
            or not _count(self.incomplete_count)
            or self.completed_count + self.incomplete_count > self.entry_count
            or not _count(self.elapsed_ns)
            or not _count(self.sample_count)
            or any(
                not _optional_count(value)
                for value in (
                    self.sampled_memory_current_max_bytes,
                    self.sampled_anon_max_bytes,
                    self.sampled_file_max_bytes,
                    self.sampled_shmem_max_bytes,
                )
            )
            or sampled != all(value is not None for value in (
                self.sampled_memory_current_max_bytes,
                self.sampled_anon_max_bytes,
                self.sampled_file_max_bytes,
                self.sampled_shmem_max_bytes,
            ))
            or self.sampled_shmem_max_bytes is not None
            and self.sampled_shmem_max_bytes > self.sampled_file_max_bytes
            or self.event_quality not in {
                "measured", "counter_reset", "unavailable",
            }
            or self.event_quality == "measured"
            and any(value is None for value in (
                self.event_start, self.event_end, self.event_delta,
            ))
            or self.event_quality == "counter_reset"
            and (
                self.event_start is None or self.event_end is None
                or self.event_delta is not None
            )
            or self.event_quality == "unavailable"
            and any(value is not None for value in (
                self.event_start, self.event_end, self.event_delta,
            ))
        ):
            raise WgsRuntimeError("runtime phase aggregate is invalid")


@dataclass(frozen=True, slots=True)
class RuntimeMetrics:
    sample_count: int
    sampled_memory_current_max_bytes: int | None
    kernel_memory_peak_bytes: int | None
    memory_limit_bytes: int | None
    sampled_anon_max_bytes: int | None
    sampled_file_max_bytes: int | None
    sampled_shmem_max_bytes: int | None
    event_start: CgroupEvents | None
    event_end: CgroupEvents | None
    event_delta: CgroupEvents | None
    event_quality: str
    cgroup_quality: str
    scratch_start_bytes: int | None
    scratch_high_water_bytes: int | None
    scratch_growth_high_water_bytes: int | None
    scratch_method: str | None
    scratch_quality: str
    phases: tuple[RuntimePhaseAggregate, ...]
    phase_coverage: str
    phase_checkpoint_sha256: str | None
    cgroup_interval_ms: int = 200
    scratch_interval_seconds: int = 60
    cgroup_sample_failures: int = 0
    cgroup_first_sample_ns: int | None = None
    cgroup_last_sample_ns: int | None = None
    scratch_scan_count: int = 0
    scratch_scan_failures: int = 0
    scratch_first_scan_ns: int | None = None
    scratch_last_scan_ns: int | None = None
    scratch_coverage: str = "unavailable"

    def __post_init__(self) -> None:
        sampled = self.sample_count > 0
        memory_values = (
            self.sampled_memory_current_max_bytes,
            self.kernel_memory_peak_bytes,
            self.sampled_anon_max_bytes,
            self.sampled_file_max_bytes,
            self.sampled_shmem_max_bytes,
        )
        scratch = self.scratch_high_water_bytes is not None
        if (
            self.cgroup_interval_ms != 200 or type(self.cgroup_interval_ms) is not int
            or self.scratch_interval_seconds != 60 or type(self.scratch_interval_seconds) is not int
            or not _count(self.cgroup_sample_failures)
            or not _count(self.scratch_scan_count) or not _count(self.scratch_scan_failures)
            or self.scratch_coverage not in {"start-and-exit", "partial", "unavailable"}
            or self.scratch_coverage == "start-and-exit" and (self.scratch_scan_count < 2 or self.scratch_scan_failures)
            or (self.scratch_scan_count > 0) != scratch
            or sampled != (self.cgroup_first_sample_ns is not None and self.cgroup_last_sample_ns is not None)
            or scratch != (self.scratch_first_scan_ns is not None and self.scratch_last_scan_ns is not None)
            or any(not _optional_count(v) for v in (self.cgroup_first_sample_ns, self.cgroup_last_sample_ns,
                self.scratch_first_scan_ns, self.scratch_last_scan_ns))
            or sampled and self.cgroup_first_sample_ns > self.cgroup_last_sample_ns
            or scratch and self.scratch_first_scan_ns > self.scratch_last_scan_ns
            or not _count(self.sample_count)
            or any(not _optional_count(value) for value in memory_values)
            or not _optional_count(self.memory_limit_bytes, 1)
            or sampled != all(value is not None for value in memory_values)
            or sampled != (self.cgroup_quality == "measured")
            or self.cgroup_quality not in {"measured", "unavailable"}
            or not sampled and self.memory_limit_bytes is not None
            or sampled
            and self.kernel_memory_peak_bytes < self.sampled_memory_current_max_bytes
            or sampled
            and self.sampled_shmem_max_bytes > self.sampled_file_max_bytes
            or self.event_quality not in {
                "measured", "counter_reset", "unavailable",
            }
            or self.event_quality == "measured"
            and any(value is None for value in (
                self.event_start, self.event_end, self.event_delta,
            ))
            or self.event_quality == "counter_reset"
            and (
                self.event_start is None or self.event_end is None
                or self.event_delta is not None
            )
            or self.event_quality == "unavailable"
            and any(value is not None for value in (
                self.event_start, self.event_end, self.event_delta,
            ))
            or not _optional_count(self.scratch_start_bytes)
            or not _optional_count(self.scratch_high_water_bytes)
            or not _optional_count(self.scratch_growth_high_water_bytes)
            or scratch != (self.scratch_start_bytes is not None)
            or scratch != (self.scratch_growth_high_water_bytes is not None)
            or scratch != (self.scratch_method is not None)
            or scratch != (self.scratch_quality == "measured")
            or self.scratch_quality not in {"measured", "unavailable"}
            or self.scratch_method not in {
                None, "allocated-blocks-recursive-sampled-high-water",
            }
            or scratch
            and (
                self.scratch_high_water_bytes < self.scratch_start_bytes
                or self.scratch_growth_high_water_bytes
                != self.scratch_high_water_bytes - self.scratch_start_bytes
            )
            or type(self.phases) is not tuple
            or len(self.phases) > len(PHASES)
            or any(type(value) is not RuntimePhaseAggregate for value in self.phases)
            or tuple(value.phase for value in self.phases)
            != tuple(sorted(
                {value.phase for value in self.phases}, key=_PHASE_ORDER.__getitem__,
            ))
            or self.phase_coverage not in {
                "complete", "incomplete", "unavailable", "invalid",
            }
            or self.phase_coverage in {"complete", "incomplete"}
            and (
                not self.phases or not _digest(self.phase_checkpoint_sha256)
                or any(
                    value.entry_count
                    != value.completed_count + value.incomplete_count
                    for value in self.phases
                )
            )
            or self.phase_coverage == "complete"
            and any(value.incomplete_count for value in self.phases)
            or self.phase_coverage == "unavailable"
            and (self.phases or self.phase_checkpoint_sha256 is not None)
            or self.phase_coverage == "invalid"
            and (self.phases or not _digest(self.phase_checkpoint_sha256))
        ):
            raise WgsRuntimeError("runtime metrics are invalid")


@dataclass(frozen=True, slots=True)
class RuntimePhaseCheckpoint:
    stage: str
    attempt_id: str
    execution_config_sha256: str
    elapsed_ns: int
    phases: tuple[RuntimePhaseAggregate, ...]
    active_phases: tuple[str, ...]
    terminal: bool
    schema: str = RUNTIME_PHASE_SCHEMA
    version: int = 2
    classification: str = "controlled"

    def __post_init__(self) -> None:
        active_counts = {
            phase: self.active_phases.count(phase) for phase in set(self.active_phases)
        }
        if (
            self.schema != RUNTIME_PHASE_SCHEMA
            or self.version != 2
            or self.classification != "controlled"
            or self.stage not in STAGES
            or not _digest(self.attempt_id)
            or not _digest(self.execution_config_sha256)
            or not _count(self.elapsed_ns)
            or type(self.phases) is not tuple
            or not self.phases
            or len(self.phases) > len(PHASES)
            or tuple(value.phase for value in self.phases)
            != tuple(sorted(
                {value.phase for value in self.phases}, key=_PHASE_ORDER.__getitem__,
            ))
            or type(self.active_phases) is not tuple
            or len(self.active_phases) > MAX_PHASE_DEPTH
            or any(value not in _PHASE_ORDER for value in self.active_phases)
            or type(self.terminal) is not bool
            or self.terminal and self.active_phases
            or any(
                value.entry_count
                != value.completed_count + value.incomplete_count
                + active_counts.get(value.phase, 0)
                for value in self.phases
            )
        ):
            raise WgsRuntimeError("runtime phase checkpoint is invalid")

    def as_bytes(self) -> bytes:
        content = _canonical(asdict(self))
        if len(content) > MAX_RUNTIME_DOCUMENT_BYTES:
            raise WgsRuntimeError("runtime phase checkpoint exceeds its bound")
        return content


@dataclass(frozen=True, slots=True)
class StageRuntimeSample:
    stage: str
    attempt_id: str
    execution_config_sha256: str
    child_exit_code: int
    elapsed_ns: int
    metrics: RuntimeMetrics
    schema: str = RUNTIME_SAMPLE_SCHEMA
    version: int = 3
    classification: str = "controlled"

    def __post_init__(self) -> None:
        if (
            self.schema != RUNTIME_SAMPLE_SCHEMA
            or self.version != 3
            or self.classification != "controlled"
            or self.stage not in STAGES
            or not _digest(self.attempt_id)
            or not _digest(self.execution_config_sha256)
            or not _count(self.child_exit_code)
            or self.child_exit_code > 255
            or not _count(self.elapsed_ns)
            or type(self.metrics) is not RuntimeMetrics
        ):
            raise WgsRuntimeError("bounded runtime sample is invalid")

    def as_bytes(self) -> bytes:
        content = _canonical(asdict(self))
        if len(content) > MAX_RUNTIME_DOCUMENT_BYTES:
            raise WgsRuntimeError("bounded runtime sample exceeds its bound")
        return content


@dataclass(frozen=True, slots=True)
class StageRuntimeEvidence:
    stage: str
    attempt_id: str
    execution_config_sha256: str
    exit_code: int
    oom_killed: bool | None
    termination: str
    metrics: RuntimeMetrics
    sample_sha256: str | None
    sample_quality: str
    schema: str = RUNTIME_EVIDENCE_SCHEMA
    version: int = 3
    classification: str = "controlled"

    def __post_init__(self) -> None:
        if (
            self.schema != RUNTIME_EVIDENCE_SCHEMA
            or self.version != 3
            or self.classification != "controlled"
            or self.stage not in STAGES
            or not _digest(self.attempt_id)
            or not _digest(self.execution_config_sha256)
            or not _count(self.exit_code)
            or self.exit_code > 255
            or self.oom_killed is not None and type(self.oom_killed) is not bool
            or self.termination not in {
                "completed", "nonzero-exit", "exit-137", "oom-killed",
            }
            or self.termination
            != (
                "oom-killed" if self.oom_killed is True else
                "exit-137" if self.exit_code == 137 else
                "completed" if self.exit_code == 0 else
                "nonzero-exit"
            )
            or type(self.metrics) is not RuntimeMetrics
            or self.sample_quality not in {"measured", "unavailable"}
            or (self.sample_quality == "measured") != _digest(self.sample_sha256)
            or self.sample_quality == "unavailable" and (
                self.metrics.cgroup_quality != "unavailable"
                or self.metrics.scratch_quality != "unavailable"
                or self.metrics.phase_coverage != "unavailable"
            )
        ):
            raise WgsRuntimeError("bounded runtime evidence is invalid")

    def as_bytes(self) -> bytes:
        content = _canonical(asdict(self))
        if len(content) > MAX_RUNTIME_DOCUMENT_BYTES:
            raise WgsRuntimeError("bounded runtime evidence exceeds its bound")
        return content


def _fields(kind) -> frozenset[str]:
    return frozenset(value.name for value in fields(kind))


def _decode_events(value) -> CgroupEvents | None:
    if value is None:
        return None
    if type(value) is not dict or set(value) != _fields(CgroupEvents):
        raise ValueError
    return CgroupEvents(**value)


def _decode_phase(value) -> RuntimePhaseAggregate:
    if type(value) is not dict or set(value) != _fields(RuntimePhaseAggregate):
        raise ValueError
    return RuntimePhaseAggregate(
        **{
            **value,
            "event_start": _decode_events(value["event_start"]),
            "event_end": _decode_events(value["event_end"]),
            "event_delta": _decode_events(value["event_delta"]),
        }
    )


def _decode_metrics(value) -> RuntimeMetrics:
    if type(value) is not dict or set(value) != _fields(RuntimeMetrics):
        raise ValueError
    phases = value["phases"]
    if type(phases) is not list or len(phases) > len(PHASES):
        raise ValueError
    return RuntimeMetrics(
        **{
            **value,
            "event_start": _decode_events(value["event_start"]),
            "event_end": _decode_events(value["event_end"]),
            "event_delta": _decode_events(value["event_delta"]),
            "phases": tuple(_decode_phase(item) for item in phases),
        }
    )


def _document(content: bytes) -> dict:
    if (
        type(content) is not bytes
        or not content
        or len(content) > MAX_RUNTIME_DOCUMENT_BYTES
    ):
        raise ValueError

    def unique(pairs):
        value = dict(pairs)
        if len(value) != len(pairs):
            raise ValueError
        return value

    value = json.loads(content, object_pairs_hook=unique)
    if type(value) is not dict or _canonical(value) != content:
        raise ValueError
    return value


def load_runtime_phase_checkpoint(
    content: bytes,
    *,
    expected_stage: str,
    expected_attempt_id: str,
    expected_execution_config_sha256: str,
) -> RuntimePhaseCheckpoint:
    try:
        value = _document(content)
        if set(value) != _fields(RuntimePhaseCheckpoint):
            raise ValueError
        phases, active = value["phases"], value["active_phases"]
        if type(phases) is not list or type(active) is not list:
            raise ValueError
        result = RuntimePhaseCheckpoint(
            **{
                **value,
                "phases": tuple(_decode_phase(item) for item in phases),
                "active_phases": tuple(active),
            }
        )
        if (
            result.stage != expected_stage
            or result.attempt_id != expected_attempt_id
            or result.execution_config_sha256
            != expected_execution_config_sha256
        ):
            raise ValueError
        return result
    except (ValueError, TypeError, KeyError, UnicodeError, json.JSONDecodeError,
            WgsRuntimeError):
        raise WgsRuntimeError("runtime phase checkpoint is invalid") from None


def load_stage_runtime_sample(
    content: bytes,
    *,
    expected_stage: str,
    expected_attempt_id: str,
    expected_execution_config_sha256: str,
    expected_exit_code: int | None = None,
) -> StageRuntimeSample:
    try:
        value = _document(content)
        if set(value) != _fields(StageRuntimeSample):
            raise ValueError
        result = StageRuntimeSample(**{
            **value,
            "metrics": _decode_metrics(value["metrics"]),
        })
        if (
            result.stage != expected_stage
            or result.attempt_id != expected_attempt_id
            or result.execution_config_sha256
            != expected_execution_config_sha256
            or expected_exit_code is not None
            and result.child_exit_code != expected_exit_code
        ):
            raise ValueError
        return result
    except (ValueError, TypeError, KeyError, UnicodeError, json.JSONDecodeError,
            WgsRuntimeError):
        raise WgsRuntimeError("bounded runtime sample is invalid") from None


def load_stage_runtime_evidence(
    content: bytes,
    *,
    expected_stage: str,
    expected_attempt_id: str | None = None,
    expected_execution_config_sha256: str | None = None,
) -> StageRuntimeEvidence:
    try:
        value = _document(content)
        if set(value) != _fields(StageRuntimeEvidence):
            raise ValueError
        result = StageRuntimeEvidence(**{
            **value,
            "metrics": _decode_metrics(value["metrics"]),
        })
        if (
            result.stage != expected_stage
            or expected_attempt_id is not None
            and result.attempt_id != expected_attempt_id
            or expected_execution_config_sha256 is not None
            and result.execution_config_sha256
            != expected_execution_config_sha256
        ):
            raise ValueError
        return result
    except (ValueError, TypeError, KeyError, UnicodeError, json.JSONDecodeError,
            WgsRuntimeError):
        raise WgsRuntimeError("bounded runtime evidence is invalid") from None


def _read_cgroup_text(path: Path) -> str:
    descriptor = -1
    try:
        descriptor = os.open(path, os.O_RDONLY | os.O_NOFOLLOW)
        content = bytearray()
        while len(content) <= MAX_CGROUP_FILE_BYTES:
            block = os.read(
                descriptor,
                min(8192, MAX_CGROUP_FILE_BYTES + 1 - len(content)),
            )
            if not block:
                break
            content.extend(block)
        if not content or len(content) > MAX_CGROUP_FILE_BYTES:
            raise OSError
        return bytes(content).decode("ascii", errors="strict")
    except (OSError, UnicodeError):
        raise WgsRuntimeError(
            "container cgroup measurement is unavailable"
        ) from None
    finally:
        if descriptor >= 0:
            os.close(descriptor)


def _named_counts(path: Path, required: tuple[str, ...]) -> dict[str, int]:
    values = {}
    try:
        for line in _read_cgroup_text(path).splitlines():
            parts = line.split()
            if len(parts) != 2 or parts[0] in values:
                raise ValueError
            value = int(parts[1])
            if value < 0:
                raise ValueError
            values[parts[0]] = value
    except (OSError, UnicodeError, ValueError):
        raise WgsRuntimeError("container cgroup measurement is unavailable") from None
    if any(name not in values for name in required):
        raise WgsRuntimeError("container cgroup measurement is unavailable")
    return values


def _read_count(path: Path, *, allow_max: bool = False) -> int | None:
    try:
        value = _read_cgroup_text(path).strip()
        if allow_max and value == "max":
            return None
        result = int(value)
        if result < 0 or allow_max and result == 0:
            raise ValueError
        return result
    except (OSError, UnicodeError, ValueError, WgsRuntimeError):
        raise WgsRuntimeError("container cgroup measurement is unavailable") from None


def read_cgroup_v2_snapshot(root: str | Path) -> CgroupSnapshot:
    """Read one bounded cgroup-v2 snapshot, ignoring unknown kernel fields."""

    root = Path(root)
    current = _read_count(root / "memory.current")
    peak = _read_count(root / "memory.peak")
    limit = _read_count(root / "memory.max", allow_max=True)
    memory = _named_counts(root / "memory.stat", ("anon", "file", "shmem"))
    events = _named_counts(root / "memory.events", _CORE_EVENT_NAMES)
    return CgroupSnapshot(
        current,
        peak,
        limit,
        memory["anon"],
        memory["file"],
        memory["shmem"],
        CgroupEvents(
            *(events[name] for name in _CORE_EVENT_NAMES),
            events.get("oom_group_kill"),
        ),
    )


def require_bounded_cgroup_limits(root: str | Path) -> None:
    """Reject a stage launched outside the fixed whole-container v5 quotas."""
    root = Path(root)
    try:
        cpu = _read_cgroup_text(root / "cpu.max").split()
        if len(cpu) != 2:
            raise ValueError
        quota, period = (int(value) for value in cpu)
        if (
            quota <= 0 or period <= 0 or quota != 12 * period
            or _read_count(root / "memory.max") != 24 * 1024**3
            or _read_count(root / "memory.swap.max") != 0
            or _read_count(root / "pids.max") != 512
        ):
            raise ValueError
    except (ValueError, TypeError, WgsRuntimeError):
        raise WgsRuntimeError("bounded container limits are incompatible") from None


def supervise_operator_run(command, *, workspace: Path, started: float,
                           maximum_seconds: float | None, kill_grace_seconds=30, lock_fd=None):
    """Supervise one complete Linux case, including descendants in new sessions.

    Start an independent interpreter. Stage recorders remain in the
    worker. The parent enforces the optional single deadline.
    """
    import ctypes
    import math
    from dataclasses import asdict
    if (not command or any(type(value) is not str or not value for value in command)
            or maximum_seconds is not None and (
            not math.isfinite(maximum_seconds) or maximum_seconds <= 0)
            or not 0 < kill_grace_seconds <= 30):
        raise WgsRuntimeError('operator supervision settings are invalid')
    # Adopt grandchildren if a worker exits, including children that called
    # setsid. A process-group kill alone cannot cover those children.
    libc = ctypes.CDLL(None, use_errno=True)
    if libc.prctl(36, 1, 0, 0, 0) != 0:  # PR_SET_CHILD_SUBREAPER
        raise WgsRuntimeError('operator child supervision is unavailable')
    cancelled = False
    def cancel(*_):
        nonlocal cancelled
        cancelled = True
    handlers = {sig: signal.signal(sig, cancel) for sig in (signal.SIGINT, signal.SIGTERM)}
    def descendants():
        parents = {}
        for entry in Path('/proc').iterdir():
            if entry.name.isdigit():
                try:
                    # comm may contain spaces or parentheses.
                    fields = (entry / 'stat').read_text().rsplit(')', 1)[1].split()
                    parents[int(entry.name)] = int(fields[1])
                except (OSError, ValueError, IndexError):
                    pass
        owned = {os.getpid()}
        while True:
            expanded = owned | {pid for pid, parent in parents.items() if parent in owned}
            if expanded == owned:
                return owned - {os.getpid()}
            owned = expanded
    existing_children = descendants()
    def owned_descendants():
        return descendants() - existing_children
    def signal_owned(sig):
        for pid in owned_descendants():
            try:
                os.kill(pid, sig)
            except ProcessLookupError:
                pass
    worker = None
    def reap_owned():
        for pid in owned_descendants() | ({worker} if worker is not None else set()):
            try:
                reaped, result = os.waitpid(pid, os.WNOHANG)
                if reaped == worker:
                    process.returncode = os.waitstatus_to_exitcode(result)
            except ChildProcessError:
                pass
    status, code, peak = 'failed', 1, None
    reason = None
    try:
        if maximum_seconds is not None and time.monotonic() - started >= maximum_seconds:
            status, code = 'timed-out', 124
            return code
        process = subprocess.Popen(command, start_new_session=True,
            pass_fds=() if lock_fd is None else (lock_fd,))
        worker = process.pid
        next_notice = started
        while True:
            waited, outcome = os.waitpid(worker, os.WNOHANG)
            if waited:
                code = os.waitstatus_to_exitcode(outcome)
                process.returncode = code
                if code < 0:
                    code = 128 - code
                status = {0: 'completed-ready', 2: 'completed-non-ready'}.get(code, 'failed')
                break
            elapsed = time.monotonic() - started
            if cancelled:
                status, code = 'cancelled', 130
                break
            if maximum_seconds is not None and elapsed >= maximum_seconds:
                status, code = 'timed-out', 124
                break
            snapshot = read_cgroup_v2_snapshot('/sys/fs/cgroup')
            peak = snapshot
            if time.monotonic() >= next_notice:
                print(f'Analysis running: {elapsed / 60:.1f} minutes elapsed; '
                      + ('no overall limit.' if maximum_seconds is None else 'configured hours limit active.'), flush=True)
                notice = operator_progress_notice(workspace.parent.parent)
                if notice:
                    print(notice, flush=True)
                next_notice = time.monotonic() + 30
            time.sleep(0.2)
    finally:
        if worker != 0:
            # Also clean up leaked descendants after a nominally successful exit.
            pending = owned_descendants()
            if pending:
                signal_owned(signal.SIGTERM)
                until = time.monotonic() + kill_grace_seconds
                while owned_descendants() and time.monotonic() < until:
                    reap_owned()
                    time.sleep(0.05)
                signal_owned(signal.SIGKILL)
                until = time.monotonic() + 2
                while owned_descendants() and time.monotonic() < until:
                    reap_owned()
                    time.sleep(0.05)
                reap_owned()
                if owned_descendants():
                    status, code, reason = 'failed', 125, 'child_shutdown_incomplete'
            for sig, previous in handlers.items():
                signal.signal(sig, previous)
            _publish_private(workspace / 'outcome.json', _canonical(dict(
                schema='genosapiens.operator-invocation-outcome', version=1,
                status=status, exit_code=code, reason=reason,
                child_exit_code=None if worker is None else process.returncode,
                elapsed_seconds=time.monotonic() - started, maximum_seconds=maximum_seconds,
                last_cgroup_sample=None if peak is None else asdict(peak))))
    return code


def _allocated_bytes(root: Path, *, excluded_roots: tuple[Path, ...] = ()) -> int:
    # Workers unlink completed temporary spools while the observer walks them.
    # A vanished child contributes no live allocation. Other traversal errors
    # must remain failures, including unreadable directories (os.walk's default
    # otherwise silently suppresses those errors).
    def walk_error(error):
        if not isinstance(error, FileNotFoundError):
            raise error

    try:
        if root.is_symlink() or not root.is_dir():
            raise ValueError
        total = root.lstat().st_blocks * 512
        for directory, directories, files in os.walk(root, followlinks=False, onerror=walk_error):
            parent = Path(directory)
            directories[:] = [name for name in directories if parent / name not in excluded_roots]
            for name in (*directories, *files):
                path = parent / name
                try:
                    metadata = path.lstat()
                except FileNotFoundError:
                    continue
                if stat.S_ISLNK(metadata.st_mode):
                    raise ValueError
                total += metadata.st_blocks * 512
        return total
    except (OSError, ValueError):
        raise WgsRuntimeError("scratch measurement is unavailable") from None


def _publish_private(path: Path, content: bytes) -> None:
    descriptor = -1
    try:
        descriptor = os.open(
            path,
            os.O_WRONLY | os.O_CREAT | os.O_EXCL | os.O_NOFOLLOW,
            0o600,
        )
        os.fchmod(descriptor, 0o600)
        view = memoryview(content)
        while view:
            written = os.write(descriptor, view)
            if written <= 0:
                raise OSError
            view = view[written:]
        os.fsync(descriptor)
    except OSError:
        raise WgsRuntimeError("runtime evidence publication failed") from None
    finally:
        if descriptor >= 0:
            os.close(descriptor)


def _replace_private(path: Path, content: bytes) -> None:
    temporary = path.with_name(path.name + ".tmp")
    try:
        if temporary.exists() or temporary.is_symlink():
            temporary.unlink()
        _publish_private(temporary, content)
        os.replace(temporary, path)
    except (OSError, WgsRuntimeError):
        try:
            temporary.unlink()
        except OSError:
            pass
        raise WgsRuntimeError("runtime checkpoint publication failed") from None


def _read_private_bounded(path: Path, maximum_bytes: int) -> bytes:
    descriptor = -1
    try:
        descriptor = os.open(path, os.O_RDONLY | os.O_NOFOLLOW)
        metadata = os.fstat(descriptor)
        if (
            not stat.S_ISREG(metadata.st_mode)
            or metadata.st_uid != os.geteuid()
            or stat.S_IMODE(metadata.st_mode) != 0o600
            or metadata.st_size <= 0
            or metadata.st_size > maximum_bytes
        ):
            raise OSError
        with os.fdopen(descriptor, "rb", closefd=True) as stream:
            descriptor = -1
            content = stream.read(maximum_bytes + 1)
        if len(content) != metadata.st_size:
            raise OSError
        return content
    except OSError:
        raise WgsRuntimeError("runtime evidence input is invalid") from None
    finally:
        if descriptor >= 0:
            os.close(descriptor)


def _empty_metrics() -> RuntimeMetrics:
    return RuntimeMetrics(
        0, None, None, None, None, None, None,
        None, None, None, "unavailable", "unavailable",
        None, None, None, None, "unavailable",
        (), "unavailable", None,
    )


# This vocabulary contains code-owned labels only, never biological identifiers.
_ANALYSIS_PROGRESS_PHASES = {
    'completed_input_reuse': {'metadata', 'staging', 'assessment_snapshot'},
    'association_discovery': {'working', 'clinical_query', 'regulatory_query', 'writing_associations', 'complete', 'failed'},
    'context_assembly': {'resource_prepare', 'candidate_construction', 'graph_support',
        'association_batches', 'baselines', 'index_closure', 'assessment_publication',
        'assessment_decisions', 'assessment_sources', 'context_merge'},
    'scientific_evaluation': {'resources', 'member_enrichment', 'genotype_facts', 'hypothesis_evidence',
        'logical_pairs', 'population', 'phenotype', 'scores', 'representation', 'eligibility',
        'ordering', 'context_publication', 'primary_seal', 'case_seal', 'evidence_finish', 'descriptor_publication'},
    'publication': {'scientific_input_staging', 'scientific_replay', 'final_scientific_replay', 'scientific_integrity',
        'reference_metadata', 'gene_resources', 'phenotype_resources', 'population_resources',
        'products', 'reports', 'final_artifacts', 'completed_science_binding', 'saved_scientific_projections',
        'saved_context_directory', 'saved_order', 'saved_identities', 'saved_selection_inputs', 'saved_scalar_facts',
        'saved_selection', 'selected_detail', 'selected_detail_resources', 'output_copy', 'publication_closure'},
    'final_inventory_check': {'required_paths'},
}
_ANALYSIS_COUNT_KEYS = frozenset({
    'metadata_completed', 'staging_files_completed', 'staging_files_total', 'staging_bytes_completed',
    'staging_bytes_total', 'source_accounting_inputs', 'authentication_scans', 'complete_domain_scans',
    'candidates_processed', 'candidates_retained', 'graph_completed', 'graph_total',
    'regulatory_inputs_completed', 'assessed_completed', 'assessed_total', 'associations_completed', 'associations_total', 'batches_committed',
    'baselines_completed', 'baselines_total', 'contexts_obtained', 'contexts_merged', 'contexts_emitted',
    'contexts_started', 'contexts_completed', 'contexts_total', 'current_context_ordinal', 'members',
    'hypotheses', 'logical_pairs_completed', 'logical_pairs_total', 'materialized_hypotheses',
    'assessment_decisions_completed', 'assessment_decisions_total',
    'assessment_sources_completed', 'assessment_sources_total', 'assessment_nodes',
    'files_completed', 'files_total', 'bytes_completed', 'bytes_total',
    'reused_scientific_contexts', 'selected_contexts_replayed', 'selected_contexts_read', 'selected_contexts_total',
    'scientific_resources_restored',
})
_ANALYSIS_WORKER_KEYS = frozenset({'submitted', 'received', 'outstanding', 'cancelled',
    'input_bytes', 'result_bytes', 'reserved_bytes', 'cpu_nanoseconds', 'busy_nanoseconds', 'idle_nanoseconds'})
_ANALYSIS_WRITER_KEYS = frozenset({'commits', 'checkpoints', 'commit_nanoseconds', 'checkpoint_nanoseconds',
    'supplied_bytes', 'pending_bytes', 'database_bytes', 'wal_bytes', 'shm_bytes'})
_ANALYSIS_RESOURCE_KEYS = frozenset({'memory_current_bytes', 'memory_peak_bytes', 'anon_bytes', 'file_bytes',
    'shmem_bytes', 'memory_limit_bytes', 'swap_limit_bytes', 'cpu_usage_usec', 'cpu_throttled_usec',
    'oom', 'oom_kill', 'disk_available_bytes'})


def operator_progress_notice(output: Path) -> str | None:
    """Display bounded existing counters, never arbitrary case text or an ETA."""
    path = output / 'analysis/analysis-work/state/analysis-progress.json'
    try:
        with path.open('rb') as stream:
            raw = stream.read(32 * 1024 + 1)
        if len(raw) > 32 * 1024:
            return None
        data = json.loads(raw)
        phase, subphase = data['phase'], data['subphase']
        if phase not in _ANALYSIS_PROGRESS_PHASES or subphase not in _ANALYSIS_PROGRESS_PHASES[phase]:
            return None
        counts = data['counts']
        parts = []
        for label, completed, total in (
                ('contexts', 'contexts_completed', 'contexts_total'),
                ('assessments', 'assessed_completed', 'assessed_total'),
                ('files', 'files_completed', 'files_total')):
            value, bound = counts.get(completed), counts.get(total)
            if type(value) is int and value >= 0:
                parts.append(f'{label}: {value}' +
                    (f'/{bound}' if type(bound) is int and bound >= value else ''))
        return f'Current stage: {phase}/{subphase}' + ('; ' + '; '.join(parts) if parts else '') + '.'
    except (OSError, ValueError, KeyError, TypeError, AttributeError):
        return None  # Observational telemetry never replaces the actual outcome.


def analysis_worker_metrics(metrics):
    """Project only existing bounded numeric worker totals."""
    return {
        'submitted':metrics.get('submitted_jobs'), 'received':metrics.get('received_jobs'),
        'outstanding':metrics.get('outstanding_jobs'), 'cancelled':metrics.get('cancelled_jobs'),
        'input_bytes':metrics.get('outstanding_input_bytes'), 'result_bytes':metrics.get('max_result_bytes'),
        'reserved_bytes':metrics.get('reserved_result_bytes'), 'cpu_nanoseconds':metrics.get('cpu_nanoseconds'),
        'busy_nanoseconds':sum(metrics.get('worker_busy_nanoseconds', ())),
        'idle_nanoseconds':sum(metrics.get('worker_idle_nanoseconds', ())),
    }


def analysis_progress(callback, phase, subphase, **values):
    """Observational callback: telemetry must never replace scientific failure."""
    if callback is not None:
        try:
            callback(phase, subphase, **values)
        except Exception:
            # No arbitrary callback exception text belongs in runtime evidence.
            owner = getattr(callback, '__self__', callback)
            if isinstance(owner, AnalysisProgressRecorder):
                owner.degraded = True


class AnalysisProgressRecorder:
    """One private bounded current snapshot, not an authority or event ledger."""
    def __init__(self, *, path, attempt_id, candidate_sha256, restart_stage,
                 clock=time.monotonic, utc=None):
        from collections import deque
        from datetime import datetime, timezone
        if (not _digest(attempt_id) or candidate_sha256 is not None and not _digest(candidate_sha256)
                or restart_stage not in {'analysis_entry', 'association_discovery', 'context_assembly', 'assessment_publication', 'scientific_publication'}):
            raise WgsRuntimeError('invalid analysis progress binding')
        self.path = Path(path)
        self.attempt_id, self.candidate_sha256 = attempt_id, candidate_sha256
        self.restart_stage = restart_stage
        self.clock = clock
        self.utc = utc or (lambda: datetime.now(timezone.utc).isoformat())
        self.started = self.phase_started = clock()
        self.last_write = self.last_observation = float('-inf')
        self.phase = self.subphase = self.scope = None
        self.counts = dict.fromkeys(_ANALYSIS_COUNT_KEYS)
        self.worker = dict.fromkeys(_ANALYSIS_WORKER_KEYS)
        self.writer = dict.fromkeys(_ANALYSIS_WRITER_KEYS)
        self.resources = dict.fromkeys(_ANALYSIS_RESOURCE_KEYS)
        self.reused_stages = ()
        self.state = 'running'
        self.last_work_at = None
        self.degraded = False
        self.observations = deque(maxlen=12)
        self.unit = self.completed = self.total = None
        self.initial_completed = None
        self.writes = 0

    @staticmethod
    def _numeric(values, allowed):
        if (type(values) is not dict or not set(values) <= allowed
                or any(value is not None and (type(value) is not int or value < 0)
                       for value in values.values())):
            raise WgsRuntimeError('invalid analysis progress counters')
        return values

    def __call__(self, phase, subphase, *, counts=None, worker=None, writer=None, resources=None,
                 unit=None, completed=None, total=None, scope=None, reused_stages=None, force=False):
        if self.state != 'running':
            return
        if (phase not in _ANALYSIS_PROGRESS_PHASES or subphase not in _ANALYSIS_PROGRESS_PHASES[phase]
                or scope is not None and (type(scope) is not int or scope < 0)
                or unit not in {None, 'bytes', 'files', 'candidates', 'contexts', 'inputs', 'associations', 'pairs', 'records'}
                or any(value is not None and (type(value) is not int or value < 0) for value in (completed, total))
                or total is not None and completed is not None and completed > total):
            raise WgsRuntimeError('invalid analysis progress scope')
        now = self.clock()
        changed = (phase, subphase, scope) != (self.phase, self.subphase, self.scope)
        if changed:
            self.phase, self.subphase, self.scope = phase, subphase, scope
            self.phase_started = now
            self.observations.clear()
            self.initial_completed = completed
            self.last_observation = float('-inf')
            self.unit = self.completed = self.total = None
        advanced = completed is not None and (self.completed is None or completed > self.completed)
        if self.completed is not None and completed is not None and completed < self.completed:
            raise WgsRuntimeError('analysis progress counter regressed')
        for name, values, allowed in (
                ('counts', counts, _ANALYSIS_COUNT_KEYS), ('worker', worker, _ANALYSIS_WORKER_KEYS),
                ('writer', writer, _ANALYSIS_WRITER_KEYS), ('resources', resources, _ANALYSIS_RESOURCE_KEYS)):
            if values is not None:
                checked = self._numeric(values, allowed)
                if name in {'counts', 'worker', 'writer'}:
                    previous = getattr(self, name)
                    work_keys = _ANALYSIS_COUNT_KEYS if name == 'counts' else {'received'} if name == 'worker' else {'supplied_bytes'}
                    checked_work = {key:value for key,value in checked.items() if key in work_keys}
                    advanced |= any(v is not None and (previous.get(k) is None or v > previous[k])
                                    for k, v in checked_work.items())
                getattr(self, name).update(checked)
        if reused_stages is not None:
            if type(reused_stages) is not tuple or not set(reused_stages) <= {'profile', 'association_discovery', 'admission'}:
                raise WgsRuntimeError('invalid reused analysis stages')
            self.reused_stages = reused_stages
        if unit is not None:
            self.unit, self.completed, self.total = unit, completed, total
        if advanced:
            self.last_work_at = self.utc()
        if now - self.last_observation >= 5 and self.completed is not None:
            self.observations.append((now, self.completed))
            self.last_observation = now
        self.flush(force=force or changed)

    def estimate(self, now):
        result = dict(quality='unknown', unit=self.unit, completed=self.completed, total=self.total,
            percent=None, remaining_seconds_range=None, window_seconds=None, reason='unknown_total',
            overall_remaining_seconds=None, missing_stages=([] if self.state == 'complete' else
                ['publication'] if self.phase in {'publication', 'final_inventory_check'} else
                ['downstream_science', 'publication']))
        if self.completed is None or self.total is None:
            return result
        result['percent'] = 100.0 if self.total == 0 else 100 * self.completed / self.total
        result['reason'] = 'insufficient_advancing_observations'
        points = tuple(self.observations)
        advances = sum(b[1] > a[1] for a, b in zip(points, points[1:]))
        if len(points) < 2 or advances < 2 or points[-1][0] - points[0][0] < 10:
            return result
        if points[-1][1] == points[-2][1] or now - points[-1][0] > 10:
            result['reason'] = 'stalled_window'
            return result
        elapsed = now - self.phase_started
        cumulative = (self.completed - (self.initial_completed or 0)) / elapsed if elapsed > 0 else 0
        window = points[-1][0] - points[0][0]
        rolling = (points[-1][1] - points[0][1]) / window
        if min(cumulative, rolling) <= 0:
            return result
        remaining = self.total - self.completed
        result.update(quality='heuristic', reason='nonuniform_workload_possible', window_seconds=window,
            remaining_seconds_range=sorted((remaining / cumulative, remaining / rolling)))
        return result

    def snapshot(self):
        now = self.clock()
        return dict(schema='genosapiens.analysis-progress', version=1, classification='controlled', llm_policy='no-llm',
            attempt_id=self.attempt_id, candidate_sha256=self.candidate_sha256, restart_stage=self.restart_stage,
            state=self.state, phase=self.phase, subphase=self.subphase, at_utc=self.utc(),
            elapsed_seconds=max(0, now-self.started), phase_elapsed_seconds=max(0, now-self.phase_started),
            last_work_at_utc=self.last_work_at, reused_stages=self.reused_stages,
            counts=self.counts, worker=self.worker, writer=self.writer, resources=self.resources,
            estimate=self.estimate(now), telemetry_quality='degraded' if self.degraded else 'available')

    def flush(self, *, force=False):
        now = self.clock()
        if not force and now - self.last_write < 5:
            return
        self.last_write = now
        temporary = self.path.with_suffix('.tmp')
        try:
            try:
                sample = read_cgroup_v2_snapshot('/sys/fs/cgroup')
                self.resources.update({key:getattr(sample,key) for key in (
                    'memory_current_bytes','memory_peak_bytes','anon_bytes','file_bytes','shmem_bytes','memory_limit_bytes')})
                self.resources.update(oom=sample.events.oom, oom_kill=sample.events.oom_kill)
            except (WgsRuntimeError, OSError):
                pass
            content = (json.dumps(self.snapshot(), sort_keys=True, separators=(',', ':'))+'\n').encode()
            if len(content) > 32 * 1024:
                raise WgsRuntimeError('analysis progress exceeds capacity')
            self.path.parent.mkdir(mode=0o700, parents=True, exist_ok=True)
            fd = os.open(temporary, os.O_WRONLY|os.O_CREAT|os.O_TRUNC|os.O_NOFOLLOW, 0o600)
            with os.fdopen(fd, 'wb') as out:
                out.write(content)
            os.replace(temporary, self.path)
            self.writes += 1
        except (OSError, ValueError, WgsRuntimeError):
            self.degraded = True
            try:
                temporary.unlink(missing_ok=True)
            except OSError:
                pass

    def finish(self, state):
        if state not in {'complete', 'failed', 'interrupted'}:
            raise WgsRuntimeError('invalid terminal analysis progress state')
        self.state = state
        self.flush(force=True)


class RuntimePhaseRecorder:
    """Fixed-key phase accumulator with one atomically replaced checkpoint."""

    def __init__(
        self,
        *,
        stage: str,
        workspace: str | Path,
        attempt_id: str,
        execution_config_sha256: str,
        cgroup_root: str | Path = "/sys/fs/cgroup",
    ) -> None:
        if (
            stage not in STAGES
            or not _digest(attempt_id)
            or not _digest(execution_config_sha256)
        ):
            raise WgsRuntimeError("runtime phase recorder binding is invalid")
        self.stage = stage
        self.root = Path(workspace)
        self.attempt_id = attempt_id
        self.execution_config_sha256 = execution_config_sha256
        self.cgroup_root = Path(cgroup_root)
        self.path = self.root / f".runtime-phases-{stage}.json"
        self.started_ns = time.monotonic_ns()
        self._values = {}
        self._stack = []
        self._terminal = False
        self._telemetry_failed = False

    def _snapshot(self):
        try:
            return read_cgroup_v2_snapshot(self.cgroup_root)
        except WgsRuntimeError:
            return None

    def _value(self, phase):
        return self._values.setdefault(phase, {
            "entry_count": 0,
            "completed_count": 0,
            "incomplete_count": 0,
            "elapsed_ns": 0,
            "sample_count": 0,
            "current": None,
            "anon": None,
            "file": None,
            "shmem": None,
            "event_start": None,
            "event_end": None,
        })

    @staticmethod
    def _sample(value, snapshot):
        if snapshot is None:
            return
        value["sample_count"] += 1
        for key, observed in (
            ("current", snapshot.memory_current_bytes),
            ("anon", snapshot.anon_bytes),
            ("file", snapshot.file_bytes),
            ("shmem", snapshot.shmem_bytes),
        ):
            value[key] = observed if value[key] is None else max(value[key], observed)
        if value["event_start"] is None:
            value["event_start"] = snapshot.events
        value["event_end"] = snapshot.events

    def _aggregates(self):
        result = []
        for phase in sorted(self._values, key=_PHASE_ORDER.__getitem__):
            value = self._values[phase]
            delta, quality = _event_binding(
                value["event_start"], value["event_end"],
            )
            result.append(RuntimePhaseAggregate(
                phase,
                value["entry_count"],
                value["completed_count"],
                value["incomplete_count"],
                value["elapsed_ns"],
                value["sample_count"],
                value["current"], value["anon"], value["file"], value["shmem"],
                value["event_start"] if quality != "unavailable" else None,
                value["event_end"] if quality != "unavailable" else None,
                delta,
                quality,
            ))
        return tuple(result)

    def record_worker_duration(self, phase: str, elapsed_ns: int) -> None:
        """Accumulate summed child work; no per-batch document or memory sample."""
        if self._terminal or phase not in {"prefilter_worker_query", "prefilter_adaptation"} or not _count(elapsed_ns):
            raise WgsRuntimeError("worker phase aggregate is invalid")
        value = self._value(phase)
        value["entry_count"] += 1
        value["completed_count"] += 1
        value["elapsed_ns"] += elapsed_ns

    def _checkpoint(self):
        if not self._values or self._telemetry_failed:
            return
        try:
            checkpoint = RuntimePhaseCheckpoint(
                self.stage,
                self.attempt_id,
                self.execution_config_sha256,
                max(0, time.monotonic_ns() - self.started_ns),
                self._aggregates(),
                tuple(value[0] for value in self._stack),
                self._terminal,
            )
            _replace_private(self.path, checkpoint.as_bytes())
        except WgsRuntimeError:
            self._telemetry_failed = True

    @contextmanager
    def phase(self, phase: str):
        if (
            self._terminal
            or phase not in _PHASE_ORDER
            or len(self._stack) >= MAX_PHASE_DEPTH
        ):
            raise WgsRuntimeError("runtime phase transition is invalid")
        value = self._value(phase)
        value["entry_count"] += 1
        started = time.monotonic_ns()
        self._sample(value, self._snapshot())
        self._stack.append((phase, started))
        self._checkpoint()
        try:
            yield
        except BaseException:
            try:
                self._end(phase, started, completed=False)
            except WgsRuntimeError:
                pass
            raise
        else:
            self._end(phase, started, completed=True)

    def _end(self, phase, started, *, completed):
        if not self._stack or self._stack[-1] != (phase, started):
            raise WgsRuntimeError("runtime phase nesting is invalid")
        self._stack.pop()
        value = self._value(phase)
        value["elapsed_ns"] += max(0, time.monotonic_ns() - started)
        value["completed_count" if completed else "incomplete_count"] += 1
        self._sample(value, self._snapshot())
        self._checkpoint()

    def finish(self) -> RuntimePhaseCheckpoint:
        if (
            self._terminal
            or self._stack
            or not self._values
            or self._telemetry_failed
        ):
            raise WgsRuntimeError("runtime phase finalization is invalid")
        self._terminal = True
        self._checkpoint()
        if self._telemetry_failed:
            raise WgsRuntimeError("runtime phase finalization is invalid")
        return load_runtime_phase_checkpoint(
            _read_private_bounded(self.path, MAX_RUNTIME_DOCUMENT_BYTES),
            expected_stage=self.stage,
            expected_attempt_id=self.attempt_id,
            expected_execution_config_sha256=self.execution_config_sha256,
        )


def runtime_phase_recorder_from_environment(
    *, expected_stage: str,
) -> RuntimePhaseRecorder | None:
    """Construct the worker-side recorder from observer-owned exact bindings."""

    names = (
        RUNTIME_STAGE_ENV,
        RUNTIME_WORKSPACE_ENV,
        RUNTIME_ATTEMPT_ENV,
        RUNTIME_CONFIG_ENV,
    )
    values = tuple(os.environ.get(name) for name in names)
    if all(value is None for value in values):
        return None
    if (
        any(value is None or value == "" for value in values)
        or values[0] != expected_stage
        or expected_stage not in STAGES
        or not Path(values[1]).is_absolute()
    ):
        raise WgsRuntimeError("runtime observer environment is invalid")
    return RuntimePhaseRecorder(
        stage=expected_stage,
        workspace=values[1],
        attempt_id=values[2],
        execution_config_sha256=values[3],
    )


@contextmanager
def managed_runtime_phase_recorder(*, expected_stage: str):
    """Finalize worker telemetry without replacing a primary stage failure."""

    recorder = runtime_phase_recorder_from_environment(
        expected_stage=expected_stage,
    )
    try:
        yield recorder
    except BaseException as primary:
        if recorder is not None:
            from genosapiens.wgs_failures import run_cleanup
            run_cleanup([('runtime_finalize', recorder.finish)], primary=primary)
        raise
    else:
        if recorder is not None:
            try:
                recorder.finish()
            except WgsRuntimeError:
                pass


@contextmanager
def runtime_phase(recorder: RuntimePhaseRecorder | None, phase: str):
    """Use one allowlisted phase when a bounded observer owns this process."""

    if recorder is None:
        yield
    else:
        with recorder.phase(phase):
            yield


@dataclass(frozen=True, slots=True)
class AuditedRuntimeDurations:
    audited_total_elapsed_ns: int | None
    source_audit_elapsed_ns: int | None
    independently_measured_default_elapsed_ns: None = None
    current_consumer_elapsed_ns: int | None = None
    original_source_elapsed_ns: int | None = None
    original_source_audit_elapsed_ns: int | None = None
    original_phase_checkpoint_sha256: str | None = None
    completed_source_input_sha256: str | None = None



def _original_source_durations(completed_source):
    from genosapiens.wgs_preparation import ResolvedCompletedSource
    if type(completed_source) is not ResolvedCompletedSource:
        raise WgsRuntimeError("original runtime source is invalid")
    descriptor = completed_source.descriptor
    if descriptor.source_runtime_checkpoint_path is None:
        raise WgsRuntimeError("original runtime checkpoint binding is unavailable")
    content = _read_private_bounded(Path(descriptor.source_runtime_checkpoint_path),
                                   MAX_RUNTIME_DOCUMENT_BYTES)
    value = _document(content)
    checkpoint = load_runtime_phase_checkpoint(content, expected_stage="preparation",
        expected_attempt_id=value.get("attempt_id"),
        expected_execution_config_sha256=descriptor.producer.execution_config_sha256)
    phases = {row.phase: row for row in checkpoint.phases}
    required = ("source_ingestion", "source_accounting", "source_audit")
    if not descriptor.audit_requested or any(name not in phases or
            phases[name].completed_count != 1 or phases[name].incomplete_count for name in required):
        raise WgsRuntimeError("original completed stage runtime is unavailable")
    return phases["source_ingestion"].elapsed_ns, phases["source_audit"].elapsed_ns, hashlib.sha256(content).hexdigest()


def audited_runtime_durations(
    sample: StageRuntimeSample,
    evidence: StageRuntimeEvidence,
    *, completed_source=None,
) -> AuditedRuntimeDurations:
    """Report actual audited total and audit phase, without inferring a default run."""
    if (
        type(sample) is not StageRuntimeSample
        or type(evidence) is not StageRuntimeEvidence
        or sample.stage != "preparation"
        or evidence.stage != "preparation"
        or sample.attempt_id != evidence.attempt_id
        or sample.execution_config_sha256 != evidence.execution_config_sha256
        or sample.child_exit_code != evidence.exit_code
        or sample.metrics != evidence.metrics
        or evidence.sample_sha256 != hashlib.sha256(sample.as_bytes()).hexdigest()
        or evidence.sample_quality != "measured"
    ):
        raise WgsRuntimeError("audited runtime binding is invalid")
    if completed_source is not None:
        if any(row.phase in {"source_ingestion", "source_accounting", "source_audit"} for row in sample.metrics.phases):
            raise WgsRuntimeError("imported source cannot claim fresh source runtime")
        source_ns, audit_ns, checkpoint_sha = _original_source_durations(completed_source)
        return AuditedRuntimeDurations(None, None, current_consumer_elapsed_ns=sample.elapsed_ns,
            original_source_elapsed_ns=source_ns, original_source_audit_elapsed_ns=audit_ns,
            original_phase_checkpoint_sha256=checkpoint_sha,
            completed_source_input_sha256=completed_source.descriptor.sha256)
    audit = next((phase for phase in sample.metrics.phases
                  if phase.phase == "source_audit"), None)
    if (audit is None or audit.completed_count != 1
            or audit.incomplete_count != 0
            or audit.elapsed_ns > sample.elapsed_ns):
        raise WgsRuntimeError("audited runtime phase is unavailable")
    return AuditedRuntimeDurations(sample.elapsed_ns, audit.elapsed_ns)


@dataclass(frozen=True, slots=True)
class PreparationScratchException:
    """Exact historical preparation exception; never an analysis launch permit."""

    receipt_sha256: str
    runtime_sha256: str
    producer_candidate_sha256: str
    original_analysis_attempt_id: str

    def __post_init__(self):
        if not all(_digest(value) for value in (
                self.receipt_sha256, self.runtime_sha256,
                self.producer_candidate_sha256, self.original_analysis_attempt_id)):
            raise WgsRuntimeError("preparation scratch exception is invalid")


def require_capacity_qualified_runtime(
    evidence: StageRuntimeEvidence,
    *,
    required_phases: tuple[str, ...] = (),
    completed_source=None,
    preparation_scratch_exception: PreparationScratchException | None = None,
) -> StageRuntimeEvidence:
    """Reject failed, incomplete, or unknown telemetry at a capacity gate."""

    if (
        type(evidence) is not StageRuntimeEvidence
        or type(required_phases) is not tuple
        or len(set(required_phases)) != len(required_phases)
        or any(value not in _PHASE_ORDER for value in required_phases)
    ):
        raise WgsRuntimeError("runtime capacity qualification request is invalid")
    if preparation_scratch_exception is not None:
        if (type(preparation_scratch_exception) is not PreparationScratchException
                or evidence.stage != "preparation"
                or hashlib.sha256(evidence.as_bytes()).hexdigest()
                    != preparation_scratch_exception.runtime_sha256
                or evidence.metrics.scratch_coverage != "partial"):
            raise WgsRuntimeError("preparation exception does not cover this runtime")
    if completed_source is not None:
        _original_source_durations(completed_source)
        if evidence.stage == "preparation":
            required_phases = tuple(name for name in required_phases
                if name not in {"source_ingestion", "source_accounting", "source_audit"})
            if any(row.phase in {"source_ingestion", "source_accounting", "source_audit"} for row in evidence.metrics.phases):
                raise WgsRuntimeError("imported source cannot claim fresh source capacity")
    metrics = evidence.metrics
    phases = {value.phase: value for value in metrics.phases}
    delta = metrics.event_delta
    if (
        evidence.exit_code != 0
        or evidence.oom_killed is not False
        or evidence.termination != "completed"
        or evidence.sample_quality != "measured"
        or metrics.cgroup_quality != "measured"
        or metrics.event_quality != "measured"
        or metrics.scratch_quality != "measured"
        or metrics.scratch_coverage != "start-and-exit" and preparation_scratch_exception is None
        or metrics.cgroup_sample_failures != 0
        or metrics.phase_coverage != "complete"
        or metrics.kernel_memory_peak_bytes is None
        or metrics.memory_limit_bytes is None
        or metrics.memory_limit_bytes != 24 * 1024**3
        or metrics.kernel_memory_peak_bytes > metrics.memory_limit_bytes
        or delta is None
        or delta.oom != 0
        or delta.oom_kill != 0
        or delta.oom_group_kill not in {None, 0}
        or any(
            phase not in phases
            or phases[phase].completed_count == 0
            or phases[phase].incomplete_count != 0
            for phase in required_phases
        )
    ):
        raise WgsRuntimeError("runtime capacity evidence is incomplete or failed")
    return evidence


class _MetricAccumulator:
    def __init__(self):
        self.count = 0
        self.current = self.peak = self.limit = None
        self.anon = self.file = self.shmem = None
        self.first_events = self.last_events = None
        self.stable = True

    def add(self, snapshot):
        if snapshot is None or not self.stable:
            return
        if self.count and snapshot.memory_limit_bytes != self.limit:
            self.stable = False
            self.count = 0
            self.current = self.peak = self.limit = None
            self.anon = self.file = self.shmem = None
            self.first_events = self.last_events = None
            return
        self.count += 1
        self.limit = snapshot.memory_limit_bytes
        for name, observed in (
            ("current", snapshot.memory_current_bytes),
            ("peak", snapshot.memory_peak_bytes),
            ("anon", snapshot.anon_bytes),
            ("file", snapshot.file_bytes),
            ("shmem", snapshot.shmem_bytes),
        ):
            previous = getattr(self, name)
            setattr(self, name, observed if previous is None else max(previous, observed))
        if self.first_events is None:
            self.first_events = snapshot.events
        self.last_events = snapshot.events


def _checkpoint_metrics(root, *, stage, attempt_id, config_sha256):
    path = root / f".runtime-phases-{stage}.json"
    if not path.exists() and not path.is_symlink():
        return (), "unavailable", None
    try:
        content = _read_private_bounded(path, MAX_RUNTIME_DOCUMENT_BYTES)
        checkpoint = load_runtime_phase_checkpoint(
            content,
            expected_stage=stage,
            expected_attempt_id=attempt_id,
            expected_execution_config_sha256=config_sha256,
        )
    except WgsRuntimeError:
        try:
            content = _read_private_bounded(path, MAX_RUNTIME_DOCUMENT_BYTES)
        except WgsRuntimeError:
            return (), "invalid", hashlib.sha256(b"invalid-runtime-checkpoint").hexdigest()
        return (), "invalid", hashlib.sha256(content).hexdigest()
    phases = checkpoint.phases
    if checkpoint.terminal:
        coverage = (
            "incomplete"
            if any(value.incomplete_count for value in phases)
            else "complete"
        )
    else:
        active = {
            phase: checkpoint.active_phases.count(phase)
            for phase in set(checkpoint.active_phases)
        }
        phases = tuple(replace(
            value,
            incomplete_count=value.incomplete_count + active.get(value.phase, 0),
        ) for value in phases)
        coverage = "incomplete"
    return phases, coverage, hashlib.sha256(content).hexdigest()


def _monitor_every(stop, interval_seconds, action):
    """Wait after each observation so slow work cannot cause catch-up scans."""
    while not stop.wait(interval_seconds):
        action()


def observe_bounded_stage(
    *,
    stage: str,
    workspace: str | Path,
    attempt_id: str,
    execution_config_sha256: str,
    command: tuple[str, ...],
    cgroup_root: str | Path = "/sys/fs/cgroup",
    cgroup_interval_ms: int = 200,
    scratch_interval_seconds: int = 60,
    excluded_roots: tuple[Path, ...] = (),
    require_bounded_limits: bool = False,
) -> tuple[int, StageRuntimeSample]:
    """Run a child while observing cgroup-v2 and bounded phase checkpoints."""

    root = Path(workspace)
    if (
        stage not in STAGES
        or not _digest(attempt_id)
        or not _digest(execution_config_sha256)
        or not command
        or any(type(value) is not str or not value for value in command)
        or type(cgroup_interval_ms) is not int or cgroup_interval_ms != 200
        or type(scratch_interval_seconds) is not int or scratch_interval_seconds != 60
        or type(excluded_roots) is not tuple
        or any(not path.is_absolute() or ".." in path.parts or path == root for path in excluded_roots)
        or type(require_bounded_limits) is not bool
    ):
        raise WgsRuntimeError("bounded stage observation request is invalid")
    if require_bounded_limits:
        require_bounded_cgroup_limits(cgroup_root)
    metrics = _MetricAccumulator()
    scratch_start = scratch_peak = None
    cgroup_first = cgroup_last = scratch_first = scratch_last = None
    scratch_count = scratch_failures = cgroup_failures = 0
    lock = threading.Lock()
    stop = threading.Event()
    scratch_stop = threading.Event()
    started = time.monotonic_ns()

    def sample_cgroup() -> None:
        nonlocal cgroup_first, cgroup_last, cgroup_failures
        try:
            snapshot = read_cgroup_v2_snapshot(cgroup_root)
        except WgsRuntimeError:
            cgroup_failures += 1
            return
        observed = max(0, time.monotonic_ns() - started)
        with lock:
            metrics.add(snapshot)
            cgroup_first = observed if cgroup_first is None else cgroup_first
            cgroup_last = observed

    def scan_scratch() -> None:
        nonlocal scratch_start, scratch_peak, scratch_count, scratch_failures, scratch_first, scratch_last
        try:
            used = _allocated_bytes(root, excluded_roots=excluded_roots)
        except WgsRuntimeError:
            scratch_failures += 1
            return
        observed = max(0, time.monotonic_ns() - started)
        scratch_count += 1
        scratch_first = observed if scratch_first is None else scratch_first
        scratch_last = observed
        scratch_start = used if scratch_start is None else scratch_start
        scratch_peak = used if scratch_peak is None else max(scratch_peak, used)

    def monitor_cgroup() -> None:
        _monitor_every(stop, cgroup_interval_ms / 1000, sample_cgroup)

    def monitor_scratch() -> None:
        # No catch-up scans: at least one full interval after each completed scan.
        _monitor_every(scratch_stop, scratch_interval_seconds, scan_scratch)

    sample_cgroup()
    thread = threading.Thread(target=monitor_cgroup, name="bounded-cgroup-observer")
    scratch_thread = threading.Thread(target=monitor_scratch, name="bounded-scratch-observer")
    thread.start()
    try:
        scan_scratch()
        environment = os.environ.copy()
        environment.update({
            RUNTIME_STAGE_ENV: stage,
            RUNTIME_WORKSPACE_ENV: str(root),
            RUNTIME_ATTEMPT_ENV: attempt_id,
            RUNTIME_CONFIG_ENV: execution_config_sha256,
        })
        try:
            process = subprocess.Popen(command, start_new_session=True, env=environment)
        except (OSError, ValueError, subprocess.SubprocessError):
            raise WgsRuntimeError("bounded stage child could not start") from None
        scratch_thread.start()
        try:
            return_code = process.wait()
        except KeyboardInterrupt:
            try:
                os.killpg(process.pid, signal.SIGINT)
            except ProcessLookupError:
                pass
            return_code = process.wait()
    finally:
        scratch_stop.set()
        if scratch_thread.ident is not None:
            scratch_thread.join()
        scan_scratch()
        stop.set()
        thread.join()
        sample_cgroup()
    exit_code = return_code if return_code >= 0 else 128 + (-return_code)
    elapsed = max(0, time.monotonic_ns() - started)
    with lock:
        delta, event_quality = _event_binding(
            metrics.first_events, metrics.last_events,
        )
        phases, phase_coverage, checkpoint_sha256 = _checkpoint_metrics(
            root,
            stage=stage,
            attempt_id=attempt_id,
            config_sha256=execution_config_sha256,
        )
        runtime_metrics = RuntimeMetrics(
            metrics.count,
            metrics.current,
            metrics.peak,
            metrics.limit if metrics.count else None,
            metrics.anon,
            metrics.file,
            metrics.shmem,
            metrics.first_events if event_quality != "unavailable" else None,
            metrics.last_events if event_quality != "unavailable" else None,
            delta,
            event_quality,
            "measured" if metrics.count else "unavailable",
            scratch_start,
            scratch_peak,
            None if scratch_peak is None else scratch_peak - scratch_start,
            None if scratch_peak is None else "allocated-blocks-recursive-sampled-high-water",
            "measured" if scratch_peak is not None else "unavailable",
            phases,
            phase_coverage,
            checkpoint_sha256,
            cgroup_interval_ms, scratch_interval_seconds,
            cgroup_failures, cgroup_first, cgroup_last, scratch_count, scratch_failures,
            scratch_first, scratch_last,
            "start-and-exit" if scratch_count >= 2 and not scratch_failures else
            "partial" if scratch_count else "unavailable",
        )
    sample_record = StageRuntimeSample(
        stage,
        attempt_id,
        execution_config_sha256,
        exit_code,
        elapsed,
        runtime_metrics,
    )
    try:
        _publish_private(
            root / f".runtime-sample-{stage}.json", sample_record.as_bytes(),
        )
    except WgsRuntimeError:
        # The stopped child's exit remains authoritative. Finalization records
        # unavailable telemetry rather than changing that primary exit.
        pass
    return exit_code, sample_record


def finalize_bounded_runtime(
    *,
    stage: str,
    workspace: str | Path,
    attempt_id: str,
    exit_code: int,
    oom_killed: bool | None,
    execution_config_sha256: str,
) -> StageRuntimeEvidence:
    """Publish terminal metrics before the stopped container is removed."""

    if (
        stage not in STAGES
        or not _digest(attempt_id)
        or not _digest(execution_config_sha256)
        or not _count(exit_code)
        or exit_code > 255
    ):
        raise WgsRuntimeError("runtime finalization binding is invalid")
    root = Path(workspace)
    sample_path = root / f".runtime-sample-{stage}.json"
    sample = None
    content = None
    if sample_path.exists() or sample_path.is_symlink():
        content = _read_private_bounded(sample_path, MAX_RUNTIME_DOCUMENT_BYTES)
        sample = load_stage_runtime_sample(
            content,
            expected_stage=stage,
            expected_attempt_id=attempt_id,
            expected_execution_config_sha256=execution_config_sha256,
            expected_exit_code=exit_code,
        )
    termination = (
        "oom-killed" if oom_killed is True else
        "exit-137" if exit_code == 137 else
        "completed" if exit_code == 0 else
        "nonzero-exit"
    )
    evidence = StageRuntimeEvidence(
        stage,
        attempt_id,
        execution_config_sha256,
        exit_code,
        oom_killed,
        termination,
        _empty_metrics() if sample is None else sample.metrics,
        None if content is None else hashlib.sha256(content).hexdigest(),
        "unavailable" if sample is None else "measured",
    )
    _publish_private(root / f"runtime-{stage}.json", evidence.as_bytes())
    return evidence
