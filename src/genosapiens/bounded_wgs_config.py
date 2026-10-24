"""Strict v5 technical limits for bounded whole-genome execution."""

from __future__ import annotations

from dataclasses import asdict, dataclass
import hashlib
import json
import re
from typing import Any, Literal


SCHEMA = "genosapiens.bounded-wgs-execution"
VERSION = 5
GIB = 1024 ** 3
_SHA256 = re.compile(r"^[0-9a-f]{64}$")
_ERRORS = (ValueError, TypeError, KeyError, UnicodeError, RecursionError, OverflowError)


class BoundedWgsConfigError(ValueError):
    """Static execution-contract error."""


def _fail(message: str = "invalid bounded WGS execution configuration") -> None:
    raise BoundedWgsConfigError(message)


def _pairs(items: list[tuple[str, Any]]) -> dict[str, Any]:
    value: dict[str, Any] = {}
    for key, item in items:
        if key in value:
            _fail()
        value[key] = item
    return value


def _constant(_value: str) -> None:
    _fail()


def _canonical(value: object) -> bytes:
    try:
        return (json.dumps(value, sort_keys=True, separators=(",", ":"),
                           ensure_ascii=True, allow_nan=False) + "\n").encode("ascii")
    except _ERRORS:
        _fail("cannot encode bounded WGS execution configuration")


def _object(value: object, fields: set[str]) -> dict[str, Any]:
    if type(value) is not dict or set(value) != fields:
        _fail()
    return value


def _positive(value: object) -> int:
    if type(value) is not int or value <= 0 or value > 2**63 - 1:
        _fail()
    return value


@dataclass(frozen=True, slots=True)
class InputLimits:
    max_record_bytes: int
    max_header_bytes: int
    strict_ascii: Literal[True]

    def __post_init__(self) -> None:
        _positive(self.max_record_bytes)
        _positive(self.max_header_bytes)
        if self.strict_ascii is not True or self.max_header_bytes < self.max_record_bytes:
            _fail()


@dataclass(frozen=True, slots=True)
class SourceShardLimits:
    max_uncompressed_rows: int
    max_uncompressed_bytes: int

    def __post_init__(self) -> None:
        _positive(self.max_uncompressed_rows)
        _positive(self.max_uncompressed_bytes)


@dataclass(frozen=True, slots=True)
class PublicQueryLimits:
    batch_rows: int
    concurrent_batches: int
    max_response_bytes: int
    timeout_seconds: int

    def __post_init__(self) -> None:
        _positive(self.batch_rows)
        _positive(self.concurrent_batches)
        _positive(self.max_response_bytes)
        _positive(self.timeout_seconds)
        if self.concurrent_batches != 10:
            _fail()


@dataclass(frozen=True, slots=True)
class AnnotationLimits:
    max_partition_requests: int
    max_partition_bytes: int
    buffer_size: int
    max_response_row_bytes: int
    max_shard_rows: int
    max_shard_bytes: int
    fork_processes: int

    def __post_init__(self) -> None:
        for value in asdict(self).values():
            _positive(value)
        if (self.buffer_size > self.max_partition_requests
                or self.max_response_row_bytes > self.max_shard_bytes
                or self.max_shard_bytes > self.max_partition_bytes
                or self.fork_processes != 4):
            _fail()


@dataclass(frozen=True, slots=True)
class ParallelismLimits:
    source_workers: int
    filter_workers: int
    max_inflight_cpu_batches: int
    max_inflight_filter_batches: int
    source_batch_rows: int
    max_source_batch_bytes: int

    def __post_init__(self) -> None:
        for value in asdict(self).values():
            _positive(value)
        if (self.source_workers != 4 or self.filter_workers != 10
                or self.max_inflight_cpu_batches != 4
                or self.max_inflight_filter_batches != 20
                or self.source_batch_rows != 256
                or self.max_source_batch_bytes != 64 * 1024 * 1024):
            _fail()


@dataclass(frozen=True, slots=True)
class MonitoringLimits:
    cgroup_sample_interval_ms: int
    scratch_scan_interval_seconds: int

    def __post_init__(self) -> None:
        _positive(self.cgroup_sample_interval_ms)
        _positive(self.scratch_scan_interval_seconds)
        if self.cgroup_sample_interval_ms != 200 or self.scratch_scan_interval_seconds != 60:
            _fail()


@dataclass(frozen=True, slots=True)
class VerificationLimits:
    source_replay_default: Literal[False]
    reauthenticate_completed_on_resume: Literal[False]

    def __post_init__(self) -> None:
        if (self.source_replay_default is not False
                or self.reauthenticate_completed_on_resume is not False):
            _fail()


@dataclass(frozen=True, slots=True)
class ChildProcessLimits:
    wall_timeout_seconds: int
    kill_grace_seconds: int
    max_total_output_bytes: int
    max_stderr_bytes: int
    pipe_chunk_bytes: int
    one_child_at_a_time: Literal[True]

    def __post_init__(self) -> None:
        for value in (
            self.wall_timeout_seconds, self.kill_grace_seconds, self.max_total_output_bytes,
            self.max_stderr_bytes, self.pipe_chunk_bytes,
        ):
            _positive(value)
        if (self.one_child_at_a_time is not True
                or self.kill_grace_seconds >= self.wall_timeout_seconds
                or self.max_stderr_bytes > self.max_total_output_bytes
                or self.pipe_chunk_bytes > self.max_stderr_bytes):
            _fail()


@dataclass(frozen=True, slots=True)
class ContextLimits:
    max_member_identities: int
    max_context_bytes: int
    pair_truncation_allowed: Literal[False]

    def __post_init__(self) -> None:
        _positive(self.max_member_identities)
        _positive(self.max_context_bytes)
        if self.pair_truncation_allowed is not False:
            _fail()


@dataclass(frozen=True, slots=True)
class SqliteLimits:
    cursor_batch_rows: int
    cache_bytes: int
    temp_store: Literal["file"]
    foreign_keys: Literal[True]
    load_extensions: Literal[False]

    def __post_init__(self) -> None:
        _positive(self.cursor_batch_rows)
        _positive(self.cache_bytes)
        if (self.temp_store != "file" or self.foreign_keys is not True
                or self.load_extensions is not False):
            _fail()


@dataclass(frozen=True, slots=True)
class RuntimeLimits:
    container_memory_bytes: int
    container_memory_swap_bytes: int
    host_available_memory_bytes: int
    container_cpus: int
    pids_limit: int
    tmpfs_bytes: int
    network: Literal["none"]
    non_root_uid_gid: Literal["65532:65532"]
    read_only_root: Literal[True]

    def __post_init__(self) -> None:
        for value in (
            self.container_memory_bytes, self.container_memory_swap_bytes,
            self.host_available_memory_bytes, self.container_cpus, self.pids_limit,
            self.tmpfs_bytes,
        ):
            _positive(value)
        if (self.container_memory_bytes != 24 * GIB
                or self.container_memory_swap_bytes != self.container_memory_bytes
                or self.host_available_memory_bytes != 26 * GIB
                or self.container_cpus != 12 or self.pids_limit != 512
                or self.tmpfs_bytes != GIB or self.network != "none"
                or self.non_root_uid_gid != "65532:65532" or self.read_only_root is not True
                or self.tmpfs_bytes >= self.container_memory_bytes):
            _fail()

    @property
    def additional_swap_bytes(self) -> int:
        return self.container_memory_swap_bytes - self.container_memory_bytes

    def operator_identity(self, uid: int, gid: int) -> str:
        """Bind the invoking operator; 65532 remains the public-test identity.

        This changes only local ownership. Resource and scientific policy
        configuration retains its qualified content identity.
        """
        if any(type(value) is not int or not 0 < value < 2**32 - 1 for value in (uid, gid)):
            _fail('operator execution requires a non-root UID and GID')
        return f'{uid}:{gid}'


@dataclass(frozen=True, slots=True)
class DiskLimits:
    """Retain legacy disk fields for checkpoint serialization, without capacity gates."""
    minimum_free_reserve_bytes: int
    scratch_safety_factor_numerator: int
    scratch_safety_factor_denominator: int
    projected_scratch_required: bool
    measure_scratch_high_water: Literal[True]

    def __post_init__(self) -> None:
        if (type(self.minimum_free_reserve_bytes) is not int or self.minimum_free_reserve_bytes < 0
                or type(self.scratch_safety_factor_numerator) is not int or self.scratch_safety_factor_numerator < 0):
            _fail()
        _positive(self.scratch_safety_factor_denominator)
        if (type(self.projected_scratch_required) is not bool
                or self.measure_scratch_high_water is not True):
            _fail()


@dataclass(frozen=True, slots=True)
class JournalLimits:
    default_scientific_attempts_per_partition: int
    durable_start_before_child: Literal[True]
    exclusive_stage_lock: Literal[True]

    def __post_init__(self) -> None:
        _positive(self.default_scientific_attempts_per_partition)
        if (self.durable_start_before_child is not True
                or self.exclusive_stage_lock is not True):
            _fail()


@dataclass(frozen=True, slots=True)
class BoundedWgsExecutionConfig:
    input: InputLimits
    source_shards: SourceShardLimits
    public_queries: PublicQueryLimits
    annotation: AnnotationLimits
    child_process: ChildProcessLimits
    contexts: ContextLimits
    sqlite: SqliteLimits
    runtime: RuntimeLimits
    disk: DiskLimits
    journal: JournalLimits
    parallelism: ParallelismLimits
    verification: VerificationLimits
    monitoring: MonitoringLimits
    schema: Literal["genosapiens.bounded-wgs-execution"] = SCHEMA
    version: Literal[5] = VERSION

    def __post_init__(self) -> None:
        if self.schema != SCHEMA or self.version != VERSION:
            _fail()
        if (self.source_shards.max_uncompressed_bytes < self.input.max_record_bytes
                or self.annotation.max_partition_bytes < self.input.max_record_bytes
                or self.child_process.max_total_output_bytes < self.annotation.max_partition_bytes
                or self.public_queries.batch_rows > self.contexts.max_member_identities
                or self.public_queries.concurrent_batches != self.parallelism.filter_workers
                or 2 * self.public_queries.concurrent_batches != self.parallelism.max_inflight_filter_batches):
            _fail()


_SECTIONS: dict[str, tuple[type, set[str]]] = {
    "monitoring": (MonitoringLimits, {"cgroup_sample_interval_ms", "scratch_scan_interval_seconds"}),
    "input": (InputLimits, {"max_record_bytes", "max_header_bytes", "strict_ascii"}),
    "source_shards": (SourceShardLimits, {"max_uncompressed_rows", "max_uncompressed_bytes"}),
    "public_queries": (PublicQueryLimits, {
        "batch_rows", "concurrent_batches", "max_response_bytes", "timeout_seconds",
    }),
    "annotation": (AnnotationLimits, {"max_partition_requests", "max_partition_bytes", "buffer_size",
                                      "max_response_row_bytes", "max_shard_rows", "max_shard_bytes",
                                      "fork_processes"}),
    "parallelism": (ParallelismLimits, {"source_workers", "filter_workers",
                                        "max_inflight_cpu_batches", "max_inflight_filter_batches", "source_batch_rows",
                                        "max_source_batch_bytes"}),
    "verification": (VerificationLimits, {"source_replay_default",
                                         "reauthenticate_completed_on_resume"}),
    "child_process": (ChildProcessLimits, {"wall_timeout_seconds", "kill_grace_seconds",
                                           "max_total_output_bytes", "max_stderr_bytes",
                                           "pipe_chunk_bytes", "one_child_at_a_time"}),
    "contexts": (ContextLimits, {"max_member_identities", "max_context_bytes",
                                 "pair_truncation_allowed"}),
    "sqlite": (SqliteLimits, {"cursor_batch_rows", "cache_bytes", "temp_store",
                              "foreign_keys", "load_extensions"}),
    "runtime": (RuntimeLimits, {"container_memory_bytes", "container_memory_swap_bytes",
                                "host_available_memory_bytes", "container_cpus", "pids_limit",
                                "tmpfs_bytes", "network", "non_root_uid_gid", "read_only_root"}),
    "journal": (JournalLimits, {"default_scientific_attempts_per_partition",
                                "durable_start_before_child", "exclusive_stage_lock"}),
}


def load_bounded_wgs_execution_config(
    content: bytes, expected_sha256: str
) -> BoundedWgsExecutionConfig:
    """Authenticate and load canonical technical bounds; no defaults are filled."""

    try:
        if (type(content) is not bytes or type(expected_sha256) is not str
                or _SHA256.fullmatch(expected_sha256) is None
                or hashlib.sha256(content).hexdigest() != expected_sha256):
            _fail()
        data = json.loads(content.decode("ascii"), object_pairs_hook=_pairs,
                          parse_constant=_constant)
        if _canonical(data) != content:
            _fail()
        top = _object(data, {"schema", "version", *_SECTIONS, "disk"})
        values: dict[str, object] = {}
        for name, (kind, fields) in _SECTIONS.items():
            values[name] = kind(**_object(top[name], fields))
        disk = _object(top["disk"], {"minimum_free_reserve_bytes", "scratch_safety_factor",
                                    "projected_scratch_required", "measure_scratch_high_water"})
        factor = _object(disk["scratch_safety_factor"], {"numerator", "denominator"})
        values["disk"] = DiskLimits(
            minimum_free_reserve_bytes=disk["minimum_free_reserve_bytes"],
            scratch_safety_factor_numerator=factor["numerator"],
            scratch_safety_factor_denominator=factor["denominator"],
            projected_scratch_required=disk["projected_scratch_required"],
            measure_scratch_high_water=disk["measure_scratch_high_water"],
        )
        return BoundedWgsExecutionConfig(
            **values, schema=top["schema"], version=top["version"]
        )
    except BoundedWgsConfigError:
        raise
    except _ERRORS:
        _fail()


def bounded_wgs_execution_config_as_json(config: BoundedWgsExecutionConfig) -> bytes:
    """Render canonical bytes matching the external v5 configuration shape."""

    if type(config) is not BoundedWgsExecutionConfig:
        _fail()
    value = asdict(config)
    disk = value["disk"]
    value["disk"] = {
        "minimum_free_reserve_bytes": disk["minimum_free_reserve_bytes"],
        "scratch_safety_factor": {
            "numerator": disk["scratch_safety_factor_numerator"],
            "denominator": disk["scratch_safety_factor_denominator"],
        },
        "projected_scratch_required": disk["projected_scratch_required"],
        "measure_scratch_high_water": disk["measure_scratch_high_water"],
    }
    return _canonical(value)
