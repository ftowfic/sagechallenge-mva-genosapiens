"""Bounded retained-allele cursors and deterministic VEP transport plans."""

from __future__ import annotations

from dataclasses import asdict, dataclass, replace
from contextlib import closing, contextmanager
import hashlib
import json
import multiprocessing
import os
from pathlib import Path
import queue
import re
import selectors
import signal
import socket
import sqlite3
import stat
import subprocess
import time
from typing import Iterable, Iterator

from genosapiens.bounded_wgs_config import (
    BoundedWgsExecutionConfig,
    VERSION as BOUNDED_WGS_CONFIG_VERSION,
    bounded_wgs_execution_config_as_json,
)
from genosapiens.evidence_model import SourcedRecord
from genosapiens.evidence_model import TranscriptConsequence, _typed, _typed_value, _validation_session
from genosapiens.evidence_serialization import _plain, _owned_canonical_bytes
from genosapiens.wgs_prepared_case import PreparedIndex, PREPARED_INDEX_VERSION
from genosapiens.vep_annotation import (
    EXECUTABLE as VEP_EXECUTABLE,
    FLAGS as VEP_FLAGS,
    RICH_FIELDS,
    ParsedVep,
    ParsedVepHeader,
    VepVcfRecord,
    parse_vep_header,
    parse_vep_vcf,
)
from genosapiens.wgs_ledger import (
    MAX_MANIFEST_BYTES,
    MAX_ROW_BYTES,
    MAX_SHARD_BYTES,
    MAX_SHARDS,
    SOURCE_MANIFEST_SCHEMA,
    WORKING_SET_MANIFEST_SCHEMA,
    Artifact,
    JournalEntry,
    LedgerError,
    LedgerManifest,
    ParentIdentity,
    StageCommitState,
    StageRunBinding,
    WorkingSetManifest,
    attest_existing_artifact,
    canonical,
    iter_artifact_lines,
    iter_ledger,
    iter_working_set_shard,
    load_model,
    publish_bytes,
    publish_stream,
    read_artifact,
    read_journal,
    stage_lock,
    verify_artifact,
)
from genosapiens.wgs_source_routing import (
    SourceRoutingContext,
    WgsSourceRoutingError,
    classify_source_occurrence,
    decode_source_occurrence,
)
from genosapiens.wgs_runtime import RuntimePhaseRecorder, runtime_phase
from genosapiens.wgs_annotation_coverage import ReferenceCoverage, REASON as UNASSESSED_REASON, POLICY as COVERAGE_POLICY
from genosapiens.wgs_annotation_coverage import ALIAS_REASON, ALIAS_POLICY, REASON_POLICIES


RETAINED_ALLELE_SCHEMA = "genosapiens.wgs-retained-allele"
VEP_REQUEST_SCHEMA = "genosapiens.wgs-vep-transport-request"
VEP_PARTITION_SCHEMA = "genosapiens.wgs-vep-transport-partition"
VEP_PLAN_SCHEMA = "genosapiens.wgs-vep-transport-plan"
ANNOTATION_CLOSURE_SCHEMA = "genosapiens.wgs-annotation-partition-closure"
ANNOTATION_STAGE_SCHEMA = "genosapiens.wgs-annotation-stage-manifest"
PARSED_PARTITION_SCHEMA = "genosapiens.wgs-parsed-annotation-partition"
VEP_COLUMNS = b"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
VEP_HEADER = b"##fileformat=VCFv4.2\n" + VEP_COLUMNS
_SHA256 = re.compile(r"^[0-9a-f]{64}$")
_NORMALIZED = re.compile(r"^allele-[0-9a-f]{24}$")
_SOURCE_ALLELE = re.compile(r"^source-allele-[0-9a-f]{24}$")
_IDENTIFIER = re.compile(r"^[A-Za-z0-9][A-Za-z0-9._:/-]{0,127}$")
_BASES = re.compile(r"^[ACGTN]+$")
_ANNOTATION_PARSE_WORKERS = 10
_ANNOTATION_PARSE_MAX_JOBS = 20
_ANNOTATION_PARSE_JOB_ROWS = 16
_ANNOTATION_PARSE_INPUT_BYTES = 16 * 1024**2
_ANNOTATION_PARSE_RESULT_BYTES = 64 * 1024**2


from genosapiens.wgs_failures import WgsFailure, static_failure
from genosapiens.wgs_query_limits import QueryBudget


class WgsAnnotationError(WgsFailure):
    """Static bounded-annotation preparation failure."""

    def __init__(self, message, *, category="malformed_evidence"):
        super().__init__(
            "annotation_child", category,
            operation="annotation_child_execute", exception_kind="project_failure",
        )
        RuntimeError.__init__(self, message)


def _digest(value: object) -> bool:
    return isinstance(value, str) and _SHA256.fullmatch(value) is not None


def _count(value: object, minimum: int = 0) -> bool:
    return (
        isinstance(value, int)
        and not isinstance(value, bool)
        and minimum <= value <= 2**63 - 1
    )


def _canonical_text(value: object, *, maximum: int = MAX_ROW_BYTES) -> str:
    content = canonical(value)
    if len(content) > maximum:
        raise WgsAnnotationError("bounded annotation row exceeds its capacity")
    return content.decode("ascii")


def _json_text(value: str, *, maximum: int = MAX_ROW_BYTES) -> dict:
    try:
        content = value.encode("ascii")
        if len(content) > maximum:
            raise ValueError
        parsed = json.loads(content)
        if type(parsed) is not dict or canonical(parsed) != content:
            raise ValueError
        return parsed
    except (UnicodeError, ValueError, TypeError, RecursionError, json.JSONDecodeError):
        raise WgsAnnotationError("invalid bounded annotation payload") from None


def _rows(cursor, size: int):
    while True:
        values = cursor.fetchmany(size)
        if not values:
            return
        yield from values


@dataclass(frozen=True, slots=True, repr=False)
class RetainedAllele:
    ordinal: int
    identity_id: str
    route: str
    source_allele_ids: tuple[str, ...]
    parent_source_manifest_sha256: str
    parent_working_set_manifest_sha256: str
    contig: str | None
    position: int | None
    reference: str | None
    alternate: str | None
    kind: str | None
    payload_json: str

    def __post_init__(self) -> None:
        sequence = self.route == "sequence"
        if (
            not _count(self.ordinal, 1)
            or self.route not in {"sequence", "native"}
            or not isinstance(self.source_allele_ids, tuple)
            or not self.source_allele_ids
            or len(self.source_allele_ids) > 65536
            or len(set(self.source_allele_ids)) != len(self.source_allele_ids)
            or any(_SOURCE_ALLELE.fullmatch(value) is None for value in self.source_allele_ids)
            or not _digest(self.parent_source_manifest_sha256)
            or not _digest(self.parent_working_set_manifest_sha256)
            or sequence != (_NORMALIZED.fullmatch(self.identity_id) is not None)
            or not sequence and _IDENTIFIER.fullmatch(self.identity_id) is None
        ):
            raise WgsAnnotationError("invalid retained-allele identity")
        geometry = (self.contig, self.position, self.reference, self.alternate, self.kind)
        if sequence:
            if (
                not isinstance(self.contig, str)
                or not self.contig
                or not _count(self.position, 1)
                or not isinstance(self.reference, str)
                or _BASES.fullmatch(self.reference) is None
                or not isinstance(self.alternate, str)
                or _BASES.fullmatch(self.alternate) is None
                or self.reference == self.alternate
                or self.kind not in {"snv", "mnv", "insertion", "deletion", "delins"}
            ):
                raise WgsAnnotationError("invalid retained sequence geometry")
        elif any(value is not None for value in geometry):
            if (
                any(value is None for value in geometry)
                or not isinstance(self.contig, str)
                or not self.contig
                or not _count(self.position, 1)
                or not isinstance(self.reference, str)
                or _BASES.fullmatch(self.reference) is None
                or not isinstance(self.alternate, str)
                or _BASES.fullmatch(self.alternate) is None
                or self.reference == self.alternate
                or self.kind not in {"snv", "mnv", "insertion", "deletion", "delins"}
            ):
                raise WgsAnnotationError("invalid retained native geometry")
        payload = _json_text(self.payload_json)
        if (
            payload.get("schema") != RETAINED_ALLELE_SCHEMA
            or payload.get("version") != 1
            or payload.get("classification") != "controlled"
            or payload.get("llm_policy") != "no-llm"
            or payload.get("identity_id") != self.identity_id
            or payload.get("route") != self.route
        ):
            raise WgsAnnotationError("invalid retained-allele payload")


@dataclass(frozen=True, slots=True, repr=False)
class VepTransportRequest:
    request_ordinal: int
    retained_ordinal: int
    request_id: str
    identity_id: str
    source_allele_ids: tuple[str, ...]
    parent_working_set_manifest_sha256: str
    record: VepVcfRecord

    def __post_init__(self) -> None:
        expected = "vep-" + hashlib.sha256(canonical({
            "identity_id": self.identity_id,
            "working_set_manifest_sha256": self.parent_working_set_manifest_sha256,
        })).hexdigest()[:24]
        if (
            not _count(self.request_ordinal, 1)
            or not _count(self.retained_ordinal, 1)
            or _IDENTIFIER.fullmatch(self.identity_id) is None
            or not isinstance(self.source_allele_ids, tuple)
            or not self.source_allele_ids
            or len(self.source_allele_ids) > 65536
            or len(set(self.source_allele_ids)) != len(self.source_allele_ids)
            or any(_SOURCE_ALLELE.fullmatch(value) is None for value in self.source_allele_ids)
            or not _digest(self.parent_working_set_manifest_sha256)
            or not isinstance(self.record, VepVcfRecord)
            or self.request_id != expected
            or self.record.record_id != self.request_id
        ):
            raise WgsAnnotationError("invalid VEP transport request")

    def as_bytes(self) -> bytes:
        content = canonical({
            "schema": VEP_REQUEST_SCHEMA,
            "version": 1,
            "classification": "controlled",
            "llm_policy": "no-llm",
            "data": {
                "request_ordinal": self.request_ordinal,
                "retained_ordinal": self.retained_ordinal,
                "request_id": self.request_id,
                "identity_id": self.identity_id,
                "source_allele_ids": self.source_allele_ids,
                "parent_working_set_manifest_sha256":
                    self.parent_working_set_manifest_sha256,
                "record": asdict(self.record),
            },
        })
        if len(content) > MAX_ROW_BYTES:
            raise WgsAnnotationError("VEP request exceeds its bounded row size")
        return content


def _load_request(content: bytes) -> VepTransportRequest:
    try:
        if len(content) > MAX_ROW_BYTES:
            raise ValueError
        value = json.loads(content)
        if canonical(value) != content or set(value) != {
            "schema", "version", "classification", "llm_policy", "data",
        }:
            raise ValueError
        if (
            value["schema"] != VEP_REQUEST_SCHEMA
            or value["version"] != 1
            or value["classification"] != "controlled"
            or value["llm_policy"] != "no-llm"
            or set(value["data"]) != {
                "request_ordinal", "retained_ordinal", "request_id", "identity_id",
                "source_allele_ids", "parent_working_set_manifest_sha256", "record",
            }
            or set(value["data"]["record"]) != {
                "record_id", "contig", "position", "reference", "alternates",
            }
        ):
            raise ValueError
        data = value["data"]
        record = data.pop("record")
        data["source_allele_ids"] = tuple(data["source_allele_ids"])
        record["alternates"] = tuple(record["alternates"])
        return VepTransportRequest(**data, record=VepVcfRecord(**record))
    except (KeyError, TypeError, ValueError, json.JSONDecodeError):
        raise WgsAnnotationError("invalid VEP request artifact") from None


@dataclass(frozen=True, slots=True, repr=False)
class VepTransportPartition:
    ordinal: int
    first_request_ordinal: int
    last_request_ordinal: int
    requests: int
    transport: Artifact
    expected_requests: Artifact
    transport_bytes: int
    request_bytes: int
    previous_partition_sha256: str | None

    def __post_init__(self) -> None:
        if (
            not _count(self.ordinal, 1)
            or self.ordinal > MAX_SHARDS
            or not _count(self.first_request_ordinal, 1)
            or not _count(self.last_request_ordinal, self.first_request_ordinal)
            or not _count(self.requests, 1)
            or self.last_request_ordinal - self.first_request_ordinal + 1 != self.requests
            or self.requests > 25_000
            or not isinstance(self.transport, Artifact)
            or not isinstance(self.expected_requests, Artifact)
            or self.transport.name != f"transport-{self.ordinal:06d}.vcf"
            or self.expected_requests.name != f"requests-{self.ordinal:06d}.jsonl"
            or self.transport.size_bytes != self.transport_bytes
            or self.expected_requests.size_bytes != self.request_bytes
            or not _count(self.transport_bytes, 1)
            or not _count(self.request_bytes, 1)
            or self.transport_bytes > MAX_SHARD_BYTES
            or self.request_bytes > MAX_SHARD_BYTES
            or (self.previous_partition_sha256 is not None
                and not _digest(self.previous_partition_sha256))
        ):
            raise WgsAnnotationError("invalid VEP transport partition")

    @property
    def sha256(self) -> str:
        return hashlib.sha256(canonical({
            "schema": VEP_PARTITION_SCHEMA,
            "version": 1,
            "classification": "controlled",
            "llm_policy": "no-llm",
            "data": asdict(self),
        })).hexdigest()


@dataclass(frozen=True, slots=True, repr=False)
class AnnotationIntakeBinding:
    manifest_sha256: str
    retained_identities: int
    retained_occurrences: int

    def __post_init__(self):
        if (not _digest(self.manifest_sha256) or not _count(self.retained_identities)
                or not _count(self.retained_occurrences, self.retained_identities)):
            raise WgsAnnotationError("annotation intake binding is invalid")


@dataclass(frozen=True, slots=True, repr=False)
class VepTransportPlan:
    parent_source_manifest_sha256: str
    parent_working_set_manifest_sha256: str
    execution_limits_sha256: str
    parents: tuple[ParentIdentity, ...]
    total_retained_identities: int
    sequence_requests: int
    native_identities: int
    partitions: tuple[VepTransportPartition, ...]
    intake_binding: AnnotationIntakeBinding | None = None

    def __post_init__(self) -> None:
        if (
            not _digest(self.parent_source_manifest_sha256)
            or not _digest(self.parent_working_set_manifest_sha256)
            or not _digest(self.execution_limits_sha256)
            or not _count(self.total_retained_identities)
            or not _count(self.sequence_requests)
            or not _count(self.native_identities)
            or self.sequence_requests + self.native_identities
            != self.total_retained_identities
            or not isinstance(self.parents, tuple)
            or not isinstance(self.partitions, tuple)
            or len(self.partitions) > MAX_SHARDS
        ):
            raise WgsAnnotationError("invalid VEP transport plan")
        expected_parents = {
            "source_manifest": (
                SOURCE_MANIFEST_SCHEMA, self.parent_source_manifest_sha256,
            ),
            "working_set_manifest": (
                WORKING_SET_MANIFEST_SCHEMA, self.parent_working_set_manifest_sha256,
            ),
            "execution_limits": (
                "genosapiens.bounded-wgs-execution", self.execution_limits_sha256,
            ),
        }
        if self.intake_binding is not None:
            if (type(self.intake_binding) is not AnnotationIntakeBinding
                    or self.total_retained_identities != self.intake_binding.retained_identities):
                raise WgsAnnotationError("VEP merged intake count is invalid")
            expected_parents["read_intake"] = (
                "genosapiens.wgs-merged-read-intake", self.intake_binding.manifest_sha256)
        if (
            len(self.parents) != len(expected_parents)
            or {parent.role for parent in self.parents} != set(expected_parents)
            or any(
                parent.parent_version != (
                    BOUNDED_WGS_CONFIG_VERSION
                    if parent.role == "execution_limits" else
                    2 if parent.role == "working_set_manifest" else 1
                )
                or (parent.parent_schema, parent.sha256) != expected_parents[parent.role]
                for parent in self.parents
            )
        ):
            raise WgsAnnotationError("invalid VEP transport parents")
        next_request = 1
        previous = None
        for ordinal, partition in enumerate(self.partitions, start=1):
            if (
                not isinstance(partition, VepTransportPartition)
                or partition.ordinal != ordinal
                or partition.first_request_ordinal != next_request
                or partition.previous_partition_sha256 != previous
            ):
                raise WgsAnnotationError("invalid VEP partition chain")
            next_request = partition.last_request_ordinal + 1
            previous = partition.sha256
        if next_request - 1 != self.sequence_requests:
            raise WgsAnnotationError("invalid VEP request accounting")
        if not self.partitions and self.sequence_requests:
            raise WgsAnnotationError("missing VEP transport partitions")

    def as_bytes(self) -> bytes:
        content = canonical({
            "schema": VEP_PLAN_SCHEMA,
            "version": 1,
            "classification": "controlled",
            "llm_policy": "no-llm",
            "data": asdict(self),
        })
        if len(content) > MAX_MANIFEST_BYTES:
            raise WgsAnnotationError("VEP transport plan exceeds manifest capacity")
        return content

    @property
    def sha256(self) -> str:
        return hashlib.sha256(self.as_bytes()).hexdigest()


@dataclass(frozen=True, slots=True, repr=False)
class VepTransportPublication:
    plan: VepTransportPlan
    manifest_artifact: Artifact


@dataclass(frozen=True, slots=True, repr=False)
class BoundedChildResult:
    exit_code: int
    stdout_sha256: str
    stdout_bytes: int
    stderr_sha256: str
    stderr_bytes: int
    total_output_bytes: int

    def __post_init__(self) -> None:
        if (
            self.exit_code != 0
            or not _digest(self.stdout_sha256)
            or not _digest(self.stderr_sha256)
            or not _count(self.stdout_bytes)
            or not _count(self.stderr_bytes)
            or self.total_output_bytes != self.stdout_bytes + self.stderr_bytes
        ):
            raise WgsAnnotationError("invalid bounded child result")


@dataclass(frozen=True, slots=True, repr=False)
class VerifiedVepRecord:
    request: VepTransportRequest
    parsed: ParsedVep
    raw_line_sha256: str
    parent_output_sha256: str

    def __post_init__(self) -> None:
        key = (self.request.request_id, 1)
        if (
            not isinstance(self.request, VepTransportRequest)
            or not isinstance(self.parsed, ParsedVep)
            or set(self.parsed.transcripts) != {key}
            or set(self.parsed.non_transcripts) != {key}
            or not RICH_FIELDS.issubset(self.parsed.fields)
            or not _digest(self.raw_line_sha256)
            or not _digest(self.parent_output_sha256)
        ):
            raise WgsAnnotationError("invalid verified VEP record")


@dataclass(frozen=True, slots=True, repr=False)
class UnassessedVepRecord:
    request: VepTransportRequest
    parent_output_sha256: str
    reference_fai_sha256: str
    reason: str = UNASSESSED_REASON

    @property
    def raw_line_sha256(self):
        return None

    def __post_init__(self):
        if (type(self.request) is not VepTransportRequest
                or not _digest(self.parent_output_sha256)
                or not _digest(self.reference_fai_sha256) or self.reason not in REASON_POLICIES):
            raise WgsAnnotationError('unassessed annotation outcome is invalid')

    def as_bytes(self):
        return canonical(dict(schema='genosapiens.wgs-unassessed-vep-record', version=1,
            policy=REASON_POLICIES[self.reason], request=json.loads(self.request.as_bytes()),
            parent_output_sha256=self.parent_output_sha256,
            reference_fai_sha256=self.reference_fai_sha256, reason=self.reason))


@dataclass(frozen=True, slots=True, repr=False)
class EncodedVerifiedVepRecord:
    """Validated writer-owned bytes; saved-store readers still decode strictly.

    Only the reconciliation owner creates these after scientific validation.
    The spool read authenticates both bytes and publication metadata before
    handing them to the coordinator's transactional shard writer.
    """

    request_id: str
    request_ordinal: int
    canonical_request: bytes
    parsed_record: bytes
    parsed_record_sha256: str
    raw_line_sha256: str | None
    parent_output_sha256: str

    def __post_init__(self):
        if (type(self.request_id) is not str or not _IDENTIFIER.fullmatch(self.request_id)
                or not _count(self.request_ordinal, 1)
                or type(self.canonical_request) is not bytes
                or not 0 < len(self.canonical_request) <= MAX_ROW_BYTES
                or type(self.parsed_record) is not bytes
                or not 0 < len(self.parsed_record) <= MAX_ROW_BYTES
                or not self.parsed_record.endswith(b"\n")
                or not all(_digest(value) for value in (
                    self.parsed_record_sha256, self.parent_output_sha256))
                or (self.raw_line_sha256 is not None if type(self) is EncodedUnassessedVepRecord
                    else not _digest(self.raw_line_sha256))):
            raise WgsAnnotationError("encoded VEP record is invalid")

    def metadata_sha256(self):
        return hashlib.sha256(canonical((self.request_id, self.request_ordinal,
            hashlib.sha256(self.canonical_request).hexdigest(), self.parsed_record_sha256,
            self.raw_line_sha256, self.parent_output_sha256))).hexdigest()


class EncodedUnassessedVepRecord(EncodedVerifiedVepRecord):
    """Explicit unassessed outcome with no fabricated VEP response digest."""
    __slots__ = ()


@dataclass(frozen=True, slots=True, repr=False)
class ParsedPartitionManifest:
    partition_ordinal: int
    partition_sha256: str
    parent_output_sha256: str
    records: int
    shards: tuple[Artifact, ...]

    def __post_init__(self) -> None:
        if (not _count(self.partition_ordinal, 1) or not _digest(self.partition_sha256)
                or not _digest(self.parent_output_sha256) or not _count(self.records, 1)
                or not isinstance(self.shards, tuple) or not self.shards
                or any(not isinstance(row, Artifact) for row in self.shards)):
            raise WgsAnnotationError("invalid parsed annotation partition")

    def as_bytes(self) -> bytes:
        return canonical({"schema": PARSED_PARTITION_SCHEMA, "version": 2,
                          "classification": "controlled", "llm_policy": "no-llm",
                          "data": asdict(self)})


@dataclass(frozen=True, slots=True, repr=False)
class AnnotationPartitionClosure:
    partition_ordinal: int
    partition_sha256: str
    parent_plan_sha256: str
    resource_session_sha256: str
    bundle_manifest_sha256: str
    execution_limits_sha256: str
    output: Artifact
    stderr: Artifact | None
    child: BoundedChildResult
    verified_responses: int
    parsed_manifest: Artifact
    unassessed_requests: int = 0
    reference_coverage: Artifact | None = None
    unassessed_alias_requests: int = 0

    @property
    def accounted_requests(self):
        return self.verified_responses + self.unassessed_requests

    @property
    def counts(self):
        return {'responses': self.verified_responses, **(
            {'unassessed': self.unassessed_requests} if self.unassessed_requests else {}), **(
            {'unassessed_alias': self.unassessed_alias_requests} if self.unassessed_alias_requests else {})}

    def __post_init__(self) -> None:
        if (
            not _count(self.partition_ordinal, 1)
            or not _digest(self.partition_sha256)
            or not _digest(self.parent_plan_sha256)
            or not _digest(self.resource_session_sha256)
            or not _digest(self.bundle_manifest_sha256)
            or not _digest(self.execution_limits_sha256)
            or not isinstance(self.output, Artifact)
            or self.stderr is not None and not isinstance(self.stderr, Artifact)
            or not isinstance(self.child, BoundedChildResult)
            or self.output.sha256 != self.child.stdout_sha256
            or self.output.size_bytes != self.child.stdout_bytes
            or self.stderr is not None and (
                self.stderr.sha256 != self.child.stderr_sha256
                or self.stderr.size_bytes != self.child.stderr_bytes
            )
            or self.stderr is None and self.child.stderr_bytes != 0
            or not _count(self.verified_responses)
            or not _count(self.unassessed_requests)
            or not _count(self.unassessed_alias_requests)
            or self.unassessed_alias_requests > self.unassessed_requests
            or not _count(self.accounted_requests, 1)
            or (self.reference_coverage is not None if not self.unassessed_requests
                else type(self.reference_coverage) is not Artifact)
            or not isinstance(self.parsed_manifest, Artifact)
        ):
            raise WgsAnnotationError("invalid annotation partition closure")

    def as_bytes(self) -> bytes:
        data = asdict(self)
        if not self.unassessed_alias_requests:
            data.pop('unassessed_alias_requests')
        if not self.unassessed_requests:
            data.pop('unassessed_requests')
            data.pop('reference_coverage')
        return canonical({
            "schema": ANNOTATION_CLOSURE_SCHEMA,
            "version": 4 if self.unassessed_alias_requests else 3 if self.unassessed_requests else 2,
            "classification": "controlled",
            "llm_policy": "no-llm",
            "data": data,
        })


@dataclass(frozen=True, slots=True, repr=False)
class AnnotationStageManifest:
    parent_plan_sha256: str
    resource_session_sha256: str
    bundle_manifest_sha256: str
    execution_limits_sha256: str
    parents: tuple[ParentIdentity, ...]
    sequence_partitions: int
    sequence_requests: int
    partition_closures: tuple[Artifact, ...]
    parsed_manifests: tuple[Artifact, ...]
    state_artifact: Artifact
    prefix_source: Artifact | None = None

    def __post_init__(self) -> None:
        if (
            not _digest(self.parent_plan_sha256)
            or not _digest(self.resource_session_sha256)
            or not _digest(self.bundle_manifest_sha256)
            or not _digest(self.execution_limits_sha256)
            or not isinstance(self.parents, tuple)
            or not _count(self.sequence_partitions)
            or not _count(self.sequence_requests)
            or not isinstance(self.partition_closures, tuple)
            or len(self.partition_closures) != self.sequence_partitions
            or not isinstance(self.parsed_manifests, tuple)
            or len(self.parsed_manifests) != self.sequence_partitions
            or any(not isinstance(value, Artifact) for value in self.partition_closures)
            or any(not isinstance(value, Artifact) for value in self.parsed_manifests)
            or not isinstance(self.state_artifact, Artifact)
            or self.state_artifact.name != "annotation.sqlite"
            or self.prefix_source is not None and (type(self.prefix_source) is not Artifact
                or self.prefix_source.name != "annotation-prefix.json")
        ):
            raise WgsAnnotationError("invalid annotation stage manifest")

    def as_bytes(self) -> bytes:
        data = asdict(self)
        if self.prefix_source is None:
            data.pop("prefix_source")
        content = canonical({
            "schema": ANNOTATION_STAGE_SCHEMA,
            "version": 2 if self.prefix_source is None else 3,
            "classification": "controlled",
            "llm_policy": "no-llm",
            "data": data,
        })
        if len(content) > MAX_MANIFEST_BYTES:
            raise WgsAnnotationError("annotation stage manifest exceeds capacity")
        return content

    @property
    def sha256(self) -> str:
        return hashlib.sha256(self.as_bytes()).hexdigest()


@dataclass(frozen=True, slots=True, repr=False)
class AnnotationStageResult:
    manifest: AnnotationStageManifest
    manifest_artifact: Artifact
    journal_head: JournalEntry
    reused_partitions: int
    prefix: object | None = None


def _result_prefix(result, root, plan, config):
    from genosapiens.wgs_annotation_prefix import (ResolvedAnnotationPrefix,
        _load_producer_annotation, _resolve_producer_annotation)
    artifact = result.manifest.prefix_source
    if artifact is None:
        if result.prefix is not None:
            raise WgsAnnotationError("unexpected annotation prefix")
        return None
    content = read_artifact(root, artifact)
    if result.prefix is None:
        descriptor = _load_producer_annotation(content, expected_sha256=artifact.sha256)
        resolved = _resolve_producer_annotation(descriptor, plan=plan, execution_config=config)
        # An immutable result owns this checked metadata/index cache. No process
        # global cache and no repeated aggregate scan for every scientific reader.
        object.__setattr__(result, "prefix", resolved)
    prefix = result.prefix
    if (type(prefix) is not ResolvedAnnotationPrefix or prefix.descriptor.as_bytes() != content
            or prefix.descriptor.plan_sha256 != plan.sha256):
        raise WgsAnnotationError("annotation prefix binding changed")
    return prefix


def _partition_source(prefix, partition, root, resource_session_sha256):
    if prefix is not None and partition.ordinal <= prefix.count:
        return prefix.partition_source(partition)
    return root, resource_session_sha256, None


@dataclass(frozen=True, slots=True, repr=False)
class BoundedAnnotationSource:
    """One completed immutable annotation source for authorized scientific replay."""

    plan_root: Path
    plan: VepTransportPlan
    output_root: Path
    result: AnnotationStageResult
    native_index: Artifact | None = None

    def __post_init__(self):
        if (not self.plan_root.is_absolute() or not self.output_root.is_absolute()
                or type(self.plan) is not VepTransportPlan or type(self.result) is not AnnotationStageResult
                or self.result.manifest.parent_plan_sha256 != self.plan.sha256
                or self.result.manifest_artifact.sha256 != self.result.manifest.sha256
                or self.result.journal_head.event_type != "stage_completed"
                or len(self.result.journal_head.artifacts) != 1
                or self.result.journal_head.artifacts[0].name
                != "annotation-stage-closure.json"):
            raise WgsAnnotationError("fixed annotation source binding is invalid")
        if self.native_index is not None and (
                type(self.native_index) is not Artifact or self.native_index.name != "native-index.jsonl"):
            raise WgsAnnotationError("fixed native annotation source binding is invalid")

    def as_bytes(self):
        parsed_binding = hashlib.sha256(canonical(tuple(asdict(row) for row in
            self.result.manifest.parsed_manifests))).hexdigest()
        return canonical({"schema": "genosapiens.wgs-fixed-annotation-source", "version": 3,
            "classification": "controlled", "llm_policy": "no-llm",
            "plan_sha256": self.plan.sha256, "stage_manifest_sha256": self.result.manifest.sha256,
            "parsed_manifests_sha256": parsed_binding,
            "journal_head_sha256": self.result.journal_head.sha256,
            "native_index": None if self.native_index is None else asdict(self.native_index)})

    @property
    def sha256(self):
        return hashlib.sha256(self.as_bytes()).hexdigest()

    def verify(self, *, expected_plan, source_manifest, working_manifest, resource_session_sha256,
               bundle_manifest_sha256, executable_stack_sha256, runtime_image_digest,
               execution_config, intake_binding=None, prepared_index=None,
               prepared_state_root=None):
        prefix = _result_prefix(self.result, self.output_root, self.plan, execution_config)
        expected_parents = _annotation_parents(expected_plan,
            execution_limits_sha256=hashlib.sha256(bounded_wgs_execution_config_as_json(execution_config)).hexdigest(),
            execution_limits_version=execution_config.version,
            resource_session_sha256=resource_session_sha256, bundle_manifest_sha256=bundle_manifest_sha256,
            executable_stack_sha256=executable_stack_sha256, runtime_image_digest=runtime_image_digest,
            prefix_sha256=None if self.result.manifest.prefix_source is None else self.result.manifest.prefix_source.sha256,
            prefix_version=1 if prefix is None else prefix.descriptor.version)
        if (self.plan.as_bytes() != expected_plan.as_bytes()
                or self.result.manifest.parents != expected_parents
                or read_artifact(self.output_root, self.result.manifest_artifact) != self.result.manifest.as_bytes()):
            raise WgsAnnotationError("fixed annotation source differs from current inputs")
        _verify_state_descriptor(self.output_root, self.result.manifest.state_artifact)
        from genosapiens.wgs_ledger import read_journal_head
        terminal = read_journal_head(self.output_root, stage="annotation", parent_sha256=self.plan.sha256,
                                    expected_head_sha256=self.result.journal_head.sha256)
        if terminal != self.result.journal_head or terminal.parents != expected_parents:
            raise WgsAnnotationError("fixed annotation source journal changed")
        seal = json.loads(read_artifact(
            self.output_root, self.result.journal_head.artifacts[0],
        ))
        if (seal.get("schema") != "genosapiens.wgs-annotation-stage-closure"
                or seal.get("version") != 1
                or seal.get("plan_sha256") != self.plan.sha256
                or seal.get("partition_closures") != [
                    item.sha256 for item in self.result.manifest.partition_closures
                ]):
            raise WgsAnnotationError("fixed annotation stage closure changed")
        load_completed_vep_plan(output_root=self.plan_root, expected_manifest_sha256=self.plan.sha256,
            source_manifest=source_manifest, working_manifest=working_manifest,
            execution_config=execution_config, intake_binding=intake_binding,
            prepared_index=prepared_index, prepared_state_root=prepared_state_root)
        prefix = _result_prefix(self.result, self.output_root, self.plan, execution_config)
        for partition, artifact, parsed in zip(self.plan.partitions,
                self.result.manifest.partition_closures,
                self.result.manifest.parsed_manifests, strict=True):
            origin, session, counts = _partition_source(prefix, partition, self.output_root, resource_session_sha256)
            closure = read_completed_partition_closure(origin, artifact, partition,
                plan=self.plan, plan_root=self.plan_root,
                resource_session_sha256=session,
                bundle_manifest_sha256=bundle_manifest_sha256,
                execution_config=execution_config)
            if closure.parsed_manifest != parsed:
                raise WgsAnnotationError("fixed annotation parsed binding is invalid")


def load_fixed_annotation_source(content, *, expected_sha256, plan_root, output_root):
    """Decode current saved annotation metadata without launching annotation."""
    from genosapiens.evidence_model import _typed
    from genosapiens.wgs_shards import read_bounded_metadata
    try:
        if (type(content) is not bytes or len(content) > MAX_MANIFEST_BYTES
                or not _digest(expected_sha256) or hashlib.sha256(content).hexdigest() != expected_sha256):
            raise ValueError
        data = json.loads(content)
        if (set(data) - {"native_index"} != {"schema", "version", "classification", "llm_policy",
                "plan_sha256", "stage_manifest_sha256", "parsed_manifests_sha256",
                "journal_head_sha256"}
                or data["schema"] != "genosapiens.wgs-fixed-annotation-source"
                or type(data["version"]) is not int or data["version"] != 3
                or data["classification"] != "controlled" or data["llm_policy"] != "no-llm"):
            raise ValueError
        plan_root, output_root = Path(plan_root), Path(output_root)
        plan_content = read_bounded_metadata(plan_root / "manifest.json")
        plan = _typed(json.loads(plan_content)["data"], VepTransportPlan, decode=True)
        if plan.as_bytes() != plan_content or plan.sha256 != data["plan_sha256"]:
            raise ValueError
        stage_content = read_bounded_metadata(output_root / "annotation-manifest.json")
        stage = _load_stage_manifest(stage_content)
        if stage.sha256 != data["stage_manifest_sha256"]:
            raise ValueError
        if hashlib.sha256(canonical(tuple(asdict(row) for row in
                stage.parsed_manifests))).hexdigest() != data["parsed_manifests_sha256"]:
            raise ValueError
        artifact = Artifact("annotation-manifest.json", stage.sha256, len(stage_content))
        from genosapiens.wgs_ledger import read_journal_head
        terminal = read_journal_head(output_root, stage="annotation", parent_sha256=plan.sha256,
                                    expected_head_sha256=data["journal_head_sha256"])
        result = BoundedAnnotationSource(plan_root, plan, output_root,
            AnnotationStageResult(stage, artifact, terminal, 0),
            None if data["native_index"] is None else Artifact(**data["native_index"]))
        if result.as_bytes() != content:
            raise ValueError
        return result
    except (ValueError, TypeError, KeyError, IndexError, OSError, LedgerError):
        raise WgsAnnotationError("saved fixed annotation source is invalid") from None


def _working_rows(root: Path, manifest: WorkingSetManifest):
    manifest_content = manifest.as_bytes()
    expected = Artifact(
        "manifest.json", hashlib.sha256(manifest_content).hexdigest(),
        len(manifest_content),
    )
    if read_artifact(root, expected) != manifest_content:
        raise WgsAnnotationError("working-set manifest artifact is invalid")
    summary_content = manifest.summary.as_bytes()
    summary = Artifact(
        "summary.json", hashlib.sha256(summary_content).hexdigest(),
        len(summary_content),
    )
    if read_artifact(root, summary) != summary_content:
        raise WgsAnnotationError("working-set summary artifact is invalid")
    for shard in manifest.shards:
        yield from iter_working_set_shard(root, shard)


def _working_binds_source(
    manifest: WorkingSetManifest, source_manifest: LedgerManifest,
) -> bool:
    parents = {parent.role: parent for parent in manifest.parents}
    parent = parents.get("source_manifest")
    filtered = parents.get("filter_manifest")
    return (
        manifest.version == 2
        and manifest.summary.version == 2
        and set(parents) == {"source_manifest", "filter_manifest"}
        and parent is not None
        and parent.parent_schema == SOURCE_MANIFEST_SCHEMA
        and parent.parent_version == 1
        and parent.sha256 == source_manifest.sha256
        and filtered is not None
        and filtered.parent_schema == "genosapiens.wgs-filter-manifest"
        and filtered.parent_version == 2
        and filtered.sha256 == manifest.parent_filter_manifest_sha256
    )


def _open_index(path: Path, config: BoundedWgsExecutionConfig) -> sqlite3.Connection:
    descriptor = -1
    try:
        descriptor = os.open(
            path, os.O_WRONLY | os.O_CREAT | os.O_EXCL | os.O_NOFOLLOW, 0o600,
        )
        os.fchmod(descriptor, 0o600)
        os.close(descriptor)
        descriptor = -1
        connection = sqlite3.connect(path)
        connection.enable_load_extension(False)
        connection.execute("PRAGMA foreign_keys=ON")
        connection.execute(f"PRAGMA cache_size=-{max(1, config.sqlite.cache_bytes // 1024)}")
        connection.execute("PRAGMA temp_store=FILE")
        connection.execute("PRAGMA mmap_size=0")
        connection.execute("PRAGMA trusted_schema=OFF")
        return connection
    except (OSError, sqlite3.Error):
        if descriptor >= 0:
            os.close(descriptor)
        raise WgsAnnotationError("cannot create retained-allele index") from None


def _projection_binding(prepared_index: PreparedIndex, routing: SourceRoutingContext,
                        config: BoundedWgsExecutionConfig) -> str:
    return hashlib.sha256(canonical({
        "schema": "genosapiens.wgs-retained-projection", "version": 1,
        "prepared_index_sha256": prepared_index.sha256,
        "source_routing_sha256": routing.contract.sha256,
        "execution_limits_sha256": hashlib.sha256(
            bounded_wgs_execution_config_as_json(config)).hexdigest(),
    })).hexdigest()


def _require_prepared_projection_parent(prepared_index: PreparedIndex,
                                        state_root: Path, source_manifest: LedgerManifest,
                                        working_manifest: WorkingSetManifest) -> None:
    if (not isinstance(prepared_index, PreparedIndex)
            or prepared_index.version != PREPARED_INDEX_VERSION
            or not _working_binds_source(working_manifest, source_manifest)
            or not state_root.is_absolute()
            or str(state_root) != prepared_index.preparation_state_root
            or prepared_index.source_manifest_sha256 != source_manifest.sha256
            or prepared_index.filter_manifest_sha256
            != working_manifest.parent_filter_manifest_sha256
            or prepared_index.working_manifest_sha256 != working_manifest.sha256):
        raise WgsAnnotationError("prepared projection parent is invalid")
    try:
        directory = state_root.lstat()
        if (not stat.S_ISDIR(directory.st_mode) or directory.st_uid != os.getuid()
                or stat.S_IMODE(directory.st_mode) != 0o700):
            raise ValueError
        for descriptor in prepared_index.states:
            if descriptor.stage not in {"source", "prefilter"}:
                continue
            path = prepared_index.state_path(descriptor)
            metadata = path.lstat()
            if (not stat.S_ISREG(metadata.st_mode) or metadata.st_uid != os.getuid()
                    or stat.S_IMODE(metadata.st_mode) != 0o600 or metadata.st_nlink != 1
                    or metadata.st_size != descriptor.size_bytes):
                raise ValueError
            with closing(sqlite3.connect(path.as_uri() + "?mode=ro&immutable=1", uri=True)) as db:
                if (db.execute("PRAGMA user_version").fetchone()[0] != descriptor.user_version
                        or hashlib.sha256(db.execute(
                            "SELECT canonical_bytes FROM run_binding WHERE singleton=1"
                        ).fetchone()[0]).hexdigest() != descriptor.run_binding_sha256):
                    raise ValueError
    except (OSError, ValueError, TypeError, sqlite3.Error):
        raise WgsAnnotationError("prepared projection state is missing or stale") from None


def open_retained_projection(path: Path, *, prepared_index: PreparedIndex,
                             routing: SourceRoutingContext,
                             execution_config: BoundedWgsExecutionConfig,
                             expected_identities: int,
                             expected_occurrences: int) -> sqlite3.Connection:
    """Open one completed case-local projection without replaying source shards."""
    connection = None
    try:
        metadata = path.lstat()
        if (not path.is_absolute() or not stat.S_ISREG(metadata.st_mode)
                or metadata.st_uid != os.getuid() or stat.S_IMODE(metadata.st_mode) != 0o600
                or metadata.st_nlink != 1 or metadata.st_size <= 0):
            raise ValueError
        connection = sqlite3.connect(path.as_uri() + "?mode=ro&immutable=1", uri=True)
        connection.enable_load_extension(False)
        connection.execute("PRAGMA query_only=ON")
        connection.execute("PRAGMA trusted_schema=OFF")
        connection.execute(f"PRAGMA cache_size=-{max(1, execution_config.sqlite.cache_bytes // 1024)}")
        if (connection.execute("PRAGMA user_version").fetchone()[0] != 1
                or connection.execute("SELECT binding_sha256,identities,occurrences,complete "
                    "FROM projection_meta WHERE singleton=1").fetchone() != (
                        _projection_binding(prepared_index, routing, execution_config),
                        expected_identities, expected_occurrences, 1)):
            raise ValueError
        return connection
    except (OSError, ValueError, TypeError, sqlite3.Error):
        if connection is not None:
            connection.close()
        raise WgsAnnotationError("retained projection is missing or incomplete") from None


def _iter_retained_projection_rows(connection: sqlite3.Connection,
                                   source_manifest: LedgerManifest,
                                   working_manifest: WorkingSetManifest,
                                   config: BoundedWgsExecutionConfig) -> Iterator[RetainedAllele]:
    cursor = connection.execute("""
        SELECT i.ordinal,i.identity_id,i.route,f.contig,f.position,f.reference,
               f.alternate,f.kind,f.payload_json
        FROM identities i JOIN facts f USING(identity_id) ORDER BY i.ordinal
    """)
    previous = 0
    occurrences = 0
    for values in _rows(cursor, config.sqlite.cursor_batch_rows):
        ordinal, identity_id, route, contig, position, reference, alternate, kind, payload = values
        if ordinal != previous + 1:
            raise WgsAnnotationError("retained projection order is invalid")
        source_ids = tuple(value[0] for value in _rows(
            connection.execute("SELECT source_allele_id FROM occurrences WHERE identity_id=? "
                               "ORDER BY occurrence_ordinal", (identity_id,)),
            config.sqlite.cursor_batch_rows))
        if not source_ids:
            raise WgsAnnotationError("retained projection occurrence is missing")
        previous = ordinal
        occurrences += len(source_ids)
        yield RetainedAllele(ordinal, identity_id, route, source_ids,
            source_manifest.sha256, working_manifest.sha256,
            contig, position, reference, alternate, kind, payload)
    if (previous != working_manifest.summary.retained_unique_identities
            or occurrences != working_manifest.summary.retained_occurrences):
        raise WgsAnnotationError("retained projection coverage is incomplete")


def iter_retained_alleles(
    *,
    source_root: str | Path,
    source_manifest: LedgerManifest,
    working_root: str | Path,
    working_manifest: WorkingSetManifest,
    source_routing_context: SourceRoutingContext,
    scratch_path: str | Path,
    execution_config: BoundedWgsExecutionConfig,
    prepared_index: PreparedIndex | None = None,
    prepared_state_root: str | Path | None = None,
) -> Iterator[RetainedAllele]:
    """Yield complete retained identities after a disk-backed source-ledger join."""

    source_root, working_root, scratch_path = (
        Path(source_root), Path(working_root), Path(scratch_path)
    )
    if (
        not isinstance(source_manifest, LedgerManifest)
        or not isinstance(working_manifest, WorkingSetManifest)
        or not isinstance(source_routing_context, SourceRoutingContext)
        or source_routing_context.contract.parent_source_manifest_sha256
        != source_manifest.sha256
        or source_routing_context.contract.reference_manifest_sha256
        != source_manifest.source.reference_manifest_sha256
        or not isinstance(execution_config, BoundedWgsExecutionConfig)
        or not source_root.is_absolute()
        or not working_root.is_absolute()
        or not scratch_path.is_absolute()
        or not _working_binds_source(working_manifest, source_manifest)
    ):
        raise WgsAnnotationError("retained-allele cursor binding is invalid")
    if prepared_index is not None:
        if prepared_state_root is None:
            raise WgsAnnotationError("prepared projection state root is missing")
        _require_prepared_projection_parent(prepared_index, Path(prepared_state_root),
                                            source_manifest, working_manifest)
        if scratch_path.exists():
            with closing(open_retained_projection(scratch_path,
                    prepared_index=prepared_index, routing=source_routing_context,
                    execution_config=execution_config,
                    expected_identities=working_manifest.summary.retained_unique_identities,
                    expected_occurrences=working_manifest.summary.retained_occurrences)) as existing:
                yield from _iter_retained_projection_rows(existing, source_manifest,
                                                          working_manifest, execution_config)
            return
    elif prepared_state_root is not None:
        raise WgsAnnotationError("unbound prepared projection state root")
    connection = _open_index(scratch_path, execution_config)
    sealed = False
    try:
        connection.executescript("""
            CREATE TABLE identities(
                identity_id TEXT PRIMARY KEY,
                ordinal INTEGER NOT NULL UNIQUE,
                route TEXT NOT NULL,
                parent_filter_sha256 TEXT NOT NULL
            );
            CREATE TABLE occurrences(
                source_allele_id TEXT PRIMARY KEY,
                identity_id TEXT NOT NULL REFERENCES identities(identity_id),
                occurrence_ordinal INTEGER NOT NULL,
                seen INTEGER NOT NULL DEFAULT 0,
                record_number INTEGER,
                allele_number INTEGER,
                UNIQUE(identity_id,occurrence_ordinal)
            );
            CREATE INDEX occurrence_identity ON occurrences(identity_id,occurrence_ordinal);
            CREATE UNIQUE INDEX occurrence_source_order ON occurrences(record_number,allele_number);
            CREATE TABLE source_payloads(record_number INTEGER PRIMARY KEY,
                payload_json TEXT NOT NULL);
            CREATE TABLE facts(
                identity_id TEXT PRIMARY KEY REFERENCES identities(identity_id),
                contig TEXT,
                position INTEGER,
                reference TEXT,
                alternate TEXT,
                kind TEXT,
                payload_json TEXT NOT NULL
            );
            CREATE TABLE projection_meta(singleton INTEGER PRIMARY KEY CHECK(singleton=1),
                binding_sha256 TEXT NOT NULL, identities INTEGER NOT NULL,
                occurrences INTEGER NOT NULL, complete INTEGER NOT NULL CHECK(complete IN (0,1)));
            PRAGMA user_version=1;
        """)
        identity_count = occurrence_count = 0
        for row in _working_rows(working_root, working_manifest):
            if row.ordinal != identity_count + 1:
                raise WgsAnnotationError("working-set identity order is invalid")
            connection.execute(
                "INSERT INTO identities VALUES (?,?,?,?)",
                (
                    row.identity_id, row.ordinal, row.route,
                    row.parent_filter_manifest_sha256,
                ),
            )
            connection.executemany(
                "INSERT INTO occurrences(source_allele_id,identity_id,occurrence_ordinal) "
                "VALUES (?,?,?)",
                (
                    (source_id, row.identity_id, ordinal)
                    for ordinal, source_id in enumerate(row.source_allele_ids, start=1)
                ),
            )
            identity_count += 1
            occurrence_count += len(row.source_allele_ids)
            if identity_count % execution_config.sqlite.cursor_batch_rows == 0:
                connection.commit()
        if (
            identity_count != working_manifest.summary.retained_unique_identities
            or occurrence_count != working_manifest.summary.retained_occurrences
        ):
            raise WgsAnnotationError("working-set cursor accounting is invalid")
        connection.commit()

        for source_row in iter_ledger(source_root, source_manifest):
            try:
                payload = json.loads(source_row.payload_json)
                outcomes = payload["normalization"]["alleles"]
                if len(outcomes) != len(source_row.source_allele_ids):
                    raise ValueError
                selected_in_record = False
                for allele_number, (source_id, outcome) in enumerate(
                    zip(source_row.source_allele_ids, outcomes, strict=True), start=1,
                ):
                    selected = connection.execute(
                        "SELECT identity_id,route,occurrence_ordinal,seen FROM occurrences "
                        "JOIN identities USING(identity_id) WHERE source_allele_id=?",
                        (source_id,),
                    ).fetchone()
                    if selected is None:
                        continue
                    if not selected_in_record:
                        connection.execute("INSERT INTO source_payloads VALUES (?,?)",
                                           (source_row.record_number, source_row.payload_json))
                        selected_in_record = True
                    identity_id, route, occurrence_ordinal, seen = selected
                    if seen:
                        raise WgsAnnotationError("duplicate retained source occurrence")
                    facts = decode_source_occurrence(
                        source_row.payload_json,
                        source_allele_id=source_id,
                        source_manifest=source_manifest,
                        contract=source_routing_context.contract,
                        expected_record_number=source_row.record_number,
                        expected_source_record_id=source_row.source_record_id,
                        expected_source_allele_ids=source_row.source_allele_ids,
                    )
                    expected_route = classify_source_occurrence(
                        facts,
                        contract=source_routing_context.contract,
                        renderer=source_routing_context.renderer,
                    )
                    normalized = outcome.get("normalized_allele")
                    if route == "sequence":
                        if (
                            expected_route.preparation_disposition != "normalized"
                            or expected_route.normalized_allele_id != identity_id
                            or type(normalized) is not dict
                            or normalized.get("normalized_allele_id") != identity_id
                        ):
                            raise WgsAnnotationError("retained normalization mapping is invalid")
                        fact = (
                            normalized["contig"], normalized["position"],
                            normalized["reference"], normalized["alternate"],
                            normalized["kind"],
                        )
                        detail = {
                            "normalized_allele": normalized,
                            "first_source_record_number": source_row.record_number,
                            "first_source_allele_number": allele_number,
                        }
                    else:
                        if (
                            expected_route.preparation_disposition != "native_class"
                            or expected_route.normalized_allele_id is not None
                            or identity_id != source_id
                            or occurrence_ordinal != 1
                        ):
                            raise WgsAnnotationError("retained native mapping is invalid")
                        fact = (
                            (
                                normalized["contig"], normalized["position"],
                                normalized["reference"], normalized["alternate"],
                                normalized["kind"],
                            )
                            if type(normalized) is dict else
                            (None, None, None, None, None)
                        )
                        detail = {
                            "source_record_id": source_row.source_record_id,
                            "source_record_number": source_row.record_number,
                            "source_allele_number": allele_number,
                            "source_payload_sha256": hashlib.sha256(
                                source_row.payload_json.encode("ascii")
                            ).hexdigest(),
                        }
                    existing = connection.execute(
                        "SELECT contig,position,reference,alternate,kind FROM facts "
                        "WHERE identity_id=?",
                        (identity_id,),
                    ).fetchone()
                    if existing is None:
                        payload_json = _canonical_text({
                            "schema": RETAINED_ALLELE_SCHEMA,
                            "version": 1,
                            "classification": "controlled",
                            "llm_policy": "no-llm",
                            "identity_id": identity_id,
                            "route": route,
                            "detail": detail,
                        })
                        connection.execute(
                            "INSERT INTO facts VALUES (?,?,?,?,?,?,?)",
                            (identity_id, *fact, payload_json),
                        )
                    elif tuple(existing) != fact:
                        raise WgsAnnotationError("retained identity geometry is inconsistent")
                    updated = connection.execute(
                        "UPDATE occurrences SET seen=1,record_number=?,allele_number=? "
                        "WHERE source_allele_id=? AND seen=0",
                        (source_row.record_number, allele_number, source_id),
                    )
                    if updated.rowcount != 1:
                        raise WgsAnnotationError("retained occurrence replay is invalid")
            except (
                KeyError, TypeError, ValueError, json.JSONDecodeError,
                sqlite3.Error, WgsSourceRoutingError,
            ):
                raise WgsAnnotationError("source-ledger retained join is invalid") from None
        connection.commit()
        unseen = connection.execute(
            "SELECT count(*) FROM occurrences WHERE seen=0 "
            "OR record_number IS NULL OR allele_number IS NULL"
        ).fetchone()[0]
        facts = connection.execute("SELECT count(*) FROM facts").fetchone()[0]
        if unseen or facts != identity_count:
            raise WgsAnnotationError("retained source coverage is incomplete")
        if prepared_index is not None:
            connection.execute("INSERT INTO projection_meta VALUES (1,?,?,?,1)", (
                _projection_binding(prepared_index, source_routing_context, execution_config),
                identity_count, occurrence_count))
            connection.commit()
        sealed = True
        yield from _iter_retained_projection_rows(connection, source_manifest,
                                                 working_manifest, execution_config)
    except sqlite3.Error:
        raise WgsAnnotationError("retained-allele index failed") from None
    finally:
        connection.close()
        if prepared_index is None:
            try:
                scratch_path.unlink()
            except (FileNotFoundError, OSError):
                pass
        if not sealed:
            # The static exception raised above remains authoritative.
            pass


def _plan_parents(
    source_manifest: LedgerManifest,
    working_manifest: WorkingSetManifest,
    execution_limits_sha256: str,
    execution_limits_version: int,
    intake_binding: AnnotationIntakeBinding | None = None,
) -> tuple[ParentIdentity, ...]:
    return (
        ParentIdentity(
            "source_manifest", SOURCE_MANIFEST_SCHEMA, 1, source_manifest.sha256,
        ),
        ParentIdentity(
            "working_set_manifest", WORKING_SET_MANIFEST_SCHEMA, working_manifest.version,
            working_manifest.sha256,
        ),
        ParentIdentity(
            "execution_limits", "genosapiens.bounded-wgs-execution",
            execution_limits_version,
            execution_limits_sha256,
        ),
        *((ParentIdentity("read_intake", "genosapiens.wgs-merged-read-intake", 1,
                          intake_binding.manifest_sha256),) if intake_binding else ()),
    )


def _request(request_ordinal: int, allele: RetainedAllele) -> VepTransportRequest:
    request_id = "vep-" + hashlib.sha256(canonical({
        "identity_id": allele.identity_id,
        "working_set_manifest_sha256": allele.parent_working_set_manifest_sha256,
    })).hexdigest()[:24]
    return VepTransportRequest(
        request_ordinal,
        allele.ordinal,
        request_id,
        allele.identity_id,
        allele.source_allele_ids,
        allele.parent_working_set_manifest_sha256,
        VepVcfRecord(
            request_id, allele.contig, allele.position, allele.reference,
            (allele.alternate,),
        ),
    )


def load_completed_vep_plan(
    *, output_root: str | Path, source_manifest: LedgerManifest,
    working_manifest: WorkingSetManifest, execution_config: BoundedWgsExecutionConfig,
    expected_manifest_sha256: str | None = None,
    intake_binding: AnnotationIntakeBinding | None = None,
    prepared_index: PreparedIndex | None = None,
    prepared_state_root: str | Path | None = None,
) -> VepTransportPublication:
    """Read complete plan metadata without opening a writer or retained cursor.

    Imported inputs require an externally bound manifest digest. Exact-version
    resume uses the same canonical/parent checks. Transport bytes are verified
    by their existing consumers, without rebuilding or scanning a retained DB.
    """
    from genosapiens.wgs_shards import read_bounded_metadata
    root = Path(output_root)
    if (not root.is_absolute() or not isinstance(execution_config, BoundedWgsExecutionConfig)
            or not _working_binds_source(working_manifest, source_manifest)
            or expected_manifest_sha256 is not None and not _digest(expected_manifest_sha256)):
        raise WgsAnnotationError("completed VEP plan binding is invalid")
    if prepared_index is not None:
        if prepared_state_root is None:
            raise WgsAnnotationError("VEP prepared projection state is missing")
        _require_prepared_projection_parent(prepared_index, Path(prepared_state_root),
                                            source_manifest, working_manifest)
    elif prepared_state_root is not None:
        raise WgsAnnotationError("VEP prepared projection state is unbound")
    try:
        content = read_bounded_metadata(root / "manifest.json")
        digest = hashlib.sha256(content).hexdigest()
        limits_sha256 = hashlib.sha256(bounded_wgs_execution_config_as_json(execution_config)).hexdigest()
        value = json.loads(content)
        plan = _typed(value["data"], VepTransportPlan, decode=True)
        if (expected_manifest_sha256 is not None and digest != expected_manifest_sha256
                or plan.as_bytes() != content
                or plan.parent_source_manifest_sha256 != source_manifest.sha256
                or plan.parent_working_set_manifest_sha256 != working_manifest.sha256
                or plan.execution_limits_sha256 != limits_sha256
                or plan.parents != _plan_parents(source_manifest, working_manifest,
                    limits_sha256, execution_config.version, intake_binding)
                or plan.intake_binding != intake_binding):
            raise ValueError
        return VepTransportPublication(plan, Artifact("manifest.json", digest, len(content)))
    except (OSError, KeyError, TypeError, ValueError, RuntimeError):
        raise WgsAnnotationError("completed VEP plan binding is invalid") from None


def write_vep_partitions(
    retained: Iterable[RetainedAllele],
    *,
    source_manifest: LedgerManifest,
    working_manifest: WorkingSetManifest,
    output_root: str | Path,
    execution_config: BoundedWgsExecutionConfig,
    intake_binding: AnnotationIntakeBinding | None = None,
    prepared_index: PreparedIndex | None = None,
    prepared_state_root: str | Path | None = None,
    resume: bool = False,
) -> VepTransportPublication:
    """Write immutable, source-ordered VEP partitions under count and byte caps."""

    root = Path(output_root)
    if (
        not isinstance(source_manifest, LedgerManifest)
        or not isinstance(working_manifest, WorkingSetManifest)
        or not isinstance(execution_config, BoundedWgsExecutionConfig)
        or type(resume) is not bool
        or not root.is_absolute()
        or execution_config.annotation.max_partition_requests > 25_000
        or execution_config.annotation.max_partition_bytes > MAX_SHARD_BYTES
        or not _working_binds_source(working_manifest, source_manifest)
    ):
        raise WgsAnnotationError("VEP partition binding is invalid")
    if prepared_index is not None:
        if prepared_state_root is None:
            raise WgsAnnotationError("VEP prepared projection state is missing")
        _require_prepared_projection_parent(prepared_index, Path(prepared_state_root),
                                            source_manifest, working_manifest)
    elif prepared_state_root is not None:
        raise WgsAnnotationError("VEP prepared projection state is unbound")
    limits_content = bounded_wgs_execution_config_as_json(execution_config)
    limits_sha256 = hashlib.sha256(limits_content).hexdigest()
    state_path = root / "partition-state.sqlite"
    if resume:
        manifest_path = root / "manifest.json"
        if manifest_path.exists() or manifest_path.is_symlink():
            return load_completed_vep_plan(output_root=root, source_manifest=source_manifest,
                working_manifest=working_manifest, execution_config=execution_config,
                intake_binding=intake_binding, prepared_index=prepared_index,
                prepared_state_root=prepared_state_root)
        if state_path.exists() or state_path.is_symlink():
            raise WgsAnnotationError("incomplete VEP plan checkpoint is unavailable")
    connection = _open_index(state_path, execution_config)
    partitions: list[VepTransportPartition] = []
    batch: list[VepTransportRequest] = []
    batch_transport_bytes = len(VEP_HEADER)
    batch_request_bytes = 0
    retained_count = sequence_count = native_count = occurrence_count = 0

    def flush() -> None:
        nonlocal batch, batch_transport_bytes, batch_request_bytes
        if not batch:
            return
        ordinal = len(partitions) + 1
        transport = publish_stream(
            root,
            f"transport-{ordinal:06d}.vcf",
            (chunk for request in batch for chunk in (
                VEP_HEADER if request is batch[0] else b"",
                request.record.line().encode("ascii"),
            ) if chunk),
            maximum_bytes=execution_config.annotation.max_partition_bytes,
        )
        expected = publish_stream(
            root,
            f"requests-{ordinal:06d}.jsonl",
            (request.as_bytes() for request in batch),
            maximum_bytes=execution_config.annotation.max_partition_bytes,
        )
        previous = partitions[-1].sha256 if partitions else None
        partitions.append(VepTransportPartition(
            ordinal,
            batch[0].request_ordinal,
            batch[-1].request_ordinal,
            len(batch),
            transport,
            expected,
            batch_transport_bytes,
            batch_request_bytes,
            previous,
        ))
        batch = []
        batch_transport_bytes = len(VEP_HEADER)
        batch_request_bytes = 0

    try:
        connection.executescript("""
            CREATE TABLE identities(identity_id TEXT PRIMARY KEY, ordinal INTEGER UNIQUE);
            CREATE TABLE sources(source_allele_id TEXT, identity_id TEXT NOT NULL, native INTEGER,
                PRIMARY KEY(source_allele_id,identity_id));
        """)
        with stage_lock(root, "vep-transport"):
            for allele in retained:
                if (
                    not isinstance(allele, RetainedAllele)
                    or allele.ordinal != retained_count + 1
                    or allele.parent_source_manifest_sha256 != source_manifest.sha256
                    or allele.parent_working_set_manifest_sha256 != working_manifest.sha256
                ):
                    raise WgsAnnotationError("retained transport order is invalid")
                try:
                    connection.execute(
                        "INSERT INTO identities VALUES (?,?)",
                        (allele.identity_id, allele.ordinal),
                    )
                    for source_id in allele.source_allele_ids:
                        previous = connection.execute("SELECT native FROM sources WHERE source_allele_id=?",
                                                      (source_id,)).fetchone()
                        native = allele.contig is None
                        # An authenticated intake can retain native locus/member
                        # views of one occurrence. Sequence requests still have
                        # exactly one owner, and standalone input cannot expand.
                        if previous is not None and (intake_binding is None or not native or not previous[0]):
                            raise WgsAnnotationError("retained transport source is duplicated")
                        connection.execute("INSERT INTO sources VALUES (?,?,?)",
                                           (source_id, allele.identity_id, native))
                except sqlite3.Error:
                    raise WgsAnnotationError("retained transport identity is duplicated") from None
                retained_count += 1
                occurrence_count += len(allele.source_allele_ids)
                if allele.contig is None:
                    native_count += 1
                    continue
                request = _request(sequence_count + 1, allele)
                line_bytes = len(request.record.line().encode("ascii"))
                request_bytes = len(request.as_bytes())
                limit = execution_config.annotation.max_partition_bytes
                if len(VEP_HEADER) + line_bytes > limit or request_bytes > limit:
                    raise WgsAnnotationError("VEP request exceeds partition capacity")
                if batch and (
                    len(batch) >= execution_config.annotation.max_partition_requests
                    or batch_transport_bytes + line_bytes > limit
                    or batch_request_bytes + request_bytes > limit
                ):
                    flush()
                batch.append(request)
                batch_transport_bytes += line_bytes
                batch_request_bytes += request_bytes
                sequence_count += 1
                if retained_count % execution_config.sqlite.cursor_batch_rows == 0:
                    connection.commit()
            flush()
            connection.commit()
            if (
                retained_count != (intake_binding.retained_identities if intake_binding
                                   else working_manifest.summary.retained_unique_identities)
                or occurrence_count != (intake_binding.retained_occurrences if intake_binding
                                       else working_manifest.summary.retained_occurrences)
                or sequence_count + native_count != retained_count
            ):
                raise WgsAnnotationError("VEP partition source coverage is incomplete")
            plan = VepTransportPlan(
                source_manifest.sha256,
                working_manifest.sha256,
                limits_sha256,
                _plan_parents(
                    source_manifest, working_manifest, limits_sha256,
                    execution_config.version, intake_binding,
                ),
                retained_count,
                sequence_count,
                native_count,
                tuple(partitions),
                intake_binding,
            )
            manifest = publish_bytes(root, "manifest.json", plan.as_bytes())
            return VepTransportPublication(plan, manifest)
    finally:
        connection.close()
        try:
            state_path.unlink()
        except (FileNotFoundError, OSError):
            pass


def verify_vep_partition_plan(
    *,
    root: str | Path,
    plan: VepTransportPlan,
    source_manifest: LedgerManifest,
    working_manifest: WorkingSetManifest,
    execution_config: BoundedWgsExecutionConfig,
    intake_binding: AnnotationIntakeBinding | None = None,
    prepared_index: PreparedIndex | None = None,
    prepared_state_root: str | Path | None = None,
) -> int:
    """Authenticate every planned request and its exact VCF transport line."""

    root = Path(root)
    limits_sha256 = hashlib.sha256(
        bounded_wgs_execution_config_as_json(execution_config)
    ).hexdigest()
    if (
        not root.is_absolute()
        or not isinstance(plan, VepTransportPlan)
        or plan.parent_source_manifest_sha256 != source_manifest.sha256
        or plan.parent_working_set_manifest_sha256 != working_manifest.sha256
        or plan.execution_limits_sha256 != limits_sha256
        or plan.parents
        != _plan_parents(
            source_manifest, working_manifest, limits_sha256,
            execution_config.version, intake_binding,
        )
        or plan.intake_binding != intake_binding
    ):
        raise WgsAnnotationError("VEP partition plan binding is invalid")
    if prepared_index is not None:
        if prepared_state_root is None:
            raise WgsAnnotationError("VEP prepared projection state is missing")
        _require_prepared_projection_parent(prepared_index, Path(prepared_state_root),
                                            source_manifest, working_manifest)
    elif prepared_state_root is not None:
        raise WgsAnnotationError("VEP prepared projection state is unbound")
    content = plan.as_bytes()
    artifact = Artifact(
        "manifest.json", hashlib.sha256(content).hexdigest(), len(content),
    )
    if read_artifact(root, artifact) != content:
        raise WgsAnnotationError("VEP partition manifest is invalid")
    observed = 0
    for partition in plan.partitions:
        verify_artifact(
            root, partition.transport,
            maximum_bytes=execution_config.annotation.max_partition_bytes,
        )
        verify_artifact(
            root, partition.expected_requests,
            maximum_bytes=execution_config.annotation.max_partition_bytes,
        )
        requests = (
            _load_request(line)
            for line in iter_artifact_lines(
                root,
                partition.expected_requests,
                maximum_bytes=execution_config.annotation.max_partition_bytes,
            )
        )
        transport = iter_artifact_lines(
            root,
            partition.transport,
            maximum_bytes=execution_config.annotation.max_partition_bytes,
        )
        try:
            if next(transport) != b"##fileformat=VCFv4.2\n" or next(transport) != VEP_COLUMNS:
                raise WgsAnnotationError("VEP transport header is invalid")
        except StopIteration:
            raise WgsAnnotationError("VEP transport header is incomplete") from None
        count = 0
        for request in requests:
            observed += 1
            count += 1
            if (
                request.request_ordinal != observed
                or request.parent_working_set_manifest_sha256 != working_manifest.sha256
            ):
                raise WgsAnnotationError("VEP request sequence is invalid")
            try:
                if next(transport) != request.record.line().encode("ascii"):
                    raise WgsAnnotationError("VEP transport request mismatch")
            except StopIteration:
                raise WgsAnnotationError("VEP transport is incomplete") from None
        try:
            next(transport)
        except StopIteration:
            pass
        else:
            raise WgsAnnotationError("VEP transport has extra records")
        if count != partition.requests:
            raise WgsAnnotationError("VEP partition request count is invalid")
    if observed != plan.sequence_requests:
        raise WgsAnnotationError("VEP plan request coverage is invalid")
    return observed


def iter_partition_requests(
    root: str | Path,
    partition: VepTransportPartition,
    *,
    execution_config: BoundedWgsExecutionConfig,
) -> Iterator[VepTransportRequest]:
    """Yield an authenticated partition sidecar under the declared row bound."""

    root = Path(root)
    if (
        not root.is_absolute()
        or not isinstance(partition, VepTransportPartition)
        or not isinstance(execution_config, BoundedWgsExecutionConfig)
    ):
        raise WgsAnnotationError("VEP request cursor binding is invalid")
    count = 0
    for content in iter_artifact_lines(
        root,
        partition.expected_requests,
        maximum_bytes=execution_config.annotation.max_partition_bytes,
        maximum_line_bytes=MAX_ROW_BYTES,
    ):
        request = _load_request(content)
        count += 1
        if request.request_ordinal != partition.first_request_ordinal + count - 1:
            raise WgsAnnotationError("VEP partition request order is invalid")
        yield request
    if count != partition.requests:
        raise WgsAnnotationError("VEP partition request coverage is invalid")


def _open_private_output(path: Path, root: Path) -> int:
    descriptor = -1
    try:
        if path.parent != root or not path.name or path.name in {".", ".."}:
            raise ValueError
        descriptor = os.open(
            path,
            os.O_WRONLY | os.O_CREAT | os.O_EXCL | os.O_NOFOLLOW,
            0o600,
        )
        os.fchmod(descriptor, 0o600)
        metadata = os.fstat(descriptor)
        if (
            not stat.S_ISREG(metadata.st_mode)
            or metadata.st_uid != os.getuid()
            or stat.S_IMODE(metadata.st_mode) != 0o600
            or metadata.st_nlink != 1
        ):
            raise ValueError
        return descriptor
    except (OSError, ValueError):
        if descriptor >= 0:
            os.close(descriptor)
        raise WgsAnnotationError("cannot create bounded child output") from None


def _write_all(descriptor: int, content: bytes) -> None:
    view = memoryview(content)
    while view:
        written = os.write(descriptor, view)
        if written <= 0:
            raise WgsAnnotationError("cannot write bounded child output")
        view = view[written:]


def _terminate_group(process: subprocess.Popen, grace_seconds: int) -> None:
    try:
        os.killpg(process.pid, signal.SIGTERM)
    except ProcessLookupError:
        pass
    try:
        process.wait(timeout=grace_seconds)
    except subprocess.TimeoutExpired:
        pass
    # The leader may have exited while a forked VEP worker remains in its group.
    try:
        os.killpg(process.pid, signal.SIGKILL)
    except ProcessLookupError:
        pass
    if process.poll() is None:
        try:
            process.wait(timeout=grace_seconds)
        except subprocess.TimeoutExpired:
            raise WgsAnnotationError("bounded child could not be reaped") from None
    # Admission owns an active subreaper. Reap only this child's process
    # group so cancellation cannot leave adopted native-query descendants.
    deadline = time.monotonic() + grace_seconds
    while True:
        try:
            child, _ = os.waitpid(-process.pid, os.WNOHANG)
        except ChildProcessError:
            break
        if child == 0:
            if time.monotonic() >= deadline:
                raise WgsAnnotationError("bounded child descendants could not be reaped")
            time.sleep(.01)


def run_bounded_child(
    command: tuple[str, ...],
    *,
    cwd: str | Path,
    stdout_path: str | Path,
    stderr_path: str | Path,
    execution_config: BoundedWgsExecutionConfig,
    environment: dict[str, str] | None = None,
    query_budget: QueryBudget | None = None,
    output_files: tuple[Path, ...] = (),
) -> BoundedChildResult:
    """Run one child while draining both pipes under wall/output limits."""

    root = Path(cwd)
    stdout_path, stderr_path = Path(stdout_path), Path(stderr_path)
    if (
        not isinstance(execution_config, BoundedWgsExecutionConfig)
        or not root.is_absolute()
        or not root.is_dir()
        or not isinstance(command, tuple)
        or not command
        or not Path(command[0]).is_absolute()
        or any(not isinstance(value, str) or not value or "\0" in value for value in command)
        or stdout_path == stderr_path
        or environment is not None and (
            type(environment) is not dict
            or any(
                not isinstance(key, str) or not key
                or not isinstance(value, str) or "\0" in value
                for key, value in environment.items()
            )
        )
    ):
        raise WgsAnnotationError("invalid bounded child configuration")
    output_limit = execution_config.child_process.max_total_output_bytes
    if query_budget is not None:
        query_budget.check()
        output_limit = min(output_limit, query_budget.limits.max_response_bytes)
    stdout_fd = _open_private_output(stdout_path, root)
    try:
        stderr_fd = _open_private_output(stderr_path, root)
    except Exception:
        os.close(stdout_fd)
        raise
    process = None
    selector = selectors.DefaultSelector()
    digests = {"stdout": hashlib.sha256(), "stderr": hashlib.sha256()}
    sizes = {"stdout": 0, "stderr": 0}
    failure = None
    failure_category = "unexpected"
    terminated = False
    try:
        process = subprocess.Popen(
            command,
            cwd=root,
            env=environment,
            stdin=subprocess.DEVNULL,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            start_new_session=True,
            close_fds=True,
        )
        streams = {
            process.stdout.fileno(): ("stdout", process.stdout, stdout_fd),
            process.stderr.fileno(): ("stderr", process.stderr, stderr_fd),
        }
        for descriptor, (name, stream, output_fd) in streams.items():
            os.set_blocking(descriptor, False)
            selector.register(descriptor, selectors.EVENT_READ, (name, stream, output_fd))
        deadline = time.monotonic() + execution_config.child_process.wall_timeout_seconds
        if query_budget is not None:
            deadline = min(deadline, query_budget.deadline)
        while selector.get_map() or process.poll() is None:
            from genosapiens.wgs_query_limits import _query_cancelled
            if failure is None and _query_cancelled():
                failure = "bounded child interrupted"
                failure_category = "interrupted"
            if failure is None and any(
                path.is_symlink() or path.exists() and path.stat().st_size > execution_config.annotation.max_shard_bytes
                for path in output_files
            ):
                failure = "bounded child file output limit exceeded"
                failure_category = "response_capacity"
            if failure is None and time.monotonic() >= deadline:
                failure = "bounded child timed out"
                failure_category = "timeout"
            if failure is not None and not terminated:
                _terminate_group(process, execution_config.child_process.kill_grace_seconds)
                terminated = True
            for key, _mask in selector.select(timeout=0.05):
                name, stream, output_fd = key.data
                try:
                    content = os.read(
                        key.fd, min(execution_config.child_process.pipe_chunk_bytes, output_limit + 1),
                    )
                except BlockingIOError:
                    continue
                if not content:
                    selector.unregister(key.fd)
                    stream.close()
                    continue
                proposed_total = sum(sizes.values()) + len(content)
                proposed_stream = sizes[name] + len(content)
                stream_limit = (
                    min(execution_config.child_process.max_stderr_bytes, output_limit)
                    if name == "stderr"
                    else output_limit
                )
                if failure is None and (
                    proposed_total > output_limit
                    or proposed_stream > stream_limit
                ):
                    failure = "bounded child output limit exceeded"
                    failure_category = "response_capacity"
                remaining_total = max(
                    0,
                    output_limit
                    - sum(sizes.values()),
                )
                remaining_stream = max(0, stream_limit - sizes[name])
                retained = content[:min(len(content), remaining_total, remaining_stream)]
                if retained:
                    _write_all(output_fd, retained)
                    digests[name].update(retained)
                    sizes[name] += len(retained)
        exit_code = process.wait()
        if failure is None and any(
            path.is_symlink() or path.exists() and path.stat().st_size > execution_config.annotation.max_shard_bytes
            for path in output_files
        ):
            failure = "bounded child file output limit exceeded"
            failure_category = "response_capacity"
        if failure is None and exit_code != 0:
            failure = "bounded child exited unsuccessfully"
        if failure is not None:
            if not terminated:
                _terminate_group(process, execution_config.child_process.kill_grace_seconds)
            if query_budget is not None:
                raise WgsFailure(
                    query_budget.origin, failure_category,
                    operation="annotation_query", exception_kind="project_failure",
                )
            raise WgsAnnotationError(failure, category=failure_category)
        os.fsync(stdout_fd)
        os.fsync(stderr_fd)
        return BoundedChildResult(
            exit_code,
            digests["stdout"].hexdigest(),
            sizes["stdout"],
            digests["stderr"].hexdigest(),
            sizes["stderr"],
            sum(sizes.values()),
        )
    except WgsAnnotationError:
        raise
    except (OSError, ValueError, subprocess.SubprocessError):
        if process is not None:
            _terminate_group(process, execution_config.child_process.kill_grace_seconds)
        raise WgsAnnotationError("bounded child execution failed") from None
    finally:
        selector.close()
        for descriptor in (stdout_fd, stderr_fd):
            try:
                os.fsync(descriptor)
            except OSError:
                pass
            os.close(descriptor)
        if process is not None:
            for stream in (process.stdout, process.stderr):
                if stream is not None and not stream.closed:
                    stream.close()


def _vep_native_record(record: VepVcfRecord) -> VepVcfRecord:
    """Render the canonical mitochondrial name in the pinned Ensembl namespace.

    The plan/request and source identities remain canonical M. Only the VEP
    boundary uses MT, matching its pinned FASTA/cache and the ordinary adapter.
    No coordinates, alleles, reference roles or other contig names change.
    """
    return replace(record, contig="MT") if record.contig == "M" else record


@contextmanager
def _vep_partition_input(partition_root, partition, output_root, config, attempt):
    name = f"vep-input-{partition.ordinal:06d}-attempt-{attempt:02d}.vcf"
    path = output_root / name
    def rows():
        for line in iter_artifact_lines(partition_root, partition.transport,
                maximum_bytes=config.annotation.max_partition_bytes):
            yield b"MT\t" + line[2:] if line.startswith(b"M\t") else line
    # publish_stream authenticates and closes the bounded source before VEP.
    # Keep the original plan immutable, including on last-partition recovery.
    artifact = publish_stream(output_root, name, rows(),
        maximum_bytes=config.annotation.max_partition_bytes)
    try:
        yield output_root / artifact.name
    finally:
        path.unlink()


def run_vep_partition_child(
    *,
    partition_root: str | Path,
    partition: VepTransportPartition,
    output_root: str | Path,
    cache_directory: str | Path,
    fasta_path: str | Path,
    resource_session_sha256: str,
    execution_config: BoundedWgsExecutionConfig,
    attempt_number: int = 1,
) -> BoundedChildResult:
    """Invoke the fixed offline VEP command for one authenticated partition."""

    partition_root, output_root = Path(partition_root), Path(output_root)
    cache_directory, fasta_path = Path(cache_directory), Path(fasta_path)
    executable = Path(VEP_EXECUTABLE)
    try:
        if (
            not partition_root.is_absolute()
            or not output_root.is_absolute()
            or not cache_directory.is_absolute()
            or not fasta_path.is_absolute()
            or not isinstance(partition, VepTransportPartition)
            or not isinstance(execution_config, BoundedWgsExecutionConfig)
            or not _digest(resource_session_sha256)
            or not _count(attempt_number, 1)
            or any(name != "lo" for _index, name in socket.if_nameindex())
            or not output_root.is_dir()
            or output_root.is_symlink()
            or not cache_directory.is_dir()
            or cache_directory.is_symlink()
            or not fasta_path.is_file()
            or fasta_path.is_symlink()
            or not executable.is_absolute()
            or not executable.is_file()
            or executable.is_symlink()
        ):
            raise ValueError
    except (OSError, ValueError):
        raise WgsAnnotationError("offline VEP partition binding is invalid") from None
    verify_artifact(
        partition_root,
        partition.transport,
        maximum_bytes=execution_config.annotation.max_partition_bytes,
    )
    verify_artifact(
        partition_root,
        partition.expected_requests,
        maximum_bytes=execution_config.annotation.max_partition_bytes,
    )
    environment = {
        "PATH": "/usr/local/bin:/usr/bin:/bin",
        "HOME": str(output_root),
        "TMPDIR": str(output_root),
        "LANG": "C",
        "LC_ALL": "C",
        "PERL5LIB": "/opt/vep/src/ensembl-vep:/opt/vep/src/ensembl-vep/modules",
        "PERL_HASH_SEED": "0",
        "PERL_PERTURB_KEYS": "0",
    }
    command = (
        str(executable),
        *VEP_FLAGS,
        "--buffer_size", str(execution_config.annotation.buffer_size),
        "--fork", str(execution_config.annotation.fork_processes),
        "--dir", str(output_root),
        "--dir_plugins", str(cache_directory),
        "--dir_cache", str(cache_directory),
        "--fasta", str(fasta_path),
        "--warning_file", "STDERR",
        "--output_file", "STDOUT",
    )
    with _vep_partition_input(partition_root, partition, output_root,
            execution_config, attempt_number) as input_path:
        return run_bounded_child(
            (*command, "--input_file", str(input_path)),
            cwd=output_root,
            stdout_path=output_root / (
                f"annotated-{partition.ordinal:06d}-attempt-{attempt_number:02d}.vcf"
            ),
            stderr_path=output_root / (
                f"vep-{partition.ordinal:06d}-attempt-{attempt_number:02d}.stderr"
            ),
            execution_config=execution_config,
            environment=environment,
        )


def _open_stable_response(path: Path):
    try:
        descriptor = os.open(path, os.O_RDONLY | os.O_NOFOLLOW)
        stream = os.fdopen(descriptor, "rb", closefd=True)
        metadata = os.fstat(stream.fileno())
        named = os.stat(path, follow_symlinks=False)
        if (
            not stat.S_ISREG(metadata.st_mode)
            or metadata.st_uid != os.getuid()
            or stat.S_IMODE(metadata.st_mode) != 0o600
            or metadata.st_nlink != 1
            or (metadata.st_dev, metadata.st_ino) != (named.st_dev, named.st_ino)
        ):
            stream.close()
            raise ValueError
        return stream, metadata
    except (OSError, ValueError):
        raise WgsAnnotationError("VEP response file is unsafe") from None


def _parse_one_response(
    header: bytes,
    raw_line: bytes,
    request: VepTransportRequest,
    *,
    release: str,
    bundle_manifest_sha256: str,
    output_sha256: str,
    parsed_header: ParsedVepHeader,
) -> ParsedVep:
    try:
        content = (header + raw_line).decode("utf-8", errors="strict")
        source = SourcedRecord(
            "Ensembl VEP",
            release,
            request.request_id,
            hashlib.sha256(content.encode("utf-8")).hexdigest(),
            bundle_manifest_sha256,
            {
                "partition_output_sha256": output_sha256,
                "transport_request_id": request.request_id,
                "working_set_manifest_sha256":
                    request.parent_working_set_manifest_sha256,
            },
        )
        return parse_vep_vcf(
            content,
            expected=(_vep_native_record(request.record),),
            source=source,
            required_fields=RICH_FIELDS,
            parsed_header=parsed_header,
        )
    except (UnicodeError, ValueError, TypeError, RuntimeError):
        raise WgsAnnotationError("VEP response row is invalid") from None


def _parsed_record_bytes(record: VerifiedVepRecord) -> bytes:
    key = (record.request.request_id, 1)
    return _owned_canonical_bytes({
        "schema": "genosapiens.wgs-parsed-vep-record", "version": 1,
        "request": json.loads(record.request.as_bytes()),
        "fields": record.parsed.fields,
        "transcripts": record.parsed.transcripts[key],
        "non_transcripts": record.parsed.non_transcripts[key],
        "raw_line_sha256": record.raw_line_sha256,
        "parent_output_sha256": record.parent_output_sha256,
    })


def _load_parsed_record(content: bytes) -> VerifiedVepRecord | UnassessedVepRecord:
    # Every ordinary saved-record decode starts an independent validation walk.
    with _validation_session(isolated=True):
        return _load_parsed_record_owned(content)


@_validation_session()
def _load_parsed_record_owned(content: bytes) -> VerifiedVepRecord | UnassessedVepRecord:
    """Strict decode within the producer's one bounded synchronous row owner."""
    try:
        value = json.loads(content)
        if type(value) is not dict:
            raise ValueError
        if value.get('schema') == 'genosapiens.wgs-unassessed-vep-record':
            record = UnassessedVepRecord(_load_request(canonical(value['request'])),
                value['parent_output_sha256'], value['reference_fai_sha256'], value['reason'])
            if record.as_bytes() != content:
                raise ValueError
            return record
        if (canonical(value) != content or set(value) != {
                "schema", "version", "request", "fields", "transcripts",
                "non_transcripts", "raw_line_sha256", "parent_output_sha256",
        } or value["schema"] != "genosapiens.wgs-parsed-vep-record"
                or value["version"] != 1):
            raise ValueError
        request = _load_request(canonical(value["request"]))
        key = (request.request_id, 1)
        parsed = ParsedVep(
            tuple(value["fields"]),
            {key: tuple(_typed_value(row, TranscriptConsequence, decode=True)
                        for row in value["transcripts"])},
            {key: tuple(_typed_value(row, SourcedRecord, decode=True)
                        for row in value["non_transcripts"])},
        )
        record = VerifiedVepRecord(request, parsed, value["raw_line_sha256"],
                                   value["parent_output_sha256"])
        if _parsed_record_bytes(record) != content:
            raise ValueError
        return record
    except (ValueError, TypeError, KeyError, RuntimeError, UnicodeError, LedgerError):
        raise WgsAnnotationError("parsed VEP record is invalid") from None


def _run_annotation_parse_worker(index, jobs, results):
    """Parse bounded messages only: no resource sessions, state or publication."""
    def interrupt(_signum, _frame):
        raise KeyboardInterrupt

    signal.signal(signal.SIGTERM, interrupt)
    generation = ordinal = None
    binding = None
    try:
        results.put(("ready", index, None, None, os.getpid()))
        while True:
            message = jobs.get()
            if message is None:
                return
            kind, generation, ordinal, payload = message
            if kind == "bind":
                binding = payload
                results.put(("bound", index, generation, None, None))
                del message, payload
                continue
            if kind == "drain":
                binding = None
                results.put(("drained", index, generation, None, None))
                del message, payload
                continue
            if kind != "parse" or binding is None:
                raise ValueError
            header, parsed_header, release, bundle_sha, output_sha, response_limit = binding
            if (type(payload) is not tuple or not 0 < len(payload) <= _ANNOTATION_PARSE_JOB_ROWS
                    or sum(len(request) + len(raw) for request, raw in payload)
                        > _ANNOTATION_PARSE_INPUT_BYTES):
                raise ValueError
            started = time.process_time_ns()
            records = []
            size = 0
            parse_calls = encode_calls = 0
            for request_bytes, raw in payload:
                if (type(request_bytes) is not bytes or len(request_bytes) > MAX_ROW_BYTES
                        or type(raw) is not bytes or len(raw) > response_limit):
                    raise ValueError
                request = _load_request(request_bytes)
                parse_calls += 1
                parsed = _parse_one_response(header, raw, request, release=release,
                    bundle_manifest_sha256=bundle_sha, output_sha256=output_sha,
                    parsed_header=parsed_header)
                record = VerifiedVepRecord(request, parsed, hashlib.sha256(raw).hexdigest(), output_sha)
                encode_calls += 1
                content = _parsed_record_bytes(record)
                if len(content) > MAX_ROW_BYTES or size + len(content) > _ANNOTATION_PARSE_RESULT_BYTES:
                    raise ValueError
                size += len(content)
                records.append(EncodedVerifiedVepRecord(request.request_id, request.request_ordinal,
                    request_bytes, content, hashlib.sha256(content).hexdigest(),
                    record.raw_line_sha256, output_sha))
                del request, parsed, record, content, request_bytes, raw
            results.put(("result", index, generation, ordinal,
                         (tuple(records), time.process_time_ns() - started, parse_calls, encode_calls)))
            del message, payload, records
    except BaseException:
        # Scientific payloads and exception text never cross the error channel.
        results.put(("error", index, generation, ordinal, None))
    finally:
        signal.signal(signal.SIGTERM, signal.SIG_IGN)


class _AnnotationParsePool:
    """Ten lazy spawned parsers, at most two jobs each and twenty in total.

    Results are consumed in arrival order into the coordinator's temporary
    SQLite store. Only its final ordered join drives publication. An idle pool
    may span VEP partitions; all jobs must be consumed before rebinding or VEP.
    """

    def __init__(self, config, *, _worker_count=_ANNOTATION_PARSE_WORKERS,
                 _worker_target=_run_annotation_parse_worker):
        if type(_worker_count) is not int or not 1 <= _worker_count <= _ANNOTATION_PARSE_WORKERS:
            raise WgsAnnotationError("annotation parser worker count is invalid")
        self.config = config
        self.worker_count = _worker_count
        context = multiprocessing.get_context("spawn")
        self.results = context.Queue(maxsize=_ANNOTATION_PARSE_MAX_JOBS)
        self.jobs = tuple(context.Queue(maxsize=2) for _ in range(_worker_count))
        self.workers = tuple(context.Process(target=_worker_target,
            args=(i, self.jobs[i], self.results), name=f"wgs-annotation-parser-{i}")
            for i in range(_worker_count))
        self.outstanding = {}
        self.loads = [0] * _worker_count
        self.started = self.closed = False
        self.generation = self.next_ordinal = self.next_worker = 0
        self.output_sha256 = None
        self.worker_pids = {}
        self.worker_rows = [0] * _worker_count
        self.parse_count = self.encode_count = self.cpu_nanoseconds = 0
        self.max_outstanding = self.max_rows = self.max_input_bytes = self.max_result_bytes = 0

    def _health(self, deadline):
        if self.closed or any(worker.exitcode is not None for worker in self.workers):
            raise WgsAnnotationError("annotation parser worker exited")
        if time.monotonic() >= deadline:
            raise WgsAnnotationError("annotation parser deadline exceeded", category="timeout")

    def _receive(self, deadline):
        while True:
            self._health(deadline)
            try:
                message = self.results.get(timeout=min(.1, max(.001, deadline - time.monotonic())))
            except queue.Empty:
                continue
            if (type(message) is not tuple or len(message) != 5
                    or type(message[1]) is not int or not 0 <= message[1] < self.worker_count):
                raise WgsAnnotationError("annotation parser message is malformed")
            if message[0] == "error":
                raise WgsAnnotationError("annotation parser failed")
            return message

    def bind(self, *, header, parsed_header, release, bundle_sha256, output_sha256):
        if self.closed or self.outstanding:
            raise WgsAnnotationError("annotation parser still has pending jobs")
        deadline = time.monotonic() + min(300, self.config.child_process.wall_timeout_seconds)
        if not self.started:
            self.started = True
            for worker in self.workers:
                worker.start()
            while len(self.worker_pids) < self.worker_count:
                kind, index, generation, ordinal, pid = self._receive(deadline)
                if (kind != "ready" or generation is not None or ordinal is not None
                        or index in self.worker_pids or pid != self.workers[index].pid):
                    raise WgsAnnotationError("annotation parser startup is invalid")
                self.worker_pids[index] = pid
        self.generation += 1
        self.output_sha256 = output_sha256
        binding = (header, parsed_header, release, bundle_sha256, output_sha256,
                   self.config.annotation.max_response_row_bytes)
        for jobs in self.jobs:
            jobs.put(("bind", self.generation, None, binding), timeout=1)
        ready = set()
        while len(ready) < self.worker_count:
            kind, index, generation, ordinal, value = self._receive(deadline)
            if (kind != "bound" or generation != self.generation or ordinal is not None
                    or value is not None or index in ready):
                raise WgsAnnotationError("annotation parser binding is invalid")
            ready.add(index)

    @property
    def full(self):
        return len(self.outstanding) >= min(_ANNOTATION_PARSE_MAX_JOBS, 2 * self.worker_count)

    def submit(self, rows):
        """Rows contain parent-authenticated request ID/ordinal/bytes and raw line."""
        size = sum(len(row[2]) + len(row[3]) for row in rows)
        if (not self.started or self.full or not 0 < len(rows) <= _ANNOTATION_PARSE_JOB_ROWS
                or size > _ANNOTATION_PARSE_INPUT_BYTES):
            raise WgsAnnotationError("annotation parser submission exceeds bounds")
        deadline = time.monotonic() + self.config.child_process.wall_timeout_seconds
        self._health(min([deadline, *(item[2] for item in self.outstanding.values())]))
        for offset in range(self.worker_count):
            index = (self.next_worker + offset) % self.worker_count
            if self.loads[index] < 2:
                break
        self.next_worker = (index + 1) % self.worker_count
        ordinal = self.next_ordinal
        self.next_ordinal += 1
        expected = tuple((identity, order, request, hashlib.sha256(raw).hexdigest())
                         for identity, order, request, raw in rows)
        self.outstanding[ordinal] = (index, expected, deadline)
        self.loads[index] += 1
        self.jobs[index].put(("parse", self.generation, ordinal,
                             tuple((row[2], row[3]) for row in rows)), timeout=1)
        self.max_outstanding = max(self.max_outstanding, len(self.outstanding))
        self.max_rows = max(self.max_rows, len(rows))
        self.max_input_bytes = max(self.max_input_bytes, size)

    def take(self):
        if not self.outstanding:
            raise WgsAnnotationError("annotation parser has no pending result")
        kind, index, generation, ordinal, value = self._receive(
            min(item[2] for item in self.outstanding.values()))
        if (kind != "result" or generation != self.generation or type(ordinal) is not int
                or ordinal not in self.outstanding or type(value) is not tuple or len(value) != 4):
            raise WgsAnnotationError("annotation parser result is unexpected")
        owner, expected, _ = self.outstanding[ordinal]
        records, cpu_ns, parse_calls, encode_calls = value
        if (owner != index or type(records) is not tuple or len(records) != len(expected)
                or type(cpu_ns) is not int or cpu_ns < 0
                or type(parse_calls) is not int or type(encode_calls) is not int
                or parse_calls != len(records) or encode_calls != len(records)):
            raise WgsAnnotationError("annotation parser result binding is invalid")
        size = 0
        for record, (identity, order, request, raw_sha) in zip(records, expected, strict=True):
            if type(record) is not EncodedVerifiedVepRecord:
                raise WgsAnnotationError("annotation parser result type is invalid")
            record.__post_init__()
            size += len(record.parsed_record)
            if (size > _ANNOTATION_PARSE_RESULT_BYTES or record.request_id != identity
                    or record.request_ordinal != order or record.canonical_request != request
                    or record.raw_line_sha256 != raw_sha or record.parent_output_sha256 != self.output_sha256
                    or hashlib.sha256(record.parsed_record).hexdigest() != record.parsed_record_sha256):
                raise WgsAnnotationError("annotation parser record binding is invalid")
        self.max_result_bytes = max(self.max_result_bytes, size)
        self.worker_rows[index] += len(records)
        self.parse_count += parse_calls
        self.encode_count += encode_calls
        self.cpu_nanoseconds += cpu_ns
        self.loads[index] -= 1
        del self.outstanding[ordinal]
        return records

    def finish_partition(self):
        """Fence every worker and reject duplicate/trailing results before closure."""
        if self.outstanding:
            raise WgsAnnotationError("annotation parser still has pending jobs")
        deadline = time.monotonic() + min(300, self.config.child_process.wall_timeout_seconds)
        for jobs in self.jobs:
            jobs.put(("drain", self.generation, None, None), timeout=1)
        drained = set()
        while len(drained) < self.worker_count:
            kind, index, generation, ordinal, value = self._receive(deadline)
            if (kind != "drained" or generation != self.generation or ordinal is not None
                    or value is not None or index in drained):
                raise WgsAnnotationError("annotation parser closure is invalid")
            drained.add(index)

    def close(self):
        if self.closed:
            return
        self.closed = True
        owned = [worker for worker in self.workers if worker.pid is not None]
        for index, worker in enumerate(self.workers):
            if worker.pid is None or not worker.is_alive():
                continue
            if self.outstanding:
                worker.terminate()
            else:
                try:
                    self.jobs[index].put_nowait(None)
                except queue.Full:
                    worker.terminate()
        deadline = time.monotonic() + self.config.child_process.kill_grace_seconds
        while any(worker.is_alive() for worker in owned) and time.monotonic() < deadline:
            try:
                self.results.get(timeout=.05)
            except queue.Empty:
                pass
        for worker in owned:
            if worker.is_alive():
                worker.kill()
        reap_deadline = time.monotonic() + 1
        for worker in owned:
            worker.join(timeout=max(0, reap_deadline - time.monotonic()))
        for channel in (*self.jobs, self.results):
            channel.cancel_join_thread()
            channel.close()
        self.outstanding.clear()
        if any(worker.is_alive() for worker in owned):
            raise WgsAnnotationError("annotation parser cleanup failed")


def iter_verified_vep_records(
    expected_requests: Iterable[VepTransportRequest],
    *,
    response_path: str | Path,
    expected_output_sha256: str,
    bundle_manifest_sha256: str,
    release: str,
    scratch_path: str | Path,
    execution_config: BoundedWgsExecutionConfig,
    _parse_pool: _AnnotationParsePool | None = None,
    reference_coverage: ReferenceCoverage | None = None,
) -> Iterator[EncodedVerifiedVepRecord]:
    """Close all raw responses, then authenticate source-ordered encoded bytes.

    First publication never reconstructs typed evidence from this private spool.
    Actual persisted annotation readers retain `_load_parsed_record` validation.
    """

    response_path, scratch_path = Path(response_path), Path(scratch_path)
    if (
        not response_path.is_absolute()
        or not scratch_path.is_absolute()
        or not _digest(expected_output_sha256)
        or not _digest(bundle_manifest_sha256)
        or not isinstance(release, str)
        or not release
        or not isinstance(execution_config, BoundedWgsExecutionConfig)
    ):
        raise WgsAnnotationError("VEP response verification binding is invalid")
    connection = _open_index(scratch_path, execution_config)
    pool = _parse_pool if _parse_pool is not None else _AnnotationParsePool(execution_config)
    stream = None
    try:
        connection.executescript("""
            CREATE TABLE expected(
                request_id TEXT PRIMARY KEY,
                request_ordinal INTEGER NOT NULL UNIQUE,
                canonical_request BLOB NOT NULL
            );
            CREATE TABLE responses(
                request_id TEXT PRIMARY KEY REFERENCES expected(request_id),
                parsed_record BLOB NOT NULL,
                parsed_record_sha256 TEXT NOT NULL,
                raw_line_sha256 TEXT,
                parent_output_sha256 TEXT NOT NULL,
                metadata_sha256 TEXT NOT NULL
            );
        """)
        expected_count = 0
        previous_ordinal = None
        for request in expected_requests:
            if (
                not isinstance(request, VepTransportRequest)
                or previous_ordinal is not None
                and request.request_ordinal != previous_ordinal + 1
            ):
                raise WgsAnnotationError("VEP expected request order is invalid")
            try:
                connection.execute(
                    "INSERT INTO expected VALUES (?,?,?)",
                    (request.request_id, request.request_ordinal, request.as_bytes()),
                )
            except sqlite3.Error:
                raise WgsAnnotationError("VEP expected request identity is duplicated") from None
            expected_count += 1
            previous_ordinal = request.request_ordinal
            if expected_count % execution_config.sqlite.cursor_batch_rows == 0:
                connection.commit()
        if expected_count == 0:
            raise WgsAnnotationError("VEP response verification requires requests")
        connection.commit()

        stream, initial = _open_stable_response(response_path)
        if initial.st_size > execution_config.child_process.max_total_output_bytes:
            raise WgsAnnotationError("VEP response exceeds total output capacity")
        digest = hashlib.sha256()
        total = 0
        header_parts: list[bytes] = []
        header = None
        parsed_header = None
        header_bytes = 0
        columns = False
        response_count = 0
        pending = []
        pending_size = 0

        def consume():
            for encoded in pool.take():
                try:
                    connection.execute("INSERT INTO responses VALUES (?,?,?,?,?,?)",
                        (encoded.request_id, encoded.parsed_record, encoded.parsed_record_sha256,
                         encoded.raw_line_sha256, encoded.parent_output_sha256,
                         encoded.metadata_sha256()))
                except sqlite3.IntegrityError:
                    raise WgsAnnotationError("VEP response identity is duplicated") from None

        def submit():
            nonlocal pending_size
            if pool.full:
                consume()
            pool.submit(pending)
            pending.clear()
            pending_size = 0

        while True:
            limit = (
                execution_config.annotation.max_response_row_bytes
                if columns else execution_config.input.max_header_bytes
            )
            line = stream.readline(limit + 1)
            if not line:
                break
            total += len(line)
            digest.update(line)
            if len(line) > limit or not line.endswith(b"\n"):
                raise WgsAnnotationError("VEP response row exceeds its bounded size")
            if not columns:
                header_bytes += len(line)
                if header_bytes > execution_config.input.max_header_bytes:
                    raise WgsAnnotationError("VEP response header exceeds its bounded size")
                header_parts.append(line)
                if line == VEP_COLUMNS:
                    columns = True
                    header = b"".join(header_parts)
                    try:
                        parsed_header = parse_vep_header(
                            header.decode("utf-8", errors="strict"), required_fields=RICH_FIELDS)
                    except (UnicodeError, ValueError, RuntimeError):
                        raise WgsAnnotationError("VEP response header is invalid") from None
                    pool.bind(header=header, parsed_header=parsed_header, release=release,
                        bundle_sha256=bundle_manifest_sha256, output_sha256=expected_output_sha256)
                elif not line.startswith(b"#"):
                    raise WgsAnnotationError("VEP response header is invalid")
                continue
            fields = line.split(b"\t", 3)
            if len(fields) != 4:
                raise WgsAnnotationError("VEP response row is malformed")
            try:
                request_id = fields[2].decode("ascii")
            except UnicodeError:
                raise WgsAnnotationError("VEP response identity is invalid") from None
            stored = connection.execute(
                "SELECT request_ordinal, canonical_request FROM expected WHERE request_id=?",
                (request_id,),
            ).fetchone()
            if stored is None:
                raise WgsAnnotationError("VEP response contains an extra identity")
            row_size = len(stored[1]) + len(line)
            if row_size > _ANNOTATION_PARSE_INPUT_BYTES:
                raise WgsAnnotationError("VEP response job exceeds capacity")
            if pending and (len(pending) >= _ANNOTATION_PARSE_JOB_ROWS
                            or pending_size + row_size > _ANNOTATION_PARSE_INPUT_BYTES):
                submit()
            pending.append((request_id, stored[0], stored[1], line))
            pending_size += row_size
            response_count += 1
            if response_count % execution_config.sqlite.cursor_batch_rows == 0:
                connection.commit()
        if pending:
            submit()
        while pool.outstanding:
            consume()
        if pool.started:
            pool.finish_partition()
        final = os.fstat(stream.fileno())
        named = os.stat(response_path, follow_symlinks=False)
        if (
            not columns
            or total != initial.st_size
            or digest.hexdigest() != expected_output_sha256
            or (initial.st_dev, initial.st_ino, initial.st_size,
                initial.st_mtime_ns, initial.st_ctime_ns)
            != (final.st_dev, final.st_ino, final.st_size,
                final.st_mtime_ns, final.st_ctime_ns)
            or (final.st_dev, final.st_ino) != (named.st_dev, named.st_ino)
        ):
            raise WgsAnnotationError("VEP response file changed during verification")
        connection.commit()
        missing = connection.execute("""
            SELECT e.canonical_request FROM expected e
            LEFT JOIN responses r USING(request_id) WHERE r.request_id IS NULL
            ORDER BY e.request_ordinal
        """)
        for (request_bytes,) in missing:
            request = _load_request(request_bytes)
            reason = (reference_coverage.unassessed_reason(request.record.contig)
                      if type(reference_coverage) is ReferenceCoverage else None)
            if reason is None:
                raise WgsAnnotationError("VEP response coverage is incomplete", category="missing_response")
            record = UnassessedVepRecord(request, expected_output_sha256, reference_coverage.sha256, reason)
            content = record.as_bytes()
            encoded = EncodedUnassessedVepRecord(request.request_id, request.request_ordinal,
                request_bytes, content, hashlib.sha256(content).hexdigest(), None, expected_output_sha256)
            connection.execute('INSERT INTO responses VALUES (?,?,?,?,?,?)',
                (encoded.request_id, content, encoded.parsed_record_sha256, None,
                 expected_output_sha256, encoded.metadata_sha256()))
            response_count += 1
        if response_count != expected_count:
            raise WgsAnnotationError("VEP response coverage is incomplete", category="missing_response")
        cursor = connection.execute("""
            SELECT e.request_id, e.request_ordinal, e.canonical_request,
                   r.parsed_record, r.parsed_record_sha256, r.raw_line_sha256,
                   r.parent_output_sha256, r.metadata_sha256 FROM expected e
            JOIN responses r USING(request_id) ORDER BY e.request_ordinal
        """)
        for values in _rows(
            cursor, execution_config.sqlite.cursor_batch_rows,
        ):
            encoded_type = EncodedUnassessedVepRecord if values[5] is None else EncodedVerifiedVepRecord
            encoded = encoded_type(*values[:-1])
            if (hashlib.sha256(encoded.parsed_record).hexdigest() != encoded.parsed_record_sha256
                    or encoded.parent_output_sha256 != expected_output_sha256
                    or encoded.metadata_sha256() != values[-1]):
                raise WgsAnnotationError("encoded VEP spool binding is invalid")
            yield encoded
    except sqlite3.Error:
        raise WgsAnnotationError("VEP response reconciliation failed") from None
    finally:
        try:
            if _parse_pool is None:
                pool.close()
        finally:
            if stream is not None:
                stream.close()
            connection.close()
            try:
                scratch_path.unlink()
            except (FileNotFoundError, OSError):
                pass


def _artifact_value(value: object) -> Artifact:
    try:
        if type(value) is not dict or set(value) != {"name", "sha256", "size_bytes"}:
            raise ValueError
        return Artifact(**value)
    except (TypeError, ValueError, LedgerError):
        raise WgsAnnotationError("invalid annotation artifact binding") from None


def _load_parsed_partition(content: bytes) -> ParsedPartitionManifest:
    try:
        value = json.loads(content)
        if (canonical(value) != content or set(value) != {
                "schema", "version", "classification", "llm_policy", "data",
        } or value["schema"] != PARSED_PARTITION_SCHEMA or value["version"] != 2
                or value["classification"] != "controlled" or value["llm_policy"] != "no-llm"
                or set(value["data"]) != {
                    "partition_ordinal", "partition_sha256", "parent_output_sha256",
                    "records", "shards",
                }):
            raise ValueError
        data = value["data"]
        data["shards"] = tuple(_artifact_value(row) for row in data["shards"])
        return ParsedPartitionManifest(**data)
    except WgsAnnotationError:
        raise
    except (ValueError, TypeError, KeyError, json.JSONDecodeError):
        raise WgsAnnotationError("invalid parsed annotation partition") from None


def _load_partition_closure(content: bytes) -> AnnotationPartitionClosure:
    try:
        value = json.loads(content)
        if (
            canonical(value) != content
            or set(value) != {
                "schema", "version", "classification", "llm_policy", "data",
            }
            or value["schema"] != ANNOTATION_CLOSURE_SCHEMA
            or value["version"] not in (2, 3, 4)
            or value["classification"] != "controlled"
            or value["llm_policy"] != "no-llm"
            or set(value["data"]) != {
                "partition_ordinal", "partition_sha256", "parent_plan_sha256",
                "resource_session_sha256", "bundle_manifest_sha256",
                "execution_limits_sha256", "output", "stderr", "child",
                "verified_responses", "parsed_manifest",
            } | ({'unassessed_requests', 'reference_coverage'} if value['version'] >= 3 else set())
              | ({'unassessed_alias_requests'} if value['version'] == 4 else set())
        ):
            raise ValueError
        data = value["data"]
        child = BoundedChildResult(**data.pop("child"))
        output = _artifact_value(data.pop("output"))
        stderr_value = data.pop("stderr")
        stderr = None if stderr_value is None else _artifact_value(stderr_value)
        parsed_manifest = _artifact_value(data.pop("parsed_manifest"))
        if value['version'] >= 3:
            if not _count(data['unassessed_requests'], 1):
                raise ValueError
            data['reference_coverage'] = _artifact_value(data['reference_coverage'])
        if value['version'] == 4 and not _count(data['unassessed_alias_requests'], 1):
            raise ValueError
        return AnnotationPartitionClosure(
            **data, output=output, stderr=stderr, child=child,
            parsed_manifest=parsed_manifest,
        )
    except WgsAnnotationError:
        raise
    except (KeyError, TypeError, ValueError, json.JSONDecodeError, LedgerError):
        raise WgsAnnotationError("invalid annotation partition closure") from None


def _load_stage_manifest(content: bytes) -> AnnotationStageManifest:
    try:
        value = json.loads(content)
        if (
            canonical(value) != content
            or set(value) != {
                "schema", "version", "classification", "llm_policy", "data",
            }
            or value["schema"] != ANNOTATION_STAGE_SCHEMA
            or value["version"] not in (2, 3)
            or value["classification"] != "controlled"
            or value["llm_policy"] != "no-llm"
            or set(value["data"]) != {
                "parent_plan_sha256", "resource_session_sha256",
                "bundle_manifest_sha256", "execution_limits_sha256", "parents",
                "sequence_partitions", "sequence_requests", "partition_closures",
                "parsed_manifests", "state_artifact",
            } | ({"prefix_source"} if value["version"] == 3 else set())
        ):
            raise ValueError
        data = value["data"]
        data["parents"] = tuple(ParentIdentity(**item) for item in data["parents"])
        data["partition_closures"] = tuple(
            _artifact_value(item) for item in data["partition_closures"]
        )
        data["parsed_manifests"] = tuple(
            _artifact_value(item) for item in data["parsed_manifests"]
        )
        data["state_artifact"] = _artifact_value(data["state_artifact"])
        if value["version"] == 3:
            data["prefix_source"] = _artifact_value(data["prefix_source"])
        return AnnotationStageManifest(**data)
    except WgsAnnotationError:
        raise
    except (KeyError, TypeError, ValueError, json.JSONDecodeError, LedgerError):
        raise WgsAnnotationError("invalid annotation stage manifest") from None


def _completed_child_artifact(path: Path, *, expected_sha256: str,
                              expected_size: int, maximum_bytes: int) -> Artifact | None:
    """Use the successful writer's identity after stable private metadata checks.

    No payload read occurs here. Stdout's necessary reconciliation authenticates
    its bytes before completion; unread stderr retains the writer's identity.
    Failure/recovery readers without successful writer facts use _existing_artifact.
    """
    if (not _digest(expected_sha256) or not _count(expected_size)
            or not _count(maximum_bytes) or maximum_bytes < 1
            or expected_size > maximum_bytes
            or expected_size == 0 and expected_sha256 != hashlib.sha256(b"").hexdigest()):
        raise WgsAnnotationError("annotation child artifact binding is invalid")
    try:
        if not stat.S_ISREG(path.lstat().st_mode):
            raise OSError
    except OSError:
        raise WgsAnnotationError("annotation child artifact is unsafe") from None
    stream, initial = _open_stable_response(path)
    try:
        final = os.fstat(stream.fileno())
        named = os.stat(path, follow_symlinks=False)
        def identity(value):
            return (value.st_dev, value.st_ino, value.st_size, value.st_mtime_ns,
                    value.st_ctime_ns, value.st_mode, value.st_uid, value.st_nlink)
        if initial.st_size != expected_size or identity(initial) != identity(final) or identity(final) != identity(named):
            raise WgsAnnotationError("annotation child artifact changed")
        return None if expected_size == 0 else Artifact(path.name, expected_sha256, expected_size)
    except OSError:
        raise WgsAnnotationError("annotation child artifact changed") from None
    finally:
        stream.close()


def _existing_artifact(path: Path, *, maximum_bytes: int) -> Artifact | None:
    if not path.exists():
        return None
    stream, initial = _open_stable_response(path)
    try:
        if initial.st_size == 0:
            return None
        if initial.st_size > maximum_bytes:
            raise WgsAnnotationError("annotation attempt artifact exceeds capacity")
        digest = hashlib.sha256()
        size = 0
        while block := stream.read(1024 * 1024):
            size += len(block)
            if size > maximum_bytes:
                raise WgsAnnotationError("annotation attempt artifact exceeds capacity")
            digest.update(block)
        final = os.fstat(stream.fileno())
        named = os.stat(path, follow_symlinks=False)
        if (
            size != initial.st_size
            or (initial.st_dev, initial.st_ino, initial.st_size,
                initial.st_mtime_ns, initial.st_ctime_ns)
            != (final.st_dev, final.st_ino, final.st_size,
                final.st_mtime_ns, final.st_ctime_ns)
            or (final.st_dev, final.st_ino) != (named.st_dev, named.st_ino)
        ):
            raise WgsAnnotationError("annotation attempt artifact changed")
        return Artifact(path.name, digest.hexdigest(), size)
    finally:
        stream.close()


def _verify_existing(root: Path, artifact: Artifact, *, maximum_bytes: int) -> None:
    observed = _existing_artifact(root / artifact.name, maximum_bytes=maximum_bytes)
    if observed != artifact:
        raise WgsAnnotationError("annotation attempt artifact is invalid")


def _annotation_parents(
    plan: VepTransportPlan,
    *,
    execution_limits_sha256: str,
    execution_limits_version: int,
    resource_session_sha256: str,
    bundle_manifest_sha256: str,
    executable_stack_sha256: str,
    runtime_image_digest: str,
    prefix_sha256: str | None = None,
    prefix_version: int = 1,
) -> tuple[ParentIdentity, ...]:
    runtime_sha256 = runtime_image_digest.removeprefix("sha256:")
    if (
        not _digest(execution_limits_sha256)
        or not _digest(resource_session_sha256)
        or not _digest(bundle_manifest_sha256)
        or not _digest(executable_stack_sha256)
        or not re.fullmatch(r"sha256:[0-9a-f]{64}", runtime_image_digest)
        or prefix_sha256 is not None and not _digest(prefix_sha256)
        or type(prefix_version) is not int or prefix_version not in (1, 2)
    ):
        raise WgsAnnotationError("annotation journal binding is invalid")
    return (
        ParentIdentity("stage_parent", VEP_PLAN_SCHEMA, 1, plan.sha256),
        ParentIdentity("partition_plan", VEP_PLAN_SCHEMA, 1, plan.sha256),
        ParentIdentity(
            "execution_limits", "genosapiens.bounded-wgs-execution",
            execution_limits_version,
            execution_limits_sha256,
        ),
        ParentIdentity(
            "vep_resource_session", "genosapiens.vep-resource-session", 1,
            resource_session_sha256,
        ),
        ParentIdentity(
            "vep_bundle_manifest", "genosapiens.vep-bundle-manifest", 1,
            bundle_manifest_sha256,
        ),
        ParentIdentity(
            "executable_stack", "genosapiens.executable-stack", 1,
            executable_stack_sha256,
        ),
        ParentIdentity(
            "runtime_image", "genosapiens.runtime-image", 1, runtime_sha256,
        ),
    ) + (() if prefix_sha256 is None else (ParentIdentity("completed_annotation_prefix",
        "genosapiens.wgs-completed-annotation-prefix", prefix_version, prefix_sha256),))


def _attempt_paths(root: Path, partition: VepTransportPartition, attempt: int):
    return (
        root / f"annotated-{partition.ordinal:06d}-attempt-{attempt:02d}.vcf",
        root / f"vep-{partition.ordinal:06d}-attempt-{attempt:02d}.stderr",
        root / f"annotation-{partition.ordinal:06d}-attempt-{attempt:02d}.json",
    )


def _attempt_evidence(
    root: Path,
    partition: VepTransportPartition,
    attempt: int,
    config: BoundedWgsExecutionConfig,
) -> tuple[Artifact, ...]:
    output, stderr, closure = _attempt_paths(root, partition, attempt)
    evidence = []
    for path, maximum in (
        (output, config.child_process.max_total_output_bytes),
        (stderr, config.child_process.max_stderr_bytes),
        (closure, MAX_MANIFEST_BYTES),
    ):
        artifact = _existing_artifact(path, maximum_bytes=maximum)
        if artifact is not None:
            evidence.append(artifact)
    return tuple(evidence)


def _read_partition_closure_metadata(
    root: Path,
    closure_artifact: Artifact,
    partition: VepTransportPartition,
    *,
    plan: VepTransportPlan,
    plan_root: Path,
    resource_session_sha256: str,
    bundle_manifest_sha256: str,
    execution_config: BoundedWgsExecutionConfig,
) -> tuple[AnnotationPartitionClosure, ParsedPartitionManifest]:
    content = read_artifact(root, closure_artifact)
    closure = _load_partition_closure(content)
    limits_sha256 = hashlib.sha256(
        bounded_wgs_execution_config_as_json(execution_config)
    ).hexdigest()
    if (
        closure.partition_ordinal != partition.ordinal
        or closure.partition_sha256 != partition.sha256
        or closure.parent_plan_sha256 != plan.sha256
        or closure.resource_session_sha256 != resource_session_sha256
        or closure.bundle_manifest_sha256 != bundle_manifest_sha256
        or closure.execution_limits_sha256 != limits_sha256
        or closure.accounted_requests != partition.requests
    ):
        raise WgsAnnotationError("annotation partition closure binding is invalid")
    parsed = _load_parsed_partition(read_artifact(root / "parsed", closure.parsed_manifest))
    if (parsed.partition_ordinal != partition.ordinal
            or parsed.partition_sha256 != partition.sha256
            or parsed.parent_output_sha256 != closure.output.sha256
            or parsed.records != closure.accounted_requests):
        raise WgsAnnotationError("parsed annotation partition binding is invalid")
    if closure.reference_coverage is not None:
        _closure_reference_coverage(root, closure)
    return closure, parsed


def _closure_reference_coverage(root, closure):
    artifact = closure.reference_coverage
    if artifact is None:
        return None
    try:
        if artifact.name != f'reference-coverage-{artifact.sha256}.fai':
            raise ValueError
        return ReferenceCoverage(read_artifact(root, artifact), artifact.sha256)
    except (ValueError, LedgerError):
        raise WgsAnnotationError('annotation reference coverage proof is invalid') from None


def _verify_partition_closure(
    root: Path,
    closure_artifact: Artifact,
    partition: VepTransportPartition,
    *,
    plan: VepTransportPlan,
    plan_root: Path,
    resource_session_sha256: str,
    bundle_manifest_sha256: str,
    execution_config: BoundedWgsExecutionConfig,
) -> AnnotationPartitionClosure:
    closure, parsed = _read_partition_closure_metadata(root, closure_artifact, partition,
        plan=plan, plan_root=plan_root, resource_session_sha256=resource_session_sha256,
        bundle_manifest_sha256=bundle_manifest_sha256, execution_config=execution_config)
    limits_sha256 = hashlib.sha256(bounded_wgs_execution_config_as_json(execution_config)).hexdigest()
    binding = _annotation_state_binding(plan.sha256, resource_session_sha256,
                                        bundle_manifest_sha256, limits_sha256)
    with closing(_open_annotation_state(root, binding, execution_config, writable=False)) as state:
        attempt_match = re.fullmatch(
            rf"annotation-{partition.ordinal:06d}-attempt-([0-9]{{2}})\.json",
            closure_artifact.name,
        )
        if attempt_match is None:
            raise WgsAnnotationError("annotation partition closure name is invalid")
        commit = state.execute(
            "SELECT cursor_bytes,outputs_bytes FROM commits WHERE partition_id=?",
            (f"partition-{partition.ordinal:06d}",),
        ).fetchone()
        counts = tuple(
            state.execute(f"SELECT count(*) FROM {table} WHERE partition_ordinal=?",
                (partition.ordinal,)).fetchone()[0] for table in ("locations", "expected", "responses"))
        if (commit is None
                or json.loads(commit[0]) != [
                    partition.ordinal,
                    int(attempt_match.group(1)),
                    closure.accounted_requests,
                ]
                or json.loads(commit[1])[-1:] != [closure_artifact.name]
                or counts != (parsed.records, parsed.records, closure.verified_responses)
                or closure.unassessed_requests and state.execute(
                    'SELECT count(*) FROM unassessed WHERE partition_ordinal=?',
                    (partition.ordinal,)).fetchone()[0] != closure.unassessed_requests):
            raise WgsAnnotationError("annotation partition checkpoint is inconsistent")
    return closure



def read_completed_partition_closure(
    root: Path,
    closure_artifact: Artifact,
    partition: VepTransportPartition,
    *,
    plan: VepTransportPlan,
    plan_root: Path,
    resource_session_sha256: str,
    bundle_manifest_sha256: str,
    execution_config: BoundedWgsExecutionConfig,
) -> AnnotationPartitionClosure:
    """Reuse a producer partition commitment; validate requested rows in readers."""
    closure, parsed = _read_partition_closure_metadata(root, closure_artifact, partition,
        plan=plan, plan_root=plan_root, resource_session_sha256=resource_session_sha256,
        bundle_manifest_sha256=bundle_manifest_sha256, execution_config=execution_config)
    limits_sha256 = hashlib.sha256(bounded_wgs_execution_config_as_json(execution_config)).hexdigest()
    binding = _annotation_state_binding(plan.sha256, resource_session_sha256, bundle_manifest_sha256, limits_sha256)
    with closing(_open_annotation_state(root, binding, execution_config, writable=False)) as state:
        commit = state.execute('SELECT cursor_bytes,outputs_bytes,counts_bytes FROM commits WHERE partition_id=?',
            (f'partition-{partition.ordinal:06d}',)).fetchone()
        attempt = re.fullmatch(rf'annotation-{partition.ordinal:06d}-attempt-([0-9]{{2}})\.json', closure_artifact.name)
        if (commit is None or attempt is None
                or json.loads(commit[0]) != [partition.ordinal,int(attempt.group(1)),closure.accounted_requests]
                or json.loads(commit[1])[-1:] != [closure_artifact.name]
                or json.loads(commit[2]) != closure.counts):
            raise WgsAnnotationError('completed annotation partition metadata differs')
    return closure


def _annotation_state_binding(plan_sha256: str, resource_sha256: str,
                              bundle_sha256: str, limits_sha256: str) -> str:
    return hashlib.sha256(canonical({"schema": "genosapiens.wgs-annotation-state",
        "version": 3, "plan_sha256": plan_sha256, "resource_session_sha256": resource_sha256,
        "bundle_manifest_sha256": bundle_sha256, "execution_limits_sha256": limits_sha256})).hexdigest()


def _verify_state_descriptor(root: Path, artifact: Artifact) -> None:
    try:
        metadata = (root.parent / "state" / artifact.name).lstat()
        if (not stat.S_ISREG(metadata.st_mode) or metadata.st_uid != os.getuid()
                or stat.S_IMODE(metadata.st_mode) != 0o600 or metadata.st_nlink != 1
                or metadata.st_size != artifact.size_bytes):
            raise ValueError
    except (OSError, ValueError):
        raise WgsAnnotationError("annotation state descriptor is invalid") from None


def _open_annotation_state(root: Path, binding: str, config: BoundedWgsExecutionConfig,
                           *, writable: bool, connection: sqlite3.Connection | None = None,
                           created: bool = False) -> sqlite3.Connection:
    path = root.parent / "state" / "annotation.sqlite"
    path.parent.mkdir(mode=0o700, exist_ok=True)
    owned = connection is None
    if writable and owned or created and not writable:
        raise WgsAnnotationError("annotation state owner is invalid")
    try:
        metadata = path.lstat()
        if (not stat.S_ISREG(metadata.st_mode) or metadata.st_uid != os.getuid()
                or stat.S_IMODE(metadata.st_mode) != 0o600 or metadata.st_nlink != 1):
            raise ValueError
        if owned:
            from genosapiens.wgs_read_cache import cache_bounded_readonly_connection
            connection = cache_bounded_readonly_connection(path)
            connection.enable_load_extension(False)
            connection.execute("PRAGMA foreign_keys=ON")
            connection.execute(f"PRAGMA cache_size=-{max(1, config.sqlite.cache_bytes // 1024)}")
            connection.execute("PRAGMA temp_store=FILE")
            connection.execute("PRAGMA mmap_size=0")
            connection.execute("PRAGMA trusted_schema=OFF")
        if created:
            connection.executescript("""
                CREATE TABLE annotation_binding(binding_sha256 TEXT PRIMARY KEY);
                CREATE TABLE locations(request_id TEXT PRIMARY KEY,
                    partition_ordinal INTEGER NOT NULL, request_ordinal INTEGER NOT NULL UNIQUE,
                    shard_name TEXT NOT NULL, row_ordinal INTEGER NOT NULL);
                CREATE TABLE expected(request_id TEXT PRIMARY KEY,
                    partition_ordinal INTEGER NOT NULL, request_ordinal INTEGER NOT NULL UNIQUE,
                    canonical_request BLOB NOT NULL);
                CREATE TABLE responses(request_id TEXT PRIMARY KEY REFERENCES expected(request_id),
                    partition_ordinal INTEGER NOT NULL, raw_line_sha256 TEXT NOT NULL,
                    parent_output_sha256 TEXT NOT NULL);
                CREATE INDEX locations_partition ON locations(partition_ordinal,request_ordinal);
                CREATE INDEX expected_partition ON expected(partition_ordinal);
                CREATE INDEX responses_partition ON responses(partition_ordinal);
            """)
            connection.execute("INSERT INTO annotation_binding VALUES (?)", (binding,))
        if (connection.execute("PRAGMA user_version").fetchone()[0] != 1
                or connection.execute("SELECT binding_sha256 FROM annotation_binding").fetchone() != (binding,)):
            raise ValueError
        return connection
    except (OSError, ValueError, sqlite3.Error):
        if owned and connection is not None:
            connection.close()
        raise WgsAnnotationError("annotation state is missing or incompatible") from None


def _publish_parsed_partition(records: Iterable[EncodedVerifiedVepRecord], *, root: Path,
                              partition: VepTransportPartition, output_sha256: str,
                              attempt: int, config: BoundedWgsExecutionConfig,
                              state: sqlite3.Connection, outcome_counts=None) -> tuple[Artifact, int]:
    parsed_root = root / "parsed"
    parsed_root.mkdir(mode=0o700, exist_ok=True)
    shards: list[Artifact] = []
    pending: list[tuple[EncodedVerifiedVepRecord, bytes]] = []
    pending_bytes = 0
    count = 0

    def flush() -> None:
        nonlocal pending_bytes
        if not pending:
            return
        name = (f"partition-{partition.ordinal:06d}-attempt-{attempt:02d}-"
                f"shard-{len(shards) + 1:06d}.jsonl")
        artifact = publish_stream(parsed_root, name, (row for _, row in pending),
                                  maximum_bytes=config.annotation.max_shard_bytes)
        for ordinal, (record, _) in enumerate(pending, 1):
            state.execute("INSERT INTO expected VALUES (?,?,?,?)", (
                record.request_id, partition.ordinal,
                record.request_ordinal, record.canonical_request))
            table = 'unassessed' if type(record) is EncodedUnassessedVepRecord else 'responses'
            state.execute(f"INSERT INTO {table} VALUES (?,?,?,?)", (
                record.request_id, partition.ordinal,
                record.raw_line_sha256, record.parent_output_sha256))
            state.execute("INSERT INTO locations VALUES (?,?,?,?,?)", (
                record.request_id, partition.ordinal,
                record.request_ordinal, artifact.name, ordinal))
        shards.append(artifact)
        pending.clear()
        pending_bytes = 0

    for record in records:
        content = record.parsed_record
        if outcome_counts is not None and type(record) is EncodedUnassessedVepRecord:
            reason = _load_parsed_record(content).reason
            outcome_counts[reason] = outcome_counts.get(reason, 0) + 1
        if (len(content) > MAX_ROW_BYTES or len(content) > config.annotation.max_shard_bytes):
            raise WgsAnnotationError("parsed VEP record exceeds shard capacity")
        if (pending and (len(pending) >= config.annotation.max_shard_rows
                         or pending_bytes + len(content) > config.annotation.max_shard_bytes)):
            flush()
        pending.append((record, content))
        pending_bytes += len(content)
        count += 1
    flush()
    manifest = ParsedPartitionManifest(partition.ordinal, partition.sha256,
                                       output_sha256, count, tuple(shards))
    name = f"partition-{partition.ordinal:06d}-attempt-{attempt:02d}-manifest.json"
    return publish_bytes(parsed_root, name, manifest.as_bytes()), count


def _stage_manifest_artifact(
    root: Path, manifest: AnnotationStageManifest,
) -> Artifact:
    content = manifest.as_bytes()
    expected = Artifact(
        "annotation-manifest.json", hashlib.sha256(content).hexdigest(), len(content),
    )
    path = root / expected.name
    if path.exists():
        if read_artifact(root, expected) != content:
            raise WgsAnnotationError("annotation stage orphan is invalid")
        return expected
    return publish_bytes(root, expected.name, content)


def run_annotation_partitions(
    *,
    plan_root: str | Path,
    plan: VepTransportPlan,
    source_manifest: LedgerManifest,
    working_manifest: WorkingSetManifest,
    output_root: str | Path,
    cache_directory: str | Path,
    fasta_path: str | Path,
    resource_session_sha256: str,
    bundle_manifest_sha256: str,
    executable_stack_sha256: str,
    runtime_image_digest: str,
    execution_config: BoundedWgsExecutionConfig,
    intake_binding: AnnotationIntakeBinding | None = None,
    prepared_index: PreparedIndex | None = None,
    prepared_state_root: str | Path | None = None,
    resume: bool = False,
    run_options_sha256: str | None = None,
    phase_recorder: RuntimePhaseRecorder | None = None,
    completed_prefix=None,
    reference_coverage: ReferenceCoverage | None = None,
) -> AnnotationStageResult:
    """Run/resume a fixed annotation plan with durable scientific attempts."""

    plan_root, root = Path(plan_root), Path(output_root)
    if reference_coverage is not None and type(reference_coverage) is not ReferenceCoverage:
        raise WgsAnnotationError('annotation reference coverage owner is invalid')
    if completed_prefix is not None:
        from genosapiens.wgs_annotation_prefix import ResolvedAnnotationPrefix
        if (type(completed_prefix) is not ResolvedAnnotationPrefix or resume
                or completed_prefix.descriptor.version != 2
                or completed_prefix.count != len(plan.partitions)
                or completed_prefix.descriptor.plan_sha256 != plan.sha256
                or completed_prefix.producer.runtime_image_digest != runtime_image_digest
                or completed_prefix.producer.binding("vep_bundle") != bundle_manifest_sha256):
            raise WgsAnnotationError("annotation prefix execution is incompatible")
    stage_state_path = root.parent / "state" / "annotation.sqlite"
    stage_head_path = root / "annotation-head.json"
    fresh_stage = not (stage_state_path.exists() or stage_state_path.is_symlink())
    if completed_prefix is not None:
        load_completed_vep_plan(output_root=plan_root, expected_manifest_sha256=plan.sha256,
            source_manifest=source_manifest, working_manifest=working_manifest,
            execution_config=execution_config, intake_binding=intake_binding,
            prepared_index=prepared_index, prepared_state_root=prepared_state_root)
    elif not resume or fresh_stage and not (
            stage_head_path.exists() or stage_head_path.is_symlink()):
        verify_vep_partition_plan(
            root=plan_root,
            plan=plan,
            source_manifest=source_manifest,
            working_manifest=working_manifest,
            execution_config=execution_config,
            intake_binding=intake_binding,
            prepared_index=prepared_index,
            prepared_state_root=prepared_state_root,
        )
    if not root.is_absolute():
        raise WgsAnnotationError("annotation stage root is invalid")
    limits_sha256 = hashlib.sha256(
        bounded_wgs_execution_config_as_json(execution_config)
    ).hexdigest()
    parents = _annotation_parents(
        plan,
        execution_limits_sha256=limits_sha256,
        execution_limits_version=execution_config.version,
        resource_session_sha256=resource_session_sha256,
        bundle_manifest_sha256=bundle_manifest_sha256,
        executable_stack_sha256=executable_stack_sha256,
        runtime_image_digest=runtime_image_digest,
        prefix_sha256=None if completed_prefix is None else completed_prefix.descriptor.sha256,
        prefix_version=1 if completed_prefix is None else completed_prefix.descriptor.version,
    )
    return _run_annotation_partitions_committed(
        plan_root=plan_root, plan=plan,
        source_manifest=source_manifest, working_manifest=working_manifest,
        root=root, cache_directory=cache_directory, fasta_path=fasta_path,
        resource_session_sha256=resource_session_sha256,
        bundle_manifest_sha256=bundle_manifest_sha256,
        executable_stack_sha256=executable_stack_sha256,
        runtime_image_digest=runtime_image_digest,
        execution_config=execution_config, limits_sha256=limits_sha256,
        parents=parents, resume=resume,
        run_options_sha256=run_options_sha256,
        phase_recorder=phase_recorder,
        completed_prefix=completed_prefix,
        reference_coverage=reference_coverage,
    )


def _annotation_attempt_limit(stage_state, partition_id, *, default, resume, coverage_repair=False):
    """Permit one new attempt after an eligible failure in explicit resume.

    Retain every prior start and terminal event. Completed partitions bypass
    this function. A coverage repair permits one old missing-response retry.
    Other response failures keep their configured attempt allowance.
    """
    if not resume:
        return default
    saved = stage_state.connection.execute(
        "SELECT event_bytes,event_sha256 FROM events WHERE partition_id=? "
        "ORDER BY ordinal DESC LIMIT 1", (partition_id,),
    ).fetchone()
    if saved is None:
        return default
    terminal = load_model(saved[0], JournalEntry)
    if (terminal.sha256 != saved[1] or terminal.partition_id != partition_id
            or terminal.stage != stage_state.binding.stage
            or terminal.stage_parent_sha256 != stage_state.binding.stage_parent_sha256
            or terminal.parents != stage_state.binding.parents):
        raise WgsAnnotationError("annotation retry journal binding is invalid")
    eligible = (
        terminal.event_type == "attempt_interrupted" and terminal.failure_code == "process_lost"
        or terminal.event_type == "attempt_failed" and terminal.failure_code in {
            "annotation_child_storage", "annotation_child_unexpected",
        }
        or coverage_repair and terminal.event_type == 'attempt_failed'
        and terminal.failure_code == 'annotation_child_missing_response'
    )
    if eligible and terminal.scientific and terminal.attempt_number < 99:
        return max(default, terminal.attempt_number + 1)
    return default


def _run_annotation_partitions_committed(
    *, plan_root, plan, source_manifest, working_manifest, root,
    cache_directory, fasta_path, resource_session_sha256,
    bundle_manifest_sha256, executable_stack_sha256,
    runtime_image_digest, execution_config, limits_sha256, parents,
    resume, run_options_sha256, phase_recorder, completed_prefix=None,
    reference_coverage=None,
) -> AnnotationStageResult:
    if (type(resume) is not bool
            or run_options_sha256 is not None and not _digest(run_options_sha256)
            or phase_recorder is not None
            and type(phase_recorder) is not RuntimePhaseRecorder):
        raise WgsAnnotationError("annotation resume binding is invalid")
    options_sha256 = run_options_sha256 or hashlib.sha256(canonical({
        "schema": "genosapiens.wgs-annotation-run-options", "version": 1,
        "mode": "default",
    })).hexdigest()
    binding = StageRunBinding(
        "annotation", plan.sha256, executable_stack_sha256,
        runtime_image_digest, limits_sha256, options_sha256,
        source_manifest.source.input_sha256,
        source_manifest.source.reference_manifest_sha256,
        working_manifest.parent_filter_manifest_sha256,
        (source_manifest.sha256, working_manifest.sha256,
         resource_session_sha256, bundle_manifest_sha256),
        parents,
    )
    state_binding = _annotation_state_binding(
        plan.sha256, resource_session_sha256,
        bundle_manifest_sha256, limits_sha256,
    )
    prefix_artifact = None if completed_prefix is None else publish_bytes(
        root, "annotation-prefix.json", completed_prefix.descriptor.as_bytes())
    state_path = root.parent / "state" / "annotation.sqlite"
    state_path.parent.mkdir(mode=0o700, exist_ok=True)
    state_exists = state_path.exists() or state_path.is_symlink()
    if state_exists and not resume:
        raise WgsAnnotationError("annotation checkpoint requires explicit resume")
    if resume and not state_exists and (
            (root / "annotation-head.json").exists()
            or (root / "annotation-head.json").is_symlink()):
        raise WgsAnnotationError("annotation checkpoint state is missing")

    def materialize_manifest(stage_state, completion):
        seal = completion.artifacts[0]
        seal_content = read_artifact(root, seal)
        seal_value = json.loads(seal_content)
        if (seal_value.get("schema") != "genosapiens.wgs-annotation-stage-closure"
                or seal_value.get("version") != 1
                or seal_value.get("plan_sha256") != plan.sha256
                or type(seal_value.get("partition_closures")) is not list
                or len(seal_value["partition_closures"]) != len(plan.partitions)
                or canonical(seal_value) != seal_content):
            raise WgsAnnotationError("annotation stage closure is invalid")
        path = root / "annotation-manifest.json"
        if path.exists() or path.is_symlink():
            from genosapiens.wgs_shards import read_bounded_metadata
            content = read_bounded_metadata(path)
            manifest = _load_stage_manifest(content)
            artifact = Artifact(path.name, hashlib.sha256(content).hexdigest(), len(content))
            if (read_artifact(root, artifact) != content
                    or manifest.parent_plan_sha256 != plan.sha256
                    or manifest.parents != parents
                    or manifest.resource_session_sha256 != resource_session_sha256
                    or manifest.bundle_manifest_sha256 != bundle_manifest_sha256
                    or manifest.execution_limits_sha256 != limits_sha256
                    or manifest.sequence_partitions != len(plan.partitions)
                    or manifest.sequence_requests != plan.sequence_requests
                    or seal_value["partition_closures"] != [
                        item.sha256 for item in manifest.partition_closures
                    ]):
                raise WgsAnnotationError("annotation completed manifest changed")
            _verify_state_descriptor(root, manifest.state_artifact)
            return manifest, artifact
        closures = [] if completed_prefix is None else list(completed_prefix.closures)
        parsed_manifests = [] if completed_prefix is None else list(completed_prefix.parsed_manifests)
        for partition in plan.partitions if completed_prefix is None else ():
            committed = stage_state.committed(f"partition-{partition.ordinal:06d}")
            if committed is None:
                raise WgsAnnotationError("annotation committed partition is missing")
            artifact = committed.outputs[-1]
            closure = _load_partition_closure(read_artifact(root, artifact))
            if (closure.partition_ordinal != partition.ordinal
                    or closure.partition_sha256 != partition.sha256
                    or closure.accounted_requests != partition.requests):
                raise WgsAnnotationError("annotation committed closure is inconsistent")
            closures.append(artifact)
            parsed_manifests.append(closure.parsed_manifest)
        if seal_value["partition_closures"] != [item.sha256 for item in closures]:
            raise WgsAnnotationError("annotation stage closure changed")
        state_artifact = _existing_artifact(state_path, maximum_bytes=2**63 - 1)
        if state_artifact is None:
            raise WgsAnnotationError("annotation completed state is missing")
        manifest = AnnotationStageManifest(
            plan.sha256, resource_session_sha256, bundle_manifest_sha256,
            limits_sha256, parents, len(plan.partitions), plan.sequence_requests,
            tuple(closures), tuple(parsed_manifests), state_artifact,
            prefix_artifact,
        )
        return manifest, _stage_manifest_artifact(root, manifest)

    try:
        with stage_lock(root, "annotation") as lock, StageCommitState(
                root, state_path, binding, lock=lock, create=not state_exists,
        ) as stage_state, closing(_AnnotationParsePool(execution_config)) as parse_pool:
            state = _open_annotation_state(
                root, state_binding, execution_config, writable=True,
                connection=stage_state.connection, created=not state_exists,
            )
            tail = state.execute(
                "SELECT event_bytes FROM events ORDER BY ordinal DESC LIMIT 1"
            ).fetchone()
            if tail is not None:
                last_event = load_model(tail[0], JournalEntry)
                if last_event.event_type == "stage_completed":
                    manifest, artifact = materialize_manifest(stage_state, last_event)
                    return AnnotationStageResult(
                        manifest, artifact, last_event, len(plan.partitions), completed_prefix,
                    )
            # Additive tables preserve every old row, binding, and commitment.
            # A completed stage returns above without changing its database.
            state.executescript('''
                CREATE TABLE IF NOT EXISTS unassessed(
                    request_id TEXT PRIMARY KEY REFERENCES expected(request_id),
                    partition_ordinal INTEGER NOT NULL, raw_line_sha256 TEXT CHECK(raw_line_sha256 IS NULL),
                    parent_output_sha256 TEXT NOT NULL);
                CREATE INDEX IF NOT EXISTS unassessed_partition ON unassessed(partition_ordinal);
                CREATE TABLE IF NOT EXISTS coverage_attempts(
                    partition_id TEXT PRIMARY KEY, attempt_number INTEGER NOT NULL,
                    reference_fai_sha256 TEXT NOT NULL, policy TEXT NOT NULL);
                CREATE TABLE IF NOT EXISTS alias_coverage_attempts(
                    partition_id TEXT PRIMARY KEY, attempt_number INTEGER NOT NULL,
                    reference_fai_sha256 TEXT NOT NULL, policy TEXT NOT NULL);
            ''')
            closures = [] if completed_prefix is None else list(completed_prefix.closures)
            parsed_manifests = [] if completed_prefix is None else list(completed_prefix.parsed_manifests)
            reused = 0 if completed_prefix is None else len(plan.partitions)
            for partition in plan.partitions if completed_prefix is None else ():
                partition_id = f"partition-{partition.ordinal:06d}"
                completed = stage_state.committed(partition_id)
                if completed is not None:
                    artifact = completed.outputs[-1]
                    closure = _verify_partition_closure(
                        root, artifact, partition, plan=plan, plan_root=plan_root,
                        resource_session_sha256=resource_session_sha256,
                        bundle_manifest_sha256=bundle_manifest_sha256,
                        execution_config=execution_config,
                    )
                    closures.append(artifact)
                    parsed_manifests.append(closure.parsed_manifest)
                    reused += 1
                    continue
                incomplete = stage_state.open_attempt(partition_id)
                if incomplete is not None:
                    output, stderr, closure_path = _attempt_paths(
                        root, partition, incomplete.attempt_number,
                    )
                    stage_state.archive_orphans(
                        names=(output.name, stderr.name, closure_path.name),
                        attempt_id=incomplete.attempt_id,
                    )
                    stage_state.fail(
                        incomplete, failure_code="process_lost", interrupted=True,
                    )
                try:
                    coverage_repair = False
                    if resume and reference_coverage is not None:
                        marked = state.execute('SELECT 1 FROM coverage_attempts WHERE partition_id=?',
                            (partition_id,)).fetchone()
                        alias_marked = state.execute('SELECT 1 FROM alias_coverage_attempts WHERE partition_id=?',
                            (partition_id,)).fetchone()
                        previous = state.execute('SELECT event_bytes FROM events WHERE partition_id=? '
                            'ORDER BY ordinal DESC LIMIT 1', (partition_id,)).fetchone()
                        if ((marked is None or alias_marked is None) and previous is not None
                                and load_model(previous[0], JournalEntry).failure_code == 'annotation_child_missing_response'):
                            coverage_repair = any(
                                reason == UNASSESSED_REASON and marked is None
                                or reason == ALIAS_REASON and alias_marked is None
                                for request in iter_partition_requests(plan_root, partition,
                                    execution_config=execution_config)
                                for reason in (reference_coverage.unassessed_reason(request.record.contig),))
                    start = stage_state.begin(
                        partition_id=partition_id, scientific=True,
                        max_attempts=_annotation_attempt_limit(stage_state, partition_id,
                            default=execution_config.journal.default_scientific_attempts_per_partition,
                            resume=resume, coverage_repair=coverage_repair),
                    )
                except LedgerError:
                    raise WgsAnnotationError("annotation attempt budget is exhausted") from None
                output_path, stderr_path, closure_path = _attempt_paths(
                    root, partition, start.attempt_number,
                )
                try:
                    if reference_coverage is not None:
                        state.execute('INSERT OR IGNORE INTO coverage_attempts VALUES (?,?,?,?)',
                            (partition_id, start.attempt_number, reference_coverage.sha256, COVERAGE_POLICY))
                        state.execute('INSERT OR IGNORE INTO alias_coverage_attempts VALUES (?,?,?,?)',
                            (partition_id, start.attempt_number, reference_coverage.sha256, ALIAS_POLICY))
                        state.commit()
                    if parse_pool.outstanding:
                        raise WgsAnnotationError("annotation parsing overlaps a VEP child")
                    child = run_vep_partition_child(
                        partition_root=plan_root, partition=partition,
                        output_root=root, cache_directory=cache_directory,
                        fasta_path=fasta_path,
                        resource_session_sha256=resource_session_sha256,
                        execution_config=execution_config,
                        attempt_number=start.attempt_number,
                    )
                    output = _completed_child_artifact(
                        output_path,
                        expected_sha256=child.stdout_sha256, expected_size=child.stdout_bytes,
                        maximum_bytes=execution_config.child_process.max_total_output_bytes,
                    )
                    stderr = _completed_child_artifact(
                        stderr_path,
                        expected_sha256=child.stderr_sha256, expected_size=child.stderr_bytes,
                        maximum_bytes=execution_config.child_process.max_stderr_bytes,
                    )
                    if output is None:
                        raise WgsAnnotationError("annotation child output binding is invalid")
                    state.execute("BEGIN IMMEDIATE")
                    outcome_counts = {}
                    with runtime_phase(phase_recorder, "annotation_response_reconciliation"):
                        parsed_manifest, verified = _publish_parsed_partition(
                            iter_verified_vep_records(
                                iter_partition_requests(
                                    plan_root, partition, execution_config=execution_config,
                                ),
                                response_path=output_path,
                                expected_output_sha256=output.sha256,
                                bundle_manifest_sha256=bundle_manifest_sha256,
                                release="116",
                                scratch_path=root / (
                                    f"responses-{partition.ordinal:06d}-"
                                    f"attempt-{start.attempt_number:02d}.sqlite"
                                ),
                                execution_config=execution_config,
                                _parse_pool=parse_pool,
                                reference_coverage=reference_coverage,
                            ),
                            root=root, partition=partition, output_sha256=output.sha256,
                            attempt=start.attempt_number, config=execution_config,
                            state=state,
                            outcome_counts=outcome_counts,
                        )
                    unassessed = state.execute('SELECT count(*) FROM unassessed WHERE partition_ordinal=?',
                        (partition.ordinal,)).fetchone()[0]
                    verified -= unassessed
                    if sum(outcome_counts.values()) != unassessed:
                        raise WgsAnnotationError('unassessed annotation counts differ')
                    coverage_artifact = None
                    if unassessed:
                        coverage_artifact = Artifact(f'reference-coverage-{reference_coverage.sha256}.fai',
                            reference_coverage.sha256, len(reference_coverage.content))
                        if (root / coverage_artifact.name).exists():
                            attest_existing_artifact(root, coverage_artifact)
                        else:
                            publish_bytes(root, coverage_artifact.name, reference_coverage.content)
                    closure = AnnotationPartitionClosure(
                        partition.ordinal, partition.sha256, plan.sha256,
                        resource_session_sha256, bundle_manifest_sha256,
                        limits_sha256, output, stderr, child, verified,
                        parsed_manifest, unassessed, coverage_artifact, outcome_counts.get(ALIAS_REASON, 0),
                    )
                    closure_artifact = publish_bytes(
                        root, closure_path.name, closure.as_bytes(),
                    )
                    stage_state.commit_partition(
                        start,
                        cursor=(partition.ordinal, start.attempt_number, closure.accounted_requests),
                        outputs=(closure_artifact,), counts=closure.counts,
                        prepared_transaction=True,
                    )
                    closures.append(closure_artifact)
                    parsed_manifests.append(parsed_manifest)
                except Exception as error:
                    if state.in_transaction:
                        state.rollback()
                    if state.execute(
                            "SELECT 1 FROM commits WHERE partition_id=?",
                            (partition_id,),
                    ).fetchone() is not None:
                        raise WgsAnnotationError(
                            "annotation completion journal is inconsistent"
                        ) from None
                    try:
                        stage_state.fail(
                            start, failure_code=static_failure(
                                error, "annotation_child",
                                operation="annotation_child_execute",
                            ).failure_code,
                        )
                    except Exception:
                        pass
                    if isinstance(error, WgsAnnotationError):
                        raise
                    raise WgsAnnotationError("annotation partition execution failed") from None
            seal_content = canonical({
                "schema": "genosapiens.wgs-annotation-stage-closure", "version": 1,
                "plan_sha256": plan.sha256,
                "partition_closures": [item.sha256 for item in closures],
            })
            seal_path = root / "annotation-stage-closure.json"
            seal = (
                attest_existing_artifact(root, Artifact(
                    seal_path.name, hashlib.sha256(seal_content).hexdigest(),
                    len(seal_content),
                ))
                if seal_path.exists() or seal_path.is_symlink()
                else publish_bytes(root, seal_path.name, seal_content)
            )
            head = stage_state.finish(
                expected_partitions=tuple(
                    f"partition-{partition.ordinal:06d}" for partition in plan.partitions
                    if completed_prefix is None
                ), output=seal, process_id="process-" + os.urandom(16).hex(),
                closure_check=lambda checked: checked.execute(
                    "SELECT (SELECT count(*) FROM responses) + (SELECT count(*) FROM unassessed)"
                ).fetchone() == (plan.sequence_requests - (0 if completed_prefix is None
                    else completed_prefix.descriptor.completed_requests),),
            )
            manifest, artifact = materialize_manifest(stage_state, head)
            return AnnotationStageResult(manifest, artifact, head, reused, completed_prefix)
    except WgsAnnotationError:
        raise
    except (LedgerError, OSError, ValueError, TypeError, sqlite3.Error):
        raise WgsAnnotationError("annotation checkpoint is missing or incompatible") from None


@dataclass(frozen=True, slots=True, repr=False)
class AnnotationStageEncodedRecord:
    """Authenticated saved row plus its original producer/index location.

    This is a downstream transport, not an annotation publication format.
    The reader authenticates the entire bounded shard before exposing bytes;
    the consuming worker owns strict typed decode and index reconciliation.
    """

    encoded: EncodedVerifiedVepRecord
    producer_root: str
    resource_session_sha256: str
    partition_ordinal: int
    shard_name: str
    row_ordinal: int

    def __post_init__(self):
        if (type(self.encoded) not in (EncodedVerifiedVepRecord, EncodedUnassessedVepRecord)
                or type(self.producer_root) is not str or not self.producer_root
                or not _digest(self.resource_session_sha256)
                or not _count(self.partition_ordinal, 1)
                or not _count(self.row_ordinal, 1)
                or type(self.shard_name) is not str
                or Path(self.shard_name).name != self.shard_name):
            raise WgsAnnotationError("saved annotation row binding is invalid")

    def decode(self) -> VerifiedVepRecord | UnassessedVepRecord:
        return self._decode(_load_parsed_record)

    def _decode_owned(self) -> VerifiedVepRecord | UnassessedVepRecord:
        return self._decode(_load_parsed_record_owned)

    def _decode(self, load_record) -> VerifiedVepRecord | UnassessedVepRecord:
        encoded = self.encoded
        if hashlib.sha256(encoded.parsed_record).hexdigest() != encoded.parsed_record_sha256:
            raise WgsAnnotationError("saved annotation row digest is invalid")
        record = load_record(encoded.parsed_record)
        if ((type(record) is UnassessedVepRecord) != (type(encoded) is EncodedUnassessedVepRecord)
                or record.request.request_id != encoded.request_id
                or record.request.request_ordinal != encoded.request_ordinal
                or record.request.as_bytes() != encoded.canonical_request
                or record.raw_line_sha256 != encoded.raw_line_sha256
                or record.parent_output_sha256 != encoded.parent_output_sha256):
            raise WgsAnnotationError("annotation parsed record order is invalid")
        return record


def iter_annotation_stage_encoded_records(
    result: AnnotationStageResult,
    *,
    plan_root: str | Path,
    plan: VepTransportPlan,
    source_manifest: LedgerManifest,
    working_manifest: WorkingSetManifest,
    output_root: str | Path,
    resource_session_sha256: str,
    bundle_manifest_sha256: str,
    execution_config: BoundedWgsExecutionConfig,
    intake_binding: AnnotationIntakeBinding | None = None,
    prepared_index: PreparedIndex | None = None,
    prepared_state_root: str | Path | None = None,
) -> Iterator[AnnotationStageEncodedRecord]:
    """Authenticate saved shards/indexes in order without rebuilding evidence.

    Producer algorithms, formats and original-prefix ownership are unchanged.
    Each consumer must decode the supplied row and reconcile its exact indexed
    request/raw/output fields via AnnotationStageEncodedRecord.decode().
    """

    plan_root, root = Path(plan_root), Path(output_root)
    if (
        not isinstance(result, AnnotationStageResult)
        or result.manifest.parent_plan_sha256 != plan.sha256
        or result.manifest.resource_session_sha256 != resource_session_sha256
        or result.manifest.bundle_manifest_sha256 != bundle_manifest_sha256
        or result.manifest.sequence_partitions != len(plan.partitions)
        or result.manifest.sequence_requests != plan.sequence_requests
        or len(result.manifest.partition_closures) != len(plan.partitions)
        or len(result.manifest.parsed_manifests) != len(plan.partitions)
    ):
        raise WgsAnnotationError("annotation stage replay binding is invalid")
    load_completed_vep_plan(
        output_root=plan_root,
        expected_manifest_sha256=plan.sha256,
        source_manifest=source_manifest,
        working_manifest=working_manifest,
        execution_config=execution_config,
        intake_binding=intake_binding,
        prepared_index=prepared_index,
        prepared_state_root=prepared_state_root,
    )
    _verify_state_descriptor(root, result.manifest.state_artifact)
    prefix = _result_prefix(result, root, plan, execution_config)
    previous = 0
    for partition, artifact, parsed_artifact in zip(
        plan.partitions,
        result.manifest.partition_closures,
        result.manifest.parsed_manifests,
        strict=True,
    ):
        origin, session, counts = _partition_source(prefix, partition, root, resource_session_sha256)
        closure = read_completed_partition_closure(
            origin,
            artifact,
            partition,
            plan=plan,
            plan_root=plan_root,
            resource_session_sha256=session,
            bundle_manifest_sha256=bundle_manifest_sha256,
            execution_config=execution_config,
        )
        if closure.parsed_manifest != parsed_artifact:
            raise WgsAnnotationError("annotation parsed stage binding is invalid")
        parsed = _load_parsed_partition(read_artifact(origin / "parsed", parsed_artifact))
        partition_count = 0
        binding = _annotation_state_binding(plan.sha256, session,
            bundle_manifest_sha256, hashlib.sha256(
                bounded_wgs_execution_config_as_json(execution_config)).hexdigest())
        state = _open_annotation_state(origin, binding, execution_config, writable=False)
        try:
            coverage = _closure_reference_coverage(origin, closure)
            outcome_table = ('(SELECT * FROM responses UNION ALL SELECT * FROM unassessed)'
                if closure.unassessed_requests else 'responses')
            indexed = state.execute("SELECT l.request_id, l.request_ordinal, l.shard_name, "
                "l.row_ordinal, e.canonical_request, r.raw_line_sha256, r.parent_output_sha256 "
                "FROM locations l JOIN expected e USING(request_id) "
                f"JOIN {outcome_table} r USING(request_id) WHERE l.partition_ordinal=? "
                "AND e.partition_ordinal=l.partition_ordinal "
                "AND r.partition_ordinal=l.partition_ordinal "
                "AND l.request_ordinal BETWEEN ? AND ? "
                "ORDER BY l.request_ordinal",
                (partition.ordinal, partition.first_request_ordinal, partition.last_request_ordinal))
            unassessed_count = 0
            unassessed_alias_count = 0
            for shard in parsed.shards:
                row_ordinal = 0
                for line in iter_artifact_lines(origin / "parsed", shard,
                        maximum_bytes=execution_config.annotation.max_shard_bytes,
                        maximum_line_bytes=MAX_ROW_BYTES):
                    row_ordinal += 1
                    expected = indexed.fetchone()
                    if (expected is None or expected[1] != previous + 1
                            or expected[2:4] != (shard.name, row_ordinal)
                            or expected[6] != closure.output.sha256):
                        raise WgsAnnotationError("annotation parsed record order is invalid")
                    encoded_type = EncodedUnassessedVepRecord if expected[5] is None else EncodedVerifiedVepRecord
                    encoded = encoded_type(
                        expected[0], expected[1], expected[4], line,
                        hashlib.sha256(line).hexdigest(), expected[5], expected[6],
                    )
                    if encoded_type is EncodedUnassessedVepRecord:
                        record = _load_parsed_record(line)
                        if (type(record) is not UnassessedVepRecord or coverage is None
                                or record.reference_fai_sha256 != coverage.sha256
                                or coverage.unassessed_reason(record.request.record.contig) != record.reason):
                            raise WgsAnnotationError('unassessed annotation reference proof differs')
                        unassessed_count += 1
                        unassessed_alias_count += record.reason == ALIAS_REASON
                    previous = encoded.request_ordinal
                    partition_count += 1
                    yield AnnotationStageEncodedRecord(encoded, str(origin), session,
                        partition.ordinal, shard.name, row_ordinal)
            if (indexed.fetchone() is not None or unassessed_count != closure.unassessed_requests
                    or unassessed_alias_count != closure.unassessed_alias_requests):
                raise WgsAnnotationError("annotation parsed index has extra records")
        except LedgerError:
            raise WgsAnnotationError("annotation parsed shard is invalid") from None
        except sqlite3.Error:
            raise WgsAnnotationError("annotation parsed index is invalid") from None
        finally:
            state.close()
        if partition_count != parsed.records:
            raise WgsAnnotationError("annotation parsed partition is incomplete")
    if previous != plan.sequence_requests:
        raise WgsAnnotationError("annotation stage replay is incomplete")


def annotation_coverage_summary(source, execution_config):
    """Read bounded committed metadata, without scanning annotation rows."""
    result, plan, root = source.result, source.plan, source.output_root
    manifest = result.manifest
    prefix = _result_prefix(result, root, plan, execution_config)
    verified = unassessed = unassessed_alias = 0
    references = set()
    for partition, artifact in zip(plan.partitions, manifest.partition_closures, strict=True):
        origin, session, _ = _partition_source(prefix, partition, root, manifest.resource_session_sha256)
        closure = read_completed_partition_closure(origin, artifact, partition,
            plan=plan, plan_root=source.plan_root, resource_session_sha256=session,
            bundle_manifest_sha256=manifest.bundle_manifest_sha256, execution_config=execution_config)
        verified += closure.verified_responses
        unassessed += closure.unassessed_requests
        unassessed_alias += closure.unassessed_alias_requests
        if closure.reference_coverage is not None:
            references.add(closure.reference_coverage.sha256)
    if verified + unassessed != manifest.sequence_requests:
        raise WgsAnnotationError('annotation coverage summary is incomplete')
    return dict(schema='genosapiens.wgs-annotation-coverage', version=1,
        policy=ALIAS_POLICY if unassessed_alias else COVERAGE_POLICY, annotation_manifest_sha256=manifest.sha256,
        sequence_requests=manifest.sequence_requests, verified_responses=verified,
        unassessed_requests=unassessed,
        reason_counts={UNASSESSED_REASON: unassessed - unassessed_alias,
                       **({ALIAS_REASON: unassessed_alias} if unassessed_alias else {})},
        reference_fai_sha256=sorted(references), annotation_root=str(root),
        detail_location='committed parsed annotation shards and admission input records',
        interpretation='Missing annotation does not establish a negative finding.')


def iter_annotation_stage_records(
    result: AnnotationStageResult,
    *,
    plan_root: str | Path,
    plan: VepTransportPlan,
    source_manifest: LedgerManifest,
    working_manifest: WorkingSetManifest,
    output_root: str | Path,
    resource_session_sha256: str,
    bundle_manifest_sha256: str,
    execution_config: BoundedWgsExecutionConfig,
    intake_binding: AnnotationIntakeBinding | None = None,
    prepared_index: PreparedIndex | None = None,
    prepared_state_root: str | Path | None = None,
) -> Iterator[VerifiedVepRecord]:
    """Strict typed consumer of the shared authenticated saved-row reader."""
    with closing(iter_annotation_stage_encoded_records(result, plan_root=plan_root,
            plan=plan, source_manifest=source_manifest, working_manifest=working_manifest,
            output_root=output_root, resource_session_sha256=resource_session_sha256,
            bundle_manifest_sha256=bundle_manifest_sha256, execution_config=execution_config,
            intake_binding=intake_binding, prepared_index=prepared_index,
            prepared_state_root=prepared_state_root)) as records:
        for saved in records:
            yield saved.decode()
