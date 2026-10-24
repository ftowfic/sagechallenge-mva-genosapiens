"""Fail-closed boundary for any local or remote model-facing operation."""

from __future__ import annotations

from dataclasses import asdict, dataclass
from enum import StrEnum
import hashlib
import json
from pathlib import Path
import re
from typing import Callable, Mapping, TypeVar

from genosapiens.bounded_json import BoundedJSONError, load_bounded_json, read_bounded_bytes
from genosapiens.vcf_header import VcfHeaderError, validate_schema_header


EXCEPTION_REGISTRY_SCHEMA = "genosapiens.llm-exception-registry"
MAX_MODEL_ARTIFACT_BYTES = 2 * 1024 * 1024
MAX_MODEL_ARTIFACTS = 1_000
IDENTIFIER = re.compile(r"^[A-Za-z0-9][A-Za-z0-9._:-]*$")
SHA256 = re.compile(r"^[0-9a-f]{64}$")
HPO_TERM = re.compile(r"^HP:\d{7}(?:[ \t]+[^\x00-\x1f\x7f]+)?$")


class LLMGuardError(RuntimeError):
    """Raised without protected content when model dispatch is denied."""

    def __init__(self, code: str) -> None:
        self.code = code
        super().__init__(f"model operation blocked ({code})")


class LLMPolicy(StrEnum):
    NO_LLM = "no-llm"
    ALLOWED = "allowed"
    SANITIZED_HEADER_EXCEPTION = "sanitized-header-exception"
    EXAMPLE_HPO_EXCEPTION = "example-hpo-exception"


class LLMArtifactKind(StrEnum):
    CHALLENGE_HPO = "challenge-hpo"
    RAW_VCF = "raw-vcf"
    RAW_VCF_INDEX = "raw-vcf-index"
    FASTQ = "fastq"
    PATIENT_DERIVED_OUTPUT = "patient-derived-output"
    SANITIZED_VCF_HEADER = "sanitized-vcf-header"
    NON_SENSITIVE_EXAMPLE_HPO = "non-sensitive-example-hpo"
    PUBLIC = "public"
    INDEPENDENTLY_SYNTHETIC = "independently-synthetic"
    NON_SENSITIVE_OPERATIONAL_METADATA = "non-sensitive-operational-metadata"


class ModelSink(StrEnum):
    PROMPT = "prompt"
    CONTEXT = "context"
    TOOL_RESULT = "tool-result"
    SUMMARY = "summary"
    EMBEDDING = "embedding"
    RETRIEVAL_INDEX = "retrieval-index"
    TRANSCRIPT = "transcript"
    REQUEST = "request"


class ModelDestination(StrEnum):
    LOCAL = "local"
    REMOTE = "remote"


class ExceptionKind(StrEnum):
    SANITIZED_VCF_HEADER = "sanitized-vcf-header"
    NON_SENSITIVE_EXAMPLE_HPO = "non-sensitive-example-hpo"


PROTECTED_KINDS = {
    LLMArtifactKind.CHALLENGE_HPO,
    LLMArtifactKind.RAW_VCF,
    LLMArtifactKind.RAW_VCF_INDEX,
    LLMArtifactKind.FASTQ,
    LLMArtifactKind.PATIENT_DERIVED_OUTPUT,
}
ALLOWED_KINDS = {
    LLMArtifactKind.PUBLIC,
    LLMArtifactKind.INDEPENDENTLY_SYNTHETIC,
    LLMArtifactKind.NON_SENSITIVE_OPERATIONAL_METADATA,
}
EXPECTED_EXCEPTION_DESIGNATIONS = {
    ExceptionKind.SANITIZED_VCF_HEADER: "validated-sanitized-vcf-header",
    ExceptionKind.NON_SENSITIVE_EXAMPLE_HPO: "explicitly-non-sensitive-example-hpo",
}


@dataclass(frozen=True, slots=True)
class LLMArtifact:
    artifact_id: str
    path: str
    sha256: str
    kind: LLMArtifactKind
    policy: LLMPolicy | None
    parent_artifact_ids: tuple[str, ...] = ()

    def __post_init__(self) -> None:
        path = Path(self.path)
        if (
            IDENTIFIER.fullmatch(self.artifact_id) is None
            or not path.is_absolute()
            or SHA256.fullmatch(self.sha256) is None
            or self.parent_artifact_ids
            != tuple(sorted(set(self.parent_artifact_ids)))
            or any(
                IDENTIFIER.fullmatch(parent) is None
                for parent in self.parent_artifact_ids
            )
            or self.artifact_id in self.parent_artifact_ids
        ):
            raise LLMGuardError("invalid-artifact-metadata")
        if self.policy is None:
            return
        expected_policy = (
            LLMPolicy.NO_LLM
            if self.kind in PROTECTED_KINDS
            else LLMPolicy.SANITIZED_HEADER_EXCEPTION
            if self.kind is LLMArtifactKind.SANITIZED_VCF_HEADER
            else LLMPolicy.EXAMPLE_HPO_EXCEPTION
            if self.kind is LLMArtifactKind.NON_SENSITIVE_EXAMPLE_HPO
            else LLMPolicy.ALLOWED
            if self.kind in ALLOWED_KINDS
            else None
        )
        if self.policy is not expected_policy:
            raise LLMGuardError("invalid-artifact-classification")


@dataclass(frozen=True, slots=True)
class LLMException:
    artifact_id: str
    path: str
    sha256: str
    kind: ExceptionKind
    designation: str

    def __post_init__(self) -> None:
        if (
            IDENTIFIER.fullmatch(self.artifact_id) is None
            or not Path(self.path).is_absolute()
            or SHA256.fullmatch(self.sha256) is None
            or self.designation != EXPECTED_EXCEPTION_DESIGNATIONS.get(self.kind)
        ):
            raise LLMGuardError("invalid-exception-registry")


@dataclass(frozen=True, slots=True)
class LLMExceptionRegistry:
    exceptions: tuple[LLMException, ...]
    schema: str = EXCEPTION_REGISTRY_SCHEMA

    def __post_init__(self) -> None:
        identities = [item.artifact_id for item in self.exceptions]
        kinds = [item.kind for item in self.exceptions]
        if (
            self.schema != EXCEPTION_REGISTRY_SCHEMA
            or len(identities) != len(set(identities))
            or len(kinds) != len(set(kinds))
            or self.exceptions
            != tuple(sorted(self.exceptions, key=lambda item: item.artifact_id))
        ):
            raise LLMGuardError("invalid-exception-registry")


@dataclass(frozen=True, slots=True)
class GuardAuditEvent:
    allowed: bool
    destination: ModelDestination
    sink: ModelSink
    artifact_count: int
    exception_count: int
    reason: str


Payload = TypeVar("Payload")
Result = TypeVar("Result")
AuditRecorder = Callable[[GuardAuditEvent], None]


def exception_registry_as_bytes(registry: LLMExceptionRegistry) -> bytes:
    return (
        json.dumps(asdict(registry), sort_keys=True, separators=(",", ":")) + "\n"
    ).encode()


def load_exception_registry(path: str | Path) -> LLMExceptionRegistry:
    try:
        data, content = load_bounded_json(path)
        if (
            not isinstance(data, Mapping)
            or set(data) != {"exceptions", "schema"}
            or data["schema"] != EXCEPTION_REGISTRY_SCHEMA
            or not isinstance(data["exceptions"], list)
            or len(data["exceptions"]) > 2
            or any(
                not isinstance(item, Mapping)
                or set(item)
                != {"artifact_id", "designation", "kind", "path", "sha256"}
                for item in data["exceptions"]
            )
        ):
            raise LLMGuardError("invalid-exception-registry")
        registry = LLMExceptionRegistry(
            exceptions=tuple(
                LLMException(
                    artifact_id=item["artifact_id"],
                    path=item["path"],
                    sha256=item["sha256"],
                    kind=ExceptionKind(item["kind"]),
                    designation=item["designation"],
                )
                for item in data["exceptions"]
            ),
            schema=data["schema"],
        )
    except (BoundedJSONError, KeyError, TypeError, ValueError) as error:
        raise LLMGuardError("invalid-exception-registry") from error
    if exception_registry_as_bytes(registry) != content:
        raise LLMGuardError("invalid-exception-registry")
    return registry


def _validate_sanitized_header(content: bytes) -> None:
    try:
        validate_schema_header(content)
    except VcfHeaderError:
        raise LLMGuardError("exception-content-invalid") from None


def _validate_example_hpo(content: bytes) -> None:
    try:
        text = content.decode("utf-8", errors="strict")
    except UnicodeDecodeError as error:
        raise LLMGuardError("exception-content-invalid") from error
    lines = tuple(line.strip() for line in text.splitlines())
    if (
        not text.endswith("\n")
        or not lines
        or any(not line or HPO_TERM.fullmatch(line) is None for line in lines)
    ):
        raise LLMGuardError("exception-content-invalid")


def _record_audit(recorder: AuditRecorder, event: GuardAuditEvent) -> None:
    try:
        recorder(event)
    except Exception as error:
        raise LLMGuardError("audit-failed") from error


def _deny(
    code: str,
    *,
    sink: ModelSink,
    destination: ModelDestination,
    artifact_count: int,
    exception_count: int,
    record_audit: AuditRecorder,
) -> None:
    _record_audit(
        record_audit,
        GuardAuditEvent(
            allowed=False,
            destination=destination,
            sink=sink,
            artifact_count=artifact_count,
            exception_count=exception_count,
            reason=code,
        ),
    )
    raise LLMGuardError(code)


def guarded_model_dispatch(
    *,
    artifact_ids: tuple[str, ...],
    artifacts: Mapping[str, LLMArtifact],
    exception_registry: LLMExceptionRegistry,
    sink: ModelSink,
    destination: ModelDestination,
    build_payload: Callable[[Mapping[str, bytes]], Payload],
    dispatch: Callable[[Payload], Result],
    record_audit: AuditRecorder,
) -> Result:
    """Authorize inputs before payload construction and dispatch to any model sink."""

    artifact_count = len(artifact_ids)
    exception_count = 0

    def deny(code: str) -> None:
        _deny(
            code,
            sink=sink,
            destination=destination,
            artifact_count=artifact_count,
            exception_count=exception_count,
            record_audit=record_audit,
        )

    if (
        not artifact_ids
        or len(artifact_ids) > MAX_MODEL_ARTIFACTS
        or len(artifact_ids) != len(set(artifact_ids))
        or len(artifacts) > MAX_MODEL_ARTIFACTS
    ):
        deny("untracked-input")

    states: dict[str, str] = {}

    def visit(artifact_id: str, *, selected: bool) -> None:
        artifact = artifacts.get(artifact_id)
        if artifact is None:
            deny("invalid-lineage")
        if artifact.policy in {
            LLMPolicy.SANITIZED_HEADER_EXCEPTION,
            LLMPolicy.EXAMPLE_HPO_EXCEPTION,
        } and not selected:
            deny("exception-not-inherited")
        state = states.get(artifact_id)
        if state == "visiting":
            deny("invalid-lineage")
        if state == "visited":
            return
        states[artifact_id] = "visiting"
        if artifact.policy is None:
            deny("unclassified-artifact")
        if artifact.policy is LLMPolicy.NO_LLM:
            deny("no-llm-artifact")
        if artifact.policy in {
            LLMPolicy.SANITIZED_HEADER_EXCEPTION,
            LLMPolicy.EXAMPLE_HPO_EXCEPTION,
        }:
            states[artifact_id] = "visited"
            return
        for parent_id in artifact.parent_artifact_ids:
            visit(parent_id, selected=False)
        states[artifact_id] = "visited"

    for artifact_id in artifact_ids:
        visit(artifact_id, selected=True)

    exceptions = {item.artifact_id: item for item in exception_registry.exceptions}
    content_by_id: dict[str, bytes] = {}
    for artifact_id in artifact_ids:
        artifact = artifacts[artifact_id]
        exception_kind: ExceptionKind | None = None
        if artifact.policy is LLMPolicy.SANITIZED_HEADER_EXCEPTION:
            exception_kind = ExceptionKind.SANITIZED_VCF_HEADER
        elif artifact.policy is LLMPolicy.EXAMPLE_HPO_EXCEPTION:
            exception_kind = ExceptionKind.NON_SENSITIVE_EXAMPLE_HPO
        if exception_kind is not None:
            exception_count += 1
            exception = exceptions.get(artifact.artifact_id)
            if (
                exception is None
                or exception.path != artifact.path
                or exception.sha256 != artifact.sha256
                or exception.kind is not exception_kind
            ):
                deny("exception-binding-invalid")
        try:
            content = read_bounded_bytes(
                artifact.path,
                maximum_bytes=MAX_MODEL_ARTIFACT_BYTES,
            )
        except BoundedJSONError:
            deny("artifact-read-invalid")
        if hashlib.sha256(content).hexdigest() != artifact.sha256:
            deny("artifact-identity-mismatch")
        try:
            if exception_kind is ExceptionKind.SANITIZED_VCF_HEADER:
                _validate_sanitized_header(content)
            elif exception_kind is ExceptionKind.NON_SENSITIVE_EXAMPLE_HPO:
                _validate_example_hpo(content)
        except LLMGuardError:
            deny("exception-content-invalid")
        content_by_id[artifact_id] = content

    _record_audit(
        record_audit,
        GuardAuditEvent(
            allowed=True,
            destination=destination,
            sink=sink,
            artifact_count=artifact_count,
            exception_count=exception_count,
            reason="authorized",
        ),
    )
    return dispatch(build_payload(content_by_id))
