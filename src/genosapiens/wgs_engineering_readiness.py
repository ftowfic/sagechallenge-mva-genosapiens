"""Freeze and independently read back bounded-WGS engineering readiness."""

from __future__ import annotations

from dataclasses import asdict, dataclass
import hashlib
import json
import os
from pathlib import Path
import re
import stat
from typing import Literal

from genosapiens.track1_public_validation_coordinator import executable_stack_bytes
from genosapiens.track1_release import EVIDENCE_INTERFACE_FILES
from genosapiens.wgs_ledger import canonical


READINESS_SCHEMA = "genosapiens.track1-bounded-wgs-engineering-readiness"
CANDIDATE_SCHEMA = "genosapiens.track1-bounded-wgs-engineering-candidate"
_SHA = re.compile(r"^[0-9a-f]{64}$")
_REVISION = re.compile(r"^[0-9a-f]{40}$")
_MAX_RECORD_BYTES = 1024 * 1024
RESOURCE_PATHS = (
    "config/public-resources/gnomad-joint-faf95-precision-v1.json",
    "config/public-resources/track1-common-allele-exception-curation-v1.json",
    "config/public-resources/track1-common-allele-exceptions-v1.json",
    "config/track1/bounded-wgs-execution-v5.json",
    "config/track1/common-allele-prefilter-v1.json",
    "config/track1/analysis-eligibility-v1.json",
    "config/track1/context-writer-runtime-v1.json",
    "config/track1/public-ranking-development-policy-family.json",
    "config/track1/submission-acceptance-profile.json",
)
REQUIRED_COMPONENT_PATHS = (
    *RESOURCE_PATHS,
    "config/track1/handoff-contract.json",
    "src/genosapiens/wgs_analysis_policy.py",
    "src/genosapiens/wgs_analysis_inputs.py",
    "src/genosapiens/wgs_admission_inputs.py",
    "src/genosapiens/wgs_association_codec.py",
    "src/genosapiens/wgs_admission_graph.py",
    "src/genosapiens/wgs_admission_screening.py",
    "src/genosapiens/wgs_allele_query_cache.py",
    "src/genosapiens/candidate_admission.py",
    "src/genosapiens/annotation_admission.py",
    "src/genosapiens/track1_methods.py",
    "src/genosapiens/wgs_context_admission.py",
    "src/genosapiens/wgs_context_payloads.py",
    "src/genosapiens/wgs_context_pipeline.py",
    "src/genosapiens/wgs_annotation_prefix.py",
    "config/track1/handoff-pin.json",
    "schemas/track1/mechanism-ready-handoff.schema.json",
    "scripts/build_common_allele_public_indexes.py",
    "scripts/build_public_wgs_unprefixed_derivative.py",
    "scripts/build_track1_common_allele_exception_bundle.py",
    "scripts/build_track1_common_allele_exception_sources.py",
    "scripts/operator_compatibility_repair/analyze_bounded_case.py",
    "scripts/operator_compatibility_repair/archive_bounded_resume.py",
    "scripts/operator_compatibility_repair/archive_failed_preparation.py",
    "scripts/operator_compatibility_repair/check_bounded_host_capacity.py",
    "scripts/operator_compatibility_repair/check_stage_prerequisite.py",
    "scripts/operator_compatibility_repair/finalize_bounded_runtime.py",
    "scripts/operator_compatibility_repair/prepare_bounded_case.py",
    "scripts/operator_compatibility_repair/run_bounded_stage.py",
    "scripts/operator_compatibility_repair/seal_bounded_stage_evidence.py",
    "scripts/prepare_track1_operator_bounded_wgs.py",
    "scripts/qualify_gnomad_faf95_precision.py",
    "src/genosapiens/bounded_wgs_config.py",
    "src/genosapiens/common_allele_exception_bundle.py",
    "src/genosapiens/common_allele_filter.py",
    "src/genosapiens/common_allele_public_indexes.py",
    "src/genosapiens/common_allele_resources.py",
    "src/genosapiens/gnomad_faf95_adapter.py",
    "src/genosapiens/population_clinical_session.py",
    "src/genosapiens/population_annotation.py",
    "src/genosapiens/evidence_serialization.py",
    "src/genosapiens/wgs_admission.py",
    "src/genosapiens/wgs_analysis.py",
    "src/genosapiens/wgs_annotation.py",
    "src/genosapiens/wgs_contexts.py",
    "src/genosapiens/wgs_clinical_session.py",
    "src/genosapiens/wgs_engineering_readiness.py",
    "src/genosapiens/wgs_failures.py",
    "src/genosapiens/wgs_filter.py",
    "src/genosapiens/wgs_handoff.py",
    "src/genosapiens/wgs_ledger.py",
    "src/genosapiens/wgs_longtr.py",
    "src/genosapiens/wgs_matching_validation.py",
    "src/genosapiens/wgs_operator_package.py",
    "src/genosapiens/wgs_preparation.py",
    "src/genosapiens/wgs_prepared_case.py",
    "src/genosapiens/wgs_provenance.py",
    "src/genosapiens/wgs_public_analysis.py",
    "src/genosapiens/wgs_public_case.py",
    "src/genosapiens/wgs_public_execution.py",
    "src/genosapiens/wgs_public_filter.py",
    "src/genosapiens/wgs_public_qualification.py",
    "src/genosapiens/wgs_public_replay.py",
    "src/genosapiens/wgs_public_runner.py",
    "src/genosapiens/wgs_public_store.py",
    "src/genosapiens/wgs_query_limits.py",
    "src/genosapiens/wgs_read_rescue.py",
    "src/genosapiens/wgs_reporting.py",
    "src/genosapiens/wgs_resume.py",
    "src/genosapiens/wgs_runtime.py",
    "src/genosapiens/wgs_scientific.py",
    "src/genosapiens/wgs_sensitivity.py",
    "src/genosapiens/wgs_shards.py",
    "src/genosapiens/wgs_source_routing.py",
    "src/genosapiens/wgs_validation_audit.py",
    "src/genosapiens/wgs_validation_authority.py",
    "src/genosapiens/wgs_validation_metrics.py",
)
# Current evidence components are explicit requirements in addition to the
# executable stack's complete Python/config inventory.
REQUIRED_COMPONENT_PATHS = tuple(dict.fromkeys((*REQUIRED_COMPONENT_PATHS,
    *("src/genosapiens/" + name for name in EVIDENCE_INTERFACE_FILES))))

BLOCKERS = (
    "affected-public-acceptance-not-qualified",
    "full-wgs-capacity-not-qualified",
    "operator-package-not-qualified",
    "public-pilot-not-qualified",
    "superseding-release-not-frozen",
)
GATE_NAMES = (
    "public_pilot",
    "full_wgs_capacity",
    "affected_public_acceptance",
    "superseding_release",
    "operator_package",
)


class WgsEngineeringReadinessError(RuntimeError):
    """Public static readiness evidence is incomplete or changed."""


def _digest(value: object) -> bool:
    return type(value) is str and _SHA.fullmatch(value) is not None


def _decode(content: bytes) -> dict:
    try:
        value = json.loads(
            content,
            object_pairs_hook=lambda pairs: (
                dict(pairs)
                if len(pairs) == len({key for key, _item in pairs})
                else (_ for _ in ()).throw(ValueError())
            ),
        )
        if type(value) is not dict:
            raise ValueError
        return value
    except (UnicodeError, ValueError, TypeError, RecursionError):
        raise WgsEngineeringReadinessError(
            "engineering readiness record is invalid"
        ) from None


def _read_bounded(path: Path) -> bytes:
    descriptor = -1
    try:
        metadata = path.lstat()
        if (
            path.is_symlink()
            or not stat.S_ISREG(metadata.st_mode)
            or metadata.st_size <= 0
            or metadata.st_size > _MAX_RECORD_BYTES
        ):
            raise OSError
        descriptor = os.open(path, os.O_RDONLY | os.O_NOFOLLOW)
        with os.fdopen(descriptor, "rb", closefd=True) as stream:
            descriptor = -1
            content = stream.read(_MAX_RECORD_BYTES + 1)
        if len(content) != metadata.st_size:
            raise OSError
        return content
    except OSError:
        raise WgsEngineeringReadinessError(
            "engineering readiness source is unavailable"
        ) from None
    finally:
        if descriptor >= 0:
            os.close(descriptor)


@dataclass(frozen=True, slots=True)
class EngineeringFileBinding:
    path: str
    sha256: str
    size_bytes: int

    def __post_init__(self) -> None:
        path = Path(self.path)
        if (
            not self.path
            or path.is_absolute()
            or path.as_posix() != self.path
            or any(part in {"", ".", ".."} for part in path.parts)
            or not _digest(self.sha256)
            or type(self.size_bytes) is not int
            or self.size_bytes <= 0
        ):
            raise WgsEngineeringReadinessError(
                "engineering file binding is invalid"
            )


@dataclass(frozen=True, slots=True)
class EngineeringGate:
    name: str
    passed: Literal[False] = False

    def __post_init__(self) -> None:
        if self.name not in GATE_NAMES or self.passed is not False:
            raise WgsEngineeringReadinessError(
                "engineering qualification gate is invalid"
            )


@dataclass(frozen=True, slots=True)
class EngineeringCandidate:
    baseline_revision: str
    readiness_path: str
    readiness_sha256: str
    executable_stack_sha256: str
    executable_stack_file_count: int
    runtime_image_digest: str
    required_component_paths: tuple[str, ...]
    resources: tuple[EngineeringFileBinding, ...]
    gates: tuple[EngineeringGate, ...]
    blockers: tuple[str, ...]
    completed_source_input_sha256: str | None = None
    schema: Literal[
        "genosapiens.track1-bounded-wgs-engineering-candidate"
    ] = CANDIDATE_SCHEMA
    version: Literal[3] = 3
    optional_qualification_status: Literal["not_qualified"] = "not_qualified"
    classification: Literal["public"] = "public"
    status: Literal[
        "engineering-ready-awaiting-public-qualification"
    ] = "engineering-ready-awaiting-public-qualification"
    controlled_use_authorized: Literal[False] = False

    def __post_init__(self) -> None:
        if (
            type(self.version) is not int or self.version != 3
            or self.optional_qualification_status != "not_qualified"
            or self.completed_source_input_sha256 is not None and not _digest(self.completed_source_input_sha256)
            or _REVISION.fullmatch(self.baseline_revision) is None
            or not _digest(self.readiness_sha256)
            or not _digest(self.executable_stack_sha256)
            or type(self.executable_stack_file_count) is not int
            or self.executable_stack_file_count <= 0
            or not self.runtime_image_digest.startswith(
                "genosapiens/track1-pipeline@sha256:"
            )
            or not _digest(self.runtime_image_digest.rsplit(":", 1)[-1])
            or self.required_component_paths != REQUIRED_COMPONENT_PATHS
            or tuple(value.path for value in self.resources) != RESOURCE_PATHS
            or tuple(value.name for value in self.gates) != GATE_NAMES
            or any(value.passed is not False for value in self.gates)
            or self.blockers != BLOCKERS
            or Path(self.readiness_path).is_absolute()
            or self.controlled_use_authorized is not False
        ):
            raise WgsEngineeringReadinessError(
                "engineering candidate is invalid"
            )

    def as_bytes(self) -> bytes:
        return canonical(asdict(self))

    @property
    def sha256(self) -> str:
        # Bind readiness bytes, excluding only their identity-derived locator.
        identity = asdict(self)
        identity.pop("readiness_path")
        return hashlib.sha256(canonical(identity)).hexdigest()


def _validate_binding(root: Path, value: dict) -> EngineeringFileBinding:
    try:
        binding = EngineeringFileBinding(
            value["path"], value["sha256"], value["size_bytes"],
        )
        if set(value) != {"path", "sha256", "size_bytes"}:
            raise ValueError
        content = _read_bounded(root / binding.path)
        if (
            len(content) != binding.size_bytes
            or hashlib.sha256(content).hexdigest() != binding.sha256
        ):
            raise ValueError
        return binding
    except (KeyError, TypeError, ValueError):
        raise WgsEngineeringReadinessError(
            "engineering file binding changed"
        ) from None


def freeze_engineering_candidate(
    repository_root: str | Path,
    *,
    readiness_path: str,
    expected_readiness_sha256: str,
) -> EngineeringCandidate:
    """Authenticate current public readiness and freeze the blocked candidate."""

    root = Path(repository_root).resolve(strict=True)
    if Path(readiness_path).is_absolute() or not _digest(expected_readiness_sha256):
        raise WgsEngineeringReadinessError("engineering readiness request is invalid")
    readiness_content = _read_bounded(root / readiness_path)
    if hashlib.sha256(readiness_content).hexdigest() != expected_readiness_sha256:
        raise WgsEngineeringReadinessError("engineering readiness identity changed")
    readiness = _decode(readiness_content)
    try:
        invocations = readiness["qualification_invocations"]
        source_review = readiness["source_boundary_review"]
        if (
            readiness["schema"] != READINESS_SCHEMA
            or readiness["version"] != 1
            or readiness["classification"] != "public"
            or readiness["qualification_status"] != "engineering-ready"
            or readiness["controlled_use_authorized"] is not False
            or readiness["controlled_input_accessed"] is not False
            or readiness["scientific_invocations"] != 0
            or readiness["full_wgs_capacity_qualified"] is not False
            or readiness["clinical_sensitivity_qualified"] is not False
            or type(invocations) is not list
            or not invocations
            or any(value.get("status") != "passed" for value in invocations)
            or sum(value["tests_run"] for value in invocations)
            != readiness["test_count"]
            or source_review["eager_whole_file_or_fetchall_calls_in_runtime_core"] != 0
            or source_review["forbidden_driver_dependencies"] != 0
            or source_review["unbounded_child_output_paths"] != 0
            or source_review["shell_help_and_syntax"] != "passed"
            or source_review["source_credential_scan_status"] != "passed"
        ):
            raise ValueError
    except (KeyError, TypeError, ValueError):
        raise WgsEngineeringReadinessError(
            "engineering readiness gates are incomplete"
        ) from None

    stack_content = executable_stack_bytes(root)
    stack_sha256 = hashlib.sha256(stack_content).hexdigest()
    stack = _decode(stack_content)
    try:
        stack_paths = tuple(value["path"] for value in stack["files"])
        if (
            stack_sha256 != readiness["executable_stack"]["sha256"]
            or len(stack_paths) != readiness["executable_stack"]["file_count"]
            or any(path not in stack_paths for path in REQUIRED_COMPONENT_PATHS)
        ):
            raise ValueError
    except (KeyError, TypeError, ValueError):
        raise WgsEngineeringReadinessError(
            "engineering executable stack is incomplete"
        ) from None

    resources_by_path = {
        value.get("path"): value
        for value in readiness.get("resource_and_config_sources", ())
        if type(value) is dict
    }
    if set(resources_by_path) != set(RESOURCE_PATHS):
        raise WgsEngineeringReadinessError(
            "engineering resource closure is incomplete"
        )
    resources = tuple(
        _validate_binding(root, resources_by_path[path]) for path in RESOURCE_PATHS
    )
    for name in ("source_entrypoints", "test_sources"):
        values = readiness.get(name)
        if type(values) is not list or not values:
            raise WgsEngineeringReadinessError(
                "engineering source closure is incomplete"
            )
        paths = set()
        for value in values:
            try:
                path = value["path"]
                expected = value["sha256"]
                if set(value) != {"path", "sha256"} or path in paths:
                    raise ValueError
                content = _read_bounded(root / path)
                if not _digest(expected) or hashlib.sha256(content).hexdigest() != expected:
                    raise ValueError
                paths.add(path)
            except (KeyError, TypeError, ValueError):
                raise WgsEngineeringReadinessError(
                    "engineering source closure changed"
                ) from None
    image = readiness.get("production_runtime_image", {}).get("digest")
    return EngineeringCandidate(
        readiness["baseline_revision"],
        readiness_path,
        expected_readiness_sha256,
        stack_sha256,
        len(stack_paths),
        image,
        REQUIRED_COMPONENT_PATHS,
        resources,
        tuple(EngineeringGate(name) for name in GATE_NAMES),
        BLOCKERS,
        readiness.get("completed_source_input_sha256"),
    )


def load_engineering_candidate(
    content: bytes,
    *,
    expected_sha256: str,
) -> EngineeringCandidate:
    """Strictly decode canonical frozen candidate bytes."""

    if (
        type(content) is not bytes
        or not _digest(expected_sha256)
    ):
        raise WgsEngineeringReadinessError("engineering candidate identity changed")
    value = _decode(content)
    try:
        resources = tuple(
            EngineeringFileBinding(**item) for item in value.pop("resources")
        )
        gates = tuple(EngineeringGate(**item) for item in value.pop("gates"))
        value["required_component_paths"] = tuple(value["required_component_paths"])
        value["blockers"] = tuple(value["blockers"])
        candidate = EngineeringCandidate(**value, resources=resources, gates=gates)
        if candidate.as_bytes() != content or candidate.sha256 != expected_sha256:
            raise ValueError
        return candidate
    except (KeyError, TypeError, ValueError):
        raise WgsEngineeringReadinessError(
            "engineering candidate is invalid"
        ) from None


def read_back_engineering_candidate(
    content: bytes,
    *,
    expected_sha256: str,
    repository_root: str | Path,
) -> EngineeringCandidate:
    """Recompute the complete closure and compare it to frozen candidate bytes."""

    candidate = load_engineering_candidate(
        content, expected_sha256=expected_sha256,
    )
    current = freeze_engineering_candidate(
        repository_root,
        readiness_path=candidate.readiness_path,
        expected_readiness_sha256=candidate.readiness_sha256,
    )
    if current != candidate:
        raise WgsEngineeringReadinessError("engineering candidate readback changed")
    return candidate
