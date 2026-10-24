"""Producer-neutral metadata and file-authentication kernels for bounded rescue.

Current delta, source-ledger replay and persisted merge are in wgs_read_rescue.
This module neither executes a producer nor stores a whole-callset document.
"""

from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass
from enum import StrEnum
import hashlib
import json
import os
from pathlib import Path, PurePosixPath
import re
import stat
from typing import Literal

from genosapiens.allele_identity import (
    AlleleKind,
    MitochondrialAllele,
    MosaicAllele,
    ReferenceRole,
    RepeatAllele,
    SequenceAllele,
    StructuralAllele,
)
from genosapiens.case_manifest import CaseFileClass
from genosapiens.evidence_model import (
    AlleleEvidence,
    EvidenceModelError,
    _Controlled,
)
from genosapiens.evidence_serialization import _plain
from genosapiens.governance import (
    ArtifactClassification,
    ArtifactRecord,
    GovernanceError,
    validate_relative_classification,
)
from genosapiens.reference_concordance import (
    RescueReferenceAuthorization,
)
from genosapiens.variant_inventory import WgsClass


ADDED_READ_INVENTORY_SCHEMA = "genosapiens.track1-added-read-inventory"
READ_IDENTITY_SCHEMA = "genosapiens.track1-read-vcf-identity"
PRODUCER_STACK_SCHEMA = "genosapiens.track1-evidence-producer-stack"
SELECTED_PLAN_SCHEMA = "genosapiens.track1-evidence-plan"

SHA256 = re.compile(r"^[0-9a-f]{64}$")
IDENTIFIER = re.compile(r"^[A-Za-z0-9][A-Za-z0-9._:-]*$")


class AuthorizedEvidenceError(RuntimeError):
    """Static failure at the controlled authorized-evidence boundary."""


def _fail(message: str = "invalid authorized evidence structure") -> None:
    raise AuthorizedEvidenceError(message)


def _sha(value: object) -> None:
    if type(value) is not str or SHA256.fullmatch(value) is None:
        _fail()


def _identifier(value: object) -> None:
    if type(value) is not str or IDENTIFIER.fullmatch(value) is None:
        _fail()


def _canonical_bytes(value: object) -> bytes:
    return (
        json.dumps(
            _plain(value),
            sort_keys=True,
            separators=(",", ":"),
            ensure_ascii=True,
            allow_nan=False,
        )
        + "\n"
    ).encode("ascii")


def _content_sha256(value: object) -> str:
    return hashlib.sha256(_canonical_bytes(value)).hexdigest()


class _BatchRecord(_Controlled):
    __slots__ = ()

    def __post_init__(self) -> None:
        try:
            _Controlled.__post_init__(self)
        except EvidenceModelError:
            _fail()


@dataclass(frozen=True, slots=True, repr=False)
class AddedReadFile(_BatchRecord):
    lane_id: str
    read_group_id: str
    file_class: CaseFileClass
    artifact: ArtifactRecord

    def __post_init__(self) -> None:
        _BatchRecord.__post_init__(self)
        _identifier(self.lane_id)
        _identifier(self.read_group_id)
        if (
            self.file_class not in {CaseFileClass.FASTQ_R1, CaseFileClass.FASTQ_R2}
            or self.artifact.classification is not ArtifactClassification.CONTROLLED
            or self.artifact.status != "present"
            or self.artifact.deleted_on is not None
            or self.artifact.retention != "delete"
            or self.artifact.kind != self.file_class.value
        ):
            _fail()


@dataclass(frozen=True, slots=True, repr=False)
class AddedReadInventory(_BatchRecord):
    parent_early_snapshot_sha256: str
    case_id: str
    sample_id: str
    reads: tuple[AddedReadFile, ...]
    schema: Literal["genosapiens.track1-added-read-inventory"] = (
        ADDED_READ_INVENTORY_SCHEMA
    )

    def __post_init__(self) -> None:
        _BatchRecord.__post_init__(self)
        _sha(self.parent_early_snapshot_sha256)
        _identifier(self.case_id)
        _identifier(self.sample_id)
        keys = tuple(
            (item.lane_id, item.read_group_id, item.file_class.value)
            for item in self.reads
        )
        pairs: dict[tuple[str, str], set[CaseFileClass]] = {}
        for item in self.reads:
            pairs.setdefault((item.lane_id, item.read_group_id), set()).add(
                item.file_class
            )
        if (
            self.schema != ADDED_READ_INVENTORY_SCHEMA
            or not self.reads
            or keys != tuple(sorted(keys))
            or len(keys) != len(set(keys))
            or len({item.artifact.artifact_id for item in self.reads})
            != len(self.reads)
            or len({item.artifact.relative_path for item in self.reads})
            != len(self.reads)
            or any(
                classes != {CaseFileClass.FASTQ_R1, CaseFileClass.FASTQ_R2}
                for classes in pairs.values()
            )
        ):
            _fail()

    @property
    def sha256(self) -> str:
        return _content_sha256(self)


class ProducerComponentRole(StrEnum):
    EXECUTABLE = "executable"
    MODEL = "model"
    CATALOG = "catalog"
    PARAMETERS = "parameters"
    IMAGE = "image"


@dataclass(frozen=True, slots=True, repr=False)
class ProducerComponent(_BatchRecord):
    role: ProducerComponentRole
    component_id: str
    sha256: str

    def __post_init__(self) -> None:
        _BatchRecord.__post_init__(self)
        _identifier(self.component_id)
        _sha(self.sha256)


@dataclass(frozen=True, slots=True, repr=False)
class ProducerStack(_BatchRecord):
    producer_id: str
    release: str
    components: tuple[ProducerComponent, ...]
    schema: Literal["genosapiens.track1-evidence-producer-stack"] = (
        PRODUCER_STACK_SCHEMA
    )

    def __post_init__(self) -> None:
        _BatchRecord.__post_init__(self)
        _identifier(self.producer_id)
        _identifier(self.release)
        keys = tuple((item.role.value, item.component_id) for item in self.components)
        if (
            self.schema != PRODUCER_STACK_SCHEMA
            or keys != tuple(sorted(keys))
            or len(keys) != len(set(keys))
            or len(self.components) != len(ProducerComponentRole)
            or {item.role for item in self.components} != set(ProducerComponentRole)
        ):
            _fail()

    @property
    def sha256(self) -> str:
        return _content_sha256(self)


@dataclass(frozen=True, slots=True, repr=False)
class PlanStage(_BatchRecord):
    stage_id: str
    variant_classes: tuple[WgsClass, ...]
    mandatory: bool

    def __post_init__(self) -> None:
        _BatchRecord.__post_init__(self)
        _identifier(self.stage_id)
        if (
            not self.variant_classes
            or self.variant_classes
            != tuple(sorted(set(self.variant_classes), key=lambda item: item.value))
        ):
            _fail()


@dataclass(frozen=True, slots=True, repr=False)
class SelectedPlan(_BatchRecord):
    plan_id: str
    stages: tuple[PlanStage, ...]
    schema: Literal["genosapiens.track1-evidence-plan"] = SELECTED_PLAN_SCHEMA

    def __post_init__(self) -> None:
        _BatchRecord.__post_init__(self)
        _identifier(self.plan_id)
        identifiers = tuple(item.stage_id for item in self.stages)
        if (
            self.schema != SELECTED_PLAN_SCHEMA
            or not self.stages
            or identifiers != tuple(sorted(identifiers))
            or len(identifiers) != len(set(identifiers))
        ):
            _fail()

    @property
    def sha256(self) -> str:
        return _content_sha256(self)


class StageStatus(StrEnum):
    SUCCEEDED = "succeeded"
    FAILED = "failed"
    INCONCLUSIVE = "inconclusive"
    NOT_RUN = "not-run"


@dataclass(frozen=True, slots=True, repr=False)
class StageState(_BatchRecord):
    stage_id: str
    status: StageStatus
    result_sha256: str | None
    reason_code: str | None

    def __post_init__(self) -> None:
        _BatchRecord.__post_init__(self)
        _identifier(self.stage_id)
        if self.result_sha256 is not None:
            _sha(self.result_sha256)
        if self.reason_code is not None:
            _identifier(self.reason_code)
        if (
            self.status is StageStatus.SUCCEEDED
            and (self.result_sha256 is None or self.reason_code is not None)
            or self.status is not StageStatus.SUCCEEDED
            and (self.result_sha256 is not None or self.reason_code is None)
        ):
            _fail()


class AssessmentOutcome(StrEnum):
    VALIDATED = "validated"
    REJECTED = "rejected"
    INCONCLUSIVE = "inconclusive"


@dataclass(frozen=True, slots=True, repr=False)
class ReadAssessment(_BatchRecord):
    assessment_id: str
    variant_class: WgsClass
    source_input_sha256: str
    sample_name: str
    allele_id: str
    native_identity_id: str
    source_record_id: str
    source_allele_id: str
    reasons: tuple[str, ...]
    scope: str
    outcome: AssessmentOutcome
    measurements: Mapping[str, object]
    measurement_artifact_sha256: str

    def __post_init__(self) -> None:
        _BatchRecord.__post_init__(self)
        for value in (
            self.assessment_id,
            self.sample_name,
            self.allele_id,
            self.native_identity_id,
            self.source_record_id,
            self.source_allele_id,
            self.scope,
        ):
            _identifier(value)
        _sha(self.source_input_sha256)
        _sha(self.measurement_artifact_sha256)
        for value in self.reasons:
            _identifier(value)
        if (
            not self.reasons
            or len(self.reasons) != len(set(self.reasons))
            or not self.measurements
            or any(IDENTIFIER.fullmatch(key) is None for key in self.measurements)
        ):
            _fail()


@dataclass(frozen=True, slots=True, repr=False)
class BatchAllele(_BatchRecord):
    variant_class: WgsClass
    evidence: AlleleEvidence

    def __post_init__(self) -> None:
        _BatchRecord.__post_init__(self)
        native = next(
            item.identity
            for item in self.evidence.identities
            if item.identity_id == self.evidence.native_identity_id
        )
        compatible = {
            AlleleKind.SEQUENCE: {WgsClass.SMALL_VARIANT, WgsClass.REGULATORY},
            AlleleKind.STRUCTURAL: {WgsClass.STRUCTURAL, WgsClass.COPY_NUMBER},
            AlleleKind.MITOCHONDRIAL: {WgsClass.MITOCHONDRIAL},
            AlleleKind.MOSAIC: {WgsClass.MOSAIC},
            AlleleKind.REPEAT: {WgsClass.REPEAT},
        }
        if self.variant_class not in compatible[native.kind]:
            _fail()


class ClassAssessmentStatus(StrEnum):
    ASSESSED = "assessed"
    ABSENT = "absent-in-assessed-scope"
    NOT_ASSESSED = "not-assessed"
    FAILED = "failed"
    INCONCLUSIVE = "inconclusive"


@dataclass(frozen=True, slots=True, repr=False)
class ClassAssessmentState(_BatchRecord):
    variant_class: WgsClass
    stage_ids: tuple[str, ...]
    status: ClassAssessmentStatus
    scope: str
    reason_code: str | None

    def __post_init__(self) -> None:
        _BatchRecord.__post_init__(self)
        for value in self.stage_ids:
            _identifier(value)
        _identifier(self.scope)
        if self.reason_code is not None:
            _identifier(self.reason_code)
        if (
            self.stage_ids != tuple(sorted(set(self.stage_ids)))
            or self.status
            in {ClassAssessmentStatus.ASSESSED, ClassAssessmentStatus.ABSENT}
            and self.reason_code is not None
            or self.status
            not in {ClassAssessmentStatus.ASSESSED, ClassAssessmentStatus.ABSENT}
            and self.reason_code is None
        ):
            _fail()


@dataclass(frozen=True, slots=True, repr=False)
class ReferenceAuthorizationBinding(_BatchRecord):
    source_input_sha256: str
    authorization: RescueReferenceAuthorization

    def __post_init__(self) -> None:
        _BatchRecord.__post_init__(self)
        authorization = self.authorization
        _sha(self.source_input_sha256)
        _identifier(authorization.reference_identity_id)
        for value in (
            authorization.reference_manifest_sha256,
            authorization.input_sha256,
            authorization.validated_input_sha256,
            authorization.concordance_report_sha256,
        ):
            _sha(value)
        if (
            authorization.authorized is not True
            or authorization.total_records <= 0
            or authorization.matching_records != authorization.total_records
            or authorization.input_sha256 != self.source_input_sha256
        ):
            _fail()


class EvidenceUse(StrEnum):
    VALIDATION_ONLY = "public-synthetic-validation-only"
    CONTROLLED_PRIMARY = "controlled-primary"


@dataclass(frozen=True, slots=True, repr=False)
class EvidenceUseBinding(_BatchRecord):
    intended_use: EvidenceUse
    candidate_stack_sha256: str
    release_freeze_sha256: str | None

    def __post_init__(self) -> None:
        _BatchRecord.__post_init__(self)
        _sha(self.candidate_stack_sha256)
        if self.release_freeze_sha256 is not None:
            _sha(self.release_freeze_sha256)
        if (
            self.intended_use is EvidenceUse.VALIDATION_ONLY
            and self.release_freeze_sha256 is not None
            or self.intended_use is EvidenceUse.CONTROLLED_PRIMARY
            and self.release_freeze_sha256 is None
        ):
            _fail()


class IdentityStatus(StrEnum):
    PASSED = "passed"
    CONFLICTING = "conflicting"
    INCONCLUSIVE = "inconclusive"


@dataclass(frozen=True, slots=True, repr=False)
class ReadVcfIdentityReport(_BatchRecord):
    parent_early_snapshot_sha256: str
    case_id: str
    sample_id: str
    source_input_sha256s: tuple[str, ...]
    validated_input_sha256s: tuple[str, ...]
    added_read_inventory_sha256: str
    identity_policy_sha256: str
    supporting_artifacts: tuple[ArtifactRecord, ...]
    informative_markers: int
    concordant_markers: int
    discordant_markers: int
    status: IdentityStatus
    schema: Literal["genosapiens.track1-read-vcf-identity"] = READ_IDENTITY_SCHEMA

    def __post_init__(self) -> None:
        _BatchRecord.__post_init__(self)
        _sha(self.parent_early_snapshot_sha256)
        _identifier(self.case_id)
        _identifier(self.sample_id)
        _sha(self.added_read_inventory_sha256)
        _sha(self.identity_policy_sha256)
        for value in (*self.source_input_sha256s, *self.validated_input_sha256s):
            _sha(value)
        expected_status = (
            IdentityStatus.INCONCLUSIVE
            if self.informative_markers == 0
            else IdentityStatus.PASSED
            if self.discordant_markers == 0
            else IdentityStatus.CONFLICTING
        )
        artifact_ids = tuple(item.artifact_id for item in self.supporting_artifacts)
        if (
            self.schema != READ_IDENTITY_SCHEMA
            or not self.source_input_sha256s
            or self.source_input_sha256s
            != tuple(sorted(set(self.source_input_sha256s)))
            or self.validated_input_sha256s
            != tuple(sorted(set(self.validated_input_sha256s)))
            or len(self.source_input_sha256s) != len(self.validated_input_sha256s)
            or any(
                type(value) is not int or isinstance(value, bool) or value < 0
                for value in (
                    self.informative_markers,
                    self.concordant_markers,
                    self.discordant_markers,
                )
            )
            or self.informative_markers
            != self.concordant_markers + self.discordant_markers
            or self.status is not expected_status
            or not self.supporting_artifacts
            or artifact_ids != tuple(sorted(artifact_ids))
            or len(artifact_ids) != len(set(artifact_ids))
            or any(
                item.classification is not ArtifactClassification.CONTROLLED
                or item.status != "present"
                or item.deleted_on is not None
                or item.retention != "delete"
                for item in self.supporting_artifacts
            )
            or not {"alignment-identity", "read-qc"}
            <= {item.kind for item in self.supporting_artifacts}
        ):
            _fail()

    @property
    def sha256(self) -> str:
        return _content_sha256(self)


@dataclass(frozen=True, slots=True, repr=False)
class ArtifactReplay:
    artifact_id: str
    path: Path

    def __post_init__(self) -> None:
        _identifier(self.artifact_id)
        if not isinstance(self.path, Path):
            _fail("invalid authorized evidence artifact replay")

    def __repr__(self) -> str:
        return "<controlled artifact replay; no-llm>"


@dataclass(frozen=True, slots=True, repr=False)
class ReferenceAuthorizationReplay:
    source_input_sha256: str
    identity_path: Path
    report_path: Path
    input_path: Path
    reference_manifest_path: Path
    validated_input_path: Path | None = None
    order_adjudication_path: Path | None = None
    source_report_path: Path | None = None
    source_reference_manifest_path: Path | None = None
    source_reference_root: Path | None = None

    def __post_init__(self) -> None:
        _sha(self.source_input_sha256)
        path_values = (
            self.identity_path,
            self.report_path,
            self.input_path,
            self.reference_manifest_path,
        )
        optional_paths = (
            self.validated_input_path,
            self.order_adjudication_path,
            self.source_report_path,
            self.source_reference_manifest_path,
            self.source_reference_root,
        )
        if any(not isinstance(value, Path) for value in path_values) or any(
            value is not None and not isinstance(value, Path) for value in optional_paths
        ):
            _fail("invalid reference authorization replay")

    def __repr__(self) -> str:
        return "<controlled reference authorization replay; no-llm>"


@dataclass(frozen=True, slots=True, repr=False)
class AuthorizedEvidenceExpectations:
    producer_stack: ProducerStack
    selected_plan: SelectedPlan
    use_binding: EvidenceUseBinding
    identity_policy_sha256: str

    def __post_init__(self) -> None:
        if (
            type(self.producer_stack) is not ProducerStack
            or type(self.selected_plan) is not SelectedPlan
            or type(self.use_binding) is not EvidenceUseBinding
        ):
            _fail("invalid authorized evidence expectations")
        _sha(self.identity_policy_sha256)

    def __repr__(self) -> str:
        return "<authorized evidence pins>"


def _assessment_target(item: ReadAssessment) -> tuple[object, ...]:
    return (
        item.source_input_sha256,
        item.sample_name,
        item.allele_id,
        item.native_identity_id,
        item.source_record_id,
        item.source_allele_id,
        item.reasons,
        item.scope,
    )


def _artifact_identity(path: Path, root: Path) -> tuple[str, int]:
    relative = path.relative_to(root)
    directory_flags = (
        os.O_RDONLY
        | getattr(os, "O_CLOEXEC", 0)
        | getattr(os, "O_DIRECTORY", 0)
        | getattr(os, "O_NOFOLLOW", 0)
    )
    file_flags = (
        os.O_RDONLY
        | getattr(os, "O_CLOEXEC", 0)
        | getattr(os, "O_NONBLOCK", 0)
        | getattr(os, "O_NOFOLLOW", 0)
    )
    directory_fd = os.open(root, directory_flags)
    try:
        for part in relative.parts[:-1]:
            child_fd = os.open(part, directory_flags, dir_fd=directory_fd)
            os.close(directory_fd)
            directory_fd = child_fd
        file_fd = os.open(relative.name, file_flags, dir_fd=directory_fd)
        try:
            before = os.fstat(file_fd)
            if not stat.S_ISREG(before.st_mode):
                raise ValueError
            digest = hashlib.sha256()
            size = 0
            while chunk := os.read(file_fd, 1024 * 1024):
                digest.update(chunk)
                size += len(chunk)
            after = os.fstat(file_fd)
            stable = lambda value: (
                value.st_dev,
                value.st_ino,
                value.st_mode,
                value.st_size,
                value.st_mtime_ns,
                value.st_ctime_ns,
            )
            if stable(before) != stable(after) or size != before.st_size:
                raise ValueError
            return digest.hexdigest(), size
        finally:
            os.close(file_fd)
    finally:
        os.close(directory_fd)


def _verify_artifacts(
    records: tuple[ArtifactRecord, ...],
    replays: tuple[ArtifactReplay, ...],
    *,
    reference_root: str | Path,
) -> None:
    try:
        root = Path(reference_root).resolve(strict=True)
        expected = {item.artifact_id: item for item in records}
        replayed = {item.artifact_id: item.path for item in replays}
        if (
            len(expected) != len(records)
            or len({item.relative_path for item in records}) != len(records)
            or set(expected) != set(replayed)
        ):
            raise ValueError
        for artifact_id, record in expected.items():
            relative = PurePosixPath(record.relative_path)
            validate_relative_classification(
                relative.as_posix(), ArtifactClassification.CONTROLLED
            )
            path = replayed[artifact_id]
            if not path.is_absolute() or path != root.joinpath(*relative.parts):
                raise ValueError
            sha256, size = _artifact_identity(path, root)
            if (sha256, size) != (record.sha256, record.size_bytes):
                raise ValueError
    except (GovernanceError, OSError, ValueError):
        raise AuthorizedEvidenceError(
            "authorized evidence artifact replay failed"
        ) from None


def _identity_reference(identity: object) -> tuple[ReferenceRole, str]:
    if type(identity) in {SequenceAllele, StructuralAllele, RepeatAllele}:
        return identity.reference_role, identity.reference_identity_id
    if type(identity) in {MitochondrialAllele, MosaicAllele}:
        return identity.sequence.reference_role, identity.sequence.reference_identity_id
    _fail()
