"""Public-only authority for bounded-WGS qualification invocations.

This authority is deliberately narrower than a Track 1 base release.  It can
authorize the named public pilot or full-WGS qualification, but it can never
authorize controlled input use, an accepted release, or operator packaging.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass
import hashlib
import json
import os
from pathlib import Path
import re
import stat
from typing import Literal

from genosapiens.wgs_preparation import CompletedSourceArtifact, COMPLETED_SOURCE_SCHEMA
from genosapiens.bounded_json import read_bounded_bytes
from genosapiens.phenotype_ingestion import MAX_HPO_LIST_BYTES
from genosapiens.bounded_wgs_config import (
    BoundedWgsExecutionConfig,
    bounded_wgs_execution_config_as_json,
)
from genosapiens.track1_release import BaseRelease


SCHEMA = "genosapiens.track1-bounded-wgs-public-qualification-authority"
VERSION = 3
_SHA = re.compile(r"^[0-9a-f]{64}$")
_IMAGE = re.compile(r"^sha256:[0-9a-f]{64}$")
_INVOCATIONS = frozenset({"bounded-public-pilot", "complete-hg002-deepvariant-wgs", "focused-fresh-subsample"})
_PUBLIC_INPUT_ROOTS = (
    Path("/mnt/reference/releases/deepvariant"),
    Path("/mnt/reference/scratch/public-build"),
)
_PUBLIC_OUTPUT_ROOT = Path(
    "/mnt/reference/scratch/public-build/track1-bounded-wgs-qualification-20260919"
)
_R11_PUBLIC_OUTPUT_ROOT = Path(
    "/mnt/reference/scratch/public-build/track1-bounded-wgs-qualification-20260922-r11"
)
_R12_PUBLIC_OUTPUT_ROOT = Path(
    "/mnt/reference/scratch/public-build/track1-bounded-wgs-qualification-20260922-r12"
)
_R12C_PUBLIC_OUTPUT_ROOT = Path(
    "/mnt/reference/scratch/public-build/track1-bounded-wgs-qualification-20260922-r12c"
)
_R12D_PUBLIC_OUTPUT_ROOT = Path(
    "/mnt/reference/scratch/public-build/track1-bounded-wgs-qualification-20260922-r12d"
)
_R12E_PUBLIC_OUTPUT_ROOT = Path(
    "/mnt/reference/scratch/public-build/track1-bounded-wgs-qualification-20260922-r12e"
)
_PUBLIC_HPO_ROOT = Path("/workspace/tests/fixtures/track1/roster")
_FORBIDDEN_PARTS = frozenset({"controlled-work", "controlled", "data"})


class WgsPublicQualificationError(RuntimeError):
    """Static refusal at a public-qualification boundary."""


def _digest(value: object) -> bool:
    return type(value) is str and _SHA.fullmatch(value) is not None


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


def _public_path(value: str, *, roots: tuple[Path, ...]) -> Path:
    try:
        path = Path(value)
        if (
            type(value) is not str
            or not path.is_absolute()
            or value != path.as_posix()
            or ".." in path.parts
            or any(part in _FORBIDDEN_PARTS for part in path.parts)
            or not any(path.is_relative_to(root) for root in roots)
        ):
            raise ValueError
        return path
    except (OSError, TypeError, ValueError):
        raise WgsPublicQualificationError(
            "public qualification path is outside its public boundary"
        ) from None


def _public_output_roots() -> tuple[Path, ...]:
    """Return the historical and exact replacement qualification roots."""

    return (
        _PUBLIC_OUTPUT_ROOT,
        _R11_PUBLIC_OUTPUT_ROOT,
        _R12_PUBLIC_OUTPUT_ROOT,
        _R12C_PUBLIC_OUTPUT_ROOT,
        _R12D_PUBLIC_OUTPUT_ROOT,
        _R12E_PUBLIC_OUTPUT_ROOT,
    )


@dataclass(frozen=True, slots=True)
class QualificationBinding:
    role: str
    sha256: str

    def __post_init__(self) -> None:
        if (
            type(self.role) is not str
            or re.fullmatch(r"[a-z][a-z0-9_]{1,95}", self.role) is None
            or not _digest(self.sha256)
        ):
            raise WgsPublicQualificationError(
                "public qualification binding is invalid"
            )


@dataclass(frozen=True, slots=True)
class PublicQualificationAuthority:
    invocation_id: Literal[
        "bounded-public-pilot", "complete-hg002-deepvariant-wgs", "focused-fresh-subsample"
    ]
    plan_sha256: str
    authorization_receipt_sha256: str
    engineering_candidate_sha256: str
    executable_stack_sha256: str
    runtime_image_digest: str
    execution_config_sha256: str
    input_path: str
    input_sha256: str
    input_size_bytes: int
    hpo_input_path: str
    hpo_input_sha256: str
    hpo_input_size_bytes: int
    output_root: str
    preparation_run_options_sha256: str
    policy_bindings: tuple[QualificationBinding, ...]
    resource_bindings: tuple[QualificationBinding, ...]
    completed_source_input: CompletedSourceArtifact | None = None
    preparation_attempts: Literal[0, 1] = 1
    analysis_attempts: Literal[1] = 1
    automatic_retries: Literal[0] = 0
    resume_attempts: Literal[0] = 0
    schema: Literal[
        "genosapiens.track1-bounded-wgs-public-qualification-authority"
    ] = SCHEMA
    version: Literal[3] = VERSION
    classification: Literal["public"] = "public"
    controlled_input_accessed: Literal[False] = False
    controlled_use_authorized: Literal[False] = False
    accepted_release_claimed: Literal[False] = False
    operator_package_authorized: Literal[False] = False

    def __post_init__(self) -> None:
        hashes = (
            self.plan_sha256,
            self.authorization_receipt_sha256,
            self.engineering_candidate_sha256,
            self.executable_stack_sha256,
            self.execution_config_sha256,
            self.input_sha256,
            self.hpo_input_sha256,
            self.preparation_run_options_sha256,
        )
        if (
            self.invocation_id not in _INVOCATIONS
            or any(not _digest(value) for value in hashes)
            or _IMAGE.fullmatch(self.runtime_image_digest) is None
            or type(self.input_size_bytes) is not int
            or self.input_size_bytes <= 0
            or type(self.hpo_input_size_bytes) is not int
            or self.hpo_input_size_bytes <= 0
            or type(self.preparation_attempts) is not int or self.preparation_attempts not in {0, 1}
            or self.analysis_attempts != 1
            or self.automatic_retries != 0
            or self.resume_attempts != 0
            or self.schema != SCHEMA
            or self.version != VERSION
            or self.classification != "public"
            or self.controlled_input_accessed is not False
            or self.controlled_use_authorized is not False
            or self.accepted_release_claimed is not False
            or self.operator_package_authorized is not False
        ):
            raise WgsPublicQualificationError(
                "public qualification authority is invalid"
            )
        if self.completed_source_input is not None:
            if (type(self.completed_source_input) is not CompletedSourceArtifact
                    or self.completed_source_input.schema != COMPLETED_SOURCE_SCHEMA):
                raise WgsPublicQualificationError("completed source authority is invalid")
            _public_path(self.completed_source_input.path, roots=_PUBLIC_INPUT_ROOTS)
        _public_path(self.input_path, roots=_PUBLIC_INPUT_ROOTS)
        _public_path(self.hpo_input_path, roots=(_PUBLIC_HPO_ROOT,))
        output_roots = (
            *_public_output_roots(),
            Path(
                "/mnt/reference/scratch/public-build/"
                f"track1-bounded-wgs-acceleration-{self.engineering_candidate_sha256}"
            ),
        )
        output_roots += (Path("/mnt/reference/scratch/public-build/"
            f"track1-bounded-wgs-prefilter-replacement-{self.engineering_candidate_sha256}"),)
        fresh_root = Path('/mnt/reference/scratch/public-build/'
            f'track1-focused-fresh-{self.engineering_candidate_sha256}')
        if self.invocation_id == 'focused-fresh-subsample':
            output_roots = (fresh_root,)
        output = _public_path(self.output_root, roots=output_roots)
        if (
            output.parent not in output_roots
            or output.name != "full-wgs"
            and re.fullmatch(r"pilot(?:-[0-9]{2})?", output.name) is None
            and not (self.invocation_id == 'focused-fresh-subsample'
                and re.fullmatch(r'case-[0-9]{2}', output.name))
        ):
            raise WgsPublicQualificationError(
                "public qualification output root is not frozen"
            )
        for values in (self.policy_bindings, self.resource_bindings):
            if (
                type(values) is not tuple
                or not values
                or any(type(value) is not QualificationBinding for value in values)
                or values != tuple(sorted(values, key=lambda value: value.role))
                or len({value.role for value in values}) != len(values)
            ):
                raise WgsPublicQualificationError(
                    "public qualification bindings are incomplete"
                )
        roles = [value.role for value in (*self.policy_bindings, *self.resource_bindings)]
        from genosapiens.wgs_work_contract import ROLES as work_roles, FRESH_INVOCATION
        work = {value.role for value in self.resource_bindings} & work_roles
        if (any(value.role in work_roles for value in self.policy_bindings)
                or work and (work != work_roles or self.invocation_id != FRESH_INVOCATION)):
            raise WgsPublicQualificationError('focused work bindings are incomplete')
        if self.invocation_id == FRESH_INVOCATION:
            required = {*work_roles, 'analysis_profile_policy', 'analysis_eligibility_policy',
                'context_writer_profile', 'fresh_input_selection', 'fresh_coverage_matrix', 'fresh_acceptance',
                'fresh_selected_sample'}
            resource_roles = {value.role for value in self.resource_bindings}
            if (not required <= resource_roles or self.completed_source_input is not None
                    or self.preparation_attempts != 1
                    or any(role.startswith(('completed_', 'public_profile_pilot_')) for role in roles)
                    or output.parent != fresh_root or re.fullmatch(r'case-[0-9]{2}', output.name) is None):
                raise WgsPublicQualificationError('fresh qualification requires original public inputs')
        completed = any(value.role == "completed_analysis_input" for value in self.resource_bindings)
        prefix = any(value.role == "completed_annotation_prefix" for value in self.resource_bindings)
        admission = any(value.role == "completed_admission_input" for value in self.resource_bindings)
        eligibility_roles={'analysis_eligibility_policy','context_writer_profile'}
        eligibility={value.role for value in self.resource_bindings if value.role in eligibility_roles}
        pilot_roles={'public_profile_pilot_acceptance','public_profile_pilot_cohort'}
        pilot={value.role for value in self.resource_bindings if value.role in pilot_roles}
        context_restart_role = 'public_profile_pilot_context_inputs'
        context_restart = any(value.role == context_restart_role for value in self.resource_bindings)
        assessment_restart_role = 'public_profile_pilot_assessment_inputs'
        assessment_restart = any(value.role == assessment_restart_role for value in self.resource_bindings)
        scientific_restart_role = 'public_profile_pilot_scientific_inputs'
        scientific_restart = any(value.role == scientific_restart_role for value in self.resource_bindings)
        primary_restart_role = 'public_profile_pilot_primary_inputs'
        primary_restart = any(value.role == primary_restart_role for value in self.resource_bindings)
        finite_pilot=(pilot==pilot_roles and completed and prefix and admission
            and self.invocation_id=='bounded-public-pilot'
            and re.fullmatch(r'pilot(?:-[0-9]{2})?',output.name) is not None)
        if (len(roles) != len(set(roles))
                or any(value.role in {"completed_analysis_input", "completed_annotation_prefix", "completed_admission_input"} for value in self.policy_bindings)
                or any(value.role in eligibility_roles for value in self.policy_bindings)
                or any(value.role in pilot_roles for value in self.policy_bindings)
                or any(value.role == context_restart_role for value in self.policy_bindings)
                or context_restart and not finite_pilot
                or any(value.role == assessment_restart_role for value in self.policy_bindings)
                or assessment_restart and not (finite_pilot and context_restart)
                or any(value.role == scientific_restart_role for value in self.policy_bindings)
                or scientific_restart and (not finite_pilot or context_restart or assessment_restart)
                or any(value.role == primary_restart_role for value in self.policy_bindings)
                or primary_restart and (not finite_pilot or context_restart or assessment_restart or scientific_restart)
                or pilot and not finite_pilot
                or any(value.role == 'analysis_profile_policy' for value in self.policy_bindings)
                or eligibility and eligibility!=eligibility_roles
                or prefix and not completed
                or admission and not (completed and prefix)
                or (self.preparation_attempts == 0) != completed
                or completed and not finite_pilot and (self.invocation_id != "complete-hg002-deepvariant-wgs"
                                  or output.name != "full-wgs")):
            raise WgsPublicQualificationError("analysis-only input authority is invalid")

    def as_bytes(self) -> bytes:
        return _canonical(asdict(self))

    @property
    def is_finite_profile_pilot(self) -> bool:
        """A declared finite consumer of completed public whole-genome inputs."""
        roles = {value.role for value in self.resource_bindings}
        return (self.invocation_id == "bounded-public-pilot"
                and self.preparation_attempts == 0
                and {"public_profile_pilot_acceptance", "public_profile_pilot_cohort"} <= roles)

    @property
    def sha256(self) -> str:
        return hashlib.sha256(self.as_bytes()).hexdigest()

    def binding(self, role: str) -> str:
        values = {
            value.role: value.sha256
            for value in (*self.policy_bindings, *self.resource_bindings)
        }
        try:
            return values[role]
        except KeyError:
            raise WgsPublicQualificationError(
                "public qualification binding is missing"
            ) from None


def load_public_qualification_authority(
    content: bytes, *, expected_sha256: str,
) -> PublicQualificationAuthority:
    """Load exact canonical authority bytes with duplicate-key rejection."""

    try:
        if (
            type(content) is not bytes
            or not _digest(expected_sha256)
            or hashlib.sha256(content).hexdigest() != expected_sha256
        ):
            raise ValueError

        def pairs(items):
            value = {}
            for key, child in items:
                if key in value:
                    raise ValueError
                value[key] = child
            return value

        data = json.loads(content, object_pairs_hook=pairs)
        expected = set(PublicQualificationAuthority.__dataclass_fields__)
        if type(data) is not dict or set(data) != expected:
            raise ValueError
        for name in ("policy_bindings", "resource_bindings"):
            data[name] = tuple(QualificationBinding(**value) for value in data[name])
        if data["completed_source_input"] is not None:
            data["completed_source_input"] = CompletedSourceArtifact(**data["completed_source_input"])
        authority = PublicQualificationAuthority(**data)
        if authority.as_bytes() != content:
            raise ValueError
        return authority
    except (KeyError, TypeError, ValueError, json.JSONDecodeError):
        raise WgsPublicQualificationError(
            "public qualification authority is invalid or mismatched"
        ) from None


def require_private_qualification_hpo_input(
    authority: PublicQualificationAuthority,
) -> None:
    """Authenticate the tiny private-mode HPO input before large resources."""

    try:
        path = _public_path(
            authority.hpo_input_path, roots=(_PUBLIC_HPO_ROOT,),
        )
        metadata = path.lstat()
        if (
            path.is_symlink()
            or not stat.S_ISREG(metadata.st_mode)
            or metadata.st_uid != os.geteuid()
            or stat.S_IMODE(metadata.st_mode) != 0o600
            or metadata.st_size != authority.hpo_input_size_bytes
        ):
            raise ValueError
        content = read_bounded_bytes(
            path,
            maximum_bytes=MAX_HPO_LIST_BYTES,
            reference_root=public_qualification_hpo_root(authority),
            expected_size_bytes=authority.hpo_input_size_bytes,
        )
        if hashlib.sha256(content).hexdigest() != authority.hpo_input_sha256:
            raise ValueError
    except (OSError, TypeError, ValueError):
        raise WgsPublicQualificationError(
            "public qualification HPO input is unsafe or mismatched"
        ) from None


def public_qualification_hpo_root(
    authority: PublicQualificationAuthority,
) -> Path:
    """Return the authenticated HPO path's exact safe-reader containment root."""

    _public_path(authority.hpo_input_path, roots=(_PUBLIC_HPO_ROOT,))
    return _PUBLIC_HPO_ROOT


def require_public_qualification_execution(
    authority: PublicQualificationAuthority,
    *,
    original_input_path: str | Path,
    hpo_input_path: str | Path | None,
    output_path: str | Path,
    input_sha256: str,
    input_size_bytes: int,
    hpo_input_sha256: str,
    execution_config: BoundedWgsExecutionConfig,
    executable_stack_sha256: str,
    runtime_image_digest: str,
    policy_sha256: str,
    population_manifest_sha256: str,
    clinical_manifest_sha256: str,
    exception_manifest_sha256: str,
    precision_manifest_sha256: str,
    preparation_run_options_sha256: str | None = None,
) -> None:
    """Match an invocation to one exact public authority without filesystem reads."""

    try:
        output = Path(output_path)
        authorized_output = Path(authority.output_root)
        hpo_matches = hpo_input_path is None or (
            Path(hpo_input_path).as_posix() == authority.hpo_input_path
        )
        if (
            Path(original_input_path).as_posix() != authority.input_path
            or not hpo_matches
            or not output.is_absolute()
            or output == authorized_output
            or not output.is_relative_to(authorized_output)
            or any(part in _FORBIDDEN_PARTS for part in output.parts)
            or input_sha256 != authority.input_sha256
            or input_size_bytes != authority.input_size_bytes
            or hpo_input_sha256 != authority.hpo_input_sha256
            or hashlib.sha256(
                bounded_wgs_execution_config_as_json(execution_config)
            ).hexdigest() != authority.execution_config_sha256
            or executable_stack_sha256 != authority.executable_stack_sha256
            or runtime_image_digest != authority.runtime_image_digest
            or policy_sha256 != authority.binding("common_allele_policy")
            or population_manifest_sha256
            != authority.binding("gnomad_joint_manifest")
            or clinical_manifest_sha256 != authority.binding("clinvar_manifest")
            or exception_manifest_sha256
            != authority.binding("common_allele_exception_manifest")
            or precision_manifest_sha256
            != authority.binding("gnomad_projection_precision_report")
            or preparation_run_options_sha256 is not None
            and preparation_run_options_sha256
            != authority.preparation_run_options_sha256
        ):
            raise ValueError
    except (OSError, TypeError, ValueError):
        raise WgsPublicQualificationError(
            "invocation differs from public qualification authority"
        ) from None


def require_public_qualification_release(
    authority: PublicQualificationAuthority,
    release: BaseRelease,
) -> None:
    """Permit only an exact candidate-only release for public qualification."""

    try:
        bindings = release.policy_bindings
        authority_sha256 = authority.sha256
        required_resources = {
            value.sha256 for value in authority.resource_bindings
        } | {
            value.sha256 for value in authority.policy_bindings
            if value.role.endswith("manifest")
            or value.role == "gnomad_projection_precision_report"
        }
        if (
            release.authorization_status != "candidate_only"
            or release.candidate_sha256 != authority.engineering_candidate_sha256
            or release.consolidated_acceptance_sha256
            != authority.authorization_receipt_sha256
            or authority_sha256 not in release.accepted_result_sha256s
            or authority.plan_sha256 not in release.accepted_result_sha256s
            or authority.authorization_receipt_sha256
            not in release.accepted_result_sha256s
            or release.validation_plan_sha256 != authority.plan_sha256
            or release.executable_stack_sha256
            != authority.executable_stack_sha256
            or release.runtime_image_digest != authority.runtime_image_digest
            or not required_resources <= set(release.resource_manifest_sha256s)
            or bindings.common_allele_policy_sha256
            != authority.binding("common_allele_policy")
            or bindings.common_allele_exception_manifest_sha256
            != authority.binding("common_allele_exception_manifest")
            or bindings.common_allele_public_index_manifest_sha256
            != authority.binding("common_allele_public_index_manifest")
            or bindings.gnomad_projection_precision_sha256
            != authority.binding("gnomad_projection_precision_report")
            or bindings.bounded_wgs_execution_sha256
            != authority.execution_config_sha256
        ):
            raise ValueError
    except (AttributeError, KeyError, TypeError, ValueError):
        raise WgsPublicQualificationError(
            "release does not authorize public qualification"
        ) from None
