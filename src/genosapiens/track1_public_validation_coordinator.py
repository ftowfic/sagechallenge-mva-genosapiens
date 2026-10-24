"""Evaluator-side coordinator for the bounded public Track 1 experiment.

The scientific runner receives only :class:`RankerCaseView` values.  This outer
coordinator retains the public synthetic truth roster, records complete attempts,
selects the finite development family, freezes a candidate, and only then releases
held-out runner views.
"""

from __future__ import annotations

from concurrent.futures import Executor, ThreadPoolExecutor, FIRST_COMPLETED, wait
from dataclasses import asdict, dataclass, fields, replace
from fractions import Fraction
import hashlib
import json
import os
from pathlib import Path
import re
import stat
import tempfile
from threading import Event
from typing import Callable, Literal
from uuid import uuid4

from genosapiens.track1_public_validation import (
    REFERENCE_MANIFEST_SHA256,
    RUNTIME_IMAGE,
    PublicPolicyFamily,
    PublicPolicyRuntime,
    PublicResourceContext,
    PublicTrack1ValidationError,
    execute_prepared_public_case,
    load_public_policy_family,
    load_public_resource_context,
    prepare_public_case,
    _authenticated_policy_loci,
    _population_template,
)
from genosapiens.track1_roster import (
    RosterCase, Track1Roster, load_submission_acceptance_profile,
    load_submission_acceptance_roster,
    require_profile_acceptance_binding,
)
from genosapiens.track1_release import BaseRelease, freeze_base_release, evidence_interface_bytes
from genosapiens.track1_public_worker import run_isolated_public_case
from genosapiens.track1_validation_plan import (
    AcceptanceCriterion,
    RankerCaseView,
    ValidationPlan,
    validation_plan_as_json,
)
from genosapiens.track1_validation_runner import (
    REQUIRED_ACCEPTANCE_RESULTS,
    AcceptanceResultReference,
    CaseEvaluation,
    ConsolidatedAcceptance,
    DevelopmentExperiment,
    PreHeldOutCandidate,
    PolicyCandidate,
    ResourceMeasurements,
    ValidationAttempt,
    PREPARATION_STAGE_PREFIX, PREPARATION_STAGES,
    EXECUTION_STAGE_PREFIX, EXECUTION_STAGES,
    authenticated_held_out_ranker_view,
    complete_development_experiment,
    development_ranker_view,
    evaluate_attempts,
    evaluation_metrics,
    model_coverage_complete,
    freeze_preheldout_candidate,
    consolidate_public_acceptance,
    record_execution_attempt,
)
from genosapiens.wgs_public_execution import OUTPUT_ROLES


SCHEMA = "genosapiens.track1-public-candidate-validation"
STOPPED_SCHEMA = "genosapiens.track1-public-validation-stopped-outcome"
GATE_FAILURES = frozenset({
    "execution-failed", "replay-failed", "required-model-coverage-missing",
    "positive-export-incomplete", "population-probe-mismatch", "singleton-violation",
})
UNSTARTED_REASONS = frozenset({
    "stage-a-gate-failed", "development-gate-failed", "held-out-gate-failed", "owned-stop", "interrupted",
})
POPULATION_STAGE_SCHEDULE = (
    ('A', 'development', ('development-population-model-state-01',
        'development-population-model-state-02', 'development-population-model-state-03',
        'development-source-backed-singleton')),
    ('B', 'development', ('development-mixed-class', 'development-repeat')),
    ('C', 'held_out', ('held-out-population-model-state-01', 'held-out-source-backed-singleton')),
)
FINAL_SCHEMA = "genosapiens.track1-public-validation-and-release"
CURRENT_ENGINEERING_RESULT_PATH = (
    "artifacts/validation/track1-qualified-current-engineering-2026-09-25-context-completion.json"
)
CURRENT_SECURITY_RESULT_PATH = (
    "artifacts/supply-chain/track1-qualified-current-security-2026-09-25-context-completion.json"
)
Progress = Callable[[ValidationAttempt], None]


class PublicTrack1CoordinatorError(RuntimeError):
    """Static evaluator/coordinator failure without source-derived payloads."""


def _bytes(value: object) -> bytes:
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


def _sha256(content: bytes) -> str:
    return hashlib.sha256(content).hexdigest()


def _fraction(value: Fraction) -> str:
    return (
        str(value.numerator)
        if value.denominator == 1
        else f"{value.numerator}/{value.denominator}"
    )


def _floor_decimal(value: Fraction, places: int = 6) -> str:
    """Conservatively project an observed exact ratio into the decimal schema."""
    scale = 10**places
    scaled = value.numerator * scale // value.denominator
    whole, remainder = divmod(scaled, scale)
    if not remainder:
        return str(whole)
    return f"{whole}.{remainder:0{places}d}".rstrip("0")


def _ceil_decimal(value: Fraction, places: int = 6) -> str:
    scale = 10**places
    scaled = -(-value.numerator * scale // value.denominator)
    whole, remainder = divmod(scaled, scale)
    if not remainder:
        return str(whole)
    return f"{whole}.{remainder:0{places}d}".rstrip("0")


def executable_stack_bytes(repository_root: str | Path) -> bytes:
    """Bind all Track 1 Python/config bytes, not merely a mutable revision label."""
    root = Path(repository_root).resolve(strict=True)
    paths = tuple(sorted(
        {
            root / 'track1',
            root / 'scripts/resolve_track1_resume.py',
            root / 'scripts/portable_paths.py',
            root / 'scripts/ensure_track1_installation.py',
            root / 'scripts/build_track1_installation.py',
            root / 'config/track1/runtime-image.json',
            root / 'docs/track1-operator-command.md',
            root / 'docs/track1-installation.md',
            *(
                path
                for path in (root / "src/genosapiens").glob("*.py")
                if path.is_file() and not path.is_symlink()
            ),
            root / "requirements/track1-pipeline.lock",
            root / "config/track1/public-ranking-development-policy-family.json",
            root / "config/track1/public-population-development-derivation.json",
            root / "config/public-resources/grch38-model-loci.json",
            root / "config/track1/submission-acceptance-profile.json",
            root / "config/track1/challenge-evaluator.json",
            root / "config/track1/epcr-policy.json",
            root / "config/track1/epcr-numerical-suitability.json",
            root / "config/track1/common-allele-prefilter-v1.json",
            root / "config/track1/early-scope-policy.json",
            root / "config/track1/vcf-evidence-policy-v1.json",
            root / "config/track1/analysis-eligibility-v1.json",
            root / "config/track1/analysis-profiles.json",
            root / "config/track1/protected-population-work-filter-v1.json",
            root / "config/track1/focused-work-contract-v1.json",
            *(root / "config/track1").glob("regulatory-index-pin.json"),
            root / "config/track1/context-writer-runtime-v1.json",
            root / "config/track1/handoff-contract.json",
            root / "config/track1/handoff-pin.json",
            root / "schemas/track1/mechanism-ready-handoff.schema.json",
            root / "config/track1/bounded-wgs-execution-v5.json",
            root / "config/public-resources/track1-common-allele-exceptions-v1.json",
            root / "config/public-resources/track1-common-allele-exception-curation-v1.json",
            root / "config/public-resources/gnomad-joint-faf95-precision-v1.json",
            root / "scripts/run_track1_public_validation.py",
            root / "scripts/launch-track1-public-validation.sh",
            root / "scripts/build_track1_common_allele_exception_sources.py",
            root / "scripts/build_track1_common_allele_exception_bundle.py",
            root / "scripts/build_common_allele_public_indexes.py",
            root / "scripts/build_public_wgs_unprefixed_derivative.py",
            root / "scripts/qualify_gnomad_faf95_precision.py",
            root / "scripts/prepare_track1_operator_bounded_wgs.py",
            root / "scripts/operator_compatibility_repair/prepare_bounded_case.py",
            root / "scripts/operator_compatibility_repair/analyze_bounded_case.py",
            root / "scripts/operator_compatibility_repair/archive_bounded_resume.py",
            root / "scripts/operator_compatibility_repair/archive_failed_preparation.py",
            root / "scripts/operator_compatibility_repair/check_bounded_host_capacity.py",
            root / "scripts/operator_compatibility_repair/check_stage_prerequisite.py",
            root / "scripts/operator_compatibility_repair/finalize_bounded_runtime.py",
            root / "scripts/operator_compatibility_repair/run_bounded_stage.py",
            root / "scripts/operator_compatibility_repair/seal_bounded_stage_evidence.py",
            *(root / "config/public-resources").glob("*.json"),
            *(root / "config/reference-identities").glob("*.json"),
        },
        key=lambda path: path.relative_to(root).as_posix(),
    ))
    if any(not path.is_file() or path.is_symlink() for path in paths):
        raise PublicTrack1CoordinatorError("executable stack is incomplete")
    return _bytes({
        "schema": "genosapiens.track1-executable-stack",
        "runtime_image": RUNTIME_IMAGE,
        "files": tuple({
            "path": path.relative_to(root).as_posix(),
            "size_bytes": path.stat().st_size,
            "sha256": _sha256(path.read_bytes()),
        } for path in paths),
    })


@dataclass(frozen=True, slots=True)
class PublicCandidateValidation:
    executable_stack_sha256: str
    development_experiment: DevelopmentExperiment
    development_derivation: dict
    candidate: PreHeldOutCandidate | None
    held_out_attempts: tuple[ValidationAttempt, ...]
    held_out_evaluations: tuple
    held_out_metrics: tuple[tuple[str, str], ...]
    determinism_replays: tuple[tuple[str, str, str, str], ...]
    status: Literal[
        "no-development-candidate",
        "held-out-complete",
    ]
    schema: Literal[
        "genosapiens.track1-public-candidate-validation"
    ] = SCHEMA
    classification: Literal["public"] = "public"
    validation_only: Literal[True] = True
    controlled_use_authorized: Literal[False] = False
    carried_forward: dict | None = None

    def __post_init__(self) -> None:
        candidate_present = self.candidate is not None
        try:
            canonical_metrics = tuple(sorted(dict(self.held_out_metrics).items()))
        except (TypeError, ValueError):
            raise PublicTrack1CoordinatorError(
                "invalid public candidate validation"
            ) from None
        if (
            self.schema != SCHEMA
            or self.classification != "public"
            or self.validation_only is not True
            or self.controlled_use_authorized is not False
            or self.status not in {
                "no-development-candidate", "held-out-complete",
            }
            or len(self.executable_stack_sha256) != 64
            or any(character not in "0123456789abcdef" for character in self.executable_stack_sha256)
            or type(self.development_experiment) is not DevelopmentExperiment
            or type(self.development_derivation) is not dict
            or any(type(item) is not ValidationAttempt for item in self.held_out_attempts)
            or any(type(item) is not CaseEvaluation for item in self.held_out_evaluations)
            or self.held_out_metrics != canonical_metrics
            or any(
                type(name) is not str or not name or type(value) is not str
                for name, value in self.held_out_metrics
            )
            or any(
                type(item) is not tuple
                or len(item) != 4
                or item[0] not in {"development", "held_out"}
                or any(type(value) is not str or not value for value in item[:3])
                or len(item[3]) != 64
                or any(character not in "0123456789abcdef" for character in item[3])
                for item in self.determinism_replays
            )
            or candidate_present != (self.status == "held-out-complete")
            or candidate_present
            and (
                type(self.candidate) is not PreHeldOutCandidate
                or not self.development_derivation
                or not self.held_out_attempts
                or not self.held_out_evaluations
                or not self.held_out_metrics
                or self.candidate.executable_stack_sha256
                != self.executable_stack_sha256
                or self.candidate.development_experiment_sha256
                != _sha256(self.development_experiment.as_bytes())
            )
            or not candidate_present
            and (
                self.development_derivation
                or self.held_out_attempts
                or self.held_out_evaluations
                or self.held_out_metrics
            )
        ):
            raise PublicTrack1CoordinatorError("invalid public candidate validation")

    def as_bytes(self) -> bytes:
        return _bytes(asdict(self))

    @property
    def sha256(self) -> str:
        return _sha256(self.as_bytes())


@dataclass(frozen=True, slots=True)
class PublicValidationStoppedOutcome:
    execution_context_sha256: str
    executable_stack_sha256: str
    profile_sha256: str
    stage_schedule: tuple[tuple[str, str, tuple[str, ...]], ...]
    declared_case_ids: tuple[str, ...]
    stop_stage: Literal['A', 'B', 'C']
    stop_kind: Literal['gate-failure', 'interrupted', 'owned-stop']
    first_failed_case_id: str | None
    first_failed_check: str | None
    attempted_case_ids: tuple[str, ...]
    replay_started_case_ids: tuple[str, ...]
    replay_completed_case_ids: tuple[str, ...]
    unstarted_cases: tuple[tuple[str, str], ...]
    unstarted_replays: tuple[tuple[str, str], ...]
    attempt_artifacts: tuple[dict, ...]
    replay_started_artifacts: tuple[dict, ...]
    replay_artifacts: tuple[dict, ...]
    worker_evidence_artifacts: tuple[dict, ...]
    gate_checks: tuple[tuple[str, tuple[tuple[str, bool | None], ...]], ...]
    candidate: PreHeldOutCandidate | None
    released_held_out_views: tuple[dict, ...]
    schema: Literal['genosapiens.track1-public-validation-stopped-outcome'] = STOPPED_SCHEMA
    experiment_id: str = 'track1-submission-acceptance-population-models-20260917'
    runtime_image: str = RUNTIME_IMAGE
    classification: Literal['public'] = 'public'
    controlled_use_authorized: Literal[False] = False
    declared_case_count: Literal[8] = 8
    declared_development_case_count: Literal[6] = 6
    declared_held_out_case_count: Literal[2] = 2
    status: Literal['validation-stopped'] = 'validation-stopped'
    carried_forward: dict | None = None

    def __post_init__(self):
        schedule = tuple(name for _, _, names in self.stage_schedule for name in names)
        declared = set(self.declared_case_ids)
        attempted = set(self.attempted_case_ids)
        started = set(self.replay_started_case_ids)
        completed = set(self.replay_completed_case_ids)
        stage_index = next((index for index, row in enumerate(self.stage_schedule)
                            if row[0] == self.stop_stage), -1)
        admitted = {name for _, _, names in self.stage_schedule[:stage_index + 1] for name in names}
        prior = {name for _, _, names in self.stage_schedule[:stage_index] for name in names}
        checks = dict(self.gate_checks)
        reason = {'A':'stage-a-gate-failed','B':'development-gate-failed','C':'held-out-gate-failed'}.get(self.stop_stage)
        if self.stop_kind != 'gate-failure':
            reason = self.stop_kind
        if (self.schema != STOPPED_SCHEMA or self.classification != 'public'
                or self.experiment_id != 'track1-submission-acceptance-population-models-20260917'
                or self.runtime_image != RUNTIME_IMAGE or self.controlled_use_authorized is not False
                or (self.declared_case_count,self.declared_development_case_count,self.declared_held_out_case_count) != (8,6,2)
                or self.status != 'validation-stopped' or len(schedule) != 8 or len(declared) != 8
                or self.declared_case_ids != schedule
                or self.stage_schedule != POPULATION_STAGE_SCHEDULE
                or tuple((stage,split,len(names)) for stage,split,names in self.stage_schedule)
                   != (('A','development',4),('B','development',2),('C','held_out',2))
                or self.stop_stage not in {'A','B','C'} or self.stop_kind not in {'gate-failure','interrupted','owned-stop'}
                or any(re.fullmatch(r'[0-9a-f]{64}', value) is None for value in (
                    self.execution_context_sha256,self.executable_stack_sha256,self.profile_sha256))
                or not completed <= started <= attempted <= declared
                or not attempted <= admitted or not prior <= completed
                or tuple(name for name in schedule if name in attempted) != self.attempted_case_ids
                or tuple(name for name in schedule if name in started) != self.replay_started_case_ids
                or tuple(name for name in schedule if name in completed) != self.replay_completed_case_ids
                or tuple(name for name, _ in self.gate_checks) != self.attempted_case_ids
                or any(len(rows) != len(GATE_FAILURES) or set(dict(rows)) != GATE_FAILURES
                       or any(passed is not None and type(passed) is not bool for _, passed in rows)
                       for _, rows in self.gate_checks)
                or any(not all(passed for _, passed in checks[name]) for name in prior)
                or any(len(values) != len(set(values)) for values in (
                    self.attempted_case_ids,self.replay_started_case_ids,self.replay_completed_case_ids))
                or self.unstarted_cases != tuple((name,reason) for name in schedule if name not in attempted)
                or self.unstarted_replays != tuple((name,reason) for name in schedule if name not in started)
                or reason not in UNSTARTED_REASONS
                or self.stop_kind == 'gate-failure' and (
                    self.first_failed_check not in GATE_FAILURES or self.first_failed_case_id not in attempted
                    or self.first_failed_case_id not in self.stage_schedule[stage_index][2]
                    or dict(checks[self.first_failed_case_id]).get(self.first_failed_check) is not False)
                or self.stop_kind != 'gate-failure' and (
                    self.first_failed_case_id is not None or self.first_failed_check is not None)
                or self.stop_stage != 'C' and (self.candidate is not None or self.released_held_out_views)
                or self.stop_stage == 'C' and type(self.candidate) is not PreHeldOutCandidate
                or self.stop_stage == 'C' and (
                    self.candidate.executable_stack_sha256 != self.executable_stack_sha256
                    or self.candidate.roster_sha256 != self.profile_sha256
                    or len(self.released_held_out_views) != 2)
                or len(self.attempt_artifacts) != len(attempted)
                or len(self.replay_artifacts) != len(completed)
                or self.stop_kind == 'gate-failure' and completed == declared):
            raise PublicTrack1CoordinatorError('invalid stopped validation outcome')
        for collection in (self.attempt_artifacts,self.replay_started_artifacts,self.replay_artifacts,
                           self.worker_evidence_artifacts,self.released_held_out_views):
            paths = []
            for binding in collection:
                if (type(binding) is not dict or set(binding) != {'path','sha256','size_bytes'}
                        or type(binding['path']) is not str or Path(binding['path']).is_absolute()
                        or '..' in Path(binding['path']).parts or binding['path'] in paths
                        or type(binding['size_bytes']) is not int or binding['size_bytes'] < 0
                        or re.fullmatch(r'[0-9a-f]{64}',binding['sha256']) is None):
                    raise PublicTrack1CoordinatorError('invalid stopped validation artifact')
                paths.append(binding['path'])

    def as_bytes(self):
        return _bytes({**asdict(self), 'accounting': {
            'attempted_case_count':len(self.attempted_case_ids),
            'replay_started_count':len(self.replay_started_case_ids),
            'replay_completed_count':len(self.replay_completed_case_ids),
            'unstarted_case_count':len(self.unstarted_cases),
            'unstarted_replay_count':len(self.unstarted_replays),
            'unfinished_started_replay_case_ids':tuple(name for name in self.replay_started_case_ids
                                                     if name not in self.replay_completed_case_ids),
            'attempted_subset_diagnostics_satisfy_full_acceptance':False,
        }})

    @property
    def sha256(self):
        return _sha256(self.as_bytes())


@dataclass(frozen=True, slots=True)
class FinalizedPublicValidation:
    candidate_validation: PublicCandidateValidation
    acceptance: ConsolidatedAcceptance
    base_release: BaseRelease | None
    status: Literal["accepted-and-frozen", "rejected"]
    schema: Literal["genosapiens.track1-public-validation-and-release"] = FINAL_SCHEMA
    classification: Literal["public"] = "public"

    def __post_init__(self) -> None:
        if (
            self.schema != FINAL_SCHEMA
            or self.classification != "public"
            or self.status not in {"accepted-and-frozen", "rejected"}
            or type(self.candidate_validation) is not PublicCandidateValidation
            or type(self.acceptance) is not ConsolidatedAcceptance
            or self.base_release is not None
            and type(self.base_release) is not BaseRelease
            or (self.base_release is not None) != (self.status == "accepted-and-frozen")
            or (self.acceptance.status == "accepted")
            != (self.status == "accepted-and-frozen")
        ):
            raise PublicTrack1CoordinatorError("invalid finalized public validation")

    def as_bytes(self) -> bytes:
        return _bytes(asdict(self))

    @property
    def sha256(self) -> str:
        return _sha256(self.as_bytes())


def _development_derivation(
    experiment: DevelopmentExperiment,
    repository_root: Path,
) -> tuple[dict, bytes, tuple[AcceptanceCriterion, ...]]:
    if experiment.selected_policy_id is None:
        raise PublicTrack1CoordinatorError("development produced no selected policy")
    evaluations = tuple(
        item
        for item in experiment.evaluations
        if item.policy_id == experiment.selected_policy_id
    )
    metrics = evaluation_metrics(evaluations)
    profile, profile_bytes = load_submission_acceptance_profile(repository_root)
    if experiment.roster_sha256 != _sha256(profile_bytes):
        raise PublicTrack1CoordinatorError("development scope differs from the declared profile")
    if not _positive_exports_complete(profile, evaluations, "development"):
        raise PublicTrack1CoordinatorError("development positive exact-export contract failed")
    threshold_specs = tuple(tuple(row) for row in profile["fixed_metric_criteria"])
    document = {
        "schema": "genosapiens.track1-development-derived-criteria",
        "classification": "public",
        "validation_only": True,
        "development_experiment_sha256": _sha256(experiment.as_bytes()),
        "selected_policy_id": experiment.selected_policy_id,
        "observed_exact_metrics": tuple(sorted(
            (name, _fraction(value)) for name, value in metrics.items()
        )),
        "held_out_thresholds": threshold_specs,
        "threshold_derivation": "fixed submission-profile criteria declared before outcomes",
        "submission_acceptance_profile_sha256": _sha256(profile_bytes),
        "clinical_calibration_claimed": False,
    }
    content = _bytes(document)
    digest = _sha256(content)
    criteria = tuple(
        AcceptanceCriterion(name, comparator, threshold, digest)
        for name, comparator, threshold in threshold_specs
    )
    return document, content, criteria


def _positive_exports_complete(profile, evaluations, split):
    values = {value.case_id: value for value in evaluations if value.split == split}
    required = [row for row in profile["cases"] if row["split"] == split and row["positive_exact_export"]]
    limit = profile["positive_export_contract"]["eligible_rank_at_most"]
    return bool(required) and all(
        row["case_id"] in values
        and values[row["case_id"]].status == "succeeded"
        and values[row["case_id"]].target_kind == row["target_kind"]
        and values[row["case_id"]].exact_target_generated
        and values[row["case_id"]].exact_target_eligible_rank is not None
        and values[row["case_id"]].exact_target_eligible_rank <= limit
        for row in required
    )


def _candidate(
    experiment: DevelopmentExperiment,
    *,
    repository_root: Path,
    resources: PublicResourceContext,
    policies: PublicPolicyFamily,
    derivation_bytes: bytes,
    criteria: tuple[AcceptanceCriterion, ...],
    stack_sha256: str,
    recovery_resolver,
) -> PreHeldOutCandidate:
    from genosapiens.track1_public_validation import authenticate_population_development_derivation
    authenticate_population_development_derivation(experiment, policies, repository_root, recovery_resolver=recovery_resolver)
    evaluator_pin_bytes = (
        repository_root / "config/track1/challenge-evaluator.json"
    ).read_bytes()
    evaluator_bytes = _bytes({
        "pin_sha256": _sha256(evaluator_pin_bytes),
        "source_sha256": _sha256(resources.evaluator_source_path.read_bytes()),
    })
    numerical_bytes = _bytes({
        "epcr_policy_sha256": policies.epcr_policy_sha256,
        "epcr_suitability_sha256": _sha256(policies.epcr_derivation_bytes),
    })
    runtime_digest = "sha256:" + RUNTIME_IMAGE.rsplit("@sha256:", 1)[1]
    return freeze_preheldout_candidate(
        experiment,
        derivation_bytes=derivation_bytes,
        numerical_policy_evidence_sha256=_sha256(numerical_bytes),
        evaluator_sha256=_sha256(evaluator_bytes),
        executable_stack_sha256=stack_sha256,
        runtime_image_digest=runtime_digest,
        resource_manifest_sha256s=tuple(sorted({
            *(item.descriptor_sha256 for item in resources.common_runner_resources),
            _sha256((repository_root / 'config/track1/analysis-profiles.json').read_bytes()),
        })),
        acceptance_criteria=criteria,
    )


def _failed_attempt(case, view, policy: PolicyCandidate, stage: str, *, error=None) -> ValidationAttempt:
    statuses = ((stage, "failed"),)
    if stage in {"public-case-preparation", "public-resource-authentication"}:
        safe_stage = (
            error.preparation_stage if type(error) is PublicTrack1ValidationError else None
        )
        if type(safe_stage) is not str or safe_stage not in PREPARATION_STAGES:
            safe_stage = "preparation"
        statuses += ((PREPARATION_STAGE_PREFIX + safe_stage, "failed"),)
    elif stage == "public-scientific-execution":
        safe_stage = (
            error.execution_stage if type(error) is PublicTrack1ValidationError else None
        )
        if type(safe_stage) is not str or safe_stage not in EXECUTION_STAGES:
            safe_stage = "execution"
        statuses += ((EXECUTION_STAGE_PREFIX + safe_stage, "failed"),)
    return ValidationAttempt(
        case.case_id,
        case.split,
        case.scope,
        policy.policy_id,
        policy.bindings.sha256,
        1,
        "failed",
        None,
        tuple(item.sha256 for item in view.inputs),
        tuple(item.descriptor_sha256 for item in view.resources),
        statuses,
        None,
        (),
        None,
        None,
        None,
        None,
        ResourceMeasurements(None, None, None, None),
    )


def _run_development_case(
    case: RosterCase,
    view: RankerCaseView,
    resources: PublicResourceContext,
    policies: PublicPolicyFamily,
    checkpoint_root: Path | None = None,
) -> tuple[tuple[ValidationAttempt, ...], tuple[tuple[str, str, str, str], ...]]:
    attempts = []
    replays = []
    try:
        prepared = prepare_public_case(
            case.case_id,
            view,
            resources=resources,
            policies=policies,
        )
    except Exception as error:
        return (
            tuple(
                _failed_attempt(
                    case,
                    view,
                    runtime.candidate,
                    "public-case-preparation",
                    error=error,
                )
                for runtime in policies.runtimes
            ),
            (),
        )
    for runtime in policies.runtimes:
        try:
            execution, measurements = execute_prepared_public_case(
                prepared,
                runtime,
                resources=resources,
                policies=policies,
            )
        except Exception as error:
            attempts.append(_failed_attempt(
                case,
                view,
                runtime.candidate,
                "public-scientific-execution",
                error=error,
            ))
            if checkpoint_root is not None:
                _persist_case_attempts(checkpoint_root, tuple(attempts[-1:]), (), None)
            continue
        if checkpoint_root is not None:
            _persist_execution(checkpoint_root, execution)
        attempts.append(record_execution_attempt(
            execution,
            roster_case=case,
            policy=runtime.candidate,
            attempt_number=1,
            measurements=measurements,
            population_selector=_population_selector(policies),
        ))
        if checkpoint_root is not None:
            _persist_case_attempts(checkpoint_root, tuple(attempts[-1:]), (), None)
        # Keep the immutable comparison bytes, rather than two scientific
        # object graphs, alive while the matching replay is constructed.
        expected_execution_bytes = execution.as_bytes()
        execution_sha256 = _sha256(expected_execution_bytes)
        del execution
        if checkpoint_root is not None:
            _persist_replay_started(checkpoint_root, attempts[-1])
        replayed, _replayed_measurements = execute_prepared_public_case(
            prepared,
            runtime,
            resources=resources,
            policies=policies,
        )
        if replayed.as_bytes() != expected_execution_bytes:
            raise PublicTrack1CoordinatorError(
                "public scientific execution is not deterministic"
            )
        del replayed
        replays.append((
            "development",
            case.case_id,
            runtime.candidate.policy_id,
            execution_sha256,
        ))
        if checkpoint_root is not None:
            _persist_case_attempts(checkpoint_root, (), tuple(replays[-1:]), None)
    return tuple(attempts), tuple(replays)


def _run_held_out_case(
    case: RosterCase,
    view: RankerCaseView,
    runtime: PublicPolicyRuntime,
    resources: PublicResourceContext,
    policies: PublicPolicyFamily,
    checkpoint_root: Path | None = None,
) -> tuple[ValidationAttempt, tuple[str, str, str, str] | None]:
    stage = "public-case-preparation"
    try:
        prepared = prepare_public_case(
            case.case_id,
            view,
            resources=resources,
            policies=policies,
        )
        stage = "public-scientific-execution"
        execution, measurements = execute_prepared_public_case(
            prepared,
            runtime,
            resources=resources,
            policies=policies,
        )
    except Exception as error:
        return _failed_attempt(
            case,
            view,
            runtime.candidate,
            stage, error=error,
        ), None
    attempt = record_execution_attempt(
        execution,
        roster_case=case,
        policy=runtime.candidate,
        attempt_number=1,
        measurements=measurements,
        population_selector=_population_selector(policies),
    )
    if checkpoint_root is not None:
        _persist_execution(checkpoint_root, execution)
        _persist_case_attempts(checkpoint_root, (attempt,), (), None)
    expected_execution_bytes = execution.as_bytes()
    execution_sha256 = _sha256(expected_execution_bytes)
    del execution
    if checkpoint_root is not None:
        _persist_replay_started(checkpoint_root, attempt)
    replayed, _replayed_measurements = execute_prepared_public_case(
        prepared,
        runtime,
        resources=resources,
        policies=policies,
    )
    if replayed.as_bytes() != expected_execution_bytes:
        raise PublicTrack1CoordinatorError(
            "public scientific execution is not deterministic"
        )
    del replayed
    replay = (
        "held_out",
        case.case_id,
        runtime.candidate.policy_id,
        execution_sha256,
    )
    if checkpoint_root is not None:
        _persist_case_attempts(checkpoint_root, (), (replay,), None)
    return attempt, replay


def _parallel_results(
    executor: Executor,
    worker,
    jobs,
    completed: Callable[[int, object], None] | None = None,
    received: Callable[[int, object], None] | None = None,
    failed: Callable[[int], object] | None = None,
    *,
    concurrency: int | None = None,
    gate: Callable[[int, object], bool] | None = None,
    interrupted: Callable[[], None] | None = None,
):
    jobs = tuple(jobs)
    limit = concurrency if concurrency is not None else getattr(executor, '_max_workers', 1)
    if type(limit) is not int or limit < 1:
        raise PublicTrack1CoordinatorError("invalid public validation concurrency")
    futures, next_index, stopped, interruption = {}, 0, False, False
    results = [None] * len(jobs)
    def submit_available():
        nonlocal next_index
        while not stopped and len(futures) < limit and next_index < len(jobs):
            index = next_index
            next_index += 1
            futures[executor.submit(worker, jobs[index])] = index
    submit_available()
    while futures:
        try:
            done, _pending = wait(tuple(futures), return_when=FIRST_COMPLETED)
        except KeyboardInterrupt:
            stopped, interruption = True, True
            if interrupted is not None:interrupted()
            continue
        for future in sorted(done, key=lambda item: futures[item]):
            index = futures.pop(future)
            try:
                results[index] = future.result()
            except KeyboardInterrupt:
                stopped, interruption = True, True
                if interrupted is not None:interrupted()
                if failed is None:raise
                results[index] = failed(index)
            except Exception:
                if failed is None:
                    raise PublicTrack1CoordinatorError("public validation worker failed") from None
                results[index] = failed(index)
            if received is not None:
                try:
                    received(index, results[index])
                except KeyboardInterrupt:
                    stopped, interruption = True, True
                    if interrupted is not None:interrupted()
                    # The worker returned real results. Finish their durable
                    # publication before draining the remaining owned workers.
                    received(index, results[index])
            if gate is not None and not gate(index, results[index]):
                stopped = True
            if completed is not None:
                try:completed(index, results[index])
                except KeyboardInterrupt:
                    stopped, interruption = True, True
                    if interrupted is not None:interrupted()
        submit_available()
    if interruption:raise KeyboardInterrupt
    return tuple(results)


def _case_gate_checks(root, roster, case, attempts, replays, profile, *, journal=None, recovery_resolver=None):
    if len(attempts) != 1:
        return (('execution-failed',False),)
    attempt = attempts[0]
    stages = dict(attempt.stage_statuses)
    # Successful attempts report the bounded execution manifest's output roles.
    # Use its contract directly so producer and acceptance cannot drift apart.
    success = (attempt.status == 'succeeded' and attempt.source_inventory_status == 'complete'
               and attempt.required_enumerations is not None
               and attempt.required_enumerations == attempt.completed_enumerations
               and stages.keys() == OUTPUT_ROLES
               and all(status in {'succeeded','succeeded-empty'} for status in stages.values()))
    if journal is not None:
        for path in (journal/'worker-lifecycle').glob('*.json'):
            lifecycle=_checkpoint_json(_read_checkpoint(journal,path.relative_to(journal).as_posix()))
            if lifecycle.get('case_id')==case.case_id:
                success = success and lifecycle.get('exit_code')==0 and lifecycle.get('container_operation_failed') is False
    checks = [('execution-failed',success)]
    checks.append(('replay-failed', len(replays)==1 and replays[0] == (
        case.split,case.case_id,attempt.policy_id,attempt.execution_sha256)))
    coverage = dict(attempt.model_coverage_counts)
    checks.append(('singleton-violation',coverage.get('unresolved-partner-singleton-violation',0)==0))
    declaration = next(row for row in profile['cases'] if row['case_id']==case.case_id)
    evaluation = evaluate_attempts(replace(roster,cases=(case,)),(attempt,),
        split=case.split,policy_id=attempt.policy_id,recovery_resolver=recovery_resolver)[0]
    coverage_complete = (model_coverage_complete((evaluation,)) if case.required_model_coverage else
        coverage.get('clinical-annotation-completed',0)==1
        and coverage.get('population-annotation-completed',0)==1)
    checks.append(('required-model-coverage-missing',coverage_complete))
    checks.append(('positive-export-incomplete',not declaration['positive_exact_export'] or (
        evaluation.exact_target_generated and evaluation.exact_target_eligible_rank is not None
        and evaluation.exact_target_eligible_rank<=profile['positive_export_contract']['eligible_rank_at_most'])))
    expected = []
    if '-population-model-state-' in case.case_id:
        if case.split == 'development':
            selection = json.loads((root/'artifacts/validation/track1-population-development-selection-2026-09-17.json').read_bytes())
            qualification = json.loads((root/'artifacts/validation/track1-population-development-qualification-2026-09-17.json').read_bytes())
            identities = selection['selected_low_identity_sha256s'] if case.case_id.endswith('01') else (
                selection['probes'][0]['identity_sha256'],selection['selected_high_identity_sha256'])
            expected = [(row['identity_sha256'],Fraction(row['frequency']),row['allele_number'])
                        for row in qualification['selected']['rows'] if row['identity_sha256'] in identities]
        else:
            held = json.loads((root/'artifacts/validation/track1-population-held-out-selection-2026-09-17.json').read_bytes())
            expected = [(row['identity_sha256'],Fraction(row['frequency']),row['allele_number'])
                        for row in held['observations']]
    observations = (() if not expected or attempt.recovery_binding is None else
        recovery_resolver.frequency_observations(attempt, tuple(key for key, _, _ in expected)))
    observed = {key:(Fraction(frequency),an) for key,frequency,an in observations}
    checks.append(('population-probe-mismatch',not expected or len(expected)==2 and all(
        observed.get(key)==(frequency,an) for key,frequency,an in expected)))
    if attempt.status != 'succeeded' and not attempt.model_coverage_counts:
        # A lost or cancelled execution has no scientific observations for these
        # checks. Null means unassessed; it never asserts a biological failure.
        checks=[(name,passed if name in {'execution-failed','replay-failed'} else None)
                for name,passed in checks]
    return tuple(checks)


def classify_profile_case_outcome(*, execution_status, profile_disposition,
                                  exact_target_generated, eligible_rank):
    """Evaluator outcome; a scope decision never asserts genotype recovery.

    Execution failure takes precedence over scope, and an advancing unresolved
    route is distinct from a completed scientific assessment with no recovery.
    """
    if execution_status not in {'succeeded', 'failed', 'timed_out', 'truncated', 'unsupported'}:
        raise PublicTrack1CoordinatorError('invalid profile case execution status')
    if profile_disposition not in {'selected', 'unresolved', 'outside_profile', 'not_applicable_upstream'}:
        raise PublicTrack1CoordinatorError('invalid profile case disposition')
    if (type(exact_target_generated) is not bool
            or eligible_rank is not None and (type(eligible_rank) is not int or eligible_rank < 1)
            or eligible_rank is not None and not exact_target_generated):
        raise PublicTrack1CoordinatorError('invalid profile case recovery observation')
    if execution_status in {'failed', 'timed_out', 'truncated'}:
        return 'failed'
    if execution_status == 'unsupported':
        return 'unsupported'
    if profile_disposition == 'outside_profile':
        if exact_target_generated or eligible_rank is not None:
            raise PublicTrack1CoordinatorError('out-of-profile target cannot establish recovery')
        return 'scoped_out'
    if profile_disposition == 'not_applicable_upstream':
        return 'unsupported'
    if exact_target_generated and eligible_rank is not None and eligible_rank <= 10:
        return 'recovered'
    if profile_disposition == 'unresolved':
        return 'unresolved'
    return 'evaluated_no_recovery'


def run_public_candidate_validation(
    *, repository_root, plan, roster, resources, policies, progress=None,
    workers=1, checkpoint_root=None, continuation=None, analysis_profile='focused',
) -> PublicCandidateValidation | PublicValidationStoppedOutcome:
    """Run one declared A/B/C campaign; stop before admitting dependent work."""
    if type(workers) is not int or not 1 <= workers <= 16:
        raise PublicTrack1CoordinatorError("invalid public validation worker count")
    root = Path(repository_root).resolve(strict=True)
    profile, profile_bytes = load_submission_acceptance_profile(root)
    require_profile_acceptance_binding(profile, plan, analysis_profile=analysis_profile)
    schedule = profile['stage_schedule']
    names = tuple(name for row in schedule for name in row['case_ids'])
    if len(policies.runtimes) != 1 or set(names) != {case.case_id for case in roster.cases}:
        raise PublicTrack1CoordinatorError('campaign declaration differs')
    stack_sha256 = _sha256(executable_stack_bytes(root))
    journal = Path(checkpoint_root) if checkpoint_root is not None else resources.scratch_root/'outcomes'/uuid4().hex
    if not journal.absolute().is_relative_to(resources.scratch_root.resolve(strict=True)):
        raise PublicTrack1CoordinatorError('invalid validation checkpoint root')
    if continuation is not None and (journal.exists() or journal.is_symlink()):
        raise PublicTrack1CoordinatorError('continuation requires a fresh checkpoint root')
    publish_validation_artifact(journal/'execution-context.json',
        _execution_context_bytes(root,plan,roster,resources,policies))
    if (journal/'validation-stopped-outcome.json').exists():
        raise PublicTrack1CoordinatorError('a stopped campaign cannot be resumed as a new attempt')
    recovery_resolver = _recovery_resolver(journal, roster, resources, policies)
    cases = {case.case_id:case for case in roster.cases}
    all_attempts, all_replays, checks = [], [], {}
    carried_forward = None
    if continuation is not None:
        carried_forward, inherited_attempts, inherited_replays = _import_passed_stage_a(
            root, journal, continuation, roster, resources, policies, profile)
        all_attempts.extend(inherited_attempts)
        all_replays.extend(inherited_replays)
        for attempt in inherited_attempts:
            case = cases[attempt.case_id]
            checks[case.case_id] = _case_gate_checks(root, roster, case, (attempt,),
                tuple(row for row in inherited_replays if row[1] == case.case_id), profile, journal=journal, recovery_resolver=recovery_resolver)
    cancellation_event=Event()
    candidate, experiment, derivation, current_stage = None, None, {}, 'A'
    executor = ThreadPoolExecutor(max_workers=workers) if workers>1 else None
    def run_stage(row):
        nonlocal current_stage
        current_stage=row['stage']
        jobs=[]
        for name in row['case_ids']:
            if carried_forward is not None and name in carried_forward['case_ids']:
                continue
            case=cases[name]
            if row['split']=='development':
                view=development_ranker_view(roster,name,runner_resources=resources.common_runner_resources)
                selected=None
            else:
                view=authenticated_held_out_ranker_view(candidate,experiment,roster,name,
                    runner_resources=resources.common_runner_resources)
                selected=candidate.selected_policy_id
                publish_validation_artifact(journal/'held-out-views'/f'{name}.json',_bytes({
                    'case_id':name,'candidate_sha256':candidate.sha256,'view':asdict(view)}))
            jobs.append((case,view,selected))
        first_failed=[]
        def received(index,result):
            attempts,replays=result
            _persist_case_attempts(journal,attempts,replays,None if candidate is None else candidate.sha256)
            all_attempts.extend(a for a in attempts if a not in all_attempts)
            all_replays.extend(r for r in replays if r not in all_replays)
            case=jobs[index][0]
            checks[case.case_id]=_case_gate_checks(root,roster,case,attempts,replays,profile,journal=journal, recovery_resolver=recovery_resolver)
            failed=next((code for code,passed in checks[case.case_id] if not passed),None)
            if failed is not None and not first_failed:first_failed.append((case.case_id,failed))
        def notify(index,result):
            if progress is not None:
                for attempt in result[0]:progress(attempt)
        def worker(payload):
            case,view,selected=payload
            return run_isolated_public_case(case,view,resources,policies,journal,policy_id=selected,
                                            cancellation_event=cancellation_event)
        def failed(index):
            case,view,_selected=jobs[index]
            return ((_failed_attempt(case,view,policies.runtimes[0].candidate,'public-worker'),),())
        if executor is None:
            results=[]
            for index,job in enumerate(jobs):
                try:result=worker(job)
                except Exception:result=failed(index)
                results.append(result)
                received(index,result)
                notify(index,result)
                if first_failed or cancellation_event.is_set():break
        else:
            _parallel_results(executor,worker,jobs,completed=notify,received=received,failed=failed,
                concurrency=workers,gate=lambda index,result:not first_failed and not cancellation_event.is_set(),
                interrupted=cancellation_event.set)
        return first_failed[0] if first_failed else None
    def stopped(failure,kind='gate-failure'):
        return _stopped_validation_outcome(root,journal,profile,current_stage,
            tuple(all_attempts),tuple(all_replays),checks,failure,candidate=candidate,stop_kind=kind)
    try:
        failure=run_stage(schedule[0])
        if cancellation_event.is_set():return stopped(None,'interrupted')
        if failure is not None:return stopped(failure)
        failure=run_stage(schedule[1])
        if cancellation_event.is_set():return stopped(None,'interrupted')
        development=tuple(a for a in all_attempts if a.split=='development')
        development_replays=tuple(r for r in all_replays if r[0]=='development')
        if len(development)!=6 or len(development_replays)!=6:
            return stopped(failure)
        experiment=complete_development_experiment(plan=plan,roster=roster,protocol=policies.protocol,
            candidates=tuple(runtime.candidate for runtime in policies.runtimes),attempts=development, recovery_resolver=recovery_resolver)
        publish_validation_artifact(journal/'development-experiment.json',experiment.as_bytes())
        if failure is not None or experiment.selected_policy_id is None:
            outcome=PublicCandidateValidation(stack_sha256,experiment,{},None,(),(),(),tuple(all_replays),
                                               'no-development-candidate', carried_forward=carried_forward)
            persist_validation_outcome(journal,outcome)
            return outcome
        derivation,derivation_bytes,criteria=_development_derivation(experiment,root)
        candidate=_candidate(experiment,repository_root=root,resources=resources,policies=policies,
            derivation_bytes=derivation_bytes,criteria=criteria,stack_sha256=stack_sha256, recovery_resolver=recovery_resolver)
        publish_validation_artifact(journal/'development-derivation.json',derivation_bytes)
        publish_validation_artifact(journal/'candidate.json',candidate.as_bytes())
        failure=run_stage(schedule[2])
        if cancellation_event.is_set():return stopped(None,'interrupted')
        held=tuple(a for a in all_attempts if a.split=='held_out')
        held_replays=tuple(r for r in all_replays if r[0]=='held_out')
        if len(held)!=2 or len(held_replays)!=2:return stopped(failure)
        evaluations=evaluate_attempts(roster,held,split='held_out',policy_id=candidate.selected_policy_id, recovery_resolver=recovery_resolver)
        metrics=tuple(sorted((name,_fraction(value)) for name,value in evaluation_metrics(evaluations).items()))
        outcome=PublicCandidateValidation(stack_sha256,experiment,derivation,candidate,held,evaluations,metrics,
                                           tuple(all_replays),'held-out-complete', carried_forward=carried_forward)
        if _sha256(executable_stack_bytes(root)) != stack_sha256:
            raise PublicTrack1CoordinatorError('validation executable changed')
        persist_validation_outcome(journal,outcome)
        return outcome
    except KeyboardInterrupt:
        # Running workers are bounded and own their terminal cleanup. Preserve
        # their durable checkpoint evidence before the interruption receipt.
        cancellation_event.set()
        if executor is not None:executor.shutdown(wait=True,cancel_futures=True)
        return stopped(None,'interrupted')
    except PublicTrack1CoordinatorError:
        raise
    except Exception:
        raise PublicTrack1CoordinatorError('public candidate validation failed') from None
    finally:
        if executor is not None:executor.shutdown(wait=True)


def _read_json(path: Path) -> tuple[dict, bytes]:
    try:
        if path.is_symlink() or not path.is_file():
            raise ValueError
        content = path.read_bytes()
        value = json.loads(content)
        if type(value) is not dict:
            raise ValueError
        return value, content
    except (OSError, ValueError, TypeError, json.JSONDecodeError):
        raise PublicTrack1CoordinatorError(
            "required public acceptance evidence is invalid"
        ) from None


def _current_stack_evidence_matches(
    record, *, schema, executable_stack_sha256, runtime_image_digest,
):
    """Accept public evidence only for its actual tested code and image."""
    return (
        isinstance(record, dict)
        and record.get("schema") == schema
        and record.get("classification") == "public"
        and record.get("status") == "passed"
        and record.get("executable_stack_sha256") == executable_stack_sha256
        and record.get("runtime_image_digest") == runtime_image_digest
    )


def _current_contract_checks(engineering, security, *, executable_stack_sha256, runtime_image_digest):
    engineering_ok = _current_stack_evidence_matches(engineering,
        schema="genosapiens.track1-current-contract-engineering-validation",
        executable_stack_sha256=executable_stack_sha256, runtime_image_digest=runtime_image_digest)
    security_ok = _current_stack_evidence_matches(security,
        schema="genosapiens.track1-security-validation",
        executable_stack_sha256=executable_stack_sha256, runtime_image_digest=runtime_image_digest)
    groups = engineering.get("groups", {})
    checks = {name: engineering_ok and groups.get(name) == "passed" for name in (
        "authorized-evidence-and-publication-contracts", "current-contract-engineering",
        "evaluator-and-numerical-policy", "public-resource-license-and-usability",
        "security-and-data-boundaries")}
    checks["security-and-data-boundaries"] = (checks["security-and-data-boundaries"]
        and security_ok
        and security.get("boundary_checks", {}).get("controlled_and_agent_context_exclusion") == "passed"
        and security.get("repository", {}).get("status") == "passed")
    return checks


def require_current_public_acceptance_evidence(repository_root, *, executable_stack_sha256):
    """Reject missing/stale finalization evidence before starting scientific work."""
    root = Path(repository_root)
    engineering, _ = _read_json(root / CURRENT_ENGINEERING_RESULT_PATH)
    security, _ = _read_json(root / CURRENT_SECURITY_RESULT_PATH)
    if not all(_current_contract_checks(engineering, security,
            executable_stack_sha256=executable_stack_sha256,
            runtime_image_digest=RUNTIME_IMAGE.split("@", 1)[1]).values()):
        raise PublicTrack1CoordinatorError("current public acceptance evidence is not qualified")


def _acceptance_evidence(
    outcome: PublicCandidateValidation,
    repository_root: Path,
    checkpoint_root: Path,
) -> tuple[AcceptanceResultReference, ...]:
    """Adjudicate the fixed, matching public evidence set by closed fields."""
    if outcome.candidate is None or outcome.status != "held-out-complete":
        raise PublicTrack1CoordinatorError("candidate validation is incomplete")
    candidate = outcome.candidate
    engineering, engineering_bytes = _read_json(
        repository_root
        / CURRENT_ENGINEERING_RESULT_PATH
    )
    empirical, empirical_bytes = _read_json(
        repository_root
        / "artifacts/validation/track1-public-empirical-current-2026-09-12.json"
    )
    from genosapiens.track1_empirical_validation import load_empirical_component_evidence
    empirical_components = load_empirical_component_evidence(repository_root=repository_root)
    profile, _profile_bytes = load_submission_acceptance_profile(repository_root)
    security, security_bytes = _read_json(
        repository_root / CURRENT_SECURITY_RESULT_PATH
    )
    all_attempts = (
        *outcome.development_experiment.attempts,
        *outcome.held_out_attempts,
    )
    replay_bindings = {
        (split, case_id, policy_id): digest
        for split, case_id, policy_id, digest in outcome.determinism_replays
    }
    determinism_ok = (
        len(replay_bindings) == len(all_attempts)
        and all(
            replay_bindings.get((item.split, item.case_id, item.policy_id))
            == item.execution_sha256
            for item in all_attempts
        )
    )
    statuses = {
        **_current_contract_checks(engineering, security,
            executable_stack_sha256=outcome.executable_stack_sha256,
            runtime_image_digest=candidate.runtime_image_digest),
        "bounded-real-tool-vcf-first": (
            determinism_ok
            and model_coverage_complete(outcome.held_out_evaluations)
            and _positive_exports_complete(profile, outcome.held_out_evaluations, "held_out")
            and bool(outcome.development_experiment.attempts)
            and all(
                item.status == "succeeded"
                for item in outcome.development_experiment.attempts
            )
            and bool(outcome.held_out_attempts)
            and all(item.status == "succeeded" for item in outcome.held_out_attempts)
            and all(item.execution_sha256 is not None for item in outcome.held_out_attempts)
        ),
        "giab-genomic-representation": (
            "giab-genomic-representation" in empirical_components.applicable_components
            and "giab-annotation" in empirical_components.applicable_components
            and
            empirical.get("schema") == "genosapiens.track1-public-empirical-validation"
            and empirical.get("giab_representation", {}).get("outcomes", {}).get(
                "status"
            ) == "accepted"
            and empirical.get("giab_annotation", {}).get("outcomes", {}).get(
                "status"
            ) == "accepted"
        ),
        "public-caller-preparation": (
            "public-caller-preparation" in empirical_components.applicable_components
            and empirical.get("caller_preparation", {}).get("status")
            == "accepted-for-bounded-early-workload"
        ),
    }
    artifacts = {
        "authorized-evidence-and-publication-contracts": engineering_bytes,
        "bounded-real-tool-vcf-first": outcome.as_bytes(),
        "current-contract-engineering": engineering_bytes,
        "evaluator-and-numerical-policy": engineering_bytes,
        "giab-genomic-representation": empirical_bytes,
        "public-caller-preparation": empirical_bytes,
        "public-resource-license-and-usability": engineering_bytes,
        "security-and-data-boundaries": security_bytes,
    }
    references = []
    for result_id in REQUIRED_ACCEPTANCE_RESULTS:
        artifact_sha256 = _sha256(artifacts[result_id])
        applicability = (
            {
                "source_report_sha256": empirical_components.source_report_sha256,
                "generating_roster_cases": empirical_components.generating_roster_cases,
                "current_roster_annotation_accepted": False,
                "applicable_components": (
                    ("giab-genomic-representation", "giab-annotation")
                    if result_id == "giab-genomic-representation"
                    else ("public-caller-preparation",)
                ),
            }
            if result_id in {"giab-genomic-representation", "public-caller-preparation"}
            else ({"tested_current_candidate": outcome.carried_forward is None,
                   "carried_forward": outcome.carried_forward}
                  if result_id == "bounded-real-tool-vcf-first" else {"tested_current_candidate": True})
        )
        dependency = _sha256(_bytes({
            "candidate_sha256": candidate.sha256,
            "executable_stack_sha256": outcome.executable_stack_sha256,
            "result_id": result_id,
            "source_artifact_sha256": artifact_sha256,
            "source_applicability": applicability,
        }))
        scope = _sha256(_bytes({
            "plan_sha256": candidate.plan_sha256,
            "roster_sha256": candidate.roster_sha256,
            "result_id": result_id,
            "validation_only": True,
            "source_applicability": applicability,
        }))
        references.append(AcceptanceResultReference(
            result_id,
            "accepted" if statuses[result_id] else "rejected",
            candidate.sha256,
            dependency,
            scope,
            artifact_sha256,
        ))
    return tuple(references)


def finalize_public_validation(
    outcome: PublicCandidateValidation,
    *,
    repository_root: str | Path,
    roster: Track1Roster,
    resources: PublicResourceContext,
    checkpoint_root: str | Path | None = None,
) -> FinalizedPublicValidation:
    """Consolidate fixed evidence and freeze only a genuinely accepted release."""
    journal = Path(checkpoint_root) if checkpoint_root is not None else (
        resources.scratch_root / "outcomes" / outcome.sha256
    )
    stage = "outcome-publication"
    try:
        root = Path(repository_root).resolve(strict=True)
        if not journal.absolute().is_relative_to(resources.scratch_root.resolve(strict=True)):
            raise ValueError
        persist_validation_outcome(journal, outcome)
        if outcome.candidate is None:
            raise ValueError
        stage = "consolidation"
        recovery_resolver = _recovery_resolver(journal, roster, resources,
            load_public_policy_family(root, model_loci=resources.model_loci))
        references = _acceptance_evidence(outcome, root, journal)
        acceptance = consolidate_public_acceptance(
            outcome.candidate,
            roster,
            outcome.held_out_attempts,
            references,
            limitations=(
                "Synthetic cases share one public HG002 background and are not independent patients.",
                "The roster is engineering validation and does not establish clinical calibration.",
                "Per-attempt memory is the cumulative process peak through that attempt, not an isolated cgroup peak.",
                "Any scientific policy change requires a new development experiment and held-out split.",
                "Public acceptance does not authorize raw-read rescue or model access to controlled inputs.",
            ),
         recovery_resolver=recovery_resolver)
        stage = "acceptance-publication"
        publish_validation_artifact(journal / "acceptance.json", acceptance.as_bytes())
        if acceptance.status != "accepted":
            result = FinalizedPublicValidation(outcome, acceptance, None, "rejected")
            stage = "final-publication"
            publish_validation_artifact(journal / "finalized.json", result.as_bytes())
            return result
        stage = "release-freeze"
        identity_path = root / (
            "config/reference-identities/challenge-grch38-masked-v2-hs38d1.json"
        )
        handoff_schema = root / "schemas/track1/mechanism-ready-handoff.schema.json"
        handoff_contract = root / "config/track1/handoff-contract.json"
        handoff = json.loads(handoff_contract.read_bytes())
        reference_bindings = _bytes({
            "identity_contract_sha256": _sha256(identity_path.read_bytes()),
            "construction_report": asdict(resources.construction_report),
            "public_identity_id": resources.reference_contract.reference_role(
                "public_annotation_benchmark"
            ).identity_id,
        })
        assessment_scope = _bytes({
            "plan_sha256": outcome.candidate.plan_sha256,
            "roster_sha256": outcome.candidate.roster_sha256,
            "required_results": REQUIRED_ACCEPTANCE_RESULTS,
            "mode": "base-supplied-call-vcf-first",
            "read_rescue": "not-authorized-not-applicable",
        })
        evidence_interface = evidence_interface_bytes(root)
        for name, content in (
            ("reference-bindings", reference_bindings),
            ("assessment-scope", assessment_scope),
            ("evidence-interface", evidence_interface),
            ("handoff-path", _bytes(handoff["canonical_artifact_path"])),
            ("executable-stack", executable_stack_bytes(root)),
        ):
            publish_validation_artifact(journal / "release-bindings" / f"{name}.json", content)
        release = freeze_base_release(
            outcome.candidate,
            acceptance,
            reference_identity_id=resources.reference_contract.reference_role(
                "public_annotation_benchmark"
            ).identity_id,
            reference_manifest_sha256=REFERENCE_MANIFEST_SHA256,
            reference_bindings_sha256=_sha256(reference_bindings),
            assessment_scope_sha256=_sha256(assessment_scope),
            evidence_interface_sha256=_sha256(evidence_interface),
            canonical_handoff_schema_sha256=_sha256(handoff_schema.read_bytes()),
            canonical_handoff_path_sha256=_sha256(_bytes(
                handoff["canonical_artifact_path"]
            )),
            path_selection_policy_sha256=_sha256(handoff_contract.read_bytes()),
        )
        stage = "release-publication"
        publish_validation_artifact(journal / "base-release.json", release.as_bytes())
        from genosapiens.track1_release import load_base_release
        load_base_release((journal / "base-release.json").read_bytes(), expected_sha256=release.sha256)
        stage = "manifest-publication"
        publish_public_release_manifest(root, resources, journal, outcome, acceptance, release)
        result = FinalizedPublicValidation(
            outcome,
            acceptance,
            release,
            "accepted-and-frozen",
        )
        stage = "final-publication"
        publish_validation_artifact(journal / "finalized.json", result.as_bytes())
        for name, content in (
            ("validation-outcome", outcome.as_bytes()),
            ("candidate", outcome.candidate.as_bytes()),
            ("acceptance", acceptance.as_bytes()),
            ("base-release", release.as_bytes()),
            ("finalized", result.as_bytes()),
        ):
            if (journal / f"{name}.json").read_bytes() != content:
                raise ValueError
        return result
    except Exception:
        record_finalization_failure(journal, stage, outcome.sha256)
        raise PublicTrack1CoordinatorError(
            f"public validation finalization failed: {stage}"
        ) from None


def _acceptance_artifact_bytes(outcome, root, journal) -> dict[str, bytes]:
    contents = [outcome.as_bytes()]
    for relative in (
        CURRENT_ENGINEERING_RESULT_PATH,
        "artifacts/validation/track1-public-empirical-current-2026-09-12.json",
        CURRENT_SECURITY_RESULT_PATH,
    ):
        _value, content = _read_json(root / relative)
        contents.append(content)
    return {_sha256(content): content for content in contents}


def publish_public_release_manifest(root, resources, journal, outcome, acceptance, release) -> Path:
    """Publish the complete accepted freeze in the public manifest namespace.

    The content-addressed directory has no mutable latest pointer. Verification
    is committed only after all sources/components have been read back exactly.
    Scientific objects retain their classifications in the execution journal.
    """
    if type(outcome) is not PublicCandidateValidation or outcome.candidate is None:
        raise PublicTrack1CoordinatorError('a complete candidate outcome is required for publication')
    persist_validation_outcome(journal,outcome)
    from genosapiens.track1_release import load_base_release
    destination = resources.reference_root / "manifests/track1/base-releases" / release.sha256
    artifacts = _acceptance_artifact_bytes(outcome, root, journal)
    if set(release.accepted_result_sha256s) != set(artifacts):
        raise ValueError
    components = {
        "base-release.json": release.as_bytes(),
        "candidate.json": outcome.candidate.as_bytes(),
        "acceptance.json": acceptance.as_bytes(),
        "validation-outcome.json": outcome.as_bytes(),
        "development-experiment.json": outcome.development_experiment.as_bytes(),
        "development-derivation.json": _bytes(outcome.development_derivation),
        **{f"accepted-results/{digest}.json": content for digest, content in artifacts.items()},
    }
    component_files = {}
    if outcome.carried_forward is not None:
        _verify_carried_forward(journal, outcome.carried_forward)
        for name in ('carried-forward.json', 'carried-forward-stack.json',
                     'carried-forward-context.json', 'carried-forward-stopped.json'):
            components[name] = _read_checkpoint(journal, name)
        for binding in outcome.carried_forward['artifacts']:
            components['carried-forward/' + binding['path']] = _read_checkpoint(journal, binding['path'])
    for path in (journal / "release-bindings").iterdir():
        components["release-bindings/" + path.name] = _read_checkpoint(journal, "release-bindings/" + path.name)
    for path in (
        root / "schemas/track1/mechanism-ready-handoff.schema.json",
        root / "config/track1/handoff-contract.json",
        root / "config/reference-identities/challenge-grch38-masked-v2-hs38d1.json",
    ):
        if path.is_symlink() or not path.is_file():
            raise ValueError
        components["contracts/" + path.name] = path.read_bytes()
    for content in resources.common_resource_bytes:
        components[f"resource-descriptors/{_sha256(content)}.json"] = content
    evaluator_pin = (root / "config/track1/challenge-evaluator.json").read_bytes()
    evaluator_source = resources.evaluator_source_path.read_bytes()
    policies = load_public_policy_family(root, model_loci=resources.model_loci)
    selected = next(item for item in policies.runtimes if item.candidate.policy_id == outcome.candidate.selected_policy_id)
    if selected.candidate.bindings != release.policy_bindings:
        raise ValueError
    components["policies/selected-fusion.json"] = selected.fusion_policy_bytes
    components["policies/evaluator-binding.json"] = _bytes({
        "pin_sha256": _sha256(evaluator_pin), "source_sha256": _sha256(evaluator_source),
    })
    components["policies/evaluator-source.py"] = evaluator_source
    components["policies/numerical-binding.json"] = _bytes({
        "epcr_policy_sha256": policies.epcr_policy_sha256,
        "epcr_suitability_sha256": _sha256(policies.epcr_derivation_bytes),
    })
    stack = _checkpoint_json(components["release-bindings/executable-stack.json"])
    for entry in stack["files"]:
        path = root / entry["path"]
        if not path.is_relative_to(root) or path.is_symlink():
            raise ValueError
        content = path.read_bytes()
        if len(content) != entry["size_bytes"] or _sha256(content) != entry["sha256"]:
            raise ValueError
        components["stack-files/" + entry["path"]] = content
    bound_hashes = {_sha256(content) for content in components.values()}
    required_hashes = {
        release.candidate_sha256, release.consolidated_acceptance_sha256,
        release.reference_manifest_sha256, release.reference_bindings_sha256,
        release.assessment_scope_sha256, release.evidence_interface_sha256,
        release.canonical_handoff_schema_sha256, release.canonical_handoff_path_sha256,
        release.path_selection_policy_sha256, release.executable_stack_sha256,
        release.evaluator_pin_sha256, release.validation_plan_sha256,
        release.validation_roster_sha256, *release.resource_manifest_sha256s,
        outcome.candidate.numerical_policy_evidence_sha256,
    }
    required_hashes.update(release.policy_bindings.component_sha256s)
    # Plan/roster remain evaluator-side; publishing them does not mount them in
    # scientific workers. Their source bytes complete the accepted release DAG.
    for relative in (
        "tests/fixtures/track1/base-release-validation-plan.synthetic.json",
        "tests/fixtures/track1/roster/roster.json",
        "config/track1/submission-acceptance-profile.json",
    ):
        components["validation-inputs/" + Path(relative).name] = (root / relative).read_bytes()
    profile, _profile_bytes = load_submission_acceptance_profile(root)
    source_binding = profile["source_roster"]
    source_content = (root / source_binding["path"]).read_bytes()
    if (len(source_content) != source_binding["size_bytes"]
            or _sha256(source_content) != source_binding["sha256"]):
        raise ValueError
    # Fixture JSON is authenticated by its exact source binding and loader;
    # its whitespace need not follow the checkpoint serialization contract.
    source_roster = json.loads(source_content)
    source_cases = {row["case_id"]: row for row in source_roster["cases"]}
    parsed_cases = {case.case_id: case for case in load_submission_acceptance_roster(repository_root=root).cases}
    format_bindings = tuple(
        {"path": relative, "sha256": digest}
        for row in profile["cases"]
        for relative, digest in zip(source_cases[row["case_id"]]["source_format_bindings"],
                                   parsed_cases[row["case_id"]].runner_view.source_format_sha256s,
                                   strict=True)
    )
    fixture_bindings = (
        source_binding, source_roster["annotation_expectations"], source_roster["causal_truth"],
        *profile["configuration_bindings"],
        *(entry for row in profile["cases"] for entry in row["runner_inputs"]),
        *format_bindings,
    )
    for entry in fixture_bindings:
        relative = Path(entry["path"])
        path = root / relative
        if relative.is_absolute() or ".." in relative.parts or path.is_symlink():
            raise ValueError
        size = path.stat().st_size
        if "size_bytes" in entry and size != entry["size_bytes"]:
            raise ValueError
        component_files["validation-inputs/source/" + relative.as_posix()] = (path, entry["sha256"], size)
    from genosapiens.track1_validation_plan import load_validation_plan_from_pin
    components["validation-inputs/plan.canonical.json"] = validation_plan_as_json(
        load_validation_plan_from_pin(repository_root=root)
    )
    bound_hashes.update(_sha256(content) for content in components.values())
    bound_hashes.update(digest for _path, digest, _size in component_files.values())
    if not required_hashes.issubset(bound_hashes):
        raise ValueError
    for relative, content in components.items():
        publish_validation_artifact(destination / relative, content)
    for relative, (path, digest, size) in component_files.items():
        publish_validation_file(destination / relative, path, expected_sha256=digest, expected_size_bytes=size)
    for relative, content in components.items():
        if _read_checkpoint(destination, relative) != content:
            raise ValueError
    for relative, (_path, digest, size) in component_files.items():
        descriptor = os.open(destination / relative, os.O_RDONLY | os.O_NOFOLLOW | os.O_NONBLOCK)
        with os.fdopen(descriptor, "rb") as stream:
            before = os.fstat(stream.fileno())
            if not stat.S_ISREG(before.st_mode) or before.st_size != size:
                raise ValueError
            actual = hashlib.file_digest(stream, "sha256").hexdigest()
            after = os.fstat(stream.fileno())
            if (actual != digest or
                    (before.st_dev, before.st_ino, before.st_size, before.st_mtime_ns) !=
                    (after.st_dev, after.st_ino, after.st_size, after.st_mtime_ns)):
                raise ValueError
    load_base_release(_read_checkpoint(destination, "base-release.json"), expected_sha256=release.sha256)
    verification = _bytes({
        "schema": "genosapiens.track1-public-release-publication-verification",
        "classification": "public", "status": "accepted-and-published",
        "base_release_sha256": release.sha256,
        "candidate_sha256": outcome.candidate.sha256,
        "acceptance_sha256": acceptance.sha256,
        "manifest_relative_path": destination.relative_to(resources.reference_root).as_posix() + "/base-release.json",
        "components": {**{relative: _sha256(content) for relative, content in components.items()},
                       **{relative: digest for relative, (_path, digest, _size) in component_files.items()}},
        "controlled_read_rescue_authorized": False,
    })
    publish_validation_artifact(destination / "verification.json", verification)
    if _read_checkpoint(destination, "verification.json") != verification:
        raise ValueError
    return destination / "base-release.json"


def publish_validation_artifact(path: str | Path, content: bytes) -> str:
    """Durably publish complete bytes without replacing an existing artifact.

    Private files also retain controlled/no-LLM scientific envelopes; only the
    evaluator ledger and static operational failure records are public metadata.
    The hard link is the atomic commit, and both file and directory are fsynced.
    """
    temporary = None
    try:
        target = Path(path).absolute()
        if type(content) is not bytes:
            raise ValueError
        for parent in (*reversed(target.parents), target):
            if parent.is_symlink():
                raise ValueError
        target.parent.mkdir(parents=True, exist_ok=True, mode=0o700)
        with tempfile.NamedTemporaryFile(dir=target.parent, delete=False) as handle:
            temporary = Path(handle.name)
            os.chmod(temporary, 0o600)
            handle.write(content)
            handle.flush()
            os.fsync(handle.fileno())
        try:
            os.link(temporary, target)
        except FileExistsError:
            if target.is_symlink() or not target.is_file() or target.read_bytes() != content:
                raise ValueError
        descriptor = os.open(target.parent, os.O_RDONLY | os.O_DIRECTORY)
        try:
            os.fsync(descriptor)
        finally:
            os.close(descriptor)
        return _sha256(content)
    except Exception:
        raise PublicTrack1CoordinatorError("validation artifact publication failed") from None
    finally:
        if temporary is not None:
            temporary.unlink(missing_ok=True)


def publish_validation_file(path: str | Path, source_path: str | Path, *, expected_sha256: str,
                            expected_size_bytes: int) -> str:
    """Copy an authenticated public input without buffering its whole callset."""
    import stat
    temporary = None
    try:
        target, source = Path(path).absolute(), Path(source_path).absolute()
        if (type(expected_size_bytes) is not int or expected_size_bytes <= 0
                or type(expected_sha256) is not str or re.fullmatch(r"[0-9a-f]{64}", expected_sha256) is None
                or any(parent.is_symlink() for value in (target, source) for parent in (value, *value.parents))):
            raise ValueError
        target.parent.mkdir(parents=True, exist_ok=True, mode=0o700)
        descriptor = os.open(source, os.O_RDONLY | os.O_NOFOLLOW | os.O_NONBLOCK)
        with os.fdopen(descriptor, "rb") as stream, tempfile.NamedTemporaryFile(dir=target.parent, delete=False) as out:
            temporary = Path(out.name)
            os.chmod(temporary, 0o600)
            before = os.fstat(stream.fileno())
            if not stat.S_ISREG(before.st_mode) or before.st_size != expected_size_bytes:
                raise ValueError
            size, digest = 0, hashlib.sha256()
            for block in iter(lambda: stream.read(65536), b""):
                size += len(block)
                if size > expected_size_bytes:
                    raise ValueError
                out.write(block)
                digest.update(block)
            after = os.fstat(stream.fileno())
            identity = lambda row: (row.st_dev, row.st_ino, row.st_size, row.st_mtime_ns, row.st_ctime_ns)
            if (size != expected_size_bytes or digest.hexdigest() != expected_sha256
                    or identity(before) != identity(after)):
                raise ValueError
            out.flush()
            os.fsync(out.fileno())
        try:
            os.link(temporary, target)
        except FileExistsError:
            if target.is_symlink() or not target.is_file() or target.stat().st_size != expected_size_bytes:
                raise ValueError
            digest = hashlib.sha256()
            with target.open("rb") as stream:
                for block in iter(lambda: stream.read(65536), b""):
                    digest.update(block)
            if digest.hexdigest() != expected_sha256:
                raise ValueError
        descriptor = os.open(target.parent, os.O_RDONLY | os.O_DIRECTORY)
        try:
            os.fsync(descriptor)
        finally:
            os.close(descriptor)
        return expected_sha256
    except Exception:
        raise PublicTrack1CoordinatorError("validation file publication failed") from None
    finally:
        if temporary is not None:
            temporary.unlink(missing_ok=True)


def _execution_context_bytes(root, plan, roster, resources, policies) -> bytes:
    if (plan.submission_acceptance_profile_sha256 != roster.roster_sha256
            or roster != load_submission_acceptance_roster(repository_root=root)):
        raise PublicTrack1CoordinatorError("execution roster differs from the declared acceptance profile")
    loci = _authenticated_policy_loci(root, resources.model_loci)
    if policies.hypothesis_policy.loci != loci:
        raise PublicTrack1CoordinatorError("execution model-locus binding differs")
    return _bytes({
        "schema": "genosapiens.track1-public-validation-execution-context",
        "classification": "public",
        "validation_only": True,
        "executable_stack_sha256": _sha256(executable_stack_bytes(root)),
        "plan_sha256": _sha256(validation_plan_as_json(plan)),
        "roster_sha256": roster.roster_sha256,
        "model_locus_resource": {
            "descriptor_sha256": resources.model_loci.descriptor_sha256,
            "manifest_sha256": resources.model_loci.manifest_sha256,
            "output_sha256": resources.model_loci.output_sha256,
        },
        "policy_family_sha256": _sha256(_bytes({
            "protocol": asdict(policies.protocol),
            "candidates": tuple(asdict(item.candidate) for item in policies.runtimes),
        })),
        "resource_descriptor_sha256s": tuple(sorted(
            item.descriptor_sha256 for item in resources.common_runner_resources
        )),
        "input_sha256s": tuple(sorted({
            item.sha256 for case in roster.cases for item in case.runner_view.inputs
        })),
    })


def _persist_execution(journal: Path, execution) -> None:
    """Retain complete bounded products before the completion marker."""
    from genosapiens.wgs_public_store import persist_public_execution
    persist_public_execution(execution, objects=journal / "scientific-objects")


def _carried_context(journal, case_id):
    path = Path(journal) / 'carried-forward.json'
    if not path.exists():
        return None
    record = _checkpoint_json(_read_checkpoint(journal, 'carried-forward.json'))
    return record if case_id in record['case_ids'] else None


def _attempt_context_sha256(journal, case_id):
    carried = _carried_context(journal, case_id)
    return (carried['execution_context_sha256'] if carried is not None
            else _sha256(_read_checkpoint(journal, 'execution-context.json')))


def _verify_carried_forward(journal, record):
    if record is None:
        if (Path(journal) / 'carried-forward.json').exists():
            raise ValueError
        return
    if (_read_checkpoint(journal, 'carried-forward.json') != _bytes(record)
            or record['current_execution_context_sha256'] != _sha256(_read_checkpoint(journal, 'execution-context.json'))
            or record['stopped_outcome_sha256'] != _sha256(_read_checkpoint(journal, 'carried-forward-stopped.json'))
            or record['execution_context_sha256'] != _sha256(_read_checkpoint(journal, 'carried-forward-context.json'))
            or record['executable_stack_sha256'] != _sha256(_read_checkpoint(journal, 'carried-forward-stack.json'))):
        raise ValueError
    for binding in record['artifacts']:
        if _journal_artifact_binding(journal, Path(journal) / binding['path']) != binding:
            raise ValueError


def _import_passed_stage_a(root, journal, continuation, roster, resources, policies, profile):
    """Carry four authenticated passes under their original stack and context.

    An explicit stopped-outcome pin and original stack descriptor are required.
    A failed or unstarted case is never eligible, and policy/resource/input
    changes cannot use this implementation-repair continuation.
    """
    try:
        if set(continuation) != {'checkpoint_root', 'outcome_sha256', 'executable_stack'}:
            raise ValueError
        prior = Path(continuation['checkpoint_root']).resolve(strict=True)
        if (prior == journal or any(path.is_symlink() for path in (prior, *prior.parents))
                or not prior.is_relative_to(resources.reference_root / 'scratch/public-build')):
            raise ValueError
        content = _read_checkpoint(prior, 'validation-stopped-outcome.json')
        if _sha256(content) != continuation['outcome_sha256']:
            raise ValueError
        stopped = _checkpoint_json(content)
        old_context = _read_checkpoint(prior, 'execution-context.json')
        old_data = _checkpoint_json(old_context)
        current = _checkpoint_json(_read_checkpoint(journal, 'execution-context.json'))
        original_stack = _bytes(continuation['executable_stack'])
        if (stopped['execution_context_sha256'] != _sha256(old_context)
                or stopped['executable_stack_sha256'] != _sha256(original_stack)
                or old_data['executable_stack_sha256'] != _sha256(original_stack)
                or {k:v for k,v in old_data.items() if k != 'executable_stack_sha256'} !=
                   {k:v for k,v in current.items() if k != 'executable_stack_sha256'}
                or stopped['candidate'] is not None or stopped['stop_stage'] != 'B'):
            raise ValueError
        names = tuple(profile['stage_schedule'][0]['case_ids'])
        if len(names) != 4 or not set(names) <= set(dict(stopped['gate_checks'])):
            raise ValueError
        checks = dict(stopped['gate_checks'])
        if any(not all(passed is True for _, passed in checks[name]) for name in names):
            raise ValueError
        attempts, replays, selected = [], [], []
        invocations = set()
        for binding in stopped['worker_evidence_artifacts']:
            if Path(binding['path']).parent.name == 'worker-lifecycle':
                lifecycle = _checkpoint_json(_read_checkpoint(prior, binding['path']))
                if lifecycle.get('case_id') in names:
                    invocations.add(Path(binding['path']).stem)
        for key in ('attempt_artifacts', 'replay_artifacts', 'replay_started_artifacts', 'worker_evidence_artifacts'):
            for binding in stopped[key]:
                if _journal_artifact_binding(prior, prior / binding['path']) != binding:
                    raise ValueError
                data = _checkpoint_json(_read_checkpoint(prior, binding['path']))
                case_id = data.get('case_id', data.get('attempt', {}).get('case_id'))
                if key == 'replay_artifacts':
                    case_id = data['replay'][1]
                parts = Path(binding['path']).parts
                if (case_id not in names and not set(parts) & set(names)
                        and not any(part.split('.')[0] in invocations for part in parts)):
                    continue
                selected.append(binding)
                if key == 'attempt_artifacts':
                    if data['execution_context_sha256'] != _sha256(old_context):
                        raise ValueError
                    attempts.append(_load_attempt(data['attempt']))
                elif key == 'replay_artifacts':
                    if data['execution_context_sha256'] != _sha256(old_context):
                        raise ValueError
                    replays.append(tuple(data['replay']))
        if ({a.case_id for a in attempts} != set(names) or len(attempts) != 4
                or len(replays) != 4):
            raise ValueError
        recovery_resolver = _recovery_resolver(prior, roster, resources, policies, original_stack=original_stack)
        for attempt in attempts:
            case = next(case for case in roster.cases if case.case_id == attempt.case_id)
            runtime = next(row for row in policies.runtimes if row.candidate.policy_id == attempt.policy_id)
            view = development_ranker_view(roster, case.case_id, runner_resources=resources.common_runner_resources)
            if (attempt.split != 'development' or attempt.attempt_number != 1
                    or attempt.policy_bindings_sha256 != runtime.candidate.bindings.sha256
                    or attempt.input_sha256s != tuple(row.sha256 for row in view.inputs)
                    or attempt.resource_descriptor_sha256s != tuple(row.descriptor_sha256 for row in view.resources)
                    or not all(passed is True for _, passed in _case_gate_checks(root, roster, case,
                        (attempt,), tuple(row for row in replays if row[1] == case.case_id), profile, journal=prior, recovery_resolver=recovery_resolver))):
                raise ValueError
        record = {'schema': 'genosapiens.track1-public-validation-carried-forward', 'version': 1,
            'classification': 'public', 'case_ids': list(names), 'prior_checkpoint_root': str(prior),
            'stopped_outcome_sha256': continuation['outcome_sha256'],
            'execution_context_sha256': _sha256(old_context),
            'executable_stack_sha256': _sha256(original_stack),
            'current_execution_context_sha256': _sha256(_read_checkpoint(journal, 'execution-context.json')),
            'artifacts': selected, 'scientific_invocations_repeated': False}
        publish_validation_artifact(journal / 'carried-forward-stopped.json', content)
        publish_validation_artifact(journal / 'carried-forward-context.json', old_context)
        publish_validation_artifact(journal / 'carried-forward-stack.json', original_stack)
        for binding in selected:
            publish_validation_artifact(journal / binding['path'], _read_checkpoint(prior, binding['path']))
        # Content-addressed objects retain their original bytes and authority.
        # Copying them does not rerun preparation, annotation, or ranking.
        for path in (prior / 'scientific-objects').iterdir():
            if path.is_symlink() or not path.is_file() or re.fullmatch(r'[0-9a-f]{64}', path.name) is None:
                raise ValueError
            publish_validation_file(journal / 'scientific-objects' / path.name, path,
                expected_sha256=path.name, expected_size_bytes=path.stat().st_size)
        publish_validation_artifact(journal / 'carried-forward.json', _bytes(record))
        return record, tuple(attempts), tuple(replays)
    except Exception:
        raise PublicTrack1CoordinatorError('passed-case continuation evidence is invalid') from None


def _verify_saved_attempt(journal, attempt, case, view, runtime, resources, policies, *, recovery_resolver=None):
    from genosapiens.wgs_public_store import read_public_execution
    from genosapiens.track1_validation_plan import ranker_case_view_as_json
    if recovery_resolver is not None:
        manifest = recovery_resolver.manifest(attempt)
        if (manifest.case_id != case.case_id or manifest.candidate_stack_sha256 != runtime.candidate.bindings.sha256
                or manifest.ranker_case_view_sha256 != _sha256(ranker_case_view_as_json(view))):
            raise ValueError
        return
    carried = _carried_context(journal, case.case_id)
    if carried is not None:
        stack = _read_checkpoint(journal, 'carried-forward-stack.json')
        if _sha256(stack) != carried['executable_stack_sha256']:
            raise ValueError
        resources = replace(resources, executable_stack_content=stack)
    rule, _transform, _rescues = _population_template(policies.policy_bytes)
    manifest, metrics = read_public_execution(journal / "scientific-objects", expected_sha256=attempt.execution_sha256,
        resources=resources, policies=policies, scratch_root=journal / "readback-scratch",
        population_selector=(rule.subset, rule.stratum, rule.ancestry))
    if (manifest.case_id != case.case_id or manifest.candidate_stack_sha256 != runtime.candidate.bindings.sha256
            or manifest.ranker_case_view_sha256 != _sha256(ranker_case_view_as_json(view))):
        raise ValueError
    reconstructed = ValidationAttempt(case_id=case.case_id, split=case.split, scope=case.scope,
        policy_id=runtime.candidate.policy_id, policy_bindings_sha256=runtime.candidate.bindings.sha256,
        attempt_number=attempt.attempt_number, execution_sha256=manifest.sha256,
        input_sha256s=manifest.runner_input_sha256s, resource_descriptor_sha256s=manifest.runner_resource_descriptor_sha256s,
        stage_statuses=tuple((key, "succeeded") for key, _ in manifest.output_bindings),
        measurements=attempt.measurements, **metrics)
    if reconstructed != attempt:
        raise ValueError


def _persist_case_attempts(journal, attempts, replays, candidate_sha256) -> None:
    context = (journal / "execution-context.json").read_bytes()
    for attempt in attempts:
        content = _bytes({
            "schema": "genosapiens.track1-public-validation-attempt-checkpoint",
            "classification": "public",
            "validation_only": True,
            "execution_context_sha256": _sha256(context),
            "attempt": asdict(attempt),
        })
        name = f"{attempt.split}.{attempt.case_id}.{attempt.policy_id}.{attempt.attempt_number}.json"
        publish_validation_artifact(journal / "attempts" / name, content)
    for replay in replays:
        content = _bytes({
            "schema": "genosapiens.track1-public-validation-determinism-checkpoint",
            "classification": "public",
            "execution_context_sha256": _sha256(context),
            "candidate_sha256": candidate_sha256,
            "replay": replay,
        })
        publish_validation_artifact(journal / "determinism" / f"{replay[0]}.{replay[1]}.{replay[2]}.json", content)


def _persist_replay_started(journal, attempt, *, request_sha256=None, source_marker_sha256=None):
    context = (journal / "execution-context.json").read_bytes()
    content = _bytes({
        "schema": "genosapiens.track1-public-validation-replay-started",
        "classification": "public", "execution_context_sha256": _sha256(context),
        "case_id": attempt.case_id, "split": attempt.split, "policy_id": attempt.policy_id,
        "attempt_number": attempt.attempt_number, "execution_sha256": attempt.execution_sha256,
        "request_sha256": request_sha256, "source_marker_sha256": source_marker_sha256,
    })
    publish_validation_artifact(journal / "replay-started" /
        f"{attempt.split}.{attempt.case_id}.{attempt.policy_id}.json", content)


def persist_validation_outcome(journal: str | Path, outcome: PublicCandidateValidation) -> str:
    journal = Path(journal)
    if type(outcome) is not PublicCandidateValidation:
        raise PublicTrack1CoordinatorError("invalid validation outcome")
    _verify_carried_forward(journal, outcome.carried_forward)
    attempts = (*outcome.development_experiment.attempts,*outcome.held_out_attempts)
    replay_by_key = {row[:3]:row[3] for row in outcome.determinism_replays}
    if ((journal/'validation-stopped-outcome.json').exists()
            or len(outcome.development_experiment.attempts)!=6
            or len(outcome.held_out_attempts)!=(2 if outcome.candidate is not None else 0)
            or len(replay_by_key)!=len(attempts)
            or len(outcome.determinism_replays)!=len(attempts)
            or any(replay_by_key.get((a.split,a.case_id,a.policy_id))!=a.execution_sha256 for a in attempts)):
        raise PublicTrack1CoordinatorError('complete validation requires every attempt and replay')
    if outcome.candidate is not None:
        publish_validation_artifact(journal / "candidate.json", outcome.candidate.as_bytes())
    return publish_validation_artifact(journal / "validation-outcome.json", outcome.as_bytes())


def _journal_artifact_binding(journal, path):
    path = Path(path)
    relative = path.relative_to(journal).as_posix()
    content = _read_checkpoint(journal, relative)
    return {'path':relative,'size_bytes':len(content),'sha256':_sha256(content)}


def _stopped_validation_outcome(root, journal, profile, stage, attempts, replays, checks,
                                first_failed, *, candidate=None, stop_kind='gate-failure', persist=True):
    schedule = tuple((row['stage'],row['split'],tuple(row['case_ids'])) for row in profile['stage_schedule'])
    names = tuple(name for _,_,ids in schedule for name in ids)
    attempted = tuple(name for name in names if any(a.case_id == name for a in attempts))
    completed = tuple(name for name in names if any(r[1] == name for r in replays))
    marker_paths = tuple(sorted((journal/'replay-started').glob('*.json')))
    marker_names = set()
    for path in marker_paths:
        value = _checkpoint_json(_read_checkpoint(journal,path.relative_to(journal).as_posix()))
        attempt = next((a for a in attempts if a.case_id == value.get('case_id')), None)
        if (attempt is None or value != {
                'schema':'genosapiens.track1-public-validation-replay-started',
                'classification':'public',
                'execution_context_sha256':_attempt_context_sha256(journal, attempt.case_id),
                'case_id':attempt.case_id,'split':attempt.split,'policy_id':attempt.policy_id,
                'attempt_number':attempt.attempt_number,'execution_sha256':attempt.execution_sha256,
                'request_sha256':value.get('request_sha256'),
                'source_marker_sha256':value.get('source_marker_sha256')}
                or path.name != f'{attempt.split}.{attempt.case_id}.{attempt.policy_id}.json'
                or attempt.case_id in marker_names
                or (value['request_sha256'] is None) != (value['source_marker_sha256'] is None)
                or any(item is not None and re.fullmatch(r'[0-9a-f]{64}',item) is None
                       for item in (value['request_sha256'],value['source_marker_sha256']))):
            raise PublicTrack1CoordinatorError('stopped replay context is invalid')
        marker_names.add(value['case_id'])
    started = tuple(name for name in names if name in marker_names or name in completed)
    reason = {'A':'stage-a-gate-failed','B':'development-gate-failed','C':'held-out-gate-failed'}[stage]
    if stop_kind != 'gate-failure':reason=stop_kind
    attempt_paths = tuple(journal/'attempts'/f'{a.split}.{a.case_id}.{a.policy_id}.{a.attempt_number}.json'
                          for a in sorted(attempts,key=lambda a:names.index(a.case_id)))
    replay_paths = tuple(journal/'determinism'/f'{r[0]}.{r[1]}.{r[2]}.json'
                         for r in sorted(replays,key=lambda r:names.index(r[1])))
    evidence = []
    for directory in ('worker-specifications','worker-inspections','worker-lifecycle','worker-observations'):
        evidence.extend(path for path in (journal/directory).rglob('*.json') if path.is_file())
    for name in attempted:
        scratch = journal/'isolated-cases'/name/'scratch'
        evidence.extend((scratch/'requests').glob('*.json'))
        for checkpoint in (scratch/'checkpoints').glob('*'):
            evidence.extend(path for path in checkpoint.rglob('*.json') if path.is_file()
                            and 'scientific-objects' not in path.parts)
    views = tuple(sorted((journal/'held-out-views').glob('*.json'))) if candidate is not None else ()
    outcome = PublicValidationStoppedOutcome(
        _sha256(_read_checkpoint(journal,'execution-context.json')),
        _sha256(executable_stack_bytes(root)),_sha256((root/'config/track1/submission-acceptance-profile.json').read_bytes()),
        schedule,names,stage,stop_kind,
        None if first_failed is None else first_failed[0],None if first_failed is None else first_failed[1],
        attempted,started,completed,tuple((name,reason) for name in names if name not in attempted),
        tuple((name,reason) for name in names if name not in started),
        tuple(_journal_artifact_binding(journal,path) for path in attempt_paths),
        tuple(_journal_artifact_binding(journal,path) for path in marker_paths),
        tuple(_journal_artifact_binding(journal,path) for path in replay_paths),
        tuple(_journal_artifact_binding(journal,path) for path in sorted(set(evidence))),
        tuple((name,tuple(checks[name])) for name in names if name in checks),candidate,
        tuple(_journal_artifact_binding(journal,path) for path in views),
        carried_forward=(_checkpoint_json(_read_checkpoint(journal, 'carried-forward.json'))
            if (journal / 'carried-forward.json').exists() else None))
    if persist:persist_stopped_validation_outcome(journal,outcome)
    return outcome


def persist_stopped_validation_outcome(journal, outcome):
    if type(outcome) is not PublicValidationStoppedOutcome:
        raise PublicTrack1CoordinatorError('invalid stopped validation outcome')
    journal = Path(journal)
    _verify_carried_forward(journal, outcome.carried_forward)
    context_content=_read_checkpoint(journal,'execution-context.json')
    context=_checkpoint_json(context_content)
    if (_sha256(context_content) != outcome.execution_context_sha256
            or context['executable_stack_sha256'] != outcome.executable_stack_sha256
            or context['roster_sha256'] != outcome.profile_sha256):
        raise PublicTrack1CoordinatorError('stopped validation context changed')
    for collection in (outcome.attempt_artifacts,outcome.replay_started_artifacts,outcome.replay_artifacts,
                       outcome.worker_evidence_artifacts,outcome.released_held_out_views):
        for binding in collection:
            if _journal_artifact_binding(journal,journal/binding['path']) != binding:
                raise PublicTrack1CoordinatorError('stopped validation artifact changed')
    if (journal/'validation-outcome.json').exists():
        raise PublicTrack1CoordinatorError('stopped validation cannot replace a complete outcome')
    return publish_validation_artifact(journal/'validation-stopped-outcome.json',outcome.as_bytes())


def load_stopped_validation_outcome_checkpoint(checkpoint_root, *, expected_sha256,
        repository_root, plan, roster, resources, policies):
    """Authenticate a stopped ledger without granting resume or finalization authority."""
    try:
        root=Path(repository_root).resolve(strict=True)
        journal=Path(checkpoint_root).absolute()
        if not journal.is_relative_to(resources.scratch_root.resolve(strict=True)):
            raise ValueError
        context=_execution_context_bytes(root,plan,roster,resources,policies)
        if _read_checkpoint(journal,'execution-context.json')!=context:
            raise ValueError
        content=_read_checkpoint(journal,'validation-stopped-outcome.json')
        if _sha256(content)!=expected_sha256 or (journal/'validation-outcome.json').exists():
            raise ValueError
        data=_checkpoint_json(content)
        if set(data)!={field.name for field in fields(PublicValidationStoppedOutcome)}|{'accounting'}:
            raise ValueError
        _verify_carried_forward(journal, data['carried_forward'])
        recovery_resolver = _recovery_resolver(journal, roster, resources, policies)
        attempts=[]
        replays=[]
        for binding in data['attempt_artifacts']:
            if _journal_artifact_binding(journal,journal/binding['path'])!=binding:raise ValueError
            envelope=_checkpoint_json(_read_checkpoint(journal,binding['path']))
            if (envelope['schema']!='genosapiens.track1-public-validation-attempt-checkpoint'
                    or envelope['classification']!='public' or envelope['validation_only'] is not True
                    or envelope['execution_context_sha256']!=_attempt_context_sha256(journal,
                        envelope.get('attempt', {}).get('case_id') or envelope.get('replay', [None, None])[1])):raise ValueError
            attempt=_load_attempt(envelope['attempt'])
            case=next(case for case in roster.cases if case.case_id==attempt.case_id)
            runtime=next(runtime for runtime in policies.runtimes if runtime.candidate.policy_id==attempt.policy_id)
            if (attempt.split!=case.split or attempt.scope!=case.scope
                    or attempt.attempt_number!=1 or attempt.policy_bindings_sha256!=runtime.candidate.bindings.sha256
                    or attempt.input_sha256s!=tuple(row.sha256 for row in case.runner_view.inputs)):raise ValueError
            attempts.append(attempt)
        if {path.relative_to(journal).as_posix() for path in (journal/'attempts').glob('*.json')} != {
                binding['path'] for binding in data['attempt_artifacts']}:raise ValueError
        for binding in data['replay_artifacts']:
            if _journal_artifact_binding(journal,journal/binding['path'])!=binding:raise ValueError
            envelope=_checkpoint_json(_read_checkpoint(journal,binding['path']))
            if (envelope['schema']!='genosapiens.track1-public-validation-determinism-checkpoint'
                    or envelope['classification']!='public'
                    or envelope['execution_context_sha256']!=_attempt_context_sha256(journal,
                        envelope.get('attempt', {}).get('case_id') or envelope.get('replay', [None, None])[1])):raise ValueError
            row=tuple(envelope['replay'])
            if not any(row==(a.split,a.case_id,a.policy_id,a.execution_sha256) for a in attempts):raise ValueError
            replays.append(row)
        if {path.relative_to(journal).as_posix() for path in (journal/'determinism').glob('*.json')} != {
                binding['path'] for binding in data['replay_artifacts']}:raise ValueError
        candidate=None
        if data['candidate'] is not None:
            development=tuple(a for a in attempts if a.split=='development')
            experiment=complete_development_experiment(plan=plan,roster=roster,protocol=policies.protocol,
                candidates=tuple(runtime.candidate for runtime in policies.runtimes),attempts=development, recovery_resolver=recovery_resolver)
            _derivation,derivation_bytes,criteria=_development_derivation(experiment,root)
            candidate=_candidate(experiment,repository_root=root,resources=resources,policies=policies,
                derivation_bytes=derivation_bytes,criteria=criteria,stack_sha256=_sha256(executable_stack_bytes(root)), recovery_resolver=recovery_resolver)
            if candidate.as_bytes()!=_bytes(data['candidate']) or _read_checkpoint(journal,'candidate.json')!=candidate.as_bytes():
                raise ValueError
            if (_read_checkpoint(journal,'development-experiment.json') != experiment.as_bytes()
                    or _read_checkpoint(journal,'development-derivation.json') != derivation_bytes):
                raise ValueError
        for attempt in attempts:
            case=next(case for case in roster.cases if case.case_id==attempt.case_id)
            view=(development_ranker_view(roster,case.case_id,runner_resources=resources.common_runner_resources)
                  if case.split=='development' else authenticated_held_out_ranker_view(
                      candidate,experiment,roster,case.case_id,runner_resources=resources.common_runner_resources))
            if attempt.resource_descriptor_sha256s != tuple(row.descriptor_sha256 for row in view.resources):
                raise ValueError
            if attempt.execution_sha256 is not None:
                runtime = next(row for row in policies.runtimes if row.candidate.policy_id == attempt.policy_id)
                _verify_saved_attempt(journal, attempt, case, view, runtime, resources, policies, recovery_resolver=recovery_resolver)

        if candidate is not None:
            expected_views=set()
            for case in roster.cases:
                if case.split != 'held_out':continue
                view=authenticated_held_out_ranker_view(candidate,experiment,roster,case.case_id,
                    runner_resources=resources.common_runner_resources)
                relative=f'held-out-views/{case.case_id}.json'
                expected_views.add(relative)
                if _read_checkpoint(journal,relative)!=_bytes({
                        'case_id':case.case_id,'candidate_sha256':candidate.sha256,'view':asdict(view)}):
                    raise ValueError
            if {path.relative_to(journal).as_posix() for path in (journal/'held-out-views').glob('*.json')}!=expected_views:
                raise ValueError
        elif (journal/'candidate.json').exists() or (journal/'held-out-views').exists():
            raise ValueError
        profile,_=load_submission_acceptance_profile(root)
        checks={}
        for attempt in attempts:
            case=next(case for case in roster.cases if case.case_id==attempt.case_id)
            checks[case.case_id]=_case_gate_checks(root,roster,case,(attempt,),
                tuple(row for row in replays if row[1]==case.case_id),profile,journal=journal, recovery_resolver=recovery_resolver)
        first_failed=None if data['first_failed_case_id'] is None else (
            data['first_failed_case_id'],data['first_failed_check'])
        if first_failed is not None and dict(checks[first_failed[0]]).get(first_failed[1]) is not False:
            raise ValueError
        outcome=_stopped_validation_outcome(root,journal,profile,data['stop_stage'],tuple(attempts),tuple(replays),
            checks,first_failed,candidate=candidate,stop_kind=data['stop_kind'],persist=False)
        if outcome.as_bytes()!=content:raise ValueError
        for collection in (outcome.attempt_artifacts,outcome.replay_started_artifacts,outcome.replay_artifacts,
                           outcome.worker_evidence_artifacts,outcome.released_held_out_views):
            for binding in collection:
                if _journal_artifact_binding(journal,journal/binding['path'])!=binding:raise ValueError
        if _execution_context_bytes(root,plan,roster,resources,policies)!=context:raise ValueError
        return outcome
    except Exception:
        raise PublicTrack1CoordinatorError('stopped validation checkpoint is invalid') from None


_FINALIZATION_STAGES = {
    "outcome-publication", "consolidation", "acceptance-publication",
    "release-freeze", "release-publication", "manifest-publication", "final-publication", "packet-publication",
}


def record_finalization_failure(journal: str | Path, stage: str, outcome_sha256: str) -> None:
    if stage not in _FINALIZATION_STAGES or re.fullmatch(r"[0-9a-f]{64}", outcome_sha256) is None:
        raise PublicTrack1CoordinatorError("invalid finalization failure binding")
    content = _bytes({
        "schema": "genosapiens.track1-public-validation-finalization-failure",
        "classification": "public",
        "validation_outcome_sha256": outcome_sha256,
        "stage": stage,
        "error_code": "track1-" + stage + "-failed",
        "finalization_status": "failed",
        "controlled_use_authorized": False,
    })
    publish_validation_artifact(Path(journal) / "finalization-failures" / f"{uuid4().hex}.json", content)


def _checkpoint_json(content: bytes) -> dict:
    def pairs(items):
        value = {}
        for key, child in items:
            if key in value:
                raise ValueError
            value[key] = child
        return value
    if type(content) is not bytes or len(content) > 64_000_000:
        raise ValueError
    value = json.loads(content, object_pairs_hook=pairs)
    if type(value) is not dict or _bytes(value) != content:
        raise ValueError
    return value


def _read_checkpoint(journal: Path, relative: str) -> bytes:
    target = journal / relative
    if not target.is_relative_to(journal) or any(
        parent.is_symlink() for parent in (target, *target.parents)
    ):
        raise ValueError
    return target.read_bytes()


def _population_selector(policies):
    rule, _, _ = _population_template(policies.policy_bytes)
    return rule.subset, rule.stratum, rule.ancestry


def _recovery_resolver(journal, roster, resources, policies, *, original_stack=None):
    from genosapiens.wgs_validation_recovery import AuthenticatedValidationRecovery
    from genosapiens.wgs_shards import read_bounded_metadata
    derivation = json.loads(read_bounded_metadata(resources.repository_root /
        "config/track1/public-population-development-derivation.json"))
    probes = {*derivation["selected_low_identity_sha256s"], derivation["selected_high_identity_sha256"],
        *(_sha256(member.encode("ascii")) for case in roster.cases for member in case.truth_members)}
    def bound_resources(attempt):
        carried = _carried_context(journal, attempt.case_id)
        if carried is None:
            return resources if original_stack is None else replace(resources, executable_stack_content=original_stack)
        stack = _read_checkpoint(journal, "carried-forward-stack.json")
        if _sha256(stack) != carried["executable_stack_sha256"]:
            raise PublicTrack1CoordinatorError("carried validation stack differs")
        return replace(resources, executable_stack_content=stack)
    return AuthenticatedValidationRecovery(objects=journal / "scientific-objects", resources=resources,
        policies=policies, scratch_root=journal / "readback-scratch", cases=roster.cases,
        population_selector=_population_selector(policies), frequency_identities=probes,
        resources_for_attempt=bound_resources)


def _load_attempt(value: dict) -> ValidationAttempt:
    from genosapiens.evidence_model import _typed
    from genosapiens.track1_validation_runner import Track1ValidationRunnerError
    try:
        return _typed(value, ValidationAttempt, decode=True)
    except (ValueError, TypeError):
        raise Track1ValidationRunnerError("invalid validation attempt checkpoint") from None


def load_validation_outcome_checkpoint(
    checkpoint_root: str | Path,
    *,
    expected_sha256: str,
    repository_root: str | Path,
    plan: ValidationPlan,
    roster: Track1Roster,
    resources: PublicResourceContext,
    policies: PublicPolicyFamily,
) -> PublicCandidateValidation:
    """Authenticate a saved outcome, current sources and independently expected identity.

    Recompute evaluator products and candidate derivation from complete attempts.
    Reauthenticate complete saved preparation, annotation and scientific shards.
    Readback launches neither annotation nor another scientific invocation.
    """
    try:
        journal = Path(checkpoint_root).absolute()
        root = Path(repository_root).resolve(strict=True)
        if (
            re.fullmatch(r"[0-9a-f]{64}", expected_sha256) is None
            or not journal.is_relative_to(resources.scratch_root.resolve(strict=True))
            or any(parent.is_symlink() for parent in (journal, *journal.parents))
        ):
            raise ValueError
        context = _execution_context_bytes(root, plan, roster, resources, policies)
        if (journal/'validation-stopped-outcome.json').exists():
            raise PublicTrack1CoordinatorError('a stopped validation cannot be finalized')
        if _read_checkpoint(journal, "execution-context.json") != context:
            raise ValueError
        content = _read_checkpoint(journal, "validation-outcome.json")
        if _sha256(content) != expected_sha256:
            raise ValueError
        data = _checkpoint_json(content)
        if set(data) != {item.name for item in fields(PublicCandidateValidation)}:
            raise ValueError
        _verify_carried_forward(journal, data['carried_forward'])
        recovery_resolver = _recovery_resolver(journal, roster, resources, policies)
        development = tuple(_load_attempt(value) for value in data["development_experiment"]["attempts"])
        held_out = tuple(_load_attempt(value) for value in data["held_out_attempts"])
        attempts = (*development, *held_out)
        candidates = tuple(runtime.candidate for runtime in policies.runtimes)
        expected_keys = {
            (case.split, case.case_id, policy.policy_id, 1)
            for case in roster.cases if case.split == "development"
            for policy in candidates
        }
        if data["candidate"] is not None:
            expected_keys.update(
                (case.split, case.case_id, data["candidate"]["selected_policy_id"], 1)
                for case in roster.cases if case.split == "held_out"
            )
        keys = {(item.split, item.case_id, item.policy_id, item.attempt_number) for item in attempts}
        if len(keys) != len(attempts) or keys != expected_keys:
            raise ValueError
        experiment = complete_development_experiment(
            plan=plan, roster=roster, protocol=policies.protocol,
            candidates=candidates, attempts=development,
         recovery_resolver=recovery_resolver)
        if experiment.as_bytes() != _bytes(data["development_experiment"]):
            raise ValueError
        candidate = None
        derivation = {}
        evaluations = ()
        metrics = ()
        profile,_ = load_submission_acceptance_profile(root)
        development_gates_pass = all(all(passed for _, passed in _case_gate_checks(
            root,roster,next(case for case in roster.cases if case.case_id==attempt.case_id),
            (attempt,),tuple(tuple(row) for row in data['determinism_replays'] if row[1]==attempt.case_id),profile,
            journal=journal, recovery_resolver=recovery_resolver))
            for attempt in development)
        if experiment.selected_policy_id is not None and development_gates_pass:
            derivation, derivation_bytes, criteria = _development_derivation(experiment, root)
            candidate = _candidate(
                experiment, repository_root=root, resources=resources, policies=policies,
                derivation_bytes=derivation_bytes, criteria=criteria,
                stack_sha256=_sha256(executable_stack_bytes(root)),
             recovery_resolver=recovery_resolver)
            if _read_checkpoint(journal, "candidate.json") != candidate.as_bytes():
                raise ValueError
            evaluations = evaluate_attempts(roster, held_out, split="held_out", policy_id=candidate.selected_policy_id, recovery_resolver=recovery_resolver)
            metrics = tuple(sorted((name, _fraction(value)) for name, value in evaluation_metrics(evaluations).items()))
        replays = tuple(tuple(item) for item in data["determinism_replays"])
        if (len(replays) != len(attempts) or len({row[:3] for row in replays}) != len(attempts)
                or any((a.split,a.case_id,a.policy_id,a.execution_sha256) not in replays for a in attempts)):
            raise ValueError
        outcome = PublicCandidateValidation(
            _sha256(executable_stack_bytes(root)), experiment, derivation, candidate,
            held_out, evaluations, metrics, replays, data["status"], carried_forward=data["carried_forward"],
        )
        if outcome.as_bytes() != content:
            raise ValueError
        expected_attempt_files = {
            f"{item.split}.{item.case_id}.{item.policy_id}.{item.attempt_number}.json"
            for item in attempts
        }
        if {path.name for path in (journal / "attempts").iterdir()} != expected_attempt_files:
            raise ValueError
        replay_by_key = {item[:3]: item for item in replays}
        if len(replay_by_key) != len(replays):
            raise ValueError
        if {path.name for path in (journal/'determinism').iterdir()} != {
                f'{row[0]}.{row[1]}.{row[2]}.json' for row in replays}:
            raise ValueError
        for attempt in attempts:
            name = f"{attempt.split}.{attempt.case_id}.{attempt.policy_id}.{attempt.attempt_number}.json"
            saved = _checkpoint_json(_read_checkpoint(journal, "attempts/" + name))
            if saved != _checkpoint_json(_bytes({
                "schema": "genosapiens.track1-public-validation-attempt-checkpoint",
                "classification": "public", "validation_only": True,
                "execution_context_sha256": _attempt_context_sha256(journal, attempt.case_id), "attempt": asdict(attempt),
            })):
                raise ValueError
            case = next(case for case in roster.cases if case.case_id == attempt.case_id)
            runtime = next(item for item in policies.runtimes if item.candidate.policy_id == attempt.policy_id)
            view = development_ranker_view(roster, case.case_id, runner_resources=resources.common_runner_resources) if case.split == "development" else authenticated_held_out_ranker_view(
                candidate, experiment, roster, case.case_id, runner_resources=resources.common_runner_resources,
            )
            if (
                attempt.input_sha256s != tuple(item.sha256 for item in view.inputs)
                or attempt.resource_descriptor_sha256s != tuple(item.descriptor_sha256 for item in view.resources)
                or attempt.policy_bindings_sha256 != runtime.candidate.bindings.sha256
                or attempt.scope != case.scope
            ):
                raise ValueError
            if attempt.execution_sha256 is not None:
                _verify_saved_attempt(journal, attempt, case, view, runtime, resources, policies, recovery_resolver=recovery_resolver)
                replay = replay_by_key.get((attempt.split, attempt.case_id, attempt.policy_id))
                if replay is not None:
                    if replay[3] != attempt.execution_sha256 or _read_checkpoint(
                        journal, f"determinism/{replay[0]}.{replay[1]}.{replay[2]}.json",
                    ) != _bytes({
                        "schema": "genosapiens.track1-public-validation-determinism-checkpoint",
                        "classification": "public", "execution_context_sha256": _attempt_context_sha256(journal, attempt.case_id),
                        "candidate_sha256": None if attempt.split == "development" else candidate.sha256,
                        "replay": replay,
                    }):
                        raise ValueError

        if _execution_context_bytes(root, plan, roster, resources, policies) != context:
            raise ValueError
        return outcome
    except Exception:
        raise PublicTrack1CoordinatorError("invalid or stale validation outcome checkpoint") from None
