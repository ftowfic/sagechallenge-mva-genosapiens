"""Static Track 1 public validation-plan and truth-isolation contract.

This module validates the plan boundary used before development execution. It
does not derive policies, execute cases, inspect held-out results, create a
candidate-stack manifest, or authorize controlled use.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass
from enum import StrEnum
import hashlib
import json
from pathlib import Path, PurePosixPath
import re
from typing import Literal
from urllib.parse import quote

from genosapiens.bounded_json import (
    BoundedJSONError,
    load_bounded_json,
    read_bounded_bytes,
)
from genosapiens.json_schema import JsonSchemaError, validate_json_schema


PLAN_SCHEMA = "genosapiens.track1-base-release-validation-plan"
PIN_SCHEMA = "genosapiens.track1-base-release-validation-plan-pin"
ARTIFACT_BUNDLE_SCHEMA = "genosapiens.track1-base-release-validation-artifact-bundle"
CANONICAL_SCHEMA_PATH = "schemas/track1/base-release-validation-plan.schema.json"
CANONICAL_SCHEMA_SHA256 = (
    "060d61cdf7bd5c232309e9d31712d52e59dfa287a9c2a062cd95cfd986bac2a4"
)
CANONICAL_ARTIFACT_BUNDLE_SCHEMA_PATH = (
    "schemas/track1/base-release-validation-artifact-bundle.schema.json"
)
CANONICAL_ARTIFACT_BUNDLE_SCHEMA_SHA256 = (
    "9779067edc5f2912999aa1011ff4c36e513d1f048d84ae01d119861f173701f4"
)
CANONICAL_PIN_PATH = "config/track1/base-release-validation-plan-synthetic-pin.json"
CANONICAL_PIN_SHA256 = (
    "f8790c77c026bc5ea7e5eccdb48c293f0d7edd81d77ded5bbbbd1f31760a9204"
)
MAX_VALIDATION_ARTIFACT_BYTES = 2 * 1024 * 1024


def completed_context_inputs_from_plan(plan, *, plan_sha256):
    """Parse only the explicit public context restart; ordinary plans stay fresh."""
    from genosapiens.wgs_profile_projection import CompletedProfileDiscoveryInputs
    from genosapiens.wgs_ledger import canonical
    from genosapiens.track1_roster import validate_context_restart_acceptance
    value = plan.get('completed_context_inputs')
    if value is None:
        if plan.get('restart_stage') == 'scientific_publication' and plan.get('completed_scientific_inputs') is not None:
            return None
        if plan.get('restart_stage') not in (None, 'analysis_entry'):
            raise ValueError('context restart inputs are missing')
        return None
    if (plan.get('restart_stage') not in ('association_discovery', 'context_assembly', 'assessment_publication')
            or value['consumer_executable_stack_sha256'] != plan['executable_stack_sha256']
            or value['acceptance_amendment'] != plan['acceptance_amendment']
            or value['producer_candidate_sha256'] == plan['candidate_sha256']):
        raise ValueError('context restart plan is incompatible')
    result = CompletedProfileDiscoveryInputs(canonical(value), plan_sha256)
    if result.profile_only != (plan.get('restart_stage') == 'association_discovery'):
        raise ValueError('profile-only restart stage binding differs')
    validate_context_restart_acceptance(result)
    return result


def completed_assessment_inputs_from_plan(plan, *, plan_sha256):
    from genosapiens.wgs_assessment_checkpoint import CompletedAssessmentInputs
    from genosapiens.wgs_ledger import canonical
    value = plan.get('completed_assessment_inputs')
    if value is None:
        if plan.get('restart_stage') == 'assessment_publication':
            raise ValueError('assessment restart inputs are missing')
        return None
    if (plan.get('restart_stage') != 'assessment_publication'
            or plan.get('completed_context_inputs') is None
            or value['consumer_executable_stack_sha256'] != plan['executable_stack_sha256']
            or value['acceptance_amendment'] != plan['acceptance_amendment']):
        raise ValueError('assessment restart plan is incompatible')
    return CompletedAssessmentInputs(canonical(value), plan_sha256)


def completed_scientific_inputs_from_plan(plan, *, plan_sha256):
    from genosapiens.wgs_scientific_checkpoint import CompletedScientificInputs
    from genosapiens.wgs_ledger import canonical
    value = plan.get('completed_scientific_inputs')
    if value is None:
        if plan.get('restart_stage') == 'scientific_publication':
            raise ValueError('scientific restart inputs are missing')
        return None
    if (plan.get('restart_stage') != 'scientific_publication'
            or plan.get('completed_context_inputs') is not None
            or plan.get('completed_assessment_inputs') is not None
            or value['consumer_executable_stack_sha256'] != plan['executable_stack_sha256']
            or value['acceptance_amendment'] != plan['acceptance_amendment']):
        raise ValueError('scientific restart plan is incompatible')
    return CompletedScientificInputs(canonical(value), plan_sha256)


SHA256 = re.compile(r"^[0-9a-f]{64}$")
IDENTIFIER = re.compile(r"^[A-Za-z0-9][A-Za-z0-9._:-]*$")
DECIMAL = re.compile(r"^(?:0|[1-9][0-9]*)(?:\.[0-9]+)?$")
HPO_ID = re.compile(r"^HP:[0-9]{7}$")
DNA_SEQUENCE = re.compile(r"^[ACGTN]+$")
PUBLIC_BACKGROUND_SOURCE = {
    "contract_sha256": "79a1a5c64a587bcdf81924f6eff530ae7d8d2b1206a9d72aaa5fd31df6e2d21e",
    "manifest_sha256": "3a7901ba0814842fc799781916e3faf1213ac15ab45567f6b4758c90dcd403cb",
    "vcf_sha256": "758456587bdcf53e7aa56ce8abfbbeb65460d1f3c80ca5f10c9472810cd49e7a",
    "provenance_sha256": "3e0820dba8474c87c572ce786a46d00d4419d09ea5552041f423542824f0e966",
}
PUBLIC_BACKGROUND_RECORDS = {
    "unphased_insertion": (
        "public-01", "1", 602113, "T", ("TGCCCA",),
        "15f598aa7240d04822262d8d0a56b4372b4e12653de482e989d0a065bbffe461",
    ),
    "phased_snv": (
        "public-02", "1", 604358, "G", ("C",),
        "838676330d71289fd571b558d4dd412e0b93e403a4a34d0dd7e779bcfe00ff7d",
    ),
    "missing_quality_snv": (
        "public-03", "1", 604360, "T", ("C",),
        "bd4c4b0e7202a2a849363a523bff0d2821a1fbda32659c89f92728a9bd2865c5",
    ),
    "filtered_deletion": (
        "public-04", "1", 631859, "CG", ("C",),
        "5f437304b94030541a9c44c4d487a07d4c59bba2144cadd1a097749fe2e9e8c9",
    ),
    "multiallelic": (
        "public-05", "1", 1735008, "G", ("GGTTT", "GTTTT"),
        "797d2fc60e07d0860c9725cb9d8bceba7db03d1f5b059c5c76107feac381638c",
    ),
}
REFERENCE_MANIFEST_SHA256 = (
    "27eaf00fde35e5fa83f945c63e5dee871960374d90488cc0797b5ebea2fb1c25"
)
SOURCE_FORMAT_FIELDS = {
    "schema",
    "input_sha256",
    "reference_role",
    "reference_identity_id",
    "reference_manifest_sha256",
    "mitochondrial_contigs",
    "repeat_identity",
}
OPAQUE_SEED_ID = re.compile(r"^seed-[0-9a-f]{16}$")
SYNTHETIC_VCF_METADATA_KEYS = {
    "fileformat", "reference", "source", "genosapiens_seed_sha256",
    "contig", "FILTER", "INFO", "FORMAT",
}
SYNTHETIC_VCF_INFO_FIELDS = {
    "SVTYPE", "END", "CIPOS", "CIEND", "REPID", "RU",
    "START", "MOTIF", "PERIOD", "INEXACT_ALLELE", "BPDIFFS",
}
SYNTHETIC_VCF_FORMAT_FIELDS = {
    "GT", "GQ", "DP", "AD", "PL", "PS", "CN", "REPCN", "GB", "Q",
}

REQUIRED_SCOPES = (
    "compound_heterozygous",
    "dominant",
    "homozygous_recessive",
    "incomplete_recessive",
    "mitochondrial",
    "mixed_class",
    "mosaic",
    "regulatory",
    "repeat",
    "sv_cnv",
    "x_linked",
)
TRUTH_SCOPE_MEMBER_COUNTS = {
    scope: 2 if scope in {"compound_heterozygous", "mixed_class"} else 1
    for scope in REQUIRED_SCOPES
}
NON_SUCCESS_STATUSES = ("failed", "timed_out", "truncated", "unsupported")
NON_SUCCESS_ACCOUNTING = "retained-in-full-roster-denominator"
STATUS_ACCOUNTING = "reported-by-status-including-non-success"
NULL_MEASUREMENT_ACCOUNTING = (
    "status-retained-measurement-null-assessed-subset-separate"
)
ASSESSED_SUBSET_ACCOUNTING = "publish-denominator-and-exclusions-separately"
REQUIRED_ANALYSIS_MODES = ("primary",)
ADMISSION_DIMENSIONS = (
    "admission_outcome",
    "admission_route",
    "case_id",
    "case_status",
    "policy_variant",
    "scope_support",
)
RANKING_DIMENSIONS = (
    "analysis_mode",
    "case_id",
    "case_status",
    "policy_variant",
    "scope_support",
)
CASE_DIMENSIONS = ("case_id", "case_status", "policy_variant", "scope_support")
RESOURCE_DIMENSIONS = ("case_id", "case_status", "policy_variant")
METRIC_CONTRACTS = {
    "admission_causal_allele_recall": (
        "fraction",
        "truth-causal-allele-identities-admitted-by-declared-route-or-union",
        "truth-causal-allele-identities-in-complete-split-roster",
        ADMISSION_DIMENSIONS,
    ),
    "admission_outcome_count": (
        "count",
        "candidate-admission-assessments-by-declared-outcome",
        "one",
        ADMISSION_DIMENSIONS,
    ),
    "candidate_burden": (
        "count",
        "admitted-candidate-identities",
        "one",
        ADMISSION_DIMENSIONS,
    ),
    "candidate_generation_recall": (
        "fraction",
        "truth-causal-allele-identities-generated",
        "truth-causal-allele-identities-in-complete-split-roster",
        CASE_DIMENSIONS,
    ),
    "challenge_f_max": (
        "fraction",
        "pinned-evaluator-best-f-score-sum",
        "cases-in-complete-split-roster",
        RANKING_DIMENSIONS,
    ),
    "challenge_rank_points": (
        "points-per-case",
        "pinned-evaluator-rank-points-sum",
        "cases-in-complete-split-roster",
        RANKING_DIMENSIONS,
    ),
    "context_recovery": (
        "fraction",
        "truth-relevant-candidate-contexts-recovered",
        "truth-relevant-candidate-contexts-in-complete-split-roster",
        ADMISSION_DIMENSIONS,
    ),
    "cpu_time": (
        "cpu-seconds",
        "measured-cpu-seconds-or-null",
        "one",
        RESOURCE_DIMENSIONS,
    ),
    "disk_usage": (
        "bytes",
        "measured-peak-disk-bytes-or-null",
        "one",
        RESOURCE_DIMENSIONS,
    ),
    "enumeration_coverage": (
        "fraction",
        "completed-required-enumeration-assessments",
        "required-enumeration-assessments-in-complete-split-roster",
        ADMISSION_DIMENSIONS,
    ),
    "exact_full_pair_recovery": (
        "fraction",
        "truth-full-pair-cases-with-exact-pair-generated",
        "truth-full-pair-cases-in-complete-split-roster",
        RANKING_DIMENSIONS,
    ),
    "exact_singleton_recovery": (
        "fraction",
        "truth-singleton-cases-with-exact-singleton-generated",
        "truth-singleton-cases-in-complete-split-roster",
        RANKING_DIMENSIONS,
    ),
    "full_pair_member_availability": (
        "fraction",
        "truth-full-pair-member-identities-generated",
        "truth-full-pair-member-identities-in-complete-split-roster",
        ADMISSION_DIMENSIONS,
    ),
    "peak_memory": (
        "bytes",
        "maximum-measured-resident-set-bytes-or-null",
        "one",
        RESOURCE_DIMENSIONS,
    ),
    "repeat_locus_recovery": (
        "fraction",
        "truth-repeat-locus-cases-with-matched-locus-hypothesis",
        "truth-repeat-locus-cases-in-complete-split-roster",
        RANKING_DIMENSIONS,
    ),
    "repeat_projected_allele_recovery": (
        "fraction",
        "truth-repeat-projected-allele-identities-generated",
        "truth-repeat-projected-allele-identities-in-complete-split-roster",
        RANKING_DIMENSIONS,
    ),
    "top_1_full_genotype_recovery": (
        "fraction",
        "truth-genotype-cases-ranked-exactly-first",
        "truth-genotype-cases-in-complete-split-roster",
        RANKING_DIMENSIONS,
    ),
    "top_10_full_genotype_recovery": (
        "fraction",
        "truth-genotype-cases-ranked-in-exact-top-ten",
        "truth-genotype-cases-in-complete-split-roster",
        RANKING_DIMENSIONS,
    ),
    "wall_time": (
        "wall-seconds",
        "measured-wall-seconds-or-null",
        "one",
        RESOURCE_DIMENSIONS,
    ),
}
STATUS_ACCOUNTING_METRICS = {"admission_outcome_count"}
NULL_MEASUREMENT_METRICS = {
    "candidate_burden",
    "cpu_time",
    "disk_usage",
    "peak_memory",
    "wall_time",
}


class ValidationPlanError(RuntimeError):
    """Raised when the public validation plan is not replayable or isolated."""


@dataclass(frozen=True, slots=True)
class ValidationPlanPin:
    plan_id: str
    canonical_plan_path: str
    canonical_plan_sha256: str
    canonical_schema_path: str
    canonical_schema_sha256: str
    canonical_artifact_bundle_path: str
    canonical_artifact_bundle_sha256: str
    canonical_artifact_bundle_schema_path: str
    canonical_artifact_bundle_schema_sha256: str
    schema: Literal["genosapiens.track1-base-release-validation-plan-pin"] = (
        PIN_SCHEMA
    )


class Split(StrEnum):
    DEVELOPMENT = "development"
    HELD_OUT = "held_out"


class InputKind(StrEnum):
    HPO_TEXT = "hpo_text"
    VCF = "vcf"
    NATIVE = "native"


class OracleRole(StrEnum):
    GENOMIC_TRUTH = "genomic_truth"
    ANNOTATION_EXPECTATION = "annotation_expectation"
    CAUSAL_TRUTH = "causal_truth"
    COMPARISON_BASELINE = "comparison_baseline"


ORACLE_SOURCES = {
    OracleRole.GENOMIC_TRUTH: "giab_hg002",
    OracleRole.ANNOTATION_EXPECTATION: "release_annotation_expectation",
    OracleRole.CAUSAL_TRUTH: "synthetic_causal_truth",
    OracleRole.COMPARISON_BASELINE: "exomiser_comparison",
}
ORACLE_RESOURCE_KINDS = {
    OracleRole.GENOMIC_TRUTH: "genomic_truth",
    OracleRole.ANNOTATION_EXPECTATION: "annotation_expectations",
    OracleRole.CAUSAL_TRUTH: "causal_truth",
    OracleRole.COMPARISON_BASELINE: "comparison_baseline",
}
RUNNER_RESOURCE_KINDS = {"runner_annotation_resources", "runner_reference"}
RESOURCE_KINDS = RUNNER_RESOURCE_KINDS | set(ORACLE_RESOURCE_KINDS.values()) | {
    "challenge_evaluator", "submission_acceptance_profile"
}


@dataclass(frozen=True, slots=True)
class GeneratorBinding:
    generator_id: str
    version: str
    manifest_sha256: str


@dataclass(frozen=True, slots=True)
class InputBinding:
    input_id: str
    kind: InputKind
    sha256: str
    generator_id: str
    seed: int


@dataclass(frozen=True, slots=True)
class ResourceBinding:
    resource_id: str
    kind: str
    descriptor_sha256: str


@dataclass(frozen=True, slots=True)
class BackgroundBinding:
    background_id: str
    subject_id: str
    kind: Literal["HG002", "independent_synthetic"]
    shared_across_cases: bool
    disclosure: str


@dataclass(frozen=True, slots=True)
class OracleBinding:
    oracle_id: str
    role: OracleRole
    source_kind: str
    resource_id: str


@dataclass(frozen=True, slots=True)
class ValidationCase:
    case_id: str
    split: Split
    recipe_group_id: str
    recipe_manifest_sha256: str
    perturbation_group_id: str
    perturbation_manifest_sha256: str
    generator_id: str
    seed: int
    background_id: str
    runner_input_ids: tuple[str, ...]
    runner_resource_ids: tuple[str, ...]
    causal_truth_oracle_id: str
    assessed_scopes: tuple[str, ...]


@dataclass(frozen=True, slots=True)
class MetricDefinition:
    metric_id: str
    unit: str
    numerator: str
    denominator: str
    reporting_dimensions: tuple[str, ...]
    non_success_statuses: tuple[str, ...]
    non_success_accounting: str


@dataclass(frozen=True, slots=True)
class AcceptanceCriterion:
    metric_id: str
    comparator: Literal[">=", "<="]
    threshold: str
    derivation_sha256: str


@dataclass(frozen=True, slots=True)
class Lifecycle:
    state: Literal["development_plan", "held_out_ready", "adaptive_revision"]
    held_out_execution: Literal["blocked", "ready", "adaptive_only"]
    independence_status: Literal["not_evaluated", "prespecified", "invalidated"]
    development_derivation_sha256: str | None
    candidate_stack_sha256: str | None
    post_held_out_policy_change_sha256: str | None
    acceptance_criteria: tuple[AcceptanceCriterion, ...]


@dataclass(frozen=True, slots=True)
class ValidationPlan:
    plan_id: str
    submission_acceptance_profile_sha256: str
    analysis_profile: str
    analysis_profile_policy_sha256: str
    profile_acceptance_sha256: str
    lifecycle: Lifecycle
    evaluator_resource_id: str
    required_scopes: tuple[str, ...]
    required_analysis_modes: tuple[str, ...]
    generators: tuple[GeneratorBinding, ...]
    inputs: tuple[InputBinding, ...]
    resources: tuple[ResourceBinding, ...]
    backgrounds: tuple[BackgroundBinding, ...]
    oracles: tuple[OracleBinding, ...]
    cases: tuple[ValidationCase, ...]
    metrics: tuple[MetricDefinition, ...]
    assessed_subset_accounting: str
    limitations: tuple[str, ...]
    schema: Literal["genosapiens.track1-base-release-validation-plan"] = PLAN_SCHEMA
    classification: Literal["public"] = "public"
    validation_only: Literal[True] = True
    controlled_use_authorized: Literal[False] = False


@dataclass(frozen=True, slots=True)
class RankerInput:
    kind: InputKind
    sha256: str


@dataclass(frozen=True, slots=True)
class RankerResource:
    kind: str
    descriptor_sha256: str


@dataclass(frozen=True, slots=True)
class RankerCaseView:
    inputs: tuple[RankerInput, ...]
    resources: tuple[RankerResource, ...]


@dataclass(frozen=True, slots=True)
class ValidationArtifact:
    artifact_id: str
    binding_id: str | None
    parent_binding_id: str | None
    audience: Literal["evaluator", "runner"]
    path: str
    size_bytes: int
    sha256: str


@dataclass(frozen=True, slots=True)
class ValidationArtifactBundle:
    plan_id: str
    plan_sha256: str
    artifacts: tuple[ValidationArtifact, ...]
    schema: Literal[
        "genosapiens.track1-base-release-validation-artifact-bundle"
    ] = ARTIFACT_BUNDLE_SCHEMA
    classification: Literal["public"] = "public"
    validation_only: Literal[True] = True


def _fail(message: str = "invalid Track 1 validation plan") -> None:
    raise ValidationPlanError(message)


def _sha(value: object) -> None:
    if type(value) is not str or SHA256.fullmatch(value) is None:
        _fail()


def _identifier(value: object) -> None:
    if type(value) is not str or IDENTIFIER.fullmatch(value) is None:
        _fail()


def _canonical_bytes(value: object) -> bytes:
    return (
        json.dumps(value, sort_keys=True, separators=(",", ":"), ensure_ascii=True, allow_nan=False)
        + "\n"
    ).encode("ascii")


def validation_plan_as_json(plan: ValidationPlan) -> bytes:
    return _canonical_bytes(asdict(plan))


def ranker_case_view_as_json(view: RankerCaseView) -> bytes:
    return _canonical_bytes(asdict(view))


def validation_artifact_bundle_as_json(bundle: ValidationArtifactBundle) -> bytes:
    return _canonical_bytes(asdict(bundle))


def _ordered_unique(values, *, key, what: str) -> None:
    keys = tuple(key(item) for item in values)
    if not values or keys != tuple(sorted(keys)) or len(keys) != len(set(keys)):
        _fail(f"validation plan {what} must be nonempty, unique, and ordered")


def _parse(document: dict[str, object]) -> ValidationPlan:
    lifecycle_document = document["lifecycle"]
    lifecycle = Lifecycle(
        state=lifecycle_document["state"],
        held_out_execution=lifecycle_document["held_out_execution"],
        independence_status=lifecycle_document["independence_status"],
        development_derivation_sha256=lifecycle_document["development_derivation_sha256"],
        candidate_stack_sha256=lifecycle_document["candidate_stack_sha256"],
        post_held_out_policy_change_sha256=lifecycle_document[
            "post_held_out_policy_change_sha256"
        ],
        acceptance_criteria=tuple(
            AcceptanceCriterion(**item)
            for item in lifecycle_document["acceptance_criteria"]
        ),
    )
    return ValidationPlan(
        plan_id=document["plan_id"],
        submission_acceptance_profile_sha256=document["submission_acceptance_profile_sha256"],
        analysis_profile=document["analysis_profile"],
        analysis_profile_policy_sha256=document["analysis_profile_policy_sha256"],
        profile_acceptance_sha256=document["profile_acceptance_sha256"],
        lifecycle=lifecycle,
        evaluator_resource_id=document["evaluator_resource_id"],
        required_scopes=tuple(document["required_scopes"]),
        required_analysis_modes=tuple(document["required_analysis_modes"]),
        generators=tuple(GeneratorBinding(**item) for item in document["generators"]),
        inputs=tuple(
            InputBinding(
                input_id=item["input_id"],
                kind=InputKind(item["kind"]),
                sha256=item["sha256"],
                generator_id=item["generator_id"],
                seed=item["seed"],
            )
            for item in document["inputs"]
        ),
        resources=tuple(ResourceBinding(**item) for item in document["resources"]),
        backgrounds=tuple(BackgroundBinding(**item) for item in document["backgrounds"]),
        oracles=tuple(
            OracleBinding(
                oracle_id=item["oracle_id"],
                role=OracleRole(item["role"]),
                source_kind=item["source_kind"],
                resource_id=item["resource_id"],
            )
            for item in document["oracles"]
        ),
        cases=tuple(
            ValidationCase(
                case_id=item["case_id"],
                split=Split(item["split"]),
                recipe_group_id=item["recipe_group_id"],
                recipe_manifest_sha256=item["recipe_manifest_sha256"],
                perturbation_group_id=item["perturbation_group_id"],
                perturbation_manifest_sha256=item[
                    "perturbation_manifest_sha256"
                ],
                generator_id=item["generator_id"],
                seed=item["seed"],
                background_id=item["background_id"],
                runner_input_ids=tuple(item["runner_input_ids"]),
                runner_resource_ids=tuple(item["runner_resource_ids"]),
                causal_truth_oracle_id=item["causal_truth_oracle_id"],
                assessed_scopes=tuple(item["assessed_scopes"]),
            )
            for item in document["cases"]
        ),
        metrics=tuple(
            MetricDefinition(
                metric_id=item["metric_id"],
                unit=item["unit"],
                numerator=item["numerator"],
                denominator=item["denominator"],
                reporting_dimensions=tuple(item["reporting_dimensions"]),
                non_success_statuses=tuple(item["non_success_statuses"]),
                non_success_accounting=item["non_success_accounting"],
            )
            for item in document["metrics"]
        ),
        assessed_subset_accounting=document["assessed_subset_accounting"],
        limitations=tuple(document["limitations"]),
        schema=document["schema"],
        classification=document["classification"],
        validation_only=document["validation_only"],
        controlled_use_authorized=document["controlled_use_authorized"],
    )


def _parse_artifact_bundle(document: dict[str, object]) -> ValidationArtifactBundle:
    return ValidationArtifactBundle(
        plan_id=document["plan_id"],
        plan_sha256=document["plan_sha256"],
        artifacts=tuple(ValidationArtifact(**item) for item in document["artifacts"]),
        schema=document["schema"],
        classification=document["classification"],
        validation_only=document["validation_only"],
    )


def _expected_artifact_bindings(
    plan: ValidationPlan,
) -> dict[str, tuple[str, Literal["evaluator", "runner"]]]:
    expected = {
        f"generator:{item.generator_id}": (item.manifest_sha256, "evaluator")
        for item in plan.generators
    }
    expected.update(
        {
            f"input:{item.input_id}": (item.sha256, "runner")
            for item in plan.inputs
        }
    )
    expected.update(
        {
            f"resource:{item.resource_id}": (
                item.descriptor_sha256,
                "runner" if item.kind in RUNNER_RESOURCE_KINDS else "evaluator",
            )
            for item in plan.resources
        }
    )
    for case in plan.cases:
        expected[f"case:{case.case_id}:recipe"] = (
            case.recipe_manifest_sha256,
            "evaluator",
        )
        expected[f"case:{case.case_id}:perturbation"] = (
            case.perturbation_manifest_sha256,
            "evaluator",
        )
    return expected


def _validate_artifact_bundle(
    bundle: ValidationArtifactBundle,
    *,
    plan: ValidationPlan,
    plan_sha256: str,
) -> None:
    if (
        type(bundle) is not ValidationArtifactBundle
        or bundle.schema != ARTIFACT_BUNDLE_SCHEMA
        or bundle.classification != "public"
        or bundle.validation_only is not True
        or bundle.plan_id != plan.plan_id
        or bundle.plan_sha256 != plan_sha256
    ):
        _fail("validation artifact bundle does not bind the canonical plan")
    _sha(bundle.plan_sha256)
    if type(bundle.artifacts) is not tuple or any(
        type(item) is not ValidationArtifact for item in bundle.artifacts
    ):
        _fail()
    _ordered_unique(
        bundle.artifacts,
        key=lambda item: item.artifact_id,
        what="artifacts",
    )

    expected = _expected_artifact_bindings(plan)
    primary: dict[str, ValidationArtifact] = {}
    paths: set[str] = set()
    hashes: set[str] = set()
    for artifact in bundle.artifacts:
        _identifier(artifact.artifact_id)
        _sha(artifact.sha256)
        if artifact.binding_id is not None:
            _identifier(artifact.binding_id)
        if artifact.parent_binding_id is not None:
            _identifier(artifact.parent_binding_id)
        relative = PurePosixPath(artifact.path) if type(artifact.path) is str else None
        if (
            artifact.audience not in {"evaluator", "runner"}
            or type(artifact.size_bytes) is not int
            or not 0 < artifact.size_bytes <= MAX_VALIDATION_ARTIFACT_BYTES
            or relative is None
            or relative.is_absolute()
            or relative == PurePosixPath(".")
            or ".." in relative.parts
            or "\\" in artifact.path
            or relative.as_posix() != artifact.path
            or artifact.path in paths
            or artifact.sha256 in hashes
            or (artifact.binding_id is None) == (artifact.parent_binding_id is None)
        ):
            _fail("validation artifact binding is invalid")
        paths.add(artifact.path)
        hashes.add(artifact.sha256)
        if artifact.binding_id is not None:
            if artifact.binding_id in primary:
                _fail("validation artifact plan binding is duplicated")
            primary[artifact.binding_id] = artifact

    if set(primary) != set(expected):
        _fail("validation artifact bundle does not resolve every plan binding")
    for binding_id, (expected_sha256, expected_audience) in expected.items():
        artifact = primary[binding_id]
        if (artifact.sha256, artifact.audience) != (
            expected_sha256,
            expected_audience,
        ):
            _fail("validation artifact differs from its plan binding")
    for artifact in bundle.artifacts:
        if artifact.parent_binding_id is not None:
            parent = primary.get(artifact.parent_binding_id)
            if parent is None or artifact.audience != parent.audience:
                _fail("validation auxiliary artifact has an invalid parent")


def _load_artifact_json(
    verified_content: dict[str, bytes],
    artifact: ValidationArtifact,
) -> dict[str, object]:
    def reject_duplicate_keys(pairs):
        value = {}
        for key, child in pairs:
            if key in value:
                _fail("validation artifact descriptor contains duplicate fields")
            value[key] = child
        return value

    try:
        document = json.loads(
            verified_content[artifact.path],
            object_pairs_hook=reject_duplicate_keys,
        )
    except (KeyError, UnicodeDecodeError, json.JSONDecodeError, RecursionError):
        _fail("validation artifact descriptor is invalid")
    if not isinstance(document, dict):
        _fail("validation artifact descriptor must be a JSON object")
    return document


def _closed_descriptor(
    document: dict[str, object],
    *,
    fields: set[str],
    schema: str,
) -> None:
    if (
        set(document) != fields
        or document.get("schema") != schema
        or document.get("classification") != "public"
        or document.get("validation_only") is not True
    ):
        _fail("validation artifact descriptor fields are invalid")


def _synthetic_vcf_truth_identities(content: bytes) -> set[str]:
    try:
        lines = content.decode("ascii").splitlines()
    except UnicodeDecodeError:
        _fail("validation synthetic VCF is not ASCII")
    if any(marker in content.upper() for marker in (b"CAUSAL", b"TRUTH", b"MODEL")):
        _fail("validation synthetic VCF publishes a truth label")
    identities = set()
    for line in lines:
        if line.startswith("##"):
            key, separator, _ = line[2:].partition("=")
            if not separator or key not in SYNTHETIC_VCF_METADATA_KEYS:
                _fail("validation synthetic VCF metadata is invalid")
            continue
        if line.startswith("#CHROM"):
            continue
        columns = line.split("\t")
        if (
            len(columns) != 10
            or (columns[2] != "." and OPAQUE_SEED_ID.fullmatch(columns[2]) is None)
        ):
            _fail("validation synthetic VCF record is invalid")
        identities.update(
            f"{columns[0]}:{columns[1]}:{columns[3]}:{alternate}"
            for alternate in columns[4].split(",")
        )
        info = {}
        if columns[7] != ".":
            for field in columns[7].split(";"):
                key, separator, value = field.partition("=")
                if not separator or not key or key in info:
                    _fail("validation synthetic VCF INFO is invalid")
                info[key] = value
        if set(info) - SYNTHETIC_VCF_INFO_FIELDS:
            _fail("validation synthetic VCF INFO is invalid")
        format_fields = columns[8].split(":")
        if (
            not format_fields
            or len(format_fields) != len(set(format_fields))
            or set(format_fields) - SYNTHETIC_VCF_FORMAT_FIELDS
        ):
            _fail("validation synthetic VCF FORMAT is invalid")
        if "REPID" in info:
            identities.add(info["REPID"])
        if {"START", "END", "MOTIF"}.issubset(info):
            identities.add(
                columns[2]
                if columns[2] != "." and ";" not in columns[2]
                else "longtr-region:"
                + quote(columns[0], safe="")
                + f":{info['START']}-{info['END']}"
            )
    return identities


def _synthetic_background_identities(
    bundle: ValidationArtifactBundle,
    verified_content: dict[str, bytes],
) -> set[str]:
    artifacts = [
        artifact
        for artifact in bundle.artifacts
        if artifact.artifact_id == "generator-public-background"
    ]
    if len(artifacts) != 1:
        _fail("validation generator public background is missing")
    document = _load_artifact_json(verified_content, artifacts[0])
    _closed_descriptor(
        document,
        fields={
            "schema",
            "classification",
            "validation_only",
            "background_id",
            "subject_id",
            "controlled_inputs_used",
            "source_fixture",
            "records",
        },
        schema="genosapiens.track1-synthetic-public-background",
    )
    source = document.get("source_fixture")
    records = document.get("records")
    if (
        document.get("background_id") != "hg002-public-background"
        or document.get("subject_id") != "HG002"
        or document.get("controlled_inputs_used") is not False
        or not isinstance(source, dict)
        or set(source)
        != {
            "resource_id",
            "release",
            "contract_sha256",
            "manifest_sha256",
            "vcf_sha256",
            "provenance_sha256",
        }
        or source.get("resource_id")
        != "challenge-shaped-public-small-variant-fixture-v2"
        or source.get("release") != "v2"
        or any(source.get(field) != value for field, value in PUBLIC_BACKGROUND_SOURCE.items())
        or not isinstance(records, list)
        or len(records) != 5
        or any(not isinstance(record, dict) for record in records)
    ):
        _fail("validation generator public background is invalid")
    identities = set()
    record_ids = set()
    profiles = set()
    for record in records:
        alts = record.get("alts")
        if (
            set(record)
            != {
                "fixture_record_id",
                "chrom",
                "position",
                "ref",
                "alts",
                "profile",
                "source_record_sha256",
            }
            or type(record.get("fixture_record_id")) is not str
            or not record["fixture_record_id"]
            or record["fixture_record_id"] in record_ids
            or type(record.get("chrom")) is not str
            or type(record.get("position")) is not int
            or record["position"] <= 0
            or type(record.get("ref")) is not str
            or DNA_SEQUENCE.fullmatch(record["ref"]) is None
            or not isinstance(alts, list)
            or not alts
            or any(
                type(alt) is not str or DNA_SEQUENCE.fullmatch(alt) is None
                for alt in alts
            )
            or type(record.get("profile")) is not str
            or not record["profile"]
            or record["profile"] in profiles
            or (
                record["fixture_record_id"],
                record["chrom"],
                record["position"],
                record["ref"],
                tuple(alts),
                record.get("source_record_sha256"),
            )
            != PUBLIC_BACKGROUND_RECORDS.get(record.get("profile"))
        ):
            _fail("validation generator public background record is invalid")
        new_identities = {
            f"{record['chrom']}:{record['position']}:{record['ref']}:{alt}"
            for alt in alts
        }
        if identities & new_identities:
            _fail("validation generator public background allele is duplicated")
        identities.update(new_identities)
        record_ids.add(record["fixture_record_id"])
        profiles.add(record["profile"])
    if profiles != {
        "filtered_deletion",
        "missing_quality_snv",
        "multiallelic",
        "phased_snv",
        "unphased_insertion",
    }:
        _fail("validation generator public background profiles are incomplete")
    return identities


def _validate_primary_descriptors(
    plan: ValidationPlan,
    bundle: ValidationArtifactBundle,
    verified_content: dict[str, bytes],
) -> None:
    primary = {
        artifact.binding_id: artifact
        for artifact in bundle.artifacts
        if artifact.binding_id is not None
    }
    from genosapiens.track1_roster import SUBMISSION_PROFILE_SHA256
    acceptance = primary["resource:submission-acceptance-profile"]
    if (acceptance.path != "config/track1/submission-acceptance-profile.json"
            or acceptance.sha256 != SUBMISSION_PROFILE_SHA256
            or hashlib.sha256(verified_content[acceptance.path]).hexdigest() != SUBMISSION_PROFILE_SHA256):
        _fail("validation submission profile is stale or substituted")
    for generator in plan.generators:
        document = _load_artifact_json(
            verified_content,
            primary[f"generator:{generator.generator_id}"],
        )
        _closed_descriptor(
            document,
            fields={
                "schema",
                "classification",
                "validation_only",
                "generator_id",
                "version",
                "materialization",
                "controlled_sources_used",
                "executable_binding",
                "public_background_binding",
                "public_background_contract_binding",
                "generated_outputs",
                "determinism",
                "limitations",
            },
            schema="genosapiens.track1-synthetic-case-generator",
        )
        if (
            document.get("generator_id") != generator.generator_id
            or document.get("version") != generator.version
            or document.get("materialization") != "executable-seed-generator"
            or document.get("controlled_sources_used") is not False
            or not isinstance(document.get("executable_binding"), dict)
            or not isinstance(document.get("public_background_binding"), dict)
            or not isinstance(document.get("public_background_contract_binding"), dict)
            or not isinstance(document.get("generated_outputs"), list)
            or not isinstance(document.get("limitations"), list)
            or not document.get("limitations")
            or any(
                type(item) is not str or not item
                for item in document.get("limitations", [])
            )
        ):
            _fail("validation generator descriptor differs from the plan")

    input_roles = {
        InputKind.HPO_TEXT: "plain-hpo-text",
        InputKind.VCF: "synthetic-overlay-vcf",
        InputKind.NATIVE: "native-input-manifest",
    }
    plan_inputs = {item.input_id: item for item in plan.inputs}
    for case in plan.cases:
        recipe = _load_artifact_json(
            verified_content,
            primary[f"case:{case.case_id}:recipe"],
        )
        _closed_descriptor(
            recipe,
            fields={
                "schema",
                "classification",
                "validation_only",
                "case_id",
                "split",
                "recipe_group_id",
                "generator_id",
                "seed",
                "background_id",
                "input_roles",
                "hpo_terms",
                "output_profile",
                "scope",
            },
            schema="genosapiens.track1-synthetic-case-recipe",
        )
        hpo_terms = recipe.get("hpo_terms")
        if (
            (
                recipe.get("case_id"),
                recipe.get("split"),
                recipe.get("recipe_group_id"),
                recipe.get("generator_id"),
                recipe.get("seed"),
                recipe.get("background_id"),
            )
            != (
                case.case_id,
                case.split,
                case.recipe_group_id,
                case.generator_id,
                case.seed,
                case.background_id,
            )
            or recipe.get("input_roles")
            != [input_roles[plan_inputs[item].kind] for item in case.runner_input_ids]
            or not isinstance(hpo_terms, list)
            or not hpo_terms
            or len(hpo_terms) != len(set(hpo_terms))
            or any(
                type(term) is not str or HPO_ID.fullmatch(term) is None
                for term in hpo_terms
            )
            or recipe.get("output_profile") != "all-class-seed-v2"
            or recipe.get("scope") != "materialization-seed-only"
        ):
            _fail("validation recipe descriptor differs from the plan")
        hpo_inputs = [
            plan_inputs[input_id]
            for input_id in case.runner_input_ids
            if plan_inputs[input_id].kind is InputKind.HPO_TEXT
        ]
        if len(hpo_inputs) != 1:
            _fail("validation recipe HPO input is invalid")
        hpo_artifact = primary[f"input:{hpo_inputs[0].input_id}"]
        expected_hpo = ("\n".join(hpo_terms) + "\n").encode("ascii")
        if verified_content[hpo_artifact.path] != expected_hpo:
            _fail("validation recipe HPO terms differ from runner input")

        perturbations = _load_artifact_json(
            verified_content,
            primary[f"case:{case.case_id}:perturbation"],
        )
        _closed_descriptor(
            perturbations,
            fields={
                "schema",
                "classification",
                "validation_only",
                "case_id",
                "split",
                "perturbation_group_id",
                "operations",
            },
            schema="genosapiens.track1-synthetic-case-perturbations",
        )
        operations = perturbations.get("operations")
        if (
            (
                perturbations.get("case_id"),
                perturbations.get("split"),
                perturbations.get("perturbation_group_id"),
            )
            != (case.case_id, case.split, case.perturbation_group_id)
            or not isinstance(operations, list)
            or not operations
            or len(operations) != len(set(operations))
            or any(type(item) is not str or not item for item in operations)
        ):
            _fail("validation perturbation descriptor differs from the plan")

    resources = {item.resource_id: item for item in plan.resources}
    descriptor_fields = {
        "annotation-expectations": {
            "schema", "classification", "validation_only", "resource_id",
            "kind", "status", "required_release", "expectations_path",
            "expectations_sha256", "roster_path", "roster_sha256",
            "validation_report_path", "validation_report_sha256", "limitation",
        },
        "exomiser-baseline": {
            "schema", "classification", "validation_only", "resource_id",
            "kind", "status", "version", "data_release", "manifest_path",
            "manifest_sha256", "validation_report_path",
            "validation_report_sha256", "limitation",
        },
        "public-annotations": {
            "schema", "classification", "validation_only", "resource_id",
            "kind", "status", "configuration_bindings",
            "installed_manifest_bindings", "output_bindings",
            "unresolved_resource_ids", "limitation",
        },
        "public-reference": {
            "schema", "classification", "validation_only", "resource_id",
            "kind", "status", "configuration_binding",
            "installed_manifest_binding", "output_bindings", "limitation",
        },
    }
    descriptor_values = {
        "annotation-expectations": {
            "status": "materialized_and_validated",
            "required_release": "Ensembl-VEP-116-GRCh38",
            "expectations_path": "tests/fixtures/track1/roster/annotation-expectations.json",
            "expectations_sha256": "1e19d9c40dfac9c331a400662cc324f2c8e6002531f279ccd6c9467d2b4d8a65",
            "roster_path": "tests/fixtures/track1/roster/roster.json",
            "roster_sha256": "a5f71d2696842dae4007a30c078a16ead6b84ed7118ec90cd07fcda34a044cf3",
            "validation_report_path": "artifacts/validation/track1-model-roster-annotation-2026-09-15.json",
            "validation_report_sha256": "e1209abaa05fe800e574b822da2d46d356f6625dfdbd06b83cc9ac9ab06ebe2f",
        },
        "exomiser-baseline": {
            "status": "provisioned_with_historical_comparison",
            "version": "15.1.0",
            "data_release": "2512",
            "manifest_path": "config/public-resources/exomiser-15.1.0-2512.json",
            "manifest_sha256": "97cfec470f84aa0aaba4191d1ce4cf0f65b97d86e7bf84a8ac2b95f124aad3e6",
            "validation_report_path": "artifacts/validation/track1-exomiser-2026-09-12.json",
            "validation_report_sha256": "a216298c07d6ccc78865d299258741e542bbc6736fefe8274a86c03e00ea8bb4",
        },
        "public-annotations": {
            "status": "materialized_and_reverified_for_acceptance",
        },
        "public-reference": {
            "status": "materialized_and_reverified_for_acceptance",
        },
    }
    for resource_id, fields in descriptor_fields.items():
        document = _load_artifact_json(
            verified_content,
            primary[f"resource:{resource_id}"],
        )
        _closed_descriptor(
            document,
            fields=fields,
            schema="genosapiens.track1-validation-resource-descriptor",
        )
        if (
            document.get("resource_id") != resource_id
            or document.get("kind") != resources[resource_id].kind
            or any(
                document.get(key) != value
                for key, value in descriptor_values[resource_id].items()
            )
            or type(document.get("limitation")) is not str
            or not document.get("limitation")
        ):
            _fail("validation resource descriptor differs from the plan")

    annotations = _load_artifact_json(
        verified_content,
        primary["resource:public-annotations"],
    )
    configuration_bindings = annotations["configuration_bindings"]
    installed_manifest_bindings = annotations["installed_manifest_bindings"]
    output_bindings = annotations["output_bindings"]
    unresolved_resource_ids = annotations["unresolved_resource_ids"]
    if (
        not isinstance(configuration_bindings, list)
        or not isinstance(installed_manifest_bindings, list)
        or not isinstance(output_bindings, list)
        or not isinstance(unresolved_resource_ids, list)
        or any(
            not isinstance(item, dict)
            or set(item) != {"resource_id", "path", "sha256"}
            for item in configuration_bindings
        )
        or any(
            not isinstance(item, dict)
            or set(item)
            != {
                "resource_id",
                "path",
                "sha256",
                "recorded_configuration_sha256",
                "matches_current_configuration",
            }
            for item in installed_manifest_bindings
        )
        or any(
            not isinstance(item, dict)
            or set(item)
            != {
                "resource_id",
                "manifest_sha256",
                "role",
                "path",
                "size_bytes",
                "sha256",
            }
            for item in output_bindings
        )
    ):
        _fail("validation annotation resource identity bindings are invalid")
    configurations = {
        item["resource_id"]: item for item in configuration_bindings
    }
    installed = {
        item["resource_id"]: item for item in installed_manifest_bindings
    }
    expected_configurations = {
        "ensembl-vep-116-grch38": (
            "config/public-resources/ensembl-vep-116-grch38.json",
            "f1c47e44cd07d8d3298f4c2b1c4b35392f57d1a5abebb03a4646276869254e1e",
        ),
        "rare-disease-evidence-2026-08-30-v6": (
            "config/public-resources/rare-disease-evidence-2026-08-30.json",
            "cded6e513755278d628a9e262336c66a3578708162b24f9356973fb8c65afa5f",
        ),
        "variant-class-resources-v1": (
            "config/public-resources/variant-class-resources.json",
            "cdea1165c93c81f7e1a1cecbe28defcc7eda2de61732f12e11c957838e3abccc",
        ),
    }
    expected_installed = {
        "broad-str-analysis-disease-catalog-grch38-d41c85c": (
            "manifests/variant-class-resources/broad-str-analysis-disease-catalog-grch38-d41c85c.json",
            "7058790cfa8b76bcdc6684f18cfd3152b70a51026412641955262c35e83a95c6",
            "cdea1165c93c81f7e1a1cecbe28defcc7eda2de61732f12e11c957838e3abccc",
        ),
        "encode-grch38-blacklist-encff356lfx": (
            "manifests/variant-class-resources/encode-grch38-blacklist-encff356lfx.json",
            "ae43f4ac6ea0623f31c97de1bed75435e0cb74a0c060c9f6ec59cb817a42cd06",
            "cdea1165c93c81f7e1a1cecbe28defcc7eda2de61732f12e11c957838e3abccc",
        ),
        "encode-re2g-registry-v4-grch38": (
            "manifests/variant-class-resources/encode-re2g-registry-v4-grch38.json",
            "11659dada45aa5930d227a8e476559e05eca3df8991a445374ac75b8138383a7",
            "cdea1165c93c81f7e1a1cecbe28defcc7eda2de61732f12e11c957838e3abccc",
        ),
        "encode-registry-v4-ccre-grch38": (
            "manifests/variant-class-resources/encode-registry-v4-ccre-grch38.json",
            "d2ca51d87593f52dc5186028b700d26122ae51a09f116046b0fadb5b107f4ea8",
            "cdea1165c93c81f7e1a1cecbe28defcc7eda2de61732f12e11c957838e3abccc",
        ),
        "ensembl-grch38-primary-assembly-fasta-116": (
            "manifests/ensembl/vep/ensembl-grch38-primary-assembly-fasta-116.f1c47e44cd07.json",
            "0b18876d150870a9013f9daef25c17da4006efcdd1a521037a1855355942a0a8",
            "f1c47e44cd07d8d3298f4c2b1c4b35392f57d1a5abebb03a4646276869254e1e",
        ),
        "ensembl-vep-116-grch38": (
            "manifests/ensembl/vep/ensembl-vep-116-grch38.f1c47e44cd07.bundle.json",
            "3071600d1ab2db89a9d0988465d137c0221378ebfdedb021b5392d4a8deb075a",
            "f1c47e44cd07d8d3298f4c2b1c4b35392f57d1a5abebb03a4646276869254e1e",
        ),
        "ensembl-vep-cache-116-grch38": (
            "manifests/ensembl/vep/ensembl-vep-cache-116-grch38.f1c47e44cd07.json",
            "4ee6f42b07a2a38aaaf30abeb3fe58b8acdf9a158df3bf8220212c49ad8ae824",
            "f1c47e44cd07d8d3298f4c2b1c4b35392f57d1a5abebb03a4646276869254e1e",
        ),
        "gnomad-exome-cnv-v4.1": (
            "manifests/variant-class-resources/gnomad-exome-cnv-v4.1.json",
            "b0122af8dc192324b2280f619d5184311da718c84f347331ed8717b8bbd5a821",
            "cdea1165c93c81f7e1a1cecbe28defcc7eda2de61732f12e11c957838e3abccc",
        ),
        "gnomad-mtdna-v3.1": (
            "manifests/variant-class-resources/gnomad-mtdna-v3.1.json",
            "6c8eb73f91b47911cc6518ad844bcbbdf7a3503a8c4b9f7340a931921d3a603e",
            "cdea1165c93c81f7e1a1cecbe28defcc7eda2de61732f12e11c957838e3abccc",
        ),
        "gnomad-sv-v4.1": (
            "manifests/variant-class-resources/gnomad-sv-v4.1.json",
            "91c7c82d4607433b7df57b666fc6ef55a2953dd7ccbd71aa2af44227d5592f84",
            "cdea1165c93c81f7e1a1cecbe28defcc7eda2de61732f12e11c957838e3abccc",
        ),
        "rare-disease-evidence-2026-08-30-v6": (
            "manifests/rare-disease-evidence/rare-disease-evidence-2026-08-30-v6.cded6e513755.bundle.json",
            "12c2bf613424137e3c8e7aec7f69ace0cf83a84a5a3bd0bcca81fc63970823d0",
            "cded6e513755278d628a9e262336c66a3578708162b24f9356973fb8c65afa5f",
        ),
        "variant-class-resources-v1": (
            "manifests/variant-class-resources/variant-class-resources-v1.bundle.json",
            "26ae62835bcda693a77623580d0383dd77768101653662a9eea6e77ed695c247",
            "cdea1165c93c81f7e1a1cecbe28defcc7eda2de61732f12e11c957838e3abccc",
        ),
    }
    if (
        len(configurations) != len(configuration_bindings)
        or tuple(configurations) != tuple(sorted(configurations))
        or {
            key: (item["path"], item["sha256"])
            for key, item in configurations.items()
        }
        != expected_configurations
        or len(installed) != len(installed_manifest_bindings)
        or tuple(installed) != tuple(sorted(installed))
        or {
            key: (
                item["path"],
                item["sha256"],
                item["recorded_configuration_sha256"],
            )
            for key, item in installed.items()
        }
        != expected_installed
        or any(
            item["matches_current_configuration"]
            is not (
                item["recorded_configuration_sha256"]
                == configurations[
                    "ensembl-vep-116-grch38"
                    if key.startswith("ensembl-")
                    else "variant-class-resources-v1"
                    if key in {
                        "broad-str-analysis-disease-catalog-grch38-d41c85c",
                        "encode-grch38-blacklist-encff356lfx",
                        "encode-re2g-registry-v4-grch38",
                        "encode-registry-v4-ccre-grch38",
                        "gnomad-exome-cnv-v4.1",
                        "gnomad-mtdna-v3.1",
                        "gnomad-sv-v4.1",
                    }
                    else key
                ]["sha256"]
            )
            for key, item in installed.items()
        )
        or unresolved_resource_ids
    ):
        _fail("validation annotation resource identities are inconsistent")
    expected_output_manifests = {
        ("ensembl-vep-116-grch38", "cache_info"): (
            "4ee6f42b07a2a38aaaf30abeb3fe58b8acdf9a158df3bf8220212c49ad8ae824"
        ),
        ("ensembl-vep-116-grch38", "cache_inventory"): (
            "4ee6f42b07a2a38aaaf30abeb3fe58b8acdf9a158df3bf8220212c49ad8ae824"
        ),
        ("ensembl-vep-116-grch38", "fai"): (
            "0b18876d150870a9013f9daef25c17da4006efcdd1a521037a1855355942a0a8"
        ),
        ("ensembl-vep-116-grch38", "fasta"): (
            "0b18876d150870a9013f9daef25c17da4006efcdd1a521037a1855355942a0a8"
        ),
    }
    expected_output_manifests.update(
        {
            ("rare-disease-evidence-2026-08-30-v6", role): (
                "12c2bf613424137e3c8e7aec7f69ace0cf83a84a5a3bd0bcca81fc63970823d0"
            )
            for role in {
                "disease",
                "disease_alias",
                "disease_inheritance",
                "disease_onset",
                "disease_prevalence",
                "disease_relation",
                "gene_assertion",
                "phenotype_assertion",
            }
        }
    )
    expected_output_manifests.update(
        {
            ("broad-str-analysis-disease-catalog-grch38-d41c85c", "json"):
                "7058790cfa8b76bcdc6684f18cfd3152b70a51026412641955262c35e83a95c6",
            ("encode-grch38-blacklist-encff356lfx", "bed"):
                "ae43f4ac6ea0623f31c97de1bed75435e0cb74a0c060c9f6ec59cb817a42cd06",
            ("encode-re2g-registry-v4-grch38", "tsv"):
                "11659dada45aa5930d227a8e476559e05eca3df8991a445374ac75b8138383a7",
            ("encode-registry-v4-ccre-grch38", "bed"):
                "d2ca51d87593f52dc5186028b700d26122ae51a09f116046b0fadb5b107f4ea8",
            ("gnomad-exome-cnv-v4.1", "vcf"):
                "b0122af8dc192324b2280f619d5184311da718c84f347331ed8717b8bbd5a821",
            ("gnomad-mtdna-v3.1", "tbi"):
                "6c8eb73f91b47911cc6518ad844bcbbdf7a3503a8c4b9f7340a931921d3a603e",
            ("gnomad-mtdna-v3.1", "vcf"):
                "6c8eb73f91b47911cc6518ad844bcbbdf7a3503a8c4b9f7340a931921d3a603e",
            ("gnomad-sv-v4.1", "tbi"):
                "91c7c82d4607433b7df57b666fc6ef55a2953dd7ccbd71aa2af44227d5592f84",
            ("gnomad-sv-v4.1", "vcf"):
                "91c7c82d4607433b7df57b666fc6ef55a2953dd7ccbd71aa2af44227d5592f84",
        }
    )
    output_keys = []
    for item in output_bindings:
        for value in (item["resource_id"], item["role"]):
            _identifier(value)
        for value in (item["manifest_sha256"], item["sha256"]):
            _sha(value)
        if (
            item["manifest_sha256"]
            != expected_output_manifests.get((item["resource_id"], item["role"]))
            or type(item["path"]) is not str
            or not item["path"]
            or type(item["size_bytes"]) is not int
            or item["size_bytes"] <= 0
        ):
            _fail("validation annotation output binding is invalid")
        output_keys.append((item["resource_id"], item["role"], item["path"]))
    if (
        output_keys != sorted(output_keys)
        or len(output_keys) != len(set(output_keys))
        or {
            (resource_id, role) for resource_id, role, _path in output_keys
        }
        != {
            ("ensembl-vep-116-grch38", role)
            for role in {"cache_info", "cache_inventory", "fai", "fasta"}
        }
        | {
            ("rare-disease-evidence-2026-08-30-v6", role)
            for role in {
                "disease",
                "disease_alias",
                "disease_inheritance",
                "disease_onset",
                "disease_prevalence",
                "disease_relation",
                "gene_assertion",
                "phenotype_assertion",
            }
        }
        | {
            ("broad-str-analysis-disease-catalog-grch38-d41c85c", "json"),
            ("encode-grch38-blacklist-encff356lfx", "bed"),
            ("encode-re2g-registry-v4-grch38", "tsv"),
            ("encode-registry-v4-ccre-grch38", "bed"),
            ("gnomad-exome-cnv-v4.1", "vcf"),
            ("gnomad-mtdna-v3.1", "tbi"),
            ("gnomad-mtdna-v3.1", "vcf"),
            ("gnomad-sv-v4.1", "tbi"),
            ("gnomad-sv-v4.1", "vcf"),
        }
    ):
        _fail("validation annotation output set is incomplete or unordered")

    reference = _load_artifact_json(
        verified_content,
        primary["resource:public-reference"],
    )
    configuration_binding = reference["configuration_binding"]
    installed_manifest_binding = reference["installed_manifest_binding"]
    reference_outputs = reference["output_bindings"]
    if (
        not isinstance(configuration_binding, dict)
        or set(configuration_binding) != {"resource_id", "path", "sha256"}
        or not isinstance(installed_manifest_binding, dict)
        or set(installed_manifest_binding) != {"resource_id", "path", "sha256"}
        or not isinstance(reference_outputs, list)
        or any(
            not isinstance(item, dict)
            or set(item)
            != {"role", "manifest_sha256", "path", "size_bytes", "sha256"}
            for item in reference_outputs
        )
        or configuration_binding
        != {
            "resource_id": "challenge-grch38-masked-v2-plus-hs38d1",
            "path": "config/reference-identities/challenge-grch38-masked-v2-hs38d1.json",
            "sha256": "a92c50f61c1b149e635656a54086e671f3f27f09b97002309ff68a9ded1d3b62",
        }
        or installed_manifest_binding
        != {
            "resource_id": "challenge-grch38-masked-v2-plus-hs38d1",
            "path": "manifests/grch38/challenge-grch38-masked-v2-plus-hs38d1.a92c50f61c1b.json",
            "sha256": REFERENCE_MANIFEST_SHA256,
        }
    ):
        _fail("validation reference resource identity bindings are invalid")
    reference_output_keys = []
    for item in reference_outputs:
        _identifier(item["role"])
        _sha(item["manifest_sha256"])
        _sha(item["sha256"])
        if (
            item["manifest_sha256"] != REFERENCE_MANIFEST_SHA256
            or type(item["path"]) is not str
            or not item["path"]
            or type(item["size_bytes"]) is not int
            or item["size_bytes"] <= 0
        ):
            _fail("validation reference output binding is invalid")
        reference_output_keys.append((item["role"], item["path"]))
    if (
        reference_output_keys != sorted(reference_output_keys)
        or len(reference_output_keys) != len(set(reference_output_keys))
        or {role for role, _path in reference_output_keys}
        != {"contig_aliases", "public_native_fai", "public_native_fasta"}
    ):
        _fail("validation reference output set is incomplete or unordered")

    evaluator = _load_artifact_json(
        verified_content,
        primary["resource:challenge-evaluator"],
    )
    if (
        set(evaluator)
        != {
            "schema", "task", "source", "inspection", "license_review",
            "inspection_metadata", "limits",
        }
        or evaluator.get("schema") != "genosapiens.challenge-evaluator-pin"
        or evaluator.get("task") != "build-track1-wgs-prioritization/11.1"
    ):
        _fail("validation challenge evaluator descriptor is invalid")

    giab = _load_artifact_json(verified_content, primary["resource:giab-truth"])
    if (
        set(giab)
        != {
            "schema", "resource_id", "release", "sample", "assembly",
            "reference_binding", "published_checksum_source",
            "autosomal_contigs", "artifacts",
        }
        or giab.get("schema") != "genosapiens.giab-hg002"
        or giab.get("resource_id") != "giab-hg002-nist-v4.2.1-grch38"
        or giab.get("sample") != "HG002"
        or giab.get("assembly") != "GRCh38"
    ):
        _fail("validation GIAB descriptor is invalid")

    truth = _load_artifact_json(
        verified_content,
        primary["resource:synthetic-causal-truth"],
    )
    _closed_descriptor(
        truth,
        fields={
            "schema", "classification", "validation_only", "audience",
            "cases", "limitations",
        },
        schema="genosapiens.track1-synthetic-causal-truth",
    )
    truth_cases = truth.get("cases")
    if (
        truth.get("audience") != "evaluator_only"
        or not isinstance(truth_cases, list)
        or len(truth_cases) != len(plan.cases)
        or any(not isinstance(item, dict) for item in truth_cases)
        or not isinstance(truth.get("limitations"), list)
        or not truth.get("limitations")
    ):
        _fail("validation causal-truth descriptor is invalid")
    plan_cases = {item.case_id: item for item in plan.cases}
    if {item.get("case_id") for item in truth_cases} != set(plan_cases):
        _fail("validation causal truth differs from the case roster")
    auxiliary = {}
    for artifact in bundle.artifacts:
        if artifact.parent_binding_id is not None:
            auxiliary.setdefault(artifact.parent_binding_id, []).append(artifact)
    plan_inputs = {item.input_id: item for item in plan.inputs}
    available_truth = {}
    for case in plan.cases:
        identities = set()
        for input_id in case.runner_input_ids:
            plan_input = plan_inputs[input_id]
            if plan_input.kind is InputKind.VCF:
                artifacts = [primary[f"input:{input_id}"]]
            elif plan_input.kind is InputKind.NATIVE:
                artifacts = [
                    artifact
                    for artifact in auxiliary.get(f"input:{input_id}", [])
                    if artifact.path.endswith(".vcf")
                ]
            else:
                artifacts = []
            for artifact in artifacts:
                identities.update(
                    _synthetic_vcf_truth_identities(
                        verified_content[artifact.path]
                    )
                )
        available_truth[case.case_id] = identities

    truth_members_by_split = {split: set() for split in Split}
    for item in truth_cases:
        hypotheses = item.get("hypotheses")
        if (
            set(item) != {"case_id", "hypotheses"}
            or not isinstance(hypotheses, list)
            or not hypotheses
            or any(
                not isinstance(hypothesis, dict)
                or set(hypothesis) != {"scope", "members"}
                or hypothesis.get("scope")
                not in plan_cases[item["case_id"]].assessed_scopes
                or not isinstance(hypothesis.get("members"), list)
                or not hypothesis.get("members")
                or len(hypothesis.get("members", []))
                != len(set(hypothesis.get("members", [])))
                or any(
                    type(member) is not str or not member
                    for member in hypothesis.get("members", [])
                )
                for hypothesis in hypotheses
            )
        ):
            _fail("validation causal truth is invalid")
        hypothesis_scopes = [hypothesis["scope"] for hypothesis in hypotheses]
        if (
            len(hypothesis_scopes) != len(set(hypothesis_scopes))
            or set(hypothesis_scopes)
            != set(plan_cases[item["case_id"]].assessed_scopes)
            or any(
                len(hypothesis["members"])
                != TRUTH_SCOPE_MEMBER_COUNTS[hypothesis["scope"]]
                for hypothesis in hypotheses
            )
        ):
            _fail("validation causal truth scope contract is invalid")
        signatures = [
            (hypothesis["scope"], tuple(hypothesis["members"]))
            for hypothesis in hypotheses
        ]
        if len(signatures) != len(set(signatures)):
            _fail("validation causal truth contains duplicate hypotheses")
        case_members = {
            member
            for hypothesis in hypotheses
            for member in hypothesis["members"]
        }
        if not case_members.issubset(available_truth[item["case_id"]]):
            _fail("validation causal truth is absent from its runner inputs")
        split = plan_cases[item["case_id"]].split
        truth_members_by_split[split].update(case_members)
    if truth_members_by_split[Split.DEVELOPMENT] & truth_members_by_split[Split.HELD_OUT]:
        _fail("validation causal identities cross development and held-out splits")
    background_identities = _synthetic_background_identities(bundle, verified_content)
    input_identities_by_split = {split: set() for split in Split}
    for case in plan.cases:
        input_identities_by_split[case.split].update(available_truth[case.case_id])
    if (
        (
            input_identities_by_split[Split.DEVELOPMENT]
            & input_identities_by_split[Split.HELD_OUT]
        )
        != background_identities
        or background_identities
        & (
            truth_members_by_split[Split.DEVELOPMENT]
            | truth_members_by_split[Split.HELD_OUT]
        )
    ):
        _fail("validation shared background identities are invalid")


def _validate_auxiliary_bindings(
    bundle: ValidationArtifactBundle,
    verified_content: dict[str, bytes],
) -> None:
    primary = {
        artifact.binding_id: artifact
        for artifact in bundle.artifacts
        if artifact.binding_id is not None
    }
    children: dict[str, list[ValidationArtifact]] = {}
    for artifact in bundle.artifacts:
        if artifact.parent_binding_id is not None:
            children.setdefault(artifact.parent_binding_id, []).append(artifact)

    supported_parents = {
        "generator:track1-synthetic-case-generator",
        "input:development-vcf",
        "input:held-out-native",
        "resource:public-annotations",
        "resource:public-reference",
        "resource:annotation-expectations",
    }
    if set(children) != supported_parents:
        _fail("validation auxiliary artifact parent set is incomplete")

    annotation_descriptor = _load_artifact_json(
        verified_content, primary["resource:annotation-expectations"],
    )
    materialization_children = children["resource:annotation-expectations"]
    materialization_bindings = {
        (annotation_descriptor[f"{role}_path"], annotation_descriptor[f"{role}_sha256"])
        for role in ("expectations", "roster", "validation_report")
    }
    if (len(materialization_children) != 3 or materialization_bindings != {
            (item.path, item.sha256) for item in materialization_children}):
        _fail("validation model-roster evidence differs from auxiliary bindings")
    materialization = _load_artifact_json(verified_content, next(
        item for item in materialization_children
        if item.path == annotation_descriptor["validation_report_path"]
    ))
    _closed_descriptor(
        materialization,
        schema="genosapiens.track1-model-roster-annotation-validation",
        fields={
            "schema", "classification", "llm_policy", "validation_only", "status",
            "controlled_sources_used", "controlled_use_authorized", "scope", "roster",
            "annotations", "original_roster", "original_cases_preserved", "additional_cases",
            "full_roster_cases", "request_count", "runtime_image_digest",
            "generating_executable_stack_sha256", "finalizer_sha256",
            "container_inspection_sha256", "annotation_binding_sha256",
            "input_evidence_sha256", "resource_bundle_sha256", "annotation_contract",
            "artifacts", "limitations",
        },
    )
    for key, role in (("roster", "roster"), ("annotations", "expectations")):
        bound = next(item for item in materialization_children
                     if item.path == annotation_descriptor[f"{role}_path"])
        if materialization[key] != {
                "path": bound.path, "size_bytes": bound.size_bytes, "sha256": bound.sha256}:
            _fail("validation model-roster annotation evidence is stale")
    counts = {"original_cases_preserved": 22, "additional_cases": 6,
              "full_roster_cases": 28, "request_count": 48}
    if (
        materialization["classification"] != "public"
        or materialization["llm_policy"] != "no-llm"
        or materialization["validation_only"] is not True
        or materialization["status"] != "passed"
        or materialization["controlled_sources_used"] is not False
        or materialization["controlled_use_authorized"] is not False
        or any(type(materialization[key]) is not int or materialization[key] != value
               for key, value in counts.items())
        or materialization["original_roster"].get("sha256")
        != "34b1ec2285c459e2e1b98fe9c8c53bf5efcbbe96f7eee58982a393b31ead03fa"
        or not materialization["limitations"]
    ):
        _fail("validation model-roster annotation evidence scope is invalid")
    for key in ("generating_executable_stack_sha256", "finalizer_sha256",
                "container_inspection_sha256", "annotation_binding_sha256",
                "input_evidence_sha256", "resource_bundle_sha256"):
        _sha(materialization[key])
    expectations = _load_artifact_json(verified_content, next(
        item for item in materialization_children
        if item.path == annotation_descriptor["expectations_path"]
    ))
    if (expectations["source_runs"]["current-model-state"]["binding_sha256"]
            != materialization["annotation_binding_sha256"]):
        _fail("validation model-roster annotation generating binding differs")

    generator_binding_id = "generator:track1-synthetic-case-generator"
    generator = _load_artifact_json(
        verified_content,
        primary[generator_binding_id],
    )
    _closed_descriptor(
        generator,
        fields={
            "schema",
            "classification",
            "validation_only",
            "generator_id",
            "version",
            "materialization",
            "controlled_sources_used",
            "executable_binding",
            "public_background_binding",
            "public_background_contract_binding",
            "generated_outputs",
            "determinism",
            "limitations",
        },
        schema="genosapiens.track1-synthetic-case-generator",
    )
    if (
        generator.get("generator_id") != "track1-synthetic-case-generator"
        or generator.get("version") != "v2"
        or generator.get("materialization") != "executable-seed-generator"
        or generator.get("controlled_sources_used") is not False
        or type(generator.get("determinism")) is not str
        or not generator.get("determinism")
        or not isinstance(generator.get("limitations"), list)
        or not generator.get("limitations")
    ):
        _fail("validation generator descriptor is invalid")
    expected_generator_references = {
        (
            generator["executable_binding"].get("path"),
            generator["executable_binding"].get("sha256"),
            generator["executable_binding"].get("size_bytes"),
        ),
        (
            generator["public_background_binding"].get("path"),
            generator["public_background_binding"].get("sha256"),
            generator["public_background_binding"].get("size_bytes"),
        ),
        (
            generator["public_background_contract_binding"].get("path"),
            generator["public_background_contract_binding"].get("sha256"),
            generator["public_background_contract_binding"].get("size_bytes"),
        ),
    }
    generator_children = {
        (artifact.path, artifact.sha256, artifact.size_bytes)
        for artifact in children[generator_binding_id]
    }
    binding_fields = {"path", "size_bytes", "sha256"}
    if (
        set(generator["executable_binding"]) != binding_fields
        or set(generator["public_background_binding"]) != binding_fields
        or set(generator["public_background_contract_binding"]) != binding_fields
        or len(expected_generator_references) != 3
        or expected_generator_references != generator_children
    ):
        _fail("validation generator input bindings are invalid")
    generated_outputs = generator["generated_outputs"]
    if (
        any(
            not isinstance(item, dict) or set(item) != binding_fields
            for item in generated_outputs
        )
        or len(generated_outputs) != len({item.get("path") for item in generated_outputs})
    ):
        _fail("validation generator output bindings are invalid")
    expected_outputs = {
        (
            artifact.path,
            artifact.sha256,
            artifact.size_bytes,
        )
        for artifact in bundle.artifacts
        if artifact.path
        in {
            "tests/fixtures/track1/materialized/causal-truth.manifest.json",
            "tests/fixtures/track1/materialized/development.hpo.txt",
            "tests/fixtures/track1/materialized/development.source-format.json",
            "tests/fixtures/track1/materialized/development.vcf",
            "tests/fixtures/track1/materialized/held-out-longtr.source-format.json",
            "tests/fixtures/track1/materialized/held-out-longtr.vcf",
            "tests/fixtures/track1/materialized/held-out-native.json",
            "tests/fixtures/track1/materialized/held-out-small-native.vcf",
            "tests/fixtures/track1/materialized/held-out-small.source-format.json",
            "tests/fixtures/track1/materialized/held-out.hpo.txt",
        }
    }
    declared_outputs = {
        (item.get("path"), item.get("sha256"), item.get("size_bytes"))
        for item in generated_outputs
    }
    if len(expected_outputs) != 10 or declared_outputs != expected_outputs:
        _fail("validation generator outputs differ from the artifact bundle")

    source_contract_artifacts = [
        artifact
        for artifact in children[generator_binding_id]
        if artifact.artifact_id == "generator-public-background-contract"
    ]
    if len(source_contract_artifacts) != 1:
        _fail("validation generator public background contract is missing")
    source_contract = _load_artifact_json(
        verified_content,
        source_contract_artifacts[0],
    )
    if (
        source_contract.get("schema") != "genosapiens.challenge-shaped-fixture"
        or source_contract.get("resource_id")
        != "challenge-shaped-public-small-variant-fixture-v2"
        or source_contract.get("release") != "v2"
        or source_contract_artifacts[0].sha256
        != PUBLIC_BACKGROUND_SOURCE["contract_sha256"]
    ):
        _fail("validation generator public background contract is invalid")

    def validate_source_format(
        artifact: ValidationArtifact,
        *,
        input_sha256: object,
        repeat_identity: str | None,
        mitochondrial_contigs: list[str],
    ) -> None:
        source_format = _load_artifact_json(verified_content, artifact)
        if (
            set(source_format) != SOURCE_FORMAT_FIELDS
            or source_format.get("schema") != "genosapiens.supplied-vcf-format"
            or source_format.get("input_sha256") != input_sha256
            or source_format.get("reference_role")
            != "public_annotation_benchmark"
            or source_format.get("reference_identity_id")
            != "grch38-no-alt-hs38d1-public-native-v1"
            or source_format.get("reference_manifest_sha256")
            != REFERENCE_MANIFEST_SHA256
            or source_format.get("repeat_identity") != repeat_identity
            or source_format.get("mitochondrial_contigs")
            != mitochondrial_contigs
        ):
            _fail("validation native source-format profile is inconsistent")

    development_source_formats = children["input:development-vcf"]
    if (
        len(development_source_formats) != 1
        or development_source_formats[0].artifact_id
        != "development-source-format"
    ):
        _fail("validation development source-format binding is invalid")
    validate_source_format(
        development_source_formats[0],
        input_sha256=primary["input:development-vcf"].sha256,
        repeat_identity="repid-ru-repcn-pos-end",
        mitochondrial_contigs=["M"],
    )

    native = _load_artifact_json(
        verified_content,
        primary["input:held-out-native"],
    )
    _closed_descriptor(
        native,
        fields={
            "schema", "classification", "validation_only", "case_id",
            "sample_id", "reference_role", "reference_identity_id",
            "reference_manifest_sha256", "members", "limitations",
        },
        schema="genosapiens.track1-synthetic-native-input",
    )
    members = native.get("members")
    if (
        native.get("case_id") != "held-out-case"
        or native.get("sample_id") != "proband"
        or native.get("reference_role") != "public_annotation_benchmark"
        or native.get("reference_identity_id")
        != "grch38-no-alt-hs38d1-public-native-v1"
        or native.get("reference_manifest_sha256")
        != REFERENCE_MANIFEST_SHA256
        or not isinstance(members, list)
        or any(not isinstance(item, dict) for item in members)
        or native.get("limitations")
        != [
            "This manifest preserves raw native members but is not yet a completed SuppliedSourceSetInputs replay artifact.",
            "Held-out execution remains blocked until development policy derivation and authenticated candidate-stack publication.",
        ]
    ):
        _fail("validation native input members are invalid")
    member_roles = {item.get("role") for item in members}
    if (
        member_roles != {"supplied-small-and-native-vcf", "supplied-longtr-vcf"}
        or len(members) != len(member_roles)
        or any(
            set(item)
            != {
                "role",
                "path",
                "size_bytes",
                "sha256",
                "sample_name",
                "source_format_binding",
            }
            or item.get("sample_name") != "SYNTHETIC"
            for item in members
        )
    ):
        _fail("validation native input member fields are invalid")
    native_references = {
        (item.get("path"), item.get("sha256"), item.get("size_bytes"))
        for item in members
    }
    source_format_bindings = [
        item.get("source_format_binding")
        for item in members
    ]
    if len(source_format_bindings) != 2 or any(
        not isinstance(binding, dict) for binding in source_format_bindings
    ):
        _fail("validation native source-format binding is invalid")
    for source_format_binding in source_format_bindings:
        if set(source_format_binding) != {
            "path",
            "size_bytes",
            "sha256",
            "repeat_identity",
        }:
            _fail("validation native source-format binding fields are invalid")
        native_references.add(
            (
                source_format_binding.get("path"),
                source_format_binding.get("sha256"),
                source_format_binding.get("size_bytes"),
            )
        )
    native_children = {
        (item.path, item.sha256, item.size_bytes)
        for item in children["input:held-out-native"]
    }
    if len(native_references) != len(members) * 2 or native_references != native_children:
        _fail("validation native input members differ from auxiliary bindings")

    expected_profiles = {
        "supplied-small-and-native-vcf": (None, ["M"]),
        "supplied-longtr-vcf": ("longtr-fd73e0d-vcf4.1/v1", []),
    }
    for member in members:
        source_format_binding = member["source_format_binding"]
        source_format_artifacts = [
            artifact
            for artifact in children["input:held-out-native"]
            if artifact.path == source_format_binding["path"]
        ]
        if (
            len(source_format_artifacts) != 1
            or source_format_binding["repeat_identity"]
            != expected_profiles[member["role"]][0]
        ):
            _fail("validation native source-format binding is invalid")
        validate_source_format(
            source_format_artifacts[0],
            input_sha256=member["sha256"],
            repeat_identity=expected_profiles[member["role"]][0],
            mitochondrial_contigs=expected_profiles[member["role"]][1],
        )

    annotations = _load_artifact_json(
        verified_content,
        primary["resource:public-annotations"],
    )
    bindings = annotations.get("configuration_bindings")
    if not isinstance(bindings, list) or any(
        not isinstance(item, dict)
        or set(item) != {"resource_id", "path", "sha256"}
        for item in bindings
    ):
        _fail("validation annotation-resource bindings are invalid")
    annotation_references = {
        (item.get("path"), item.get("sha256")) for item in bindings
    }
    annotation_children = {
        (item.path, item.sha256)
        for item in children["resource:public-annotations"]
    }
    if (
        len(bindings) != len(annotation_references)
        or annotation_references != annotation_children
    ):
        _fail("validation annotation resources differ from auxiliary bindings")
    annotation_configs = {
        item["resource_id"]: _load_artifact_json(
            verified_content,
            next(
                artifact
                for artifact in children["resource:public-annotations"]
                if artifact.path == item["path"]
            ),
        )
        for item in bindings
    }
    vep = annotation_configs["ensembl-vep-116-grch38"]
    rare_disease = annotation_configs["rare-disease-evidence-2026-08-30-v6"]
    variant_class = annotation_configs["variant-class-resources-v1"]
    expected_annotation_outputs = {
        ("ensembl-vep-116-grch38", "cache_info"): vep["cache"]["info_path"],
        ("ensembl-vep-116-grch38", "cache_inventory"): (
            vep["cache"]["output_root"] + "/cache.inventory.json"
        ),
        ("ensembl-vep-116-grch38", "fai"): vep["fasta"]["output_path"]
        + ".fai",
        ("ensembl-vep-116-grch38", "fasta"): vep["fasta"]["output_path"],
    }
    expected_annotation_outputs.update(
        {
            ("rare-disease-evidence-2026-08-30-v6", table): (
                rare_disease["projection"]["output_prefix"] + f"/{table}.parquet"
            )
            for table in rare_disease["projection"]["expected_row_counts"]
        }
    )
    expected_annotation_outputs.update(
        {
            (resource["resource_id"], artifact["role"]): artifact["path"]
            for resource in variant_class["resources"]
            for artifact in resource["artifacts"]
        }
    )
    actual_annotation_outputs = {
        (item["resource_id"], item["role"]): item["path"]
        for item in annotations["output_bindings"]
    }
    vep_outputs = {
        item["role"]: item
        for item in annotations["output_bindings"]
        if item["resource_id"] == "ensembl-vep-116-grch38"
    }
    if (
        vep.get("bundle_id") != "ensembl-vep-116-grch38"
        or rare_disease.get("bundle_id")
        != "rare-disease-evidence-2026-08-30-v6"
        or variant_class.get("bundle_id") != "variant-class-resources-v1"
        or actual_annotation_outputs != expected_annotation_outputs
        or annotations["unresolved_resource_ids"]
        or (
            vep_outputs["cache_info"]["size_bytes"],
            vep_outputs["cache_info"]["sha256"],
        )
        != (3365, vep["cache"]["info_sha256"])
        or (
            vep_outputs["fasta"]["size_bytes"],
            vep_outputs["fasta"]["sha256"],
        )
        != (vep["fasta"]["output_size_bytes"], vep["fasta"]["output_sha256"])
        or (
            vep_outputs["fai"]["size_bytes"],
            vep_outputs["fai"]["sha256"],
        )
        != (vep["fasta"]["fai_size_bytes"], vep["fasta"]["fai_sha256"])
    ):
        _fail("validation annotation identities differ from current configurations")

    reference = _load_artifact_json(
        verified_content,
        primary["resource:public-reference"],
    )
    binding = reference.get("configuration_binding")
    if (
        not isinstance(binding, dict)
        or set(binding) != {"resource_id", "path", "sha256"}
    ):
        _fail("validation reference-resource binding is invalid")
    reference_references = {(binding.get("path"), binding.get("sha256"))}
    reference_children = {
        (item.path, item.sha256)
        for item in children["resource:public-reference"]
    }
    if (
        len(children["resource:public-reference"]) != 1
        or reference_references != reference_children
        or native.get("reference_manifest_sha256")
        != reference["installed_manifest_binding"]["sha256"]
    ):
        _fail("validation reference resource differs from its auxiliary binding")
    reference_config = _load_artifact_json(
        verified_content,
        children["resource:public-reference"][0],
    )
    expected_reference_outputs = {
        item["role"]: item["path"]
        for item in reference_config["expected_outputs"]
        if item["role"]
        in {"contig_aliases", "public_native_fai", "public_native_fasta"}
    }
    if (
        reference_config.get("resource_id") != binding["resource_id"]
        or reference["installed_manifest_binding"]["resource_id"]
        != binding["resource_id"]
        or {
            item["role"]: item["path"] for item in reference["output_bindings"]
        }
        != expected_reference_outputs
    ):
        _fail("validation reference identities differ from current configuration")


def _validate(plan: ValidationPlan) -> None:
    if type(plan) is not ValidationPlan or type(plan.lifecycle) is not Lifecycle:
        _fail()
    if (
        plan.schema != PLAN_SCHEMA
        or plan.classification != "public"
        or plan.validation_only is not True
        or plan.controlled_use_authorized is not False
    ):
        _fail()
    _identifier(plan.plan_id)
    from genosapiens.track1_roster import SUBMISSION_PROFILE_SHA256, PROFILE_ACCEPTANCE_SHA256
    from genosapiens.track1_analysis_profiles import POLICY_SHA256 as PROFILE_POLICY_SHA256
    if (plan.analysis_profile not in {'focused', 'exhaustive'}
            or plan.analysis_profile_policy_sha256 != PROFILE_POLICY_SHA256
            or plan.profile_acceptance_sha256 != PROFILE_ACCEPTANCE_SHA256):
        _fail('validation plan does not bind the declared analysis profile acceptance')
    if (plan.submission_acceptance_profile_sha256 != SUBMISSION_PROFILE_SHA256
            or tuple((r.resource_id, r.descriptor_sha256) for r in plan.resources
                if r.kind == "submission_acceptance_profile")
            != (("submission-acceptance-profile", SUBMISSION_PROFILE_SHA256),)):
        _fail("validation plan does not bind the declared submission acceptance profile")
    if plan.required_scopes != REQUIRED_SCOPES:
        _fail("validation plan required scopes are incomplete or unordered")
    if plan.required_analysis_modes != REQUIRED_ANALYSIS_MODES:
        _fail("validation plan required analysis modes are incomplete or unordered")
    if plan.assessed_subset_accounting != ASSESSED_SUBSET_ACCOUNTING:
        _fail("validation assessed-subset accounting is not fixed")
    if not plan.limitations or any(type(item) is not str or not item for item in plan.limitations):
        _fail("validation plan limitations are required")
    if any(
        type(value) is not tuple
        for value in (
            plan.required_scopes,
            plan.required_analysis_modes,
            plan.generators,
            plan.inputs,
            plan.resources,
            plan.backgrounds,
            plan.oracles,
            plan.cases,
            plan.metrics,
            plan.limitations,
            plan.lifecycle.acceptance_criteria,
        )
    ):
        _fail()
    if any(type(item) is not GeneratorBinding for item in plan.generators):
        _fail()
    if any(
        type(item) is not InputBinding or type(item.kind) is not InputKind
        for item in plan.inputs
    ):
        _fail()
    if any(type(item) is not ResourceBinding for item in plan.resources):
        _fail()
    if any(type(item) is not BackgroundBinding for item in plan.backgrounds):
        _fail()
    if any(
        type(item) is not OracleBinding or type(item.role) is not OracleRole
        for item in plan.oracles
    ):
        _fail()
    if any(
        type(item) is not ValidationCase or type(item.split) is not Split
        for item in plan.cases
    ):
        _fail()
    if any(type(item) is not MetricDefinition for item in plan.metrics):
        _fail()
    if any(
        type(item) is not AcceptanceCriterion
        for item in plan.lifecycle.acceptance_criteria
    ):
        _fail()

    collections = (
        (plan.generators, lambda item: item.generator_id, "generators"),
        (plan.inputs, lambda item: item.input_id, "inputs"),
        (plan.resources, lambda item: item.resource_id, "resources"),
        (plan.backgrounds, lambda item: item.background_id, "backgrounds"),
        (plan.oracles, lambda item: item.oracle_id, "oracles"),
        (plan.cases, lambda item: item.case_id, "cases"),
        (plan.metrics, lambda item: item.metric_id, "metrics"),
    )
    for values, key, what in collections:
        _ordered_unique(values, key=key, what=what)

    generators = {item.generator_id: item for item in plan.generators}
    inputs = {item.input_id: item for item in plan.inputs}
    resources = {item.resource_id: item for item in plan.resources}
    backgrounds = {item.background_id: item for item in plan.backgrounds}
    oracles = {item.oracle_id: item for item in plan.oracles}
    _identifier(plan.evaluator_resource_id)
    evaluator = resources.get(plan.evaluator_resource_id)
    if evaluator is None or evaluator.kind != "challenge_evaluator":
        _fail("validation plan challenge evaluator is not bound")
    for item in plan.generators:
        _identifier(item.generator_id)
        _sha(item.manifest_sha256)
    for item in plan.resources:
        _identifier(item.resource_id)
        _sha(item.descriptor_sha256)
    for item in plan.generators:
        _identifier(item.version)
    if any(
        type(item.kind) is not str or item.kind not in RESOURCE_KINDS
        for item in plan.resources
    ):
        _fail("validation resource kind is unsupported")
    resource_hashes = tuple(item.descriptor_sha256 for item in plan.resources)
    if len(resource_hashes) != len(set(resource_hashes)):
        _fail("validation resources cannot alias the same manifest bytes")
    for item in plan.inputs:
        _identifier(item.input_id)
        _sha(item.sha256)
        if item.generator_id not in generators or type(item.seed) is not int or item.seed < 0:
            _fail("validation input has an invalid generator or seed")
    for item in plan.backgrounds:
        _identifier(item.background_id)
        _identifier(item.subject_id)
        if (item.kind == "HG002") is not (item.subject_id == "HG002"):
            _fail("HG002 background subject identity must remain canonical")
        if (
            type(item.kind) is not str
            or item.kind not in {"HG002", "independent_synthetic"}
            or type(item.shared_across_cases) is not bool
            or type(item.disclosure) is not str
            or not item.disclosure
        ):
            _fail("validation background disclosure is required")
    for item in plan.oracles:
        _identifier(item.oracle_id)
        if (
            item.resource_id not in resources
            or item.source_kind != ORACLE_SOURCES[item.role]
            or resources[item.resource_id].kind != ORACLE_RESOURCE_KINDS[item.role]
        ):
            _fail("validation oracle role and source are inconsistent")
    if len(plan.oracles) != len(OracleRole) or {item.role for item in plan.oracles} != set(OracleRole):
        _fail("validation plan must keep all oracle roles separate")

    oracle_resource_ids = {item.resource_id for item in plan.oracles}
    group_splits: dict[tuple[str, str], set[Split]] = {}
    group_hashes: dict[tuple[str, str], set[str]] = {}
    group_manifest_splits: dict[str, set[Split]] = {}
    input_hash_splits: dict[str, set[Split]] = {}
    generator_seed_splits: dict[tuple[str, int], set[Split]] = {}
    split_scopes = {split: set() for split in Split}
    background_uses: dict[str, int] = {}
    background_subject_uses: dict[str, int] = {}
    for case in plan.cases:
        for value in (
            case.case_id,
            case.recipe_group_id,
            case.perturbation_group_id,
            case.generator_id,
            case.background_id,
        ):
            _identifier(value)
        _sha(case.recipe_manifest_sha256)
        _sha(case.perturbation_manifest_sha256)
        if (
            case.generator_id not in generators
            or case.background_id not in backgrounds
            or type(case.seed) is not int
            or case.seed < 0
            or not case.runner_input_ids
            or case.runner_input_ids != tuple(sorted(set(case.runner_input_ids)))
            or not case.runner_resource_ids
            or case.runner_resource_ids != tuple(sorted(set(case.runner_resource_ids)))
            or case.assessed_scopes != tuple(sorted(set(case.assessed_scopes)))
            or not case.assessed_scopes
            or set(case.assessed_scopes) - set(REQUIRED_SCOPES)
        ):
            _fail("validation case bindings are incomplete or unordered")
        selected_inputs = [inputs.get(key) for key in case.runner_input_ids]
        if (
            any(item is None for item in selected_inputs)
            or len({item.sha256 for item in selected_inputs}) != len(selected_inputs)
            or any(
                item.generator_id != case.generator_id or item.seed != case.seed
                for item in selected_inputs
            )
            or {item.kind for item in selected_inputs}.isdisjoint({InputKind.VCF, InputKind.NATIVE})
            or InputKind.HPO_TEXT not in {item.kind for item in selected_inputs}
        ):
            _fail("validation case inputs do not bind one generator/seed and required input kinds")
        if (
            any(key not in resources for key in case.runner_resource_ids)
            or set(case.runner_resource_ids) & oracle_resource_ids
            or any(resources[key].kind not in RUNNER_RESOURCE_KINDS for key in case.runner_resource_ids)
        ):
            _fail("ranker resource view contains missing or evaluator-only oracle evidence")
        causal = oracles.get(case.causal_truth_oracle_id)
        if causal is None or causal.role is not OracleRole.CAUSAL_TRUTH:
            _fail("validation case causal truth is not evaluator-only")
        for kind, group in (
            ("recipe", case.recipe_group_id),
            ("perturbation", case.perturbation_group_id),
        ):
            group_splits.setdefault((f"{kind}_id", group), set()).add(case.split)
        for kind, group, manifest_sha256 in (
            ("recipe", case.recipe_group_id, case.recipe_manifest_sha256),
            (
                "perturbation",
                case.perturbation_group_id,
                case.perturbation_manifest_sha256,
            ),
        ):
            group_hashes.setdefault((kind, group), set()).add(manifest_sha256)
            group_manifest_splits.setdefault(manifest_sha256, set()).add(case.split)
            group_splits.setdefault((f"{kind}_manifest", manifest_sha256), set()).add(
                case.split
            )
        for item in selected_inputs:
            input_hash_splits.setdefault(item.sha256, set()).add(case.split)
        generator_seed_splits.setdefault((case.generator_id, case.seed), set()).add(
            case.split
        )
        split_scopes[case.split].update(case.assessed_scopes)
        background_uses[case.background_id] = background_uses.get(case.background_id, 0) + 1
        subject_id = backgrounds[case.background_id].subject_id
        background_subject_uses[subject_id] = background_subject_uses.get(subject_id, 0) + 1
    if any(len(splits) != 1 for splits in group_splits.values()):
        _fail("related validation recipes or perturbations cross splits")
    if any(len(hashes) != 1 for hashes in group_hashes.values()):
        _fail("validation recipe or perturbation identity changed manifest")
    if any(len(splits) != 1 for splits in group_manifest_splits.values()):
        _fail("related validation recipe or perturbation manifests cross splits")
    if any(len(splits) != 1 for splits in input_hash_splits.values()):
        _fail("validation input bytes cross development and held-out splits")
    if any(len(splits) != 1 for splits in generator_seed_splits.values()):
        _fail("validation generator seed crosses development and held-out splits")
    if any(scopes != set(REQUIRED_SCOPES) for scopes in split_scopes.values()):
        _fail("each validation split must cover every required assessed scope")
    for background_id in background_uses:
        background = backgrounds[background_id]
        subject_uses = background_subject_uses[background.subject_id]
        if background.shared_across_cases is not (subject_uses > 1):
            _fail("validation background sharing declaration is inconsistent")
        if subject_uses > 1 and background.kind == "HG002" and not background.disclosure:
            _fail("shared HG002 background is not disclosed")
    used_inputs = {key for case in plan.cases for key in case.runner_input_ids}
    used_resources = (
        {key for case in plan.cases for key in case.runner_resource_ids}
        | oracle_resource_ids
        | {plan.evaluator_resource_id}
        | {"submission-acceptance-profile"}
    )
    used_generators = {case.generator_id for case in plan.cases}
    if (
        used_inputs != set(inputs)
        or used_resources != set(resources)
        or used_generators != set(generators)
        or set(background_uses) != set(backgrounds)
    ):
        _fail("validation plan contains unbound inputs, resources, generators, or backgrounds")

    metrics = {item.metric_id: item for item in plan.metrics}
    if set(metrics) != set(METRIC_CONTRACTS):
        _fail("validation metric set is incomplete")
    for metric_id, expected in METRIC_CONTRACTS.items():
        metric = metrics[metric_id]
        expected_non_success = NON_SUCCESS_ACCOUNTING
        if metric_id in STATUS_ACCOUNTING_METRICS:
            expected_non_success = STATUS_ACCOUNTING
        elif metric_id in NULL_MEASUREMENT_METRICS:
            expected_non_success = NULL_MEASUREMENT_ACCOUNTING
        if (
            (
                metric.unit,
                metric.numerator,
                metric.denominator,
                metric.reporting_dimensions,
            )
            != expected
            or metric.non_success_statuses != NON_SUCCESS_STATUSES
            or metric.non_success_accounting != expected_non_success
        ):
            _fail("validation metric accounting differs from the fixed contract")

    lifecycle = plan.lifecycle
    criteria = lifecycle.acceptance_criteria
    criterion_ids = tuple(item.metric_id for item in criteria)
    if criterion_ids != tuple(sorted(criterion_ids)) or len(criterion_ids) != len(set(criterion_ids)):
        _fail("validation acceptance criteria must be unique and ordered")
    for criterion in criteria:
        if (
            criterion.metric_id not in metrics
            or criterion.comparator not in {">=", "<="}
            or DECIMAL.fullmatch(criterion.threshold) is None
        ):
            _fail("validation acceptance criterion is invalid")
        _sha(criterion.derivation_sha256)
    for value in (
        lifecycle.development_derivation_sha256,
        lifecycle.candidate_stack_sha256,
        lifecycle.post_held_out_policy_change_sha256,
    ):
        if value is not None:
            _sha(value)
    valid = (
        lifecycle.state == "development_plan"
        and lifecycle.held_out_execution == "blocked"
        and lifecycle.independence_status == "not_evaluated"
        and lifecycle.development_derivation_sha256 is None
        and lifecycle.candidate_stack_sha256 is None
        and lifecycle.post_held_out_policy_change_sha256 is None
        and not criteria
    )
    if not valid:
        _fail("static validation-plan contract cannot authorize held-out execution")


def load_validation_plan(
    plan_path: str | Path,
    *,
    repository_root: str | Path,
    expected_sha256: str,
) -> ValidationPlan:
    """Load exact hash-bound public plan bytes and validate cross-field isolation."""
    try:
        _sha(expected_sha256)
        root = Path(repository_root).resolve(strict=True)
        schema_path = _resolve_repository_file(root, CANONICAL_SCHEMA_PATH)
        source = Path(plan_path)
        try:
            source.relative_to(root)
            source_root: Path | None = root
        except ValueError:
            source_root = None
        document, content = load_bounded_json(
            source,
            reference_root=source_root,
        )
        schema, schema_content = load_bounded_json(
            schema_path,
            reference_root=root,
        )
        if hashlib.sha256(content).hexdigest() != expected_sha256:
            _fail("validation plan bytes do not match their independent binding")
        if hashlib.sha256(schema_content).hexdigest() != CANONICAL_SCHEMA_SHA256:
            _fail("canonical validation-plan schema changed")
        if not isinstance(document, dict) or not isinstance(schema, dict):
            _fail()
        validate_json_schema(document, schema)
        plan = _parse(document)
        _validate(plan)
        return plan
    except ValidationPlanError:
        raise
    except (BoundedJSONError, JsonSchemaError, KeyError, OSError, TypeError, ValueError):
        _fail()


def _resolve_repository_file(root: Path, relative_value: str) -> Path:
    if type(relative_value) is not str:
        _fail("validation-plan pin path is invalid")
    relative = PurePosixPath(relative_value)
    if (
        relative.is_absolute()
        or relative == PurePosixPath(".")
        or ".." in relative.parts
        or "\\" in relative_value
        or relative.as_posix() != relative_value
    ):
        _fail("validation-plan pin path is invalid")
    return root.joinpath(*relative.parts)


def _load_validation_plan_pin(root: Path) -> ValidationPlanPin:
    pin_path = _resolve_repository_file(root, CANONICAL_PIN_PATH)
    document, pin_content = load_bounded_json(pin_path, reference_root=root)
    if hashlib.sha256(pin_content).hexdigest() != CANONICAL_PIN_SHA256:
        _fail("canonical validation-plan pin changed")
    required = {
        "schema",
        "plan_id",
        "canonical_plan_path",
        "canonical_plan_sha256",
        "canonical_schema_path",
        "canonical_schema_sha256",
        "canonical_artifact_bundle_path",
        "canonical_artifact_bundle_sha256",
        "canonical_artifact_bundle_schema_path",
        "canonical_artifact_bundle_schema_sha256",
    }
    if not isinstance(document, dict) or set(document) != required:
        _fail("validation-plan pin fields are invalid")
    pin = ValidationPlanPin(**document)
    if pin.schema != PIN_SCHEMA:
        _fail("validation-plan pin schema is unsupported")
    _identifier(pin.plan_id)
    for value in (
        pin.canonical_plan_sha256,
        pin.canonical_schema_sha256,
        pin.canonical_artifact_bundle_sha256,
        pin.canonical_artifact_bundle_schema_sha256,
    ):
        _sha(value)
    if (
        pin.canonical_schema_path != CANONICAL_SCHEMA_PATH
        or pin.canonical_schema_sha256 != CANONICAL_SCHEMA_SHA256
        or pin.canonical_artifact_bundle_schema_path
        != CANONICAL_ARTIFACT_BUNDLE_SCHEMA_PATH
        or pin.canonical_artifact_bundle_schema_sha256
        != CANONICAL_ARTIFACT_BUNDLE_SCHEMA_SHA256
    ):
        _fail("validation-plan pin does not bind the canonical schemas")
    for relative_path, expected_sha256 in (
        (pin.canonical_plan_path, pin.canonical_plan_sha256),
        (pin.canonical_schema_path, pin.canonical_schema_sha256),
        (
            pin.canonical_artifact_bundle_path,
            pin.canonical_artifact_bundle_sha256,
        ),
        (
            pin.canonical_artifact_bundle_schema_path,
            pin.canonical_artifact_bundle_schema_sha256,
        ),
    ):
        content = read_bounded_bytes(
            _resolve_repository_file(root, relative_path),
            maximum_bytes=MAX_VALIDATION_ARTIFACT_BYTES,
            reference_root=root,
        )
        if hashlib.sha256(content).hexdigest() != expected_sha256:
            _fail("validation-plan pin target changed")
    return pin


def load_validation_plan_from_pin(
    *,
    repository_root: str | Path,
) -> ValidationPlan:
    """Load a plan through its committed, closed repository-relative pin."""
    try:
        root = Path(repository_root).resolve(strict=True)
        pin = _load_validation_plan_pin(root)
        plan = load_validation_plan(
            _resolve_repository_file(root, pin.canonical_plan_path),
            repository_root=root,
            expected_sha256=pin.canonical_plan_sha256,
        )
        if plan.plan_id != pin.plan_id:
            _fail("validation-plan pin identifies another plan")
        return plan
    except ValidationPlanError:
        raise
    except (BoundedJSONError, OSError, TypeError, ValueError):
        _fail("validation-plan pin is invalid")


def load_validation_artifact_bundle_from_pin(
    *,
    repository_root: str | Path,
) -> ValidationArtifactBundle:
    """Verify every canonical plan artifact for evaluator-side setup."""
    try:
        root = Path(repository_root).resolve(strict=True)
        pin = _load_validation_plan_pin(root)
        plan = load_validation_plan(
            _resolve_repository_file(root, pin.canonical_plan_path),
            repository_root=root,
            expected_sha256=pin.canonical_plan_sha256,
        )
        schema_path = _resolve_repository_file(
            root,
            pin.canonical_artifact_bundle_schema_path,
        )
        bundle_path = _resolve_repository_file(
            root,
            pin.canonical_artifact_bundle_path,
        )
        schema, schema_content = load_bounded_json(
            schema_path,
            reference_root=root,
        )
        document, bundle_content = load_bounded_json(
            bundle_path,
            reference_root=root,
        )
        if (
            hashlib.sha256(schema_content).hexdigest()
            != pin.canonical_artifact_bundle_schema_sha256
            or hashlib.sha256(bundle_content).hexdigest()
            != pin.canonical_artifact_bundle_sha256
            or not isinstance(schema, dict)
            or not isinstance(document, dict)
        ):
            _fail("validation artifact bundle changed")
        validate_json_schema(document, schema)
        bundle = _parse_artifact_bundle(document)
        _validate_artifact_bundle(
            bundle,
            plan=plan,
            plan_sha256=pin.canonical_plan_sha256,
        )
        verified_content = {}
        for artifact in bundle.artifacts:
            content = read_bounded_bytes(
                root / PurePosixPath(artifact.path),
                maximum_bytes=MAX_VALIDATION_ARTIFACT_BYTES,
                reference_root=root,
                expected_size_bytes=artifact.size_bytes,
            )
            if hashlib.sha256(content).hexdigest() != artifact.sha256:
                _fail("validation artifact bytes differ from their bundle binding")
            verified_content[artifact.path] = content
        _validate_primary_descriptors(plan, bundle, verified_content)
        _validate_auxiliary_bindings(bundle, verified_content)
        return bundle
    except ValidationPlanError:
        raise
    except (BoundedJSONError, JsonSchemaError, KeyError, OSError, TypeError, ValueError):
        _fail("validation artifact bundle is invalid")


def create_ranker_case_view(plan: ValidationPlan, case_id: str) -> RankerCaseView:
    """Return only runner-safe input/resource bindings, never evaluator truth."""
    _validate(plan)
    matches = tuple(item for item in plan.cases if item.case_id == case_id)
    if len(matches) != 1:
        _fail("validation ranker case is unavailable or ambiguous")
    case = matches[0]
    if case.split is Split.HELD_OUT:
        _fail("held-out ranker view is blocked by the static development plan")
    inputs = {item.input_id: item for item in plan.inputs}
    resources = {item.resource_id: item for item in plan.resources}
    return RankerCaseView(
        inputs=tuple(
            RankerInput(item.kind, item.sha256)
            for item in (inputs[key] for key in case.runner_input_ids)
        ),
        resources=tuple(
            RankerResource(item.kind, item.descriptor_sha256)
            for item in (resources[key] for key in case.runner_resource_ids)
        ),
    )
