"""Closed immutable Track 1 base release freeze and consumer loader.

The release is the first record that can authorize controlled VCF-first use.
It binds an already authenticated candidate plus consolidated accepted public
results; neither ancestor contains the release hash, keeping the graph acyclic.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass
import hashlib
import json
from pathlib import Path
import re
from typing import TYPE_CHECKING, Literal

from genosapiens.track1_policy_bindings import (
    PolicyBindings,
    decode_policy_bindings,
)

if TYPE_CHECKING:
    from genosapiens.track1_validation_runner import (
        ConsolidatedAcceptance,
        PreHeldOutCandidate,
    )


SCHEMA = "genosapiens.track1-base-release"
SHA256 = re.compile(r"^[0-9a-f]{64}$")
IDENTIFIER = re.compile(r"^[A-Za-z0-9][A-Za-z0-9._:-]*$")

EVIDENCE_INTERFACE_FILES = (
    'evidence_model.py',
    'evidence_serialization.py',
    'track1_artifacts.py',
    'track1_analysis_profiles.py',
    'wgs_profile_projection.py',
    'regulatory_index.py',
    'wgs_admission_projection.py',
    'wgs_admission_inputs.py',
    'wgs_association_codec.py',
    'wgs_association_projection.py',
    'wgs_analysis_inputs.py',
    'wgs_analysis_policy.py',
    'wgs_analysis_eligibility.py',
    'wgs_annotation_prefix.py',
    'wgs_annotation.py',
    'wgs_case_context.py',
    'wgs_case_evidence.py',
    'wgs_case_identities.py',
    'wgs_case_order.py',
    'wgs_case_resources.py',
    'wgs_case_science_facts.py',
    'wgs_completed_science.py',
    'wgs_context_evidence.py',
    'wgs_context_order.py',
    'wgs_context_population.py',
    'wgs_context_scores.py',
    'wgs_diagnosis_derivative.py',
    'wgs_eligibility_evidence.py',
    'wgs_evidence_replay.py',
    'wgs_evidence_store.py',
    'wgs_failures.py',
    'wgs_gene_evidence.py',
    'wgs_genotype_evaluation.py',
    'wgs_genotype_facts.py',
    'wgs_handoff.py',
    'wgs_hypothesis_evidence.py',
    'wgs_hypothesis_verification.py',
    'wgs_ic_evidence.py',
    'wgs_ledger.py',
    'wgs_locator_tape.py',
    'wgs_logical_evidence.py',
    'wgs_packed_outcomes.py',
    'wgs_pair_decisions.py',
    'wgs_pair_domain.py',
    'wgs_pair_evidence.py',
    'wgs_pair_outcomes.py',
    'wgs_pair_scores.py',
    'wgs_pair_verification.py',
    'wgs_phase_evaluation.py',
    'wgs_phase_patterns.py',
    'wgs_phenotype_evidence.py',
    'wgs_population_evidence.py',
    'wgs_population_input.py',
    'wgs_population_work_filter.py',
    'population_filter_candidate.py',
    'wgs_work_contract.py',
    'wgs_prepared_case.py',
    'wgs_provenance.py',
    'wgs_public_case.py',
    'wgs_public_execution.py',
    'wgs_public_replay.py',
    'wgs_read_rescue.py',
    'wgs_reference_evidence.py',
    'wgs_reporting.py',
    'wgs_representation_evidence.py',
    'wgs_resume.py',
    'wgs_scientific_artifact.py',
    'wgs_scientific_products.py',
    'wgs_scientific_sensitivities.py',
    'wgs_score_decisions.py',
    'wgs_score_sources.py',
    'wgs_score_views.py',
    'wgs_selection_evidence.py',
    'wgs_selection_inputs.py',
    'wgs_sensitivity_evidence.py',
    'wgs_shards.py',
    'wgs_source_closure.py',
    'wgs_validation_audit.py',
    'wgs_validation_authority.py',
    'wgs_validation_projection.py',
    'wgs_validation_recovery.py',
)


def evidence_interface_bytes(repository_root: str | Path) -> bytes:
    """Bind the single current ledger/provenance/intake/handoff interface."""
    root = Path(repository_root) / "src/genosapiens"
    return _bytes({name: hashlib.sha256((root / name).read_bytes()).hexdigest()
                   for name in EVIDENCE_INTERFACE_FILES})


class Track1ReleaseError(ValueError):
    """Static release diagnostics without controlled payloads."""


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


def _sha(value: object) -> None:
    if type(value) is not str or SHA256.fullmatch(value) is None:
        raise Track1ReleaseError("invalid release binding")


@dataclass(frozen=True, slots=True)
class BaseRelease:
    release_id: str
    authorization_status: Literal["candidate_only", "accepted"]
    candidate_sha256: str
    consolidated_acceptance_sha256: str
    accepted_result_sha256s: tuple[str, ...]
    reference_identity_id: str
    reference_manifest_sha256: str
    reference_bindings_sha256: str
    assessment_scope_sha256: str
    policy_bindings: PolicyBindings
    evaluator_pin_sha256: str
    validation_plan_sha256: str
    validation_roster_sha256: str
    executable_stack_sha256: str
    runtime_image_digest: str
    resource_manifest_sha256s: tuple[str, ...]
    evidence_interface_sha256: str
    canonical_handoff_schema_sha256: str
    canonical_handoff_path_sha256: str
    path_selection_policy_sha256: str
    release_graph_parent_sha256s: tuple[str, ...]
    schema: Literal["genosapiens.track1-base-release"] = SCHEMA
    classification: Literal["public"] = "public"

    def __post_init__(self) -> None:
        if type(self.release_id) is not str or IDENTIFIER.fullmatch(self.release_id) is None:
            raise Track1ReleaseError("invalid release identity")
        if type(self.reference_identity_id) is not str or IDENTIFIER.fullmatch(self.reference_identity_id) is None:
            raise Track1ReleaseError("invalid release reference identity")
        hashes = (
            self.candidate_sha256,
            self.consolidated_acceptance_sha256,
            *self.accepted_result_sha256s,
            self.reference_manifest_sha256,
            self.reference_bindings_sha256,
            self.assessment_scope_sha256,
            self.evaluator_pin_sha256,
            self.validation_plan_sha256,
            self.validation_roster_sha256,
            self.executable_stack_sha256,
            *self.resource_manifest_sha256s,
            self.evidence_interface_sha256,
            self.canonical_handoff_schema_sha256,
            self.canonical_handoff_path_sha256,
            self.path_selection_policy_sha256,
            *self.release_graph_parent_sha256s,
        )
        for value in hashes:
            _sha(value)
        expected_parents = tuple(sorted({
            self.candidate_sha256,
            self.consolidated_acceptance_sha256,
            *self.accepted_result_sha256s,
        }))
        if (
            self.schema != SCHEMA
            or self.classification != "public"
            or self.authorization_status not in {"candidate_only", "accepted"}
            or type(self.policy_bindings) is not PolicyBindings
            or not self.accepted_result_sha256s
            or self.accepted_result_sha256s
            != tuple(sorted(set(self.accepted_result_sha256s)))
            or not self.resource_manifest_sha256s
            or self.resource_manifest_sha256s
            != tuple(sorted(set(self.resource_manifest_sha256s)))
            or self.release_graph_parent_sha256s != expected_parents
            or not re.fullmatch(r"sha256:[0-9a-f]{64}", self.runtime_image_digest)
        ):
            raise Track1ReleaseError("release freeze is incomplete")

    @property
    def ranking_policy_sha256(self) -> str:
        return self.policy_bindings.fusion_policy_sha256

    @property
    def epcr_policy_sha256(self) -> str:
        return self.policy_bindings.epcr_policy_sha256

    def require_analysis_profile(self, selection) -> None:
        """Current consumer scope uses existing resource/plan commitments.

        Loading immutable producer releases does not add fields or rewrite bytes.
        A current consumer must separately require the current policy resource.
        """
        from genosapiens.track1_analysis_profiles import AnalysisProfileSelection
        if (type(selection) is not AnalysisProfileSelection
                or selection.policy.sha256 not in self.resource_manifest_sha256s):
            raise Track1ReleaseError('release analysis profile binding is missing or stale')

    def require_focused_work(self, contract) -> None:
        """Authenticate current execution without altering old release bytes."""
        from genosapiens.wgs_work_contract import FocusedWorkContract
        if (type(contract) is not FocusedWorkContract
                or not set(contract.bindings().values()) <= set(self.resource_manifest_sha256s)):
            raise Track1ReleaseError('release focused work contract is missing or stale')

    def as_bytes(self) -> bytes:
        return _bytes(asdict(self))

    @property
    def sha256(self) -> str:
        return hashlib.sha256(self.as_bytes()).hexdigest()


def freeze_base_release(
    candidate: PreHeldOutCandidate,
    acceptance: ConsolidatedAcceptance,
    *,
    reference_identity_id: str,
    reference_manifest_sha256: str,
    reference_bindings_sha256: str,
    assessment_scope_sha256: str,
    evidence_interface_sha256: str,
    canonical_handoff_schema_sha256: str,
    canonical_handoff_path_sha256: str,
    path_selection_policy_sha256: str,
    analysis_profile: str = 'focused',
) -> BaseRelease:
    """Freeze an accepted release without adding a descendant hash upstream."""
    from genosapiens.track1_validation_runner import (
        ConsolidatedAcceptance,
        PreHeldOutCandidate,
    )

    if (
        type(candidate) is not PreHeldOutCandidate
        or type(acceptance) is not ConsolidatedAcceptance
        or acceptance.status != "accepted"
        or acceptance.candidate_sha256 != candidate.sha256
    ):
        raise Track1ReleaseError("public acceptance does not authorize a release")
    from genosapiens.track1_analysis_profiles import authenticate_analysis_profile
    scope = authenticate_analysis_profile(analysis_profile)
    if scope.policy.sha256 not in candidate.resource_manifest_sha256s:
        raise Track1ReleaseError('candidate analysis profile policy is unbound')
    accepted_results = tuple(sorted({
        item.artifact_sha256 for item in acceptance.result_references
    }))
    release_id = "track1-base-release-" + hashlib.sha256(_bytes({
        "candidate": candidate.sha256,
        "acceptance": acceptance.sha256,
        "results": accepted_results,
    })).hexdigest()[:16]
    parents = tuple(sorted({candidate.sha256, acceptance.sha256, *accepted_results}))
    return BaseRelease(
        release_id,
        "accepted",
        candidate.sha256,
        acceptance.sha256,
        accepted_results,
        reference_identity_id,
        reference_manifest_sha256,
        reference_bindings_sha256,
        assessment_scope_sha256,
        candidate.policy_bindings,
        candidate.evaluator_sha256,
        candidate.plan_sha256,
        candidate.roster_sha256,
        candidate.executable_stack_sha256,
        candidate.runtime_image_digest,
        candidate.resource_manifest_sha256s,
        evidence_interface_sha256,
        canonical_handoff_schema_sha256,
        canonical_handoff_path_sha256,
        path_selection_policy_sha256,
        parents,
    )


def load_base_release(content: bytes, *, expected_sha256: str) -> BaseRelease:
    """Load exact canonical current-contract bytes with duplicate-key rejection."""
    try:
        _sha(expected_sha256)
        if type(content) is not bytes or hashlib.sha256(content).hexdigest() != expected_sha256:
            raise ValueError

        def pairs(items):
            value = {}
            for key, child in items:
                if key in value:
                    raise ValueError
                value[key] = child
            return value

        data = json.loads(content, object_pairs_hook=pairs)
        fields = {
            "release_id",
            "authorization_status",
            "candidate_sha256",
            "consolidated_acceptance_sha256",
            "accepted_result_sha256s",
            "reference_identity_id",
            "reference_manifest_sha256",
            "reference_bindings_sha256",
            "assessment_scope_sha256",
            "policy_bindings",
            "evaluator_pin_sha256",
            "validation_plan_sha256",
            "validation_roster_sha256",
            "executable_stack_sha256",
            "runtime_image_digest",
            "resource_manifest_sha256s",
            "evidence_interface_sha256",
            "canonical_handoff_schema_sha256",
            "canonical_handoff_path_sha256",
            "path_selection_policy_sha256",
            "release_graph_parent_sha256s",
            "schema",
            "classification",
        }
        if type(data) is not dict or set(data) != fields:
            raise ValueError
        release = BaseRelease(
            **{
                **data,
                "policy_bindings": decode_policy_bindings(data["policy_bindings"]),
                "accepted_result_sha256s": tuple(data["accepted_result_sha256s"]),
                "resource_manifest_sha256s": tuple(data["resource_manifest_sha256s"]),
                "release_graph_parent_sha256s": tuple(data["release_graph_parent_sha256s"]),
            }
        )
        if release.as_bytes() != content:
            raise ValueError
        return release
    except (ValueError, TypeError, KeyError, AttributeError, json.JSONDecodeError):
        raise Track1ReleaseError("invalid or mismatched base release") from None
