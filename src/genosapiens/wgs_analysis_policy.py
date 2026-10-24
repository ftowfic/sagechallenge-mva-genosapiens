"""Current downstream policy and runtime-profile selection for public analysis."""
import hashlib
from pathlib import Path

from genosapiens.bounded_json import read_bounded_bytes
from genosapiens.track1_release import evidence_interface_bytes
from genosapiens.track1_analysis_profiles import (
    AnalysisProfileError, load_analysis_profile_policy, resolve_analysis_profile,
)

POLICY_ROLES = ("diagnosis_prior_policy", "fusion_policy", "hypothesis_policy")
ROOT = Path(__file__).resolve().parents[2]


class AnalysisPolicyError(ValueError):
    """Static failure of the public downstream policy boundary."""


def _require(value):
    if not value:
        raise AnalysisPolicyError("current analysis policy binding is invalid")


def _sha(content):
    return hashlib.sha256(content).hexdigest()


def select_analysis_profile(*, exhaustive=False, repository_root=None):
    """Authenticate the current scope before any input or resource work.

    Phase 1 recognizes exhaustive but refuses it before even reading policy.
    """
    if type(exhaustive) is not bool:
        raise AnalysisPolicyError('analysis profile selection is invalid')
    if exhaustive:
        raise AnalysisPolicyError('profile_not_qualified')
    root = ROOT if repository_root is None else Path(repository_root)
    try:
        policy = load_analysis_profile_policy(root / 'config/track1/analysis-profiles.json')
        profile = resolve_analysis_profile(policy, exhaustive=exhaustive)
    except AnalysisProfileError as error:
        raise AnalysisPolicyError(str(error)) from None
    return profile, policy


def current_interfaces(*, repository_root=None):
    root = ROOT if repository_root is None else Path(repository_root)
    return dict(evidence_interface_sha256=_sha(evidence_interface_bytes(root)),
                canonical_handoff_schema_sha256=_sha((root / "schemas/track1/mechanism-ready-handoff.schema.json").read_bytes()))


def load_analysis_selectors(authority, *, repository_root=None):
    """Load the sole current screening policy and disposable writer profile."""
    from genosapiens.wgs_analysis_eligibility import load_analysis_eligibility_policy
    from genosapiens.wgs_context_runtime import load_context_writer_profile
    root=ROOT if repository_root is None else Path(repository_root)
    values={row.role:row.sha256 for row in authority.resource_bindings}
    required={'analysis_eligibility_policy','context_writer_profile','analysis_profile_policy'}
    if not required<=values.keys():
        raise AnalysisPolicyError('current analysis selectors are missing')
    if {row.role for row in (*authority.resource_bindings, *authority.policy_bindings)} & {'completed_annotation_policy_transition',
            'completed_admission_eligibility_transition'}:
        raise AnalysisPolicyError('retired analysis selector is not supported')
    _, scope = select_analysis_profile(repository_root=root)
    if values['analysis_profile_policy'] != scope.sha256:
        raise AnalysisPolicyError('current analysis profile binding is stale')
    policy=load_analysis_eligibility_policy(read_bounded_bytes(root/'config/track1/analysis-eligibility-v1.json',
        maximum_bytes=65536),values['analysis_eligibility_policy'])
    profile=load_context_writer_profile(read_bounded_bytes(root/'config/track1/context-writer-runtime-v1.json',
        maximum_bytes=65536),values['context_writer_profile'])
    if policy.common_filter_policy_sha256!=authority.binding('common_allele_policy'):
        raise AnalysisPolicyError('screening and original common policies differ')
    return policy,profile


def select_analysis_runtime(authority, release, policies, *, repository_root=None, exhaustive=False):
    """The same exact runtime selection for preflight and real scientific start."""
    from genosapiens.wgs_public_qualification import require_public_qualification_release
    select_analysis_profile(exhaustive=exhaustive, repository_root=repository_root)
    require_public_qualification_release(authority, release)
    from genosapiens.track1_analysis_profiles import authenticate_analysis_profile
    release.require_analysis_profile(authenticate_analysis_profile())
    load_analysis_selectors(authority,repository_root=repository_root)
    work = load_analysis_work_contract(authority, repository_root=repository_root)
    if work is not None:
        release.require_focused_work(work)
    for name, expected in current_interfaces(repository_root=repository_root).items():
        _require(getattr(release,name) == expected)
    matches = tuple(r for r in policies.runtimes if r.candidate.bindings == release.policy_bindings)
    if len(matches) != 1:
        raise AnalysisPolicyError("release ranking policy is unavailable")
    declared = {r.role: r.sha256 for r in authority.policy_bindings}
    if "public_ranking_policy_family" in declared:
        _require(declared["public_ranking_policy_family"] == policies.policy_sha256)
    for role in POLICY_ROLES:
        if role in declared:
            _require(declared[role] == getattr(release.policy_bindings, role + "_sha256"))
    return matches[0]


def load_analysis_work_contract(authority, *, repository_root=None):
    """Legacy receipts remain readable; fresh qualification requires D41 work."""
    from genosapiens.wgs_work_contract import (
        ROLES, FRESH_INVOCATION, load_focused_work_contract, require_work_bindings,
    )
    values = {row.role: row.sha256 for row in authority.resource_bindings}
    if not ROLES.intersection(values) and authority.invocation_id != FRESH_INVOCATION:
        return None
    work = load_focused_work_contract(repository_root)
    require_work_bindings(values, work)
    return work


def preflight_analysis_policy(authority, release, *, repository_root, reference_root, exhaustive=False):
    """Authenticate the small pinned model-locus source and select real policy.

No large analysis resources, input parsing, VEP, admission or science are run.
"""
    from genosapiens.track1_public_validation import load_public_model_loci, load_public_policy_family
    mode, scope_policy = select_analysis_profile(exhaustive=exhaustive, repository_root=repository_root)
    policies = load_public_policy_family(repository_root,
        model_loci=load_public_model_loci(repository_root, reference_root))
    runtime = select_analysis_runtime(authority, release, policies, repository_root=repository_root,
        exhaustive=exhaustive)
    policy,profile=load_analysis_selectors(authority,repository_root=repository_root)
    work=load_analysis_work_contract(authority,repository_root=repository_root)
    return dict(status="passed", analysis_profile=mode, analysis_profile_policy_sha256=scope_policy.sha256,
        policy_family_sha256=policies.policy_sha256,
        policy_bindings_sha256=runtime.candidate.bindings.sha256,
        analysis_eligibility_policy_sha256=policy.sha256,context_writer_profile_sha256=profile.sha256,
        **({} if work is None else work.bindings()))
