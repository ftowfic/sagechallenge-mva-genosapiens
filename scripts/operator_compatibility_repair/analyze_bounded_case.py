#!/usr/bin/env python3
"""Publish the prepared bounded Track 1 primary; optional analyses run separately."""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import Mapping

from genosapiens.wgs_analysis import BoundedAnalysisResult, analyze_bounded_case
from genosapiens.wgs_runtime import managed_runtime_phase_recorder
from genosapiens.track1_analysis_profiles import POLICY_SHA256, authenticate_analysis_profile


SHARED_CONTRACT_FIELDS = frozenset({
    "policy_sha256",
    "population_manifest_sha256",
    "clinical_manifest_sha256",
    "exception_manifest_sha256",
    "precision_manifest_sha256",
    "ontology_manifest_sha256",
    "execution_config",
    "base_release_bytes",
    "expected_base_release_sha256",
    "executable_stack_sha256",
    "runtime_image_digest",
})


def execute_analysis(arguments: Mapping[str, object]) -> BoundedAnalysisResult:
    if type(arguments) is not dict:
        raise ValueError("bounded analysis arguments are required")
    selection = authenticate_analysis_profile(arguments.get('analysis_profile', 'focused'),
        expected_policy_sha256=arguments.get('analysis_profile_policy_sha256', POLICY_SHA256))
    if not SHARED_CONTRACT_FIELDS.issubset(arguments):
        raise ValueError("bounded analysis shared contract is incomplete")
    return analyze_bounded_case(**{**arguments, 'analysis_profile': selection.profile,
        'analysis_profile_policy_sha256': selection.policy.sha256})


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--resume", action="store_true")
    parser.add_argument('--exhaustive', action='store_true',
        help='analysis-only broader profile; requires matching qualification (default: focused)')
    parser.add_argument('--case-binding', type=Path)
    parser.add_argument('--case-binding-sha256')
    parsed = parser.parse_args(argv)
    authenticate_analysis_profile('exhaustive' if parsed.exhaustive else 'focused')
    if parsed.case_binding is None or parsed.case_binding_sha256 is None:
        parser.error('an immutable package and exact private case binding are required')
    from genosapiens.wgs_operator_invocation import run_operator_stage
    with managed_runtime_phase_recorder(expected_stage='analysis') as recorder:
        return run_operator_stage('analysis', package_root=Path(__file__).resolve().parent,
            parsed=parsed, execute=execute_analysis, phase_recorder=recorder)


if __name__ == "__main__":
    raise SystemExit(main())
