#!/usr/bin/env python3
"""Repository-owned entry point for bounded Track 1 preparation.

The immutable operator package authenticates one private case descriptor and
dispatches the shared preparation assembly.
"""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import Mapping

from genosapiens.wgs_preparation import (
    BoundedPreparationResult,
    PreparationRunOptions,
    prepare_bounded_case,
)
from genosapiens.wgs_runtime import managed_runtime_phase_recorder


SHARED_CONTRACT_FIELDS = frozenset({
    "policy_sha256",
    "scope_policy_sha256",
    "scope_manifest_sha256",
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
RUN_OPTION_BINDING_FIELDS = frozenset({
    "default_preparation_run_options_sha256",
    "audited_preparation_run_options_sha256",
})


def execute_preparation(arguments: Mapping[str, object]) -> BoundedPreparationResult:
    """Execute one package-bound request without exposing controlled details."""

    if type(arguments) is not dict:
        raise ValueError("bounded preparation arguments are required")
    if not (SHARED_CONTRACT_FIELDS | RUN_OPTION_BINDING_FIELDS).issubset(arguments):
        raise ValueError("bounded preparation shared contract is incomplete")
    if type(arguments.get("run_options")) is not PreparationRunOptions:
        raise ValueError("bounded preparation run options are incomplete")
    if (arguments["default_preparation_run_options_sha256"]
            != PreparationRunOptions(False).sha256
            or arguments["audited_preparation_run_options_sha256"]
            != PreparationRunOptions(True).sha256):
        raise ValueError("bounded preparation option bindings are invalid")
    selected = ("audited_preparation_run_options_sha256"
                if arguments["run_options"].verify_source_ledger
                else "default_preparation_run_options_sha256")
    if arguments[selected] != arguments["run_options"].sha256:
        raise ValueError("bounded preparation audit mode changed")
    bound = {key: value for key, value in arguments.items()
             if key not in RUN_OPTION_BINDING_FIELDS}
    if bound.get("completed_source_input") is not None:
        from genosapiens.wgs_preparation import require_completed_source_readonly
        require_completed_source_readonly(bound["completed_source_input"])
    return prepare_bounded_case(**bound)


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--resume", action="store_true")
    parser.add_argument("--verify-source-ledger", action="store_true")
    parser.add_argument('--case-binding', type=Path)
    parser.add_argument('--case-binding-sha256')
    parsed = parser.parse_args(argv)
    if parsed.case_binding is None or parsed.case_binding_sha256 is None:
        parser.error('an immutable package and exact private case binding are required')
    from genosapiens.wgs_operator_invocation import run_operator_stage
    with managed_runtime_phase_recorder(expected_stage='preparation') as recorder:
        return run_operator_stage('preparation', package_root=Path(__file__).resolve().parent,
            parsed=parsed, execute=execute_preparation, phase_recorder=recorder)


if __name__ == "__main__":
    raise SystemExit(main())
