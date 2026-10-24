#!/usr/bin/env python3
"""Build the current primary bounded-WGS package from a qualified release."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

from genosapiens.wgs_operator_package import _source, build_bounded_operator_package
from genosapiens.wgs_analysis_policy import select_analysis_profile


def main(argv=None, *, completed_source_input=None) -> int:
    """Accept a completed input only from the immutable caller binding."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repository-root", type=Path, required=True)
    parser.add_argument("--destination-parent", type=Path, required=True)
    parser.add_argument("--base-release", type=Path, required=True)
    parser.add_argument("--base-release-sha256", required=True)
    parser.add_argument("--execution-config", type=Path, required=True)
    parser.add_argument("--execution-config-sha256", required=True)
    parser.add_argument("--common-allele-policy", type=Path, required=True)
    parser.add_argument("--common-allele-policy-sha256", required=True)
    parser.add_argument("--population-manifest-sha256", required=True)
    parser.add_argument("--clinical-manifest-sha256", required=True)
    parser.add_argument("--ontology-manifest-sha256", required=True)
    parser.add_argument("--preparation-runtime-evidence", type=Path, required=True)
    parser.add_argument("--analysis-runtime-evidence", type=Path, required=True)
    parser.add_argument('--exhaustive', action='store_true',
        help='package broader analysis only after matching qualification (default: focused)')
    arguments = parser.parse_args(argv)
    select_analysis_profile(exhaustive=arguments.exhaustive, repository_root=arguments.repository_root)
    package = build_bounded_operator_package(
        repository_root=arguments.repository_root,
        completed_source_input=completed_source_input,
        destination_parent=arguments.destination_parent,
        base_release_bytes=_source(arguments.base_release),
        expected_base_release_sha256=arguments.base_release_sha256,
        execution_config_bytes=_source(arguments.execution_config),
        expected_execution_config_sha256=arguments.execution_config_sha256,
        common_allele_policy_bytes=_source(arguments.common_allele_policy),
        expected_common_allele_policy_sha256=
            arguments.common_allele_policy_sha256,
        population_manifest_sha256=arguments.population_manifest_sha256,
        clinical_manifest_sha256=arguments.clinical_manifest_sha256,
        ontology_manifest_sha256=arguments.ontology_manifest_sha256,
        preparation_runtime_evidence_bytes=_source(arguments.preparation_runtime_evidence),
        analysis_runtime_evidence_bytes=_source(arguments.analysis_runtime_evidence),
        exhaustive=arguments.exhaustive,
    )
    print(json.dumps({
        "schema": "genosapiens.bounded-wgs-operator-package-static-status",
        "version": 1,
        "status": "prepared",
        "package_name": package.root.name,
        "package_bindings_sha256": package.bindings_sha256,
        "primary_required_roles": package.contract.primary_required_roles,
        "optional_analysis_required": package.contract.optional_analysis_required,
        'analysis_profile': package.contract.analysis_profile,
        'analysis_profile_policy_sha256': package.contract.analysis_profile_policy_sha256,
        "controlled_use_authorized": False,
    }, sort_keys=True, separators=(",", ":")))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
