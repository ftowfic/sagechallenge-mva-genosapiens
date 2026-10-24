#!/usr/bin/env python3
"""Qualify and publish the immutable public common-allele exception bundle."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

from genosapiens.common_allele_exception_bundle import build_exception_bundle


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--reference-root", required=True, type=Path)
    parser.add_argument("--repository-root", required=True, type=Path)
    parser.add_argument(
        "--descriptor",
        default="config/public-resources/track1-common-allele-exceptions-v1.json",
    )
    parser.add_argument(
        "--curation",
        default="config/public-resources/track1-common-allele-exception-curation-v1.json",
    )
    arguments = parser.parse_args()
    repository = arguments.repository_root.resolve()
    result = build_exception_bundle(
        repository_root=repository,
        reference_root=arguments.reference_root,
        descriptor_path=repository / arguments.descriptor,
        curation_path=repository / arguments.curation,
    )
    print(json.dumps({field: getattr(result, field) for field in result.__dataclass_fields__},
                     sort_keys=True, separators=(",", ":")))


if __name__ == "__main__":
    main()
