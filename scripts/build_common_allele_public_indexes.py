#!/usr/bin/env python3
"""Build the complete pinned public common-allele protection indexes."""

from __future__ import annotations

import argparse
from dataclasses import asdict
import json
from pathlib import Path

from genosapiens.common_allele_public_indexes import build_common_allele_public_indexes


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--repository-root", type=Path, required=True)
    parser.add_argument("--reference-root", type=Path, required=True)
    arguments = parser.parse_args()
    from genosapiens.bounded_json import load_bounded_json
    parameters, _ = load_bounded_json(arguments.repository_root / 'config/public-resources/track1-retention-index.json')
    result = build_common_allele_public_indexes(repository_root=arguments.repository_root,
        reference_root=arguments.reference_root, **parameters)
    print(json.dumps(asdict(result), sort_keys=True, separators=(",", ":")))


if __name__ == "__main__":
    main()
