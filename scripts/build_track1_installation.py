#!/usr/bin/env python3
"""Prepare Track 1 from installed public resources without case inputs."""
import argparse
import json
import os
from pathlib import Path

from genosapiens.track1_setup import prepare_installation, prepare_public_index, prepare_early_scope, public_versions
from genosapiens.wgs_installation import build_package


def main():
    os.umask(0o077)
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--repository-root', type=Path, required=True)
    parser.add_argument('--reference-root', type=Path, required=True)
    parser.add_argument('--package-store', type=Path)
    parser.add_argument('--destination', type=Path)
    args = parser.parse_args()
    if bool(args.package_store) == bool(args.destination):
        parser.error('Supply either --package-store or --destination.')
    repository, reference = args.repository_root.resolve(), args.reference_root.resolve()
    store = (args.package_store or repository / '.track1-packages').resolve()
    store.mkdir(mode=0o700, parents=True, exist_ok=True)
    try:
        if args.destination:
            pin = prepare_public_index(repository, reference, store)
            scope_pin = prepare_early_scope(repository, reference, store, pin)
            versions, resources, bindings = public_versions(repository, reference, pin, scope_pin)
            sha = build_package(repository_root=repository, destination=args.destination.resolve(),
                reference_versions=versions, resource_versions=resources, policy_bindings=bindings,
                public_index_pin=pin, scope_pin=scope_pin)
            print(sha)
        else:
            print(json.dumps(prepare_installation(repository, reference, store), sort_keys=True))
    except (ValueError, RuntimeError, OSError) as error:
        parser.exit(2, f'Track 1 preparation failed: {error}\n')


if __name__ == '__main__':
    main()
