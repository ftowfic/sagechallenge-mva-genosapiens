#!/usr/bin/env python3
"""Explicit public preparation. This command accepts no case inputs or credentials."""
import argparse
import fcntl
import json
import os
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(ROOT), str(ROOT/'src')]


def main(argv=None):
    declaration = json.loads((ROOT/'config/installation/dependencies.json').read_bytes())
    resources = {row['id']:row for row in declaration['resources']}
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('resource', choices=[*resources, 'discovery'])
    p.add_argument('--reference-root', type=Path, default=Path('/mnt/reference'))
    p.add_argument('--source', type=Path, help='The declared local compressed LINCS matrix archive.')
    p.add_argument('--admit-copied', action='store_true', help='Authenticate copied sources and rebuild affected local metadata.')
    p.add_argument('--describe', action='store_true')
    p.add_argument('--memory-gib', type=int, default=4)
    args = p.parse_args(argv)
    if not 1 <= args.memory_gib <= 65536:
        p.error('--memory-gib requires an integer from 1 through 65536.')
    if args.describe:
        print(json.dumps(resources.get(args.resource, {'producer':'track2 prepare-discovery-resources'}), indent=2))
        return 0
    supported = {'lincs-metadata', 'lincs-matrix', 'fda', 'open-targets-drugs', 'discovery', 'public-index'}
    if args.resource not in supported:
        entry = resources[args.resource]
        p.exit(2, 'Installation limit: '+entry['availability']['limit']+'\n')
    if args.source is not None and args.resource != 'lincs-matrix':
        p.error('--source is accepted only for the compressed LINCS matrix archive.')
    if args.resource == 'lincs-matrix' and not args.admit_copied and args.source is None:
        p.error('Supply --source ARCHIVE or use --admit-copied with the original receipts.')
    from scripts.portable_paths import real_directory
    from genosapiens.portable_resources import (prepare_metadata, prepare_matrix, admit_copied_matrix,
        prepare_open_targets_metadata, checked_path)
    root = real_directory(args.reference_root)
    locks = checked_path(root, 'manifests/installation'); locks.mkdir(parents=True, exist_ok=True)
    descriptor = os.open(locks/'.prepare.lock', os.O_CREAT|os.O_RDWR|os.O_NOFOLLOW, 0o600)
    with os.fdopen(descriptor, 'a') as lock:
        fcntl.flock(lock, fcntl.LOCK_EX|fcntl.LOCK_NB)
        pins = json.loads((ROOT/'config/installation/source-pins.json').read_bytes())
        if args.resource in {'lincs-metadata', 'fda'}:
            result = prepare_metadata(root, args.resource, pins)
        elif args.resource == 'lincs-matrix':
            result = (admit_copied_matrix(root, pins['lincs_archive']) if args.admit_copied
                else prepare_matrix(root, args.source, pins['lincs_archive']))
        elif args.resource == 'public-index':
            from scripts.ensure_track1_installation import ensure
            ensure(ROOT, root)
            result = {'status':'complete', 'installation':'.track1-installation'}
        else:
            contract = json.loads((ROOT/'config/public-resources/track2-open-targets-drugs-26.06.json').read_bytes())
            prepare_open_targets_metadata(root, contract)
            from genosapiens.track2_discovery_resources import acquire_sources
            result = acquire_sources(contract, root)
            if args.resource == 'discovery':
                if args.admit_copied:
                    admit_copied_matrix(root, pins['lincs_archive'])
                    prepare_metadata(root, 'lincs-metadata', pins)
                    prepare_metadata(root, 'fda', pins)
                from genosapiens.track2_discovery_resources import admit_discovery_inputs
                source_contract = json.loads((ROOT/'config/track2/discovery-sources.json').read_bytes())
                admit_discovery_inputs(source_contract, root)
                from scripts.prepare_track2_discovery import ensure_discovery_resources
                result = ensure_discovery_resources(ROOT, root, memory_gib=args.memory_gib)
        print(json.dumps(result, sort_keys=True))
    return 0


if __name__ == '__main__':
    try: raise SystemExit(main())
    except (OSError, ValueError, RuntimeError) as error:
        print(str(error), file=sys.stderr)
        raise SystemExit(2)
