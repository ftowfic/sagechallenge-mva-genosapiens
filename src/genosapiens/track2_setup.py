"""Automatic public-only package preparation using Track 1's index owner."""
import argparse
from contextlib import redirect_stdout
import fcntl
import json
import os
from pathlib import Path
import shutil
import sys

from genosapiens.track1_setup import prepare_public_index, atomic_json
from genosapiens.public_disease_lookup import prepare_disease_lookup
from genosapiens.track2_operator_inputs import require, small_json
from genosapiens.track2_published_genotype import canonical, sha256
from genosapiens.track2_report_gene_metadata import installation_inputs


def source_files(repository):
    from genosapiens.track2_package_files import runtime_files
    yield from runtime_files(repository)


def prepare_installation(repository, reference, store, *, index_store=None,discovery_resources=None):
    """Read public metadata only. The host never mounts case paths here."""
    from genosapiens.track2_operator_resources import export_installed
    from genosapiens.track2_operator_install import build_package, validate_package
    repository, reference, store = (Path(p).resolve(strict=True) for p in (repository, reference, store))
    code = {path.relative_to(repository).as_posix(): sha256(path.read_bytes()) for path in source_files(repository)}
    companion_inputs = installation_inputs(repository, reference)
    index_store = Path(index_store or store)
    descriptor_path=repository/'config/track2/installed-public-evidence.json'
    descriptor = small_json(descriptor_path) if descriptor_path.exists() else None
    if discovery_resources is not None and (descriptor is None or not (reference/descriptor['path']/'evidence.json').is_file()):
        descriptor=None
    if descriptor is None:
        require(discovery_resources is not None,'The public study export is missing.')
        evidence=envelope=None
    else:
        evidence,envelope=_study_paths(descriptor,reference)
    discovery=None
    if discovery_resources is not None:
        from genosapiens.track2_discovery_resources import load_discovery_resources
        load_discovery_resources(discovery_resources,reference,matrix_required=False).close()
        discovery=small_json(discovery_resources)
    pin = prepare_public_index(repository, reference, index_store)
    disease_lookup = prepare_disease_lookup(reference, pin, index_store)
    request = sha256(canonical({'code': code, 'index': pin, 'disease_lookup': disease_lookup, 'evidence': descriptor,
        **({'discovery':discovery} if discovery is not None else {}),'reference_root': str(reference),
        'gene_metadata_inputs': companion_inputs}))
    saved = store/'current-package.json'
    if saved.exists():
        record = small_json(saved)
        if record.get('request_sha256') == request:
            require(record.get('package_name') == 'package-'+request, 'The current package location is invalid.')
            validate_package(store/record['package_name'], record['package_sha256'])
            print('Reuse the completed Track 2 package and public index.', file=sys.stderr, flush=True)
            return record
    target = store/('package-'+request)
    atomic_json(store/'preparation-progress.json', {'status': 'running', 'request_sha256': request})
    try:
        if target.exists() and not (target/'package.json').exists():
            require(not target.is_symlink() and not any(p.is_symlink() for p in target.rglob('*')),
                'Incomplete package preparation contains a symbolic link.')
            shutil.rmtree(target)
        if not target.exists():
            resources = store/'resources.json'
            resources.unlink(missing_ok=True)
            with redirect_stdout(sys.stderr):
                export_installed(repository, reference, resources, retention_index=pin,
                    retention_disease_lookup=disease_lookup)
                result = build_package(repository, resources, envelope, target,
                    public_evidence=evidence,**({'discovery_resources':discovery_resources} if discovery is not None else {}))
        else:
            digest = sha256((target/'package.json').read_bytes())
            manifest = validate_package(target, digest)
            result = {'package_sha256': digest, 'image': manifest.get('image',manifest.get('images',{}).get('diagnosis')),
                **({'images':manifest['images']} if 'images' in manifest else {})}
        record = {'request_sha256': request, 'package_name': target.name,
            'package_sha256': result['package_sha256'], 'image': result['image'],
            **({'images':result['images']} if 'images' in result else {})}
        atomic_json(saved, record)
        atomic_json(store/'preparation-progress.json', {'status': 'complete', **record})
        return record
    except Exception:
        atomic_json(store/'preparation-progress.json', {'status': 'incomplete', 'request_sha256': request})
        raise


def _study_paths(descriptor,reference):
    require(descriptor.get('schema') == 'genosapiens.track2-installed-public-evidence/v1',
        'The installed public evidence descriptor is invalid.')
    evidence = reference/descriptor['path']
    require(evidence.is_relative_to(reference) and evidence.resolve() == evidence.absolute()
        and (evidence/'evidence.json').is_file(),
        'Frozen Track 2 public evidence is missing. Run ./track2 install-public-evidence --source PATH.')
    require(sha256((evidence/'evidence.json').read_bytes()) == descriptor['evidence_sha256'],
        'The frozen public evidence version differs. Restore the declared public export.')
    envelope = reference/descriptor['envelope_path'] if 'envelope_path' in descriptor else evidence/'envelope.json'
    require(envelope.is_relative_to(reference) and envelope.resolve() == envelope.absolute(),
        'The frozen public envelope path is invalid.')
    if 'envelope_sha256' in descriptor:
        require(sha256(envelope.read_bytes()) == descriptor['envelope_sha256'],
            'The frozen public envelope version differs.')
    return evidence,envelope


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ('repository', 'reference', 'store', 'index-store'):
        parser.add_argument('--'+name, required=True, type=Path)
    parser.add_argument('--discovery-resources',type=Path)
    args = parser.parse_args(argv)
    for root in (args.store, args.index_store):
        require(root.is_dir() and root.resolve() == root.absolute(), 'Public preparation storage is missing or linked.')
    # Track 1's host uses this same lock for the shared preparation owner.
    descriptor = os.open(args.index_store/'.setup.lock', os.O_CREAT | os.O_RDWR | os.O_NOFOLLOW, 0o600)
    try:
        fcntl.flock(descriptor, fcntl.LOCK_EX)
        record = prepare_installation(args.repository, args.reference, args.store, index_store=args.index_store,
            discovery_resources=args.discovery_resources)
    finally:
        os.close(descriptor)
    print(json.dumps(record, sort_keys=True))


if __name__ == '__main__':
    main()
