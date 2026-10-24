"""Immutable local software packaging and explicit operator setup, without case reads."""

import argparse
import json
import os
from pathlib import Path
import re
import shutil

from genosapiens.track2_operator_inputs import require, small_json, write_new, private_path
from genosapiens.track2_published_genotype import canonical, sha256


SCHEMA = 'genosapiens.track2-operator-package/v1'
SETUP_SCHEMA = 'genosapiens.track2-operator-installation/v1'
NATIVE_SCHEMA = 'genosapiens.track2-operator-package/v2'
NATIVE_SETUP_SCHEMA = 'genosapiens.track2-operator-installation/v2'


def build_package(repository, resources, envelope, output, *, public_evidence=None, allow_synthetic=False,discovery_resources=None):
    """Copy only current code/config/schema and public input records; no case paths."""
    from genosapiens.track2_freeze import load_frozen_public_envelope
    from genosapiens.track2_operator_resources import SCHEMA as RESOURCE_SCHEMA, validate_resource_roots
    root, target = Path(repository).resolve(), Path(output).absolute()
    source = small_json(resources, maximum=128*1024*1024)
    require(source.get('schema') == RESOURCE_SCHEMA and source.get('classification') == 'public',
            'Package construction requires generated public resource metadata.')
    validate_resource_roots(source)
    frozen = load_frozen_public_envelope(envelope) if envelope is not None else None
    require(frozen is not None or discovery_resources is not None,'A package requires native discovery or a frozen study.')
    image = small_json(root/'config/track2/diagnosis-runtime.json')['image']
    require(re.fullmatch(r'genosapiens/[a-z0-9-]+@sha256:[0-9a-f]{64}', image),
            'A pinned project runtime image is required.')
    require(target.resolve() == target and not target.exists(), 'Choose a new real package directory.')
    target.mkdir(parents=True, mode=0o700)
    from genosapiens.track2_package_files import runtime_files, SCRIPTS
    for path in runtime_files(root):
        destination = target/path.relative_to(root)
        destination.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(path, destination)
        if path.relative_to(root).as_posix() in SCRIPTS:
            destination.chmod(0o700)
    # Bind configuration paths to the copied code/configuration, not the working tree.
    def relocate(value):
        if isinstance(value, str) and value.startswith(str(root)+'/config/'):
            return str(target)+value[len(str(root)):]
        if isinstance(value, dict):
            return {k:relocate(v) for k,v in value.items()}
        if isinstance(value, list):
            return [relocate(v) for v in value]
        return value
    write_new(target/'resources.json', relocate(source))
    from genosapiens.track2_report_gene_metadata import prepare_gene_metadata, FILENAME
    from genosapiens.track2_codec import decode
    from genosapiens.admission_links import PinnedPublicSource
    if source.get('gene_source') is not None:
        prepare_gene_metadata(decode(PinnedPublicSource, source['gene_source']), target/FILENAME)
    if source.get('retention_index') is not None:
        pin=target/'config/track1/installed-public-index.json'
        pin.parent.mkdir(parents=True,exist_ok=True)
        pin.write_bytes(canonical(source['retention_index']))
    if frozen is not None:write_new(target/'envelope.json', frozen)
    if discovery_resources is not None:
        write_new(target/'discovery-resources.json',small_json(discovery_resources))
    from genosapiens.track2_public_matching import load_public_evidence, empty_public_evidence
    if public_evidence is None:
        empty_public_evidence(target/'public-evidence')
    else:
        public = load_public_evidence(public_evidence, allow_synthetic=allow_synthetic)
        if public['summary'] is not None and frozen is not None:
            require(sha256(Path(envelope).read_bytes()) ==
                small_json(Path(public_evidence)/'report/manifest.json')['frozen_envelope_sha256'],
                'Installed public evidence and therapeutic envelope disagree.')
        shutil.copytree(public_evidence, target/'public-evidence')
    files = {p.relative_to(target).as_posix():sha256(p.read_bytes())
             for p in sorted(target.rglob('*')) if p.is_file()}
    images={'diagnosis':image,'discovery':small_json(root/'config/track2/discovery-runtime.json')['image']}
    manifest = {'schema':NATIVE_SCHEMA,'images':images,'files':files} if discovery_resources is not None else {'schema':SCHEMA, 'image':image, 'files':files}
    digest = write_new(target/'package.json', manifest)
    return {'image':image, **({'images':images} if discovery_resources is not None else {}),'package':str(target), 'package_sha256':digest}


def validate_package(root, expected):
    root = Path(root)
    raw = (root/'package.json').read_bytes()
    manifest = json.loads(raw)
    native=manifest.get('schema')==NATIVE_SCHEMA
    require(sha256(raw)==expected and manifest['schema'] in {SCHEMA,NATIVE_SCHEMA}
            and set(manifest)==({'schema','images','files'} if native else {'schema','image','files'}), 'Installed package binding changed.')
    if native:
        require(set(manifest['images'])=={'diagnosis','discovery'} and all(re.fullmatch(r'.+@sha256:[0-9a-f]{64}',v)
            for v in manifest['images'].values()),'Native package runtime images are invalid.')
    actual = {p.relative_to(root).as_posix() for p in root.rglob('*') if p.is_file()}
    require(actual == set(manifest['files'])|{'package.json'}, 'Installed package file closure changed.')
    for name,digest in manifest['files'].items():
        path = root/name
        require(path.resolve()==path.absolute() and sha256(path.read_bytes())==digest,
                'Installed package content changed.')
    return manifest


def load_installation(path):
    settings=small_json(private_path(path))
    require(settings.get('schema') in {SETUP_SCHEMA,NATIVE_SETUP_SCHEMA},'The operator installation is invalid.')
    # Registration authenticates the participant and incident environment.
    # Current software comes from automatic public preparation, as in Track 1.
    settings={**settings, **{key:os.environ.get('GENOSAPIENS_TRACK2_'+key.upper())
        for key in ('package','package_sha256','image')}}
    require(all(settings[key] for key in ('package','package_sha256','image')),
        'Use the ordinary Track 2 launcher to prepare the current package.')
    manifest=validate_package(settings['package'],settings['package_sha256'])
    images=manifest.get('images',{'diagnosis':manifest.get('image')})
    require(settings['image'] in images.values(),'The installed runtime image changed.')
    settings.update(images=images,image=images['diagnosis'])
    return settings,manifest


def configure(package, *, reference_root, governance, verification_key, rules,
              participant_id, incident_root, environment_id, inventory, case_id, destination):
    """Register existing private authorization; never create or sign attestations."""
    package = Path(package).absolute()
    raw = (package/'package.json').read_bytes()
    manifest = validate_package(package, sha256(raw))
    paths = {k:str(Path(v).absolute()) for k,v in {
        'reference_root':reference_root, 'governance':governance,
        'verification_key':verification_key, 'rules':rules, 'incident_root':incident_root,
        'inventory':inventory}.items()}
    require(all('\n' not in v and '\r' not in v for v in (*paths.values(),str(package))),
            'Installation paths must not contain line breaks.')
    native=manifest['schema']==NATIVE_SCHEMA
    settings = {'schema':NATIVE_SETUP_SCHEMA if native else SETUP_SCHEMA, **paths, 'package':str(package),
        'package_sha256':sha256(raw), **({'images':manifest['images']} if native else {'image':manifest['image']}),
        'verification_key_sha256':sha256(Path(verification_key).read_bytes()),
        'rules_sha256':sha256(Path(rules).read_bytes()),
        'participant_id':participant_id, 'environment_id':environment_id, 'case_id':case_id}
    # This validates actual existing signed access before installing the run entry.
    from genosapiens.track2_operator import authorize
    authorize(settings)
    destination = Path(destination).absolute()
    destination.mkdir(mode=0o700, parents=True, exist_ok=True)
    settings_path = destination/'track2-settings.json'
    write_new(settings_path, settings)
    return settings_path



def main(argv=None):
    parser = argparse.ArgumentParser(description='Prepare the local Track 2 operator installation.')
    commands = parser.add_subparsers(dest='command',required=True)
    resources = commands.add_parser('resources',help='export the installed public resource bindings')
    for name in ('repository','reference-root','destination'):
        resources.add_argument('--'+name,required=True,type=Path)
    build = commands.add_parser('package')
    for name in ('repository','resources','envelope','output'):
        build.add_argument('--'+name,required=True,type=Path)
    build.add_argument('--public-evidence', type=Path,
        help='Frozen public result export with reviewed identity mapping or an explicit mapping gap.')
    evidence = commands.add_parser('evidence', help='freeze a completed public analysis for exact offline matching')
    for name in ('bundle','plan','envelope','report','output'):
        evidence.add_argument('--'+name, required=True, type=Path)
    for name in ('mapping','mapping-review'):
        evidence.add_argument('--'+name, type=Path)
    setup = commands.add_parser('configure')
    for name in ('package','reference-root','governance','verification-key','rules','incident-root','inventory','destination'):
        setup.add_argument('--'+name,required=True,type=Path)
    for name in ('participant-id','environment-id','case-id'):
        setup.add_argument('--'+name,required=True)
    args=vars(parser.parse_args(argv));command=args.pop('command')
    if command=='resources':
        from genosapiens.track2_operator_resources import export_installed
        export_installed(**args)
        print('The installed public resource bindings were exported.')
    elif command=='evidence':
        from genosapiens.track2_public_matching import export_public_evidence
        export_public_evidence(**args)
        print('The public evidence export is frozen. Missing reviewed mappings remain explicit.')
    elif command=='package':
        print(json.dumps(build_package(**args)))
    elif command=='configure':
        configure(**args)
        print('The local Track 2 command is installed.')


if __name__ == '__main__':
    main()
