"""Preserve and bind an interrupted annotation publication before a resume."""
from pathlib import Path
import shutil
from uuid import uuid4

from genosapiens.track2_annotation_checkpoint import dependency_binding
from genosapiens.track2_operator_inputs import require, small_json, write_new
from genosapiens.track2_published_genotype import sha256
from genosapiens.bounded_json import read_bounded_bytes


def prepare_recovery(root, *, package, package_store, image):
    """Caller owns the run lock. Never accept changed producer dependencies.

    Preserve the original run record and diagnostics before the launcher
    replaces them. The worker still checks exact requests and replays admission.
    """
    root, package = Path(root), Path(package)
    checkpoint = root/'annotation-checkpoint'
    annotation = checkpoint/'annotation'
    if not annotation.exists() or (checkpoint/'annotation-completion.json').exists():
        return
    require(annotation.is_dir() and annotation.resolve() == annotation.absolute(),
        'Annotation recovery directory is missing or linked.')
    resources = small_json(package/'resources.json', maximum=128*1024*1024)
    current = dependency_binding((), resources['vep'], {'image': image}, package)
    current.pop('requests_sha256')
    destination = checkpoint/'annotation-recovery.json'
    if destination.exists():
        record = small_json(destination)
        require(record.get('schema') == 'genosapiens.track2-annotation-recovery/v1',
            'Annotation recovery dependencies changed. Original annotation was preserved.')
        if record.get('producer_binding') != current:
            # A record made by the original host code used Python's default AST
            # spelling. Accept only an exact reproduction of that old binding
            # from the current source, resources, and runtime. No dependencies
            # may be dropped or changed during this format correction.
            legacy = dependency_binding((), resources['vep'], {'image': image}, package,
                explicit_empty_ast=False)
            legacy.pop('requests_sha256')
            require(record.get('producer_binding') == legacy,
                'Annotation recovery dependencies changed. Original annotation was preserved.')
            raw = read_bounded_bytes(destination, maximum_bytes=16*1024*1024)
            backup = checkpoint/('annotation-recovery-before-ast-'+sha256(raw)+'.json')
            if backup.exists():
                require(read_bounded_bytes(backup, maximum_bytes=16*1024*1024) == raw,
                    'Preserved annotation recovery metadata changed.')
            else:
                write_new(backup, raw)
            from genosapiens.track2_run_state import atomic_record
            atomic_record(destination, {**record, 'producer_binding': current})
            print('Updated annotation recovery fingerprints for the worker Python version.', flush=True)
        return
    saved = small_json(root/'current-run.json')
    versions = saved['source_versions']
    digest = versions['package_sha256']
    import re
    require(isinstance(digest, str) and re.fullmatch('[0-9a-f]{64}', digest),
        'The original annotation package identity is invalid.')
    command = small_json(root/'current-diagnostics/construct.command.json')
    require(isinstance(command, list), 'The original annotation command is unavailable.')
    paths = [v.split('=', 1)[1] for v in command
        if isinstance(v, str) and v.startswith('GENOSAPIENS_TRACK2_PACKAGE=')]
    require(len(paths) == 1, 'The original annotation package location is unavailable.')
    original = Path(paths[0])
    require(original.is_absolute() and original.resolve() == original
        and original.parent == Path(package_store).resolve(),
        'The original annotation package location is invalid.')
    from genosapiens.track2_operator_install import validate_package
    manifest = validate_package(original, digest)
    original_image = manifest.get('image', manifest.get('images', {}).get('diagnosis'))
    require(original_image == versions['image'], 'The original annotation runtime identity changed.')
    old_resources = small_json(original/'resources.json', maximum=128*1024*1024)
    previous = dependency_binding((), old_resources['vep'], {'image': original_image}, original)
    previous.pop('requests_sha256')
    require(previous == current,
        'Annotation producer dependencies changed. Original annotation was preserved. Explicit reconstruction is required.')
    binding = read_bounded_bytes(annotation/'binding.json', maximum_bytes=64*1024*1024)
    evidence = checkpoint/('recovery-evidence-'+uuid4().hex)
    evidence.mkdir(mode=0o700)
    shutil.copyfile(root/'current-run.json', evidence/'current-run.json')
    diagnostics = root/'current-diagnostics'
    if diagnostics.exists():
        require(diagnostics.resolve() == diagnostics.absolute()
            and not any(p.is_symlink() for p in diagnostics.rglob('*')),
            'Recovery diagnostics contain a symbolic link.')
        shutil.copytree(diagnostics, evidence/'current-diagnostics')
    write_new(destination, {'schema': 'genosapiens.track2-annotation-recovery/v1',
        'producer_binding': previous, 'binding_sha256': sha256(binding),
        'source_package_sha256': digest, 'preserved_evidence': evidence.name})


def validate_recovery(artifacts, *, requests, index_bytes, image):
    from genosapiens.annotation_admission import validate_annotation_admission_files
    from genosapiens.vep_native_annotation import _strict_json
    binding = _strict_json(artifacts['binding.json'].decode('utf-8'))
    require(binding['runtime_image'] == image, 'Annotation recovery runtime changed.')
    # Replay authenticates all output hashes, complete ordered request evidence,
    # reference metadata, and original consequence claims before publication.
    validate_annotation_admission_files({**artifacts, 'reference.fai': index_bytes},
        expected_binding_sha256=sha256(artifacts['binding.json']), requests=requests)
