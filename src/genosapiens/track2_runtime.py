"""Two-stage offline execution for a frozen, reviewed diagnosis builder/consumer.

The program is project/operator-owned code, hashed with the runtime overlay.
It calls the existing scientific APIs; this launcher adds no filtering policy.
"""

import json
import math
import os
from pathlib import Path
import re
import shutil
import subprocess
import uuid

from genosapiens.governance import ArtifactClassification, validate_relative_classification
from genosapiens.incident_response import assert_processing_allowed
from genosapiens.participant_governance import authorize_controlled_access
from genosapiens.track2_published_genotype import canonical, sha256
from genosapiens.track2_storage import output_storage


class Track2RuntimeError(ValueError):
    """Static runtime error; exact paths and container logs remain classified."""


def freeze_runtime(repository, program, output, *, image=None):
    """Freeze the current code/config overlay and explicit reviewed stage program."""
    if image is None:
        image = json.loads((Path(repository) / 'config/track2/diagnosis-runtime.json').read_bytes())['image']
    if re.fullmatch(r'genosapiens/[a-z0-9-]+@sha256:[0-9a-f]{64}', image) is None:
        raise Track2RuntimeError('a pinned project runtime image is required')
    root, target = Path(repository), Path(output)
    program = Path(program)
    if program.is_symlink() or not program.is_file():
        raise Track2RuntimeError('stage program is invalid')
    target.mkdir(parents=True, exist_ok=False, mode=0o700)
    for relative in ('src/genosapiens', 'config/track2', 'schemas/track2'):
        for path in sorted((root / relative).rglob('*')):
            if path.suffix in {'.py', '.json'} and path.is_file():
                if path.is_symlink():
                    raise Track2RuntimeError('runtime source links are not permitted')
                destination = target / path.relative_to(root)
                destination.parent.mkdir(parents=True, exist_ok=True)
                shutil.copyfile(path, destination)
    relative = 'config/track1/analysis-profiles.json'
    (target / relative).parent.mkdir(parents=True, exist_ok=True)
    shutil.copyfile(root / relative, target / relative)
    shutil.copyfile(program, target / 'program.py')
    files = {path.relative_to(target).as_posix(): sha256(path.read_bytes())
             for path in sorted(target.rglob('*')) if path.is_file()}
    manifest = {'image': image, 'files': files}
    (target / 'runtime-manifest.json').write_bytes(canonical(manifest))
    return {'image': image, 'code_revision': 'sha256:' + sha256(canonical(manifest))}


def _runtime_mounts(snapshot, runtime):
    snapshot = Path(snapshot)
    raw = (snapshot / 'runtime-manifest.json').read_bytes()
    manifest = json.loads(raw)
    if (set(manifest) != {'image', 'files'} or runtime != {
            'image': manifest['image'], 'code_revision': 'sha256:' + sha256(raw)}
            or raw != canonical(manifest)):
        raise Track2RuntimeError('frozen runtime identity changed')
    if re.fullmatch(r'genosapiens/[a-z0-9-]+@sha256:[0-9a-f]{64}', runtime['image']) is None:
        raise Track2RuntimeError('a pinned project runtime image is required')
    actual = {path.relative_to(snapshot).as_posix(): path for path in snapshot.rglob('*') if path.is_file()}
    if set(actual) != {*manifest['files'], 'runtime-manifest.json'}:
        raise Track2RuntimeError('frozen runtime contains undeclared files')
    for name, digest in manifest['files'].items():
        path = actual[name]
        if path.resolve() != path.absolute() or sha256(path.read_bytes()) != digest:
            raise Track2RuntimeError('frozen runtime bytes changed')
    return [(snapshot / 'src', '/opt/genosapiens/src'),
            (snapshot / 'config/track2', '/opt/genosapiens/config/track2'),
            (snapshot / 'config/track1/analysis-profiles.json', '/opt/genosapiens/config/track1/analysis-profiles.json'),
            (snapshot / 'schemas/track2', '/opt/genosapiens/schemas/track2'),
            (snapshot / 'program.py', '/program.py')]


def run_offline_stage(stage, output, *, runtime_snapshot, runtime, classification,
                      reference_root, public_inputs, raw_vcf=None, hpo=None, package=None,
                      package_sha256=None, review_sha256=None, private_inputs=None,
                      controlled_store=None, governance=None, participant_id=None,
                      incident_root=None, environment_id=None, allow_synthetic=False,
                      timeout_seconds=120):
    """Run one fixed program with a closed set of read-only input mounts.

    Public inputs must live in the shared public roots. Additional private inputs
    are JSON files (job, plan, plan_review, snapshot); directories are forbidden.
    Only construction can mount raw VCF/HPO. Synthesis mounts an exact package.
    The host owns the inventory, audit and access controls; none is mounted into
    the worker. The worker's stdout/stderr are retained in classified storage.
    """
    if stage not in {'construct', 'synthesize'} or os.getuid() == 0:
        raise Track2RuntimeError('an offline stage and non-root host account are required')
    if timeout_seconds is not None and (isinstance(timeout_seconds, bool)
            or not isinstance(timeout_seconds, (int, float))
            or not math.isfinite(timeout_seconds) or timeout_seconds <= 0):
        raise Track2RuntimeError('stage timeout must be positive finite seconds or absent')
    root = Path(reference_root).absolute()
    if root.resolve() != root:
        raise Track2RuntimeError('reference root is invalid')
    if classification == 'controlled':
        if controlled_store is None or controlled_store.reference_root != root:
            raise Track2RuntimeError('controlled runtime requires its retention store')
        authorize_controlled_access(governance, participant_id, os.getuid())
        assert_processing_allowed(incident_root, environment_id)
    elif classification != 'public' or controlled_store is not None:
        raise Track2RuntimeError('runtime classification is invalid')
    mounts = _runtime_mounts(Path(runtime_snapshot).absolute(), runtime)

    def public_path(value):
        path = Path(value).absolute()
        try:
            validate_relative_classification(path.relative_to(root).as_posix(), ArtifactClassification.PUBLIC)
            if path.resolve() != path or not path.exists():
                raise ValueError
        except (ValueError, OSError):
            raise Track2RuntimeError('public runtime input is invalid') from None
        return path

    def case_path(value):
        path = controlled_store.path(value) if classification == 'controlled' else public_path(value)
        if not path.is_file():
            raise Track2RuntimeError('runtime input must be a regular file')
        if classification == 'controlled':
            relative = path.relative_to(root).as_posix()
            if not any(r.relative_path == relative and r.status == 'present'
                       and r.classification is ArtifactClassification.CONTROLLED
                       for r in controlled_store.inventory().artifacts):
                raise Track2RuntimeError('controlled runtime input is not inventoried')
        return path

    for name, value in public_inputs.items():
        if re.fullmatch(r'[a-z][a-z0-9_-]*', name) is None:
            raise Track2RuntimeError('public mount name is invalid')
        mounts.append((public_path(value), '/public/' + name))
    for name, value in (private_inputs or {}).items():
        if name not in {'job', 'plan', 'plan_review', 'snapshot'}:
            raise Track2RuntimeError('private runtime input role is invalid')
        path = case_path(value)
        if type(json.loads(path.read_bytes())) is not dict:
            raise Track2RuntimeError('private runtime input must be a structured record')
        mounts.append((path, '/inputs/' + name + '.json'))
    if stage == 'construct':
        if raw_vcf is None or package is not None:
            raise Track2RuntimeError('construction requires the explicit raw VCF input')
        mounts.append((case_path(raw_vcf), '/raw/input.vcf'))
        if hpo is not None:
            mounts.append((case_path(hpo), '/raw/hpo.json'))
    else:
        from genosapiens.track2_diagnosis_package import load_diagnosis_package
        if raw_vcf is not None or hpo is not None or package is None or not package_sha256 or not review_sha256:
            raise Track2RuntimeError('synthesis accepts only an immutable minimal package')
        package = (controlled_store.path(package) if classification == 'controlled' else public_path(package))
        bundle = load_diagnosis_package(package, expected_input_sha256=package_sha256,
            expected_review_sha256=review_sha256, allow_synthetic=allow_synthetic)
        names = {'package.json', 'review.json', 'acceptance.json', 'filter-exclusions', 'deferred-alleles',
                 *[r['path'] for r in bundle.data['artifacts'].values()]}
        if bundle.data['classification'] != classification or {p.name for p in package.iterdir()} != names:
            raise Track2RuntimeError('synthesis package classification or file closure is invalid')
        if classification == 'controlled':
            for name in names:
                if name in {'filter-exclusions','deferred-alleles'}:
                    for member in ('records.jsonl','manifest.json'):case_path(package/name/member)
                else:
                    case_path(package / name)
        mounts.append((package, '/package'))
    with output_storage(output, classification=classification, controlled_store=controlled_store,
                        kind='track2-offline-stage') as target:
        target.mkdir(parents=True, exist_ok=False, mode=0o700)
        results = target / 'results'
        results.mkdir(mode=0o700)
        container = 'track2-offline-' + uuid.uuid4().hex
        command = ['docker', 'create', '--name', container, '--pull', 'never', '--network', 'none',
            '--read-only', '--cap-drop', 'ALL', '--security-opt', 'no-new-privileges',
            '--user', f'{os.getuid()}:{os.getgid()}', '--cpus', '2', '--memory', '4g', '--pids-limit', '256',
            '--tmpfs', '/tmp:rw,nosuid,nodev,noexec,size=256m', '--env', 'PYTHONDONTWRITEBYTECODE=1',
            '--env', 'PYTHONPATH=/opt/genosapiens/src', '--workdir', '/opt/genosapiens', '--entrypoint', 'python']
        for source, destination in mounts:
            if ',' in str(source):
                raise Track2RuntimeError('runtime mount path is invalid')
            command += ['--mount', f'type=bind,src={source},dst={destination},readonly']
        command += ['--mount', f'type=bind,src={results},dst=/output', runtime['image'], '/program.py', stage]
        (target / 'command.json').write_bytes(canonical({'argv': command, 'runtime': runtime}))
        created = False
        status = {'status': 'incomplete', 'stage': stage, 'classification': classification}
        try:
            subprocess.run(command, check=True, capture_output=True, timeout=30)
            created = True
            inspection = json.loads(subprocess.check_output(['docker', 'inspect', container], timeout=30))[0]
            host = inspection['HostConfig']
            actual = {(m['Destination'], m['Source'], m['RW']) for m in inspection['Mounts']}
            expected = {(destination, str(source), False) for source, destination in mounts} | {('/output', str(results), True)}
            if (actual != expected or host['NetworkMode'] != 'none' or not host['ReadonlyRootfs']
                    or host['Privileged'] or host['CapDrop'] != ['ALL']
                    or inspection['Config']['User'] != f'{os.getuid()}:{os.getgid()}'):
                raise Track2RuntimeError('actual container isolation differs from the selected runtime')
            (target / 'isolation.json').write_bytes(canonical({'image_id': inspection['Image'],
                'network': host['NetworkMode'], 'readonly_root': host['ReadonlyRootfs'],
                'user': inspection['Config']['User'], 'mounts': inspection['Mounts'],
                'cap_drop': host['CapDrop'], 'security_options': host['SecurityOpt']}))
            with (target / 'worker.log').open('xb') as log:
                completed = subprocess.run(['docker', 'start', '--attach', container], stdout=log,
                    stderr=subprocess.STDOUT, timeout=timeout_seconds)
            state = json.loads(subprocess.check_output(['docker', 'inspect', container], timeout=30))[0]['State']
            (target / 'exit.json').write_bytes(canonical(state))
            if completed.returncode or state['ExitCode'] or state['OOMKilled']:
                raise Track2RuntimeError('offline stage failed; classified receipt retained')
            status['status'] = 'completed'
        except subprocess.TimeoutExpired:
            status['status'] = 'timed_out'
            raise Track2RuntimeError('offline stage timed out; classified receipt retained') from None
        except KeyboardInterrupt:
            status['status'] = 'cancelled'
            raise
        except (subprocess.SubprocessError, OSError):
            raise Track2RuntimeError('offline stage failed; classified receipt retained') from None
        finally:
            try:
                if created:
                    subprocess.run(['docker', 'rm', '--force', container], check=True, capture_output=True, timeout=30)
            finally:
                (target / 'execution.json').write_bytes(canonical(status))
    return status
