#!/usr/bin/env python3
"""Host-owned Track 2 Docker orchestration. Standard library only, like Track 1."""
import argparse
import csv
import fcntl
import io
import json
import math
import os
from pathlib import Path
import re
import shutil
import signal
import stat
import subprocess
import sys
import tempfile
import time
import uuid

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))
WORKER_INPUTS = {
    'vcf': Path('/track2-inputs/source.vcf'),
    'hpo': Path('/track2-inputs/phenotypes.txt'),
    'concordance': Path('/track2-inputs/concordance.json'),
    'review': Path('/track2-inputs/review.tsv'),
    'mechanism_review': Path('/track2-inputs/mechanism-review.tsv'),
}


def mount(source, target=None, readonly=True):
    stream = io.StringIO()
    csv.writer(stream, lineterminator='').writerow(['type=bind', 'source='+str(source),
        'target='+str(source if target is None else target)] + (['readonly'] if readonly else []))
    return stream.getvalue()


def private(path, *, directory=False):
    path = Path(path).absolute()
    if path.resolve() != path or any(p.is_symlink() for p in (path, *path.parents)):
        raise ValueError('Use real paths without symbolic links.')
    info = path.stat()
    if info.st_uid != os.getuid() or stat.S_IMODE(info.st_mode) & 0o077:
        raise ValueError('Private paths must belong to you and deny access to other users.')
    if not (stat.S_ISDIR(info.st_mode) if directory else stat.S_ISREG(info.st_mode)):
        raise ValueError('The selected private path has the wrong type.')
    return path


def read_json(path):
    with Path(path).open('rb') as stream:
        raw = stream.read(32*1024*1024+1)
    if len(raw) > 32*1024*1024:
        raise ValueError('Local metadata exceeds its supported size.')
    return json.loads(raw)


def atomic(path, data):
    path = Path(path)
    if path.is_symlink():
        raise ValueError('Current state must not be a symbolic link.')
    fd, name = tempfile.mkstemp(prefix='.current-', dir=path.parent)
    try:
        with os.fdopen(fd, 'w') as stream:
            json.dump(data, stream, sort_keys=True, separators=(',', ':'))
            stream.flush(); os.fsync(stream.fileno())
        os.replace(name, path)
    finally:
        Path(name).unlink(missing_ok=True)


def positive_hours(value):
    hours = float(value)
    if not math.isfinite(hours) or hours <= 0 or not math.isfinite(hours*3600):
        raise argparse.ArgumentTypeError('--max-hours requires positive finite hours.')
    return hours


def parser():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('command', choices=('run', 'setup', 'install-public-evidence','prepare-discovery-resources'))
    p.add_argument('--output', type=Path)
    for name in ('vcf', 'hpo', 'concordance', 'review', 'mechanism-review', 'source'):
        p.add_argument('--'+name, type=Path)
    p.add_argument('--rules', type=Path,
        help='Official rules file required for first setup. Valid registration reuses its retained copy.')
    p.add_argument('--diagnosis'); p.add_argument('--sample')
    p.add_argument('--resume', action='store_true')
    p.add_argument('--require-review', action='store_true')
    p.add_argument('--max-hours', type=positive_hours,
        help='Positive finite hours for the whole invocation. Omission means no overall deadline.')
    p.add_argument('--memory-gib', type=int, default=4)
    p.add_argument('--lincs-workers', type=int, choices=(1,2), default=2)
    p.add_argument('--figure-workers', type=int, choices=(1,2), default=2)
    p.add_argument('--reference-root', type=Path,
        help='Host reference directory. Fresh default: /mnt/reference. Resume restores the saved root.')
    p.add_argument('--settings', type=Path, help='Existing authenticated operator registration.')
    return p


class Docker:
    def __init__(self, image, memory):
        self.image, self.memory, self.name = image, memory, None

    def base(self):
        return ['--init', '--network', 'none', '--read-only', '--cpus', '2',
            '--memory', str(self.memory)+'g', '--memory-swap', str(self.memory)+'g', '--pids-limit', '256',
            '--user', f'{os.getuid()}:{os.getgid()}', '--cap-drop', 'ALL',
            '--security-opt', 'no-new-privileges', '--tmpfs', '/tmp:rw,nosuid,nodev,noexec,size=256m,mode=1777',
            '--env', 'PYTHONDONTWRITEBYTECODE=1', '--entrypoint', 'python3']

    def cleanup(self):
        if self.name:
            subprocess.run(['docker', 'rm', '--force', self.name], stdout=subprocess.DEVNULL,
                stderr=subprocess.DEVNULL, timeout=30, check=False)
            self.name = None

    def run(self, options, arguments, *, log=None, capture=False):
        self.name = 'track2-'+str(os.getuid())+'-'+uuid.uuid4().hex
        command = ['docker', 'run', '--pull', 'never', '--name', self.name,
            *self.base(), *options, self.image, *arguments]
        if log is not None:
            atomic(log.with_suffix('.command.json'), command)
        try:
            if capture:
                result = subprocess.run(command, text=True, stdout=subprocess.PIPE, check=True)
                return result.stdout
            with (log.open('w') if log is not None else open(os.devnull, 'w')) as stream:
                process = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
                try:
                    for line in process.stdout:
                        print(line, end='', flush=True); stream.write(line); stream.flush()
                    code=process.wait()
                    inspected=subprocess.run(['docker','inspect','--format','{{json .State}}',self.name],
                        text=True,capture_output=True,check=False,timeout=30)
                    state=json.loads(inspected.stdout) if inspected.returncode==0 else {'ExitCode':code}
                    if log is not None: atomic(log.with_suffix('.exit.json'),state)
                    if state.get('OOMKilled'):
                        print('The worker exceeded its memory limit. Resume with a suitable --memory-gib value.',flush=True)
                    elif code:
                        print('The worker returned status '+str(code)+'.',flush=True)
                    return code
                finally:
                    if process.poll() is None:
                        process.terminate()
                        process.wait(timeout=30)
        finally:
            self.cleanup()


def preparation_options(repository, reference, store, index_store):
    """Same public-only mount pattern as ensure_track1_installation.py."""
    options = ['--env', 'PYTHONPATH='+str(repository/'src'), '--workdir', str(repository)]
    for name in ('src', 'config', 'schemas', 'scripts', 'track2'):
        options += ['--mount', mount(repository/name)]
    for name in ('releases', 'manifests', 'derived', 'downloads'):
        if not (reference/name).is_dir() or (reference/name).is_symlink():
            raise ValueError('Installed public resource directory is missing: '+str(reference/name))
        options += ['--mount', mount(reference/name)]
    for name in ('manifests', 'derived'):
        path = reference/name/'common-allele-public-indexes'
        path.mkdir(exist_ok=True)
        if path.is_symlink():
            raise ValueError('Public index output must not be a symbolic link.')
        options += ['--mount', mount(path, readonly=False)]
    for path in (store, index_store):
        options += ['--mount', mount(path, readonly=False)]
    return options


def prepare(docker, reference):
    from scripts.prepare_track2_discovery import ensure_discovery_resources
    discovery=ensure_discovery_resources(ROOT,reference,memory_gib=docker.memory)
    store, index_store = ROOT/'.track2-packages', ROOT/'.track1-packages'
    for path in (store, index_store):
        path.mkdir(mode=0o700, exist_ok=True); private(path, directory=True)
    fd=os.open(store/'.setup.lock',os.O_CREAT|os.O_RDWR|os.O_NOFOLLOW,0o600)
    with os.fdopen(fd,'a') as lock:
        fcntl.flock(lock, fcntl.LOCK_EX)
        print('Prepare or reuse the public index and local Track 2 package.', flush=True)
        record = json.loads(docker.run(preparation_options(ROOT, reference, store, index_store),
            ['-m', 'genosapiens.track2_setup', '--repository', str(ROOT), '--reference', str(reference),
             '--store', str(store), '--index-store', str(index_store),
             '--discovery-resources',str(reference/discovery['path'])], capture=True))
        if not re.fullmatch('package-[0-9a-f]{64}', record['package_name']):
            raise ValueError('Prepared package location is invalid.')
        return store/record['package_name'], record['package_sha256']


def setup_rules(args):
    """Select rules without creating registration or public preparation state."""
    state = args.reference_root/'controlled-work'/('.track2-operator-'+str(os.getuid()))
    settings_path = args.settings or state/'command/track2-settings.json'
    if settings_path.exists():
        settings = read_json(private(settings_path))
        rules = Path(settings['rules']) if args.rules is None else args.rules.absolute()
        return state, settings_path, rules
    if args.settings is not None:
        raise ValueError('The supplied registration is missing. Preserve its records for repair.')
    if state.exists():
        private(state, directory=True)
        if any(p.name != '.setup.lock' for p in state.iterdir()):
            raise ValueError('Setup is incomplete. Preserve its records for repair.')
    if args.rules is None:
        raise ValueError('Supply --rules PATH with the official rules from '
            'https://huggingface.co/datasets/SageBio/mva-hackathon-2026-data. '
            'Public preparation has not started.')
    rules = args.rules.absolute()
    if rules.resolve() != rules or not rules.is_file():
        raise ValueError('Supply a readable official rules file without symbolic links using --rules PATH.')
    with rules.open('rb') as stream:
        content = stream.read(1024*1024+1)
    if not content.strip() or len(content) > 1024*1024:
        raise ValueError('The official rules file is empty or exceeds 1 MiB.')
    return state, None, rules


def check_setup_registration(docker, settings_path, reference, rules):
    """Use the existing authorization checks in a read-only container."""
    settings = read_json(private(settings_path))
    options = ['--mount', mount(ROOT/'src', '/setup-code'), '--env', 'PYTHONPATH=/setup-code',
        '--mount', mount(settings_path)]
    for key in ('governance', 'verification_key', 'rules', 'incident_root'):
        options += ['--mount', mount(settings[key])]
    options += ['--mount', mount(rules, '/supplied-rules.txt')]
    docker.run(options, ['-c',
        'import sys; from genosapiens.track2_operator_setup import validate_registration; '
        'validate_registration(sys.argv[1], sys.argv[2], "/supplied-rules.txt")',
        str(settings_path), str(reference)], capture=True)


def resolve(args):
    if args.output is None:
        raise ValueError('--output is required.')
    root = args.output.absolute()
    if root.resolve() != root:
        raise ValueError('Output must not contain symbolic links.')
    names = ('vcf', 'hpo', 'concordance', 'sample', 'diagnosis')
    values = {name: getattr(args, name) for name in names}
    for name in ('vcf', 'hpo', 'concordance'):
        if values[name] is not None:
            values[name] = str(Path(values[name]).absolute())
    if args.resume:
        private(root, directory=True)
        record = read_json(private(root/'current-run.json'))
        if record.get('schema') != 'genosapiens.track2-current-run/v1':
            raise ValueError('The current input manifest is incompatible. Supply the original inputs to rebuild intake.')
        for name, value in values.items():
            if value is not None and value != record['inputs'][name]:
                raise ValueError('An explicit input conflicts with the current run manifest.')
        values = record['inputs']
    elif root.exists() and any(root.iterdir()):
        raise ValueError('Existing output requires --resume.')
    if not values['vcf'] or not values['diagnosis']:
        raise ValueError('Supply --vcf and --diagnosis for a new run.')
    for name, value in values.items():
        setattr(args, name, Path(value) if name in ('vcf', 'hpo', 'concordance') and value is not None else value)
    for name in ('vcf', 'hpo', 'concordance', 'review', 'mechanism_review'):
        path = getattr(args, name)
        if path is not None:
            path=private(path)
            setattr(args,name,path)
            if name in ('vcf', 'hpo', 'concordance') and path.absolute().is_relative_to(root):
                raise ValueError('Keep source inputs outside the run output directory.')
    if args.resume:
        for name, identity in record.get('source_identities', {}).items():
            if source_identity(getattr(args,name)) != identity:
                raise ValueError('A source input changed since this run started.')
    args.output = root
    return values


def resolve_reference(args):
    from scripts.portable_paths import reference_selection, track2_legacy_reference
    saved = None
    if args.command == 'run' and args.resume:
        if args.output is None:
            raise ValueError('--output is required.')
        record = read_json(private(args.output/'current-run.json'))
        saved = record.get('reference_root') or track2_legacy_reference(record, ROOT)
        if saved is None:
            if args.reference_root is None:
                raise ValueError('This run has no saved reference root. Supply --reference-root PATH '
                    'and restore the original authenticated package for compatibility checks.')
            raise ValueError('Restore the original authenticated package in .track2-packages '
                'to check the explicit reference root before preparation.')
    args.reference_root = reference_selection(args.reference_root, saved)


def source_identity(path):
    if path is None: return None
    row=Path(path).stat()
    return [row.st_dev,row.st_ino,row.st_size,row.st_mtime_ns,row.st_ctime_ns]


def discovery_mounts(package,reference,*,matrix=False):
    descriptor=read_json(Path(package)/'discovery-resources.json')
    names=[str(Path(descriptor['catalog']['path']).parent)]
    if descriptor.get('schema')=='genosapiens.track2-discovery-resources/v2':
        names.append(str(Path(descriptor['tissue_companion']['path']).parent))
    if matrix:names += [descriptor['matrix']['path'],descriptor['matrix']['producer_receipt_path']]
    options=[]
    for name in names:
        path=Path(reference)/name
        if Path(name).is_absolute() or '..' in Path(name).parts or path.resolve()!=path.absolute():
            raise ValueError('A discovery resource path is invalid.')
        options += ['--mount',mount(path,Path('/discovery-resources')/name)]
    return options


def run_processing_worker(docker, options, arguments, root, *, log=None):
    """Keep processing scratch private and remove it after worker cleanup."""
    root=private(root,directory=True)
    with tempfile.TemporaryDirectory(prefix='.track2-tmp-',dir=root) as scratch:
        options=[*options, '--mount', mount(scratch,'/track2-scratch',False),
            '--env','TMPDIR=/track2-scratch','--env','SQLITE_TMPDIR=/track2-scratch']
        return docker.run(options,arguments,log=log)


def run_discovery_stage(docker,package,reference,plan,output,*,log=None,invocation_id=None,lincs_workers=2):
    """The discovery worker receives only its frozen plan, public assets, and output."""
    output.mkdir(mode=0o700,parents=True,exist_ok=True)
    output=private(output,directory=True)
    retained=output.with_name(output.name+'-retained')
    retained.mkdir(mode=0o700,exist_ok=True)
    retained=private(retained,directory=True)
    # Create the nested mountpoint as the operator, before Docker starts.
    (retained/'current').mkdir(mode=0o700,exist_ok=True)
    private(retained/'current',directory=True)
    # Discovery preserves old revisions beside its current directory. Mount a
    # dedicated persistent parent, not the run root or temporary scratch. Keep
    # the existing host output in place and outside the archive storage.
    worker_output='/discovery-work/current'
    options=['--mount',mount(package,'/software'),'--env','PYTHONPATH=/software/src',
        '--mount',mount(plan,'/plan'),
        '--mount',mount(retained,'/discovery-work',False),
        '--mount',mount(output,worker_output,False),
        *discovery_mounts(package,reference,matrix=True)]
    invocation_id=invocation_id or uuid.uuid4().hex
    return run_processing_worker(docker,options,['-m','genosapiens','track2','discover','/plan/plan.json',
        '/discovery-resources',worker_output,'--invocation-id',invocation_id,
        '--lincs-workers',str(lincs_workers)],output.parent,log=log)


def check_discovery_result(package,output,plan,invocation_id,runtime):
    """Check current worker authority and bytes using only standard-library code."""
    sys.path.insert(0,str(package/'src'))
    try:
        from genosapiens.track2_discovery_outcome import checked_outcome
        return checked_outcome(output,plan_path=plan/'plan.json',invocation_id=invocation_id,runtime=runtime)
    finally:sys.path.pop(0)


def record_incomplete_finalization(package, output, code):
    """Record a failed worker exit when it left no terminal diagnostic."""
    sys.path.insert(0,str(package/'src'))
    try:
        try: diagnostic=read_json(output/'current-diagnostics.json')
        except (OSError,ValueError): diagnostic={}
        if (isinstance(diagnostic,dict) and diagnostic.get('stage')=='finalize'
                and diagnostic.get('state') in {'failed','interrupted'}):
            return
        from genosapiens.track2_report_finalization import record_finalization_failure
        record_finalization_failure(output,subprocess.CalledProcessError(code,['track2','finalize']))
    except BaseException:
        print('Finalization failed. Its diagnostic record is unavailable.',file=sys.stderr)
    finally:sys.path.pop(0)


def run_discovery_figures(docker,package,discovery,output,*,log=None,figure_workers=2):
    if type(figure_workers) is not int or figure_workers not in (1,2):
        raise ValueError('Figure workers must be 1 or 2.')
    output.mkdir(mode=0o700,parents=True,exist_ok=True)
    attempts=output.parent/(output.name+'-attempts')
    attempts.mkdir(mode=0o700,parents=True,exist_ok=True)
    # Publication renames between sibling output and attempt directories.
    # Separate bind mounts reject that rename even on the same host filesystem.
    worker_output=str(Path('/figure-work')/output.name)
    options=['--mount',mount(package,'/software'),'--env','PYTHONPATH=/software/src',
        '--mount',mount(discovery,'/results'),
        '--mount',mount(output.parent,'/figure-work',False)]
    return run_processing_worker(docker,options,['-c','import sys; from genosapiens.track2_lincs_report import render_native_jobs; render_native_jobs("/results",sys.argv[1],figure_workers='+str(figure_workers)+')',worker_output],output.parent,log=log)


def run_report_overviews(docker,package,root,*,log=None):
    output=root/'synthesis/overviews'
    output.mkdir(mode=0o700,parents=True,exist_ok=True)
    options=['--mount',mount(package,'/software'),'--env','PYTHONPATH=/software/src',
        '--mount',mount(root/'synthesis/controlled-work/report/presentation.json','/presentation.json'),
        '--mount',mount(output,'/output',False)]
    return run_processing_worker(docker,options,['-c','from genosapiens.track2_report_overviews import render_overviews; from genosapiens.track2_scientific_reader import read_scientific; render_overviews("/output",read_scientific("/presentation.json"))'],root,log=log)


def prepare_annotation_resume(package, root, image):
    """Keep the immutable installed package free of interpreter cache files."""
    previous = sys.dont_write_bytecode
    sys.dont_write_bytecode = True
    sys.path.insert(0, str(package/'src'))
    try:
        from genosapiens.track2_annotation_recovery import prepare_recovery
        prepare_recovery(root, package=package, package_store=ROOT/'.track2-packages', image=image)
    finally:
        sys.path.pop(0)
        sys.dont_write_bytecode = previous


def execute(args, docker, package, digest, settings_path):
    settings = read_json(private(settings_path))
    if settings['reference_root']!=str(args.reference_root):
        raise ValueError('The registration uses a different reference root.')
    root = args.output
    native=(package/'discovery-resources.json').is_file()
    registration = []
    for key in ('governance', 'verification_key', 'rules', 'incident_root'):
        registration += ['--mount', mount(settings[key])]
    base = ['--mount', mount(package), '--mount', mount(settings_path), '--env', 'PYTHONPATH='+str(package/'src'),
        '--env', 'GENOSAPIENS_TRACK2_PACKAGE='+str(package), '--env', 'GENOSAPIENS_TRACK2_PACKAGE_SHA256='+digest,
        '--env', 'GENOSAPIENS_TRACK2_IMAGE='+docker.image, '--env', 'GENOSAPIENS_TRACK2_SETTINGS='+str(settings_path)]
    # Host paths can cross symlinks that exist only in the image (such as
    # /data). Keep worker input locations independent of the host layout.
    options = ['--output', str(root), '--vcf', str(WORKER_INPUTS['vcf']), '--diagnosis', args.diagnosis]
    for key in ('sample', 'hpo', 'concordance', 'review', 'mechanism_review'):
        if getattr(args, key) is not None:
            options += ['--'+key.replace('_', '-'), str(WORKER_INPUTS[key] if key in WORKER_INPUTS else getattr(args, key))]
    if args.resume: options += ['--resume']
    if args.require_review: options += ['--require-review']
    diagnostics = root/'current-diagnostics'
    if diagnostics.exists():
        private(diagnostics, directory=True)
        if any(p.is_symlink() for p in diagnostics.rglob('*')):
            raise ValueError('Current diagnostics contain a symbolic link.')
        shutil.rmtree(diagnostics)
    diagnostics.mkdir(mode=0o700)
    inventory = Path(settings['inventory'])
    private(inventory.parent, directory=True)
    fd=os.open(inventory.parent/'.track2.lock',os.O_CREAT|os.O_RDWR|os.O_NOFOLLOW,0o600)
    with os.fdopen(fd,'a') as lock:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        authorized = False
        primary_code = None
        finalized = False
        active_phase = None
        def report_location():
            report = root/'synthesis/controlled-work/report/report.md'
            if report.is_file():
                print(('Report: ' if finalized else 'Unfinalized report: ')+str(report), flush=True)
        def finish(code):
            nonlocal primary_code
            primary_code=code
            if active_phase in {'finalize','prepare-report'} and code not in (0,2):
                record_incomplete_finalization(package,root,code)
            report_location()
            return code
        def phase(name):
            nonlocal active_phase
            active_phase = name
            mounts = [*base, '--env', 'GENOSAPIENS_TRACK2_PHASE='+name,
                '--env', 'GENOSAPIENS_TRACK2_OUTCOME='+str(diagnostics/(name+'-outcome.json'))]
            if name == 'synthesize':
                mounts += ['--env','GENOSAPIENS_TRACK2_REPORT_MEMORY_BYTES='+str(args.memory_gib*1024**3)]
                mounts += ['--mount', mount(root/'preparation/staged/package', '/package'),
                    '--mount', mount(root/'plan', '/plan'),
                    '--mount', mount(root/'preparation/snapshot.json', '/inputs/snapshot.json'),
                    '--mount', mount(root/'current-routing.json', '/inputs/production.json'),
                    '--mount', mount(root/'synthesis', '/output', False),
                    '--env', 'GENOSAPIENS_TRACK2_OUTCOME=/output/terminal-outcome.json']
                if native:
                    mounts += [*discovery_mounts(package,args.reference_root),
                        '--mount',mount(root/'discovery','/discovery-results'),
                        '--mount',mount(root/'discovery-figures','/discovery-figures')]
                elif (root/'study-science').is_dir():
                    mounts += ['--mount',mount(root/'study-science','/study-science')]
            else:
                mounts += [*registration, '--mount', mount(root, readonly=False)]
                if name != 'construct':
                    mounts += ['--mount', mount(inventory.parent, readonly=False)]
                if name == 'construct':
                    resources = read_json(package/'resources.json')
                    for path in resources['public_mounts']:
                        mounts += ['--mount', mount(path)]
                    for key in ('vcf', 'hpo', 'concordance'):
                        if getattr(args, key) is not None:
                            mounts += ['--mount', mount(getattr(args, key), WORKER_INPUTS[key]),
                                '--env', 'GENOSAPIENS_TRACK2_HOST_'+key.upper()+'='+str(getattr(args, key))]
                if name == 'review':
                    if native:mounts += discovery_mounts(package,args.reference_root)
                    for key in ('review', 'mechanism_review'):
                        if getattr(args, key) is not None:
                            mounts += ['--mount', mount(getattr(args, key), WORKER_INPUTS[key]),
                                '--env', 'GENOSAPIENS_TRACK2_HOST_'+key.upper()+'='+str(getattr(args, key))]
            started=time.monotonic()
            print('Stage: '+name+'.', flush=True)
            arguments=['-m','genosapiens','track2','run',*options]
            log=None if name=='inventory' else diagnostics/(name+'.log')
            if name=='synthesize':
                reporting=Docker(read_json(package/'config/track2/discovery-runtime.json')['image'],args.memory_gib)
                try:result=run_processing_worker(reporting,mounts,arguments,root,log=log)
                finally:reporting.cleanup()
            elif name in {'construct','review','prepare-report','complete-report','finalize'}:
                result=run_processing_worker(docker,mounts,arguments,root,log=log)
            else:
                result=docker.run(mounts,arguments,log=log)
            print(f'Stage {name} took {time.monotonic()-started:.1f} seconds.',flush=True)
            return result
        try:
            code = phase('authorize')
            if code: return finish(code)
            authorized = True
            for name in ('construct', 'review'):
                code = phase(name)
                if code: return finish(code)
            if native:
                worker=Docker(read_json(package/'package.json')['images']['discovery'],args.memory_gib)
                try:
                    invocation_id=uuid.uuid4().hex
                    code=run_discovery_stage(worker,package,args.reference_root,root/'plan',root/'discovery',
                        log=diagnostics/'discovery.log',invocation_id=invocation_id,lincs_workers=args.lincs_workers)
                    if code in (0,2):
                        try:
                            check_discovery_result(package,root/'discovery',root/'plan',invocation_id,worker.image)
                        except (ValueError,OSError,KeyError,TypeError):
                            print('Discovery has no checked current result. Dependent work stopped.',file=sys.stderr)
                            return finish(3)
                        figure_code=run_discovery_figures(worker,package,root/'discovery',root/'discovery-figures',
                            log=diagnostics/'discovery-figures.log',figure_workers=args.figure_workers)
                        if figure_code:return finish(figure_code)
                finally:worker.cleanup()
                if code not in (0,2):return finish(code)
                code=phase('prepare-report')
                if code:return finish(code)
            route = read_json(root/'current-routing.json')
            if not route['report_complete']:
                if not route.get('synthesis_complete'):
                    code = phase('synthesize')
                    if code not in (0, 2): return finish(code)
                code = phase('finalize')
            else:
                code = route['finalization_exit_code']
            finalized = code in (0,2) and read_json(root/'current-routing.json').get('report_complete', False)
            return finish(code)
        except BaseException as error:
            if active_phase in {'finalize','prepare-report'}:
                error._track2_finalization_diagnosed = True
                sys.path.insert(0,str(package/'src'))
                try:
                    from genosapiens.track2_report_finalization import record_finalization_failure
                    record_finalization_failure(root,error)
                except BaseException:
                    print('Finalization failed. Its diagnostic record is unavailable.',file=sys.stderr)
                finally:sys.path.pop(0)
            report_location()
            raise
        finally:
            if authorized:
                primary_failed=sys.exc_info()[0] is not None or primary_code not in (None,0)
                try:
                    result = phase('inventory')
                    if result:raise RuntimeError('The controlled output inventory update failed.')
                except BaseException:
                    print('The controlled output inventory update failed.',file=sys.stderr)
                    if not primary_failed:raise


def main(argv=None):
    os.umask(0o077)
    arguments=sys.argv[1:] if argv is None else argv
    if arguments and arguments[0]=='report':
        from scripts.report_track2_run import main as report_main
        return report_main(arguments[1:])
    args = parser().parse_args(argv)
    if args.command=='run':
        from scripts.report_track2_run import historical_resume
        revised=historical_resume(args)
        if revised is not None:return revised
    if args.command=='prepare-discovery-resources' and any(getattr(args,key) is not None for key in
            ('vcf','hpo','concordance','review','mechanism_review','output','sample','diagnosis','settings')):
        raise ValueError('Public discovery preparation accepts no case inputs or operator registration.')
    resolve_reference(args)
    if args.settings is not None: args.settings=args.settings.absolute()
    registration = setup_rules(args) if args.command == 'setup' else None
    if os.getuid() == 0 or os.getgid() == 0:
        raise ValueError('Run as your ordinary non-root account.')
    if args.memory_gib < 1 or args.memory_gib > 65536:
        raise ValueError('--memory-gib requires an integer from 1 through 65536.')
    image = read_json(ROOT/'config/track2/diagnosis-runtime.json')['image']
    if not re.fullmatch(r'genosapiens/[a-z0-9-]+@sha256:[0-9a-f]{64}', image):
        raise ValueError('The project runtime pin is invalid.')
    docker = Docker(image, args.memory_gib)
    def stop(signum, frame):
        docker.cleanup()
        raise TimeoutError() if signum == signal.SIGALRM else KeyboardInterrupt()
    for sig in (signal.SIGINT, signal.SIGTERM, signal.SIGALRM): signal.signal(sig, stop)
    if args.max_hours is not None: signal.setitimer(signal.ITIMER_REAL, args.max_hours*3600)
    subprocess.run(['docker', 'image', 'inspect', image], stdout=subprocess.DEVNULL, check=True)
    if subprocess.check_output(['docker', 'info', '--format', '{{.CgroupVersion}}'], text=True).strip() != '2':
        raise ValueError('Docker with cgroups v2 is required.')
    if args.command == 'install-public-evidence':
        if args.source is None:
            raise ValueError('Supply --source with the frozen public evidence export.')
        descriptor=read_json(ROOT/'config/track2/installed-public-evidence.json')
        destination=args.reference_root/descriptor['path']
        source=args.source.absolute()
        # The public consumer checks classification, source review and every
        # retained binding before the host copies this public resource.
        code='''import sys
from pathlib import Path
from genosapiens.track2_public_matching import load_public_evidence
r=load_public_evidence(sys.argv[1])
if r['sha256']!=sys.argv[2]: raise ValueError('Public evidence version differs from the installed descriptor.')
'''
        docker.run(['--mount', mount(ROOT/'src', '/code'), '--mount', mount(source),
            '--env', 'PYTHONPATH=/code'], ['-c', code, str(source), descriptor['evidence_sha256']], capture=True)
        if not destination.exists():
            destination.parent.mkdir(parents=True, exist_ok=True)
            pending=destination.parent/'.installing'
            if pending.exists(): shutil.rmtree(pending)
            shutil.copytree(source, pending); pending.rename(destination)
        print('The declared public evidence export is installed.')
        return 0
    if args.command == 'prepare-discovery-resources':
        from scripts.prepare_track2_discovery import ensure_discovery_resources
        print(json.dumps(ensure_discovery_resources(ROOT,args.reference_root,memory_gib=args.memory_gib)))
        return 0
    discovery_image=read_json(ROOT/'config/track2/discovery-runtime.json')['image']
    subprocess.run(['docker','image','inspect',discovery_image],stdout=subprocess.DEVNULL,check=True)
    if args.command == 'setup':
        state, settings_path, rules = registration
        if settings_path is not None:
            check_setup_registration(docker, settings_path, args.reference_root, rules)
        package, digest = prepare(docker, args.reference_root)
        if settings_path is not None:
            print('Reusing the existing registration and incident state.')
            return 0
        state.mkdir(mode=0o700,parents=True,exist_ok=True); private(state,directory=True)
        fd=os.open(state/'.setup.lock',os.O_CREAT|os.O_RDWR|os.O_NOFOLLOW,0o600)
        with os.fdopen(fd,'a') as lock:
            fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
            docker.name='track2-'+str(os.getuid())+'-'+uuid.uuid4().hex
            command=['docker','run','--rm','--pull','never','--name',docker.name,'--interactive','--tty',*docker.base(),
                '--mount',mount(package),'--mount',mount(rules,'/setup-rules.txt'),
                '--mount',mount(state,readonly=False),'--env','PYTHONPATH='+str(package/'src'),docker.image,
                '-m','genosapiens.track2_operator_setup','--package',str(package),'--rules','/setup-rules.txt',
                '--reference-root',str(args.reference_root),'--state',str(state),'--image',docker.image]
            try:
                return subprocess.run(command,check=False).returncode
            finally:
                docker.cleanup()
    values = resolve(args)
    root = args.output
    settings = args.settings or args.reference_root/'controlled-work'/('.track2-operator-'+str(os.getuid()))/'command/track2-settings.json'
    if not settings.is_file():
        raise ValueError('Authenticated operator registration is missing. Run ./track2 setup.')
    if read_json(private(settings))['reference_root'] != str(args.reference_root):
        raise ValueError('The registration uses a different reference root.')
    # As in Track 1, public preparation precedes the case workspace. An
    # interrupted public build therefore leaves no half-created case manifest.
    package, digest = prepare(docker, args.reference_root)
    root.mkdir(mode=0o700, parents=True, exist_ok=True); private(root, directory=True)
    fd=os.open(root/'.run.lock',os.O_CREAT|os.O_RDWR|os.O_NOFOLLOW,0o600)
    with os.fdopen(fd,'a') as lock:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        if args.resume: resolve(args)
        if not args.resume and any(p.name!='.run.lock' for p in root.iterdir()):
            raise ValueError('Existing output requires --resume.')
        if not args.resume:
            atomic(root/'current-run.json', {'schema': 'genosapiens.track2-current-run/v1',
                'inputs': values, 'reference_root': str(args.reference_root),
                'source_versions': {'package_sha256': digest, 'image': image},
                'source_identities':{name:source_identity(getattr(args,name))
                    for name in ('vcf','hpo','concordance')}})
        else:
            prepare_annotation_resume(package, root, image)
            current=read_json(root/'current-run.json')
            current['reference_root'] = str(args.reference_root)
            current['source_versions']={'package_sha256':digest,'image':image}
            atomic(root/'current-run.json',current)
        return execute(args, docker, package, digest, settings)


if __name__ == '__main__':
    try:
        raise SystemExit(main())
    except TimeoutError:
        print('The invocation time limit expired. Completed stages remain available.', file=sys.stderr)
        raise SystemExit(124)
    except KeyboardInterrupt:
        print('Execution cancelled. Completed stages remain available.', file=sys.stderr)
        raise SystemExit(130)
    except (OSError, ValueError, RuntimeError, subprocess.SubprocessError) as error:
        print('Finalization failed. Read the private current diagnostics.'
            if getattr(error,'_track2_finalization_diagnosed',False) else str(error), file=sys.stderr)
        raise SystemExit(1)
