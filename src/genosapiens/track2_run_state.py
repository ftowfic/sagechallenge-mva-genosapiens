"""Current run inputs and dependency-bound stages, without attempt history."""
from contextlib import contextmanager
import fcntl
import os
from pathlib import Path
import shutil

from genosapiens.track2_operator_inputs import OperatorError, private_path, require, small_json
from genosapiens.track2_published_genotype import canonical, sha256

SCHEMA = 'genosapiens.track2-current-run/v1'
STAGE_SCHEMA = 'genosapiens.track2-current-stage/v1'
INPUTS = ('vcf', 'diagnosis', 'sample', 'hpo', 'concordance')
PATH_INPUTS = frozenset({'vcf', 'hpo', 'concordance'})
REFERENCED_STAGES = frozenset({'assessment','review','report'})


def atomic_record(path, value):
    """Use Track 1's durable, bounded current-document writer."""
    from genosapiens.wgs_run_manifest import atomic_current
    path = Path(path)
    require(path.resolve() == path.absolute(), 'Current run state contains a symbolic link.')
    path.parent.mkdir(mode=0o700, parents=True, exist_ok=True)
    atomic_current(path, value)


def input_values(args):
    return {name: (str(Path(value).absolute()) if name in PATH_INPUTS and value is not None else value)
        for name in INPUTS for value in (getattr(args, name),)}


def resolve_inputs(args):
    """Read and compare before writes; omitted resume inputs come from the run."""
    root = Path(args.output).absolute()
    if root.exists():
        private_path(root, directory=True)
    values = input_values(args)
    if args.resume:
        record = small_json(root/'current-run.json')
        require(record.get('schema') == SCHEMA and set(record.get('inputs', {})) == set(INPUTS),
            'The current input manifest is incompatible. Reconstruct input intake in this output directory.')
        for name, value in values.items():
            require(value is None or value == record['inputs'][name],
                'An explicit input conflicts with the current run manifest.')
        values = record['inputs']
    else:
        require(not root.exists() or not any(root.iterdir()),
            'Existing output requires --resume.')
    require(values['vcf'] is not None and values['diagnosis'] is not None,
        'Supply --vcf and --diagnosis for a new run.')
    for name, value in values.items():
        setattr(args, name, Path(value) if name in PATH_INPUTS and value is not None else value)
    return values


@contextmanager
def output_lock(args):
    """Serialize output writers; refuse conflicting inputs before lock creation."""
    resolve_inputs(args)
    root = Path(args.output).absolute()
    root.mkdir(mode=0o700, parents=True, exist_ok=True)
    path = root/'.run.lock'
    require(not path.is_symlink(), 'The output lock contains a symbolic link.')
    descriptor = os.open(path, os.O_CREAT | os.O_RDWR | os.O_NOFOLLOW, 0o600)
    try:
        try:
            fcntl.flock(descriptor, fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError:
            raise OperatorError('Another invocation holds this output lock.') from None
        # A competing writer could publish between the initial read and lock.
        if args.resume:
            resolve_inputs(args)
        else:
            require(not any(p != path for p in root.iterdir()), 'Existing output requires --resume.')
        yield RunState(root)
    finally:
        os.close(descriptor)


def files(root):
    root = Path(root)
    require(root.resolve() == root.absolute(), 'A stage output path contains a symbolic link.')
    result = {}
    from genosapiens.governance import _file_identity
    for path in sorted(root.rglob('*')):
        require(not path.is_symlink(), 'A stage output contains a symbolic link.')
        if path.is_file():
            result[path.relative_to(root).as_posix()] = _file_identity(path)[0]
        else:
            require(path.is_dir(), 'A stage output contains an unsupported file type.')
    return result


class RunState:
    """Caller holds the output lock for every mutation.

    Each stage binds its actual upstream scientific products. No global code
    hash invalidates unaffected stages. Owners choose those dependency records.
    """
    def __init__(self, root, *, direct_stage=None):
        self.root = Path(root).absolute()
        self.direct_stage = direct_stage

    def inputs(self, values, *, source_versions):
        atomic_record(self.root/'current-run.json', {
            'schema': SCHEMA, 'inputs': values, 'source_versions': source_versions})

    def stage_paths(self, stage):
        require(stage and all(c in 'abcdefghijklmnopqrstuvwxyz-' for c in stage),
            'Invalid current stage name.')
        if stage == self.direct_stage:
            return self.root, self.root/'completion.json'
        names = {'construct': 'construction', 'assessment': 'plan', 'review': 'preparation',
            'report': 'synthesis', 'concordance': 'concordance'}
        return self.root/names.get(stage, 'stages/'+stage), self.root/'completions'/(stage+'.json')

    def outputs(self, stage, target):
        rows = files(target)
        if stage == self.direct_stage:
            rows = {name:value for name,value in rows.items()
                if name not in {'completion.json','.discovery.lock','current-diagnostics.json'}}
        if stage == 'construct':
            rows = {name: value for name, value in rows.items()
                if name not in {'input-dossier.md', 'package/review.json'}}
        if stage == 'assessment':
            rows = {name: value for name, value in rows.items() if name != 'plan-review.json'}
        return rows

    def completed(self, stage, binding):
        target, path = self.stage_paths(stage)
        if not path.exists():
            return None, 'completion_missing'
        if stage == self.direct_stage or stage in REFERENCED_STAGES:
            from genosapiens.track2_scientific_completion import check_completion, SCHEMA, CONTROL_BYTES
            try:
                header = small_json(path, maximum=CONTROL_BYTES)
                if header.get('schema') != SCHEMA:
                    return None, 'completion_contract_changed'
                if header.get('binding') != binding:
                    return None, 'dependencies_changed'
                checked=check_completion(target, stage=stage, binding=binding)
                require(header==checked,'The stage completion reference changed.')
                return checked, None
            except (ValueError, OSError, KeyError, TypeError):
                return None, 'scientific_completion_invalid_or_changed'
        try:
            record = small_json(path)
        except OperatorError:
            return None, 'completion_invalid'
        if (set(record) != {'schema', 'stage', 'binding', 'outputs', 'result'}
                or record['schema'] != STAGE_SCHEMA or record['stage'] != stage):
            return None, 'completion_contract_changed'
        from genosapiens.track2_construction_dependencies import construction_binding_matches
        if not (construction_binding_matches(record['binding'], binding) if stage == 'construct'
                else record['binding'] == binding):
            return None, 'dependencies_changed'
        if not target.is_dir() or not record['outputs'] or self.outputs(stage, target) != record['outputs']:
            return None, 'output_incomplete_or_changed'
        return record, None

    def preserve_accepted_report(self, stage):
        """Preserve accepted report bytes and their reviews before replacement.

        Use the report owner's retained human plan review. Automatic technical
        acceptance cannot supply an actual human review.
        """
        target, _ = self.stage_paths(stage)
        if stage != 'report' or not target.exists():
            return
        closure = files(target)
        reviewed = any((record := small_json(target/name)).get('decision') == 'accepted'
            and record.get('reviewer_type') == 'human' for name in closure
            if Path(name).name == 'plan-review-source.json')
        if not reviewed:
            return
        identity = sha256(canonical(closure))
        destination = self.root/'accepted-reports'/identity
        if destination.exists():
            require(files(destination) == closure, 'Preserved accepted report bytes changed.')
        else:
            destination.parent.mkdir(parents=True, mode=0o700, exist_ok=True)
            pending = destination.parent/'.preserving'
            if pending.exists():
                files(pending)
                shutil.rmtree(pending)
            shutil.copytree(target, pending)
            require(files(pending) == closure, 'Accepted report preservation is incomplete.')
            # Make the retained copy durable before removing accepted source bytes.
            for name in closure:
                with (pending/name).open('rb') as stream:os.fsync(stream.fileno())
            for folder, _, _ in os.walk(pending,topdown=False):
                descriptor=os.open(folder,os.O_RDONLY|os.O_DIRECTORY|os.O_NOFOLLOW)
                try:os.fsync(descriptor)
                finally:os.close(descriptor)
            pending.rename(destination)
        descriptor=os.open(destination.parent,os.O_RDONLY|os.O_DIRECTORY|os.O_NOFOLLOW)
        try:os.fsync(descriptor)
        finally:os.close(descriptor)

    def restart(self, stage, reason):
        target, completion = self.stage_paths(stage)
        self.preserve_accepted_report(stage)
        # Unpublish before removing incomplete/stale derived output. An
        # interruption here cannot leave a valid completion for partial bytes.
        completion.unlink(missing_ok=True)
        if stage == self.direct_stage:
            # Completed scientific jobs retain their own dependency bindings.
            # A changed aggregate scope can reuse only matching completed jobs.
            target.mkdir(parents=True, mode=0o700, exist_ok=True)
            files(target)
            for path in target.iterdir():
                if stage=='figures' and path.is_dir() and path.name!='artifact-index-parts':
                    # Each figure job must re-admit its source-bound receipt.
                    # Retain these targets until the coordinator checks them.
                    continue
                if path.name in {'jobs','.discovery.lock','extractions','reuse-mappings.json','recovery-plan.json'}:
                    continue
                shutil.rmtree(path) if path.is_dir() else path.unlink()
            atomic_record(self.root/'current-diagnostics.json',
                {'stage':stage,'state':'running','restart_reason':reason})
            return target
        if target.exists():
            files(target)  # Reject links before a recursive removal.
            shutil.rmtree(target)
        target.mkdir(parents=True, mode=0o700)
        atomic_record(self.root/'current-diagnostics.json',
            {'stage': stage, 'state': 'running', 'restart_reason': reason})
        return target

    def publish(self, stage, *, binding, result):
        target, completion = self.stage_paths(stage)
        if stage == self.direct_stage or stage in REFERENCED_STAGES:
            from genosapiens.track2_scientific_completion import publish_completion
            record = publish_completion(target, stage=stage, binding=binding, result=result)
            if completion != target/'completion.json':atomic_record(completion,record)
            atomic_record(self.root/'current-diagnostics.json', {'stage': stage, 'state': 'complete'})
            return record
        outputs = self.outputs(stage, target)
        require(outputs, 'A stage requires output evidence before completion.')
        record = {'schema': STAGE_SCHEMA, 'stage': stage, 'binding': binding,
            'outputs': outputs, 'result': result}
        atomic_record(completion, record)
        atomic_record(self.root/'current-diagnostics.json', {'stage': stage, 'state': 'complete'})
        return record
