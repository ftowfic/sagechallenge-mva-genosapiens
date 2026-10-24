"""Current discovery authority and safe operational diagnostics. Standard library only."""
import errno
import json
import os
from pathlib import Path
import subprocess
import sqlite3
import sys

from genosapiens.track2_operator_inputs import require, small_json
from genosapiens.track2_run_state import atomic_record
from genosapiens.track2_scientific_storage import file_identity, json_identity, ScientificCapacityError
from genosapiens.track2_scientific_completion import reference, check_reference, CONTROL_BYTES

SCHEMA = 'genosapiens.track2-discovery-outcome/v1'


def safe_diagnostic(stage, role, error):
    sqlite_status=getattr(error,'sqlite_errorcode',None) if isinstance(error,sqlite3.Error) else None
    sqlite_category=(sqlite_status & 255) if type(sqlite_status) is int else None
    category = ('interrupted' if isinstance(error, (KeyboardInterrupt, SystemExit)) else
        'memory_capacity' if isinstance(error, MemoryError) or sqlite_category==sqlite3.SQLITE_NOMEM else
        'storage_capacity' if isinstance(error,ScientificCapacityError) or sqlite_category==sqlite3.SQLITE_FULL or isinstance(error, OSError) and error.errno in {errno.ENOSPC, errno.EDQUOT, errno.EFBIG} else
        'storage_io' if isinstance(error, OSError) or sqlite_category==sqlite3.SQLITE_IOERR else
        'worker_execution' if isinstance(error, subprocess.CalledProcessError) else
        'artifact_integrity' if isinstance(error, (ValueError, KeyError, TypeError)) else 'execution_failure')
    locations = []
    frame = error.__traceback__
    while frame is not None:
        module = frame.tb_frame.f_globals.get('__name__', '')
        if module.startswith('genosapiens.'):
            locations.append(dict(module=module, function=frame.tb_frame.f_code.co_name, line=frame.tb_lineno))
        frame = frame.tb_next
    return dict(stage=stage, artifact_role=role, category=category, locations=locations[-8:],
        errno=error.errno if isinstance(error, OSError) and type(error.errno) is int else None,
        worker_status=error.returncode if isinstance(error,subprocess.CalledProcessError) else None,
        sqlite_status=sqlite_status)


def record_failure(root, *, stage, role, error, invocation_id=None):
    """Caller owns the lock. Secondary diagnostic failure cannot replace the error."""
    diagnostic = safe_diagnostic(stage, role, error)
    state = 'interrupted' if diagnostic['category'] == 'interrupted' else 'failed'
    for name, record in [('current-diagnostics.json', dict(state=state, **diagnostic)),
            ('current-outcome.json', dict(schema=SCHEMA, state=state, invocation_id=invocation_id,
                diagnostic=diagnostic))]:
        try:
            atomic_record(Path(root)/name, record)
        except BaseException as secondary:
            print(json.dumps(dict(primary=diagnostic,
                diagnostic_write=safe_diagnostic(stage, 'diagnostic-record', secondary)), sort_keys=True), file=sys.stderr)
    print(json.dumps(diagnostic, sort_keys=True), file=sys.stderr)


def invalidate_outcome(root):
    path = Path(root)/'current-outcome.json'
    require(not path.is_symlink(), 'The discovery outcome contains a symbolic link.')
    path.unlink(missing_ok=True)
    descriptor = os.open(path.parent, os.O_RDONLY | os.O_DIRECTORY | os.O_NOFOLLOW)
    try: os.fsync(descriptor)
    finally: os.close(descriptor)


def publish_outcome(root, plan, binding, *, complete, invocation_id=None):
    from genosapiens.track2_discovery import load_completed_discovery
    from genosapiens.track2_scientific_completion import publish_completion
    completion_name='completion.json' if complete else 'partial-completion.json'
    if not complete:
        publish_completion(root,stage='discovery-partial',binding=binding,result={},filename=completion_name)
    # Replay current evidence before publishing authority for dependent work.
    load_completed_discovery(root, binding, allow_partial=not complete, check_outcome=False)
    dependency = json_identity(binding)['sha256']
    artifacts = [reference(root, name, role=role, schema=schema, dependencies=dependency)
        for name, role, schema in [
            ('request.json', 'discovery-request', 'genosapiens.track2-discovery-request/v1'),
            ('mechanism-results.json', 'mechanism-results', 'genosapiens.track2-mechanism-results/v1'),
            ('discovery-results.json', 'discovery-result', 'genosapiens.track2-discovery-results/v2'
                if plan['exploratory_discovery']['schema']=='genosapiens.track2-exploratory-plan/v2'
                else 'genosapiens.track2-discovery-results/v1')]]
    record = dict(schema=SCHEMA, state='complete' if complete else 'partial', invocation_id=invocation_id,
        plan_sha256=json_identity(plan)['sha256'], dependency_sha256=dependency, runtime=binding['runtime'],
        artifacts=artifacts, completion_sha256=file_identity(Path(root)/completion_name)['sha256'])
    atomic_record(Path(root)/'current-outcome.json', record)
    return record


def checked_outcome(root, *, plan_path=None, invocation_id=None, runtime=None, binding=None):
    root = Path(root)
    outcome = small_json(root/'current-outcome.json', maximum=CONTROL_BYTES)
    require(set(outcome) == {'schema', 'state', 'invocation_id', 'plan_sha256', 'dependency_sha256',
        'runtime', 'artifacts', 'completion_sha256'} and outcome['schema'] == SCHEMA
        and outcome['state'] in {'complete', 'partial'}, 'Discovery has no current consumable outcome.')
    require(invocation_id is None or outcome['invocation_id'] == invocation_id,
        'The discovery outcome belongs to another invocation.')
    require(plan_path is None or outcome['plan_sha256'] == file_identity(plan_path)['sha256'],
        'The discovery outcome binds a different plan.')
    require(runtime is None or outcome['runtime'] == runtime, 'The discovery runtime changed.')
    require(binding is None or outcome['dependency_sha256'] == json_identity(binding)['sha256'],
        'The discovery outcome dependencies changed.')
    from genosapiens.track2_scientific_reader import read_scientific
    result_schema=read_scientific(root/'discovery-results.json')['schema']
    require(result_schema in {'genosapiens.track2-discovery-results/v1','genosapiens.track2-discovery-results/v2'},
        'The discovery result schema is unsupported.')
    expected = [
        ('request.json', 'discovery-request', 'genosapiens.track2-discovery-request/v1'),
        ('mechanism-results.json', 'mechanism-results', 'genosapiens.track2-mechanism-results/v1'),
        ('discovery-results.json', 'discovery-result', result_schema)]
    require(len(outcome['artifacts']) == len(expected), 'Discovery outcome artifacts are incomplete.')
    for record, (name, role, schema) in zip(outcome['artifacts'], expected):
        require(record['path'] == name, 'Discovery outcome artifact roles changed.')
        check_reference(root, record, role=role, schema=schema, dependencies=outcome['dependency_sha256'])
    completion = root/'completion.json'
    if outcome['state'] == 'complete':
        require(file_identity(completion)['sha256'] == outcome['completion_sha256'], 'Discovery completion changed.')
    else:
        require(not completion.exists(),
            'Partial discovery has inconsistent completion evidence.')
        completion=root/'partial-completion.json'
        require(file_identity(completion)['sha256']==outcome['completion_sha256'], 'Partial discovery evidence changed.')
    from genosapiens.track2_scientific_completion import check_completion
    header=small_json(completion,maximum=CONTROL_BYTES)
    require(json_identity(header['binding'])['sha256']==outcome['dependency_sha256'],
        'Discovery completion dependencies changed.')
    check_completion(root,stage='discovery' if outcome['state']=='complete' else 'discovery-partial',
        binding=header['binding'],filename=completion.name)
    return outcome
