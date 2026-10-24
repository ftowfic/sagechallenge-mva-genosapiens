"""Explicit, exact-package reuse of predecessor tissue products."""
import json
from pathlib import Path
import shutil

from genosapiens.track2_operator_inputs import require
from genosapiens.track2_tissue_records import identity, validate_tissue_record
from genosapiens.track2_discovery_resources import file_digest

RECORD = 'speedup-compatibility-v1.json'
SCHEMA = 'genosapiens.track2-speedup-compatibility/v1'
MAPPINGS = 'genosapiens.track2-speedup-reuse/v1'


def consumer_contract(root):
    """Bind all runtime inputs except the record which contains this digest."""
    from genosapiens.track2_package_files import runtime_files
    root=Path(root).resolve()
    return identity({p.relative_to(root).as_posix():file_digest(p) for p in runtime_files(root)
        if p != root/'config'/'track2'/RECORD})


def compatibility(producer_code, *, root=None):
    from genosapiens.track2_discovery import scientific_identity
    root=Path(root or Path(__file__).resolve().parents[2])
    record=json.loads((root/'config'/'track2'/RECORD).read_bytes())
    contracts=record.get('consumer',{}).get('runtime_files_sha256')
    require(isinstance(contracts,list) and contracts and all(isinstance(v,str) and len(v)==64
        and all(c in '0123456789abcdef' for c in v) for v in contracts),
        'The consumer runtime inventory is invalid.')
    require(record['schema']==SCHEMA and record['status']=='validated'
        and record['producer']['scientific_code_sha256']==producer_code
        and record['consumer']['scientific_code_sha256']==scientific_identity()
        and consumer_contract(root) in contracts,
        'This package transition has no matching compatibility evidence.')
    require(record['evidence'] and all(r['status']=='passed' and len(r['sha256'])==64 for r in record['evidence']),
        'Package recovery evidence is incomplete.')
    return record


def normalized_plan(plan):
    """Replace only deterministic query and job identifiers, preserving all scope."""
    from genosapiens.track2_scientific_storage import json_chunks
    plan=json.loads(b''.join(json_chunks(plan)))
    names={}
    for query in plan['tissue_queries']:
        names[query['query_job_id']]='query:'+identity({k:v for k,v in query.items() if k!='query_job_id'})
    for job in plan['genetic_selections']:
        names[job['job_id']]='job:'+identity({k:(names[v] if k=='query_job_id' else v)
            for k,v in job.items() if k!='job_id'})
    def replace(value):
        if isinstance(value,dict):return {k:replace(v) for k,v in value.items()}
        if isinstance(value,list):return [replace(v) for v in value]
        return names.get(value,value) if isinstance(value,str) else value
    plan=replace(plan)
    for field in ('tissue_queries','genetic_selections','scenario_job_bindings'):
        plan[field]=sorted(plan[field],key=identity)
    return plan


def validate_prior_plan(plan, expected, resources):
    prior=resources.get('speedup_predecessor_request')
    if prior is None:return False
    compatibility(prior['bindings']['scientific_code_sha256'])
    return (plan==prior['plan']['exploratory_discovery']
        and normalized_plan(plan)==normalized_plan(expected))


def selection_key(selection):
    # Query linkage is checked separately against the mapped genetic query.
    return identity({k:v for k,v in selection.items() if k not in {'job_id','query_job_id'}})


def read_mapping(output, binding=None):
    from genosapiens.track2_discovery import _read
    path=Path(output)/'reuse-mappings.json'
    if not path.exists():return None
    mapping=_read(path)
    require(mapping['schema']==MAPPINGS,'The recovery mapping schema changed.')
    record=compatibility(mapping['producer_code'])
    require(mapping['compatibility_sha256']==identity(record),'The recovery compatibility record changed.')
    if binding is None:binding=_read(Path(output)/'request.json')['bindings']
    require(mapping['consumer_binding']==binding,'The recovery mapping binds another discovery request.')
    return mapping


def checked_product(path, entry, selection, binding, stage, *, revision_root):
    from genosapiens.track2_discovery import _read
    from genosapiens.track2_scientific_completion import check_completion
    from genosapiens.track2_tissue_execution import validate_tissue_job
    path=Path(path)
    require(path.resolve()==path.absolute(),'Recovery products cannot contain symbolic links.')
    require(file_digest(path/'completion.json')==entry['completion_sha256']
        and file_digest(path/'result.json')==entry['result_sha256'],'A producer artifact changed.')
    receipt=_read(path/'completion.json');saved=_read(path/'selection.json')
    require(receipt['binding']['code']==entry['producer_code']
        and receipt['binding']['runtime']==binding['runtime']
        and receipt['binding']['resources']==binding['resources']
        and receipt['binding']['selection_sha256']==identity(saved)
        and selection_key(saved)==selection_key(selection)
        and entry['selection_sha256']==identity(selection),'A producer dependency changed.')
    check_completion(path,stage=stage,binding=receipt['binding'])
    result=validate_tissue_job(path,revision_root=revision_root) if stage=='tissue-drug' else _read(path/'result.json')
    key='job_id' if stage=='tissue-drug' else 'query_job_id'
    require(result[key]==saved[key]==entry['producer_id'],'The producer attribution changed.')
    query=_read(path/'query.json')
    validate_tissue_record(query,'query')
    require(query['scientific_code_sha256']==entry['producer_code']
        and query['source_sha256']==binding['resources']['matrix_sha256']
        and query['policy_sha256']==binding['policy_sha256'], 'The producer query contract changed.')
    if stage=='tissue-query':require(result['query']==query,'The producer query result changed.')
    return result


def mapped_result(output, key, selection, mapping, *, stage):
    if mapping is None or key not in mapping['products']:return None
    entry=mapping['products'][key]
    require(entry['producer_code']==mapping['producer_code'],'A recovery mapping changed producer.')
    return checked_product(Path(output)/'jobs'/key,entry,selection,mapping['consumer_binding'],stage,revision_root=output)


def preserve_recovery_scope(output):
    from genosapiens.track2_run_state import files
    output=Path(output)
    closure=files(output);retained=output.parent/(output.name+'-retained-'+identity(closure))
    if not retained.exists():shutil.copytree(output,retained)
    else:require(files(retained)==closure,'The preserved discovery revision changed.')
    return dict(schema='genosapiens.track2-speedup-recovery-plan/v1',preserved_revision=str(retained),
        execution_started=False,recompute=['discovery'],preparation='requires_separate_compatibility_check',
        annotation='requires_separate_compatibility_check')


def reject_transition(output):
    from genosapiens.track2_scientific_storage import publish_json
    report=preserve_recovery_scope(output)
    publish_json(Path(output)/'recovery-plan.json',dict(report,status='unsupported'))


def prepare_transition(output, plan, binding):
    """Preserve first, report scope, then reuse only authenticated complete products."""
    from genosapiens.track2_discovery import _read
    from genosapiens.track2_run_state import files
    from genosapiens.track2_scientific_storage import publish_json
    output=Path(output)
    if not (output/'request.json').exists():return read_mapping(output,binding)
    prior=_read(output/'request.json')
    if prior['bindings']==binding:return read_mapping(output)
    report=preserve_recovery_scope(output);retained=Path(report['preserved_revision'])
    # The existing legacy-to-tissue transition starts a distinct revision and
    # reuses no legacy scientific products. Preserve that route's behavior.
    if (prior['plan']['exploratory_discovery']['schema']=='genosapiens.track2-exploratory-plan/v1'
            and prior['bindings']['scientific_code_sha256']==binding['scientific_code_sha256']):
        return None,dict(report,status='recompute',reuse=[])
    try:
        record=compatibility(prior['bindings']['scientific_code_sha256'])
        require(all(prior['bindings'][k]==binding[k] for k in ('resources','policy_sha256','runtime')),
            'The source, policy, or runtime changed.')
        require(prior['plan']['mechanism_assessment']==plan['mechanism_assessment']
            and normalized_plan(prior['plan']['exploratory_discovery'])==normalized_plan(plan['exploratory_discovery']),
            'The scientific plan changed.')
    except (OSError,ValueError,KeyError,TypeError) as error:
        publish_json(output/'recovery-plan.json',dict(report,status='unsupported'))
        raise ValueError('Unsupported package transition. The original revision is preserved. '
            'Use a separate discovery revision. Read recovery-plan.json before execution.') from error
    mapping=dict(schema=MAPPINGS,producer_code=prior['bindings']['scientific_code_sha256'],
        producer=record['producer'],compatibility_sha256=identity(record),consumer_binding=binding,products={})
    before=prior['plan']['exploratory_discovery'];after=plan['exploratory_discovery']
    for field,key,stage in [('tissue_queries','query_job_id','tissue-query'),('genetic_selections','job_id','tissue-drug')]:
        candidates={selection_key(r):r for r in before[field]}
        for selection in after[field]:
            old=candidates[selection_key(selection)];name=selection[key];source=retained/'jobs'/old[key]
            try:
                entry=dict(producer_id=old[key],producer_code=mapping['producer_code'],
                    selection_sha256=identity(selection),completion_sha256=file_digest(source/'completion.json'),
                    result_sha256=file_digest(source/'result.json'))
                result=checked_product(source,entry,selection,binding,stage,revision_root=retained)
                if stage=='tissue-drug':
                    q=mapping['products'].get(selection['query_job_id'])
                    require(q is not None and q['producer_id']==old['query_job_id'],
                        'The completed drug job has no compatible genetic query.')
                    query=_read(retained/'jobs'/q['producer_id']/'query.json')
                    require(result['query_sha256']==query['query_sha256'],'The producer query linkage changed.')
            except (OSError,ValueError,KeyError,TypeError):continue
            target=output/'jobs'/name
            if name!=old[key]:
                if not target.exists() or files(target)!=files(source):
                    from genosapiens.track2_tissue_execution import _begin
                    attempt=target.with_name(target.name+'.reuse-attempt')
                    _begin(attempt,stage,{},force_restart=True)
                    shutil.copytree(source,attempt,dirs_exist_ok=True)
                    (attempt/'progress.json').unlink(missing_ok=True)
                    # Restore the producer progress file as part of its exact inventory.
                    if (source/'progress.json').exists():shutil.copy2(source/'progress.json',attempt/'progress.json')
                    _begin(target,stage,{},force_restart=True);shutil.rmtree(target)
                    attempt.rename(target)
            mapping['products'][name]=entry
    report.update(status='supported',reuse=sorted(mapping['products']),recompute=[
        r[k] for field,k in [('tissue_queries','query_job_id'),('genetic_selections','job_id')]
        for r in after[field] if r[k] not in mapping['products']]+['discovery_aggregation'])
    # The caller publishes this mapping after its ordinary restart boundary.
    return mapping,report
