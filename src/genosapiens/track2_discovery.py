"""Offline native scientific jobs with retained queries and atomic completion."""
import csv
from contextlib import contextmanager
import fcntl
from functools import lru_cache
import json
import os
from pathlib import Path
import shutil
import subprocess

from genosapiens.track2_operator_inputs import require
from genosapiens.track2_published_genotype import canonical,sha256
from genosapiens.track2_run_state import RunState,atomic_record,files
from genosapiens.track2_discovery_resources import load_discovery_resources,resource_bindings,public_path,file_digest
from genosapiens.track2_scientific_storage import Records,json_identity, publish_json, publish_chunks, json_chunks
from genosapiens.track2_scientific_completion import check_completion, publish_completion


@lru_cache(maxsize=1)
def scientific_identity():
    from genosapiens.track2_construction_dependencies import python_dependencies,STORAGE_DEPENDENCIES
    root=Path(__file__).resolve().parents[2]
    return sha256(canonical(python_dependencies(root,[('genosapiens.track2_lincs_scoring','refine_exploratory_query'),
        ('genosapiens.track2_lincs_scoring','score_cell'),('genosapiens.track2_lincs_selection','select_genetic_signatures'),
        ('genosapiens.track2_lincs_selection','select_drug_conditions'),
        ('genosapiens.track2_tissue_queries','refine_tissue_query'),
        ('genosapiens.track2_lincs_selection','select_tissue_genetic_signatures'),
        ('genosapiens.track2_tissue_scoring','score_tissue_condition'),
        ('genosapiens.track2_tissue_execution','execute_tissue_query'),
        ('genosapiens.track2_tissue_execution','execute_tissue_drug_job'),
        ('genosapiens.track2_discovery','execute_lincs_job')],ignored=STORAGE_DEPENDENCIES) | {
            'reader':file_digest(root/'scripts/read_track2_lincs.R')}))


def discovery_bindings(plan,resources,runtime):
    from genosapiens.track2_construction_dependencies import python_dependencies
    return dict(plan_sha256=json_identity({k:plan[k] for k in ['mechanism_assessment','exploratory_discovery']})['sha256'],resources=resource_bindings(resources['descriptor']),
        policy_sha256=sha256(canonical(resources['policy'])),runtime=runtime,scientific_code_sha256=scientific_identity(),
        aggregation_code_sha256=sha256(canonical(python_dependencies(Path(__file__).resolve().parents[2],
            [('genosapiens.track2_drug_mechanisms','nominate_for_objective'),
             ('genosapiens.track2_discovery_union','reconcile_results'),
             ('genosapiens.track2_discovery','execute_discovery')]))))


def _read(path):
    from genosapiens.track2_scientific_reader import read_scientific
    return read_scientific(path)


def _extract(row_ids,signature_ids,resources,output):
    descriptor=resources['descriptor'];matrix=public_path(resources['reference_root'],descriptor['matrix']['path'])
    request=dict(schema='genosapiens.track2-lincs-discovery-source/v1',matrix=str(matrix),row_ids=row_ids,
        signature_ids=signature_ids,expected_dimensions=resources['dimensions'],
        block_size=resources['policy']['lincs']['column_block_size'],output=str(output))
    from genosapiens.json_schema import validate_json_schema
    root=Path(__file__).resolve().parents[2]
    validate_json_schema(request,_read(root/'schemas/track2/lincs-discovery-source.schema.json'))
    cache=resources.get('extraction_cache',{})
    cache_key=sha256(canonical(dict(rows=row_ids,columns=signature_ids,matrix=descriptor['matrix']['content_sha256'])))
    request_path=output.with_suffix('.request.json');publish_json(request_path,request)
    if cache_key in cache:
        output=cache[cache_key]
    elif resources.get('reader_session_enabled'):
        from genosapiens.track2_reader_session import ReaderSession
        if resources.get('reader_session') is None:
            checked=load_discovery_resources(descriptor,resources['reference_root'])
            checked.close()
            resources['reader_session']=ReaderSession(dict(matrix=str(matrix),
                expected_dimensions=resources['dimensions'],matrix_identity=descriptor['matrix']),
                resources['reader_directory'],timeout=resources.get('reader_timeout',900),
                counters=resources.setdefault('execution_counters',{}))
        blocks=resources['reader_session'].extract(request)
        require(_read(output/'completion.json')['blocks']==blocks, 'The reader artifact inventory changed.')
    else:
        subprocess.run(['Rscript',str(root/'scripts/read_track2_lincs.R'),'--indexed',str(request_path)],
            check=True,stdout=subprocess.DEVNULL,stderr=subprocess.DEVNULL)
    completion=_read(output/'completion.json')
    require(completion['status']=='complete','Matrix extraction is incomplete.')
    seen=[]
    for block in completion['blocks']:
        with (output/block['path']).open() as stream:rows=list(csv.DictReader(stream))
        require([r['gene_id'] for r in rows]==row_ids and set(rows[0])=={'gene_id',*block['signature_ids']},
            'Extracted matrix row or column identity differs.')
        require(block['columns']<=request['block_size'],'Matrix extraction exceeded its column bound.')
        seen.extend(block['signature_ids'])
        yield {name:[float(r[name]) if r[name] else None for r in rows] for name in block['signature_ids']}
    require(seen==signature_ids,'Extracted signature inventory differs.')
    cache[cache_key]=output


def execute_lincs_job(job,resources,output):
    from genosapiens.track2_lincs_selection import select_drug_conditions
    from genosapiens.track2_lincs_scoring import refine_exploratory_query,score_cell,summarize_conditions
    output=Path(output)
    connection=None
    try:
        binding=dict(selection_sha256=json_identity(job)['sha256'],resources=resource_bindings(resources['descriptor']),
            code=scientific_identity(),runtime=resources['runtime'])
        completion=output/'completion.json'
        if completion.exists():
            try:
                check_completion(output,stage='lincs-job',binding=binding)
            except (ValueError,OSError,KeyError,TypeError):pass
            else:return _read(output/'result.json')
        if output.exists():
            files(output);shutil.rmtree(output)
        output.mkdir(parents=True,mode=0o700)
        publish_json(output/'selection.json',job)
        atomic_record(output/'progress.json',dict(status='running',binding=binding))
        connection=load_discovery_resources(resources['descriptor'],resources['reference_root'])
        features=[dict(r) for r in connection.execute('SELECT pr_gene_id,symbol,measured_landmark FROM lincs_features WHERE measured_landmark=1 ORDER BY pr_gene_id')]
        ids=[r['pr_gene_id'] for r in features]
        genetic={}
        for block in _extract(ids,[job['selected'][d]['sig_id'] for d in sorted(job['selected'])],resources,output/'genetic'):
            genetic.update(block)
        query=refine_exploratory_query(job,genetic,features,source_sha256=resources['descriptor']['matrix']['content_sha256'],policy=resources['policy'])
        publish_json(output/'query.json',query)  # Freeze before any drug values.
        conditions=select_drug_conditions(connection,job['cell_id'],None,resources['policy'])
        require([r['sig_id'] for r in conditions['conditions']]==job['drug_signature_ids'],'Prospective drug metadata scope changed.')
        metadata={r['sig_id']:r for r in conditions['conditions']}
        summaries=[]
        def condition_records():
            if query['state']=='ready':
                for block in _extract(ids,job['drug_signature_ids'],resources,output/'drugs'):
                    cell=dict(cell_id=job['cell_id'],genes=ids,exclusions=[],conditions=[metadata[s]|dict(values=v) for s,v in block.items()])
                    scored=score_cell(query,cell)
                    for row in scored['conditions']:
                        summaries.append({k:v for k,v in row.items() if k not in {'values','gene_ids','medication'}})
                        yield from json_chunks(row)
        publish_chunks(output/'conditions.jsonl',condition_records())
        result=dict(schema='genosapiens.track2-lincs-result/v2',job_id=job['job_id'],gene_id=job['gene_id'],
            cell_id=job['cell_id'],method=job['method'],query_sha256=query['query_sha256'],gene_count=query['gene_count'],
            state='completed' if query['state']=='ready' else 'unavailable',reason=query['reason'],
            medications=summarize_conditions(summaries),exclusions=conditions['exclusions'],conditions_path='conditions.jsonl',
            query_path='query.json',source_sha256=resources['descriptor']['matrix']['content_sha256'])
        publish_json(output/'result.json',result)
        atomic_record(output/'progress.json',dict(status='complete',binding=binding))
        publish_completion(output,stage='lincs-job',binding=binding,result={})
        return _read(output/'result.json')
    except BaseException as error:
        from genosapiens.track2_discovery_outcome import record_failure
        record_failure(output,stage='lincs-job',role='job-artifacts',error=error)
        raise
    finally:
        if connection is not None:connection.close()


def derive_route_outcomes(nested,mechanisms,jobs,bindings):
    outcomes=[]
    for route in nested['routes']:
        if route['route']=='qualified_target_drug_linkage':
            local=next(r for r in mechanisms if r['objective_id']==route['objective_id'])
            outcome='completed_with_candidates' if local['nominations'] else 'completed_no_candidates'
            refs=['mechanism-results.json'];executed=True;reasons=[]
        else:
            selected=[jobs[b['job_id']] for b in nested['scenario_job_bindings'] if b['objective_id']==route['objective_id']]
            executed=any(r['state']=='completed' for r in selected)
            outcome='failed' if any(r['state']=='failed' for r in selected) else ('completed_with_candidates' if any(m['state']=='pass' for r in selected for m in r['medications'])
                else 'completed_no_candidates') if executed else 'unavailable'
            refs=['jobs/'+r['job_id']+'/result.json' for r in selected]
            reasons=route['reason_codes'] if not selected else sorted({r['reason'] for r in selected if r['reason']})
        outcomes.append(dict(**route,state=outcome,executed=executed,query_ids=[route['request_id']],
            source_bindings=bindings,reason_codes_result=reasons,result_references=refs))
    return outcomes


@contextmanager
def discovery_lock(output):
    output=Path(output).absolute()
    require(output.resolve()==output,'Discovery output contains a symbolic link.')
    output.mkdir(parents=True,mode=0o700,exist_ok=True)
    descriptor=os.open(output/'.discovery.lock',os.O_CREAT|os.O_RDWR|os.O_NOFOLLOW,0o600)
    try:
        fcntl.flock(descriptor,fcntl.LOCK_EX|fcntl.LOCK_NB)
        yield output
    finally:os.close(descriptor)


def discover_file(plan_path,reference_root,output,*,invocation_id=None,lincs_workers=2):
    """Include plan and resource admission in the locked stage failure boundary."""
    from genosapiens.track2_exploratory_planning import load_plan_resources,validate_exploratory_plan
    from genosapiens.track2_discovery_outcome import invalidate_outcome,record_failure
    with discovery_lock(output) as output:
        role='discovery-request'
        context=None
        try:
            invalidate_outcome(output)
            plan=_read(plan_path)
            role='resource-admission'
            context=load_plan_resources(plan_path,plan['exploratory_discovery'],reference_root)
            if (output/'request.json').exists():
                context['speedup_predecessor_request']=_read(output/'request.json')
            if plan.get('input_origin')=='published_patient_genotype':
                from genosapiens.track2_public_discovery import check_public_plan
                check_public_plan(plan,context)
            else:
                validate_exploratory_plan(plan['exploratory_discovery'],plan['mechanism_assessment'],context)
            context['connection'].close();context['connection']=None
            context.update(reference_root=reference_root,lincs_workers=lincs_workers,
                dimensions=_read(Path(__file__).resolve().parents[2]/'config/track2/discovery-sources.json')['lincs']['dimensions'])
        except BaseException as error:
            if context is not None and context.get('speedup_predecessor_request') is not None:
                prior=context['speedup_predecessor_request']
                if prior['bindings']['scientific_code_sha256']!=scientific_identity():
                    from genosapiens.track2_speedup_recovery import reject_transition
                    reject_transition(output)
            record_failure(output,stage='discovery',role=role,error=error,invocation_id=invocation_id)
            raise
        finally:
            if context is not None and context.get('connection') is not None:context['connection'].close()
        return _execute_discovery(plan,context,output,plan['exploratory_discovery']['runtime_images']['discovery'],invocation_id)


def execute_discovery(plan,resources,output,runtime,*,invocation_id=None):
    with discovery_lock(output) as output:
        return _execute_discovery(plan,resources,output,runtime,invocation_id)


def _execute_discovery(plan,resources,output,runtime,invocation_id):
    if plan['exploratory_discovery']['schema']=='genosapiens.track2-exploratory-plan/v2':
        from genosapiens.track2_tissue_discovery import execute_tissue_discovery
        return execute_tissue_discovery(plan,resources,output,runtime,invocation_id)
    from genosapiens.track2_drug_mechanisms import lookup_gene_mechanisms,nominate_for_objective
    from genosapiens.track2_discovery_union import reconcile_results
    from genosapiens.track2_discovery_outcome import invalidate_outcome,publish_outcome,record_failure
    state=RunState(output,direct_stage='discovery')
    connection=None;role='discovery-request'
    try:
        invalidate_outcome(output)
        binding=discovery_bindings(plan,resources,runtime)
        completed,reason=state.completed('discovery',binding)
        if completed:
            result=load_completed_discovery(output,binding,check_outcome=False)
            publish_outcome(output,plan,binding,complete=True,invocation_id=invocation_id)
            return result
        state.restart('discovery',reason)
        publish_json(output/'request.json',dict(schema='genosapiens.track2-discovery-request/v1',bindings=binding,plan=plan))
        nested=plan['exploratory_discovery'];jobs={}
        resources={**resources,'extraction_cache':{}}
        role='resource-admission'
        connection=load_discovery_resources(resources['descriptor'],resources['reference_root'])
        role='mechanism-results'
        def mechanism_records():
            for objective in nested['objectives']:
                rows=lookup_gene_mechanisms(connection,objective['gene_id'])
                yield nominate_for_objective(objective,rows,connection,resources['policy'])
        publish_json(output/'mechanism-results.json',Records(mechanism_records()))
        results=_read(output/'mechanism-results.json')
        role='scientific-jobs'
        for job in nested['genetic_selections']:
            try:
                jobs[job['job_id']]=execute_lincs_job(job,{**resources,'runtime':runtime},output/'jobs'/job['job_id'])
            except Exception:
                jobs[job['job_id']]=dict(job_id=job['job_id'],gene_id=job['gene_id'],cell_id=job['cell_id'],
                    method=job['method'],state='failed',reason='scientific_job_failed',medications=[])
        role='discovery-result'
        outcomes=derive_route_outcomes(nested,results,jobs,binding['resources'])
        union=reconcile_results(nested,results,jobs,connection,resources['policy'],outcomes)
        complete=not any(r['state']=='failed' for r in jobs.values())
        result=dict(schema='genosapiens.track2-discovery-results/v1',bindings=binding,route_outcomes=outcomes,
            jobs=[dict(job_id=k,path='jobs/'+k+'/result.json',
                completion_sha256=file_digest(output/'jobs'/k/'completion.json') if jobs[k]['state']!='failed' and (output/'jobs'/k/'completion.json').exists() else None,
                state=jobs[k]['state']) for k in sorted(jobs)],
            **union,coverage=nested['coverage'],completion_state='complete' if complete else 'incomplete')
        publish_json(output/'discovery-results.json',result)
        result=_read(output/'discovery-results.json')
        role='discovery-completion'
        if complete:state.publish('discovery',binding=binding,result=dict(result_sha256=file_digest(output/'discovery-results.json')))
        else:atomic_record(output/'current-diagnostics.json',dict(stage='discovery',state='incomplete'))
        role='discovery-outcome'
        publish_outcome(output,plan,binding,complete=complete,invocation_id=invocation_id)
        return result
    except BaseException as error:
        # Completion alone cannot authorize a failed invocation.
        record_failure(output,stage='discovery',role=role,error=error,invocation_id=invocation_id)
        raise
    finally:
        if connection is not None:connection.close()


def validate_job(path, *, revision_root=None, shared_checked=None):
    """Replay aggregates from retained condition scores without reading the matrix."""
    from genosapiens.track2_lincs_scoring import summarize_conditions,validate_query
    path=Path(path);receipt=_read(path/'completion.json')
    if receipt['stage']=='tissue-drug':
        from genosapiens.track2_tissue_execution import validate_tissue_job
        return validate_tissue_job(path,revision_root=revision_root,shared_checked=shared_checked)
    check_completion(path,stage='lincs-job',binding=receipt['binding'])
    result=_read(path/'result.json');query=_read(path/'query.json')
    validate_query(query,dict(cell_id=query['cell_id'],genes=[g['gene_id'] for g in query['genes']]))
    from genosapiens.track2_figure_stream import condition_rows
    # Replay scores without retaining every expression vector and source record.
    conditions=[{k:v for k,v in row.items() if k not in {'values','gene_ids','medication'}}
        for row in condition_rows(path/'conditions.jsonl')]
    require(result['medications']==summarize_conditions(conditions), 'Medication summaries differ from retained conditions.')
    require(result['query_sha256']==query['query_sha256'] and result['gene_count']==query['gene_count']
        and all(r['query_sha256']==query['query_sha256'] for r in conditions),'A scientific query binding differs.')
    selected=_read(path/'selection.json')
    require(receipt['binding']['selection_sha256']==json_identity(selected)['sha256'],
        'The scientific job selection changed.')
    require(query['method']==selected['method'] and query['cell_id']==selected['cell_id']
        and sorted(r['sig_id'] for r in conditions)==(sorted(selected['drug_signature_ids']) if query['state']=='ready' else []),
        'A scientific job omitted or added drug conditions.')
    return result


def load_completed_discovery(path,expected_bindings,*,allow_partial=False,check_outcome=True):
    if check_outcome:
        from genosapiens.track2_discovery_outcome import checked_outcome
        checked_outcome(path,binding=expected_bindings)
    state=RunState(path,direct_stage='discovery')
    complete,reason=state.completed('discovery',expected_bindings)
    require(complete is not None or allow_partial and not (Path(path)/'completion.json').exists(),
        'Required discovery is incomplete or changed: '+str(reason))
    result=_read(Path(path)/'discovery-results.json')
    if result['schema']=='genosapiens.track2-discovery-results/v2':
        from genosapiens.track2_tissue_discovery import validate_tissue_discovery
        return validate_tissue_discovery(path,expected_bindings,result,complete=complete,allow_partial=allow_partial)
    require(result['bindings']==expected_bindings and result['completion_state']==('complete' if complete else 'incomplete'), 'Discovery result bindings differ.')
    request=_read(Path(path)/'request.json');nested=request['plan']['exploratory_discovery']
    require(request['bindings']==expected_bindings and result['coverage']==nested['coverage'],'Discovery request bindings differ.')
    require({(r['objective_id'],r['route']) for r in result['route_outcomes']}==
        {(r['objective_id'],r['route']) for r in nested['routes']},'Discovery route coverage differs.')
    require({r['job_id'] for r in result['jobs']}=={r['job_id'] for r in nested['genetic_selections']},'Discovery job coverage differs.')
    jobs={}
    for job in result['jobs']:
        if job['completion_sha256'] is None:
            require(allow_partial and job['state']=='failed','A required job completion is missing.')
            jobs[job['job_id']]=dict(job_id=job['job_id'],state='failed',reason='scientific_job_failed',medications=[])
            continue
        target=Path(path)/'jobs'/job['job_id']
        require(job['completion_sha256']==file_digest(target/'completion.json'),'A job completion binding changed.')
        retained=validate_job(target)
        require(retained['job_id']==job['job_id'] and retained['state']==job['state'],'The scientific job identity changed.')
        jobs[job['job_id']]=retained
    mechanisms=_read(Path(path)/'mechanism-results.json')
    require(result['route_outcomes']==derive_route_outcomes(nested,mechanisms,jobs,expected_bindings['resources']),
        'Discovery route completion differs from retained jobs.')
    return result
