"""Streaming native tissue execution with complete planned scope and recovery."""
from pathlib import Path
import shutil
import time

from genosapiens.track2_discovery_resources import file_digest, load_discovery_resources
from genosapiens.track2_operator_inputs import require
from genosapiens.track2_run_state import RunState, files
from genosapiens.track2_scientific_completion import check_completion
from genosapiens.track2_scientific_storage import publish_json, Records
from genosapiens.track2_tissue_records import identity, validate_tissue_record


def tissue_outcomes(nested,mechanisms,jobs,queries,union,bindings):
    outcomes=[]
    for route in nested['routes']:
        oid=route['objective_id']
        if route['route']=='qualified_target_drug_linkage':
            local=next(r for r in mechanisms if r['objective_id']==oid)
            state='completed_with_candidates' if local['nominations'] else 'completed_no_candidates'
            executed=True; refs=['mechanism-results.json']; reasons=[]
        else:
            inventory=next(r for r in nested['coverage'] if r['objective_id']==oid)['tissues']
            ids=[j for t in inventory for j in t['drug_job_ids']]
            failed=any(jobs[j]['state']=='failed' for j in ids) or any(queries[t['query_job_id']]['state']=='failed' for t in inventory)
            executed=any(jobs[j]['state']=='completed' for j in ids)
            nominated=any(r['objective_id']==oid and r['lincs_nominated'] for r in union['concordance'])
            state='failed' if failed else ('completed_with_candidates' if nominated else 'completed_no_candidates') if executed else 'unavailable'
            reasons=['incomplete_scientific_jobs'] if failed else [] if executed else ['no_evaluable_tissue_query_or_drug_coverage']
            refs=['jobs/'+j+'/result.json' for j in ids]
        outcomes.append(dict(route,state=state,executed=executed,query_ids=[route['request_id']],
            source_bindings=bindings,reason_codes_result=reasons,result_references=refs))
    return outcomes


def execute_tissue_discovery(plan,resources,output,runtime,invocation_id):
    from genosapiens.track2_discovery import _read, discovery_bindings, load_completed_discovery
    from genosapiens.track2_discovery_outcome import invalidate_outcome, publish_outcome, record_failure
    from genosapiens.track2_drug_mechanisms import lookup_gene_mechanisms, nominate_for_objective
    from genosapiens.track2_tissue_execution import execute_tissue_query, execute_tissue_drug_job, admit_tissue_drug_job
    from genosapiens.track2_shared_blocks import plan_shared_blocks, publish_shared_block
    from genosapiens.track2_scoring_workers import execute_scoring_jobs
    from genosapiens.track2_tissue_union import reconcile_tissue_results
    from genosapiens.track2_tissue_admission import TissueAdmission
    from genosapiens.track2_speedup_recovery import prepare_transition, mapped_result, read_mapping
    output=Path(output); state=RunState(output,direct_stage='discovery'); connection=None
    started=time.monotonic()
    try:
        binding=discovery_bindings(plan,resources,runtime)
        complete,reason=state.completed('discovery',binding)
        if complete:
            result=load_completed_discovery(output,binding,check_outcome=False)
            from genosapiens.track2_run_state import atomic_record
            atomic_record(output/'current-execution.json',dict(requested_workers=resources.get('lincs_workers',2),
                effective_workers=0,reused_discovery=True,elapsed_seconds=time.monotonic()-started))
            publish_outcome(output,plan,binding,complete=True,invocation_id=invocation_id)
            return result
        transition=prepare_transition(output,plan,binding)
        mapping,report=transition if isinstance(transition,tuple) else (transition,None)
        invalidate_outcome(output)
        state.restart('discovery',reason)
        if mapping is not None:publish_json(output/'reuse-mappings.json',mapping)
        if report is not None:publish_json(output/'recovery-plan.json',report)
        publish_json(output/'request.json',dict(schema='genosapiens.track2-discovery-request/v1',bindings=binding,plan=plan))
        # Atomic replacement invalidates disk-backed views from the prior attempt.
        if mapping is not None:mapping=read_mapping(output)
        nested=plan['exploratory_discovery']; jobs={}; queries={}
        resources=dict(resources,runtime=runtime,extraction_cache={},reader_session_enabled=True,
            reader_session=None,reader_directory=output/'reader-session',execution_counters={})
        connection=load_discovery_resources(resources['descriptor'],resources['reference_root'])
        resources['tissue_admission']=TissueAdmission(resources['descriptor'],resources['policy'],connection,
            counters=resources.setdefault('execution_counters',{}))
        publish_json(output/'mechanism-results.json',Records(nominate_for_objective(o,
            lookup_gene_mechanisms(connection,o['gene_id']),connection,resources['policy']) for o in nested['objectives']))
        mechanisms=_read(output/'mechanism-results.json')
        # Complete all prospective genetic queries before reading any drug values.
        for selection in nested['tissue_queries']:
            key=selection['query_job_id']
            try:
                queries[key]=mapped_result(output,key,selection,mapping,stage='tissue-query') or execute_tissue_query(selection,resources,output/'jobs'/key)
            except Exception:
                queries[key]=dict(query_job_id=key,state='failed',reason='genetic_query_job_failed')
        admitted=[]
        for job in nested['genetic_selections']:
            q=queries[job['query_job_id']]
            if q['state']=='completed':
                try:
                    admit_tissue_drug_job(job,q['query'],resources,resources['tissue_admission'])
                    saved=mapped_result(output,job['job_id'],job,mapping,stage='tissue-drug')
                    if saved is not None:
                        require(saved['query_sha256']==q['query']['query_sha256'],'The reused drug query changed.')
                        jobs[job['job_id']]=saved
                    else:admitted.append(job)
                except Exception as error:
                    record_failure(output/'jobs'/job['job_id'],stage='tissue-drug',role='admission',error=error)
                    jobs[job['job_id']]=dict(job_id=job['job_id'],state='failed',reason='drug_admission_failed',medications=[])
        resources['revision_root']=output
        resources['shared_blocks_by_job']={job['job_id']:[] for job in admitted}
        for block in plan_shared_blocks(admitted,queries,resources):
            try:
                ref=publish_shared_block(output,block,resources)
                for key in block['consumers']:resources['shared_blocks_by_job'][key].append(ref)
            except Exception:
                for key in block['consumers']:
                    jobs[key]=dict(job_id=key,state='failed',reason='shared_extraction_failed',medications=[])
        for job in nested['genetic_selections']:
            if queries[job['query_job_id']]['state']=='failed':
                jobs[job['job_id']]=dict(job_id=job['job_id'],state='failed',reason='genetic_query_job_failed',medications=[])
        jobs.update(execute_scoring_jobs([j for j in admitted if j['job_id'] not in jobs],queries,resources,output,
            requested=resources.get('lincs_workers',2)))
        union=reconcile_tissue_results(nested,mechanisms,jobs,queries,connection,resources['policy'])
        outcomes=tissue_outcomes(nested,mechanisms,jobs,queries,union,binding['resources'])
        complete=all(r['state']!='failed' for r in [*jobs.values(),*queries.values()])
        if resources.get('reader_session') is not None:
            resources['reader_session'].close()
        # Recheck the admitted descriptor and retained matrix identity before publication.
        checked=load_discovery_resources(resources['descriptor'],resources['reference_root'])
        checked.close()
        from genosapiens.track2_run_state import atomic_record
        atomic_record(output/'current-execution.json',dict(resources['execution_counters'],elapsed_seconds=time.monotonic()-started))
        def references(records,key):
            return [dict(**{key:k},path='jobs/'+k+'/result.json',state=r['state'],
                completion_sha256=file_digest(output/'jobs'/k/'completion.json') if r['state']!='failed' else None)
                for k,r in sorted(records.items())]
        result=dict(schema='genosapiens.track2-discovery-results/v2',bindings=binding,route_outcomes=outcomes,
            jobs=references(jobs,'job_id'),queries=references(queries,'query_job_id'),**union,
            coverage=nested['coverage'],completion_state='complete' if complete else 'incomplete')
        publish_json(output/'discovery-results.json',result)
        if complete: state.publish('discovery',binding=binding,result=dict(result_sha256=file_digest(output/'discovery-results.json')))
        publish_outcome(output,plan,binding,complete=complete,invocation_id=invocation_id)
        return _read(output/'discovery-results.json')
    except BaseException as error:
        record_failure(output,stage='discovery',role='tissue-discovery',error=error,invocation_id=invocation_id)
        raise
    finally:
        if resources.get('reader_session') is not None:
            resources['reader_session'].abort()
        if connection is not None: connection.close()


def retained_tissue_jobs(path,result):
    from genosapiens.track2_discovery import _read, validate_job
    from genosapiens.track2_speedup_recovery import read_mapping,mapped_result
    path=Path(path); jobs={};queries={};shared_checked=set();mapping=read_mapping(path)
    nested=_read(path/'request.json')['plan']['exploratory_discovery'] if mapping is not None else None
    for entries,key,records in [(result['jobs'],'job_id',jobs),(result['queries'],'query_job_id',queries)]:
        for entry in entries:
            name=entry[key]; target=path/'jobs'/name
            require(entry['path']=='jobs/'+name+'/result.json','The tissue result reference changed.')
            if entry['completion_sha256'] is None:
                require(entry['state']=='failed','A completed tissue job lacks its completion.')
                records[name]=dict(**{key:name},state='failed',medications=[])
                continue
            require(file_digest(target/'completion.json')==entry['completion_sha256'],'A tissue job completion changed.')
            if mapping is not None and name in mapping['products']:
                selection=next(r for r in nested['genetic_selections' if key=='job_id' else 'tissue_queries'] if r[key]==name)
                records[name]=mapped_result(path,name,selection,mapping,stage='tissue-drug' if key=='job_id' else 'tissue-query')
                require(records[name]['state']==entry['state'],'A reused tissue state changed.')
                continue
            if key=='job_id':
                from genosapiens.track2_tissue_execution import validate_tissue_job
                records[name]=validate_tissue_job(target,revision_root=path,shared_checked=shared_checked)
            else:
                receipt=_read(target/'completion.json')
                check_completion(target,stage='tissue-query',binding=receipt['binding'])
                records[name]=_read(target/'result.json')
                validate_tissue_record(records[name]['query'],'query')
                require(receipt['binding']['selection_sha256']==identity(_read(target/'selection.json')),
                    'The retained tissue contributor selection changed.')
            require(records[name][key]==name and records[name]['state']==entry['state'],'A tissue job identity changed.')
    return jobs,queries


def validate_tissue_discovery(path,binding,result,*,complete,allow_partial):
    from genosapiens.track2_discovery import _read
    path=Path(path);request=_read(path/'request.json');nested=request['plan']['exploratory_discovery']
    require(result['bindings']==binding==request['bindings'] and result['coverage']==nested['coverage']
        and result['completion_state']==('complete' if complete else 'incomplete'), 'Tissue discovery bindings changed.')
    require({r['job_id'] for r in result['jobs']}=={r['job_id'] for r in nested['genetic_selections']}
        and {r['query_job_id'] for r in result['queries']}=={r['query_job_id'] for r in nested['tissue_queries']},
        'Tissue discovery omitted planned jobs or queries.')
    jobs,queries=retained_tissue_jobs(path,result)
    require(allow_partial or all(r['state']!='failed' for r in [*jobs.values(),*queries.values()]),
        'Required tissue discovery is incomplete.')
    for job in nested['genetic_selections']:
        local=jobs[job['job_id']];query=queries[job['query_job_id']]
        if local['state']=='failed': continue
        require(query['state']=='completed' and local['query_sha256']==query['query']['query_sha256'],
            'A drug-cell result binds a different pooled query.')
    expected=tissue_outcomes(nested,_read(path/'mechanism-results.json'),jobs,queries,result,binding['resources'])
    require(result['route_outcomes']==expected,'Tissue route outcomes differ from retained evidence.')
    return result
