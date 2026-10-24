"""Bounded spawned scoring processes with coordinator-owned publication."""
import multiprocessing
from multiprocessing.connection import wait
import os
from pathlib import Path
import shutil

from genosapiens.track2_operator_inputs import require


def effective_workers(requested, ready, *, capacity=None):
    require(type(requested) is int and requested in (1,2), 'LINCS workers must be 1 or 2.')
    if capacity is None:
        capacity=len(os.sched_getaffinity(0)) if hasattr(os,'sched_getaffinity') else os.cpu_count() or 1
        try:
            quota,period=Path('/sys/fs/cgroup/cpu.max').read_text().split()
            if quota!='max':capacity=min(capacity,max(1,int(quota)//int(period)))
        except (OSError,ValueError):
            pass
    return min(requested,max(1,int(capacity)),ready)


def _worker(channel, work):
    from genosapiens.track2_discovery_outcome import safe_diagnostic
    try:
        while True:
            assignment=channel.recv()
            if assignment is None:return
            try:
                result=work(assignment)
                channel.send(dict(state='completed',result=result))
            except BaseException as error:
                channel.send(dict(state='failed',diagnostic=safe_diagnostic('tissue-drug','scoring-worker',error)))
    except (EOFError,BrokenPipeError):
        return
    finally:
        channel.close()


def _stop_worker(process, channel):
    if process.is_alive():process.terminate()
    process.join(timeout=2)
    if process.is_alive():process.kill();process.join(timeout=2)
    channel.close()


def run_bounded_scoring(assignments, workers, *, work, finish, counters):
    """At most one job per process and one decoded block per scoring stream."""
    require(workers in (1,2), 'LINCS workers must be 1 or 2.')
    context=multiprocessing.get_context('spawn')
    slots=[]
    def start():
        parent,child=context.Pipe()
        process=context.Process(target=_worker,args=(child,work))
        process.start();child.close()
        return process,parent
    try:
        for _ in range(min(workers,len(assignments))):slots.append(start())
        for offset in range(0,len(assignments),workers):
            batch=assignments[offset:offset+workers]
            pending={};outcomes={}
            for index,assignment in enumerate(batch):
                process,channel=slots[index]
                if not process.is_alive():
                    channel.close();process.join();slots[index]=start()
                    process,channel=slots[index]
                try:
                    channel.send(assignment.get('worker_payload',assignment))
                    pending[channel]=index
                except (BrokenPipeError, EOFError, ConnectionResetError):
                    outcomes[index]=dict(state='failed',reason='scoring_worker_exited')
                    _stop_worker(process,channel)
            counters['max_pending_scoring_jobs']=max(counters.get('max_pending_scoring_jobs',0),len(pending))
            while pending:
                for channel in wait(list(pending),timeout=.1):
                    index=pending.pop(channel)
                    try:outcomes[index]=channel.recv()
                    except (EOFError, ConnectionResetError):
                        outcomes[index]=dict(state='failed',reason='scoring_worker_exited')
                        # EOF can precede the process exit observed by is_alive().
                        # Reap it before the next batch uses this slot.
                        _stop_worker(slots[index][0],channel)
                for channel,index in list(pending.items()):
                    if not slots[index][0].is_alive():
                        outcomes[index]=dict(state='failed',reason='scoring_worker_exited')
                        del pending[channel]
            for index,assignment in enumerate(batch):
                finish(assignment,outcomes[index])
    finally:
        for process,channel in slots:
            _stop_worker(process,channel)


def score_assignment(assignment):
    from genosapiens.track2_discovery import _read
    from genosapiens.track2_tissue_execution import execute_tissue_drug_job
    assignment=_read(assignment['assignment_path'])
    resources=dict(assignment['resources'],execution_counters={})
    result=execute_tissue_drug_job(assignment['job'],assignment['query'],resources,assignment['attempt'])
    return dict(state=result['state'],counters=resources.get('execution_counters',{}))


def execute_scoring_jobs(jobs, queries, resources, output, *, requested=2):
    from genosapiens.track2_discovery import _read
    from genosapiens.track2_tissue_execution import tissue_drug_binding, validate_tissue_job, _begin
    from genosapiens.track2_scientific_completion import check_completion
    from genosapiens.track2_discovery_outcome import record_failure
    from genosapiens.track2_scientific_storage import publish_json
    require(all(q['state'] in {'completed','failed'} for q in queries.values()),
        'Scoring requires terminal genetic query attempts.')
    output=Path(output);results={};assignments=[]
    for job in sorted(jobs,key=lambda j:j['job_id']):
        key=job['job_id'];query=queries[job['query_job_id']]['query'];target=output/'jobs'/key
        binding=tissue_drug_binding(job,query,resources)
        try:
            check_completion(target,stage='tissue-drug',binding=binding)
            results[key]=validate_tissue_job(target,revision_root=output)
            continue
        except (OSError,ValueError,KeyError,TypeError):pass
        local={k:resources[k] for k in ('descriptor','reference_root','policy','runtime','revision_root')}
        for name in ('reference_root','revision_root'):local[name]=str(local[name])
        local.update(shared_blocks_by_job={key:resources['shared_blocks_by_job'][key]},execution_counters={})
        assignment=dict(job=job,query=query,resources=local,target=str(target),
            attempt=str(output/'jobs'/(key+'.attempt')),binding=binding)
        assignment_path=output/'scoring-assignments'/(key+'.json')
        publish_json(assignment_path,assignment)
        assignment['worker_payload']=dict(assignment_path=str(assignment_path))
        assignments.append(assignment)
    workers=effective_workers(requested,len(assignments))
    counters=resources.setdefault('execution_counters',{})
    counters.update(requested_workers=requested,effective_workers=workers)
    def finish(assignment,outcome):
        job=assignment['job'];key=job['job_id'];attempt=Path(assignment['attempt']);target=Path(assignment['target'])
        if outcome['state']=='completed':
            validate_tissue_job(attempt,revision_root=output)
            # The coordinator alone publishes final job directories in stable order.
            _begin(target,'tissue-drug',assignment['binding'],force_restart=True)
            shutil.rmtree(target)
            attempt.rename(target)
            results[key]=_read(target/'result.json')
            for name,value in outcome['result']['counters'].items():
                counters[name]=counters.get(name,0)+value
        else:
            record_failure(attempt,stage='tissue-drug',role='scoring-worker',error=RuntimeError('Scoring worker failed.'))
            results[key]=dict(job_id=key,state='failed',reason=outcome.get('reason','scientific_job_failed'),medications=[])
    if assignments:run_bounded_scoring(assignments,workers,work=score_assignment,finish=finish,counters=counters)
    return results
