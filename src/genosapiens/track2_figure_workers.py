"""Bounded R process groups with private attempts and per-job publication."""
import itertools
import os
from pathlib import Path
import signal
import subprocess
import time
import uuid

from genosapiens.track2_operator_inputs import require
from genosapiens.track2_scientific_storage import file_identity, json_identity, publish_json
from genosapiens.track2_scientific_reader import read_scientific


def effective_workers(requested, *, capacity=None):
    require(type(requested) is int and requested in (1,2),'Figure workers must be 1 or 2.')
    from genosapiens.track2_scoring_workers import effective_workers as capacity_workers
    return capacity_workers(requested,2,capacity=capacity)


def _inventory(folder):
    from genosapiens.track2_challenge_revision import inventory
    values=inventory(folder,omit={'render-receipt.json'})
    expected={'figure-data.jsonl'}
    import json
    with (Path(folder)/'figure-data.jsonl').open() as stream:
        for line in stream:
            identifier=json.loads(line)['id']
            require(isinstance(identifier,str) and identifier not in {'.','..'}
                and '/' not in identifier and '\\' not in identifier,'A figure identifier is not local.')
            expected.update(identifier+suffix for suffix in ('.png','.svg'))
    require(set(values)==expected,'A figure job asset inventory differs.')
    return values


def _stop(processes):
    for process in processes:
        if process.poll() is None:
            try:os.killpg(process.pid,signal.SIGTERM)
            except ProcessLookupError:pass
    deadline=time.monotonic()+2
    while any(p.poll() is None for p in processes) and time.monotonic()<deadline:time.sleep(.02)
    for process in processes:
        # Kill the group even if its leader exited before a descendant.
        try:os.killpg(process.pid,signal.SIGKILL)
        except ProcessLookupError:pass
        process.wait()


def render_jobs(jobs, attempts, *, workers=2, capacity=None, command=None):
    """Prepare at most two jobs. Publish each terminal batch in input order.

    Each job supplies target, renderer, runtime, prepare(folder), and check(folder).
    Preparation regenerates or admits page semantics before receipt-based reuse.
    No persisted receipt alone authorizes scientific or page reuse.
    """
    from genosapiens.track2_challenge_overviews import pinned_runtime
    requested_count=workers
    workers=effective_workers(workers,capacity=capacity)
    attempts=Path(attempts).absolute()
    require(attempts.resolve()==attempts,'The figure attempt path contains a symbolic link.')
    attempts.mkdir(mode=0o700,parents=True,exist_ok=True)
    from concurrent.futures import ThreadPoolExecutor,wait
    import threading
    source=iter(jobs);first=list(itertools.islice(source,workers))
    effective=min(workers,len(first))
    outcomes=[];active=[];logs=[];handlers={};cancelled=threading.Event();lock=threading.Lock()
    pool=ThreadPoolExecutor(max_workers=max(1,effective),thread_name_prefix='track2-figure')
    invocation=attempts/('invocation-'+uuid.uuid4().hex+'.json')
    operational=dict(state='incomplete',requested_workers=requested_count,effective_workers=effective)
    publish_json(invocation,operational)
    def execute(argv,log):
        with lock:
            if cancelled.is_set():raise RuntimeError('Figure assignment cancelled.')
            process=subprocess.Popen(argv,stdout=log,stderr=subprocess.STDOUT,start_new_session=True)
            active.append(process)
        process.wait()
        return process
    def interrupted(signum,frame):raise KeyboardInterrupt('Figure rendering interrupted.')
    try:
        import threading
        if threading.current_thread() is threading.main_thread():
            for sig in (signal.SIGTERM,signal.SIGINT):
                handlers[sig]=signal.signal(sig,interrupted)
        while True:
            batch=first if first is not None else list(itertools.islice(source,workers))
            first=None
            if not batch:break
            pending=[]
            for job in batch:
                target=Path(job['target']).absolute()
                require(target.resolve()==target,'The figure target contains a symbolic link.')
                require(job['runtime']==pinned_runtime(),'The figure job runtime differs from the pinned runtime.')
                attempt=attempts/uuid.uuid4().hex;attempt.mkdir(mode=0o700)
                folder=attempt/'artifacts';folder.mkdir(mode=0o700)
                publish_json(attempt/'status.json',dict(state='incomplete',target=str(target)))
                source_binding=job['prepare'](folder)
                binding=dict(source=source_binding,data=file_identity(folder/'figure-data.jsonl'),
                    renderer=file_identity(job['renderer']),runtime=job['runtime'])
                reused=False
                if (target/'render-receipt.json').is_file():
                    try:
                        receipt=read_scientific(target/'render-receipt.json')
                        require(receipt==dict(schema='genosapiens.track2-figure-job/v1',binding=binding,outputs=_inventory(target)),
                            'The retained figure job differs.')
                        job['check'](target)
                        reused=True
                    except (ValueError,OSError,KeyError,TypeError):pass
                if reused:
                    publish_json(attempt/'status.json',dict(state='reused',target=str(target)))
                    outcomes.append(dict(target=str(target),reused=True))
                    continue
                argv=(command(job,folder) if command else
                    ['Rscript',str(job['renderer']),str(folder/'figure-data.jsonl'),str(folder),'--jsonl'])
                log=(attempt/'renderer.log').open('wb');logs.append(log)
                future=pool.submit(execute,argv,log)
                pending.append((job,attempt,folder,target,binding,future))
            futures=[row[-1] for row in pending]
            while any(not future.done() for future in futures):wait(futures,timeout=.05)
            failed=False
            for job,attempt,folder,target,binding,future in pending:
                process=future.result()
                if process.returncode:
                    publish_json(attempt/'status.json',dict(state='failed',exit_code=process.returncode,target=str(target)))
                    failed=True;continue
                job['check'](folder)
                publish_json(folder/'render-receipt.json',dict(schema='genosapiens.track2-figure-job/v1',binding=binding,outputs=_inventory(folder)))
                target.parent.mkdir(parents=True,mode=0o700,exist_ok=True)
                if target.exists():
                    _inventory(target)  # Reject links before preserving old artifacts.
                    target.rename(attempt/'previous')
                folder.rename(target)
                publish_json(attempt/'status.json',dict(state='complete',target=str(target)))
                outcomes.append(dict(target=str(target),reused=False))
            _stop(active)
            active.clear()
            for log in logs:log.close()
            logs.clear()
            require(not failed,'A figure renderer failed. Its private attempt remains incomplete.')
        publish_json(invocation,dict(operational,state='complete',jobs=outcomes))
        return dict(effective_workers=effective,jobs=outcomes)
    finally:
        with lock:
            cancelled.set()
            processes=list(active)
        _stop(processes)
        pool.shutdown(wait=True,cancel_futures=True)
        for log in logs:log.close()
        for sig,handler in handlers.items():signal.signal(sig,handler)
