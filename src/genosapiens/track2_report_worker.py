"""Invocation-owned, bounded calibration subprocesses on Linux."""
import json
import math
import os
from pathlib import Path
import resource
import signal
import subprocess
import sys
import threading
import time
import uuid

from genosapiens.track2_operator_inputs import require
from genosapiens.track2_scientific_storage import publish_json


def limits(*, cpus=2, memory_bytes=4*1024**3, deadline_seconds=None):
    require(type(cpus) is int and 1 <= cpus <= 64, 'Invalid calibration CPU limit.')
    require(type(memory_bytes) is int and memory_bytes >= 128*1024**2,
            'Invalid calibration memory limit.')
    require(deadline_seconds is None or isinstance(deadline_seconds,(float,int)) and math.isfinite(deadline_seconds)
            and deadline_seconds > 0, 'Invalid calibration deadline.')
    return dict(cpus=cpus,memory_bytes=memory_bytes,deadline_seconds=deadline_seconds)


def _terminate(worker):
    # This process group belongs only to this invocation. Kill descendants even
    # if the immediate worker already exited or ignored graceful termination.
    try:os.killpg(worker.pid,signal.SIGTERM)
    except ProcessLookupError:return
    try:worker.wait(timeout=1)
    except subprocess.TimeoutExpired:pass
    try:os.killpg(worker.pid,signal.SIGKILL)
    except ProcessLookupError:pass
    worker.wait()


def run_owned_worker(command, attempt, *, controls, cancelled=None):
    attempt=Path(attempt);start=time.monotonic();worker=None
    publish_json(attempt/'worker-started.json',dict(status='incomplete',requested=controls))
    reason=None
    previous=None
    if threading.current_thread() is threading.main_thread():
        previous=signal.getsignal(signal.SIGTERM)
        def terminated(signum, frame):
            raise InterruptedError('The reporting invocation was terminated.')
        signal.signal(signal.SIGTERM,terminated)
    try:
        require(cancelled is None or not cancelled(), 'Calibration was cancelled.')
        with (attempt/'worker.log').open('xb') as log:
            os.chmod(log.name,0o600)
            worker=subprocess.Popen(command,stdout=log,stderr=subprocess.STDOUT,start_new_session=True)
            while worker.poll() is None:
                if cancelled is not None and cancelled():
                    reason='cancelled';raise InterruptedError('Calibration was cancelled.')
                if controls['deadline_seconds'] is not None and time.monotonic()-start >= controls['deadline_seconds']:
                    reason='deadline';raise TimeoutError('Calibration exceeded its invocation deadline.')
                time.sleep(0.05)
            require(worker.returncode==0, 'The calibration worker failed. Inspect its retained operation log.')
    except BaseException as error:
        if worker is not None:_terminate(worker)
        publish_json(attempt/'worker-outcome.json',dict(status='incomplete',reason=reason or type(error).__name__,
            elapsed_seconds=time.monotonic()-start,returncode=worker.returncode if worker else None))
        raise
    finally:
        if worker is not None:_terminate(worker)
        if previous is not None:signal.signal(signal.SIGTERM,previous)
    publish_json(attempt/'worker-outcome.json',dict(status='complete',elapsed_seconds=time.monotonic()-start,
        returncode=worker.returncode))


def prepare_bounded_calibration(operation, inputs, *, classification, source_manifest_sha256,
                                controls=None, cancelled=None, retained_calibration=()):
    operation=Path(operation).absolute()
    require(operation.resolve()==operation, 'Calibration operation links are unsupported.')
    operation.mkdir(mode=0o700,parents=True,exist_ok=True)
    require(operation.stat().st_uid==os.geteuid() and not operation.stat().st_mode & 0o077,
            'Calibration state requires a private owned directory.')
    attempts=operation/'workers';attempts.mkdir(mode=0o700,exist_ok=True)
    require(attempts.resolve()==attempts, 'Calibration worker links are unsupported.')
    attempt=attempts/uuid.uuid4().hex;attempt.mkdir(mode=0o700)
    controls=limits(**(controls or {}))
    publish_json(attempt/'request.json',dict(operation=str(operation),inputs=inputs,classification=classification,
        source_manifest_sha256=source_manifest_sha256,controls=controls,retained_calibration=retained_calibration))
    run_owned_worker([sys.executable,'-m','genosapiens.track2_report_worker',str(attempt)],attempt,
        controls=controls,cancelled=cancelled)
    result=json.loads((attempt/'result.json').read_bytes())
    result['bank']=Path(result['bank']);result['resource_contract']=Path(result['resource_contract'])
    return result


def main():
    attempt=Path(sys.argv[1]);request=json.loads((attempt/'request.json').read_bytes())
    controls=limits(**request.pop('controls'))
    available=sorted(os.sched_getaffinity(0));effective=available[:controls['cpus']]
    os.sched_setaffinity(0,effective)
    resource.setrlimit(resource.RLIMIT_AS,(controls['memory_bytes'],controls['memory_bytes']))
    # CPU time is an additional fail-safe. The parent enforces elapsed time.
    cpu_seconds=None if controls['deadline_seconds'] is None else max(1,math.ceil(controls['deadline_seconds']*len(effective)))
    if cpu_seconds is not None:resource.setrlimit(resource.RLIMIT_CPU,(cpu_seconds,cpu_seconds))
    publish_json(attempt/'effective-limits.json',dict(cpu_affinity=effective,
        memory_bytes=resource.getrlimit(resource.RLIMIT_AS)[0],cpu_seconds=cpu_seconds,
        deadline_seconds=controls['deadline_seconds']))
    if request.pop('kind','calibration')=='render':
        script=Path(__file__).resolve().parents[2]/'scripts/render_track2_challenge_overviews.R'
        subprocess.run(['Rscript',str(script),request['data'],request['destination']],check=True)
        result=dict(operations=dict(status='complete'))
    else:
        from genosapiens.track2_report_calibration_stage import prepare_report_calibration
        result=prepare_report_calibration(**request)
        result['bank']=str(result['bank']);result['resource_contract']=str(result['resource_contract'])
    result['operations']['peak_memory_bytes']=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss*1024
    result['operations']['requested_limits']=controls
    result['operations']['effective_cpu_count']=len(effective)
    publish_json(attempt/'result.json',result)


def bounded_renderer(operation, *, controls=None, cancelled=None):
    """Return the existing renderer callback with invocation-owned supervision."""
    operation=Path(operation).absolute();controls=limits(**(controls or {}))
    def render(script,data,destination):
        expected=Path(__file__).resolve().parents[2]/'scripts/render_track2_challenge_overviews.R'
        require(Path(script)==expected, 'An unsupported overview renderer was requested.')
        require(operation.resolve()==operation, 'Figure operation links are unsupported.')
        operation.mkdir(mode=0o700,parents=True,exist_ok=True)
        attempt=operation/uuid.uuid4().hex;attempt.mkdir(mode=0o700)
        publish_json(attempt/'request.json',dict(kind='render',data=str(data),destination=str(destination),controls=controls))
        run_owned_worker([sys.executable,'-m','genosapiens.track2_report_worker',str(attempt)],attempt,
            controls=controls,cancelled=cancelled)
    return render


if __name__=='__main__':main()
