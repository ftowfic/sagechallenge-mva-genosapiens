"""Persistent Linux workers with one bounded request per process.

The parent freezes prepared inputs before forking. Initializers replace all
used database connections and query sessions before accepting work. Pipes carry
only task descriptors and bounded completion metadata, never resource graphs.
"""
from collections import deque
from contextlib import ExitStack
import json
import gc
import multiprocessing
from multiprocessing.connection import wait
import os
import signal
import time
from pathlib import Path


def process_metrics():
    import resource
    usage = resource.getrusage(resource.RUSAGE_SELF)
    result = dict(user_seconds=usage.ru_utime, system_seconds=usage.ru_stime,
        peak_rss_bytes=usage.ru_maxrss*1024)
    try:
        memory = {row.split(':')[0]: int(row.split()[1])*1024
            for row in Path('/proc/self/smaps_rollup').read_text().splitlines()[1:] if ':' in row}
        result['private_bytes'] = memory.get('Private_Clean', 0)+memory.get('Private_Dirty', 0)
        result['pss_bytes'] = memory.get('Pss')
        result['io'] = {key: int(value) for key, value in
            (row.split(':') for row in Path('/proc/self/io').read_text().splitlines())}
    except (OSError, ValueError):
        result['private_memory_and_io_available'] = False
    return result

MAX_MESSAGE_BYTES = 64 * 1024


def open_worker_query_cache(path, *, binding, clinical, execution_config):
    """Bind a new empty cache to this worker's independently owned raw records."""
    from genosapiens.wgs_allele_query_cache import CaseAlleleQueryCache, clinical_source_bindings
    from genosapiens.wgs_ledger import canonical
    parents = json.loads(binding)
    parents['clinical_sources'] = clinical_source_bindings(clinical)
    return CaseAlleleQueryCache(path, binding=canonical(parents), execution_config=execution_config)


def open_worker_reader(reader, execution_config):
    from genosapiens.wgs_admission_projection import AdmissionProjectionReader
    if type(reader) is not AdmissionProjectionReader:
        raise ValueError('scientific worker requires an authenticated admission reader')
    old = reader._owner
    arguments = dict(expected_binding=old.binding, expected_completion_sha256=old.completion.sha256,
        execution_config=execution_config)
    path = old.path
    old.close()
    return AdmissionProjectionReader(path, **arguments)


def _send(channel, value):
    raw = json.dumps(value, separators=(',', ':'), allow_nan=False).encode()
    if len(raw) > MAX_MESSAGE_BYTES:
        raise ValueError('scientific worker message exceeds capacity')
    channel.send_bytes(raw)


def _receive(channel):
    return json.loads(channel.recv_bytes(MAX_MESSAGE_BYTES))


def _worker(channel, inherited_channels, initialize):
    # Own the process group so bounded cancellation also stops query children.
    os.setsid()
    for other in inherited_channels:
        other.close()
    signal.setitimer(signal.ITIMER_REAL, 0)
    def cancelled(*_):
        raise SystemExit(143)
    signal.signal(signal.SIGTERM, cancelled)
    try:
        with ExitStack() as stack:
            started = time.monotonic()
            evaluate = initialize(stack)
            _send(channel, dict(kind='ready', pid=os.getpid(), seconds=time.monotonic()-started,
                phases=getattr(evaluate, 'initialization_metrics', {}), resources=process_metrics()))
            while True:
                job = _receive(channel)
                if job is None:
                    break
                result = evaluate(job)
                result['worker_resources'] = process_metrics()
                _send(channel, dict(kind='result', job=job, result=result))
    except BaseException as error:
        from genosapiens.wgs_science_failures import observe_setup_failure
        observe_setup_failure(error)
        try:
            _send(channel, dict(kind='fatal', error_type=type(error).__name__))
        except (OSError, ValueError):
            pass
    finally:
        channel.close()


class PersistentScienceWorkers:
    """Bound active tasks and reorder storage to the configured worker count."""
    def __init__(self, workers, initialize, *, startup_timeout=600):
        if type(workers) is not int or not 1 <= workers <= 8:
            raise ValueError('invalid scientific worker count')
        self.channels, self.processes = [], []
        self.pending, self.results = deque(), {}
        self.jobs, self.deadlines = {}, {}
        self.fatal = None
        started = time.monotonic()
        self.metrics = dict(workers=workers, startup_seconds=[], startup_resources=[], max_outstanding=0,
            startup_phases=[], max_result_bytes=0, wait_seconds=0.0, start_method='fork_read_only_resources')
        ctx = multiprocessing.get_context('fork')
        freeze = gc.get_freeze_count() == 0
        try:
            # Avoid a child's cyclic collector dirtying every shared resource
            # page. The parent restores its ordinary collector after spawning.
            if freeze:
                gc.collect()
                gc.freeze()
            for _ in range(workers):
                parent, child = ctx.Pipe()
                process = ctx.Process(target=_worker, args=(child, (*self.channels, parent), initialize))
                process.start()
                child.close()
                self.channels.append(parent)
                self.processes.append(process)
            if freeze:
                gc.unfreeze()
            remaining = set(self.channels)
            deadline = time.monotonic() + startup_timeout
            while remaining:
                ready = wait(list(remaining), timeout=max(0, min(0.2, deadline-time.monotonic())))
                if time.monotonic() >= deadline or any(p.exitcode is not None for p in self.processes):
                    raise RuntimeError('scientific worker initialization failed')
                for channel in ready:
                    message = _receive(channel)
                    if message.get('kind') != 'ready':
                        raise RuntimeError('scientific worker initialization failed')
                    self.metrics['startup_seconds'].append(message['seconds'])
                    self.metrics['startup_phases'].append(message['phases'])
                    self.metrics['startup_resources'].append(message['resources'])
                    remaining.remove(channel)
            self.metrics['startup_wall_seconds'] = time.monotonic()-started
        except BaseException:
            if freeze:
                gc.unfreeze()
            self.close()
            raise

    def submit(self, job, *, timeout):
        self.poll()
        if self.fatal is not None or len(self.pending) >= len(self.channels):
            raise RuntimeError('scientific worker queue unavailable')
        channel = next(c for c in self.channels if c not in self.jobs)
        key = job['ordinal']
        if key in self.pending:
            raise ValueError('duplicate scientific worker job')
        self.pending.append(key)
        self.jobs[channel] = job
        self.deadlines[channel] = time.monotonic() + timeout
        _send(channel, job)
        self.metrics['max_outstanding'] = max(self.metrics['max_outstanding'], len(self.pending))

    def poll(self, timeout=0):
        active = [c for c in self.jobs if self.jobs[c]['ordinal'] not in self.results]
        for channel in wait(active, timeout=timeout) if active else ():
            try:
                message = _receive(channel)
                if message.get('kind') != 'result' or message['job'] != self.jobs[channel]:
                    raise ValueError('scientific worker reply differs')
                result = message['result']
                self.results[message['job']['ordinal']] = result
                self.metrics['max_result_bytes'] = max(self.metrics['max_result_bytes'], len(json.dumps(message)))
                if result.get('failure', {}).get('scope') in ('shared', 'unknown'):
                    self.fatal = 'worker_failure'
            except (EOFError, OSError, ValueError, KeyError):
                self.fatal = 'worker_exit_or_invalid_reply'
        for channel in active:
            if self.jobs[channel]['ordinal'] in self.results:
                continue
            if time.monotonic() >= self.deadlines[channel]:
                self.fatal = 'worker_timeout'
        if any(p.exitcode is not None for p in self.processes):
            self.fatal = self.fatal or 'worker_exit'

    def take(self):
        key = self.pending[0]
        started = time.monotonic()
        while key not in self.results and self.fatal is None:
            self.poll(0.1)
        self.metrics['wait_seconds'] += time.monotonic()-started
        if key not in self.results:
            raise RuntimeError(self.fatal)
        self.pending.popleft()
        channel = next(c for c, job in self.jobs.items() if job['ordinal'] == key)
        del self.jobs[channel], self.deadlines[channel]
        return key, self.results.pop(key)

    def close(self):
        for channel in self.channels:
            try:
                _send(channel, None)
            except (OSError, ValueError):
                pass
        stop = time.monotonic()+2
        for process in self.processes:
            process.join(max(0, stop-time.monotonic()))
        descendants = {}
        pending = [p.pid for p in self.processes if p.is_alive()]
        while pending and len(descendants) < 512:
            pid = pending.pop()
            try:
                children = Path(f'/proc/{pid}/task/{pid}/children').read_text().split()
            except OSError:
                continue
            for child in children:
                child = int(child)
                if child in descendants:
                    continue
                try:
                    signature = Path(f'/proc/{child}/stat').read_text().rsplit(')', 1)[1].split()[19]
                except (OSError, IndexError):
                    continue
                descendants[child] = signature
                pending.append(child)
        def signal_descendants(signum):
            for pid, expected in descendants.items():
                try:
                    actual = Path(f'/proc/{pid}/stat').read_text().rsplit(')', 1)[1].split()[19]
                    if actual == expected:
                        os.kill(pid, signum)
                except (OSError, IndexError):
                    pass
        signal_descendants(signal.SIGTERM)
        for process in self.processes:
            if process.is_alive():
                try:
                    os.killpg(process.pid, signal.SIGTERM)
                except ProcessLookupError:
                    pass
        stop = time.monotonic()+3
        for process in self.processes:
            process.join(max(0, stop-time.monotonic()))
            if process.is_alive():
                try:
                    os.killpg(process.pid, signal.SIGKILL)
                except ProcessLookupError:
                    pass
                process.join(1)
        signal_descendants(signal.SIGKILL)
        for channel in self.channels:
            channel.close()

    def __enter__(self):
        return self

    def __exit__(self, *_):
        self.close()
