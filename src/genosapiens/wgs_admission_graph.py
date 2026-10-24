"""Bounded process execution of the unchanged categorical support kernel.

Workers receive only a target's source fields and bounded context jobs. They
never load a resource bundle, inherit a case graph, or own publication state.
Pickle is private parent/child transport, never an accepted file format.
"""
from collections import defaultdict
from dataclasses import dataclass
from datetime import datetime, timezone
import hashlib
import io
import json
import multiprocessing
import os
from pathlib import Path
import pickle
import queue
import sys

from genosapiens.wgs_failures import (WgsFailure, static_failure, worker_failure_message,
    checked_worker_failure, run_cleanup, failure_boundary)
import time


WORKERS = 10
MAX_JOB_CONTEXTS = 64
MAX_JOB_BYTES = 8 * 1024**2
MAX_RESULT_BYTES = 64 * 1024**2
_SOURCE_FIELDS = ('source_dataset', 'source_disease_id', 'is_ontology_propagated',
                  'disease_specificity', 'association_scope', 'publications')


class GraphSupportError(WgsFailure):
    """Only explicit static fields survive; supplied messages are discarded."""
    def __init__(self, message=None, *, category="unexpected", operation="graph_result_validate",
                 exception_kind="project_failure"):
        super().__init__("bounded_products", category, operation=operation, exception_kind=exception_kind)


def live_worker_count(known_items=None):
    """Size execution only; subsequent jobs still consume every source row."""
    if known_items is None:
        return WORKERS
    if type(known_items) is not int or known_items < 0:
        raise GraphSupportError(category='malformed_evidence', operation='graph_job_encode')
    return min(WORKERS, max(1, known_items))


class _LimitedBuffer(io.BytesIO):
    def __init__(self, limit, operation='graph_job_encode'):
        super().__init__()
        self.limit = limit
        self.operation = operation

    def write(self, data):
        if self.tell() + len(data) > self.limit:
            raise GraphSupportError(category='response_capacity', operation=self.operation)
        return super().write(data)


def _pack(value, limit, *, operation='graph_job_encode'):
    with failure_boundary('bounded_products', operation), _LimitedBuffer(limit, operation) as stream:
        pickle.Pickler(stream, protocol=5).dump(value)
        return stream.getvalue()


@dataclass(frozen=True, slots=True)
class _Record:
    source_fields: dict


@dataclass(frozen=True, slots=True)
class _Store:
    records: dict
    views: dict


@dataclass(frozen=True, slots=True)
class _Evidence:
    evidence_id: str
    public_manifest_sha256: str
    store: _Store


def _project_sources(indexed, reverse_views):
    result = []
    for public, assertions in indexed:
        records, views = {}, defaultdict(set)
        for assertion in assertions:
            key = assertion.record_id
            fields = public.store.records[key].source_fields
            records[key] = _Record({k: tuple(fields[k]) if k == 'publications'
                and fields[k] is not None else fields[k]
                for k in _SOURCE_FIELDS if k in fields})
            for disease in reverse_views[public.evidence_id].get(key, ()):
                views[disease].add(key)
        result.append((_Evidence(public.evidence_id, public.public_manifest_sha256,
            _Store(records, {k: frozenset(v) for k, v in views.items()})), tuple(assertions)))
    return tuple(result)


def _jobs(candidates, support_index, public_by_id, reach_cache):
    """Group by target, retaining every assertion and context alternative."""
    from genosapiens.disease_gene_evidence import assertion_native_disease_id
    reverse_views = {}
    for public in public_by_id.values():
        reverse = defaultdict(list)
        for disease, records in public.store.views.items():
            for record in records:
                reverse[record].append(disease)
        reverse_views[public.evidence_id] = reverse
    targets = defaultdict(list)
    for key, (context, namespace) in candidates.items():
        targets[namespace, context.locus_id or context.gene_id].append(key)
    for target in sorted(targets):
        indexed = _project_sources(tuple((public_by_id[p], assertions)
            for p, assertions in support_index[target].items()), reverse_views)
        native_ids = {assertion_native_disease_id(public.store, a)
                      for public, assertions in indexed for a in assertions}
        batch = []

        def encode(rows):
            try:
                yield _pack((indexed, tuple(rows)), MAX_JOB_BYTES)
            except GraphSupportError:
                if len(rows) == 1:
                    raise
                middle = len(rows) // 2
                yield from encode(rows[:middle])
                yield from encode(rows[middle:])

        for key in sorted(targets[target]):
            context, namespace = candidates[key]
            reached = reach_cache[context.native_disease_id]
            # _support only queries native IDs of this exact target's assertions.
            relevant = {native: reached[native] for native in native_ids if native in reached}
            batch.append((key, context, namespace, relevant))
            if len(batch) == MAX_JOB_CONTEXTS:
                yield from encode(batch)
                batch = []
        if batch:
            yield from encode(batch)


def _job_packet(packet):
    if (type(packet) is not tuple or len(packet) != 2 or type(packet[0]) is not int
            or packet[0] <= 0 or type(packet[1]) is not bytes):
        raise GraphSupportError(category='malformed_evidence', operation='graph_worker_execute')
    return packet


def _worker(index, jobs, results):
    from genosapiens.candidate_admission import _support
    job = digest = None
    operation = 'graph_worker_execute'
    try:
        results.put(('ready', index, os.getpid()))
        while True:
            idle_started=time.monotonic_ns()
            packet = jobs.get()
            idle=time.monotonic_ns()-idle_started
            if packet is None:
                results.put(('stopped', index))
                return
            job, payload = _job_packet(packet)
            digest = hashlib.sha256(payload).hexdigest()
            operation = 'graph_worker_execute'
            if type(payload) is not bytes or len(payload) > MAX_JOB_BYTES:
                raise GraphSupportError('graph input is invalid')
            indexed, rows = pickle.loads(payload)
            if type(rows) is not tuple or not 0 < len(rows) <= MAX_JOB_CONTEXTS:
                raise GraphSupportError('graph input is invalid')
            started = time.process_time_ns()
            wall_started=time.monotonic_ns()
            output = tuple((key, _support(context, namespace, indexed, reached))
                           for key, context, namespace, reached in rows)
            operation = 'graph_result_validate'
            message = _pack(output, MAX_RESULT_BYTES, operation=operation)
            results.put(('result', index, hashlib.sha256(payload).hexdigest(),
                         message, time.process_time_ns() - started,time.monotonic_ns()-wall_started,idle,job))
            del payload, indexed, rows, output, message
    except BaseException as error:
        results.put(worker_failure_message(error, worker=index, job=job, digest=digest, operation=operation))


class _Pool:
    """One outstanding job per worker; receive credits include queued results."""
    def __init__(self, *, startup_seconds=300, job_seconds=7200, shutdown_seconds=30,
                 worker_target=None, worker_count=WORKERS):
        if type(worker_count) is not int or not 1 <= worker_count <= WORKERS:
            raise GraphSupportError(category="malformed_evidence", operation="graph_job_encode")
        self.worker_count = worker_count
        self.startup_seconds = min(300, startup_seconds)
        self.job_seconds = min(7200, job_seconds)
        self.shutdown_seconds = min(30, shutdown_seconds)
        if min(self.startup_seconds, self.job_seconds, self.shutdown_seconds) <= 0:
            raise GraphSupportError(category='malformed_evidence', operation='graph_job_encode')
        ctx = multiprocessing.get_context('spawn')
        self.results = ctx.Queue(maxsize=self.worker_count)
        self.jobs = [ctx.Queue(maxsize=1) for _ in range(self.worker_count)]
        self.workers = [ctx.Process(target=_worker if worker_target is None else worker_target,
                                   args=(i, self.jobs[i], self.results),
                                   name=f'wgs-graph-{i}') for i in range(self.worker_count)]
        self.pending = {}
        self.job_ids = {}
        self.next_job_id = 0
        self.worker_operation = "graph_worker_execute"
        self.result_operation = "graph_result_validate"
        self.error_operations = frozenset({"graph_worker_execute", "graph_result_validate", "graph_job_encode"})
        self.metrics = dict(worker_count=self.worker_count, worker_pids=[], worker_contexts=[0] * self.worker_count,
            max_outstanding=0, max_job_bytes=0, max_result_bytes=0,
            max_reserved_input_bytes=0, max_reserved_result_bytes=0,
            cpu_nanoseconds=0, completed_jobs=0,submitted_jobs=0,received_jobs=0,cancelled_jobs=0,
            outstanding_jobs=0,outstanding_rows=0,outstanding_input_bytes=0,reserved_result_bytes=0,
            worker_busy_nanoseconds=[0]*self.worker_count,worker_idle_nanoseconds=[0]*self.worker_count,
            worker_cpu_nanoseconds=[0]*self.worker_count)

    def _reservations(self):
        self.metrics.update(outstanding_jobs=len(self.pending),
            outstanding_rows=sum(len(row[1]) for row in self.pending.values()),
            outstanding_input_bytes=sum(row[3] for row in self.pending.values()),
            reserved_result_bytes=len(self.pending)*MAX_RESULT_BYTES)

    def _received(self,index,cpu,wall,idle):
        self.metrics['received_jobs']+=1
        self.metrics['worker_cpu_nanoseconds'][index]+=cpu
        self.metrics['worker_busy_nanoseconds'][index]+=wall
        self.metrics['worker_idle_nanoseconds'][index]+=idle
        self._reservations()

    def _failure(self, category, *, index=None, operation=None, kind='project_failure'):
        operation = self.result_operation if operation is None else operation
        pending = self.pending.get(index)
        self.metrics['failure'] = dict(worker_index=index, job_number=self.job_ids.get(index),
            input_bytes=None if pending is None else pending[3],
            category=category, operation=operation, exception_kind=kind)
        return GraphSupportError(category=category, operation=operation, exception_kind=kind)

    def _receive(self, deadline, *, stopping=False):
        while time.monotonic() < deadline:
            heartbeat=getattr(self,'progress_callback',None)
            if heartbeat is not None:heartbeat()
            try:
                message = self.results.get(timeout=min(.1, max(.001, deadline-time.monotonic())))
            except queue.Empty:
                dead = next((i for i,w in enumerate(self.workers) if w.exitcode is not None
                    and (not stopping or w.exitcode != 0)), None)
                if dead is not None:
                    raise self._failure('interrupted', index=dead, operation=self.worker_operation)
                continue
            if (type(message) is not tuple or len(message) < 2 or type(message[1]) is not int
                    or not 0 <= message[1] < self.worker_count):
                raise self._failure('malformed_evidence')
            if message[0] in ('error', 'result'):
                index = message[1]
                if index not in self.pending or index not in self.job_ids:
                    raise self._failure('authentication', index=index)
                if message[0] == 'error':
                    try:
                        failure = checked_worker_failure(message, worker=index, job=self.job_ids[index],
                            digest=self.pending[index][0], operations=self.error_operations)
                    except ValueError:
                        raise self._failure('malformed_evidence') from None
                    raise self._failure(failure.category, index=index, operation=failure.operation,
                                        kind=failure.exception_kind) from None
                if len(message) != 8 or type(message[7]) is not int or message[7] != self.job_ids[index]:
                    raise self._failure('authentication', index=index)
                message = message[:7]
            return message
        raise self._failure('timeout', index=min(self.pending, key=lambda i: self.pending[i][2], default=None),
                            operation=self.worker_operation)

    def start(self):
        deadline = time.monotonic() + self.startup_seconds
        for worker in self.workers:
            worker.start()
        ready = set()
        while len(ready) < self.worker_count:
            message = self._receive(deadline)
            if (len(message) != 3 or message[0] != 'ready' or message[1] in ready
                    or message[2] != self.workers[message[1]].pid):
                raise GraphSupportError(category='malformed_evidence', operation=self.result_operation)
            ready.add(message[1])
        self.metrics['worker_pids'] = [w.pid for w in self.workers]

    def submit(self, index, payload, *, keys=None):
        if index in self.pending or not 0 <= index < self.worker_count or type(payload) is not bytes or len(payload) > MAX_JOB_BYTES:
            raise GraphSupportError(category='response_capacity', operation=self.result_operation)
        if keys is None:
            _, rows = pickle.loads(payload)
            keys = tuple(row[0] for row in rows)
        if not 0 < len(keys) <= MAX_JOB_CONTEXTS or len(set(keys)) != len(keys):
            raise GraphSupportError(category='malformed_evidence', operation=self.result_operation)
        self.pending[index] = (hashlib.sha256(payload).hexdigest(), keys,
                               time.monotonic()+self.job_seconds, len(payload))
        self.next_job_id += 1
        self.job_ids[index] = self.next_job_id
        self.jobs[index].put((self.job_ids[index], payload), timeout=1)
        self.metrics['submitted_jobs']+=1
        self._reservations()
        self.metrics['max_outstanding'] = max(self.metrics['max_outstanding'], len(self.pending))
        self.metrics['max_job_bytes'] = max(self.metrics['max_job_bytes'], len(payload))
        self.metrics['max_reserved_input_bytes'] = max(self.metrics['max_reserved_input_bytes'],
            sum(row[3] for row in self.pending.values()))
        self.metrics['max_reserved_result_bytes'] = max(self.metrics['max_reserved_result_bytes'],
            len(self.pending) * MAX_RESULT_BYTES)

    @failure_boundary('bounded_products', 'graph_result_validate',
                      known_errors=(pickle.UnpicklingError, ValueError, TypeError, EOFError))
    def receive(self):
        from genosapiens.candidate_admission import DiseaseSupportEntry
        message = self._receive(min(row[2] for row in self.pending.values()))
        if len(message) != 7 or message[0] != 'result' or message[1] not in self.pending:
            raise GraphSupportError(category='malformed_evidence', operation=self.result_operation)
        _, index, digest, content, cpu,wall,idle = message
        expected, keys, _, _ = self.pending[index]
        if (digest != expected or type(content) is not bytes or len(content) > MAX_RESULT_BYTES
                or any(type(value) is not int or value<0 for value in (cpu,wall,idle))):
            raise GraphSupportError(category='malformed_evidence', operation=self.result_operation)
        rows = pickle.loads(content)
        if (type(rows) is not tuple or len(rows) != len(keys)
                or any(type(row) is not tuple or len(row) != 2 or row[0] != key
                    or type(row[1]) is not tuple or any(type(e) is not DiseaseSupportEntry for e in row[1])
                    for row, key in zip(rows, keys))):
            raise GraphSupportError(category='malformed_evidence', operation=self.result_operation)
        self.metrics['max_result_bytes'] = max(self.metrics['max_result_bytes'], len(content))
        self.metrics['worker_contexts'][index] += len(rows)
        self.metrics['cpu_nanoseconds'] += cpu
        self.metrics['completed_jobs'] += 1
        # The caller consumes this result before submitting another job.
        del self.pending[index]
        self._received(index,cpu,wall,idle)
        return index, rows

    def finish(self):
        if self.pending:
            raise GraphSupportError(category='malformed_evidence', operation=self.result_operation)
        deadline = time.monotonic()+self.shutdown_seconds
        for jobs in self.jobs:
            jobs.put(None, timeout=1)
        stopped = set()
        while len(stopped) < self.worker_count:
            message = self._receive(deadline, stopping=True)
            if len(message) != 2 or message[0] != 'stopped' or message[1] in stopped:
                raise GraphSupportError(category='malformed_evidence', operation=self.result_operation)
            stopped.add(message[1])
        for worker in self.workers:
            worker.join(timeout=max(0, deadline-time.monotonic()))
            if worker.exitcode != 0:
                raise GraphSupportError(category='malformed_evidence', operation=self.result_operation)
        try:
            self.results.get_nowait()
        except queue.Empty:
            return
        raise GraphSupportError(category='malformed_evidence', operation=self.result_operation)

    def close(self):
        # Workers run pure Python and launch no descendants. Use a shared
        # cleanup deadline, including interruption and partial startup failures.
        deadline = time.monotonic()+self.shutdown_seconds
        def terminate(worker):
            if worker.pid is not None and worker.is_alive(): worker.terminate()
        def join(worker, divisor=1):
            if worker.pid is not None: worker.join(timeout=max(0, (deadline-time.monotonic())/divisor))
        def kill(worker):
            if worker.pid is not None and worker.is_alive(): worker.kill()
        def stopped(worker):
            if worker.pid is not None and worker.is_alive():
                raise GraphSupportError(category='interrupted', operation='worker_shutdown')
        actions = []
        for action in (terminate, lambda w: join(w, 2), kill, join, stopped):
            actions.extend(('worker_shutdown', lambda w=w, action=action: action(w)) for w in self.workers)
        for channel in (*self.jobs, self.results):
            actions.extend([('worker_shutdown', channel.cancel_join_thread), ('worker_shutdown', channel.close)])
        run_cleanup(actions, primary=sys.exc_info()[1])
        self.job_ids.clear()
        self.metrics['cancelled_jobs']+=len(self.pending)
        self.pending.clear()
        self._reservations()


from genosapiens.wgs_runtime import analysis_progress, analysis_worker_metrics


class GraphProgress:
    def __init__(self, path, total, progress=None):
        self.progress = progress
        self.path = None if path is None else Path(path)
        self.started = time.monotonic()
        self.last = 0
        self.total = total
        self.completed = 0
        self.metrics = {}

    def write(self, phase, *, force=False):
        analysis_progress(self.progress, 'context_assembly', 'graph_support',
            counts={'graph_completed':self.completed, 'graph_total':self.total},
            worker=analysis_worker_metrics(self.metrics), unit='contexts', completed=self.completed, total=self.total)
        now = time.monotonic()
        if self.path is None or not force and now-self.last < 1:
            return
        record = dict(schema='genosapiens.graph-support-progress', version=1,
            at_utc=datetime.now(timezone.utc).isoformat(), phase=phase,
            total_contexts=self.total, completed_contexts=self.completed,
            elapsed_seconds=now-self.started, workers=self.metrics.get("worker_count", 0), **self.metrics)
        temporary = self.path.with_suffix('.tmp')
        fd = os.open(temporary, os.O_WRONLY | os.O_CREAT | os.O_TRUNC | os.O_NOFOLLOW, 0o600)
        with os.fdopen(fd, 'w') as stream:
            json.dump(record, stream, sort_keys=True)
            stream.write('\n')
        os.replace(temporary, self.path)
        self.last = now


def parallel_support(candidates, support_index, public_by_id, reach_cache, *, progress):
    pool = _Pool(worker_count=live_worker_count(len(candidates)))
    pool.progress_callback = lambda: progress.write('support')
    jobs = _jobs(candidates, support_index, public_by_id, reach_cache)
    progress.write('worker_startup', force=True)
    try:
        pool.start()
        progress.metrics = pool.metrics
        for index in range(pool.worker_count):
            payload = next(jobs, None)
            if payload is None:
                break
            pool.submit(index, payload)
        progress.write('support', force=True)
        while pool.pending:
            index, rows = pool.receive()
            yield from rows
            progress.completed += len(rows)
            progress.write('support')
            payload = next(jobs, None)
            if payload is not None:
                pool.submit(index, payload)
        if progress.completed != len(candidates):
            raise GraphSupportError(category='malformed_evidence', operation=pool.result_operation)
        pool.finish()
        progress.write('assembling_graph', force=True)
    except BaseException as error:
        primary = error if isinstance(error, (KeyboardInterrupt, GeneratorExit)) else static_failure(error, 'bounded_products', operation='graph_worker_execute')
        run_cleanup([('diagnostic_publish', lambda: progress.write('failed', force=True))], primary=primary)
        raise primary from None
    finally:
        run_cleanup([('worker_shutdown', jobs.close), ('worker_shutdown', pool.close)], primary=sys.exc_info()[1])
