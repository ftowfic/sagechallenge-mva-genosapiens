"""Two-window screening handoff; sessions and SQLite stay in one thread.

This is a live operation, not a checkpoint or another computation pool. Only
immutable setup, decisions and witness bytes cross the thread boundary.
"""
from collections import OrderedDict
from contextlib import contextmanager
from dataclasses import dataclass
import queue
import sys
import threading
import time
from types import MappingProxyType

from genosapiens.wgs_analysis_eligibility import AnalysisEligibilitySession, AnalysisEligibilityError, EligibilityBatch
from genosapiens.wgs_query_limits import _cancellable_queries, _check_query_cancellation


@dataclass(frozen=True,slots=True)
class ScreeningSetup:
    screening: object
    source_sha256: str
    resource_bindings: tuple


class AdmissionScreeningThread:
    def __init__(self, factory, *, startup_seconds=300, shutdown_seconds=30):
        if not 0<startup_seconds<=300 or not 0<shutdown_seconds<=30:
            raise AnalysisEligibilityError('invalid screening thread deadline')
        self._jobs=queue.Queue(maxsize=2)
        self._results=queue.Queue(maxsize=2)
        self._cancelled=threading.Event()
        self._lock=threading.RLock()
        self._error=None
        self._pending=OrderedDict()
        self._window_sizes={}
        self._next=0
        self.closed=self.finished=False
        self.shutdown_seconds=shutdown_seconds
        self.metrics=dict(submitted_windows=0,screened_windows=0,committed_windows=0,
            submitted_inputs=0,screened_inputs=0,committed_inputs=0,max_windows=0,
            submitted_occurrences=0,screened_occurrences=0,committed_occurrences=0,
            logical_proof_references=0,window_witnesses=0,witness_bytes=0,max_witness_bytes=0,
            cancelled_windows=0,screening_cpu_nanoseconds=0,screening_wall_nanoseconds=0,
            public_query_batches=0,thread_id=None,thread_stopped=False,cleanup_failed=False,
            cancelled_states={'received':0,'screening':0,'ready_to_write':0})
        self.thread=threading.Thread(target=self._run,args=(factory,),name='wgs-admission-screening',daemon=True)
        self.thread.start()
        try:
            kind,self.setup=self._get(deadline=time.monotonic()+startup_seconds)
            if kind!='ready' or type(self.setup) is not ScreeningSetup:
                raise AnalysisEligibilityError('screening startup reply differs')
        except BaseException as error:
            self.fail(error)
            self.close()
            raise self._error

    def fail(self,error):
        with self._lock:
            if self._error is None:self._error=error
        self._cancelled.set()

    def check(self):
        with self._lock:error=self._error
        if error is not None:
            from genosapiens.wgs_admission_diagnostics import current, tag
            owner = current()
            # Coordinator failures also cancel this thread. Preserve their owner.
            if owner is None or owner.primary_error is not error:
                tag(error, 'admission_screening_receive')
            raise error

    def _get(self,*,deadline=None,check_wait=None):
        while True:
            self.check()
            if check_wait is not None:check_wait()
            if deadline is not None and time.monotonic()>=deadline:
                raise AnalysisEligibilityError('screening reply deadline exceeded')
            try:return self._results.get(timeout=.1)
            except queue.Empty:
                if not self.thread.is_alive():
                    self.check()
                    raise AnalysisEligibilityError('screening thread ended without closure')

    def _send(self,value):
        while not self._cancelled.is_set():
            try:
                self._results.put(value,timeout=.1)
                return
            except queue.Full:pass
        _check_query_cancellation()

    def _observe(self,owner):
        with self._lock:
            self.metrics['public_query_batches']=owner.public_query_batches
            if owner.case_cache is not None:self.metrics['cache']=dict(owner.case_cache.metrics)
            if owner.public_session is not None:
                self.metrics['population']=dict(owner.public_session._population.work_counts)

    @contextmanager
    def _preserve_origin(self,owner):
        try:yield
        except BaseException as error:
            self.fail(error)
            raise
        finally:self._observe(owner)

    def _run(self,factory):
        try:
            with _cancellable_queries(self._cancelled),factory() as owner,self._preserve_origin(owner):
                if type(owner) is not AnalysisEligibilitySession or owner.finished:
                    raise AnalysisEligibilityError('fresh screening thread owner is required')
                with self._lock:self.metrics['thread_id']=threading.get_ident()
                self._send(('ready',ScreeningSetup(owner.screening,owner.source_sha256,owner.resource_bindings)))
                while True:
                    _check_query_cancellation()
                    try:job=self._jobs.get(timeout=.1)
                    except queue.Empty:continue
                    if job[0]=='finish':
                        result=owner.finish(expected_inputs=job[1],expected_occurrences=job[2])
                        break
                    _,token,rows=job
                    with self._lock:
                        count,state=self._pending[token]
                        if state!='received':raise AnalysisEligibilityError('screening window state differs')
                        self._pending[token]=(count,'screening')
                    cpu=time.thread_time_ns();started=time.monotonic_ns()
                    try:batch=owner.evaluate_window(rows)
                    finally:
                        with self._lock:
                            self.metrics['screening_cpu_nanoseconds']+=time.thread_time_ns()-cpu
                            self.metrics['screening_wall_nanoseconds']+=time.monotonic_ns()-started
                    _check_query_cancellation()
                    witnesses=MappingProxyType(dict(owner.witnesses))
                    if type(batch) is not EligibilityBatch or any(type(value) is not bytes for value in witnesses.values()):
                        raise AnalysisEligibilityError('screening result ownership differs')
                    with self._lock:
                        self.metrics['screened_windows']+=1
                        self.metrics['screened_inputs']+=len(rows)
                        self.metrics['screened_occurrences']+=sum(len(row.source_allele_ids) for row in batch.decisions)
                        self.metrics['logical_proof_references']+=sum(len(row.proof_references) for row in batch.decisions)
                        self.metrics['window_witnesses']+=len(witnesses)
                        witness_bytes=sum(len(value) for value in witnesses.values())
                        self.metrics['witness_bytes']+=witness_bytes
                        self.metrics['max_witness_bytes']=max(self.metrics['max_witness_bytes'],witness_bytes)
                        self._observe(owner)
                        self._pending[token]=(count,'ready_to_write')
                    self._send(('screened',token,batch,witnesses))
                    del rows,batch,witnesses,job
            # Resource/cache closure must succeed before the completion reply.
            self._send(('finished',result))
        except BaseException as error:
            with self._lock:
                primary = self._error
            if primary is not None and primary is not error:
                from genosapiens.wgs_failures import run_cleanup
                from genosapiens.wgs_admission_diagnostics import failure
                self.cleanup_failure = failure('screening_close', error)
                run_cleanup((('resource_session_close', lambda: (_ for _ in ()).throw(error)),), primary=primary)
            self.fail(error)
        finally:
            with self._lock:self.metrics['thread_stopped']=True

    def submit(self,rows):
        from genosapiens.wgs_admission_projection import EncodedAdmissionRow
        self.check()
        if (self.closed or self.finished or type(rows) is not tuple or not 0<len(rows)<=64
                or any(not isinstance(row,EncodedAdmissionRow) for row in rows)
                or sum(row.input_size+row.proof_size+len(row.screening_facts) for row in rows)>64*1024**2):
            raise AnalysisEligibilityError('invalid screening submission')
        with self._lock:
            if len(self._pending)>=2:raise AnalysisEligibilityError('screening windows exceed capacity')
            token=self._next;self._next+=1
            self._pending[token]=(len(rows),'received')
            occurrences=sum(row.source_occurrences for row in rows)
            self._window_sizes[token]=(occurrences,sum(row.input_size+row.proof_size+len(row.screening_facts) for row in rows))
            self.metrics['submitted_windows']+=1;self.metrics['submitted_inputs']+=len(rows)
            self.metrics['submitted_occurrences']+=occurrences
            self.metrics['max_windows']=max(self.metrics['max_windows'],len(self._pending))
        self._jobs.put_nowait(('window',token,rows))
        return token

    def receive(self,*,check_wait=None):
        message=self._get(check_wait=check_wait)
        if (len(message)!=4 or message[0]!='screened' or not self._pending
                or message[1]!=next(iter(self._pending))):
            raise AnalysisEligibilityError('screening result order differs')
        return message[1:]

    def ready_to_write(self):
        """Whether the oldest held window can retire after both writer commits."""
        self.check()
        with self._lock:
            return bool(self._pending) and next(iter(self._pending.values()))[1]=='ready_to_write'

    def committed(self,token):
        with self._lock:
            if not self._pending or token!=next(iter(self._pending)) or self._pending[token][1]!='ready_to_write':
                raise AnalysisEligibilityError('screening retirement is out of order')
            count,_=self._pending.pop(token)
            occurrences,_=self._window_sizes.pop(token)
            self.metrics['committed_windows']+=1;self.metrics['committed_inputs']+=count
            self.metrics['committed_occurrences']+=occurrences
        self.check()

    def snapshot(self):
        with self._lock:
            return {**self.metrics,'windows':tuple((key,*value) for key,value in self._pending.items()),
                'outstanding_inputs':sum(value[0] for value in self._pending.values()),
                'outstanding_occurrences':sum(value[0] for value in self._window_sizes.values()),
                'outstanding_expanded_bytes':sum(value[1] for value in self._window_sizes.values())}

    def finish(self,*,expected_inputs,expected_occurrences):
        self.check()
        if self.closed or self.finished or self._pending:
            raise AnalysisEligibilityError('screening completion is premature')
        self._jobs.put_nowait(('finish',expected_inputs,expected_occurrences))
        kind,result=self._get()
        if kind!='finished':raise AnalysisEligibilityError('screening completion reply differs')
        self.thread.join(timeout=self.shutdown_seconds)
        if self.thread.is_alive():raise AnalysisEligibilityError('screening completion did not stop')
        self.check()
        self.finished=True
        return result

    def close(self):
        if self.closed:return
        primary=sys.exc_info()[1]
        previous=self._error
        if primary is not None:self.fail(primary)
        self._cancelled.set()
        self.thread.join(timeout=self.shutdown_seconds)
        self.closed=True
        with self._lock:
            self.metrics['cancelled_windows']+=len(self._pending)
            for _,state in self._pending.values():self.metrics['cancelled_states'][state]+=1
            self.metrics['cleanup_failed']=self.thread.is_alive()
            self._pending.clear()
            self._window_sizes.clear()
        for channel in (self._jobs,self._results):
            while True:
                try:channel.get_nowait()
                except queue.Empty:break
        if self.thread.is_alive() and self._error is None:
            raise AnalysisEligibilityError('screening thread cleanup failed')
        if (primary is None and self._error is not None
                and not (previous is None and getattr(self._error,'category',None)=='interrupted')):
            raise self._error
