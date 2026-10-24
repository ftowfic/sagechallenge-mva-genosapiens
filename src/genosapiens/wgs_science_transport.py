"""Scientific transport acceleration with the original upstream executor type."""
from copyreg import __newobj__
import pickle
from types import MappingProxyType, MethodType

from genosapiens.evidence_model import EvidenceReference
from genosapiens import wgs_context_admission as ordinary
from genosapiens.wgs_context_admission import admission, transport, MAX_JOB_BYTES, MAX_RESULT_BYTES
from genosapiens.wgs_failures import failure_boundary


def _reference_state(value):
    return (__newobj__, (EvidenceReference,),
            [value.evidence_id, value.record_id, value.public_manifest_sha256], None, None)


def _pack(value, limit, *, operation='context_job_encode'):
    with failure_boundary('bounded_products', operation), transport._LimitedBuffer(limit, operation) as stream:
        writer = pickle.Pickler(stream, protocol=5)
        writer.dispatch_table = {MappingProxyType: ordinary._mapping, EvidenceReference: _reference_state}
        writer.dump(value)
        return stream.getvalue()


def scientific_context_executor(*args, **kwargs):
    """Keep all ordinary owner checks and bind acceleration to this owner only."""
    executor = ordinary.ContextAdmissionExecutor(*args, **kwargs)
    try:
        executor._association_jobs = MethodType(_association_jobs, executor)
        return executor
    except BaseException:
        executor.close()
        raise


def _association_jobs(self, rows, key, eligible, hashes, binding_sha, by_id, context_keys=None):
    projection = self._projection(rows,hashes,binding_sha,by_id=by_id,
        include_contexts=eligible,context_keys=context_keys)
    split = self._association_output_bound(projection) > MAX_RESULT_BYTES//2
    payload = None
    if not split:
        try:
            payload = _pack(('associations',((key,projection),)),MAX_JOB_BYTES)
        except transport.GraphSupportError as error:
            if error.category != 'response_capacity': raise
            split = True
    if split and len(rows)>1:
        half=len(rows)//2
        self.pool.metrics['association_partition_splits'] += 1
        yield from self._association_jobs(rows[:half],key+'a',eligible,hashes,binding_sha,by_id)
        yield from self._association_jobs(rows[half:],key+'b',eligible,hashes,binding_sha,by_id)
        return
    if split:
        association, = rows
        graph, reach = projection[3], projection[2].reach_cache
        reached = tuple(sorted(key for key,(context,namespace) in graph.candidates.items()
            if association.state == 'match' and admission._association_reaches_context(
                association,context,namespace,reach)))
        if len(reached)>1:
            half=len(reached)//2
            self.pool.metrics['context_partition_splits'] += 1
            yield from self._association_jobs(rows,key+'a',eligible,hashes,binding_sha,by_id,frozenset(reached[:half]))
            yield from self._association_jobs(rows,key+'b',eligible,hashes,binding_sha,by_id,frozenset(reached[half:]))
            return
        # A conservative bound can exceed the actual output. An indivisible
        # association/context still uses every real input/result/expanded
        # limit. It is never retried after a worker capacity failure.
        projection = self._projection(rows,hashes,binding_sha,by_id=by_id,
            include_contexts=eligible,context_keys=frozenset(reached))
        payload = _pack(('associations',((key,projection),)),MAX_JOB_BYTES)
    yield payload, ('associations',((key,projection),))
