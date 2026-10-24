"""Once-only phenotype seeds and bounded process evaluation of context additions.

Private worker messages are not a file/checkpoint format. Workers see only an
exact association-specific projection, never the prepared reference graph.
Canonical output rows belong to the coordinator until global closure succeeds.
"""
from collections import OrderedDict
from collections.abc import Mapping
from contextlib import closing
from dataclasses import dataclass
from datetime import datetime, timezone
import hashlib
import json
import os
from pathlib import Path
import pickle
import sys
import time
from types import MappingProxyType, SimpleNamespace

from genosapiens import candidate_admission as admission
from genosapiens import wgs_admission_graph as transport
from genosapiens.evidence_model import _validation_session
from genosapiens.wgs_runtime import analysis_progress, analysis_worker_metrics
from genosapiens.wgs_context_payloads import Encoder, header, unpack, ContextPayloadError
from genosapiens.wgs_failures import WgsFailure, static_failure, worker_failure_message, run_cleanup, failure_boundary

WORKERS = 10
ASSOCIATIONS_PER_JOB = 64
BASELINES_PER_JOB = 16
MAX_OUTPUT_COST_BYTES = 1024**2
MAX_JOB_BYTES = transport.MAX_JOB_BYTES
MAX_RESULT_BYTES = transport.MAX_RESULT_BYTES
MAX_DISCOVERY_INPUT_BYTES = 64 * 1024**2
DISCOVERY_FRAGMENT_BYTES = MAX_JOB_BYTES - 8192
_SEAL = object()


def _discovery_fragments(payload,deadline):
    """Only a dense singleton may use these private numbered packets."""
    count=(len(payload)+DISCOVERY_FRAGMENT_BYTES-1)//DISCOVERY_FRAGMENT_BYTES
    digest=hashlib.sha256(payload).hexdigest()
    for index in range(count):
        yield _pack(('discovery_fragment',(digest,index,count,len(payload),deadline,
            payload[index*DISCOVERY_FRAGMENT_BYTES:(index+1)*DISCOVERY_FRAGMENT_BYTES])),MAX_JOB_BYTES)


def _assemble_discovery(payload,jobs,job_id=None):
    import queue
    assembled=bytearray();expected=None;index=0
    while True:
        if type(payload) is not bytes or len(payload)>MAX_JOB_BYTES:
            raise ContextAdmissionError(category='response_capacity', operation='context_worker_execute')
        message=pickle.loads(payload)
        if type(message) is not tuple or len(message)!=2 or message[0]!='discovery_fragment':
            raise ContextAdmissionError(category='malformed_evidence', operation='context_worker_execute')
        digest,number,count,total,deadline,fragment=message[1]
        if (type(digest) is not str or len(digest)!=64 or type(number) is not int or number!=index
                or type(count) is not int or type(total) is not int or not MAX_JOB_BYTES<total<=MAX_DISCOVERY_INPUT_BYTES
                or count!=(total+DISCOVERY_FRAGMENT_BYTES-1)//DISCOVERY_FRAGMENT_BYTES
                or type(deadline) not in (int,float) or not time.monotonic()<deadline<=time.monotonic()+7200
                or type(fragment) is not bytes or len(fragment)!=min(DISCOVERY_FRAGMENT_BYTES,total-len(assembled))):
            raise ContextAdmissionError(category='malformed_evidence', operation='context_worker_execute')
        identity=(digest,count,total,deadline)
        if expected is None:expected=identity
        if identity!=expected:raise ContextAdmissionError(category='authentication', operation='context_worker_execute')
        assembled.extend(fragment);index+=1
        if index==count:break
        try:
            payload=jobs.get(timeout=max(.001,deadline-time.monotonic()))
            if job_id is not None:
                next_id,payload=transport._job_packet(payload)
                if next_id != job_id:
                    raise ContextAdmissionError(category='authentication', operation='context_worker_execute')
        except queue.Empty:raise ContextAdmissionError(category='timeout', operation='context_worker_execute') from None
    result=bytes(assembled)
    if hashlib.sha256(result).hexdigest()!=digest:
        raise ContextAdmissionError(category='malformed_evidence', operation='context_worker_execute')
    kind,rows=pickle.loads(result)
    if kind!='discovery' or len(rows)!=1 or len(rows[0][1][1])!=1:
        raise ContextAdmissionError(category='malformed_evidence', operation='context_worker_execute')
    return result


class ContextAdmissionError(WgsFailure):
    """Static failure; arbitrary constructor text never becomes diagnostics."""
    def __init__(self, message=None, *, category="unexpected", operation="context_result_validate",
                 exception_kind="project_failure"):
        super().__init__("bounded_products", category, operation=operation, exception_kind=exception_kind)


def _mapping(value):
    # Existing source-field owners are immutable in the parent. A worker gets
    # its own bounded projection and validates constructed outputs normally.
    return dict, (dict(value),)


def _pack(value, limit, *, operation='context_job_encode'):
    with failure_boundary('bounded_products', operation), transport._LimitedBuffer(limit, operation) as stream:
        writer = pickle.Pickler(stream, protocol=5)
        writer.dispatch_table = {MappingProxyType: _mapping}
        writer.dump(value)
        return stream.getvalue()


def _seed_row(seed, encoder=None):
    from genosapiens.wgs_contexts import _grouping_bytes
    return (seed.context.context_id, _grouping_bytes(seed), (admission._bytes(seed) if encoder is None else encoder.record(seed)),
            tuple((r.allele_id, r.identity_id) for r in seed.allele_refs))


def _execute(kind, value):
    if kind == 'discovery':
        from genosapiens.wgs_context_pipeline import _build_discovery_job
        return _build_discovery_job(value)
    if kind == 'decode':
        from genosapiens.wgs_admission_projection import _reader_row
        binding, values = value
        row = _reader_row(values, hashlib.sha256(binding).hexdigest())
        decoded, proofs = row.decode(binding=binding)
        # No process-local seal crosses IPC. The receiver creates its owner
        # only after checking this live reply against the exact submitted row.
        return row, decoded, proofs
    with _validation_session():
        if kind == 'baseline':
            context, namespace, state, proof_ids, ledger = value
            seeds = (admission.CandidateContextSeed(context, namespace,
                ('phenotype_first',), (), state, proof_ids, (), ledger),)
            decisions = sources = coverage = ()
        elif kind == 'associations':
            inventory, associations, prepared, graph, binding_sha, hashes = value
            with failure_boundary('bounded_products', 'context_graph_evaluate'):
                seeds, decisions, _, _ = admission._evaluate_graph(
                    inventory, associations, prepared, graph, None,
                    include_phenotype=False, subject_hashes=hashes, binding_sha256=binding_sha)
            with failure_boundary('bounded_products', 'context_evidence_attach'):
                seeds, decisions, sources = admission._attach_public_association_refs(
                    seeds, decisions, associations)
            coverage = tuple(sorted({(d.subject_id, d.rule_id, d.adapter) for d in decisions}))
            if coverage != tuple(sorted({(p.subject, p.rule, p.adapter) for p in associations})):
                raise ContextAdmissionError(category='malformed_evidence', operation='context_result_validate')
        else:
            raise ContextAdmissionError(category='malformed_evidence', operation='context_result_validate')
        with failure_boundary('bounded_products', 'context_payload_encode'):
            encoder = Encoder()
            output = (tuple(_seed_row(s, encoder) for s in seeds),
                      tuple((d.decision_id, encoder.record(d)) for d in decisions),
                      tuple((s.evidence_id, encoder.record(s)) for s in sources), coverage)
            return (*output, encoder.rows())


def _worker(index, jobs, results):
    job = digest = None
    operation = 'context_worker_execute'
    try:
        results.put(('ready', index, os.getpid()))
        while True:
            idle_started=time.monotonic_ns()
            packet = jobs.get()
            idle=time.monotonic_ns()-idle_started
            if packet is None:
                results.put(('stopped', index))
                return
            job, payload = transport._job_packet(packet)
            digest = hashlib.sha256(payload).hexdigest()
            operation = 'context_worker_execute'
            if type(payload) is not bytes or len(payload) > MAX_JOB_BYTES:
                raise ContextAdmissionError(category='response_capacity', operation='context_worker_execute')
            kind, rows = pickle.loads(payload)
            if kind == 'discovery_fragment':
                # The submitted job commits the assembled payload, including
                # when transport fails before all fragments have arrived.
                if type(rows) is tuple and rows and type(rows[0]) is str:
                    digest = rows[0]
                payload=_assemble_discovery(payload,jobs,job)
                digest=hashlib.sha256(payload).hexdigest()
                kind,rows=pickle.loads(payload)
            if type(rows) is not tuple or not 0 < len(rows) <= (8 if kind == 'decode' else BASELINES_PER_JOB):
                raise ContextAdmissionError(category='malformed_evidence', operation='context_worker_execute')
            if kind == 'decode' and sum(value[1][5]+value[1][6] for _, value in rows) > MAX_JOB_BYTES:
                raise ContextAdmissionError(category='response_capacity', operation='context_worker_execute')
            started = time.process_time_ns()
            wall_started=time.monotonic_ns()
            with admission._record_byte_scope():
                output = tuple((key, _execute(kind, value)) for key, value in rows)
            operation = 'context_result_encode'
            content = _pack(output, MAX_RESULT_BYTES, operation=operation)
            results.put(('result', index, hashlib.sha256(payload).hexdigest(),
                         content, time.process_time_ns()-started,time.monotonic_ns()-wall_started,idle,job))
            del payload, rows, output, content
    except BaseException as error:
        results.put(worker_failure_message(error, worker=index, job=job, digest=digest, operation=operation))


@dataclass(frozen=True, slots=True)
class EncodedContextPart:
    seeds: tuple
    decisions: tuple
    public_sources: tuple
    coverage: tuple
    fragments: tuple
    _seal: object

    def __post_init__(self):
        if self._seal is not _SEAL:
            raise ContextAdmissionError(category='authentication', operation='context_result_validate')


@failure_boundary('bounded_products', 'context_result_validate')
def _part(value, *, allowed_keys=None, binding_sha=None):
    try:
        return _checked_part(value, allowed_keys=allowed_keys, binding_sha=binding_sha)
    except (ValueError, TypeError, KeyError, IndexError, AttributeError, UnicodeError, ContextPayloadError):
        raise ContextAdmissionError(category='malformed_evidence', operation='context_result_validate') from None


def _decoded_part(value, expected):
    from genosapiens import wgs_admission_projection as projection
    try:
        binding, columns = expected
        if type(value) is not tuple or len(value) != 3:
            raise ValueError
        row, decoded, proofs = value
        if (type(row) is not projection._ReaderRow or type(decoded) is not projection.LedgerAdmissionInput
                or type(proofs) is not tuple or any(type(p) is not projection.AnnotationAssociationProof for p in proofs)):
            raise ValueError
        actual = (row.retained_ordinal, row.identity_id, row.allele_id, row.input_blob, row.proof_blob,
            row.input_size, row.proof_size, row.input_sha256, row.proof_sha256, row.proof_count,
            row.source_occurrences, row.sequence_requests, projection.canonical(row.gene_ids), row.screening_facts)
        if (actual != columns or any(type(a) is not type(b) for a,b in zip(actual, columns, strict=True))
                or type(row.gene_ids) is not tuple or row.binding_sha256 != hashlib.sha256(binding).hexdigest()
                or type(row.input_bytes) is not bytes or type(row.proof_bytes) is not bytes
                or len(row.input_bytes) != row.input_size or len(row.proof_bytes) != row.proof_size
                or hashlib.sha256(row.input_bytes).hexdigest() != row.input_sha256
                or hashlib.sha256(row.proof_bytes).hexdigest() != row.proof_sha256
                or type(decoded.retained_ordinal) is not int
                or decoded.retained_ordinal != row.retained_ordinal or decoded.identity_id != row.identity_id
                or (None if decoded.allele is None else decoded.allele.allele_id) != row.allele_id
                or len(proofs) != row.proof_count):
            raise ValueError
        return projection.AdmissionProjectionRecord(row, decoded, proofs, projection._VALIDATED_ROW_SEAL)
    except (ValueError, TypeError, AttributeError, KeyError):
        raise ContextAdmissionError(category='malformed_evidence', operation='context_result_validate') from None


def _checked_part(value, *, allowed_keys, binding_sha):
    if type(value) is not tuple or len(value) != 5 or any(type(v) is not tuple for v in value):
        raise ContextAdmissionError(category='malformed_evidence', operation='context_result_validate')
    seeds, decisions, sources, coverage, fragments = value
    shared = {}
    expanded_bytes = sum(len(row[2]) for row in seeds if type(row) is tuple and len(row) == 4 and type(row[2]) is bytes)
    expanded_bytes += sum(len(row[1]) for rows in (decisions, sources) for row in rows
                          if type(row) is tuple and len(row) == 2 and type(row[1]) is bytes)
    for row in fragments:
        if type(row) is not tuple or len(row) != 3 or type(row[0]) is not str or row[0] in shared:
            raise ContextAdmissionError(category='malformed_evidence', operation='context_result_validate')
        digest, size, blob = row
        if type(size) is not int or size <= 0:
            raise ContextAdmissionError(category='response_capacity', operation='context_result_validate')
        expanded_bytes += size
        if expanded_bytes > MAX_RESULT_BYTES:
            raise ContextAdmissionError(category='response_capacity', operation='context_result_validate')
        content = unpack(blob, size, digest)
        shared[digest] = (size, content)
    used = set()
    def metadata(content):
        record = header(content)
        refs = record['shared']
        if any(type(d) is not str or d not in shared for d in refs.values()):
            raise ContextAdmissionError(category='malformed_evidence', operation='context_result_validate')
        if len(content) + sum(shared[d][0] for d in refs.values()) > MAX_RESULT_BYTES:
            raise ContextAdmissionError(category='response_capacity', operation='context_result_validate')
        used.update(refs.values())
        data = record['data']
        if 'context' in refs:
            digest = refs['context']; _, content = shared[digest]
            data['context'] = json.loads(content)
        return data
    if any(type(r) is not tuple or len(r) != 3 or any(type(s) is not str for s in r) for r in coverage):
        raise ContextAdmissionError(category='malformed_evidence', operation='context_result_validate')
    covered = set()
    subjects = {row[0] for row in coverage}
    for rows, width in ((seeds, 4), (decisions, 2), (sources, 2)):
        if any(type(r) is not tuple or len(r) != width or type(r[0]) is not str for r in rows):
            raise ContextAdmissionError(category='malformed_evidence', operation='context_result_validate')
    for key, grouping, content, refs in seeds:
        if (type(grouping) is not bytes or type(content) is not bytes
                or not 0 < len(content) <= MAX_RESULT_BYTES or not 0 < len(grouping) <= MAX_RESULT_BYTES
                or type(refs) is not tuple or any(type(r) is not tuple or len(r) != 2
                    or any(type(s) is not str for s in r) for r in refs)):
            raise ContextAdmissionError(category='malformed_evidence', operation='context_result_validate')
        data, group = metadata(content), json.loads(grouping)
        if (data['context']['context_id'] != key or group['context_id'] != key
                or refs != tuple((r['allele_id'], r['identity_id']) for r in data['allele_refs'])
                or allowed_keys is not None and key not in allowed_keys
                or any(ref[0] not in subjects for ref in refs)):
            raise ContextAdmissionError(category='malformed_evidence', operation='context_result_validate')
    for rows, field in ((decisions, 'decision_id'), (sources, 'evidence_id')):
        for key, content in rows:
            if type(content) is not bytes or not 0 < len(content) <= MAX_RESULT_BYTES:
                raise ContextAdmissionError(category='malformed_evidence', operation='context_result_validate')
            data = metadata(content)
            if data[field] != key:
                raise ContextAdmissionError(category='malformed_evidence', operation='context_result_validate')
            if field == 'decision_id':
                covered.add((data['subject_id'],data['rule_id'],data['adapter']))
                if (binding_sha is not None and data['binding_sha256'] != binding_sha
                        or allowed_keys is not None and data['evaluated_context_id'] is not None
                        and data['evaluated_context_id'] not in allowed_keys
                        or allowed_keys is not None and any(k not in allowed_keys for k in data['context_ids'])):
                    raise ContextAdmissionError(category='authentication', operation='context_result_validate')
    if tuple(sorted(covered)) != coverage:
        raise ContextAdmissionError(category='malformed_evidence', operation='context_result_validate')
    if used != set(shared):
        raise ContextAdmissionError(category='malformed_evidence', operation='context_result_validate')
    return EncodedContextPart(seeds, decisions, sources, coverage, fragments, _SEAL)


class _Pool(transport._Pool):
    def __init__(self, **kwargs):
        super().__init__(worker_target=_worker, **kwargs)
        self.worker_operation = 'context_worker_execute'
        self.result_operation = 'context_result_validate'
        self.error_operations = frozenset({'context_worker_execute', 'context_result_encode',
            'context_result_validate', 'context_job_encode', 'graph_worker_execute',
            'context_graph_evaluate', 'context_evidence_attach', 'context_payload_encode',
            'context_fragment_encode', 'context_record_encode'})
        self.metrics['worker_associations'] = [0]*self.worker_count
        self.metrics.update(worker_decoded_rows=[0]*self.worker_count, decoded_inputs=0,
                            decode_jobs=0, decode_cpu_nanoseconds=0,
                            fragmented_build_jobs=0,max_fragment_packet_bytes=0,max_assembled_build_bytes=0)
        self.expected = {}

    def _failure(self, category, *, index=None, operation=None, kind='project_failure'):
        failure = super()._failure(category,index=index,operation=operation,kind=kind)
        expected = self.expected.get(index)
        if expected is not None and expected[0] in {'associations','baseline'}:
            rows = expected[1]
            self.metrics['failure'].update(association_count=sum(row[2] for row in rows),
                candidate_context_count=sum(len(row[3]) for row in rows))
        return failure

    def submit(self, index, payload, *, metadata=None):
        if type(payload) is not bytes or len(payload) > MAX_DISCOVERY_INPUT_BYTES:
            raise ContextAdmissionError(category='response_capacity', operation='context_result_validate')
        if len(payload)>MAX_JOB_BYTES and (metadata is None or metadata[0]!='discovery'):
            raise ContextAdmissionError(category='malformed_evidence', operation='context_result_validate')
        kind, rows = pickle.loads(payload) if metadata is None else metadata
        if kind not in ('baseline', 'associations', 'decode', 'discovery'):
            raise ContextAdmissionError(category='malformed_evidence', operation='context_result_validate')
        if kind == 'decode' and (type(rows) is not tuple or not 0 < len(rows) <= 8):
            raise ContextAdmissionError(category='response_capacity', operation='context_result_validate')
        if kind == 'discovery' and (type(rows) is not tuple or len(rows)!=1):
            raise ContextAdmissionError(category='malformed_evidence', operation='context_result_validate')
        expected = rows if kind in ('decode','discovery') else tuple((key, (() if kind == 'baseline' else
            tuple(sorted({(p.subject, p.rule, p.adapter) for p in value[1]}))),
            0 if kind == 'baseline' else len(value[1]),
            frozenset((key,)) if kind == 'baseline' else frozenset(value[3].candidates),
            None if kind == 'baseline' else value[4]) for key, value in rows)
        if len(payload)>MAX_JOB_BYTES:
            if (kind!='discovery' or len(rows)!=1 or len(rows[0][1][1])!=1
                    or index in self.pending or not 0<=index<self.worker_count):
                raise ContextAdmissionError(category='malformed_evidence', operation='context_result_validate')
            deadline=time.monotonic()+self.job_seconds
            self.pending[index]=(hashlib.sha256(payload).hexdigest(),tuple(row[0] for row in rows),deadline,len(payload))
            self.expected[index]=(kind,expected)
            self.metrics['submitted_jobs']+=1
            self._reservations()
            self.metrics['max_outstanding']=max(self.metrics['max_outstanding'],len(self.pending))
            self.metrics['max_job_bytes']=max(self.metrics['max_job_bytes'],len(payload))
            self.metrics['max_reserved_input_bytes']=max(self.metrics['max_reserved_input_bytes'],
                sum(r[3] for r in self.pending.values()))
            self.metrics['max_reserved_result_bytes']=max(self.metrics['max_reserved_result_bytes'],len(self.pending)*MAX_RESULT_BYTES)
            self.metrics['fragmented_build_jobs']+=1
            self.metrics['max_assembled_build_bytes']=max(self.metrics['max_assembled_build_bytes'],len(payload))
            self.next_job_id += 1
            self.job_ids[index] = self.next_job_id
            for fragment in _discovery_fragments(payload,deadline):
                self.metrics['max_fragment_packet_bytes']=max(self.metrics['max_fragment_packet_bytes'],len(fragment))
                self.jobs[index].put((self.job_ids[index],fragment),timeout=min(1,max(.001,deadline-time.monotonic())))
        else:
            super().submit(index, payload, keys=tuple(row[0] for row in rows))
        self.expected[index] = (kind, expected)

    @failure_boundary('bounded_products', 'context_result_validate',
                      known_errors=(pickle.UnpicklingError, ValueError, TypeError, EOFError))
    def receive(self):
        message = self._receive(min(r[2] for r in self.pending.values()))
        if len(message) != 7 or message[0] != 'result' or message[1] not in self.pending:
            raise ContextAdmissionError(category='malformed_evidence', operation='context_result_validate')
        _, index, digest, content, cpu,wall,idle = message
        expected_digest, keys, _, _ = self.pending[index]
        kind, expected = self.expected[index]
        if (digest != expected_digest or type(content) is not bytes or len(content) > MAX_RESULT_BYTES
                or any(type(value) is not int or value<0 for value in (cpu,wall,idle))):
            raise ContextAdmissionError(category='authentication', operation='context_result_validate')
        rows = pickle.loads(content)
        if (type(rows) is not tuple or len(rows) != len(keys)
                or any(type(r) is not tuple or len(r) != 2 or r[0] != key for r, key in zip(rows, keys))):
            raise ContextAdmissionError(category='malformed_evidence', operation='context_result_validate')
        parts = []
        if kind == 'discovery':
            from genosapiens.wgs_context_pipeline import _checked_discovery_reply
            parts=tuple(_checked_discovery_reply(value,original,reply_bytes=len(content))
                for (_,value),(_,original) in zip(rows,expected,strict=True))
            self.metrics['cpu_nanoseconds']+=cpu
            self.metrics['max_result_bytes']=max(self.metrics['max_result_bytes'],len(content))
            self.metrics['completed_jobs']+=1
            del self.pending[index],self.expected[index]
            self._received(index,cpu,wall,idle)
            return index,parts
        if kind == 'decode':
            parts = tuple(_decoded_part(value, expected_value)
                          for (_, value), (_, expected_value) in zip(rows, expected, strict=True))
            self.metrics['worker_decoded_rows'][index] += len(parts)
            self.metrics['decoded_inputs'] += len(parts)
            self.metrics['decode_jobs'] += 1
            self.metrics['decode_cpu_nanoseconds'] += cpu
            self.metrics['cpu_nanoseconds'] += cpu
            self.metrics['max_result_bytes'] = max(self.metrics['max_result_bytes'], len(content))
            del self.pending[index], self.expected[index]
            self.metrics['completed_jobs']+=1
            self._received(index,cpu,wall,idle)
            return index, parts
        for (key, value), (_, coverage, _, allowed_keys, binding_sha) in zip(rows, expected, strict=True):
            part = _part(value, allowed_keys=allowed_keys, binding_sha=binding_sha)
            if part.coverage != coverage or (kind == 'baseline' and (
                    len(part.seeds) != 1 or part.seeds[0][0] != key or part.seeds[0][3]
                    or part.decisions or part.public_sources)):
                raise ContextAdmissionError(category='malformed_evidence', operation='context_result_validate')
            parts.append(part)
        self.metrics['worker_contexts'][index] += sum(len(p.seeds) for p in parts)
        self.metrics['worker_associations'][index] += sum(row[2] for row in expected)
        self.metrics['cpu_nanoseconds'] += cpu
        self.metrics['completed_jobs'] += 1
        self.metrics['max_result_bytes'] = max(self.metrics['max_result_bytes'], len(content))
        del self.pending[index], self.expected[index]
        self._received(index,cpu,wall,idle)
        return index, tuple(parts)


def _decode_admission_window(pool, rows, binding, *, max_rows, progress):
    """Strictly decode one bounded window with the owning stage pool."""
    if type(rows) is not tuple or not 0 < len(rows) <= max_rows:
        raise ContextAdmissionError(category='malformed_evidence', operation='context_result_validate')
    size = sum(row[5]+row[6] for row in rows)
    if not 0 < size <= 32*1024**2:
        raise ContextAdmissionError(category='response_capacity', operation='context_result_validate')
    target = max(1, min(MAX_JOB_BYTES, size//(2*pool.worker_count)))

    def encode(values):
        job = ('decode', tuple((v[0], (binding, v)) for v in values))
        try:
            yield _pack(job, MAX_JOB_BYTES), job
        except transport.GraphSupportError:
            if len(values) == 1:
                raise
            middle = len(values)//2
            yield from encode(values[:middle])
            yield from encode(values[middle:])

    def jobs():
        pending, size = [], 0
        for row in rows:
            cost = row[5]+row[6]
            if pending and (len(pending) >= 8 or size+cost > target):
                yield from encode(pending)
                pending, size = [], 0
            pending.append(row)
            size += cost
        if pending:
            yield from encode(pending)
    output = []
    with closing(jobs()) as source:
        for index in range(pool.worker_count):
            job = next(source, None)
            if job is None:
                break
            pool.submit(index, job[0], metadata=job[1])
        while pool.pending:
            index, decoded = pool.receive()
            output.extend(decoded)
            job = next(source, None)
            if job is not None:
                pool.submit(index, job[0], metadata=job[1])
            progress()
    output.sort(key=lambda value: value.encoded.retained_ordinal)
    if tuple(value.encoded.retained_ordinal for value in output) != tuple(row[0] for row in rows):
        raise ContextAdmissionError(category='malformed_evidence', operation='context_result_validate')
    return tuple(output)


@dataclass(frozen=True, slots=True, repr=False)
class SparseContextBatch:
    subjects: tuple
    binding: Mapping
    parts: object
    _seal: object
    association_owner: object = None

    def __post_init__(self):
        from genosapiens.wgs_association_projection import ContextSubjectProjection
        if (self._seal is not _SEAL or type(self.subjects) is not tuple or not self.subjects
                or any(type(value) is not ContextSubjectProjection for value in self.subjects)):
            raise ContextAdmissionError(category='authentication', operation='context_result_validate')
        self.require_source()

    def require_source(self):
        """Retain discovery ownership until context/audit retirement.

        Independent strict full-row callers have no discovery batch. Production
        supplies the exact compact owner; neither path reconstructs rich rows.
        """
        if self.association_owner is not None:
            from genosapiens.wgs_association_projection import AssociationBatch
            owner = self.association_owner
            if type(owner) is not AssociationBatch or self.subjects is not owner.subjects:
                raise ContextAdmissionError(category='authentication', operation='context_result_validate')
            owner.require(owner.resources)


class ContextAdmissionExecutor:
    def __init__(self, prepared, *, execution_config, total_inputs, progress_path=None, profile_reader=None, progress=None):
        from genosapiens.bounded_wgs_config import BoundedWgsExecutionConfig
        if (type(prepared) is not admission._PreparedShardedAdmission
                or prepared.seal is not admission._PREPARED_SHARDED_SEAL
                or type(execution_config) is not BoundedWgsExecutionConfig
                or execution_config.public_queries.concurrent_batches != WORKERS):
            raise ContextAdmissionError(category='authentication', operation='context_result_validate')
        self.progress = progress
        self.progress_subphase = 'association_batches'
        self.prepared, self.config = prepared, execution_config
        self.projections = OrderedDict()
        self.projection_bytes = 0
        self.output_costs = OrderedDict()
        self.output_cost_bytes = 0
        self.total_inputs = total_inputs
        if profile_reader is not None:
            from genosapiens.wgs_profile_projection import ProfileProjectionReader
            if type(profile_reader) is not ProfileProjectionReader:
                raise ContextAdmissionError(category='authentication', operation='context_result_validate')
            profile_reader.require_stable()
        self.profile_reader=profile_reader
        self.expected_assessed_inputs=None
        self.phenotype_keys = frozenset(key for p in prepared.graph.pheno for key in p.context_ids)
        self.baseline_keys = frozenset()
        self.context_keys = set()
        self.pool = _Pool(worker_count=transport.live_worker_count(max(total_inputs, len(self.phenotype_keys))))
        self.pool.metrics.update(association_partition_splits=0,context_partition_splits=0,baseline_partition_splits=0)
        self.path = None if progress_path is None else Path(progress_path)
        self.started = time.monotonic()
        self.last = 0
        self.baselines = self.inputs = self.assessed_inputs = self.batches = self.associations = 0
        self.screening_complete=False; self.eligibility_mode=True; self.eligibility_policy_sha256=None
        self.dispositions={'assess':0,'exclude_common':0,'defer_refcall':0}
        self.occurrence_dispositions=dict.fromkeys(self.dispositions,0)
        self.context_writer_observation=None
        self.accounting_binding=None
        self.accounting_origin=None
        self.admission_completion_sha256=None
        self.discovery_binding_sha256=None
        self.assessed_ordinals=hashlib.sha256()
        self.previous_assessed_ordinal=0
        self.admission_owner=None
        self.active = self.finished = self.closed = self.baseline_done = self.batch_open = False
        self.pool.progress_callback = lambda: self.write('waiting')
        self.write('worker_startup', force=True)
        try:
            self.pool.start()
        except BaseException as error:
            primary = error if isinstance(error, (KeyboardInterrupt, GeneratorExit)) else static_failure(
                error, 'bounded_products', operation='context_worker_execute')
            run_cleanup([('worker_shutdown', self.close)], primary=primary)
            raise primary from None

    def write(self, phase, *, force=False):
        if phase == 'baseline_contexts': self.progress_subphase = 'baselines'
        elif phase == 'index_closure': self.progress_subphase = 'index_closure'
        elif phase == 'evidence_batches': self.progress_subphase = 'association_batches'
        writer = None
        if self.context_writer_observation is not None:
            observed = self.context_writer_observation
            writer = {key:observed.get(key) for key in ('commits','checkpoints','supplied_bytes')}
            writer.update(pending_bytes=observed.get('pending_supplied_bytes'),
                commit_nanoseconds=int(observed.get('commit_seconds',0)*1e9),
                checkpoint_nanoseconds=int(observed.get('checkpoint_seconds',0)*1e9))
            writer.update({key:observed.get('sidecar_bytes',{}).get(source)
                for key,source in (('database_bytes','db'),('wal_bytes','wal'),('shm_bytes','shm'))})
        baseline = self.progress_subphase == 'baselines'
        analysis_progress(self.progress, 'context_assembly', self.progress_subphase,
            counts={'assessed_completed':self.assessed_inputs, 'assessed_total':self.expected_assessed_inputs,
                'associations_completed':self.associations, 'batches_committed':self.batches,
                'baselines_completed':self.baselines, 'baselines_total':len(self.baseline_keys)},
            worker=analysis_worker_metrics(self.pool.metrics), writer=writer,
            unit='contexts' if baseline else 'inputs',
            completed=self.baselines if baseline else self.assessed_inputs,
            total=len(self.baseline_keys) if baseline else self.expected_assessed_inputs)
        now = time.monotonic()
        if self.path is None or not force and now-self.last < 1:
            return
        record = dict(schema='genosapiens.context-formation-progress', version=1,
            at_utc=datetime.now(timezone.utc).isoformat(), phase=phase,
            total_inputs=self.total_inputs, completed_inputs=self.inputs,
            completed_inputs_scope='source_accounting', source_accounting_origin=self.accounting_origin,
            source_accounting_completion_sha256=self.admission_completion_sha256,
            completed_batches=self.batches, completed_associations=self.associations,
            total_baseline_contexts=len(self.baseline_keys), completed_baseline_contexts=self.baselines,
            elapsed_seconds=now-self.started, workers=self.pool.worker_count, **self.pool.metrics)
        if self.eligibility_mode:
            record.update(analysis_eligibility_policy_sha256=self.eligibility_policy_sha256,
                completed_assessed_inputs=self.assessed_inputs,
                total_assessed_inputs=self.expected_assessed_inputs,
                upstream_assessed_inputs=self.dispositions['assess'],
                profile_manifest_sha256=None if self.profile_reader is None else self.profile_reader.manifest_sha256,
                completed_excluded_inputs=self.dispositions['exclude_common'],
                completed_deferred_inputs=self.dispositions['defer_refcall'],
                completed_input_occurrences=sum(self.occurrence_dispositions.values()),
                occurrence_dispositions=dict(self.occurrence_dispositions),
                assessed_associations=self.associations)
        if self.context_writer_observation is not None:
            record['context_writer']=self.context_writer_observation
        if self.admission_owner is not None:
            record['admission_reads']=dict(self.admission_owner.metrics)
        temporary = self.path.with_suffix('.tmp')
        fd = os.open(temporary, os.O_WRONLY|os.O_CREAT|os.O_TRUNC|os.O_NOFOLLOW, 0o600)
        with os.fdopen(fd, 'w') as out:
            json.dump(record, out, sort_keys=True)
            out.write('\n')
        os.replace(temporary, self.path)
        self.last = now

    def observe_writer(self, observation):
        self.context_writer_observation=observation
        self.write('context_writing' if not self.finished else 'index_closure')

    def _run(self, jobs):
        if self.active or self.finished or self.closed:
            raise ContextAdmissionError(category='malformed_evidence', operation='context_result_validate')
        self.active = True
        try:
            for index in range(self.pool.worker_count):
                payload = next(jobs, None)
                if payload is None:
                    break
                self.pool.submit(index, payload[0], metadata=payload[1])
            while self.pool.pending:
                index, parts = self.pool.receive()
                for part in parts:
                    self.context_keys.update(row[0] for row in part.seeds)
                yield from parts
                self.write('baseline_contexts' if self.baselines < len(self.baseline_keys) else 'evidence_batches')
                payload = next(jobs, None)
                if payload is not None:
                    self.pool.submit(index, payload[0], metadata=payload[1])
        except BaseException as error:
            primary = error if isinstance(error, (KeyboardInterrupt, GeneratorExit)) else static_failure(error, 'bounded_products', operation='context_worker_execute')
            run_cleanup([('diagnostic_publish', lambda: self.write('failed', force=True)),
                         ('worker_shutdown', self.close)], primary=primary)
            raise primary from None
        finally:
            self.active = False
            run_cleanup([('worker_shutdown', jobs.close)], primary=sys.exc_info()[1])

    def reuse_completed_accounting(self, reader):
        """Reuse the producer's full-source totals; process only actual subjects."""
        from genosapiens.wgs_admission_projection import AdmissionProjectionReader
        if (type(reader) is not AdmissionProjectionReader or self.screening_complete
                or self.inputs or self.batch_open):
            raise ContextAdmissionError(category='authentication', operation='context_result_validate')
        reader._assert_stable()
        manifest = reader.eligibility_manifest
        if manifest.input_rows != self.total_inputs or reader.completion.rows != self.total_inputs:
            raise ContextAdmissionError(category='malformed_evidence', operation='context_result_validate')
        self.admission_owner = reader._owner
        self.inputs = manifest.input_rows
        self.dispositions = dict(manifest.input_dispositions)
        self.occurrence_dispositions = dict(manifest.occurrence_dispositions)
        self.accounting_binding = reader.binding_sha256
        self.accounting_origin = 'reused_producer_completion'
        self._complete_screening(reader, manifest)

    def _complete_screening(self, reader, manifest):
        self.eligibility_policy_sha256=manifest.eligibility_policy_sha256
        self.admission_completion_sha256=reader.completion.sha256
        from genosapiens.wgs_association_projection import _discovery_domain
        self.expected_assessed_inputs,self.expected_assessed_ordinals_sha256,_,_=\
            _discovery_domain(reader._owner,self.profile_reader)
        self.screening_complete=True
        self.write('evidence_batches',force=True)

    def baseline_parts(self):
        if self.baseline_done:
            raise ContextAdmissionError(category='malformed_evidence', operation='context_result_validate')
        self.baseline_keys=self.phenotype_keys & self.context_keys
        self.progress_subphase = 'baselines'
        self.write('baseline_contexts', force=True)
        def encode(rows):
            try:
                if len(rows)>1 and sum(len(admission._bytes(value))+4096 for _,value in rows)>MAX_RESULT_BYTES//2:
                    raise transport.GraphSupportError(category='response_capacity',operation='context_job_encode')
                yield _pack(('baseline', tuple(rows)), MAX_JOB_BYTES), ('baseline', tuple(rows))
            except transport.GraphSupportError as error:
                if error.category != 'response_capacity': raise
                if len(rows) == 1:
                    raise
                self.pool.metrics['baseline_partition_splits'] += 1
                half = len(rows)//2
                yield from encode(rows[:half])
                yield from encode(rows[half:])
        def jobs():
            rows = []
            for key in sorted(self.baseline_keys):
                rows.append((key, self.prepared.graph.seed_contexts[key]))
                if len(rows) == BASELINES_PER_JOB:
                    yield from encode(rows)
                    rows = []
            if rows:
                yield from encode(rows)
        with closing(self._run(jobs())) as parts:
            for part in parts:
                yield part
                self.baselines += len(part.seeds)
        if self.baselines != len(self.baseline_keys):
            raise ContextAdmissionError(category='malformed_evidence', operation='context_result_validate')
        self.baseline_done = True
        self.write('evidence_batches', force=True)

    def _projection(self, rows, hashes, binding_sha, *, by_id, include_contexts=True, context_keys=None):
        graph, prepared = self.prepared.graph, self.prepared.phenotype
        keys, targets, conditions = set(), {}, {}
        needed_reach = {}
        for p in rows if include_contexts else ():
            if p.state != 'match':
                continue
            condition = admission._normalise_curie(p.condition)
            stable = p.target.split('.')[0] if p.namespace == 'Ensembl' and p.target else p.target
            if p.target:
                target = p.namespace, stable
                possible = graph.candidate_by_target.get(target, ())
                if context_keys is not None: possible = tuple(key for key in possible if key in context_keys)
                targets[target] = possible
            else:
                possible = graph.candidate_by_condition.get(condition, ())
                if context_keys is not None: possible = tuple(key for key in possible if key in context_keys)
                conditions[condition] = possible
            keys.update(possible)
            for key in possible:
                native = graph.candidates[key][0].native_disease_id
                reached = prepared.reach_cache[native]
                selected = needed_reach.setdefault(native, {})
                if condition is not None and condition in reached:
                    selected[condition] = reached[condition]
        subjects = {p.subject for p in rows}
        mini_inventory = SimpleNamespace(alleles=tuple(SimpleNamespace(
            allele_id=a.allele_id, native_identity_id=a.native_identity_id)
            for a in (by_id[s] for s in sorted(subjects) if s in by_id)), source_alleles=())
        ordered_keys = tuple(sorted(keys))
        projected = self.projections.get(ordered_keys)
        if projected is None:
            projected = {name: {key: getattr(graph, name)[key] for key in ordered_keys}
                for name in ('candidates','ledgers','supported','phenotype_state','seed_contexts')}
            # These maps reference the immutable prepared owner's existing
            # values/keys; only dictionary and tuple storage is newly retained.
            cost = sys.getsizeof(ordered_keys) + sys.getsizeof(projected)
            cost += sum(sys.getsizeof(value) for value in projected.values()) + 256
            if cost <= MAX_JOB_BYTES:
                while self.projections and self.projection_bytes + cost > MAX_JOB_BYTES:
                    _, (_, old_cost) = self.projections.popitem(last=False)
                    self.projection_bytes -= old_cost
                self.projections[ordered_keys] = (projected, cost)
                self.projection_bytes += cost
        else:
            self.projections.move_to_end(ordered_keys)
            projected = projected[0]
        mini_graph = SimpleNamespace(**projected,
            candidate_by_target=targets, candidate_by_condition=conditions, pheno=())
        return (mini_inventory, tuple(rows), SimpleNamespace(reach_cache=needed_reach, limits=()),
                mini_graph, binding_sha, {s: hashes[s] for s in subjects})

    def _context_output_cost(self, key):
        cached = self.output_costs.get(key)
        if cached is not None:
            self.output_costs.move_to_end(key)
            return cached[:2]
        graph = self.prepared.graph
        context = graph.candidates[key][0]
        # Whole canonical fields are counted independently. Sharing in pickle
        # or the output fragment codec never lowers this conservative bound.
        ledger = len(admission._bytes(graph.ledgers[key]))
        seed = len(admission._bytes(graph.seed_contexts[key])) + 4096
        decision = len(admission._bytes(context.evidence_refs)) + 2*ledger + 4*len(key) + 4096
        cost = 512 + len(key)
        while self.output_costs and self.output_cost_bytes + cost > MAX_OUTPUT_COST_BYTES:
            _, previous = self.output_costs.popitem(last=False)
            self.output_cost_bytes -= previous[2]
        self.output_costs[key] = (seed, decision, cost)
        self.output_cost_bytes += cost
        return seed, decision

    def _association_output_bound(self, projection):
        inventory, rows, prepared, graph, _, _ = projection
        by_id = {a.allele_id:a for a in inventory.alleles}
        total, seen = 4096, set()
        for association in rows:
            size = len(admission._bytes(association))
            public_refs = 0
            for record in association.source_refs:
                if record.public_manifest_sha256 is not None:
                    identifier = len(admission._bytes(record.source_record_id))
                    public_refs += identifier + 512
                    total += len(admission._record_bytes(record)) + identifier + 1024
            stable = (association.target.split('.')[0] if association.namespace == 'Ensembl'
                      and association.target else association.target)
            condition = admission._normalise_curie(association.condition)
            possible = (() if association.state != 'match' else
                graph.candidate_by_target.get((association.namespace,stable),()) if association.target else
                graph.candidate_by_condition.get(condition,()))
            if not possible:
                total += size + public_refs + 4096
            allele = by_id.get(association.subject)
            allele_bytes = 0 if allele is None else len(admission._bytes((allele.allele_id,allele.native_identity_id)))
            for key in possible:
                seed, decision = self._context_output_cost(key)
                if key not in seen:
                    total += seed
                    seen.add(key)
                native = graph.candidates[key][0].native_disease_id
                reach = prepared.reach_cache.get(native,{}).get(condition)
                relation_bytes = 0 if reach is None else len(admission._bytes(reach[1]))
                # One decision and its worst-case attached seed contribution,
                # including member and decision references. No biological
                # outcome is predicted or omitted by this transport bound.
                total += 2*(size + public_refs + decision + relation_bytes + allele_bytes + 128)
        return total

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


    def batch(self, *, subjects, associations, prepared, binding, association_owner=None):
        if prepared is not self.prepared or self.active or self.closed or self.finished or self.batch_open:
            raise ContextAdmissionError(category='authentication', operation='context_result_validate')
        if self.baseline_done or not self.screening_complete:
            raise ContextAdmissionError(category='malformed_evidence', operation='context_result_validate')
        if self.profile_reader is not None and association_owner is None:
            raise ContextAdmissionError(category='authentication', operation='context_result_validate')
        from genosapiens.wgs_association_projection import ContextSubjectProjection
        if (type(subjects) is not tuple or not 0<len(subjects)<=self.config.public_queries.batch_rows
                or any(type(value) is not ContextSubjectProjection for value in subjects)):
            raise ContextAdmissionError(category='malformed_evidence', operation='context_result_validate')
        if association_owner is not None:
            from genosapiens.wgs_association_projection import AssociationBatch
            if (type(association_owner) is not AssociationBatch
                    or subjects is not association_owner.subjects
                    or association_owner.manifest.profile_manifest_sha256 != (
                        None if self.profile_reader is None else self.profile_reader.manifest_sha256)):
                raise ContextAdmissionError(category='authentication', operation='context_result_validate')
            owned, _, _ = association_owner.require(prepared.resources)
            if len(owned) != len(associations) or any(a is not b for a,b in zip(owned,associations)):
                raise ContextAdmissionError(category='authentication', operation='context_result_validate')
        for value in subjects:
            if (value.admission_completion_sha256!=self.admission_completion_sha256
                    or value.retained_ordinal<=self.previous_assessed_ordinal
                    or value.admission_ref[0]!=self.accounting_binding
                    or self.discovery_binding_sha256 not in (None,value.discovery_binding_sha256)):
                raise ContextAdmissionError(category='malformed_evidence', operation='context_result_validate')
            self.discovery_binding_sha256=value.discovery_binding_sha256
            self.previous_assessed_ordinal=value.retained_ordinal
            self.assessed_ordinals.update((str(value.retained_ordinal)+'\n').encode())
        if admission._hash([value.admission_ref[4] for value in subjects])!=binding.get('input_rows_sha256'):
            raise ContextAdmissionError(category='authentication', operation='context_result_validate')
        expected={(value.subject_id,rule,adapter) for value in subjects for rule,adapter in value.required_coverage}
        if {(p.subject,p.rule,p.adapter) for p in associations}!=expected:
            raise ContextAdmissionError(category='malformed_evidence', operation='context_result_validate')
        hashes={value.subject_id:value.subject_sha256 for value in subjects}
        binding_sha = admission._hash(binding)
        applicability=(True,)*len(associations) if association_owner is None else tuple(
            eligible for row in association_owner.rows
            for eligible in ((True,)*len(row.associations) if row.profile_applicability is None
                             else row.profile_applicability))
        scoped=tuple(sorted(zip(associations,applicability,strict=True),key=lambda pair:(pair[1],
            pair[0].namespace or '',pair[0].target or '',pair[0].condition or '',pair[0].state,
            pair[0].subject,pair[0].rule,pair[0].adapter)))
        self.batch_open = True
        by_id = {value.allele_id:value for value in subjects if value.allele_id is not None}
        def encode(rows, key, eligible):
            yield from self._association_jobs(rows,key,eligible,hashes,binding_sha,by_id)
        def jobs():
            rows=[];eligibility=None;first=0
            for proof,eligible in scoped:
                if rows and (len(rows)==ASSOCIATIONS_PER_JOB or eligibility!=eligible):
                    yield from encode(tuple(rows),str(first),eligibility)
                    first+=len(rows);rows=[]
                rows.append(proof);eligibility=eligible
            if rows:yield from encode(tuple(rows),str(first),eligibility)
        def parts():
            coverage = set()
            with closing(self._run(jobs())) as results:
                for part in results:
                    coverage.update(part.coverage)
                    yield part
            if coverage != expected:
                raise ContextAdmissionError(category='malformed_evidence', operation='context_result_validate')
            self.assessed_inputs += len(subjects)
            self.batches += 1
            self.associations += len(associations)
            self.batch_open = False
            self.write('evidence_batches', force=True)
        return SparseContextBatch(subjects, MappingProxyType(binding), parts(), _SEAL, association_owner)

    def finish(self):
        if (self.active or self.batch_open or self.finished or not self.baseline_done
                or not self.screening_complete or self.inputs != self.total_inputs
                or self.assessed_inputs != self.expected_assessed_inputs or self.baselines != len(self.baseline_keys)):
            raise ContextAdmissionError(category='malformed_evidence', operation='context_result_validate')
        if self.eligibility_mode and sum(self.dispositions.values())!=self.total_inputs:
            raise ContextAdmissionError(category='malformed_evidence', operation='context_result_validate')
        # Full-source callers close membership with the producer's assessed
        # count, strict order and each requested scalar's assess disposition.
        # Profile callers additionally bind the exact selected ordinal digest.
        expected = self.expected_assessed_ordinals_sha256 or self.admission_owner.assessed_ordinals_sha256
        if expected is not None and self.assessed_ordinals.hexdigest()!=expected:
            raise ContextAdmissionError(category='malformed_evidence', operation='context_result_validate')
        if self.profile_reader is not None:self.profile_reader.require_stable()
        self.pool.finish()
        self.finished = True
        self.write('index_closure', force=True)

    def release_preparation(self):
        if not self.finished or self.active or self.batch_open:
            raise ContextAdmissionError(category='malformed_evidence', operation='context_result_validate')
        self.projections.clear()
        self.output_costs.clear()
        self.output_cost_bytes = 0
        self.projection_bytes = 0
        self.prepared = None

    def close(self):
        self.projections.clear()
        self.output_costs.clear()
        self.output_cost_bytes = 0
        self.projection_bytes = 0
        if not self.closed:
            try:
                self.pool.close()
            finally:
                self.closed = True
