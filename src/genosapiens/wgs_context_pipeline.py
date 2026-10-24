"""Bounded admission windows and shared-worker decoding/discovery construction.

Discovery owns ten persistent workers through complete source closure.
Bounded windows retain at most two windows and one lookahead row; no disk trust
marker exists.
"""
from contextlib import closing
from collections import defaultdict
from dataclasses import dataclass,replace,asdict
import hashlib
import json
import time
from types import SimpleNamespace
import uuid

MAX_WINDOW_BYTES = 32 * 1024**2


def _discovery_source_sha(sources):
    from genosapiens import candidate_admission as admission
    bindings=admission._source_bindings(admission.AdmissionSources(None,sources.encode,
        sources.open_targets,sources.clinvar,sources.reference,sources.repeat_catalog))
    bindings['annotation']=sources.annotation_binding_sha256
    bindings['regulatory_index_manifest_sha256']=sources.regulatory_index_manifest_sha256
    bindings['encode_paths']=None if sources.encode is None else tuple(
        tuple(str(getattr(pin,key)) for key in ('contract_path','root','manifest_path'))
        for pin in sources.encode)
    bindings['paths']=tuple(None if pin is None else tuple(str(getattr(pin,key)) for key in
        ('contract_path','root','manifest_path')) for pin in (sources.open_targets,sources.clinvar,sources.reference))
    return hashlib.sha256(admission._bytes(bindings)).hexdigest()


def _encoded_windows(source, max_rows, max_bytes):
    rows, size = [], 0
    for value in source:
        if (type(value) is not tuple or len(value) != 14
                or any(type(v) is not int or not 0 < v <= 4*1024**2 for v in value[5:7])
                or any(type(v) is not bytes or len(v) > 4*1024**2+65536 for v in value[3:5])
                or type(value[13]) is not bytes or not 0 < len(value[13]) <= 4*1024**2):
            raise RuntimeError('admission encoded window is invalid')
        cost = value[5]+value[6]+len(value[13])
        if cost > max_bytes:
            raise RuntimeError('admission row exceeds preparation window')
        if rows and (len(rows) >= max_rows or size+cost > max_bytes):
            yield tuple(rows)
            rows, size = [], 0
        rows.append(value)
        size += cost
    if rows:
        yield tuple(rows)


def admission_windows(reader_factory, *, max_rows, decode_window, max_bytes=MAX_WINDOW_BYTES,
                      profile_reader=None):
    """One production discovery path; the decoder drains each bounded window."""
    if type(max_rows) is not int or max_rows <= 0 or not 0 < max_bytes <= MAX_WINDOW_BYTES or not callable(decode_window):
        raise ValueError('invalid admission discovery window')
    def decode_rows(source, binding):
        with closing(_encoded_windows(source, max_rows, max_bytes)) as windows:
            for window in windows:
                yield from decode_window(window, binding)
    rows, size = [], 0
    with closing(reader_factory()) as reader:
        selected=(reader.iter_assessed_rows(_decode_rows=decode_rows,purpose='discovery')
                  if profile_reader is None else reader.iter_selected_rows(profile_reader,
                      _decode_rows=decode_rows,purpose='discovery'))
        with closing(selected) as source:
            for row in source:
                cost = row.encoded.input_size+row.encoded.proof_size+len(row.encoded.screening_facts)
                if rows and (len(rows) >= max_rows or size+cost > max_bytes):
                    yield tuple(rows)
                    rows, size = [], 0
                rows.append(row)
                size += cost
            if rows:
                yield tuple(rows)


class AdmissionDiscoveryExecutor:
    """One pool for strict decode, then query-fed association construction."""
    def __init__(self, *, execution_config, progress, total_inputs=None):
        from genosapiens.wgs_context_admission import _Pool, WORKERS
        from genosapiens.bounded_wgs_config import BoundedWgsExecutionConfig
        if (type(execution_config) is not BoundedWgsExecutionConfig
                or execution_config.public_queries.concurrent_batches != WORKERS):
            raise ValueError('invalid discovery decoder configuration')
        self.max_rows = execution_config.public_queries.batch_rows
        self.progress = progress
        from genosapiens.wgs_admission_graph import live_worker_count
        self.pool = _Pool(worker_count=live_worker_count(total_inputs))
        self.pool.progress_callback=progress
        self.phase='worker_startup'
        self.active = self.finished = self.closed = False
        self.metrics = self.pool.metrics
        self.metrics.update(decode_windows=0, max_window_rows=0, max_window_bytes=0,
                            workers_stopped=False,build_windows=0,build_rows=0,
                            build_jobs=0,worker_build_rows=[0]*self.pool.worker_count,
                            max_unretired_build_jobs=0,unretired_build_jobs=0,
                            max_reorder_jobs=0,canonical_encodes=0,logical_bytes=0,
                            physical_bytes=0,physical_parts=0,
                            compressed_bytes=0,worker_build_associations=[0]*self.pool.worker_count,
                            worker_build_nanoseconds=0,encode_hash_nanoseconds=0,compression_nanoseconds=0,
                            decode_window_nanoseconds=0,build_window_nanoseconds=0,max_owned_compact_window_bytes=0,
                            cancelled_build_jobs=0,
                            unretired_input_bytes=0,reserved_build_result_bytes=0,reorder_bytes=0,
                            max_unretired_input_bytes=0,max_reserved_build_result_bytes=0,max_reorder_bytes=0)
        self.generation=uuid.uuid4().hex
        try:
            self.pool.start()
        except BaseException:
            try:self.close()
            except BaseException:pass
            raise

    def decode_window(self, rows, binding):
        from genosapiens.wgs_context_admission import _decode_admission_window
        if self.active or self.finished or self.closed or self.pool.pending:
            raise RuntimeError('discovery decoder is unavailable')
        self.active = True
        self.phase='decoding_inputs';started=time.monotonic_ns()
        try:
            output = _decode_admission_window(self.pool, rows, binding,
                max_rows=self.max_rows, progress=self.progress)
            self.metrics['decode_windows'] += 1
            self.metrics['max_window_rows'] = max(self.metrics['max_window_rows'],len(rows))
            self.metrics['max_window_bytes'] = max(self.metrics['max_window_bytes'],
                sum(row[5]+row[6]+len(row[13]) for row in rows))
            return output
        except BaseException:
            try:self.close()
            except BaseException:pass
            raise
        finally:
            self.metrics['decode_window_nanoseconds']+=time.monotonic_ns()-started
            self.active = False

    def finish(self, expected_rows):
        if (self.closed or self.finished or self.active or self.pool.pending
                or self.metrics['decoded_inputs'] != expected_rows
                or self.metrics['build_rows'] != expected_rows
                or self.metrics['unretired_build_jobs']):
            raise RuntimeError('discovery decoder coverage is incomplete')
        self.pool.finish()
        self.metrics['workers_stopped'] = True
        self.finished = True

    def close(self):
        if not self.closed:
            self.closed = True
            self.pool.close()

    def build_window(self,rows,*,clinical,resources,admission_owner,binding_sha256,
            admission_completion_sha256,reference_role,reference_identity_id,reference_manifest_sha256,
            regulatory=None,profile_rows=None):
        """Yield in source order; resume only after the caller's durable write.

        A received result still owns its slot. No new job uses that worker until
        every earlier result and this result have been durably retired.
        """
        from genosapiens import candidate_admission as admission
        from genosapiens.wgs_context_admission import _pack,MAX_JOB_BYTES,MAX_RESULT_BYTES,MAX_DISCOVERY_INPUT_BYTES
        from genosapiens.wgs_admission_graph import GraphSupportError
        if self.active or self.finished or self.closed or self.pool.pending or not rows:
            raise RuntimeError('discovery build is unavailable')
        admission_owner._assert_stable()
        if len(rows)>self.max_rows or sum(r.encoded.input_size+r.encoded.proof_size+
                len(r.encoded.screening_facts) for r in rows)>MAX_WINDOW_BYTES:
            raise RuntimeError('discovery build window exceeds capacity')
        if resources.sources.encode is not None and (
                resources.sources.regulatory_index_manifest_sha256 is None or regulatory is None):
            raise RuntimeError('discovery regulatory index results are missing')
        if (regulatory is not None and type(regulatory) is not admission.ShardedRegulatoryAssociations
                or (resources.sources.regulatory_index_manifest_sha256 is None) != (regulatory is None)):
            raise RuntimeError('discovery regulatory results are invalid')
        if clinical is None and any(s is not None for s in (resources.sources.open_targets,resources.sources.clinvar)):
            raise RuntimeError('discovery clinical results are missing')
        if clinical is not None and type(clinical) is not admission.ShardedClinicalAssociations:
            raise RuntimeError('discovery clinical results are invalid')
        subjects={r.value.subject_id for r in rows}
        if clinical is not None and any(p.proof.allele_id not in subjects for p in clinical.proofs):
            raise RuntimeError('discovery clinical result is foreign')
        if regulatory is not None and any(p.allele_id not in subjects for p in regulatory.proofs):
            raise RuntimeError('discovery regulatory result is foreign')
        if profile_rows is not None:
            from genosapiens.wgs_profile_projection import ProfileDecisionRow
            if (type(profile_rows) is not tuple or len(profile_rows)!=len(rows)
                    or any(type(profile) is not ProfileDecisionRow or profile.original_ref!=row.ref
                           for profile,row in zip(profile_rows,rows,strict=True))):
                raise RuntimeError('discovery profile target scope is foreign')
        profiles={} if profile_rows is None else {p.original_ref[1]:p for p in profile_rows}
        self.active=True
        self.phase='building_associations';started=time.monotonic_ns()
        self.metrics['max_owned_compact_window_bytes']=max(self.metrics['max_owned_compact_window_bytes'],
            sum(r.compact_bytes for r in rows))
        self.metrics['build_windows']+=1
        window=self.metrics['build_windows']
        source_sha=_discovery_source_sha(resources.sources)
        query_sha=hashlib.sha256(admission._bytes((clinical,regulatory))).hexdigest()
        common=dict(generation=self.generation,window=window,binding_sha256=binding_sha256,
            admission_binding_sha256=admission_owner.binding_sha256,
            admission_completion_sha256=admission_completion_sha256,source_sha256=source_sha,
            query_sha256=query_sha,reference_role=reference_role,reference_identity_id=reference_identity_id,
            reference_manifest_sha256=reference_manifest_sha256,
            required_order=tuple({(r.value.subject_id,rule) for r in rows for rule in r.value.required_rule_ids}))
        def encode(group,key):
            snapshots=[]
            for row in group:
                row.require_input(row.value)
                compact=row.compact if row.compact is not None else admission_owner._compact(row.value.retained_ordinal)
                header=dict(ref=row.ref,input_sha256=row.encoded.input_sha256,
                    proof_sha256=row.encoded.proof_sha256,proof_size=row.encoded.proof_size,
                    input_size=row.encoded.input_size,proof_count=row.encoded.proof_count)
                snapshots.append((header,row.value,row.proofs,compact))
            assigned={r.value.subject_id for r in group}
            subset=None if clinical is None else replace(clinical,
                proofs=tuple(p for p in clinical.proofs if p.proof.allele_id in assigned))
            meta=dict(common,job=key,refs=tuple(r.ref for r in group),
                proof_sha256s=tuple(r.encoded.proof_sha256 for r in group))
            regulatory_subset=None if regulatory is None else replace(regulatory,
                proofs=tuple(p for p in regulatory.proofs if p.allele_id in assigned))
            body=(meta,tuple(snapshots),resources.sources,subset,regulatory_subset,
                  None if profile_rows is None else tuple(profiles[r.value.retained_ordinal] for r in group))
            job=('discovery',((key,body),))
            try:
                payload=_pack(job,MAX_JOB_BYTES)
                # Reservation includes potential association expansion; actual
                # replies independently obey the unchanged 64-MiB bound.
                if len(group)>1 and len(payload)*len(group)>MAX_RESULT_BYTES//2:
                    raise GraphSupportError('discovery reply reservation exceeds capacity')
                yield payload,job,group
            except GraphSupportError:
                if len(group)==1:
                    yield _pack(job,MAX_DISCOVERY_INPUT_BYTES),job,group
                    return
                middle=len(group)//2
                yield from encode(group[:middle],key+'a')
                yield from encode(group[middle:],key+'b')
        def jobs():
            width=max(1,min(8,(len(rows)+self.pool.worker_count-1)//self.pool.worker_count))
            for first in range(0,len(rows),width):
                yield from encode(rows[first:first+width],str(first))
        source=jobs();held={};ready={};next_submit=next_retire=0
        def reservations():
            if len(held)>self.pool.worker_count:raise RuntimeError('discovery retirement credits exceeded')
            current=dict(unretired_build_jobs=len(held),unretired_input_bytes=sum(v[2] for v in held.values()),
                reserved_build_result_bytes=len(held)*MAX_RESULT_BYTES,reorder_bytes=sum(v.reply_bytes for v in ready.values()))
            self.metrics.update(current)
            for key,value in current.items():
                maximum='max_'+key
                self.metrics[maximum]=max(self.metrics.get(maximum,0),value)
        def submit(index,job):
            nonlocal next_submit
            self.pool.submit(index,job[0],metadata=job[1])
            held[next_submit]=(index,job[2],len(job[0]));next_submit+=1
            reservations()
        try:
            for index in range(self.pool.worker_count):
                job=next(source,None)
                if job is None:break
                submit(index,job)
            while held:
                if next_retire not in ready:
                    index,outputs=self.pool.receive()
                    sequence=next((n for n,v in held.items() if v[0]==index),None)
                    if sequence is None or sequence in ready or len(outputs)!=1:
                        raise RuntimeError('discovery result ownership differs')
                    ready[sequence]=outputs[0]
                    reservations()
                    self.metrics['max_reorder_jobs']=max(self.metrics['max_reorder_jobs'],len(ready))
                    self.progress()
                    continue
                built=ready.pop(next_retire)
                reservations()
                index,original,_=held[next_retire]
                if len(built.rows)!=len(original):raise RuntimeError('discovery reply row coverage differs')
                for output,row in zip(built.rows,original,strict=True):
                    output.subject.row.require_record(row,reference_role=reference_role,
                        reference_identity_id=reference_identity_id)
                    compact=row.compact if row.compact is not None else admission_owner._compact(row.value.retained_ordinal)
                    output.subject.row.require_compact(compact)
                yield built
                # The caller commits parts, subjects and targets before asking
                # for another result. This is the only retirement point.
                del held[next_retire];next_retire+=1
                reservations()
                admission_owner._assert_stable()
                self.metrics['build_jobs']+=1;self.metrics['build_rows']+=len(original)
                self.metrics['worker_build_rows'][index]+=len(original)
                self.metrics['worker_build_associations'][index]+=built.metrics['association_count']
                self.metrics['worker_build_nanoseconds']+=built.metrics['build_nanoseconds']
                for name in ('canonical_encodes','logical_bytes','physical_bytes','physical_parts',
                        'compressed_bytes','encode_hash_nanoseconds','compression_nanoseconds'):
                    self.metrics[name]+=built.metrics[name]
                job=next(source,None)
                if job is not None:submit(index,job)
                self.progress()
        except BaseException:
            try:self.close()
            except BaseException:pass
            self.metrics['cancelled_build_jobs']+=len(held)
            held.clear();ready.clear();reservations()
            raise
        finally:
            self.metrics['build_window_nanoseconds']+=time.monotonic_ns()-started
            source.close();self.active=False


@dataclass(frozen=True,slots=True)
class BuiltDiscoveryRow:
    parts: tuple
    subject: object
    association_count: int
    logical_bytes: int


@dataclass(frozen=True,slots=True)
class BuiltDiscovery:
    rows: tuple
    targets: tuple
    metrics: dict
    reply_bytes: int


def _build_discovery_job(value):
    """Worker-only immutable snapshot transfer; no DB/session or rich decode."""
    from genosapiens import candidate_admission as admission
    from genosapiens import wgs_admission_projection as ap
    from genosapiens import wgs_association_projection as projection
    from genosapiens.wgs_association_codec import encode_association_row,iter_association_parts,MAX_LOGICAL_BYTES
    from genosapiens.wgs_shards import encode_row
    meta,snapshots,sources,clinical,regulatory,profiles=value
    if not 0<len(snapshots)<=8 or len(meta['refs'])!=len(snapshots):
        raise RuntimeError('discovery job coverage is invalid')
    if (sources.encode is not None and (sources.regulatory_index_manifest_sha256 is None or regulatory is None)
            or (sources.regulatory_index_manifest_sha256 is None) != (regulatory is None)
            or clinical is None and any(s is not None for s in (sources.open_targets,sources.clinvar))):
        raise RuntimeError('discovery worker query ownership differs')
    if _discovery_source_sha(sources)!=meta['source_sha256']:
        raise RuntimeError('discovery source commitment differs')
    assigned={ref[3] for ref in meta['refs'] if ref[3] is not None}
    if clinical is not None and any(p.proof.allele_id not in assigned for p in clinical.proofs):
        raise RuntimeError('discovery clinical subset is foreign')
    if regulatory is not None and (type(regulatory) is not admission.ShardedRegulatoryAssociations
            or any(p.allele_id not in assigned for p in regulatory.proofs)):
        raise RuntimeError('discovery regulatory subset is foreign')
    resources=admission.prepare_sharded_resources(sources)
    rows=[]
    for (header,decoded,proofs,compact),ref,digest in zip(snapshots,meta['refs'],meta['proof_sha256s'],strict=True):
        if (header['ref']!=ref or header['input_sha256']!=ref[4] or header['input_size']!=ref[5]
                or header['proof_sha256']!=digest or len(proofs)!=header['proof_count']
                or decoded.retained_ordinal!=ref[1] or decoded.identity_id!=ref[2]
                or (None if decoded.allele is None else decoded.allele.allele_id)!=ref[3]
                or compact.decision.admission_binding_sha256!=meta['admission_binding_sha256']
                or compact.decision.input_sha256!=ref[4] or compact.decision.proof_sha256!=digest):
            raise RuntimeError('discovery snapshot commitment differs')
        rows.append(ap.AdmissionProjectionRecord(SimpleNamespace(**header),decoded,proofs,ap._VALIDATED_ROW_SEAL))
    before=time.monotonic_ns()
    _,associations,links,repeats=admission._load_sharded_associations(tuple(r.value for r in rows),
        tuple(p for r in rows for p in r.proofs),sources,reference_role=meta['reference_role'],
        reference_identity_id=meta['reference_identity_id'],reference_manifest_sha256=meta['reference_manifest_sha256'],
        clinical_associations=clinical,regulatory_associations=regulatory,resource_preparation=resources)
    # The unchanged kernel enumerates missing sources from a set. Spawned
    # processes have independent hash seeds. Restore the original coordinator
    # window's exact traversal order for that suffix, leaving every actual
    # source proof and its multiplicity/order untouched.
    prefix=sum(len(r.proofs) for r in rows)+len(links)+len(repeats)
    required_order={key:index for index,key in enumerate(meta['required_order'])}
    associations=(*associations[:prefix],*sorted(associations[prefix:],key=lambda p:(
        required_order[p.subject,p.rule],admission._RULE_ADAPTERS[p.rule].index(p.adapter))))
    built=time.monotonic_ns()
    grouped,link_groups,repeat_groups=defaultdict(list),defaultdict(list),defaultdict(list)
    for p in associations:grouped[p.subject].append(p)
    for p in links:link_groups[p.allele_id].append(p)
    for p in repeats:repeat_groups[p.allele_ref.allele_id].append(p)
    result=[];logical_bytes=physical_bytes=physical_parts=compressed_bytes=0
    encode_hash_nanoseconds=compression_nanoseconds=0
    for row,(_,_,_,compact) in zip(rows,snapshots,strict=True):
        v=row.value
        applicability=decision_sha=scope_summary=None
        if profiles is not None:
            from genosapiens.wgs_profile_projection import profile_association_applicability,ProfileDecisionRow
            from genosapiens.track1_association_scope import (
                scope_associations, AssociationScopeSummary, curated_associations,
            )
            from genosapiens.wgs_ledger import canonical
            profile=profiles[len(result)]
            if type(profile) is not ProfileDecisionRow or profile.original_ref!=row.ref:
                raise RuntimeError('discovery profile row differs')
            if admission.inv.RULE_IDS[4] in v.required_rule_ids:
                grouped[v.subject_id] = [a for a in grouped[v.subject_id] if a.adapter != 'curated']
                grouped[v.subject_id].extend(curated_associations(profile))
            scoped=scope_associations(profile,tuple(grouped[v.subject_id]))
            grouped[v.subject_id]=scoped.associations
            scope_summary=AssociationScopeSummary.from_result(scoped)
            applicability=profile_association_applicability(profile,scoped.associations)
            decision_sha=hashlib.sha256(canonical(asdict(profile))).hexdigest()
        logical=projection.AssociationRow(v.retained_ordinal,v.identity_id,v.subject_id,row.encoded.input_sha256,
            row.encoded.proof_sha256,tuple(grouped[v.subject_id]),
            # Scoped associations own the selected link sources/provenance.
            # Keeping the original rich link tuple would reintroduce every
            # excluded regulatory biosample into this publication domain.
            tuple(link_groups[v.subject_id]) if profiles is None else (),tuple(repeat_groups[v.subject_id]),
            profile_applicability=applicability,profile_decision_sha256=decision_sha,scope_summary=scope_summary)
        logical.require(row)
        encoded_started=time.monotonic_ns()
        owned=encode_association_row(logical)
        encode_hash_nanoseconds+=time.monotonic_ns()-encoded_started
        subject=projection._context_subject(row,logical,compact,binding_sha256=meta['binding_sha256'],
            admission_completion_sha256=meta['admission_completion_sha256'],reference_role=meta['reference_role'],
            reference_identity_id=meta['reference_identity_id'],owned_association=owned)
        subject_bytes=subject.as_bytes()
        if owned.size_bytes+len(subject_bytes)>MAX_LOGICAL_BYTES:
            raise RuntimeError('discovery logical subject exceeds capacity')
        compress_started=time.monotonic_ns()
        parts=[]
        for part in iter_association_parts(owned):
            # Ship each physical byte string once. Sending its parsed payload
            # too would double base64 storage and exceed the 64-MiB reply bound
            # for an otherwise valid 32-MiB logical singleton.
            parts.append(part.content)
            payload=part.row.payload
            compressed_bytes+=len(payload)//4*3-(len(payload)-len(payload.rstrip('=')))
        parts=tuple(parts)
        compression_nanoseconds+=time.monotonic_ns()-compress_started
        result.append((parts,subject,subject_bytes,len(logical.associations),owned.size_bytes))
        logical_bytes+=owned.size_bytes;physical_bytes+=sum(len(p) for p in parts);physical_parts+=len(parts)
    advancing_associations=associations if profiles is None else tuple(
        proof for row,profile in zip(rows,profiles,strict=True)
        for proof,eligible in zip(grouped[row.value.subject_id],
            profile_association_applicability(profile,tuple(grouped[row.value.subject_id])),strict=True) if eligible)
    targets=tuple(sorted({admission._bytes(t) for t in admission.association_targets(advancing_associations)}))
    metrics=dict(canonical_encodes=len(rows),logical_bytes=logical_bytes,physical_bytes=physical_bytes,
        physical_parts=physical_parts,association_count=sum(item[3] for item in result),build_nanoseconds=built-before,
        compressed_bytes=compressed_bytes,encode_hash_nanoseconds=encode_hash_nanoseconds,
        compression_nanoseconds=compression_nanoseconds)
    return meta,tuple(result),targets,metrics


def _checked_discovery_reply(value,expected,*,reply_bytes):
    """Authenticate live worker reply and retain its exact bounded writer bytes."""
    from genosapiens import candidate_admission as admission
    from genosapiens import wgs_association_projection as projection
    from genosapiens.wgs_association_codec import AssociationPart,MAX_PART_BYTES,MAX_PARTS,MAX_LOGICAL_BYTES
    from genosapiens.wgs_shards import CanonicalRow,_CANONICAL_ROW_SEAL
    meta,snapshots,_,_,_,_=expected
    if type(value) is not tuple or len(value)!=4 or value[0]!=meta:
        raise RuntimeError('discovery reply commitment differs')
    _,rows,targets,metrics=value
    if type(rows) is not tuple or len(rows)!=len(snapshots) or type(targets) is not tuple or targets!=tuple(sorted(set(targets))):
        raise RuntimeError('discovery reply coverage differs')
    result=[]
    for output,ref,proof in zip(rows,meta['refs'],meta['proof_sha256s'],strict=True):
        parts,subject,content,count,size=output
        if (type(subject) is not projection.ContextSubjectProjection or subject.admission_ref!=ref
                or subject.proof_sha256!=proof or subject.discovery_binding_sha256!=meta['binding_sha256']
                or subject.admission_completion_sha256!=meta['admission_completion_sha256']
                or type(content) is not bytes or content!=subject.as_bytes()
                or type(size) is not int or not 0<size or size+len(content)>MAX_LOGICAL_BYTES
                or type(count) is not int or count<0 or type(parts) is not tuple or not 0<len(parts)<=MAX_PARTS):
            raise RuntimeError('discovery reply subject differs')
        canonical=[]
        for index,part_bytes in enumerate(parts):
            if type(part_bytes) is not bytes or not 0<len(part_bytes)<=MAX_PART_BYTES:
                raise RuntimeError('discovery reply part exceeds capacity')
            part=admission._typed(json.loads(part_bytes),AssociationPart,decode=True)
            if (part.retained_ordinal!=ref[1]
                    or part.logical_sha256!=subject.association_sha256 or part.logical_bytes!=size
                    or part.part_index!=index or part.part_count!=len(parts)
                    or type(part_bytes) is not bytes or not 0<len(part_bytes)<=MAX_PART_BYTES):
                raise RuntimeError('discovery reply part differs')
            part.__post_init__()
            if admission._bytes(part)!=part_bytes:
                raise RuntimeError('discovery part bytes differ from header')
            canonical.append(CanonicalRow(part,part_bytes,_CANONICAL_ROW_SEAL))
        result.append(BuiltDiscoveryRow(tuple(canonical),CanonicalRow(subject,content,_CANONICAL_ROW_SEAL),count,size))
    for content in targets:
        target=admission._typed(json.loads(content),admission.AssociationTarget,decode=True)
        if admission._bytes(target)!=content:raise RuntimeError('discovery target is noncanonical')
    expected_counts=dict(canonical_encodes=len(rows),logical_bytes=sum(r.logical_bytes for r in result),
        physical_bytes=sum(len(p.content) for r in result for p in r.parts),
        physical_parts=sum(len(r.parts) for r in result),association_count=sum(r.association_count for r in result),
        compressed_bytes=sum(len(p.row.payload)//4*3-(len(p.row.payload)-len(p.row.payload.rstrip('='))) for r in result for p in r.parts))
    if (type(metrics) is not dict or any(metrics.get(k)!=v for k,v in expected_counts.items())
            or any(type(metrics.get(k)) is not int or metrics[k]<0 for k in
                ('build_nanoseconds','encode_hash_nanoseconds','compression_nanoseconds'))):
        raise RuntimeError('discovery reply counters differ')
    return BuiltDiscovery(tuple(result),targets,metrics,reply_bytes)
