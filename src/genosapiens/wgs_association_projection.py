"""Bounded, source-bound association discovery before final graph closure.

Every assessed input has a durable row, including unresolved and unlinked
associations. Public queries run once during discovery. Context formation reads
the same rows against the authenticated admission input/proof commitments.
"""
from collections import defaultdict
from contextlib import closing
from dataclasses import dataclass,asdict,replace
import hashlib
import json
import os
from pathlib import Path
import time

from genosapiens import candidate_admission as admission
from genosapiens.evidence_model import _sha, _typed
from genosapiens.allele_identity import ReferenceRole
from genosapiens.wgs_shards import ScientificShard, _ShardWriter, read_bounded_metadata, WgsScientificError
from genosapiens.wgs_scientific import _bytes, _load_manifest, _read_rows, _write_manifest
from genosapiens.wgs_evidence_store import EvidenceManifest,EvidenceWriter,EvidenceReader,EvidenceStoreError
from genosapiens.track1_association_scope import AssociationScopeSummary, POLICY_SHA256 as ASSOCIATION_SCOPE_SHA256
from genosapiens.wgs_association_codec import (
    CODEC, MAX_LOGICAL_BYTES, MAX_PART_BYTES, MAX_PARTS, AssociationPart, OwnedAssociationRow,
    encode_association_row, iter_association_parts, read_association_rows,
)

_SEAL=object()
MAX_ROW_BYTES=4*1024**2
MAX_TARGET_BYTES=8*1024**2
MAX_DISCOVERY_MANIFEST_BYTES=16*1024**2
MAX_DESCRIPTOR_BYTES=8*1024**2


def _discovery_domain(owner, profile_reader):
    """Bind discovery to a complete profile or the isolated expanded oracle."""
    from genosapiens.wgs_admission_projection import AuthenticatedAdmissionOwner
    if type(owner) is not AuthenticatedAdmissionOwner:
        raise admission.CandidateAdmissionError('authenticated discovery owner required')
    owner._assert_stable()
    if profile_reader is None:
        return (dict(owner.eligibility_manifest.input_dispositions)['assess'],
            owner.assessed_ordinals_sha256,
            dict(owner.eligibility_manifest.occurrence_dispositions)['assess'], None)
    from genosapiens.wgs_profile_projection import ProfileProjectionReader
    if (type(profile_reader) is not ProfileProjectionReader
            or profile_reader.identity['parent']['binding_sha256'] != owner.binding_sha256
            or profile_reader.identity['parent']['completion_sha256'] != owner.completion.sha256
            or profile_reader.identity['parent']['rows'] != owner.row_count
            or profile_reader.identity['parent']['occurrences'] != owner.occurrence_count):
        raise admission.CandidateAdmissionError('foreign discovery profile')
    profile_reader.require_stable()
    counts=profile_reader.manifest['counts']
    return (counts['advancing_rows'],counts['advancing_ordinals_sha256'],
        counts['advancing_occurrences'],profile_reader.manifest_sha256)


@dataclass(frozen=True,slots=True)
class AssociationRow:
    retained_ordinal: int
    identity_id: str
    subject_id: str
    input_sha256: str
    proof_sha256: str
    associations: tuple[admission._Association,...]
    link_proofs: tuple[admission.links.AssociationProof,...]
    repeat_proofs: tuple[admission.repeat.RepeatAssociationProof,...]
    profile_applicability: tuple[bool,...]|None=None
    profile_decision_sha256: str|None=None
    scope_summary: AssociationScopeSummary|None=None

    def __post_init__(self):
        _sha(self.input_sha256);_sha(self.proof_sha256)
        if self.profile_decision_sha256 is not None:_sha(self.profile_decision_sha256)
        if (type(self.retained_ordinal) is not int or self.retained_ordinal<1
                or not self.identity_id or not self.subject_id
                or any(type(values) is not tuple for values in
                       (self.associations,self.link_proofs,self.repeat_proofs))
                or (self.profile_applicability is None)!=(self.profile_decision_sha256 is None)
                or (self.profile_applicability is None)!=(self.scope_summary is None)
                or self.scope_summary is not None and (
                    type(self.scope_summary) is not AssociationScopeSummary
                    or self.scope_summary.grouped_count+self.scope_summary.exclusion_summaries!=len(self.associations)
                    or self.link_proofs)
                or self.profile_applicability is not None and (
                    type(self.profile_applicability) is not tuple
                    or len(self.profile_applicability)!=len(self.associations)
                    or any(type(value) is not bool for value in self.profile_applicability))
                or any(value.subject!=self.subject_id for value in self.associations)
                or any(value.allele_id!=self.subject_id for value in self.link_proofs)
                or any(value.allele_ref.allele_id!=self.subject_id for value in self.repeat_proofs)):
            raise admission.CandidateAdmissionError('invalid discovered association row')

    def require(self,row):
        row.require_input(row.value)
        if ((self.retained_ordinal,self.identity_id,self.subject_id,self.input_sha256,self.proof_sha256)
                != (row.value.retained_ordinal,row.value.identity_id,row.value.subject_id,
                    row.encoded.input_sha256,row.encoded.proof_sha256)
                or {(value.rule,value.adapter) for value in self.associations}
                != {(rule,adapter) for rule in row.value.required_rule_ids
                    for adapter in admission._RULE_ADAPTERS[rule]}):
            raise admission.CandidateAdmissionError('association/admission linkage differs')


@dataclass(frozen=True,slots=True,repr=False)
class ContextSubjectProjection:
    admission_ref: tuple[str,int,str,str|None,str,int]
    proof_sha256: str
    proof_size: int
    subject_id: str
    source_allele_ids: tuple[str,...]
    parent_source_manifest_sha256: str
    parent_working_set_manifest_sha256: str
    native_source_contract_sha256: str
    required_coverage: tuple[tuple[str,str],...]
    eligibility_decision_sha256: str
    native_identity_id: str|None
    subject_sha256: str
    association_sha256: str
    admission_completion_sha256: str
    discovery_binding_sha256: str
    disposition: str='assess'
    schema: str='genosapiens.wgs-context-subject-projection'
    version: int=1

    def __post_init__(self):
        from genosapiens.wgs_admission_projection import AdmissionAccountingRow
        # Reuse the exact original-reference/occurrence scalar contract.
        AdmissionAccountingRow(self.admission_ref,self.source_allele_ids,self.proof_size,0,0,self.disposition)
        for value in (self.proof_sha256,self.parent_source_manifest_sha256,
                self.parent_working_set_manifest_sha256,self.native_source_contract_sha256,
                self.eligibility_decision_sha256,self.subject_sha256,self.association_sha256,
                self.admission_completion_sha256,self.discovery_binding_sha256):_sha(value)
        if (self.schema!='genosapiens.wgs-context-subject-projection' or type(self.version) is not int or self.version!=1
                or self.disposition!='assess' or type(self.subject_id) is not str or not self.subject_id
                or self.native_identity_id is not None and (type(self.native_identity_id) is not str or not self.native_identity_id)
                or (self.allele_id is None)!=(self.native_identity_id is None)
                or self.allele_id is not None and self.subject_id!=self.allele_id
                or self.allele_id is None and self.subject_id not in self.source_allele_ids
                or type(self.required_coverage) is not tuple
                or any(type(pair) is not tuple or len(pair)!=2 or pair[0] not in admission._RULE_ADAPTERS
                    or pair[1] not in admission._RULE_ADAPTERS[pair[0]] for pair in self.required_coverage)
                or self.required_coverage!=tuple((rule,adapter) for rule in self.required_rule_ids
                    for adapter in admission._RULE_ADAPTERS[rule])):
            raise admission.CandidateAdmissionError('invalid context subject projection')

    @property
    def retained_ordinal(self):return self.admission_ref[1]

    @property
    def identity_id(self):return self.admission_ref[2]

    @property
    def allele_id(self):return self.admission_ref[3]

    @property
    def required_rule_ids(self):return tuple(dict.fromkeys(rule for rule,adapter in self.required_coverage))

    def as_bytes(self):return _bytes(self)

    def require_compact(self,compact):
        facts,decision=compact.facts,compact.decision
        if (self.admission_ref!=compact.ref or self.proof_sha256!=decision.proof_sha256
                or self.proof_size!=compact.proof_size or self.subject_id!=facts.subject_id
                or self.source_allele_ids!=facts.source_allele_ids
                or (self.parent_source_manifest_sha256,self.parent_working_set_manifest_sha256,
                    self.native_source_contract_sha256)!=(facts.source_manifest_sha256,
                        facts.working_manifest_sha256,facts.source_routing_sha256)
                or self.required_coverage!=tuple((rule,adapter) for rule in facts.required_rule_ids
                    for adapter in admission._RULE_ADAPTERS[rule])
                or decision.disposition!='assess'
                or self.eligibility_decision_sha256!=hashlib.sha256(decision.as_bytes()).hexdigest()):
            raise admission.CandidateAdmissionError('context subject admission linkage differs')

    def require_association(self,association,*,owned=None):
        if owned is None:
            owned=encode_association_row(association)
        if type(owned) is not OwnedAssociationRow:
            raise admission.CandidateAdmissionError('association byte owner required')
        owned.require(association)
        if ((self.retained_ordinal,self.identity_id,self.subject_id,self.admission_ref[4],self.proof_sha256)
                !=(association.retained_ordinal,association.identity_id,association.subject_id,
                    association.input_sha256,association.proof_sha256)
                or self.association_sha256!=owned.sha256
                or set(self.required_coverage)!={(p.rule,p.adapter) for p in association.associations}):
            raise admission.CandidateAdmissionError('context subject association linkage differs')

    def require_record(self,row,*,reference_role,reference_identity_id):
        row.require_input(row.value)
        value=row.value
        expected_hash=admission._context_subject_hash(allele=value.allele,source=value.unresolved_source,
            reference_role=reference_role,reference_identity_id=reference_identity_id)
        if (self.admission_ref!=row.ref or self.proof_sha256!=row.encoded.proof_sha256
                or self.proof_size!=row.encoded.proof_size or self.subject_id!=value.subject_id
                or self.source_allele_ids!=value.source_allele_ids
                or self.native_identity_id!=(None if value.allele is None else value.allele.native_identity_id)
                or self.required_rule_ids!=value.required_rule_ids or self.subject_sha256!=expected_hash
                or (self.parent_source_manifest_sha256,self.parent_working_set_manifest_sha256,
                    self.native_source_contract_sha256)!=(value.parent_source_manifest_sha256,
                        value.parent_working_set_manifest_sha256,value.native_source_contract_sha256)):
            raise admission.CandidateAdmissionError('context subject original input differs')


def _context_subject(row,association,compact,*,binding_sha256,admission_completion_sha256,
        reference_role,reference_identity_id,owned_association=None):
    owned=encode_association_row(association) if owned_association is None else owned_association
    owned.require(association)
    value=row.value
    result=ContextSubjectProjection(row.ref,row.encoded.proof_sha256,row.encoded.proof_size,
        value.subject_id,value.source_allele_ids,value.parent_source_manifest_sha256,
        value.parent_working_set_manifest_sha256,value.native_source_contract_sha256,
        tuple((rule,adapter) for rule in value.required_rule_ids for adapter in admission._RULE_ADAPTERS[rule]),
        hashlib.sha256(compact.decision.as_bytes()).hexdigest(),
        None if value.allele is None else value.allele.native_identity_id,
        admission._context_subject_hash(allele=value.allele,source=value.unresolved_source,
            reference_role=reference_role,reference_identity_id=reference_identity_id),
        owned.sha256,admission_completion_sha256,binding_sha256)
    result.require_compact(compact);result.require_association(association,owned=owned)
    result.require_record(row,reference_role=reference_role,reference_identity_id=reference_identity_id)
    return result


def load_context_subject_projection(content):
    try:
        if type(content) is not bytes or not 0<len(content)<=MAX_ROW_BYTES:raise ValueError
        from genosapiens.wgs_analysis_eligibility import _unique
        row=_typed(json.loads(content,object_pairs_hook=_unique),ContextSubjectProjection,decode=True)
        if row.as_bytes()!=content:raise ValueError
        return row
    except (ValueError,TypeError,KeyError,RuntimeError,UnicodeError):
        raise admission.CandidateAdmissionError('invalid serialized context subject projection') from None


@dataclass(frozen=True,slots=True)
class AssociationProjection:
    binding_sha256: str
    admission_completion_sha256: str
    clinical_index_manifest_sha256: str
    row_count: int
    association_count: int
    shards: tuple[ScientificShard,...]
    stream_sha256: str
    targets: tuple[admission.AssociationTarget,...]
    subject_count: int
    subject_occurrences: int
    subject_shards: tuple[ScientificShard,...]
    subject_stream_sha256: str
    assessed_ordinals_sha256: str
    eligibility_manifest_sha256: str
    original_inputs: int
    reference_role: ReferenceRole
    reference_identity_id: str
    reference_manifest_sha256: str
    part_count: int
    logical_stream_sha256: str
    profile_manifest_sha256: str|None=None
    regulatory_index_manifest_sha256: str|None=None
    association_scope_sha256: str|None=None
    codec: str=CODEC
    schema: str='genosapiens.wgs-association-projection'
    version: int=5

    def __post_init__(self):
        if self.profile_manifest_sha256 is not None:
            _sha(self.profile_manifest_sha256)
        if self.association_scope_sha256 != (None if self.profile_manifest_sha256 is None else ASSOCIATION_SCOPE_SHA256):
            raise admission.CandidateAdmissionError('association scope policy differs')
        if self.regulatory_index_manifest_sha256 is not None:
            _sha(self.regulatory_index_manifest_sha256)
        target_bytes=tuple(admission._bytes(value) for value in self.targets)
        if (sum(map(len,target_bytes))>MAX_TARGET_BYTES
                or any(sum(len(_bytes(value)) for value in descriptors)>MAX_DESCRIPTOR_BYTES
                    for descriptors in (self.shards,self.subject_shards))):
            raise WgsScientificError('discovery manifest component exceeds metadata capacity',category='response_capacity')
        for value in (self.binding_sha256,self.admission_completion_sha256,
                self.clinical_index_manifest_sha256,self.stream_sha256,self.subject_stream_sha256,
                self.assessed_ordinals_sha256,self.eligibility_manifest_sha256,self.reference_manifest_sha256,
                self.logical_stream_sha256): _sha(value)
        if (self.schema!='genosapiens.wgs-association-projection' or type(self.version) is not int or self.version!=5
                or self.codec!=CODEC or type(self.part_count) is not int or self.part_count<self.row_count
                or type(self.row_count) is not int or self.row_count<0
                or self.part_count>self.row_count*MAX_PARTS
                or type(self.association_count) is not int or self.association_count<0
                or sum(value.row_count for value in self.shards)!=self.part_count
                or any(value.kind!='associations' for value in self.shards)
                or type(self.subject_count) is not int or self.subject_count!=self.row_count
                or type(self.subject_occurrences) is not int or self.subject_occurrences<self.subject_count
                or type(self.original_inputs) is not int or self.original_inputs<self.subject_count
                or sum(value.row_count for value in self.subject_shards)!=self.subject_count
                or any(value.kind!='context_subjects' for value in self.subject_shards)
                or type(self.reference_role) is not ReferenceRole
                or type(self.reference_identity_id) is not str or not self.reference_identity_id
                or target_bytes != tuple(sorted(set(target_bytes)))):
            raise admission.CandidateAdmissionError('invalid association projection coverage')

    def as_bytes(self):
        return _bytes(dict(classification='controlled',llm_policy='no-llm',manifest=self))

    @property
    def sha256(self): return hashlib.sha256(self.as_bytes()).hexdigest()


def load_association_projection(content, *, expected_sha256):
    """Strict complete discovery metadata, with its separately approved ceiling."""
    return _load_manifest(content,expected_sha256=expected_sha256,
        result_type=AssociationProjection,result_key='manifest',max_bytes=MAX_DISCOVERY_MANIFEST_BYTES)


def _completed_association_metadata(*, inputs, admission_owner,
        profile_reader, resources, clinical_index_manifest_sha256, reference_role,
        reference_identity_id, reference_manifest_sha256):
    """Check both manifest closures before any staging payload is copied."""
    from genosapiens.wgs_profile_projection import (
        CompletedProfileDiscoveryInputs, _require, _signature)
    _require(type(inputs) is CompletedProfileDiscoveryInputs and not inputs.profile_only)
    identity = inputs.identity
    descriptor = identity['association_manifest']
    source = Path(descriptor['path']).parent
    _require(_signature(source/'manifest.json')[2] == descriptor['size_bytes'])
    arguments = dict(expected_sha256=descriptor['sha256'],
        binding_sha256=identity['discovery_binding_sha256'],
        admission_completion_sha256=admission_owner.completion.sha256,
        resources=resources, admission_owner=admission_owner, profile_reader=profile_reader)
    reader = AssociationProjectionReader(source, **arguments)
    try:
        value = reader.manifest
        _require(value.clinical_index_manifest_sha256 == clinical_index_manifest_sha256
            and value.reference_role == reference_role and value.reference_identity_id == reference_identity_id
            and value.reference_manifest_sha256 == reference_manifest_sha256)
        paths = tuple(d.relative_path for d in (*value.shards, *value.subject_shards))
        _require(len(set(paths)) == len(paths) and all(Path(p).name == p
            and p not in {'manifest.json', '.', '..'} for p in paths))
        _require(all(_signature(source/d.relative_path)[2] == d.size_bytes
            for d in (*value.shards, *value.subject_shards)))
    finally:
        reader.close()
    return value


def stage_completed_association_projection(destination, *, progress=None, **arguments):
    """Stage completed discovery, excluding all old assessments and mutable state."""
    from genosapiens.wgs_profile_projection import _stage_completed_files
    value = _completed_association_metadata(**arguments)
    descriptor = arguments['inputs'].identity['association_manifest']
    _stage_completed_files(Path(descriptor['path']).parent, destination,
        manifest_sha256=descriptor['sha256'], shards=(*value.shards, *value.subject_shards), progress=progress)
    return value


@dataclass(frozen=True,slots=True)
class AssociationBatch:
    rows: tuple[AssociationRow,...]
    subjects: tuple[ContextSubjectProjection,...]
    resources: admission.PreparedShardedResources
    manifest: AssociationProjection
    admission_owner: object
    encoded: tuple[OwnedAssociationRow,...]
    seal: object

    @property
    def clinical_index_manifest_sha256(self):return self.manifest.clinical_index_manifest_sha256

    def require(self,resources):
        if (self.seal is not _SEAL or self.resources is not resources
                or type(self.manifest) is not AssociationProjection or not self.rows
                or len(self.subjects)!=len(self.rows) or len(self.encoded)!=len(self.rows)
                or self.admission_owner.completion.sha256!=self.manifest.admission_completion_sha256):
            raise admission.CandidateAdmissionError('foreign association batch owner')
        self.admission_owner._assert_stable()
        previous=size=0
        for value,subject,encoded in zip(self.rows,self.subjects,self.encoded,strict=True):
            if (value.profile_applicability is None)!=(self.manifest.profile_manifest_sha256 is None):
                raise admission.CandidateAdmissionError('association target scope is missing or foreign')
            subject.require_association(value,owned=encoded)
            size+=encoded.size_bytes+len(subject.as_bytes())
            if size>MAX_LOGICAL_BYTES:
                raise WgsScientificError('association batch exceeds logical capacity',category='response_capacity')
            subject.require_compact(self.admission_owner._compact(subject.retained_ordinal))
            if (subject.retained_ordinal<=previous
                    or subject.admission_completion_sha256!=self.manifest.admission_completion_sha256
                    or subject.discovery_binding_sha256!=self.manifest.binding_sha256):
                raise admission.CandidateAdmissionError('foreign association subject')
            previous=subject.retained_ordinal
        return tuple(value for row in self.rows for value in row.associations),tuple(
            value for row in self.rows for value in row.link_proofs),tuple(
            value for row in self.rows for value in row.repeat_proofs)

    def require_records(self,projection_rows,resources):
        """Explicit strict full-row reference comparison, outside preparation."""
        result=self.require(resources)
        if type(projection_rows) is not tuple or len(projection_rows)!=len(self.rows):
            raise admission.CandidateAdmissionError('foreign association reference rows')
        for value,subject,row in zip(self.rows,self.subjects,projection_rows,strict=True):
            value.require(row)
            subject.require_record(row,reference_role=self.manifest.reference_role,
                reference_identity_id=self.manifest.reference_identity_id)
        return result


class AssociationProjectionReader:
    def __init__(self,root,*,expected_sha256,binding_sha256,admission_completion_sha256,resources,admission_owner,
            profile_reader=None):
        from genosapiens.wgs_admission_projection import AuthenticatedAdmissionOwner
        self.root=Path(root)
        self.manifest=load_association_projection(read_bounded_metadata(self.root/'manifest.json',
            max_bytes=MAX_DISCOVERY_MANIFEST_BYTES),expected_sha256=expected_sha256)
        rows,ordinals,occurrences,profile_sha256=_discovery_domain(admission_owner,profile_reader)
        if (self.manifest.binding_sha256!=binding_sha256
                or self.manifest.admission_completion_sha256!=admission_completion_sha256
                or type(resources) is not admission.PreparedShardedResources
                or type(admission_owner) is not AuthenticatedAdmissionOwner
                or admission_owner.completion.sha256!=admission_completion_sha256
                or admission_owner.eligibility_manifest.sha256!=self.manifest.eligibility_manifest_sha256
                or admission_owner.row_count!=self.manifest.original_inputs
                or ordinals is not None and ordinals!=self.manifest.assessed_ordinals_sha256
                or rows!=self.manifest.row_count or occurrences!=self.manifest.subject_occurrences
                or profile_sha256!=self.manifest.profile_manifest_sha256
                or resources.sources.regulatory_index_manifest_sha256!=self.manifest.regulatory_index_manifest_sha256):
            raise admission.CandidateAdmissionError('foreign association projection')
        admission_owner._assert_stable()
        self.admission_owner=admission_owner
        self.profile_reader=profile_reader
        self.profile_rows=None if profile_reader is None else (
            row for row in profile_reader.iter_rows() if row.decision['disposition'] in {'selected','unresolved'})
        self.resources=resources
        self.parts=_read_rows(self.root,self.manifest.shards,'associations',AssociationPart,
            expected_stream_sha256=self.manifest.stream_sha256,row_key=lambda value:str(value.retained_ordinal),
            max_row_bytes=MAX_PART_BYTES)
        self.rows=read_association_rows(part for _,_,part in self.parts)
        self.subjects=_read_rows(self.root,self.manifest.subject_shards,'context_subjects',ContextSubjectProjection,
            expected_stream_sha256=self.manifest.subject_stream_sha256,row_key=lambda value:str(value.retained_ordinal))
        self.count=self.associations=self.previous=self.occurrences=0
        self.ordinals=hashlib.sha256()
        self.logical_stream=hashlib.sha256()

    def _next(self):
        try:
            encoded=next(self.rows)
            value=encoded.row
            _,_,subject=next(self.subjects)
        except StopIteration:
            raise admission.CandidateAdmissionError('discovery projection coverage is incomplete') from None
        subject.require_association(value,owned=encoded)
        if self.profile_rows is not None:
            from genosapiens.wgs_profile_projection import profile_association_applicability
            from genosapiens.wgs_ledger import canonical
            profile=next(self.profile_rows,None)
            if (profile is None or profile.original_ref!=subject.admission_ref
                    or value.profile_decision_sha256!=hashlib.sha256(canonical(asdict(profile))).hexdigest()
                    or value.profile_applicability!=profile_association_applicability(profile,value.associations)):
                raise admission.CandidateAdmissionError('association target scope differs from profile')
            from genosapiens.track1_association_scope import validate_scoped_associations
            validate_scoped_associations(profile,value.associations,value.scope_summary)
        elif value.profile_applicability is not None:
            raise admission.CandidateAdmissionError('foreign association target scope')
        if encoded.size_bytes+len(subject.as_bytes())>MAX_LOGICAL_BYTES:
            raise WgsScientificError("discovery subject exceeds logical capacity",category='response_capacity')
        self.logical_stream.update((str(value.retained_ordinal)+":"+encoded.sha256+"\n").encode())
        self.encoded_row=encoded
        compact=self.admission_owner._compact(subject.retained_ordinal)
        subject.require_compact(compact)
        self.admission_row_bytes=subject.admission_ref[5]+subject.proof_size+len(compact.facts.as_bytes())
        if (subject.discovery_binding_sha256!=self.manifest.binding_sha256
                or subject.admission_completion_sha256!=self.manifest.admission_completion_sha256
                or value.retained_ordinal<=self.previous):
            raise admission.CandidateAdmissionError('discovery subject order or parents differ')
        self.previous=value.retained_ordinal
        self.count+=1;self.associations+=len(value.associations);self.occurrences+=len(subject.source_allele_ids)
        self.ordinals.update((str(value.retained_ordinal)+'\n').encode())
        return value,subject

    def batch(self,rows):
        values=[];subjects=[];encoded=[];size=0
        for row in rows:
            value,subject=self._next()
            value.require(row)
            subject.require_record(row,reference_role=self.manifest.reference_role,
                reference_identity_id=self.manifest.reference_identity_id)
            size+=self.encoded_row.size_bytes+len(subject.as_bytes())
            if size>MAX_LOGICAL_BYTES or len(values)>=self.admission_owner.config.public_queries.batch_rows:
                raise WgsScientificError('association batch exceeds logical capacity',category='response_capacity')
            values.append(value);subjects.append(subject);encoded.append(self.encoded_row)
        return AssociationBatch(tuple(values),tuple(subjects),self.resources,self.manifest,self.admission_owner,tuple(encoded),_SEAL)

    def batches(self,*,max_rows,max_bytes=32*1024**2):
        """Stream the sealed projection without another rich admission pass."""
        if (type(max_rows) is not int or not 0<max_rows<=self.admission_owner.config.public_queries.batch_rows
                or type(max_bytes) is not int or not 0<max_bytes<=32*1024**2 or self.count):
            raise admission.CandidateAdmissionError('invalid discovery batch capacity or position')
        values=[];subjects=[];encoded=[];size=0;original_size=0
        while self.count<self.manifest.row_count:
            value,subject=self._next()
            cost=self.encoded_row.size_bytes+len(subject.as_bytes())
            original_cost=self.admission_row_bytes
            if cost>max_bytes or original_cost>max_bytes:
                raise admission.CandidateAdmissionError('discovery subject exceeds context batch capacity')
            if values and (len(values)>=max_rows or size+cost>max_bytes or original_size+original_cost>max_bytes):
                yield AssociationBatch(tuple(values),tuple(subjects),self.resources,self.manifest,self.admission_owner,tuple(encoded),_SEAL)
                values=[];subjects=[];encoded=[];size=0;original_size=0
            values.append(value);subjects.append(subject);encoded.append(self.encoded_row);size+=cost;original_size+=original_cost
        if values:
            yield AssociationBatch(tuple(values),tuple(subjects),self.resources,self.manifest,self.admission_owner,tuple(encoded),_SEAL)
        self.finish()

    def finish(self):
        if (next(self.rows,None) is not None or self.count!=self.manifest.row_count
                or next(self.subjects,None) is not None or self.count!=self.manifest.subject_count
                or self.occurrences!=self.manifest.subject_occurrences
                or self.ordinals.hexdigest()!=self.manifest.assessed_ordinals_sha256
                or self.associations!=self.manifest.association_count
                or self.logical_stream.hexdigest()!=self.manifest.logical_stream_sha256):
            raise admission.CandidateAdmissionError('incomplete association projection consumption')
        self.admission_owner._assert_stable()
        if self.profile_rows is not None and next(self.profile_rows,None) is not None:
            raise admission.CandidateAdmissionError('incomplete association profile coverage')
        _discovery_domain(self.admission_owner,self.profile_reader)

    def close(self):
        if self.profile_rows is not None:self.profile_rows.close()
        self.rows.close();self.parts.close();self.subjects.close()


def association_at(root,manifest,retained_ordinal):
    """Read original-ordinal evidence, including an input with no graph links."""
    if type(manifest) is not AssociationProjection or type(retained_ordinal) is not int or retained_ordinal<1:
        raise admission.CandidateAdmissionError('invalid association evidence address')
    candidates=tuple(shard for shard in manifest.shards
        if int(shard.first_key)<=retained_ordinal<=int(shard.last_key))
    if not candidates:
        raise admission.CandidateAdmissionError('association evidence address is absent')
    result=None
    with closing(_read_rows(Path(root),candidates,'associations',AssociationPart,
            row_key=lambda value:str(value.retained_ordinal),max_row_bytes=MAX_PART_BYTES)) as parts, closing(read_association_rows(
                part for _,_,part in parts if part.retained_ordinal==retained_ordinal)) as rows:
        for owned in rows:
            if result is not None:
                raise admission.CandidateAdmissionError('association evidence address is duplicated')
            result=owned.row
    if result is None:
        raise admission.CandidateAdmissionError('association evidence address is absent')
    return result


@dataclass(frozen=True,slots=True)
class AssessmentProjection:
    association_projection_sha256: str
    admission_completion_sha256: str
    eligibility_manifest_sha256: str
    policy_sha256: str
    original_inputs: int
    assessed_inputs: int
    decision_count: int
    unlinked_decision_count: int
    public_source_count: int
    decisions_sha256: str
    public_sources_sha256: str
    evidence: EvidenceManifest
    schema: str='genosapiens.wgs-assessment-projection'
    version: int=1

    def __post_init__(self):
        for value in (self.association_projection_sha256,self.admission_completion_sha256,
                self.eligibility_manifest_sha256,self.policy_sha256,self.decisions_sha256,self.public_sources_sha256): _sha(value)
        if (self.schema!='genosapiens.wgs-assessment-projection' or type(self.version) is not int or self.version!=1
                or any(type(n) is not int or n<0 for n in (self.original_inputs,self.assessed_inputs,self.decision_count,
                    self.unlinked_decision_count,self.public_source_count))
                or self.assessed_inputs>self.original_inputs
                or self.unlinked_decision_count>self.decision_count):
            raise admission.CandidateAdmissionError('invalid assessment projection coverage')

    def as_bytes(self):
        return _bytes(dict(classification='controlled',llm_policy='no-llm',manifest=self))

    @property
    def sha256(self): return hashlib.sha256(self.as_bytes()).hexdigest()


class AssessmentProjectionOwner:
    def __init__(self,root,associations,*,policy_sha256,profile_reader=None,progress=None,publication_counts=None):
        self.root,self.associations,self.policy_sha256=Path(root),associations,policy_sha256
        if type(associations) is not AssociationProjection:
            raise admission.CandidateAdmissionError('invalid assessment source owner')
        if (None if profile_reader is None else profile_reader.manifest_sha256)!=associations.profile_manifest_sha256:
            raise admission.CandidateAdmissionError('assessment profile owner differs')
        if profile_reader is not None:
            from genosapiens.wgs_profile_projection import ProfileProjectionReader
            if type(profile_reader) is not ProfileProjectionReader:
                raise admission.CandidateAdmissionError('assessment profile owner is invalid')
            profile_reader.require_stable()
        self.profile_reader=profile_reader
        self.manifest=None
        self.progress=progress
        self.publication_counts=publication_counts

    def publish(self,connection,payloads,admission_reader):
        """Retain every source decision before the disposable context index closes."""
        from genosapiens.evidence_model import _typed,PublicSourcedEvidence
        if read_bounded_metadata(self.root/'manifest.json',max_bytes=MAX_DISCOVERY_MANIFEST_BYTES) != self.associations.as_bytes():
            raise admission.CandidateAdmissionError('assessment discovery manifest changed')
        rows,expected_ordinals,_,profile_sha256=_discovery_domain(admission_reader._owner,self.profile_reader)
        if (self.manifest is not None or admission_reader.completion.sha256!=self.associations.admission_completion_sha256
                or rows!=self.associations.row_count
                or expected_ordinals is not None and expected_ordinals!=self.associations.assessed_ordinals_sha256
                or profile_sha256!=self.associations.profile_manifest_sha256
                or admission_reader.eligibility_manifest.sha256!=self.associations.eligibility_manifest_sha256
                or admission_reader.row_count!=self.associations.original_inputs):
            raise admission.CandidateAdmissionError('assessment projection source closure differs')
        ordinals=hashlib.sha256();subjects=0
        for ordinal, in connection.execute('SELECT retained_ordinal FROM assessment_subjects ORDER BY retained_ordinal'):
            ordinals.update((str(ordinal)+'\n').encode());subjects+=1
        if (subjects!=self.associations.subject_count
                or ordinals.hexdigest()!=self.associations.assessed_ordinals_sha256):
            raise admission.CandidateAdmissionError('assessment subject source coverage differs')
        admission_reader._assert_stable()
        directory=self.root/'assessments';directory.mkdir(mode=0o700)
        counts=dict(decisions=0,unlinked=0,sources=0)
        sequence='decisions'
        nodes=0
        def report(node_count=None, *, force=False):
            nonlocal nodes
            if node_count is not None:nodes=node_count
            from genosapiens.wgs_runtime import analysis_progress
            total=None if self.publication_counts is None else self.publication_counts[
                'decisions' if sequence=='decisions' else 'public_sources']
            analysis_progress(self.progress,'context_assembly','assessment_'+sequence,
                counts={'assessment_'+sequence+'_completed':counts[sequence],
                    'assessment_'+sequence+'_total':total,'assessment_nodes':nodes},
                unit='records',completed=counts[sequence],total=total,force=force)
        def decisions():
            for batch,key,digest in connection.execute(
                    'SELECT batch_ordinal,decision_id,payload_sha256 FROM decisions ORDER BY batch_ordinal,decision_id,payload_sha256'):
                value=_typed(json.loads(payloads.expanded(digest)),admission.AdmissionDecision,decode=True)
                ordinal=connection.execute('SELECT retained_ordinal FROM assessment_subjects WHERE subject_id=?',
                    (value.subject_id,)).fetchone()
                if value.decision_id!=key or ordinal is None:
                    raise admission.CandidateAdmissionError('assessment decision address differs')
                yield dict(retained_ordinal=ordinal[0],batch_ordinal=batch,decision=value)
                counts['decisions']+=1;counts['unlinked']+=int(not value.context_ids)
                if counts['decisions']%32==0:report()
        def sources():
            for batch,key,content in connection.execute('SELECT s.batch_ordinal,s.evidence_id,p.canonical_source '
                    'FROM public_sources s JOIN public_payloads p ON p.sha256=s.source_sha256 ORDER BY s.batch_ordinal,s.evidence_id'):
                from genosapiens.wgs_context_payloads import expand
                value=_typed(json.loads(expand(content,payloads.get)),PublicSourcedEvidence,decode=True)
                if value.evidence_id!=key:
                    raise admission.CandidateAdmissionError('assessment public source address differs')
                yield dict(batch_ordinal=batch,source=value)
                counts['sources']+=1
                if counts['sources']%32==0:report()
        with EvidenceWriter(directory,scratch_directory=directory,progress=report) as writer:
            report(force=True)
            decision_sha=writer.add_sequence('decisions','source-assessment-decisions-v1',decisions())
            report(writer.node_count,force=True)
            sequence='sources'
            report(force=True)
            source_sha=writer.add_sequence('public-sources','source-assessment-public-evidence-v1',sources())
            report(writer.node_count,force=True)
            evidence=writer.finish()
        self.manifest=AssessmentProjection(self.associations.sha256,self.associations.admission_completion_sha256,
            admission_reader.eligibility_manifest.sha256,self.policy_sha256,admission_reader.row_count,
            self.associations.row_count,counts['decisions'],counts['unlinked'],counts['sources'],
            decision_sha,source_sha,evidence)
        _write_manifest(directory,'manifest.json',self.manifest.as_bytes())
        return self.manifest


def assessment_parents(associations,assessments):
    if (type(associations) is not AssociationProjection or type(assessments) is not AssessmentProjection
            or assessments.association_projection_sha256!=associations.sha256):
        raise admission.CandidateAdmissionError('foreign assessment roots')
    return dict(association_projection=associations.sha256,assessment_projection=assessments.sha256)


def verify_assessment_projections(scope,parents,*,full=True,scratch_directory=None):
    """Authenticate retained domain evidence using the scientific parent hashes."""
    keys={'association_projection','assessment_projection'}
    if not keys & parents.keys():
        if any(key.startswith('analysis_profile_') for key in parents):
            raise admission.CandidateAdmissionError('profile assessment parents are missing')
        return None
    if not keys<=parents.keys():
        raise admission.CandidateAdmissionError('incomplete assessment parents')
    root=Path(scope)/'associations'
    associations=load_association_projection(read_bounded_metadata(root/'manifest.json',
        max_bytes=MAX_DISCOVERY_MANIFEST_BYTES),expected_sha256=parents['association_projection'])
    assessments=_load_manifest(read_bounded_metadata(root/'assessments/manifest.json'),expected_sha256=parents['assessment_projection'],
        result_type=AssessmentProjection,result_key='manifest')
    if (assessment_parents(associations,assessments)!={key:parents[key] for key in keys}
            or associations.admission_completion_sha256!=assessments.admission_completion_sha256
            or associations.row_count!=assessments.assessed_inputs
            or associations.original_inputs!=assessments.original_inputs
            or associations.eligibility_manifest_sha256!=assessments.eligibility_manifest_sha256
            or parents.get('analysis_eligibility_manifest')!=assessments.eligibility_manifest_sha256
            or parents.get('analysis_profile_manifest')!=associations.profile_manifest_sha256):
        raise admission.CandidateAdmissionError('assessment source parents differ')
    if full:
        # Discovery and context consumption already validate the rich rows.
        # Scientific readback authenticates their exact durable bytes; rich
        # source proofs are decoded only by an explicit evidence lookup.
        for shards,expected_count,expected_stream in (
                (associations.shards,associations.part_count,associations.stream_sha256),
                (associations.subject_shards,associations.subject_count,associations.subject_stream_sha256)):
            stream=hashlib.sha256();count=0
            for shard in shards:
                from genosapiens.wgs_scientific import _scientific_shard_lines
                shard_count=shard_size=0;digest=hashlib.sha256()
                with closing(_scientific_shard_lines(root,shard,max_row_bytes=MAX_PART_BYTES)) as lines:
                    for line in lines:
                        stream.update(line);digest.update(line);shard_count+=1;shard_size+=len(line)
                if (shard_count!=shard.row_count or shard_size!=shard.size_bytes
                        or digest.hexdigest()!=shard.sha256):
                    raise admission.CandidateAdmissionError('assessment source shard changed')
                count+=shard_count
            if count!=expected_count or stream.hexdigest()!=expected_stream:
                raise admission.CandidateAdmissionError('assessment source coverage differs')
        with EvidenceReader(assessments.evidence,root/'assessments',scratch_directory=scratch_directory) as reader:
            for name,kind,digest,count in (('decisions','source-assessment-decisions-v1',assessments.decisions_sha256,assessments.decision_count),
                    ('public-sources','source-assessment-public-evidence-v1',assessments.public_sources_sha256,assessments.public_source_count)):
                if (reader.reference(name,kind=kind).sha256!=digest or reader._lookup(digest)[0]!='array'
                        or reader._lookup(digest)[2]!=count):
                    raise EvidenceStoreError('assessment evidence root differs')
    return associations,assessments


def assessment_retained_files(scope,parents):
    manifests=verify_assessment_projections(scope,parents,full=False)
    if manifests is None: return
    associations,assessments=manifests;root=Path(scope)/'associations'
    for relative,digest in (('manifest.json',associations.sha256),('assessments/manifest.json',assessments.sha256)):
        path=root/relative
        yield path,'associations/'+relative,digest,path.stat().st_size
    for shard in (*associations.shards,*associations.subject_shards):
        yield root/shard.relative_path,'associations/'+shard.relative_path,shard.sha256,shard.size_bytes
    for shard in (*assessments.evidence.node_shards,*assessments.evidence.root_shards):
        relative='assessments/'+shard.relative_path
        yield root/relative,'associations/'+relative,shard.sha256,shard.size_bytes


def iter_assessment_decisions(scope,parents,*,retained_ordinal=None,scratch_directory=None):
    """Explicit typed replay of retained decisions, including missing links."""
    from genosapiens.evidence_model import _typed
    _,manifest=verify_assessment_projections(scope,parents,full=False)
    if retained_ordinal is not None and (type(retained_ordinal) is not int or not 1<=retained_ordinal<=manifest.original_inputs):
        raise admission.CandidateAdmissionError('invalid assessment decision address')
    count=unlinked=0
    with EvidenceReader(manifest.evidence,Path(scope)/'associations/assessments',scratch_directory=scratch_directory) as reader:
        root=reader.reference('decisions',kind='source-assessment-decisions-v1')
        if root.sha256!=manifest.decisions_sha256:
            raise admission.CandidateAdmissionError('assessment decision root changed')
        for row in reader.iter_array(root,max_item_bytes=64*1024**2):
            if (set(row)!={'retained_ordinal','batch_ordinal','decision'}
                    or type(row['retained_ordinal']) is not int or not 1<=row['retained_ordinal']<=manifest.original_inputs
                    or type(row['batch_ordinal']) is not int or row['batch_ordinal']<1):
                raise admission.CandidateAdmissionError('invalid assessment decision address')
            value=_typed(row['decision'],admission.AdmissionDecision,decode=True)
            count+=1;unlinked+=int(not value.context_ids)
            if retained_ordinal is None or row['retained_ordinal']==retained_ordinal:
                yield row['retained_ordinal'],value
        if (count,unlinked)!=(manifest.decision_count,manifest.unlinked_decision_count):
            raise admission.CandidateAdmissionError('assessment decision coverage changed')


def _write_discovery_batch(batch, *, writer, subjects, index, previous,
        logical_stream, ordinals, timing):
    """Durably publish an owned batch before its executor credit is retired."""
    count=association_count=subject_occurrences=part_count=0
    write_started=time.monotonic_ns()
    for row in batch.rows:
        if row.logical_bytes+len(row.subject.content)>MAX_LOGICAL_BYTES:
            raise WgsScientificError('discovery subject exceeds logical capacity',category='response_capacity')
        subject=row.subject.row
        ordinal=subject.retained_ordinal
        if ordinal<=previous:
            raise admission.CandidateAdmissionError('association discovery input order differs')
        previous=ordinal
        for part in row.parts:
            writer.write(part.row,key=str(ordinal),canonical_bytes=part)
            part_count+=1
        logical_stream.update((str(ordinal)+":"+subject.association_sha256+"\n").encode())
        subjects.write(subject,key=str(ordinal),canonical_bytes=row.subject)
        subject_occurrences+=len(subject.source_allele_ids)
        ordinals.update((str(ordinal)+'\n').encode())
        count+=1;association_count+=row.association_count
    index.executemany('INSERT OR IGNORE INTO targets VALUES (?)',((v,) for v in batch.targets))
    writer.sync();subjects.sync()
    commit_started=time.monotonic_ns()
    try:index.commit()
    finally:timing['index_commit_nanoseconds']+=time.monotonic_ns()-commit_started
    timing['index_commits']+=1
    timing['write_nanoseconds']+=time.monotonic_ns()-write_started
    return previous,count,association_count,subject_occurrences,part_count


def _regulatory_batch(values, *, sources, regulatory_index, with_size=False, profile_rows=None):
    """Query one live authenticated index; transfer only complete bounded proofs."""
    from genosapiens.regulatory_index import RegulatoryIndex,indexed_regulatory_report
    from genosapiens.wgs_context_admission import MAX_DISCOVERY_INPUT_BYTES
    if sources.regulatory_index_manifest_sha256 is None:
        if sources.encode is not None or regulatory_index is not None:
            raise admission.CandidateAdmissionError('discovery regulatory index binding is missing')
        return (None,0) if with_size else None
    if (type(regulatory_index) is not RegulatoryIndex
            or regulatory_index.manifest_sha256!=sources.regulatory_index_manifest_sha256
            or sources.encode is None or sources.reference is None
            or sorted((pin.contract_sha256,pin.manifest_sha256) for pin in sources.encode)!=sorted(
                (row['contract_sha256'],row['manifest_sha256']) for row in regulatory_index.manifest['sources'])
            or sources.reference.manifest_sha256!=regulatory_index.manifest['reference']['reference_manifest_sha256']
            or sources.reference.contract_sha256!=regulatory_index.manifest['reference']['reference_contract_sha256']):
        raise admission.CandidateAdmissionError('discovery regulatory index binding differs')
    regulatory_index.require_stable()
    alleles=tuple(value.allele for value in values if value.allele is not None
                  and admission.inv.RULE_IDS[2] in value.required_rule_ids)
    selected_targets=None
    if profile_rows is not None:
        from genosapiens.wgs_profile_projection import ProfileDecisionRow
        if (type(profile_rows) is not tuple or len(profile_rows)!=len(values)
                or any(type(p) is not ProfileDecisionRow or p.original_ref[1]!=v.retained_ordinal
                    or p.original_ref[2]!=v.identity_id for p,v in zip(profile_rows,values,strict=True))):
            raise admission.CandidateAdmissionError('regulatory discovery profile binding differs')
        by_subject={v.subject_id:p for v,p in zip(values,profile_rows,strict=True)}
        selected_targets={a.allele_id:tuple(sorted({r['target_id'] for r in by_subject[a.allele_id].decision['routes']
            if r['route']=='regulatory' and r['outcome']=='selected' and r['target_id'] is not None})) for a in alleles}
    proofs=[];size=0
    for proof in indexed_regulatory_report(alleles,index=regulatory_index,reference=sources.reference,
            selected_targets=selected_targets):
        size+=len(admission._bytes(proof))
        if size>MAX_DISCOVERY_INPUT_BYTES:
            raise WgsScientificError('regulatory window exceeds discovery capacity',category='response_capacity')
        proofs.append(proof)
    regulatory_index.require_stable()
    result=admission.ShardedRegulatoryAssociations(regulatory_index.manifest_sha256,
        tuple(pin.manifest_sha256 for pin in sources.encode),sources.reference.manifest_sha256,tuple(proofs))
    return (result,size) if with_size else result


def _regulatory_windows(values, *, sources, regulatory_index, on_input=None, profile_rows=None):
    """Split on complete input boundaries before aggregate proofs exceed capacity.

    Each input is queried once. A dense individual input still fails at the
    original limit; no proof is dropped, and no query is retried after overflow.
    Only the current bounded window and one bounded lookahead input are held.
    """
    from genosapiens.wgs_context_admission import MAX_DISCOVERY_INPUT_BYTES
    if sources.regulatory_index_manifest_sha256 is None:
        bound=_regulatory_batch(values,sources=sources,regulatory_index=regulatory_index)
        if on_input is not None:on_input(len(values))
        yield 0,len(values),bound
        return
    start=size=0;proofs=[];bound=None
    for position,value in enumerate(values):
        bound,cost=_regulatory_batch((value,),sources=sources,regulatory_index=regulatory_index,with_size=True,
            profile_rows=None if profile_rows is None else profile_rows[position:position+1])
        if on_input is not None:on_input(1)
        if position>start and size+cost>MAX_DISCOVERY_INPUT_BYTES:
            yield start,position,replace(bound,proofs=tuple(proofs))
            start=position;size=0;proofs=[]
        proofs.extend(bound.proofs);size+=cost
    if start<len(values):
        yield start,len(values),replace(bound,proofs=tuple(proofs))


def discover_associations(*,root,binding_sha256,admission_completion_sha256,expected_rows,
        resources,clinical_session,reference_role,reference_identity_id,reference_manifest_sha256,execution_config,
        admission_owner,profile_reader=None,regulatory_index=None,progress_callback=None):
    """Produce one authenticated projection; partial work is never resumable."""
    from genosapiens.wgs_admission_projection import AuthenticatedAdmissionOwner
    if (type(admission_owner) is not AuthenticatedAdmissionOwner
            or admission_owner.completion.sha256!=admission_completion_sha256):
        raise admission.CandidateAdmissionError('discovery admission owner differs')
    domain_rows,domain_ordinals,domain_occurrences,profile_sha256=_discovery_domain(admission_owner,profile_reader)
    if expected_rows!=domain_rows:
        raise admission.CandidateAdmissionError('discovery profile row count differs')
    # Empty domains still require the installed source binding; absence must
    # never masquerade as a complete zero-match report.
    _regulatory_batch((),sources=resources.sources,regulatory_index=regulatory_index)
    if (profile_reader is not None and resources.sources.regulatory_index_manifest_sha256 is not None
            and profile_reader.identity['regulatory_index_sha256']!=resources.sources.regulatory_index_manifest_sha256):
        raise admission.CandidateAdmissionError('discovery profile regulatory binding differs')
    admission_owner._assert_stable()
    root=Path(root)
    root.mkdir(mode=0o700,parents=True,exist_ok=False)
    from genosapiens.wgs_contexts import _open_index
    index=_open_index(root/'targets.sqlite',execution_config)
    index.execute('CREATE TABLE targets(content BLOB PRIMARY KEY)')
    writer=_ShardWriter(root,'associations',max_rows=execution_config.source_shards.max_uncompressed_rows,
        max_bytes=execution_config.source_shards.max_uncompressed_bytes)
    subjects=_ShardWriter(root,'context_subjects',max_rows=execution_config.source_shards.max_uncompressed_rows,
        max_bytes=min(execution_config.source_shards.max_uncompressed_bytes,MAX_ROW_BYTES))
    count=association_count=previous=subject_occurrences=part_count=regulatory_inputs=0
    logical_stream=hashlib.sha256()
    ordinals=hashlib.sha256()
    batches=0;started=time.monotonic();cpu_started=time.process_time();last_progress=0.0
    timing=dict(index_commits=0,index_commit_nanoseconds=0,seal_nanoseconds=0,
        clinical_query_nanoseconds=0,clinical_batches=0,write_nanoseconds=0,
        regulatory_query_nanoseconds=0,regulatory_batches=0)
    decoder=None
    profile_source=None
    failed=False
    rich_before={key:value for key,value in admission_owner.metrics.items() if key.endswith('_discovery')}
    def progress(phase,*,force=False,target_count=None):
        nonlocal last_progress
        now=time.monotonic()
        if not force and now-last_progress<1: return
        last_progress=now
        content=_bytes(dict(schema='genosapiens.wgs-association-progress',version=1,phase=phase,
            expected_assessed_inputs=expected_rows,discovered_inputs=count,associations=association_count,
            logical_subjects=count,physical_parts=part_count,assessed_occurrences=subject_occurrences,
            discovery_batches=batches,regulatory_inputs_completed=regulatory_inputs,
            reached_targets=target_count,elapsed_seconds=now-started,
            coordinator_cpu_seconds=time.process_time()-cpu_started,
            timing=timing,admission_owner=dict(admission_owner.metrics),
            query_cache=None if getattr(clinical_session,'query_cache',None) is None else dict(clinical_session.query_cache.metrics),
            association_shards=dict(writer.metrics),subject_shards=dict(subjects.metrics),
            admission_completion_sha256=admission_completion_sha256,
            profile_manifest_sha256=profile_sha256,
            decoder=None if decoder is None else decoder.metrics,
            admission_reads={key:value-rich_before.get(key,0) for key,value in admission_owner.metrics.items()
                if key.endswith('_discovery')}))
        temporary=root/'progress.partial'
        fd=os.open(temporary,os.O_WRONLY|os.O_CREAT|os.O_TRUNC|os.O_NOFOLLOW,0o600)
        with os.fdopen(fd,'wb') as handle: handle.write(content)
        os.replace(temporary,root/'progress.json')
        from genosapiens.wgs_runtime import analysis_progress
        subphase=phase if phase in {'clinical_query','regulatory_query','writing_associations','complete','failed'} else 'working'
        analysis_progress(progress_callback, 'association_discovery', subphase,
            counts={'assessed_completed':count,'assessed_total':expected_rows,
                'associations_completed':association_count,'regulatory_inputs_completed':regulatory_inputs},
            unit='inputs',completed=regulatory_inputs if phase=='regulatory_query' else count,
            total=expected_rows,force=force)
    def regulatory_input_completed(number):
        nonlocal regulatory_inputs
        regulatory_inputs+=number
        progress('regulatory_query')
    try:
        progress('discovery',force=True)
        profile_source=None if profile_reader is None else (
            row for row in profile_reader.iter_rows() if row.decision['disposition'] in {'selected','unresolved'})
        from genosapiens.wgs_context_pipeline import AdmissionDiscoveryExecutor,admission_windows
        decoder=AdmissionDiscoveryExecutor(execution_config=execution_config, total_inputs=expected_rows,
            progress=lambda: progress('worker_startup' if decoder is None else decoder.phase))
        with closing(decoder),closing(admission_windows(admission_owner.reader,
                max_rows=execution_config.public_queries.batch_rows,decode_window=decoder.decode_window,
                profile_reader=profile_reader)) as windows:
            for rows in windows:
                values=tuple(row.value for row in rows)
                profile_rows=None if profile_source is None else tuple(next(profile_source,None) for _ in rows)
                progress('clinical_query')
                query_started=time.monotonic_ns()
                try:clinical=clinical_session.association_batch(tuple(value.allele for value in values if value.allele is not None))
                finally:timing['clinical_query_nanoseconds']+=time.monotonic_ns()-query_started
                timing['clinical_batches']+=1
                # Validate the complete response before partitioning it. A
                # foreign proof must fail, not disappear through the slice.
                if clinical is not None:
                    if type(clinical) is not admission.ShardedClinicalAssociations:
                        raise RuntimeError('discovery clinical results are invalid')
                    window_subjects={value.subject_id for value in values}
                    if any(p.proof.allele_id not in window_subjects for p in clinical.proofs):
                        raise RuntimeError('discovery clinical result is foreign')
                with closing(_regulatory_windows(values,sources=resources.sources,
                        regulatory_index=regulatory_index,on_input=regulatory_input_completed,
                        profile_rows=profile_rows)) as regulatory_windows:
                    while True:
                        progress('regulatory_query')
                        query_started=time.monotonic_ns()
                        try:window=next(regulatory_windows,None)
                        finally:timing['regulatory_query_nanoseconds']+=time.monotonic_ns()-query_started
                        if window is None:break
                        first,last,regulatory=window
                        if regulatory is not None:timing['regulatory_batches']+=1
                        selected=rows[first:last]
                        subjects_in_window={row.value.subject_id for row in selected}
                        selected_clinical=None if clinical is None else replace(clinical,
                            proofs=tuple(p for p in clinical.proofs if p.proof.allele_id in subjects_in_window))
                        with closing(decoder.build_window(selected,clinical=selected_clinical,resources=resources,
                                admission_owner=admission_owner,binding_sha256=binding_sha256,
                                admission_completion_sha256=admission_completion_sha256,
                                reference_role=reference_role,reference_identity_id=reference_identity_id,
                                reference_manifest_sha256=reference_manifest_sha256,regulatory=regulatory,
                                profile_rows=None if profile_rows is None else profile_rows[first:last])) as built:
                            for batch in built:
                                progress('writing_associations')
                                previous,n,a,o,p=_write_discovery_batch(batch,writer=writer,subjects=subjects,
                                    index=index,previous=previous,logical_stream=logical_stream,ordinals=ordinals,timing=timing)
                                count+=n;association_count+=a;subject_occurrences+=o;part_count+=p
                        batches+=1;progress('discovery')
            decoder.finish(expected_rows)
        if profile_source is not None and next(profile_source,None) is not None:
            raise admission.CandidateAdmissionError('discovery profile target tail is incomplete')
        # The expanded oracle has no producer ordinal digest. Its necessary
        # admission traversal checks each assessed record and derives the digest
        # above; opening a completed owner must not pre-scan that same domain.
        _,domain_ordinals,_,_=_discovery_domain(admission_owner,profile_reader)
        if (count!=expected_rows or domain_ordinals is None or ordinals.hexdigest()!=domain_ordinals
                or subject_occurrences!=domain_occurrences):
            raise admission.CandidateAdmissionError('association discovery assessed coverage differs')
        from genosapiens.evidence_model import _typed
        targets=[];size=0
        for content, in index.execute('SELECT content FROM targets ORDER BY content'):
            size+=len(content)
            if size>MAX_TARGET_BYTES:
                raise admission.CandidateAdmissionError('association target closure exceeds metadata capacity')
            targets.append(_typed(json.loads(content),admission.AssociationTarget,decode=True))
        seal_started=time.monotonic_ns()
        manifest=AssociationProjection(binding_sha256,admission_completion_sha256,
            clinical_session.indexes.manifest_sha256,count,association_count,writer.close(),
            writer.stream.hexdigest(),tuple(targets),count,subject_occurrences,subjects.close(),
            subjects.stream.hexdigest(),ordinals.hexdigest(),admission_owner.eligibility_manifest.sha256,
            admission_owner.row_count,reference_role,reference_identity_id,reference_manifest_sha256,
            part_count,logical_stream.hexdigest(),profile_manifest_sha256=profile_sha256,
            regulatory_index_manifest_sha256=resources.sources.regulatory_index_manifest_sha256,
            association_scope_sha256=None if profile_sha256 is None else ASSOCIATION_SCOPE_SHA256)
        if len(manifest.as_bytes())>MAX_DISCOVERY_MANIFEST_BYTES:
            raise WgsScientificError("complete discovery manifest exceeds metadata capacity",category='response_capacity')
        admission_owner._assert_stable()
        _discovery_domain(admission_owner,profile_reader)
        _write_manifest(root,'manifest.json',manifest.as_bytes(),max_bytes=MAX_DISCOVERY_MANIFEST_BYTES)
        timing['seal_nanoseconds']+=time.monotonic_ns()-seal_started
        progress('complete',force=True,target_count=len(targets))
        return manifest
    except BaseException:
        failed=True
        try:(root/'manifest.json').unlink(missing_ok=True)
        except OSError:pass
        try:progress('failed',force=True)
        except BaseException:pass
        raise
    finally:
        if profile_source is not None:profile_source.close()
        for cleanup in (writer.close,subjects.close,index.close,lambda:(root/'targets.sqlite').unlink(missing_ok=True)):
            try:cleanup()
            except BaseException:
                if not failed:raise
