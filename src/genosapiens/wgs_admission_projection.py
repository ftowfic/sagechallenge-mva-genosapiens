"""Private bounded admission bytes and storage; published evidence stays expanded.

The row owner validates/encodes once. SQLite and IPC consumers authenticate
lengths, hashes, bindings and scalar coverage before using canonical bytes.
No predecessor partial store is a supported input or migration source.
"""
from __future__ import annotations

from genosapiens.wgs_read_cache import ReadFileCache

from collections import OrderedDict, deque
from contextlib import closing, contextmanager
from dataclasses import asdict, dataclass, field, fields, replace
import hashlib
import json
import os
from pathlib import Path
import queue
import sqlite3
import stat
import threading
import time
from types import MappingProxyType, MethodType
import zlib

from genosapiens.annotation_admission import AnnotationAssociationProof
from genosapiens.bounded_wgs_config import BoundedWgsExecutionConfig
from genosapiens.evidence_model import _typed, _typed_value, _validation_session, SourcedRecord, IdentityMapping
from genosapiens.evidence_serialization import _plain, _owned_canonical_bytes
from genosapiens.wgs_admission import (
    LedgerAdmissionInput, WgsAdmissionError, load_ledger_admission_input,
    _owned_admission_input_bytes,
)
from genosapiens.wgs_ledger import MAX_ROW_BYTES, canonical
from genosapiens.wgs_analysis_eligibility import (
    AdmissionScreeningContext, METHOD_SHA256, load_admission_screening_context,
    AdmissionScreeningFacts, AnalysisEligibilityDecision, AnalysisEligibilityManifest,
    load_admission_screening_facts, load_analysis_eligibility_decision,
    load_analysis_eligibility_manifest,
    _screening_facts_from_mapping, _screening_facts_mapping,
    _eligibility_decision_from_mapping, _eligibility_decision_mapping, _unique,
)


class AdmissionProjectionError(RuntimeError):
    """Static downstream cache diagnostics, never source values or paths."""


@_validation_session(isolated=True)
def _admission_proofs_bytes(proofs: tuple[AnnotationAssociationProof, ...]) -> bytes:
    if (type(proofs) is not tuple
            or any(type(value) is not AnnotationAssociationProof for value in proofs)):
        raise AdmissionProjectionError("admission proof projection is invalid")
    return _owned_admission_proofs_bytes(_typed(proofs, tuple[AnnotationAssociationProof, ...]))


def _owned_admission_proofs_bytes(proofs):
    """Intern only the record owner's fully validated proof copies."""
    # A CSQ with several consequence terms produces several distinct proofs of
    # the same source record. Store that source and the common native identity
    # and replay provenance once, without discarding any proof or raising the
    # existing per-row byte limit.
    tables = {name: [] for name in ("records", "identities")}
    keys = {name: {} for name in tables}
    seen = {name: {} for name in tables}

    def intern(name, value):
        cached = seen[name].get(id(value))
        if cached is not None and cached[0] is value:
            return cached[1]
        key = _owned_canonical_bytes(value)
        index = keys[name].get(key)
        if index is None:
            index = len(tables[name])
            keys[name][key] = index
            tables[name].append(value)
        seen[name][id(value)] = (value, index)
        return index

    rows = []
    for proof in proofs:
        row = {field.name: getattr(proof, field.name) for field in fields(proof)
               if field.name not in {"native", "source_refs", "provenance"}}
        row.update(native=intern("identities", proof.native),
            source_refs=[intern("records", source) for source in proof.source_refs],
            provenance=intern("records", proof.provenance))
        rows.append(row)
    content = _owned_canonical_bytes({"schema": "genosapiens.wgs-enriched-admission-proofs",
        "version": 2, "classification": "controlled", "llm_policy": "no-llm",
        **tables, "proofs": rows})
    if len(content) > MAX_ROW_BYTES:
        raise AdmissionProjectionError("admission proof projection exceeds capacity")
    return content


@_validation_session(isolated=True)
def _load_admission_proofs(content: bytes) -> tuple[AnnotationAssociationProof, ...]:
    try:
        if type(content) is not bytes or len(content) > MAX_ROW_BYTES:
            raise ValueError
        value = json.loads(content)
        if (canonical(value) != content or set(value) != {
                "schema", "version", "classification", "llm_policy", "proofs", "records", "identities"}
                or value["schema"] != "genosapiens.wgs-enriched-admission-proofs"
                or type(value["version"]) is not int or value["version"] != 2 or value["classification"] != "controlled"
                or value["llm_policy"] != "no-llm"
                or any(type(value[name]) is not list for name in ("proofs", "records", "identities"))):
            raise ValueError
        # One isolated external decode owns all children. Validate even unused
        # table entries, then let the strict walker construct complete proofs
        # from these same raw owners without restarting validation per field.
        records, identities = value["records"], value["identities"]
        for row in records:
            _typed_value(row, SourcedRecord, decode=True)
        for row in identities:
            _typed_value(row, IdentityMapping, decode=True)
        def reference(table, index):
            if type(index) is not int or not 0 <= index < len(table):
                raise ValueError
            return table[index]
        names = {field.name for field in fields(AnnotationAssociationProof)}
        decoded = []
        for row in value["proofs"]:
            if type(row) is not dict or set(row) != names or type(row["source_refs"]) is not list:
                raise ValueError
            arguments = {**row, 'native': reference(identities, row["native"]),
                'source_refs': [reference(records, index) for index in row["source_refs"]],
                'provenance': reference(records, row["provenance"])}
            decoded.append(_typed_value(arguments, AnnotationAssociationProof, decode=True))
        proofs = tuple(decoded)
        if _owned_admission_proofs_bytes(proofs) != content:
            raise ValueError
        return proofs
    except (ValueError, TypeError, KeyError, RuntimeError, UnicodeError):
        raise AdmissionProjectionError("admission proof projection is invalid") from None



def _pack_admission_blob(content: bytes) -> bytes:
    """Keep repeated publisher fields off disk without changing scientific bytes."""
    if type(content) is not bytes or not 0 < len(content) <= MAX_ROW_BYTES:
        raise AdmissionProjectionError("admission cache row exceeds capacity")
    return zlib.compress(content, level=1)


def _unpack_admission_blob(content: bytes) -> bytes:
    try:
        # DEFLATE framing can add a small overhead to incompressible rows. The
        # scientific byte limit applies to the uncompressed value, unchanged.
        if type(content) is not bytes or not 0 < len(content) <= MAX_ROW_BYTES + 65536:
            raise ValueError
        decoder = zlib.decompressobj()
        value = decoder.decompress(content, MAX_ROW_BYTES + 1)
        if (not 0 < len(value) <= MAX_ROW_BYTES or not decoder.eof
                or decoder.unused_data or decoder.unconsumed_tail):
            raise ValueError
        return value
    except (ValueError, zlib.error):
        raise AdmissionProjectionError("admission cache row is invalid") from None



_BINDING_FIELDS = {
    'schema', 'version', 'prepared_index_sha256', 'source_manifest_sha256',
    'working_manifest_sha256', 'source_routing_sha256', 'plan_sha256',
    'annotation_stage_sha256', 'resource_session_sha256', 'bundle_manifest_sha256',
    'execution_limits_sha256', 'read_intake_sha256', 'consumer_executable_sha256',
    'completed_analysis_input_sha256', 'completed_annotation_prefix_sha256',
    'codec', 'identities', 'occurrences', 'sequence_requests',
    'read_intake_source_bindings_sha256',
    'screening_context_sha256', 'eligibility_policy_sha256', 'analysis_method_sha256',
}
_NULLABLE_BINDINGS = {'read_intake_sha256', 'completed_analysis_input_sha256',
    'completed_annotation_prefix_sha256', 'read_intake_source_bindings_sha256'}


def _sha(value):
    return type(value) is str and len(value) == 64 and all(c in '0123456789abcdef' for c in value)


def _count(value, minimum=0):
    return type(value) is int and value >= minimum


def _binding(content):
    try:
        if type(content) is not bytes or not 0 < len(content) <= MAX_ROW_BYTES:
            raise ValueError
        value = json.loads(content)
        if (type(content) is not bytes or canonical(value) != content
                or set(value) != _BINDING_FIELDS
                or value['schema'] != 'genosapiens.wgs-enriched-admission-projection'
                or type(value['version']) is not int or value['version'] != 5
                or value['codec'] != 'canonical-input-v1/proofs-v2/facts-v1/zlib-1/layout-v5'
                or value['analysis_method_sha256'] != METHOD_SHA256
                or any(not _count(value[k]) for k in ('identities','occurrences','sequence_requests'))
                or value['sequence_requests'] > value['identities']
                or value['occurrences'] < value['identities']):
            raise ValueError
        for key in _BINDING_FIELDS:
            if key.endswith('_sha256') and not (key in _NULLABLE_BINDINGS and value[key] is None) and not _sha(value[key]):
                raise ValueError
        if (value['read_intake_sha256'] is None) != (value['read_intake_source_bindings_sha256'] is None):
            raise ValueError
        return value
    except (TypeError, ValueError, KeyError):
        raise AdmissionProjectionError('admission projection binding is invalid') from None


def _row_parents(parents):
    intake = parents['read_intake_sha256']
    if intake is not None:
        return intake,intake,parents['read_intake_source_bindings_sha256']
    return (parents['source_manifest_sha256'],parents['working_manifest_sha256'],
        parents['source_routing_sha256'])


def _require_fact_parents(facts, parents):
    if ((facts.source_manifest_sha256, facts.working_manifest_sha256,
            facts.source_routing_sha256) != _row_parents(parents)
            or facts.merged_intake_sha256 != parents['read_intake_sha256']
            or facts.eligibility_policy_sha256 != parents['eligibility_policy_sha256']):
        raise AdmissionProjectionError('admission screening parents are invalid')


def _require_screening_context(screening, parents):
    if (type(screening) is not AdmissionScreeningContext
            or screening.sha256 != parents['screening_context_sha256']
            or screening.policy.sha256 != parents['eligibility_policy_sha256']
            or (None if screening.merged_intake is None else screening.merged_intake.sha256)
                != parents['read_intake_sha256']):
        raise AdmissionProjectionError('admission screening context changed')


def _require_execution_config(config, parents):
    from genosapiens.bounded_wgs_config import bounded_wgs_execution_config_as_json
    if (type(config) is not BoundedWgsExecutionConfig
            or hashlib.sha256(bounded_wgs_execution_config_as_json(config)).hexdigest()
                != parents['execution_limits_sha256']):
        raise AdmissionProjectionError('admission execution limits changed')


@dataclass(frozen=True, slots=True, repr=False)
class EncodedAdmissionRow:
    """One validated canonical row; no decoded graph is retained by the writer."""
    retained_ordinal: int
    identity_id: str
    allele_id: str | None
    input_blob: bytes
    proof_blob: bytes
    input_size: int
    proof_size: int
    input_sha256: str
    proof_sha256: str
    proof_count: int
    source_occurrences: int
    sequence_requests: int
    gene_ids: tuple[str, ...]
    binding_sha256: str
    screening_facts: bytes

    def __post_init__(self):
        self.authenticated_bytes()

    def authenticated_bytes(self):
        content, proofs, _ = self._authenticated_parts()
        return content, proofs

    def _authenticated_parts(self, *, binding=None):
        """One synchronous authentication of blobs and their compact linkage."""
        content = _unpack_admission_blob(self.input_blob)
        proofs = _unpack_admission_blob(self.proof_blob)
        if (not _count(self.retained_ordinal, 1)
                or type(self.identity_id) is not str or not self.identity_id
                or self.allele_id is not None and (type(self.allele_id) is not str or not self.allele_id)
                or not _count(self.input_size,1) or not _count(self.proof_size,1)
                or len(content) != self.input_size or len(proofs) != self.proof_size
                or not all(_sha(s) for s in (self.input_sha256,self.proof_sha256,self.binding_sha256))
                or hashlib.sha256(content).hexdigest() != self.input_sha256
                or hashlib.sha256(proofs).hexdigest() != self.proof_sha256
                or not _count(self.proof_count) or not _count(self.source_occurrences,1)
                or type(self.sequence_requests) is not int or self.sequence_requests not in (0,1)
                or type(self.gene_ids) is not tuple
                or any(type(g) is not str or not g.isascii() or not g for g in self.gene_ids)
                or tuple(sorted(set(self.gene_ids))) != self.gene_ids):
            raise AdmissionProjectionError('encoded admission row is invalid')
        facts = self.authenticated_facts(binding=binding)
        return content, proofs, facts

    def authenticated_facts(self, *, binding=None):
        facts = load_admission_screening_facts(self.screening_facts)
        if (facts.retained_ordinal != self.retained_ordinal or facts.identity_id != self.identity_id
                or facts.allele_id != self.allele_id or facts.gene_ids != self.gene_ids
                or len(facts.source_allele_ids) != self.source_occurrences
                or int(facts.transport == 'sequence_vep') != self.sequence_requests):
            raise AdmissionProjectionError('admission screening row linkage is invalid')
        if binding is not None:
            parents = _binding(binding)
            if self.binding_sha256 != hashlib.sha256(binding).hexdigest():
                raise AdmissionProjectionError('admission screening parents are invalid')
            _require_fact_parents(facts, parents)
        return facts

    @property
    def input_bytes(self):
        return _unpack_admission_blob(self.input_blob)

    @property
    def proof_bytes(self):
        return _unpack_admission_blob(self.proof_blob)

    @classmethod
    def encode(cls, value, proofs, *, binding, screening):
        # Ordinary callers always begin a fresh validation operation. Workers
        # call the private encoder within their synchronous construction owner.
        with _validation_session(isolated=True):
            return cls._encode_owned(value,proofs,binding=binding,screening=screening)

    @classmethod
    def _encode_owned(cls, value, proofs, *, binding, screening):
        parents = _binding(binding)
        _require_screening_context(screening, parents)
        source,working,routing = _row_parents(parents)
        if (type(value) is not LedgerAdmissionInput
                or value.parent_source_manifest_sha256 != source
                or value.parent_working_set_manifest_sha256 != working
                or value.native_source_contract_sha256 != routing):
            raise AdmissionProjectionError('admission row parents are invalid')
        # The only final encoder owns both previous capacity checks. Copies
        # validated during this operation cannot leak into subsequent rows.
        with _validation_session():
            value = _typed(value, LedgerAdmissionInput)
            proofs = _typed(proofs, tuple[AnnotationAssociationProof, ...])
            content = _owned_admission_input_bytes(value)
            proof_content = _owned_admission_proofs_bytes(proofs)
            facts = screening.derive(value)
        genes = tuple(sorted({t.ensembl_gene_id for t in
            (() if value.allele is None else value.allele.transcript_consequences or ())
            if t.ensembl_gene_id is not None}))
        return cls(value.retained_ordinal, value.identity_id,
            None if value.allele is None else value.allele.allele_id,
            _pack_admission_blob(content), _pack_admission_blob(proof_content),len(content),len(proof_content),
            hashlib.sha256(content).hexdigest(), hashlib.sha256(proof_content).hexdigest(),
            len(proofs), len(value.source_allele_ids), int(value.transport == 'sequence_vep'),
            genes, hashlib.sha256(binding).hexdigest(), facts.as_bytes())

    @property
    def ref(self):
        return (self.binding_sha256,self.retained_ordinal,self.identity_id,
            self.allele_id,self.input_sha256,self.input_size)

    def decode(self, *, binding):
        value = self.decode_input(binding=binding)
        proofs = _load_admission_proofs(self.proof_bytes)
        if len(proofs) != self.proof_count:
            raise AdmissionProjectionError('admission row proof count is invalid')
        return value,proofs

    def decode_input(self, *, binding):
        parents = _binding(binding)
        value, _ = self._decode_input_with_content(parents, hashlib.sha256(binding).hexdigest())
        return value

    def _decode_input_with_content(self, parents, binding_sha256):
        """Carry canonical bytes from strict decoding into one private owner."""
        content, _, facts = self._authenticated_parts()
        _require_fact_parents(facts, parents)
        source,working,routing = _row_parents(parents)
        if binding_sha256 != self.binding_sha256:
            raise AdmissionProjectionError('admission row binding is invalid')
        value = load_ledger_admission_input(content)
        genes = tuple(sorted({t.ensembl_gene_id for t in
            (() if value.allele is None else value.allele.transcript_consequences or ())
            if t.ensembl_gene_id is not None}))
        if (value.retained_ordinal != self.retained_ordinal or value.identity_id != self.identity_id
                or (None if value.allele is None else value.allele.allele_id) != self.allele_id
                or len(value.source_allele_ids) != self.source_occurrences
                or int(value.transport == 'sequence_vep') != self.sequence_requests or genes != self.gene_ids
                or value.parent_source_manifest_sha256 != source
                or value.parent_working_set_manifest_sha256 != working
                or value.native_source_contract_sha256 != routing):
            raise AdmissionProjectionError('admission row metadata is invalid')
        return value, content


@dataclass(frozen=True, slots=True, repr=False)
class _ReaderRow(EncodedAdmissionRow):
    """Fresh disk row's authenticated bytes, held only in a bounded window.

    The ordinary producer/IPC row remains unchanged and always reauthenticates.
    This class is constructed only from the reader's scalar SQLite projection.
    """
    _input_content: bytes = field(init=False, compare=False)
    _proof_content: bytes = field(init=False, compare=False)

    def __post_init__(self):
        content, proofs = self.authenticated_bytes()
        object.__setattr__(self, '_input_content', content)
        object.__setattr__(self, '_proof_content', proofs)

    def __eq__(self, other):
        if not isinstance(other, EncodedAdmissionRow):
            return NotImplemented
        return all(getattr(self, f.name) == getattr(other, f.name) for f in fields(EncodedAdmissionRow))

    @property
    def input_bytes(self):
        return self._input_content

    @property
    def proof_bytes(self):
        return self._proof_content


def _reader_row(values, binding_sha256):
    """The same strict disk-row boundary for local and live worker decoding."""
    try:
        ordinal,identity,allele,content,proofs,size,proof_size,digest,proof_digest,count,occurrences,requests,genes,facts = values
        gene_ids = json.loads(genes)
        if type(gene_ids) is not list or canonical(gene_ids) != genes:
            raise ValueError
        return _ReaderRow(ordinal,identity,allele,content,proofs,size,proof_size,digest,proof_digest,
            count,occurrences,requests,tuple(gene_ids),binding_sha256,facts)
    except (TypeError, ValueError):
        raise AdmissionProjectionError('admission projection row is invalid') from None


def _settings(db, config, *, readonly):
    db.enable_load_extension(False)
    db.execute('PRAGMA trusted_schema=OFF')
    db.execute('PRAGMA mmap_size=0')
    db.execute('PRAGMA temp_store=FILE')
    db.execute('PRAGMA foreign_keys=ON')
    db.execute(f'PRAGMA cache_size=-{max(1,config.sqlite.cache_bytes//1024)}')
    if readonly:
        db.execute('PRAGMA query_only=ON')
    else:
        db.execute('PRAGMA journal_mode=DELETE')
        db.execute('PRAGMA synchronous=FULL')


def _private_path(path, *, exists):
    path = Path(path)
    if not path.is_absolute() or any(p.is_symlink() for p in (path,*path.parents)):
        raise AdmissionProjectionError('admission projection path is invalid')
    parent = path.parent.stat()
    if (not stat.S_ISDIR(parent.st_mode) or parent.st_uid != os.getuid()
            or stat.S_IMODE(parent.st_mode) != 0o700):
        raise AdmissionProjectionError('admission projection parent is not private')
    if exists:
        meta = path.lstat()
        if (not stat.S_ISREG(meta.st_mode) or meta.st_uid != os.getuid()
                or stat.S_IMODE(meta.st_mode) != 0o600 or meta.st_nlink != 1 or meta.st_size <= 0):
            raise AdmissionProjectionError('admission projection file is not private')
    return path


_COLUMNS = ('ordinal,identity_id,allele_id,admission_bytes,proofs_bytes,input_size,proof_size,'
    'input_sha256,proof_sha256,proof_count,occurrences,requests,gene_ids,screening_facts')
_INPUT_COLUMNS = ('ordinal,identity_id,allele_id,admission_bytes,input_size,input_sha256,'
    'occurrences,requests,gene_ids')
_READER_CACHE_BYTES = 64*1024**2
_READER_CACHE_ENTRIES = 1024
_READ_PURPOSES = ('profile','discovery','context','selected','audit')
_COMPACT_READ_BYTES = 8*1024**2
_COMPACT_READ_ROWS = 128


@contextmanager
def _compact_read_ahead(db, progress, metrics, *, cache_release=None):
    """Temporarily transfer the idle connection to one bounded I/O thread.

    Only this thread uses SQLite until it has joined. Validation stays in the
    caller. One queued batch, one producer batch, one consumer batch and one
    overflow row together retain at most 32 MiB of compact bytes. No rich column
    is selected, and exhaustion still has to authenticate the complete closure.
    """
    pending=queue.Queue(maxsize=1);stop=threading.Event();done=threading.Event();errors=[]
    def put(batch,size):
        while not stop.is_set():
            try:
                pending.put(tuple(batch),timeout=.1)
                metrics['compact_read_batches']+=1
                metrics['max_compact_read_batch_bytes']=max(metrics['max_compact_read_batch_bytes'],size)
                if cache_release is not None:cache_release(force=True)
                return True
            except queue.Full:pass
        return False
    def read():
        try:
            db.set_progress_handler(lambda:int(stop.is_set()),100000)
            with closing(db.execute('SELECT ordinal,canonical_bytes,sha256,disposition,accounting_bytes '
                    'FROM screening ORDER BY ordinal')) as cursor:
                batch=[];size=0
                for row in cursor:
                    if stop.is_set():break
                    if (type(row[1]) is not bytes or not 0<len(row[1])<=MAX_ROW_BYTES
                            or type(row[4]) is not bytes or not 0<len(row[4])<=MAX_ROW_BYTES
                            or type(row[2]) is not str or len(row[2])!=64
                            or row[3] not in ('assess','exclude_common','defer_refcall')):
                        raise AdmissionProjectionError('compact read exceeds row capacity')
                    cost=len(row[1])+len(row[4])
                    if batch and (len(batch)>=_COMPACT_READ_ROWS or size+cost>_COMPACT_READ_BYTES):
                        if not put(batch,size):return
                        batch=[];size=0
                    batch.append(row);size+=cost;metrics['compact_read_rows']+=1
                if batch and not stop.is_set():put(batch,size)
        except BaseException as error:
            if not stop.is_set():errors.append(error.with_traceback(None))
        finally:
            try:db.set_progress_handler(None,0)
            except BaseException as error:
                if not stop.is_set():errors.append(error.with_traceback(None))
            finally:done.set()
    def rows():
        while True:
            try:batch=pending.get(timeout=.1)
            except queue.Empty:
                # Progress remains coordinator-owned, including disk stalls.
                progress()
                if not done.is_set():continue
                if errors:raise errors[0]
                return
            yield from batch
            del batch
    worker=threading.Thread(target=read,name='admission-compact-read')
    worker.start()
    try:
        with closing(rows()) as source:yield source
    finally:
        stop.set();db.interrupt();worker.join(timeout=30)
        if worker.is_alive():
            raise AdmissionProjectionError('compact read thread did not stop')
        metrics['compact_read_stopped']=True
        db.set_progress_handler(lambda:(progress(),0)[1],100000)

_VALIDATED_ROW_SEAL = object()


@dataclass(frozen=True, slots=True, repr=False)
class AdmissionProjectionRecord:
    """Strictly decoded current row and its original canonical byte owner.

    Retain at most the caller's bounded batch/context. This process-local seal
    is never serialized, accepted from workers, or used to authorize disk data.
    """
    encoded: EncodedAdmissionRow
    value: LedgerAdmissionInput
    proofs: tuple[AnnotationAssociationProof, ...]
    seal: object
    compact: AdmissionScreeningRecord | None = field(default=None,compare=False,repr=False)
    compact_bytes: int = field(default=0,compare=False,repr=False)

    @property
    def ref(self):
        return self.encoded.ref

    def require_input(self, value):
        if self.seal is not _VALIDATED_ROW_SEAL or self.value is not value:
            raise AdmissionProjectionError('validated admission row does not own input')


@dataclass(frozen=True, slots=True, repr=False)
class AdmissionScreeningRecord:
    """Compact producer attestation; authority comes from the anchored completion."""
    facts: AdmissionScreeningFacts
    decision: AnalysisEligibilityDecision
    input_size: int
    proof_size: int
    proof_count: int

    def __post_init__(self):
        self._validate()
        self.facts.__post_init__(); self.decision.__post_init__()

    def _validate(self):
        if (type(self.facts) is not AdmissionScreeningFacts
                or type(self.decision) is not AnalysisEligibilityDecision
                or not _count(self.input_size, 1) or self.input_size > MAX_ROW_BYTES
                or not _count(self.proof_size, 1) or self.proof_size > MAX_ROW_BYTES
                or not _count(self.proof_count)):
            raise AdmissionProjectionError('invalid admission screening record')
        facts, decision = self.facts, self.decision
        if ((facts.retained_ordinal, facts.identity_id, facts.allele_id, facts.source_allele_ids,
                facts.source_manifest_sha256, facts.eligibility_policy_sha256) !=
                (decision.retained_ordinal, decision.identity_id, decision.allele_id, decision.source_allele_ids,
                    decision.source_manifest_sha256, decision.eligibility_policy_sha256)
                or facts.reference_manifest_sha256s and
                    decision.reference_manifest_sha256 not in facts.reference_manifest_sha256s
                or facts.merged_intake_sha256 is not None and decision.disposition != 'assess'
                or facts.refcall_can_defer != (decision.disposition == 'defer_refcall')):
            raise AdmissionProjectionError('screening decision does not own facts')

    @classmethod
    def from_encoded(cls, row, decision):
        if not isinstance(row, EncodedAdmissionRow):
            raise AdmissionProjectionError('screening decision does not own blobs')
        return cls._from_authenticated_facts(row, decision, row.authenticated_facts())

    @classmethod
    def _from_authenticated_facts(cls, row, decision, facts):
        # Only the current synchronous blob/fact authentication passes this
        # private value. The constructor still checks facts, decision and all
        # their semantic links; nothing is retained across writer operations.
        if (not isinstance(row, EncodedAdmissionRow)
                or type(decision) is not AnalysisEligibilityDecision
                or (row.binding_sha256, row.input_sha256, row.proof_sha256) !=
                    (decision.admission_binding_sha256, decision.input_sha256, decision.proof_sha256)):
            raise AdmissionProjectionError('screening decision does not own blobs')
        return cls(facts, decision, row.input_size, row.proof_size, row.proof_count)

    @property
    def ref(self):
        value = self.decision
        return (value.admission_binding_sha256, value.retained_ordinal, value.identity_id,
            value.allele_id, value.input_sha256, self.input_size)

    def as_bytes(self):
        content = canonical(dict(facts=_screening_facts_mapping(self.facts,json.loads(self.facts.source_call_facts)),
            decision=_eligibility_decision_mapping(self.decision), input_size=self.input_size,
            proof_size=self.proof_size, proof_count=self.proof_count))
        if len(content) > MAX_ROW_BYTES:
            raise AdmissionProjectionError('screening record exceeds row capacity')
        return content


@_validation_session(isolated=True)
def load_admission_screening_record(content):
    try:
        if type(content) is not bytes or not 0 < len(content) <= MAX_ROW_BYTES:
            raise ValueError
        data = json.loads(content, object_pairs_hook=_unique)
        if type(data) is not dict or set(data) != set(AdmissionScreeningRecord.__dataclass_fields__):
            raise ValueError
        facts = _screening_facts_from_mapping(data['facts'])
        decision = _eligibility_decision_from_mapping(data['decision'])
        # Keep the nested public decision byte bound as well as the outer bound.
        if len(decision.as_bytes()) > 1024*1024: raise ValueError
        record = object.__new__(AdmissionScreeningRecord)
        for name,item in dict(data,facts=facts,decision=decision).items():
            object.__setattr__(record,name,item)
        record._validate()
        if canonical(dict(facts=_screening_facts_mapping(facts,data['facts']['source_call_facts']),
                decision=_eligibility_decision_mapping(decision),input_size=record.input_size,
                proof_size=record.proof_size,proof_count=record.proof_count)) != content:
            raise ValueError
        return record
    except (TypeError, ValueError, KeyError, UnicodeError):
        raise AdmissionProjectionError('invalid serialized screening record') from None


@dataclass(frozen=True, slots=True, repr=False)
class AdmissionAccountingRow:
    """Original membership and counts, authenticated against the compact row."""
    ref: tuple
    source_allele_ids: tuple[str, ...]
    proof_size: int
    proof_count: int
    sequence_requests: int
    disposition: str

    def __post_init__(self):
        if (type(self.ref) is not tuple or len(self.ref)!=6 or not _sha(self.ref[0])
                or not _count(self.ref[1],1) or type(self.ref[2]) is not str or not self.ref[2]
                or self.ref[3] is not None and (type(self.ref[3]) is not str or not self.ref[3])
                or not _sha(self.ref[4]) or not _count(self.ref[5],1) or self.ref[5]>MAX_ROW_BYTES
                or type(self.source_allele_ids) is not tuple or not self.source_allele_ids
                or any(type(key) is not str or not key for key in self.source_allele_ids)
                or len(set(self.source_allele_ids))!=len(self.source_allele_ids)
                or not _count(self.proof_size,1) or self.proof_size>MAX_ROW_BYTES or not _count(self.proof_count)
                or type(self.sequence_requests) is not int or self.sequence_requests not in (0,1)
                or type(self.disposition) is not str or self.disposition not in ('assess','exclude_common','defer_refcall')):
            raise AdmissionProjectionError('invalid scalar admission accounting')

    @classmethod
    def from_record(cls,record):
        return cls(record.ref,record.facts.source_allele_ids,record.proof_size,record.proof_count,
            int(record.facts.transport=='sequence_vep'),record.decision.disposition)

    def as_bytes(self):
        return canonical({name: getattr(self, name) for name in AdmissionAccountingRow.__dataclass_fields__})


def _load_admission_accounting(content):
    try:
        if type(content) is not bytes or not 0<len(content)<=MAX_ROW_BYTES: raise ValueError
        value=json.loads(content,object_pairs_hook=_unique)
        if type(value) is not dict or set(value)!=set(AdmissionAccountingRow.__dataclass_fields__): raise ValueError
        for field in ('ref','source_allele_ids'):
            if type(value[field]) is not list: raise ValueError
            value[field]=tuple(value[field])
        row=AdmissionAccountingRow(**value)
        if row.as_bytes()!=content: raise ValueError
        return row
    except (TypeError,ValueError,KeyError,UnicodeError):
        raise AdmissionProjectionError('invalid serialized scalar admission accounting') from None


@dataclass(frozen=True, slots=True)
class AdmissionProjectionCompletion:
    binding_sha256: str
    rows: int
    proofs: int
    occurrences: int
    requests: int
    input_bytes: int
    proof_bytes: int
    genes: int
    gene_sha256: str
    input_rows_sha256: str
    ordered_input_hashes_sha256: str
    ordered_proof_hashes_sha256: str
    screening_sha256: str
    eligibility_manifest_sha256: str
    schema: str = 'genosapiens.wgs-admission-completion'
    version: int = 1

    def __post_init__(self):
        if (self.schema != 'genosapiens.wgs-admission-completion'
                or type(self.version) is not int or self.version != 1
                or any(not _sha(getattr(self, f.name)) for f in fields(self) if f.name.endswith('_sha256'))
                or any(not _count(getattr(self, key)) for key in (
                    'rows','proofs','occurrences','requests','input_bytes','proof_bytes','genes'))
                or self.occurrences < self.rows or self.requests > self.rows
                or self.input_bytes < self.rows or self.proof_bytes < self.rows):
            raise AdmissionProjectionError('invalid admission completion')

    def as_bytes(self):
        return canonical(asdict(self))

    @property
    def sha256(self):
        return hashlib.sha256(self.as_bytes()).hexdigest()


def load_admission_projection_completion(content, *, expected_sha256):
    try:
        if (type(content) is not bytes or not 0 < len(content) <= 8192 or not _sha(expected_sha256)
                or hashlib.sha256(content).hexdigest() != expected_sha256):
            raise ValueError
        data = json.loads(content)
        if set(data) != set(AdmissionProjectionCompletion.__dataclass_fields__):
            raise ValueError
        completion = AdmissionProjectionCompletion(**data)
        if completion.as_bytes() != content:
            raise ValueError
        return completion
    except (ValueError, TypeError, KeyError):
        raise AdmissionProjectionError('admission completion is not authenticated') from None


def publish_admission_projection_completion(completion, *, state, start):
    """Anchor A after E/A closure in the existing authenticated stage journal."""
    from genosapiens.wgs_ledger import StageCommitState, publish_bytes
    if (type(completion) is not AdmissionProjectionCompletion or type(state) is not StageCommitState
            or state.binding.stage != 'admission_projection'
            or state.binding.stage_parent_sha256 != completion.binding_sha256):
        raise AdmissionProjectionError('admission completion stage binding changed')
    artifact=publish_bytes(state.stage_root,'admission-completion.json',completion.as_bytes())
    counts=dict(rows=completion.rows,occurrences=completion.occurrences,proofs=completion.proofs,requests=completion.requests)
    state.commit_partition(start,cursor=(completion.rows,),outputs=(artifact,),counts=counts)
    def closed(_connection):
        committed=state.committed('admission')
        return committed is not None and committed.counts == counts and committed.outputs == (artifact,)
    state.finish(expected_partitions=('admission',),output=artifact,
        process_id=start.process_id,closure_check=closed)
    return artifact


def read_admission_projection_completion(state):
    """Read only a completed current stage under its caller-supplied binding."""
    from genosapiens.wgs_ledger import StageCommitState, JournalEntry, load_model, read_artifact
    if type(state) is not StageCommitState or state.binding.stage != 'admission_projection':
        raise AdmissionProjectionError('admission completion stage is invalid')
    committed=state.committed('admission')
    tail=state.connection.execute('SELECT event_bytes FROM events ORDER BY ordinal DESC LIMIT 1').fetchone()
    if committed is None or tail is None:
        raise AdmissionProjectionError('admission stage is incomplete')
    event=load_model(tail[0],JournalEntry)
    if (event.event_type != 'stage_completed' or len(committed.outputs) != 1
            or event.artifacts != committed.outputs or committed.outputs[0].name != 'admission-completion.json'):
        raise AdmissionProjectionError('admission stage completion is missing')
    artifact=committed.outputs[0]
    completion=load_admission_projection_completion(read_artifact(state.stage_root,artifact),expected_sha256=artifact.sha256)
    if (completion.binding_sha256 != state.binding.stage_parent_sha256
            or committed.counts != dict(rows=completion.rows,occurrences=completion.occurrences,
                proofs=completion.proofs,requests=completion.requests)):
        raise AdmissionProjectionError('admission stage completion coverage changed')
    return completion


def _projection_operation(function):
    from functools import wraps
    import time
    @wraps(function)
    def operation(self,*args,**kwargs):
        if self.failed:
            raise AdmissionProjectionError('admission writer previously failed')
        started=time.monotonic_ns()
        try:
            return function(self,*args,**kwargs)
        except BaseException:
            self.failed=True
            raise
        finally:
            key='seal_nanoseconds' if function.__name__=='seal' else 'write_nanoseconds'
            self.metrics[key]+=time.monotonic_ns()-started
    return operation


class AdmissionProjectionWriter:
    """Sole coordinator writer. A failed or interrupted store is never resumed."""
    def __init__(self, path, *, binding, execution_config):
        self.parents = MappingProxyType(_binding(binding))
        _require_execution_config(execution_config,self.parents)
        self.binding = binding
        self.binding_sha256 = hashlib.sha256(binding).hexdigest()
        self.config = execution_config
        self.count = self.proofs = self.occurrences = self.requests = self.windows = 0
        self.input_size = self.proof_size = 0
        self.inputs_digest = hashlib.sha256()
        self.input_hashes = hashlib.sha256()
        self.proof_hashes = hashlib.sha256()
        self.screening_digest = hashlib.sha256()
        self.finished = self.failed = False
        self.connection = None
        self.metrics=dict(commits=0,commit_nanoseconds=0,write_nanoseconds=0,
            seal_nanoseconds=0,committed_inputs=0,committed_occurrences=0)
        try:
            self.path = _private_path(path, exists=False)
            fd = os.open(self.path, os.O_WRONLY|os.O_CREAT|os.O_EXCL|os.O_NOFOLLOW,0o600)
            os.close(fd)
            self.connection = sqlite3.connect(self.path)
            db = self.connection
            _settings(db, execution_config, readonly=False)
            db.executescript('''
                PRAGMA user_version=5;
                CREATE TABLE run_binding(singleton INTEGER PRIMARY KEY CHECK(singleton=1), canonical_bytes BLOB NOT NULL);
                CREATE TABLE enriched(ordinal INTEGER PRIMARY KEY,identity_id TEXT NOT NULL UNIQUE,allele_id TEXT,
                    admission_bytes BLOB NOT NULL,proofs_bytes BLOB NOT NULL,input_size INTEGER NOT NULL,
                    proof_size INTEGER NOT NULL,input_sha256 TEXT NOT NULL,proof_sha256 TEXT NOT NULL,
                    proof_count INTEGER NOT NULL,occurrences INTEGER NOT NULL,requests INTEGER NOT NULL,
                    gene_ids BLOB NOT NULL,screening_facts BLOB NOT NULL);
                CREATE INDEX enriched_allele ON enriched(allele_id,ordinal);
                CREATE INDEX enriched_counts ON enriched(ordinal,proof_count,occurrences,requests);
                CREATE TABLE genes(gene_id TEXT PRIMARY KEY);
                CREATE TABLE screening(ordinal INTEGER PRIMARY KEY REFERENCES enriched(ordinal),
                    canonical_bytes BLOB NOT NULL,sha256 TEXT NOT NULL,disposition TEXT NOT NULL,
                    accounting_bytes BLOB NOT NULL);
                CREATE INDEX screening_assessed ON screening(disposition,ordinal);
                CREATE TABLE source_occurrences(source_id TEXT PRIMARY KEY,ordinal INTEGER NOT NULL REFERENCES enriched(ordinal));
                CREATE TABLE eligibility(singleton INTEGER PRIMARY KEY CHECK(singleton=1),canonical_bytes BLOB NOT NULL);
                CREATE TABLE completion(singleton INTEGER PRIMARY KEY CHECK(singleton=1),canonical_bytes BLOB NOT NULL);
            ''')
            db.execute('INSERT INTO run_binding VALUES (1,?)',(binding,))
            self._commit()
        except (OSError, sqlite3.Error):
            self.close()
            raise AdmissionProjectionError('cannot create admission projection') from None

    @_projection_operation
    def add(self, row, decision):
        self._add(row,decision)
        if self.count % self.config.sqlite.cursor_batch_rows == 0:
            self._commit()

    @_projection_operation
    def add_window(self, rows, batch):
        self._require_window(rows,batch)
        for row,decision in zip(rows,batch.decisions,strict=True):
            self._add(row,decision)
        self._commit()
        self.windows += 1

    def _require_window(self, rows, batch):
        from genosapiens.wgs_analysis_eligibility import EligibilityBatch
        if (type(rows) is not tuple or not rows or type(batch) is not EligibilityBatch
                or len(rows) != len(batch.decisions)
                or len(rows) > min(_ADMISSION_WINDOW_ROWS,self.config.public_queries.batch_rows,
                                  self.config.sqlite.cursor_batch_rows)):
            raise AdmissionProjectionError('admission publication window is invalid')

    def add_screened_window(self, rows, batch, *, eligibility_writer, witnesses):
        """Authenticate one owned row for both stores within this synchronous call.

        No prepared row or validation authority escapes or survives this call.
        Standalone writer APIs still authenticate their own external inputs.
        """
        import time
        from genosapiens.wgs_analysis_eligibility import AnalysisEligibilityWriter, EligibilityBatch
        started=time.monotonic_ns();eligibility_ns=0
        def eligibility_write(function,*args,**kwargs):
            nonlocal eligibility_ns
            before=time.monotonic_ns()
            try:return function(*args,**kwargs)
            finally:
                elapsed=time.monotonic_ns()-before
                eligibility_ns+=elapsed
                eligibility_writer.metrics['write_nanoseconds']+=elapsed
        try:
            self._require_window(rows,batch)
            if (type(eligibility_writer) is not AnalysisEligibilityWriter
                    or self.finished or self.failed or eligibility_writer.finished
                    or getattr(eligibility_writer,'failed',False)
                    or eligibility_writer.rows!=self.count
                    or eligibility_writer.bindings[0]!=self.binding_sha256
                    or eligibility_writer.config!=self.config):
                raise AdmissionProjectionError('joint admission writer ownership differs')
            # Decision fields are strictly checked immutable scalar/tuple values.
            # Copy their owning objects before invoking either store.
            batch.__post_init__()
            batch=EligibilityBatch(tuple(replace(decision) for decision in batch.decisions))
            eligibility_write(eligibility_writer._begin_window,rows,batch,witnesses=witnesses)
            for external,decision in zip(rows,batch.decisions,strict=True):
                if type(external) is not EncodedAdmissionRow:
                    raise AdmissionProjectionError('joint admission row type is invalid')
                # Snapshot bytes/scalars without calling the constructor twice;
                # _add below performs the complete fresh blob/fact authentication.
                row=object.__new__(EncodedAdmissionRow)
                for name in EncodedAdmissionRow.__dataclass_fields__:
                    object.__setattr__(row,name,getattr(external,name))
                record,compact=self._add(row,decision)
                eligibility_write(eligibility_writer._add_authenticated,row,record,compact)
                del row,record,compact
            self._commit();self.windows+=1
            eligibility_write(eligibility_writer._commit_window)
        except BaseException:
            self.failed=True
            if type(eligibility_writer) is AnalysisEligibilityWriter:eligibility_writer.failed=True
            raise
        finally:
            self.metrics['write_nanoseconds']+=time.monotonic_ns()-started-eligibility_ns

    def _commit(self):
        import time
        started=time.monotonic_ns()
        try:self.connection.commit()
        finally:self.metrics['commit_nanoseconds']+=time.monotonic_ns()-started
        self.metrics['commits']+=1
        self.metrics['committed_inputs']=self.count
        self.metrics['committed_occurrences']=self.occurrences

    def _add(self, row, decision):
        if (self.finished or self.failed or type(row) is not EncodedAdmissionRow or row.binding_sha256 != self.binding_sha256
                or row.retained_ordinal != self.count+1 or self.count >= self.parents['identities']):
            raise AdmissionProjectionError('admission projection row order is invalid')
        content, _, facts = row._authenticated_parts()
        _require_fact_parents(facts, self.parents)
        record = AdmissionScreeningRecord._from_authenticated_facts(row, decision, facts)
        compact = record.as_bytes()
        try:
            self.connection.execute('INSERT INTO enriched VALUES ('+','.join('?' for _ in range(14))+')',
                (row.retained_ordinal,row.identity_id,row.allele_id,row.input_blob,
                 row.proof_blob,row.input_size,row.proof_size,
                 row.input_sha256,row.proof_sha256,row.proof_count,row.source_occurrences,
                 row.sequence_requests,canonical(row.gene_ids),row.screening_facts))
            self.connection.execute('INSERT INTO screening VALUES (?,?,?,?,?)',
                (row.retained_ordinal, compact, hashlib.sha256(compact).hexdigest(), decision.disposition,
                    AdmissionAccountingRow.from_record(record).as_bytes()))
            self.connection.executemany('INSERT INTO source_occurrences VALUES (?,?)',
                ((key,row.retained_ordinal) for key in record.facts.source_allele_ids))
            if decision.disposition == 'assess':
                self.connection.executemany('INSERT OR IGNORE INTO genes VALUES (?)',((g,) for g in row.gene_ids))
            self.count += 1
            self.proofs += row.proof_count
            self.occurrences += row.source_occurrences
            self.requests += row.sequence_requests
            self.inputs_digest.update(content)
            self.input_hashes.update((row.input_sha256+'\n').encode())
            self.proof_hashes.update((row.proof_sha256+'\n').encode())
            self.screening_digest.update(compact)
            self.input_size += row.input_size
            self.proof_size += row.proof_size
            return record,compact
        except sqlite3.Error:
            self.failed = True
            raise AdmissionProjectionError('admission projection write failed') from None

    @_projection_operation
    def seal(self, *, eligibility_manifest):
        if (self.finished or self.failed or self.count != self.parents['identities'] or self.occurrences != self.parents['occurrences']
                or self.requests != self.parents['sequence_requests']):
            raise AdmissionProjectionError('admission projection coverage is incomplete')
        if (type(eligibility_manifest) is not AnalysisEligibilityManifest
                or (eligibility_manifest.admission_binding_sha256, eligibility_manifest.eligibility_policy_sha256,
                    eligibility_manifest.input_rows, eligibility_manifest.source_occurrences,
                    eligibility_manifest.ordered_original_inputs_sha256, eligibility_manifest.ordered_original_proofs_sha256,
                    eligibility_manifest.ordered_screening_sha256) !=
                    (self.binding_sha256,self.parents['eligibility_policy_sha256'],self.count,self.occurrences,
                        self.input_hashes.hexdigest(),self.proof_hashes.hexdigest(),self.screening_digest.hexdigest())):
            raise AdmissionProjectionError('admission eligibility closure is incomplete')
        stored=hashlib.sha256(); stored_rows=0
        for ordinal,content,row_sha in self.connection.execute(
                'SELECT ordinal,canonical_bytes,sha256 FROM screening ORDER BY ordinal'):
            if ordinal != stored_rows+1 or hashlib.sha256(content).hexdigest() != row_sha:
                raise AdmissionProjectionError('admission compact rows changed before sealing')
            stored.update(content);stored_rows+=1
        if (stored_rows != self.count or stored.hexdigest() != self.screening_digest.hexdigest()
                or self.connection.execute('SELECT count(*) FROM source_occurrences').fetchone() != (self.occurrences,)):
            raise AdmissionProjectionError('admission compact coverage changed before sealing')
        digest = hashlib.sha256()
        count = 0
        for (gene,) in self.connection.execute('SELECT gene_id FROM genes ORDER BY gene_id'):
            digest.update(gene.encode('ascii')+b'\n')
            count += 1
        try:
            if self.connection.execute('PRAGMA foreign_key_check').fetchall():
                raise AdmissionProjectionError('admission source coverage is invalid')
            completion = AdmissionProjectionCompletion(self.binding_sha256,self.count,self.proofs,self.occurrences,
                self.requests,self.input_size,self.proof_size,count,digest.hexdigest(),self.inputs_digest.hexdigest(),
                self.input_hashes.hexdigest(),self.proof_hashes.hexdigest(),self.screening_digest.hexdigest(),eligibility_manifest.sha256)
            self.connection.execute('INSERT INTO eligibility VALUES (1,?)',(eligibility_manifest.as_bytes(),))
            self.connection.execute('INSERT INTO completion VALUES (1,?)',(completion.as_bytes(),))
            self._commit()
            self.finished = True
            return completion
        except sqlite3.Error:
            self.failed = True
            raise AdmissionProjectionError('admission projection seal failed') from None

    def close(self):
        if self.connection is not None:
            self.connection.close()
            self.connection = None


class AuthenticatedAdmissionOwner:
    """Reuse anchored producer completion and validate evidence when consumed.

    Independent processes check the same completion and bounded metadata. They
    trust producer validation of unread payloads; borrowed views share the file
    lifetime, requested-record checks and bounded caches, never a trust flag.
    """
    def __init__(self, path, *, expected_binding, expected_completion_sha256, execution_config, progress=None):
        if progress is not None and not callable(progress):
            raise AdmissionProjectionError('admission progress callback is invalid')
        self._progress=progress
        self._progress_started=time.monotonic()
        self._progress_last=0.0
        self._progress_phase='opening'
        self.parents = _binding(expected_binding)
        _require_execution_config(execution_config,self.parents)
        self.binding = expected_binding
        self.binding_sha256 = hashlib.sha256(expected_binding).hexdigest()
        self.config = execution_config
        self._connections = []
        self._signature_fd = None
        self._file_cache = None
        self._connection_local = threading.local()
        self._lock = threading.RLock()
        self._closed = self._failed = self._authenticated = False
        self._cache = OrderedDict()
        self.cache_bytes = 0
        self.accounting_sha256 = self.assessed_ordinals_sha256 = None
        self.metrics = dict(compact_rows=0,rich_inputs=0,rich_proofs=0,
            authentication_scans=0,gene_authentications=0,connections=0,borrowed_views=0,scalar_rows=0,
            compact_cache_hits=0,input_cache_hits=0,scalar_cache_hits=0,max_cache_bytes=0,max_cache_entries=0)
        self.metrics['authentication_nanoseconds']=0
        self.metrics.update(completion_reuses=0,metadata_check_nanoseconds=0,
            complete_domain_traversals=0)
        self.metrics.update(compact_read_rows=0,compact_read_batches=0,
            max_compact_read_batch_bytes=0,compact_read_stopped=False)
        self.metrics.update((f'{kind}_{purpose}',0) for kind in ('rich_inputs','rich_proofs') for purpose in _READ_PURPOSES)
        try:
            self.path = _private_path(path, exists=True)
            self._signature_fd = os.open(self.path,os.O_RDONLY|os.O_NOFOLLOW)
            self._file_cache = ReadFileCache(self._signature_fd)
            self._file_signature = self._signature()
            db = self.connection
            db.set_progress_handler(lambda:(self._notify_progress(),0)[1],100000)
            self._notify_progress(force=True)
            if (db.execute('PRAGMA user_version').fetchone() != (5,)
                    or db.execute('SELECT canonical_bytes FROM run_binding WHERE singleton=1').fetchone() != (expected_binding,)):
                raise ValueError
            completion = db.execute('SELECT canonical_bytes FROM completion WHERE singleton=1').fetchone()
            if completion is None:
                raise ValueError
            self.completion = load_admission_projection_completion(completion[0],expected_sha256=expected_completion_sha256)
            if self.completion.binding_sha256 != self.binding_sha256:
                raise ValueError
            complete=self.completion
            self.row_count,self.proof_count,self.occurrence_count,self.request_count = (
                complete.rows,complete.proofs,complete.occurrences,complete.requests)
            self.gene_count,self.gene_sha256,self.input_rows_sha256 = (
                complete.genes,complete.gene_sha256,complete.input_rows_sha256)
            eligibility = db.execute('SELECT canonical_bytes FROM eligibility WHERE singleton=1').fetchone()
            if eligibility is None:
                raise ValueError
            self.eligibility_manifest = load_analysis_eligibility_manifest(eligibility[0],complete.eligibility_manifest_sha256)
            if (self.eligibility_manifest.admission_binding_sha256 != self.binding_sha256
                    or self.eligibility_manifest.ordered_screening_sha256 != complete.screening_sha256
                    or self.eligibility_manifest.input_rows != complete.rows
                    or self.eligibility_manifest.source_occurrences != complete.occurrences):
                raise ValueError
            if (self.row_count != self.parents['identities'] or self.occurrence_count != self.parents['occurrences']
                    or self.request_count != self.parents['sequence_requests']):
                raise ValueError
            if (self.eligibility_manifest.ordered_original_inputs_sha256 != complete.ordered_input_hashes_sha256
                    or self.eligibility_manifest.ordered_original_proofs_sha256 != complete.ordered_proof_hashes_sha256
                    or self.eligibility_manifest.eligibility_policy_sha256 != self.parents['eligibility_policy_sha256']
                    or self.eligibility_manifest.source_manifest_sha256 != self.parents['source_manifest_sha256']):
                raise ValueError
            tables={r[0] for r in db.execute("SELECT name FROM sqlite_schema WHERE type='table'")}
            indexes={r[0] for r in db.execute("SELECT name FROM sqlite_schema WHERE type='index'")}
            if not {'run_binding','completion','eligibility','enriched','screening','genes','source_occurrences'} <= tables \
                    or not {'enriched_counts','screening_assessed'} <= indexes:
                raise ValueError
            self._assert_stable()
            self._authenticated = True
            db.set_progress_handler(lambda:(self._file_cache.release(),0)[1],100000)
            self.metrics['completion_reuses']=1
            self.metrics['metadata_check_nanoseconds']=int((time.monotonic()-self._progress_started)*1e9)
            self._progress_phase='completed_input_reused'
            self._notify_progress(force=True)
        except (OSError,sqlite3.Error,ValueError,TypeError):
            try:self.close()
            except (OSError,sqlite3.Error):pass
            raise AdmissionProjectionError('admission projection is missing or incomplete') from None
        except BaseException:
            try:self.close()
            except (OSError,sqlite3.Error):pass
            raise

    def _notify_progress(self,*,force=False):
        self._file_cache.release(force=force)
        self.metrics.update(file_cache_releases=self._file_cache.releases,
            file_cache_release_nanoseconds=self._file_cache.release_nanoseconds)
        if self._progress is None:return
        now=time.monotonic()
        if not force and now-self._progress_last<1:return
        self._progress_last=now
        self._progress(dict(phase=self._progress_phase,elapsed_seconds=now-self._progress_started,
            metrics=dict(self.metrics),cache_bytes=self.cache_bytes,cache_entries=len(self._cache)))

    @property
    def connection(self):
        self._assert_stable()
        with self._lock:
            connection = getattr(self._connection_local, 'connection', None)
            if connection is None:
                # Sessions use persistent stage threads. Bound even accidental
                # thread churn; each thread gets its own read-only connection.
                if len(self._connections) >= 12:
                    raise AdmissionProjectionError('admission connection owner exceeds capacity')
                connection = sqlite3.connect(self.path.as_uri()+'?mode=ro&immutable=1',uri=True,
                    check_same_thread=False)
                try:
                    _settings(connection,self.config,readonly=True)
                    connection.set_progress_handler(lambda:(self._file_cache.release(),0)[1],100000)
                except BaseException:
                    connection.close()
                    raise
                self._connections.append(connection)
                self._connection_local.connection = connection
                self.metrics['connections'] += 1
            return connection

    def reader(self):
        self._assert_stable()
        if not self._authenticated:
            raise AdmissionProjectionError('admission authentication has not completed')
        return AdmissionProjectionReader._borrow(self)

    @property
    def cache_entries(self):
        return len(self._cache)

    def _cached(self,key):
        item=self._cache.get(key)
        if item is None:return None
        self._cache.move_to_end(key)
        self.metrics[{'compact':'compact_cache_hits','input':'input_cache_hits','scalar':'scalar_cache_hits'}[key[0]]]+=1
        return item[0]

    def _remember(self,key,value,cost):
        if cost>_READER_CACHE_BYTES:return
        old=self._cache.pop(key,None)
        if old is not None:self.cache_bytes-=old[1]
        while self._cache and (self.cache_bytes+cost>_READER_CACHE_BYTES or len(self._cache)>=_READER_CACHE_ENTRIES):
            _,old=self._cache.popitem(last=False);self.cache_bytes-=old[1]
        self._cache[key]=(value,cost);self.cache_bytes+=cost
        self.metrics['max_cache_bytes']=max(self.metrics['max_cache_bytes'],self.cache_bytes)
        self.metrics['max_cache_entries']=max(self.metrics['max_cache_entries'],len(self._cache))

    def _count_read(self,purpose,*,proofs):
        if purpose not in _READ_PURPOSES:raise AdmissionProjectionError('invalid admission read purpose')
        self.metrics['rich_inputs']+=1;self.metrics[f'rich_inputs_{purpose}']+=1
        if proofs:
            self.metrics['rich_proofs']+=1;self.metrics[f'rich_proofs_{purpose}']+=1

    def _row(self, values):
        return _reader_row(values, self.binding_sha256)

    def _signature(self):
        value=self.path.lstat()
        if not stat.S_ISREG(value.st_mode) or self.path.is_symlink():
            raise AdmissionProjectionError('admission projection path changed')
        # Fast consecutive commits can share a filesystem timestamp tick. The
        # pinned SQLite header also binds its change counter/schema/version-valid
        # fields, without asking an immutable SQLite pager to refresh its cache.
        header=os.pread(self._signature_fd,100,0)
        if len(header)!=100: raise AdmissionProjectionError('admission projection header changed')
        return (value.st_dev,value.st_ino,value.st_size,value.st_mtime_ns,value.st_ctime_ns,
            value.st_mode,value.st_uid,value.st_nlink,header)

    def _assert_stable(self):
        if self._closed or self._failed:
            raise AdmissionProjectionError('admission owner is closed or failed')
        try:
            if self._signature()!=self._file_signature:
                raise AdmissionProjectionError('admission projection changed during read')
            self._file_cache.release()
        except (OSError,AdmissionProjectionError):
            self._failed=True
            raise AdmissionProjectionError('admission projection changed during read') from None

    def _ensure_screening(self):
        """Require live completion authority; never initiate a payload traversal."""
        self._assert_stable()
        if not self._authenticated:
            raise AdmissionProjectionError('admission completion is not available')

    def _validate_compact(self, record, ordinal, disposition):
        facts,decision=record.facts,record.decision
        if (not _count(ordinal,1) or ordinal>self.row_count
                or facts.retained_ordinal!=ordinal or decision.admission_binding_sha256!=self.binding_sha256
                or decision.disposition!=disposition
                or facts.eligibility_policy_sha256!=self.parents['eligibility_policy_sha256']
                or (facts.source_manifest_sha256,facts.working_manifest_sha256,facts.source_routing_sha256)
                    !=_row_parents(self.parents)
                or facts.merged_intake_sha256!=self.parents['read_intake_sha256']
                or decision.resource_bindings!=self.eligibility_manifest.resource_bindings):
            raise AdmissionProjectionError('requested compact parents or identity changed')

    def _compact(self,ordinal,*,with_size=False):
        with self._lock:
            self._ensure_screening()
            cached=self._cached(('compact',ordinal))
            if cached is not None:
                return (cached,self._cache['compact',ordinal][1]) if with_size else cached
            saved=self.connection.execute('SELECT canonical_bytes,sha256,disposition FROM screening WHERE ordinal=?',(ordinal,)).fetchone()
            if saved is None or hashlib.sha256(saved[0]).hexdigest()!=saved[1]:
                raise AdmissionProjectionError('requested compact input is missing or changed')
            record=load_admission_screening_record(saved[0])
            self._validate_compact(record,ordinal,saved[2])
            self.metrics['compact_rows']+=1
            self._remember(('compact',ordinal),record,len(saved[0]))
            return (record,len(saved[0])) if with_size else record

    def iter_screening_records(self, *, original_ordinals=None):
        """Explicit compact reads under an already authenticated live owner."""
        self._ensure_screening()
        if original_ordinals is not None:
            if (type(original_ordinals) is not tuple or not original_ordinals
                    or len(original_ordinals) > 50000
                    or any(type(v) is not int or not 1 <= v <= self.row_count for v in original_ordinals)
                    or original_ordinals != tuple(sorted(set(original_ordinals)))):
                raise AdmissionProjectionError('finite original compact domain is invalid')
            for ordinal in original_ordinals:
                yield self._compact(ordinal)
            self._assert_stable()
            return
        yield from self._iter_complete_screening()

    def _iter_complete_screening(self):
        """Validate the full domain while it feeds actual full-domain processing."""
        self._ensure_screening()
        self.metrics['complete_domain_traversals']+=1
        self._progress_phase='processing_complete_domain'
        self._notify_progress(force=True)
        read_start=self.metrics['compact_read_rows']
        digest=hashlib.sha256(); inputs=hashlib.sha256(); proofs=hashlib.sha256()
        decisions=hashlib.sha256()
        accounting=hashlib.sha256()
        assessed_ordinals=hashlib.sha256()
        count=occurrences=proof_count=requests=input_size=proof_size=0
        totals={key:0 for key in ('assess','exclude_common','defer_refcall')}
        occurrence_totals=dict(totals)
        try:
            with _compact_read_ahead(self.connection,self._notify_progress,self.metrics,
                    cache_release=self._file_cache.release) as cursor:
                for ordinal,content,expected,disposition,scalar in cursor:
                    if ordinal != count+1 or hashlib.sha256(content).hexdigest() != expected:
                        raise AdmissionProjectionError('screening compact coverage changed')
                    record=load_admission_screening_record(content)
                    if AdmissionAccountingRow.from_record(record).as_bytes()!=scalar:
                        raise AdmissionProjectionError('scalar accounting does not own compact record')
                    accounting.update(scalar)
                    facts,decision=record.facts,record.decision
                    self._validate_compact(record,ordinal,disposition)
                    digest.update(content); count+=1; occurrences+=len(facts.source_allele_ids)
                    proof_count+=record.proof_count;requests+=int(facts.transport=='sequence_vep')
                    input_size+=record.input_size;proof_size+=record.proof_size
                    inputs.update((decision.input_sha256+'\n').encode());proofs.update((decision.proof_sha256+'\n').encode())
                    decisions.update(decision.as_bytes())
                    totals[disposition]+=1;occurrence_totals[disposition]+=len(facts.source_allele_ids)
                    if disposition=='assess':assessed_ordinals.update((str(ordinal)+'\n').encode())
                    self.metrics['compact_rows']+=1
                    self._notify_progress()
                    self._remember(('compact',ordinal),record,len(content))
                    yield record
            complete=self.completion
            if (self.metrics['compact_read_rows']-read_start!=count or not self.metrics['compact_read_stopped']
                    or (count,occurrences,proof_count,requests,input_size,proof_size,digest.hexdigest(),inputs.hexdigest(),proofs.hexdigest()) !=
                    (complete.rows,complete.occurrences,complete.proofs,complete.requests,complete.input_bytes,complete.proof_bytes,
                        complete.screening_sha256,complete.ordered_input_hashes_sha256,complete.ordered_proof_hashes_sha256)
                    or tuple(sorted(totals.items())) != self.eligibility_manifest.input_dispositions
                    or tuple(sorted(occurrence_totals.items())) != self.eligibility_manifest.occurrence_dispositions
                    or decisions.hexdigest() != self.eligibility_manifest.ordered_decisions_sha256):
                raise AdmissionProjectionError('screening compact closure changed')
            self._assert_stable()
            self.accounting_sha256=accounting.hexdigest()
            self.assessed_ordinals_sha256=assessed_ordinals.hexdigest()
        except (sqlite3.Error,TypeError,ValueError):
            raise AdmissionProjectionError('cannot authenticate compact screening') from None

    def iter_accounting(self):
        """Consume and close original accounting in one pass, without rich reads."""
        with closing(self.iter_screening_records()) as records:
            for record in records:
                self.metrics['scalar_rows']+=1
                yield AdmissionAccountingRow.from_record(record)

    def _scalar(self,ordinal):
        with self._lock:
            self._ensure_screening()
            row=self._cached(('scalar',ordinal))
            if row is not None:return row
            saved=self.connection.execute('SELECT accounting_bytes FROM screening WHERE ordinal=?',(ordinal,)).fetchone()
            if saved is None:raise AdmissionProjectionError('requested scalar admission is missing')
            row=_load_admission_accounting(saved[0])
            if row!=AdmissionAccountingRow.from_record(self._compact(ordinal)):
                raise AdmissionProjectionError('requested scalar admission changed')
            self.metrics['scalar_rows']+=1
            self._remember(('scalar',ordinal),row,len(saved[0]))
            return row

    def iter_selected_rows(self, profile_reader, *, _decode_rows=None, purpose='discovery'):
        """Strict rich rows in the independently sealed advancing profile domain."""
        from genosapiens.wgs_profile_projection import ProfileProjectionReader
        if (type(profile_reader) is not ProfileProjectionReader
                or profile_reader.identity['parent']['binding_sha256'] != self.binding_sha256
                or profile_reader.identity['parent']['completion_sha256'] != self.completion.sha256):
            raise AdmissionProjectionError('selected profile admission owner differs')
        profile_reader.require_stable()
        yield from self.iter_assessed_rows(_decode_rows=_decode_rows,purpose=purpose,_profile=profile_reader)
        profile_reader.require_stable()

    def iter_assessed_rows(self, *, _decode_rows=None, purpose='discovery', _profile=None):
        """Read only advancing evidence and close its domain during consumption."""
        self._ensure_screening()
        if purpose not in _READ_PURPOSES:raise AdmissionProjectionError('invalid admission read purpose')
        pending=deque()
        decoded_count=occurrences=previous=0
        ordinals=hashlib.sha256()
        expected_count=(dict(self.eligibility_manifest.input_dispositions)['assess'] if _profile is None
                        else _profile.manifest['counts']['advancing_rows'])
        def references():
            if _profile is None:
                with closing(self.connection.execute("SELECT ordinal FROM screening INDEXED BY screening_assessed WHERE disposition='assess' ORDER BY ordinal")) as records:
                    for ordinal, in records:yield ordinal,None
            else:
                for selected in _profile.iter_rows():
                    if selected.decision['disposition'] in ('selected','unresolved'):
                        yield selected.original_ref[1],selected.original_ref
        def raw_rows():
            with closing(references()) as records:
                for ordinal,ref in records:
                    compact,compact_size=self._compact(ordinal,with_size=True)
                    if compact.decision.disposition!='assess' or ref is not None and compact.ref!=ref:
                        raise AdmissionProjectionError('assessed index changed')
                    values=self.connection.execute('SELECT '+_COLUMNS+' FROM enriched WHERE ordinal=?',
                        (compact.facts.retained_ordinal,)).fetchone()
                    if values is None:
                        raise AdmissionProjectionError('assessed rich input is missing')
                    self._count_read(purpose,proofs=True)
                    pending.append((compact,compact_size))
                    if len(pending)>max(self.config.public_queries.batch_rows,self.config.sqlite.cursor_batch_rows)+1:
                        raise AdmissionProjectionError('assessed decode window exceeds capacity')
                    yield values
        def serial(source):
            for values in source:
                row=self._row(values)
                value,proofs=row.decode(binding=self.binding)
                yield AdmissionProjectionRecord(row,value,proofs,_VALIDATED_ROW_SEAL)
        with closing(raw_rows()) as source, closing(serial(source) if _decode_rows is None else
                _decode_rows(source,self.binding)) as decoded:
            for record in decoded:
                if (not pending or type(record) is not AdmissionProjectionRecord
                        or record.seal is not _VALIDATED_ROW_SEAL):
                    raise AdmissionProjectionError('assessed input lacks its live decode owner')
                compact,compact_size=pending.popleft()
                ordinal=compact.facts.retained_ordinal
                if ordinal<=previous:
                    raise AdmissionProjectionError('assessed input order changed')
                previous=ordinal;occurrences+=len(compact.facts.source_allele_ids)
                ordinals.update((str(ordinal)+'\n').encode())
                decoded_count+=1
                yield self._attach_compact_record(record,compact,compact_size)
            expected_occurrences=(dict(self.eligibility_manifest.occurrence_dispositions)['assess'] if _profile is None
                else _profile.manifest['counts']['advancing_occurrences'])
            if pending or decoded_count!=expected_count or occurrences!=expected_occurrences:
                raise AdmissionProjectionError('assessed input coverage is incomplete')
            if _profile is None:
                self.assessed_ordinals_sha256=ordinals.hexdigest()
            elif ordinals.hexdigest()!=_profile.manifest['counts']['advancing_ordinals_sha256']:
                raise AdmissionProjectionError('selected input coverage changed')
        self._assert_stable()

    def _attach_compact_record(self, record, compact, compact_size):
        """Link this owner's validated compact row to its live decoded reply.

        Both sides already passed their strict decoder. Comparing the exact
        fact bytes and every scalar commitment retains the full linkage check
        without parsing and reconstructing those same facts in the coordinator.
        Callers hold compact rows from this owner in their bounded read window.
        """
        self._assert_stable()
        if (type(record) is not AdmissionProjectionRecord or record.seal is not _VALIDATED_ROW_SEAL
                or type(record.encoded) is not _ReaderRow or type(compact) is not AdmissionScreeningRecord):
            raise AdmissionProjectionError('assessed input lacks its live decode owner')
        row=record.encoded
        if (row.ref!=compact.ref or row.binding_sha256!=self.binding_sha256
                or row.proof_sha256!=compact.decision.proof_sha256
                or (row.proof_size,row.proof_count,row.source_occurrences,row.sequence_requests,row.gene_ids) !=
                    (compact.proof_size,compact.proof_count,len(compact.facts.source_allele_ids),
                        int(compact.facts.transport=='sequence_vep'),compact.facts.gene_ids)
                or row.screening_facts!=compact.facts.as_bytes()):
            raise AdmissionProjectionError('assessed rich input does not own screening record')
        return replace(record,compact=compact,compact_bytes=compact_size)

    def iter_gene_ids(self):
        self._ensure_screening()
        digest = hashlib.sha256()
        count = 0
        try:
            with closing(self.connection.execute('SELECT gene_id FROM genes ORDER BY gene_id')) as cursor:
                for (gene,) in cursor:
                    if type(gene) is not str or not gene or not gene.isascii():
                        raise ValueError
                    digest.update(gene.encode('ascii')+b'\n')
                    count += 1
                    yield gene
            if count != self.gene_count or digest.hexdigest() != self.gene_sha256:
                raise ValueError
            self._assert_stable()
        except (sqlite3.Error,ValueError):
            raise AdmissionProjectionError('admission gene projection is incomplete') from None

    def _bounded_refs(self, refs):
        self._ensure_screening()
        bounded = []
        seen = set()
        total = 0
        for ref in refs:
            if (type(ref) is not tuple or len(ref) != 6 or ref[0] != self.binding_sha256
                    or not _count(ref[1],1) or type(ref[2]) is not str
                    or ref[3] is not None and type(ref[3]) is not str
                    or not _sha(ref[4]) or not _count(ref[5],1) or ref[5]>MAX_ROW_BYTES):
                raise AdmissionProjectionError('admission input reference is invalid')
            total += ref[5]
            if (len(bounded) >= self.config.contexts.max_member_identities or ref[1] in seen
                    or total > self.config.contexts.max_context_bytes):
                raise AdmissionProjectionError('admission input references exceed capacity or repeat')
            seen.add(ref[1])
            bounded.append(ref)
        return tuple(bounded)

    def get_inputs(self, refs, *, purpose='context'):
        with self._lock:
            if purpose not in _READ_PURPOSES:raise AdmissionProjectionError('invalid admission read purpose')
            return self._get_inputs(self._bounded_refs(refs),purpose)

    def _get_inputs(self, bounded, purpose):
        result = []
        try:
            for ref in bounded:
                cached=self._cached(('input',ref))
                if cached is not None:
                    result.append(cached)
                    continue
                compact=self._compact(ref[1])
                if compact.ref!=ref:
                    raise AdmissionProjectionError('requested compact reference changed')
                values = self.connection.execute('SELECT '+_INPUT_COLUMNS+' FROM enriched WHERE ordinal=?',(ref[1],)).fetchone()
                if values is None:
                    raise AdmissionProjectionError('admission input reference is missing')
                self._count_read(purpose,proofs=False)
                ordinal,identity,allele,blob,size,digest,occurrences,requests,genes=values
                facts=compact.facts
                if ((self.binding_sha256,ordinal,identity,allele,digest,size)!=ref
                        or occurrences!=len(facts.source_allele_ids) or requests!=int(facts.transport=='sequence_vep')
                        or genes!=canonical(facts.gene_ids)):
                    raise AdmissionProjectionError('requested input metadata changed')
                content=_unpack_admission_blob(blob)
                if len(content)!=size or hashlib.sha256(content).hexdigest()!=digest:
                    raise AdmissionProjectionError('requested input bytes changed')
                value=load_ledger_admission_input(content)
                actual_genes=tuple(sorted({t.ensembl_gene_id for t in
                    (() if value.allele is None else value.allele.transcript_consequences or ())
                    if t.ensembl_gene_id is not None}))
                if ((value.retained_ordinal,value.identity_id,None if value.allele is None else value.allele.allele_id)
                            !=(ordinal,identity,allele)
                        or value.source_allele_ids!=facts.source_allele_ids or value.transport!=facts.transport
                        or actual_genes!=facts.gene_ids
                        or (value.parent_source_manifest_sha256,value.parent_working_set_manifest_sha256,
                            value.native_source_contract_sha256)!=_row_parents(self.parents)):
                    raise AdmissionProjectionError('requested input does not own compact record')
                result.append(value)
                self._remember(('input',ref),value,len(content))
            self._assert_stable()
            return tuple(result)
        except sqlite3.Error:
            raise AdmissionProjectionError('admission input reference lookup failed') from None

    def get_records(self, refs, *, purpose='audit'):
        """Explicit full evidence reads always authenticate requested proof blobs."""
        with self._lock:
            if purpose not in _READ_PURPOSES:raise AdmissionProjectionError('invalid admission read purpose')
            result=[];size=0
            for ref in self._bounded_refs(refs):
                compact=self._compact(ref[1])
                if compact.ref!=ref:raise AdmissionProjectionError('requested compact reference changed')
                size+=compact.input_size+compact.proof_size
                if size>self.config.contexts.max_context_bytes:
                    raise AdmissionProjectionError('requested full evidence exceeds capacity')
                values=self.connection.execute('SELECT '+_COLUMNS+' FROM enriched WHERE ordinal=?',(ref[1],)).fetchone()
                if values is None:raise AdmissionProjectionError('requested full evidence is missing')
                self._count_read(purpose,proofs=True)
                row=self._row(values)
                if AdmissionScreeningRecord.from_encoded(row,compact.decision)!=compact:
                    raise AdmissionProjectionError('requested full evidence does not own compact record')
                value,proofs=row.decode(binding=self.binding)
                result.append(AdmissionProjectionRecord(row,value,proofs,_VALIDATED_ROW_SEAL))
            self._assert_stable()
            return tuple(result)

    def refs_for_inputs(self, inputs, *, expected_rows_sha256):
        """Reconcile one admitted batch against the sealed scalar projection.

        Candidate admission already bound its original row hashes. This lookup
        carries those hashes without encoding the cloned candidate inputs again.
        Context members are subsequently decoded from the authoritative store.
        """
        if (type(inputs) is not tuple or not 0 < len(inputs) <= self.config.public_queries.batch_rows
                or not _sha(expected_rows_sha256)):
            raise AdmissionProjectionError('admission reference batch is invalid')
        refs = []
        source,working,routing = _row_parents(self.parents)
        try:
            ordinals=tuple(value.retained_ordinal for value in inputs)
            if ordinals != tuple(sorted(set(ordinals))):
                raise AdmissionProjectionError('admission reference batch repeats or reorders inputs')
            self._ensure_screening()
            for value in inputs:
                scalar=self._scalar(value.retained_ordinal)
                if (value.identity_id!=scalar.ref[2] or value.source_allele_ids!=scalar.source_allele_ids
                        or (None if value.allele is None else value.allele.allele_id)!=scalar.ref[3]
                        or (value.parent_source_manifest_sha256,value.parent_working_set_manifest_sha256,
                            value.native_source_contract_sha256)!=(source,working,routing)
                        or scalar.disposition!='assess'):
                    raise AdmissionProjectionError('admission reference batch changed')
                refs.append(scalar.ref)
            if hashlib.sha256(canonical([ref[4] for ref in refs])).hexdigest() != expected_rows_sha256:
                raise AdmissionProjectionError('admission reference batch hash changed')
            return tuple(refs)
        except (sqlite3.Error,ValueError):
            raise AdmissionProjectionError('admission reference batch lookup failed') from None

    def refs_for_subjects(self, subjects, *, expected_rows_sha256):
        """Reconcile discovery subjects without allocating rich input objects."""
        from genosapiens.wgs_association_projection import ContextSubjectProjection
        if (type(subjects) is not tuple or not 0<len(subjects)<=self.config.public_queries.batch_rows
                or not _sha(expected_rows_sha256)):
            raise AdmissionProjectionError('context subject reference batch is invalid')
        self._ensure_screening()
        refs=[];previous=0
        for value in subjects:
            if (type(value) is not ContextSubjectProjection or value.retained_ordinal<=previous
                    or value.admission_completion_sha256!=self.completion.sha256):
                raise AdmissionProjectionError('context subject admission owner differs')
            value.require_compact(self._compact(value.retained_ordinal))
            scalar=self._scalar(value.retained_ordinal)
            if (scalar.disposition!='assess' or scalar.ref!=value.admission_ref
                    or scalar.source_allele_ids!=value.source_allele_ids):
                raise AdmissionProjectionError('context subject scalar accounting differs')
            refs.append(scalar.ref);previous=value.retained_ordinal
        if hashlib.sha256(canonical([ref[4] for ref in refs])).hexdigest()!=expected_rows_sha256:
            raise AdmissionProjectionError('context subject original byte binding differs')
        self._assert_stable()
        return tuple(refs)

    def iter_allele_inputs(self, allele_id):
        self._ensure_screening()
        try:
            with closing(self.connection.execute('SELECT ordinal FROM enriched WHERE allele_id=? ORDER BY ordinal',(allele_id,))) as cursor:
                for ordinal, in cursor:
                    scalar=self._scalar(ordinal)
                    value,=self.get_inputs((scalar.ref,),purpose='selected')
                    if value.allele is None or value.allele.allele_id!=allele_id:
                        raise AdmissionProjectionError('selected admission lookup is invalid')
                    yield value
        except sqlite3.Error:
            raise AdmissionProjectionError('selected admission lookup failed') from None

    def close(self):
        with self._lock:
            self._closed=True
            self._authenticated=False
            self._cache.clear()
            self.cache_bytes = 0
            for connection in self._connections:
                connection.close()
            self._connections.clear()
            if self._signature_fd is not None:
                try:
                    if self._file_cache is not None:self._file_cache.release(force=True)
                finally:
                    os.close(self._signature_fd)
                    self._signature_fd=None


class AdmissionProjectionReader:
    """Borrowed view, or an ordinary independent read boundary owning its session."""
    def __init__(self, path, *, expected_binding, expected_completion_sha256, execution_config):
        self._owner=AuthenticatedAdmissionOwner(path,expected_binding=expected_binding,
            expected_completion_sha256=expected_completion_sha256,execution_config=execution_config)
        self._owns_owner=True
        self._closed=False

    @classmethod
    def _borrow(cls, owner):
        if type(owner) is not AuthenticatedAdmissionOwner or not owner._authenticated:
            raise AdmissionProjectionError('invalid admission reader owner')
        owner._assert_stable()
        result=object.__new__(cls)
        result._owner=owner;result._owns_owner=False;result._closed=False
        owner.metrics['borrowed_views']+=1
        return result

    def __getattr__(self, name):
        owner=object.__getattribute__(self,'_owner')
        if name=='connection' and object.__getattribute__(self,'_closed'):
            raise AdmissionProjectionError('admission reader view is closed')
        value=getattr(owner,name)
        if isinstance(value,MethodType):
            def borrowed(*args, **kwargs):
                if self._closed:
                    raise AdmissionProjectionError('admission reader view is closed')
                self._owner._assert_stable()
                return value(*args,**kwargs)
            return borrowed
        return value

    def close(self):
        if not self._closed:
            self._closed=True
            if self._owns_owner:self._owner.close()


# Private IPC carries canonical bytes and scalar metadata only. A full reservation
# remains charged until the ordered consumer returns, including reorder holds.
_ADMISSION_JOB_BYTES = 16 * 1024 * 1024
_ADMISSION_RESULT_BYTES = 64 * 1024 * 1024
# Deferred jobs additionally retain complete saved annotation until screening.
# Charge that reservation only to the jobs that carry it.
_ADMISSION_JOB_ROWS = min(16, (_ADMISSION_RESULT_BYTES - 1024*1024) // (3*(MAX_ROW_BYTES+65536)))
_DEFERRED_ADMISSION_JOB_ROWS = min(16, (_ADMISSION_RESULT_BYTES - 1024*1024) // (4*(MAX_ROW_BYTES+65536)))
_ADMISSION_MAX_JOBS = 20
_ADMISSION_WINDOW_ROWS = 64
_ADMISSION_WINDOW_BYTES = 64 * 1024 * 1024
_ADMISSION_MAX_WINDOWS = 2


@dataclass(frozen=True, slots=True, repr=False)
class AdmissionProjectionJob:
    worker: int
    ordinal: int
    binding_sha256: str
    setup_sha256: str
    # (retained ordinal, identity, authenticated saved annotation, optional
    # canonical input from the existing complete-locus/read-intake adapter).
    rows: tuple


@dataclass(frozen=True, slots=True, repr=False)
class PendingNativeAdmission:
    retained_ordinal: int
    identity_id: str
    input_bytes: bytes
    input_sha256: str
    binding_sha256: str
    source_occurrences: int

    def decode(self, binding):
        if (type(self.input_bytes) is not bytes or len(self.input_bytes) > MAX_ROW_BYTES
                or hashlib.sha256(self.input_bytes).hexdigest() != self.input_sha256
                or self.binding_sha256 != hashlib.sha256(binding).hexdigest()):
            raise AdmissionProjectionError('native admission transport is invalid')
        value = load_ledger_admission_input(self.input_bytes)
        if (value.retained_ordinal != self.retained_ordinal or value.identity_id != self.identity_id
                or value.transport != 'native_vep' or len(value.source_allele_ids)!=self.source_occurrences):
            raise AdmissionProjectionError('native admission transport changed')
        return value


@dataclass(frozen=True, slots=True, repr=False)
class ValidatedSequenceAdmission:
    """Encoded input plus original annotation; no broad association proof graph."""
    row: EncodedAdmissionRow
    annotation: object

    @property
    def retained_ordinal(self): return self.row.retained_ordinal

    @property
    def identity_id(self): return self.row.identity_id

    @property
    def binding_sha256(self): return self.row.binding_sha256

    @property
    def source_occurrences(self): return self.row.source_occurrences

    def authenticated_bytes(self):
        from genosapiens.wgs_annotation import AnnotationStageEncodedRecord
        if (type(self.row) is not EncodedAdmissionRow
                or type(self.annotation) is not AnnotationStageEncodedRecord
                or self.row.proof_count != 0
                or self.row.proof_sha256 != hashlib.sha256(_owned_admission_proofs_bytes(())).hexdigest()):
            raise AdmissionProjectionError('validated sequence intake is invalid')
        self.row.authenticated_bytes()
        self.annotation.__post_init__()
        self.annotation.encoded.__post_init__()

    def decode(self, binding):
        self.authenticated_bytes()
        value=self.row.decode_input(binding=binding)
        verified=self.annotation.decode()
        parents=_binding(binding)
        self._require_decoded(value, verified, parents)
        return value,verified

    def _require_decoded(self, value, verified, parents):
        if (value.transport != 'sequence_vep' or value.annotation_status != 'annotated'
                or (verified.request.retained_ordinal,verified.request.identity_id,verified.request.source_allele_ids) !=
                   (value.retained_ordinal,value.identity_id,value.source_allele_ids)
                or verified.request.parent_working_set_manifest_sha256 != parents['working_manifest_sha256']):
            raise AdmissionProjectionError('validated sequence intake identity differs')
        transcripts=verified.parsed.transcripts[(verified.request.request_id,1)]
        if transcripts and (value.allele.transcript_consequences or ())[-len(transcripts):] != transcripts:
            raise AdmissionProjectionError('validated sequence intake transcript content differs')


class _ProfileRecordContext:
    """Immutable admission parents for one intake owner, never a persisted flag."""
    def __init__(self, binding, screening):
        self.binding = binding
        self.parents = MappingProxyType(_binding(binding))
        _require_screening_context(screening, self.parents)
        self.binding_sha256 = hashlib.sha256(binding).hexdigest()
        self.policy_sha256 = screening.policy.sha256
        self.empty_proofs_sha256 = hashlib.sha256(_owned_admission_proofs_bytes(())).hexdigest()


class _ProfileRecord:
    """One private decoded input, annotation and canonical-byte lifetime.

    Strict external decoders own all nested data. Profile evaluation borrows
    that copy and its authenticated digest. Closing releases every reference.
    External codecs always validate independently, even while this owner lives.
    """
    def __init__(self, intake, context):
        self.active = False
        self.context = context
        sequence = type(intake) is ValidatedSequenceAdmission
        row = intake.row if sequence else intake
        if (type(row) is not EncodedAdmissionRow or row.proof_count != 0
                or row.proof_sha256 != context.empty_proofs_sha256):
            raise AdmissionProjectionError('profile record input is invalid')
        with _validation_session(isolated=True):
            value, content = row._decode_input_with_content(context.parents, context.binding_sha256)
            verified = None
            if sequence:
                from genosapiens.wgs_annotation import AnnotationStageEncodedRecord
                if type(intake.annotation) is not AnnotationStageEncodedRecord:
                    raise AdmissionProjectionError('profile annotation owner is invalid')
                intake.annotation.__post_init__()
                intake.annotation.encoded.__post_init__()
                verified = intake.annotation._decode_owned()
                intake._require_decoded(value, verified, context.parents)
        self.value, self.verified, self.content = value, verified, content
        self.sha256 = row.input_sha256
        self.active = True

    def input_digest(self, value):
        if not self.active or self.value is not value:
            raise AdmissionProjectionError('profile record owner differs')
        return self.sha256

    def close(self):
        self.active = False
        self.value = self.verified = self.content = self.context = None

    def __enter__(self):
        if not self.active:
            raise AdmissionProjectionError('profile record owner is closed')
        return self

    def __exit__(self, *_):
        self.close()


@dataclass(frozen=True, slots=True, repr=False)
class AdmissionProjectionResult:
    worker: int
    ordinal: int
    binding_sha256: str
    setup_sha256: str
    job_sha256: str
    job_size: int
    rows: tuple
    encode_count: int
    cpu_nanoseconds: int
    wall_nanoseconds: int
    idle_nanoseconds: int


@dataclass(frozen=True, slots=True, repr=False)
class _AdmissionResultSlice:
    result: AdmissionProjectionResult
    start: int
    stop: int


@dataclass(frozen=True, slots=True, repr=False)
class AdmissionProjectionWindow:
    """Live acquisition token; only its originating pool can retire it.

    Slices retain the complete encoded job until its final row is published.
    They never retain expanded evidence graphs or release the IPC reservation.
    """
    slices: tuple[_AdmissionResultSlice, ...]
    expanded_bytes: int

    @property
    def rows(self):
        return tuple(row for part in self.slices
                     for row in part.result.rows[part.start:part.stop])


def admission_worker_setup(context, binding, *, screening, defer_sequence_proofs=False):
    """Small authenticated metadata; no live DB, resource session or row graph."""
    from genosapiens.bounded_wgs_config import bounded_wgs_execution_config_as_json
    from genosapiens.longtr_admission import FORMAT as LONGTR_FORMAT
    range_mode = (context.read_intake is None
        and context.source_routing_context.contract.repeat_identity_format != LONGTR_FORMAT)
    _require_screening_context(screening, _binding(binding))
    if type(defer_sequence_proofs) is not bool:
        raise AdmissionProjectionError('admission intake mode is invalid')
    data = {'binding': json.loads(binding), 'screening': json.loads(screening.as_bytes()),
        'config': json.loads(bounded_wgs_execution_config_as_json(context.execution_config)),
        'projection': None, 'defer_sequence_proofs': defer_sequence_proofs}
    if range_mode:
        data['projection'] = {
            'source_root': str(context.source_root), 'source': _plain(context.source_manifest),
            'working': _plain(context.working_manifest), 'routing': _plain(context.source_routing_context),
            'prepared': _plain(context.prepared_index),
            'state_root': str(context.prepared.root.parent / 'state'),
            'retained_path': str(context.retained_projection_path)}
    content = canonical(data)
    if len(content) > _ADMISSION_JOB_BYTES:
        raise AdmissionProjectionError('admission worker metadata exceeds capacity')
    return content


def _open_admission_worker_inputs(content):
    from genosapiens.bounded_wgs_config import load_bounded_wgs_execution_config
    from genosapiens.wgs_annotation import open_retained_projection, _require_prepared_projection_parent
    from genosapiens.wgs_prepared_case import PreparedIndex
    from genosapiens.wgs_source_routing import SourceRoutingContext
    from genosapiens.wgs_ledger import LedgerManifest, WorkingSetManifest
    from genosapiens.wgs_admission import _source_header
    if type(content) is not bytes or len(content) > _ADMISSION_JOB_BYTES:
        raise AdmissionProjectionError('admission worker setup is invalid')
    data = json.loads(content)
    if (canonical(data) != content or set(data) != {'binding','config','projection','screening','defer_sequence_proofs'}
            or type(data['defer_sequence_proofs']) is not bool):
        raise AdmissionProjectionError('admission worker setup is invalid')
    binding = canonical(data['binding'])
    parents = _binding(binding)
    screening = load_admission_screening_context(canonical(data['screening']))
    _require_screening_context(screening, parents)
    config_bytes = canonical(data['config'])
    config = load_bounded_wgs_execution_config(config_bytes,parents['execution_limits_sha256'])
    projection = data['projection']
    if projection is None:
        return binding, config, None, None, screening
    if set(projection) != {'source_root','source','working','routing','prepared','state_root','retained_path'}:
        raise AdmissionProjectionError('admission worker source binding is invalid')
    source = _typed(projection['source'],LedgerManifest,decode=True)
    working = _typed(projection['working'],WorkingSetManifest,decode=True)
    routing = _typed(projection['routing'],SourceRoutingContext,decode=True)
    prepared = _typed(projection['prepared'],PreparedIndex,decode=True)
    if (source.sha256 != parents['source_manifest_sha256']
            or working.sha256 != parents['working_manifest_sha256']
            or routing.contract.sha256 != parents['source_routing_sha256']
            or prepared.sha256 != parents['prepared_index_sha256']):
        raise AdmissionProjectionError('admission worker source parents changed')
    _require_prepared_projection_parent(prepared,Path(projection['state_root']),source,working)
    header = _source_header(Path(projection['source_root']), source)
    if screening.header != header or screening.merged_intake is not None:
        raise AdmissionProjectionError('admission screening source header changed')
    arguments = dict(source_root=Path(projection['source_root']),source_manifest=source,
        working_manifest=working,source_routing_context=routing,header=header,
        definitions={(ns,name):(number,kind) for ns,name,number,kind in header.class_field_definitions},
        contig_lengths={v.contig:v.length for v in header.declared_contigs},
        execution_config=config,from_retained_projection=True)
    db = open_retained_projection(Path(projection['retained_path']),prepared_index=prepared,
        routing=routing,execution_config=config,
        expected_identities=working.summary.retained_unique_identities,
        expected_occurrences=working.summary.retained_occurrences)
    _settings(db,config,readonly=True)
    return binding, config, db, arguments, screening


@_validation_session(isolated=True)
def _produce_admission_record(inputs, entry, *, from_projection, parents, binding, screening,
                              defer_sequence_proofs=False):
    """One synchronous owner from reconstruction through encoded IPC bytes.

    next(inputs) runs inside the owner, so source constructors and annotation
    proof constructors can reuse copied immutable children. Nothing typed or
    validation-related crosses its return boundary; success and failure discard
    all walk entries. External loaders/encoders keep their isolated sessions.
    """
    from genosapiens.wgs_admission import replay_sequence_annotation_admission, validated_sequence_annotation_intake
    from genosapiens.wgs_annotation import AnnotationStageEncodedRecord, UnassessedVepRecord
    value = next(inputs)
    ordinal,identity,saved,transported = entry
    if (value.retained_ordinal != ordinal or value.identity_id != identity
            or from_projection and transported is not None):
        raise AdmissionProjectionError('admission record ownership changed')
    proofs = ()
    if value.transport == 'sequence_vep':
        if type(saved) is not AnnotationStageEncodedRecord:
            raise AdmissionProjectionError('admission annotation is missing')
        saved.__post_init__()
        saved.encoded.__post_init__()
        record = saved._decode_owned()
        if (record.request.identity_id != identity or record.request.retained_ordinal != ordinal
                or record.request.parent_working_set_manifest_sha256 != parents['working_manifest_sha256']):
            raise AdmissionProjectionError('admission annotation identity changed')
        if type(record) is UnassessedVepRecord:
            if (value.annotation_status != 'pending'
                    or record.request.source_allele_ids != value.source_allele_ids):
                raise AdmissionProjectionError('unassessed annotation input state differs')
            value = replace(value, annotation_status='unassessed', unresolved_reasons=(record.reason,))
            return EncodedAdmissionRow._encode_owned(value,(),binding=binding,screening=screening)
        if defer_sequence_proofs:
            value,_ = validated_sequence_annotation_intake(value,record,
                annotation_stage_manifest_sha256=parents['annotation_stage_sha256'],
                bundle_manifest_sha256=parents['bundle_manifest_sha256'])
            return ValidatedSequenceAdmission(
                EncodedAdmissionRow._encode_owned(value,(),binding=binding,screening=screening),saved)
        value,proofs = replay_sequence_annotation_admission(value,record,
            annotation_stage_manifest_sha256=parents['annotation_stage_sha256'],
            bundle_manifest_sha256=parents['bundle_manifest_sha256'])
    elif saved is not None:
        raise AdmissionProjectionError('admission annotation is unexpected')
    if value.transport == 'native_vep':
        content = transported if transported is not None else value.as_bytes()
        return PendingNativeAdmission(ordinal,identity,content,
            hashlib.sha256(content).hexdigest(),hashlib.sha256(binding).hexdigest(),len(value.source_allele_ids))
    return EncodedAdmissionRow._encode_owned(value,proofs,binding=binding,screening=screening)


@_validation_session(isolated=True)
def _produce_native_admission_record(row, session, *, binding, screening):
    value,proofs = session.enrich(row.decode(binding))
    return EncodedAdmissionRow._encode_owned(value,proofs,binding=binding,screening=screening)


@_validation_session(isolated=True)
def produce_profile_sequence_admission(value, verified, *, binding, screening,
        upstream_decision, selection, catalog, route_sources, clinical_session, regulatory_index, reference):
    """Fresh validated intake: construct broad sequence proofs only if advancing.

    Source/annotation input bytes are kept for every row. Only proof content
    changes after screening; the input identity/digest and upstream outcome do
    not. The coordinator persists the returned scope facts beside its ledger.
    """
    from genosapiens.wgs_admission import validated_sequence_annotation_intake
    parents = _binding(binding)
    _require_screening_context(screening, parents)
    enriched,provenance = validated_sequence_annotation_intake(value,verified,
        annotation_stage_manifest_sha256=parents['annotation_stage_sha256'],
        bundle_manifest_sha256=parents['bundle_manifest_sha256'])
    return _finish_profile_sequence_admission(enriched,verified,provenance,
        binding=binding,screening=screening,upstream_decision=upstream_decision,
        selection=selection,catalog=catalog,route_sources=route_sources,clinical_session=clinical_session,
        regulatory_index=regulatory_index,reference=reference)


@_validation_session(isolated=True)
def finalize_profile_sequence_admission(intake, *, binding, screening, upstream_decision,
        selection, catalog, route_sources, clinical_session, regulatory_index, reference):
    """Finalize a screened worker input without re-encoding its original bytes."""
    from genosapiens.wgs_admission import sequence_annotation_intake_provenance
    if type(intake) is not ValidatedSequenceAdmission:
        raise AdmissionProjectionError('validated sequence intake owner is required')
    value,verified=intake.decode(binding)
    parents=_binding(binding)
    provenance=sequence_annotation_intake_provenance(verified,
        annotation_stage_manifest_sha256=parents['annotation_stage_sha256'],
        bundle_manifest_sha256=parents['bundle_manifest_sha256'])
    return _finish_profile_sequence_admission(value,verified,provenance,
        binding=binding,screening=screening,upstream_decision=upstream_decision,
        selection=selection,catalog=catalog,route_sources=route_sources,clinical_session=clinical_session,
        regulatory_index=regulatory_index,reference=reference,original_row=intake.row)


def _finish_profile_sequence_admission(enriched, verified, provenance, *, binding, screening,
        upstream_decision, selection, catalog, route_sources, clinical_session, regulatory_index, reference,
        original_row=None, _clinical_window=None, _record_owner=None, _route_owner=None):
    from genosapiens.wgs_admission import sequence_annotation_intake_proofs
    from genosapiens.wgs_analysis_eligibility import evaluate_profile_intake
    if _record_owner is None:
        context = _ProfileRecordContext(binding, screening)
    else:
        if (type(_record_owner) is not _ProfileRecord or not _record_owner.active
                or _record_owner.value is not enriched or _record_owner.verified is not verified
                or _record_owner.context.binding != binding):
            raise AdmissionProjectionError('profile record owner differs')
        context = _record_owner.context
    parents = context.parents
    if (upstream_decision.admission_binding_sha256 != context.binding_sha256
            or upstream_decision.eligibility_policy_sha256 != context.policy_sha256
            or route_sources['molecular'] != parents['annotation_stage_sha256']
            or verified.request.parent_working_set_manifest_sha256 != parents['working_manifest_sha256']):
        raise AdmissionProjectionError('profile intake producer/screening binding differs')
    key = verified.request.request_id, 1
    molecular_sources = (*(t.source for t in enriched.allele.transcript_consequences or ()),
                         *verified.parsed.non_transcripts[key])
    facts,decision = evaluate_profile_intake(enriched,upstream_decision=upstream_decision,
        molecular_sources=molecular_sources,selection=selection,catalog=catalog,route_sources=route_sources,
        clinical_session=clinical_session,regulatory_index=regulatory_index,reference=reference,
        _clinical_window=_clinical_window, _record_owner=_record_owner, _route_owner=_route_owner)
    proofs = (sequence_annotation_intake_proofs(enriched,verified,provenance)
              if decision.disposition in ('selected','unresolved') else ())
    if original_row is None:
        row = EncodedAdmissionRow._encode_owned(enriched,proofs,binding=binding,screening=screening)
    else:
        proof_content = _owned_admission_proofs_bytes(proofs)
        row = replace(original_row,proof_blob=_pack_admission_blob(proof_content),proof_size=len(proof_content),
                      proof_sha256=hashlib.sha256(proof_content).hexdigest(),proof_count=len(proofs))
    if row.input_sha256 != upstream_decision.input_sha256:
        raise AdmissionProjectionError('profile intake original input changed')
    return row,replace(upstream_decision,proof_sha256=row.proof_sha256),facts,decision


def _run_admission_projection_worker(index, setup, jobs, results):
    """One spawned process owns strict decode, enrichment and final compression."""
    import pickle
    import signal
    import time
    from genosapiens.wgs_admission import _iter_admission_projection_rows
    def interrupted(_signal, _frame):
        raise KeyboardInterrupt
    db = None
    job = None
    try:
        os.setsid()
        signal.signal(signal.SIGTERM,interrupted)
        binding,config,db,arguments,screening = _open_admission_worker_inputs(setup)
        defer_sequence_proofs=json.loads(setup)['defer_sequence_proofs']
        job_rows = _DEFERRED_ADMISSION_JOB_ROWS if defer_sequence_proofs else _ADMISSION_JOB_ROWS
        parents = _binding(binding)
        binding_sha = hashlib.sha256(binding).hexdigest()
        setup_sha = hashlib.sha256(setup).hexdigest()
        results.put(('ready',index,setup_sha,binding_sha,os.getpid()))
        while True:
            idle_started=time.monotonic_ns()
            payload = jobs.get()
            idle=time.monotonic_ns()-idle_started
            if payload is None:
                results.put(('stopped',index,setup_sha,binding_sha,None))
                return
            if type(payload) is not bytes or len(payload) > _ADMISSION_JOB_BYTES:
                raise ValueError
            job = pickle.loads(payload)
            if (type(job) is not AdmissionProjectionJob or job.worker != index
                    or job.binding_sha256 != binding_sha or job.setup_sha256 != setup_sha
                    or type(job.rows) is not tuple or not 0 < len(job.rows) <= job_rows):
                raise ValueError
            started = time.process_time_ns()
            wall_started=time.monotonic_ns()
            if db is not None:
                inputs = _iter_admission_projection_rows(db,**arguments,
                    first_ordinal=job.rows[0][0],last_ordinal=job.rows[-1][0])
            else:
                inputs = (load_ledger_admission_input(item[3]) for item in job.rows)
            output = []
            encoded = 0
            with closing(inputs):
                for entry in job.rows:
                    row = _produce_admission_record(inputs,entry,from_projection=db is not None,
                        parents=parents,binding=binding,screening=screening,
                        defer_sequence_proofs=defer_sequence_proofs)
                    output.append(row)
                    encoded += type(row) in (EncodedAdmissionRow,ValidatedSequenceAdmission)
                # Exhaust the original generator's final coverage checks without
                # retaining a validation session across its yields or job rows.
                exhausted = object()
                if next(inputs,exhausted) is not exhausted:
                    raise AdmissionProjectionError('admission input coverage is extra')
            result = AdmissionProjectionResult(index,job.ordinal,binding_sha,setup_sha,
                hashlib.sha256(payload).hexdigest(),len(payload),tuple(output),encoded,
                time.process_time_ns()-started,time.monotonic_ns()-wall_started,idle)
            message = pickle.dumps(result,protocol=5)
            if len(message) > _ADMISSION_RESULT_BYTES:
                raise ValueError
            results.put(('result',index,setup_sha,binding_sha,message))
            del result,message,output,payload,job,inputs
            job = None
    except BaseException as error:
        from genosapiens.wgs_admission_diagnostics import worker_envelope
        results.put(('error',index,None,None,worker_envelope(error)))
    finally:
        # Cleanup may deliver another termination while multiprocessing flushes
        # its private queues during interpreter shutdown.
        signal.signal(signal.SIGTERM,signal.SIG_IGN)
        if db is not None:
            db.close()


class _AdmissionProjectionPool:
    """At most ten workers and one ordered consumer with bounded credits."""
    def __init__(self, config, *, setup, _worker_target=_run_admission_projection_worker,
                 _startup_seconds=300, _job_seconds=7200, _shutdown_seconds=30, total_inputs=None):
        import multiprocessing
        self.config,self.setup = config,setup
        self.job_rows = (_DEFERRED_ADMISSION_JOB_ROWS if json.loads(setup)['defer_sequence_proofs']
                         else _ADMISSION_JOB_ROWS)
        self.binding = canonical(json.loads(setup)['binding'])
        _binding(self.binding)
        self.binding_sha = hashlib.sha256(self.binding).hexdigest()
        self.setup_sha = hashlib.sha256(setup).hexdigest()
        self.worker_count = config.public_queries.concurrent_batches
        if self.worker_count != 10:
            raise AdmissionProjectionError('admission worker count must be ten')
        from genosapiens.wgs_admission_graph import live_worker_count
        self.worker_count = live_worker_count(total_inputs)
        self.startup_seconds = min(300,_startup_seconds)
        self.job_seconds = min(7200,_job_seconds)
        self.shutdown_seconds = min(30,_shutdown_seconds)
        context = multiprocessing.get_context('spawn')
        self.max_jobs = min(_ADMISSION_MAX_JOBS, 2*self.worker_count)
        self.results = context.Queue(maxsize=self.max_jobs)
        self.jobs = tuple(context.Queue(maxsize=2) for _ in range(self.worker_count))
        self.workers = tuple(context.Process(target=_worker_target,
            args=(i,setup,self.jobs[i],self.results),name=f'wgs-admission-{i}')
            for i in range(self.worker_count))
        self.outstanding,self.held,self.worker_pids = {},{},{}
        self.loads = [0]*self.worker_count
        self.next_job = self.next_consume = self.next_acquire = self.next_worker = 0
        self.acquire_offset = 0
        self.acquired = deque()
        self.started = self.closed = False
        self.subreaper_previous = None
        self.cleanup_failed = False
        self.abort_check=lambda:None
        self.held_bytes={}
        self.metrics = dict(worker_count=self.worker_count,max_outstanding=0,max_reorder_held=0,max_input_bytes=0,
            max_result_bytes=0,max_job_rows=0,max_reserved_input_bytes=0,
            max_reserved_result_bytes=0,encode_count=0,cpu_nanoseconds=0,
            worker_rows=[0]*self.worker_count,submitted_jobs=0,received_jobs=0,
            retired_jobs=0,acquired_windows=0,retired_windows=0,acquired_rows=0,writer_priority_yields=0,
            retired_rows=0,max_window_rows=0,max_window_bytes=0,max_acquired_windows=0,
            outstanding_jobs=0,reserved_input_bytes=0,reserved_result_bytes=0,
            acquired_window_bytes=0,received_rows=0,received_occurrences=0,
            received_native_inputs=0,received_native_occurrences=0,
            acquired_occurrences=0,retired_occurrences=0,cancelled_jobs=0,cancelled_rows=0,
            worker_busy_nanoseconds=[0]*self.worker_count,worker_idle_nanoseconds=[0]*self.worker_count,
            worker_cpu_nanoseconds=[0]*self.worker_count,worker_validated_rows=[0]*self.worker_count)

    def _receive(self, deadline, *, stopped=(), yield_when=None):
        import queue
        import time
        while True:
            self.abort_check()
            if time.monotonic() >= deadline:
                raise AdmissionProjectionError('admission worker deadline exceeded')
            self._check_worker_exit(stopped=stopped)
            try:
                if yield_when is not None and yield_when():
                    # Drain already available messages, but let a ready screened
                    # window reach its writers before waiting on a future job.
                    try:message=self.results.get_nowait()
                    except queue.Empty:return None
                else:
                    message = self.results.get(timeout=min(.1,max(.001,deadline-time.monotonic())))
            except queue.Empty:
                continue
            self._check_message(message)
            return message

    def _check_message(self, message):
        if (type(message) is not tuple or len(message) != 5 or type(message[1]) is not int
                or not 0 <= message[1] < self.worker_count):
            raise AdmissionProjectionError('admission worker failed')
        if message[0] == 'error':
            from genosapiens.wgs_admission_diagnostics import worker_error
            value = message[4] if message[2:4] == (None, None) else None
            raise worker_error(value) from None

    def _check_worker_exit(self, *, stopped=()):
        import queue
        if not any(w.exitcode is not None for i,w in enumerate(self.workers) if i not in stopped):
            return
        # Exit remains fatal. Inspect at most one queue capacity for error detail.
        # Other queued results cannot establish successful worker completion.
        for _ in range(min(_ADMISSION_MAX_JOBS, 2*self.worker_count)):
            try:
                message = self.results.get_nowait()
            except queue.Empty:
                break
            self._check_message(message)
        raise AdmissionProjectionError('admission worker exited')

    def check_running(self):
        from genosapiens.wgs_admission_diagnostics import tag
        try:
            if self.closed:
                raise AdmissionProjectionError('admission worker exited before closure')
            self._check_worker_exit()
        except BaseException as error:
            tag(error, 'admission_result_collect')
            raise

    def start(self):
        import ctypes
        import time
        if self.started or self.closed:
            raise AdmissionProjectionError('admission workers already started')
        self.started = True
        deadline = time.monotonic()+self.startup_seconds
        libc = ctypes.CDLL(None,use_errno=True)
        previous = ctypes.c_int()
        if libc.prctl(37,ctypes.byref(previous),0,0,0) or libc.prctl(36,1,0,0,0):
            raise AdmissionProjectionError('admission descendant ownership unavailable')
        self.subreaper_previous = previous.value
        for worker in self.workers:
            worker.start()
        while len(self.worker_pids) < self.worker_count:
            kind,index,setup,binding,pid = self._receive(deadline)
            if (kind != 'ready' or setup != self.setup_sha or binding != self.binding_sha
                    or index in self.worker_pids or pid != self.workers[index].pid):
                raise AdmissionProjectionError('admission worker startup binding changed')
            self.worker_pids[index] = pid

    @property
    def full(self):
        return len(self.outstanding) >= self.max_jobs

    def _message(self, rows, index):
        import pickle
        return pickle.dumps(AdmissionProjectionJob(index,self.next_job,self.binding_sha,
            self.setup_sha,tuple(rows)),protocol=5)

    def fits(self, rows):
        # Worker IDs 0..9 have identical pickle integer width.
        return (0 < len(rows) <= self.job_rows
            and len(self._message(rows,0)) <= _ADMISSION_JOB_BYTES)

    def submit(self, rows):
        import time
        if not self.started or self.closed or self.full or not self.fits(rows):
            raise AdmissionProjectionError('admission job exceeds capacity')
        for offset in range(self.worker_count):
            index = (self.next_worker+offset)%self.worker_count
            if self.loads[index] < 2:
                break
        else:
            raise AdmissionProjectionError("admission worker credit exhausted")
        self.next_worker = (index+1)%self.worker_count
        payload = self._message(rows,index)
        expected = tuple((entry[0],entry[1]) for entry in rows)
        if (any(type(o) is not int or o < 1 or type(i) is not str or not i for o,i in expected)
                or any(b[0] != a[0]+1 for a,b in zip(expected,expected[1:]))):
            raise AdmissionProjectionError('admission job keys are invalid')
        self.outstanding[self.next_job] = (index,expected,hashlib.sha256(payload).hexdigest(),
            len(payload),time.monotonic()+self.job_seconds)
        self.loads[index] += 1
        self.jobs[index].put(payload,timeout=1)
        self.next_job += 1
        metrics = self.metrics
        metrics['max_outstanding'] = max(metrics['max_outstanding'],len(self.outstanding))
        metrics['max_input_bytes'] = max(metrics['max_input_bytes'],len(payload))
        metrics['max_job_rows'] = max(metrics['max_job_rows'],len(rows))
        metrics['max_reserved_input_bytes'] = max(metrics['max_reserved_input_bytes'],len(self.outstanding)*_ADMISSION_JOB_BYTES)
        metrics['max_reserved_result_bytes'] = max(metrics['max_reserved_result_bytes'],len(self.outstanding)*_ADMISSION_RESULT_BYTES)
        metrics['submitted_jobs'] += 1
        self._reservations()

    def _reservations(self):
        self.metrics.update(outstanding_jobs=len(self.outstanding),
            outstanding_rows=sum(len(value[1]) for value in self.outstanding.values()),
            outstanding_input_bytes=sum(value[3] for value in self.outstanding.values()),
            reorder_jobs=len(self.held),reorder_rows=sum(len(value.rows) for value in self.held.values()),
            reorder_result_bytes=sum(self.held_bytes.values()),
            acquired_window_rows=sum(len(window.rows) for window in self.acquired),
            reserved_input_bytes=len(self.outstanding)*_ADMISSION_JOB_BYTES,
            reserved_result_bytes=len(self.outstanding)*_ADMISSION_RESULT_BYTES,
            acquired_window_bytes=sum(window.expanded_bytes for window in self.acquired))

    def _ordered_result(self, requested_ordinal, *, yield_when=None):
        import pickle
        if not self.started or self.closed or requested_ordinal not in self.outstanding:
            raise AdmissionProjectionError('admission result has no owner')
        while requested_ordinal not in self.held:
            pending = [v[4] for k,v in self.outstanding.items() if k not in self.held]
            message = self._receive(min(pending),yield_when=yield_when)
            if message is None:return None
            kind,index,setup,binding,payload = message
            if (kind != 'result' or setup != self.setup_sha or binding != self.binding_sha
                    or type(payload) is not bytes or len(payload) > _ADMISSION_RESULT_BYTES):
                raise AdmissionProjectionError('admission result envelope is invalid')
            result = pickle.loads(payload)
            if (type(result) is not AdmissionProjectionResult or result.ordinal not in self.outstanding
                    or result.ordinal in self.held):
                raise AdmissionProjectionError('admission result is unexpected')
            owner,keys,digest,size,_ = self.outstanding[result.ordinal]
            if (result.worker != index or owner != index or result.binding_sha256 != binding
                    or result.setup_sha256 != setup or result.job_sha256 != digest or result.job_size != size
                    or type(result.rows) is not tuple or len(result.rows) != len(keys)
                    or not _count(result.cpu_nanoseconds) or not _count(result.encode_count)
                    or not _count(result.wall_nanoseconds) or not _count(result.idle_nanoseconds)):
                raise AdmissionProjectionError('admission result binding is invalid')
            encoded = 0
            for row,(ordinal,identity) in zip(result.rows,keys,strict=True):
                if (type(row) not in (EncodedAdmissionRow,PendingNativeAdmission,ValidatedSequenceAdmission)
                        or row.retained_ordinal != ordinal or row.identity_id != identity
                        or row.binding_sha256 != binding or not _count(row.source_occurrences) or row.source_occurrences<1):
                    raise AdmissionProjectionError('admission result coverage is invalid')
                if type(row) in (EncodedAdmissionRow,ValidatedSequenceAdmission):
                    row.authenticated_bytes()
                elif (type(row.input_bytes) is not bytes or not 0 < len(row.input_bytes) <= MAX_ROW_BYTES
                        or hashlib.sha256(row.input_bytes).hexdigest() != row.input_sha256):
                    raise AdmissionProjectionError('native admission transport is invalid')
                encoded += type(row) in (EncodedAdmissionRow,ValidatedSequenceAdmission)
            if encoded != result.encode_count:
                raise AdmissionProjectionError('admission encode accounting is invalid')
            self.held[result.ordinal] = result
            self.held_bytes[result.ordinal]=len(payload)
            self.metrics['max_reorder_held'] = max(self.metrics['max_reorder_held'],len(self.held))
            self.metrics['max_result_bytes'] = max(self.metrics['max_result_bytes'],len(payload))
            self.metrics['received_jobs'] += 1
            self.metrics['received_rows']+=len(result.rows)
            self.metrics['received_occurrences']+=sum(row.source_occurrences for row in result.rows)
            native=tuple(row for row in result.rows if type(row) is PendingNativeAdmission)
            self.metrics['received_native_inputs']+=len(native)
            self.metrics['received_native_occurrences']+=sum(row.source_occurrences for row in native)
            self.metrics['worker_busy_nanoseconds'][index]+=result.wall_nanoseconds
            self.metrics['worker_idle_nanoseconds'][index]+=result.idle_nanoseconds
            self.metrics['worker_cpu_nanoseconds'][index]+=result.cpu_nanoseconds
            self.metrics['worker_validated_rows'][index]+=len(result.rows)
            self._reservations()
        return self.held[requested_ordinal]

    @staticmethod
    def _expanded_row_bytes(row):
        if type(row) is ValidatedSequenceAdmission:
            return (row.row.input_size + len(row.row.screening_facts) + MAX_ROW_BYTES +
                    len(row.annotation.encoded.canonical_request) + len(row.annotation.encoded.parsed_record))
        if type(row) is PendingNativeAdmission:
            # Native enrichment is coordinator-owned. Reserve its maximum final
            # input/proof/facts size before invoking the external native child.
            return 3 * MAX_ROW_BYTES
        return row.input_size + row.proof_size + len(row.screening_facts)

    def take_window(self, max_rows=_ADMISSION_WINDOW_ROWS, max_bytes=_ADMISSION_WINDOW_BYTES, *, yield_when=None):
        """Acquire an ordered bounded window without retiring any job credit.

        Drain available submitted work even below the cap (including EOF or a
        full twenty-job budget). Whole jobs stay together whenever they fit;
        an oversized job streams bounded row slices with one retained credit.
        A ready writer can interrupt a future-result wait. Keep received jobs
        held for coalescing, restore the acquisition cursor, and return None;
        never retire any credit or force a tiny query window on that yield.
        """
        if (type(max_rows) is not int or not 1 <= max_rows <= _ADMISSION_WINDOW_ROWS
                or type(max_bytes) is not int or not 1 <= max_bytes <= _ADMISSION_WINDOW_BYTES
                or len(self.acquired) >= _ADMISSION_MAX_WINDOWS):
            raise AdmissionProjectionError('admission window exceeds capacity')
        parts = []
        rows = size = 0
        original_cursor=(self.next_acquire,self.acquire_offset)
        while self.next_acquire in self.outstanding:
            result = self._ordered_result(self.next_acquire,yield_when=yield_when)
            if result is None:
                self.metrics['writer_priority_yields']+=1
                self.next_acquire,self.acquire_offset=original_cursor
                return None
            start = self.acquire_offset
            sizes = tuple(self._expanded_row_bytes(row) for row in result.rows[start:])
            if parts and (rows + len(sizes) > max_rows or size + sum(sizes) > max_bytes):
                break
            stop = start
            for row_size in sizes:
                if rows == max_rows or size + row_size > max_bytes:
                    break
                size += row_size
                rows += 1
                stop += 1
            if stop == start:
                raise AdmissionProjectionError('admission row exceeds window capacity')
            parts.append(_AdmissionResultSlice(result,start,stop))
            if stop == len(result.rows):
                self.next_acquire += 1
                self.acquire_offset = 0
            else:
                self.acquire_offset = stop
                break
            if rows == max_rows or size == max_bytes:
                break
        if not parts:
            if yield_when is not None:return None
            raise AdmissionProjectionError('admission window has no owner')
        window = AdmissionProjectionWindow(tuple(parts),size)
        self.acquired.append(window)
        metrics = self.metrics
        metrics['acquired_windows'] += 1
        metrics['acquired_rows'] += rows
        metrics['acquired_occurrences']+=sum(row.source_occurrences for row in window.rows)
        metrics['max_window_rows'] = max(metrics['max_window_rows'],rows)
        metrics['max_window_bytes'] = max(metrics['max_window_bytes'],size)
        metrics['max_acquired_windows'] = max(metrics['max_acquired_windows'],len(self.acquired))
        self._reservations()
        return window

    def consumed_window(self, window):
        # Object identity protects against forged/copied/partial tokens. The
        # deque fences ordered publication, including slices of the same job.
        if self.closed or not self.acquired or self.acquired[0] is not window:
            raise AdmissionProjectionError('admission window consumption is invalid')
        for part in window.slices:
            if part.result is not self.held.get(self.next_consume):
                raise AdmissionProjectionError('admission window job ownership changed')
            if part.stop == len(part.result.rows):
                self._retire(part.result)
        self.acquired.popleft()
        self.metrics['retired_windows'] += 1
        self.metrics['retired_rows'] += len(window.rows)
        self.metrics['retired_occurrences']+=sum(row.source_occurrences for row in window.rows)
        self._reservations()

    def _retire(self, result):
        # Caller writes every row before releasing this credit. Merely receiving
        # a result never admits another job into its reserved memory budget.
        if self.held.get(self.next_consume) is not result:
            raise AdmissionProjectionError('admission result consumption is invalid')
        index = self.outstanding[self.next_consume][0]
        self.loads[index] -= 1
        del self.outstanding[self.next_consume],self.held[self.next_consume]
        del self.held_bytes[self.next_consume]
        self.next_consume += 1
        self.metrics['encode_count'] += result.encode_count
        self.metrics['cpu_nanoseconds'] += result.cpu_nanoseconds
        self.metrics['worker_rows'][index] += len(result.rows)
        self.metrics['retired_jobs'] += 1

    def finish(self):
        """Fence every worker before completion, rejecting late/extra results."""
        import time
        if self.outstanding or self.held or self.acquired:
            raise AdmissionProjectionError('admission workers still have pending jobs')
        deadline = time.monotonic()+self.shutdown_seconds
        for jobs in self.jobs:
            jobs.put(None,timeout=1)
        stopped = set()
        while len(stopped) != self.worker_count:
            # A clean worker can exit immediately after queuing its fence. Read
            # available messages first; only a nonzero exit is an early failure.
            import queue
            try:
                message = self.results.get(timeout=max(.001,min(.1,deadline-time.monotonic())))
            except queue.Empty:
                if (time.monotonic() >= deadline
                        or any(w.exitcode not in (None,0) for w in self.workers)):
                    raise AdmissionProjectionError('admission worker shutdown failed')
                continue
            if (type(message) is not tuple or len(message) != 5
                    or message[0] != 'stopped' or message[1] not in self.worker_pids
                    or message[1] in stopped or message[2:] != (self.setup_sha,self.binding_sha,None)):
                raise AdmissionProjectionError('admission worker completion fence changed')
            stopped.add(message[1])
        for worker in self.workers:
            worker.join(timeout=max(0,deadline-time.monotonic()))
            if worker.exitcode != 0:
                raise AdmissionProjectionError('admission worker did not stop')
        try:
            self.results.get_nowait()
        except queue.Empty:
            return
        raise AdmissionProjectionError('admission worker emitted after completion')

    def close(self):
        import ctypes
        import signal
        import time
        if self.closed:
            return
        self.closed = True
        # Shared deadline for all children, with escalation halfway through. Each
        # child owns a process group; adopt/reap only its orphaned descendants.
        # Never launch replacement children or wait on unrelated process groups.
        deadline = time.monotonic()+self.shutdown_seconds
        for worker in self.workers:
            if worker.pid is not None:
                try:
                    if os.getpgid(worker.pid) == worker.pid:
                        os.killpg(worker.pid,signal.SIGTERM)
                    else:
                        worker.terminate()
                except ProcessLookupError:
                    pass
        soft = time.monotonic()+self.shutdown_seconds/2
        for worker in self.workers:
            if worker.pid is not None:
                worker.join(timeout=max(0,soft-time.monotonic()))
        for worker in self.workers:
            if worker.pid is not None:
                try:
                    os.killpg(worker.pid,signal.SIGKILL)
                except ProcessLookupError:
                    if worker.is_alive():
                        worker.kill()
        for worker in self.workers:
            if worker.pid is not None:
                worker.join(timeout=max(0,deadline-time.monotonic()))
                if worker.is_alive():
                    self.cleanup_failed = True
                while True:
                    try:
                        child,_ = os.waitpid(-worker.pid,os.WNOHANG)
                    except ChildProcessError:
                        break
                    if child == 0:
                        if time.monotonic() >= deadline:
                            self.cleanup_failed = True
                            break
                        time.sleep(.01)
        if self.subreaper_previous is not None:
            if ctypes.CDLL(None).prctl(36,self.subreaper_previous,0,0,0):
                self.cleanup_failed = True
        for channel in (*self.jobs,self.results):
            channel.cancel_join_thread()
            channel.close()
        self.metrics['cancelled_jobs']+=len(self.outstanding)
        self.metrics['cancelled_rows']+=self.metrics['acquired_rows']-self.metrics['retired_rows']
        self.outstanding.clear()
        self.held.clear()
        self.held_bytes.clear()
        self.acquired.clear()
        self.loads[:] = [0]*self.worker_count
        self._reservations()

    def __enter__(self):
        from genosapiens.wgs_admission_diagnostics import boundary
        try:
            with boundary('admission_result_collect', blocking=True):
                self.start()
            return self
        except BaseException as error:
            self.__exit__(type(error),error,None)
            raise

    def __exit__(self, kind, value, *_):
        from genosapiens.wgs_admission_diagnostics import close_preserving
        def close():
            self.close()
            if self.cleanup_failed:
                raise AdmissionProjectionError('admission worker cleanup failed')
        close_preserving(close, 'worker_close', primary=value)
