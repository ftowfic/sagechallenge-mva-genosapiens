"""Bounded enrichment of every member of a globally closed context.

The disposable tables are query projections in the normalized evidence owner's
own database. The published member stream retains both the admitted input and
the enriched value, plus every population/clinical batch report and public record.
No per-context rich tuple is created, and no biological filtering occurs here.
"""
from collections.abc import Mapping
from dataclasses import fields, is_dataclass, replace
from enum import Enum
import hashlib
import json
import sqlite3
import sys

from genosapiens.bounded_wgs_config import BoundedWgsExecutionConfig
from genosapiens.evidence_model import IdentityStage, PublicSourcedEvidence, SourcedRecord, _typed
from genosapiens.evidence_serialization import _canonical, _owned_canonical_bytes
from genosapiens.population_annotation import PopulationRequest, PopulationResourceSession
from genosapiens.wgs_admission import LedgerAdmissionInput
from genosapiens.wgs_contexts import ClosedContext, ClosedPublicRecord
from genosapiens.wgs_evidence_store import EvidenceStoreError
from genosapiens.wgs_member_json import load_member_json

MAX_BATCH_BYTES = 16 * 1024**2
MAX_MEMBER_BYTES = 64 * 1024**2
MAX_RECORD_BYTES = 64 * 1024**2
MAX_MEMBER_CACHE_BYTES = 8 * 1024**2
MAX_MEMBER_CACHE_ENTRIES = 128


def _member_cache_size(content, value, *, limit=None):
    """Conservatively account private decoded data, stopping at the cache cap."""
    limit = MAX_MEMBER_CACHE_BYTES if limit is None else limit
    size = 2048 + sys.getsizeof(content)
    pending,seen = [value],set()
    while pending and size <= limit:
        current = pending.pop()
        if id(current) in seen:
            continue
        seen.add(id(current))
        size += sys.getsizeof(current)
        if size > limit:
            break
        if isinstance(current,Enum) or current is None or type(current) in (str,int,float,bool,bytes):
            continue
        if is_dataclass(current) and not isinstance(current,type):
            pending.extend(getattr(current,item.name) for item in fields(current))
        elif type(current) is tuple:
            pending.extend(current)
        elif isinstance(current,Mapping):
            # A mapping proxy's getsizeof excludes its private backing dict.
            # Count table/entry overhead conservatively, in addition to leaves.
            size += 1024 + 128 * len(current)
            if size > limit:
                break
            for key,child in current.items():
                pending.extend((key,child))
        else:
            return limit + 1
    return size


def _member_bytes(member):
    """Internal context member encoding; strict producer admission row bounds stay fixed."""
    content = _owned_canonical_bytes(member)
    if len(content) > MAX_MEMBER_BYTES:
        raise EvidenceStoreError("enriched context member exceeds capacity",capacity=True)
    return content


def _member_from_bytes(content):
    if not isinstance(content,bytes) or len(content) > MAX_MEMBER_BYTES:
        raise EvidenceStoreError("invalid enriched context member")
    value = _typed(load_member_json(content),LedgerAdmissionInput,decode=True)
    if _member_bytes(value) != content:
        raise EvidenceStoreError("noncanonical enriched context member")
    return value


class EnrichedContextEvidence:
    def __init__(self, writer, closed, *, population_session, execution_config,
                 clinical_annotator=None, query_cache=None, name=None, selection=None, failure_observer=None):
        def observe(stage):
            if failure_observer is not None:
                failure_observer.stage(stage)
        observe('member_enrichment')
        if (type(closed) is not ClosedContext or type(population_session) is not PopulationResourceSession
                or type(execution_config) is not BoundedWgsExecutionConfig or writer.connection is None
                or clinical_annotator is not None and not callable(clinical_annotator)):
            raise EvidenceStoreError("closed context enrichment inputs are invalid")
        self.writer, self.connection = writer, writer.connection
        self.context_id, self.ordinal, self.seed = closed.context_id, closed.ordinal, closed.seed
        self.member_count, self.input_bytes = closed.member_count, closed.input_bytes
        self.work_filter = None
        self._work_selection = selection
        if selection is not None:
            from genosapiens.wgs_population_work_filter import ContextPopulationWorkFilter
            if type(selection) is not ContextPopulationWorkFilter:
                raise EvidenceStoreError('authenticated population work selection required')
            selection.require_owner(writer, closed)
            self.seed = selection.seed
            self.member_count, self.input_bytes = selection.retained_count, selection.input_bytes
            self.work_filter = selection.reference
            self.work_filter_header = selection.header
        self.admitted = closed
        self.name = name or "context." + hashlib.sha256(self.context_id.encode()).hexdigest()
        self.prefix = "enriched_" + hashlib.sha256(self.context_id.encode()).hexdigest()
        self.tables, self.closed = [], False
        self._member_cache,self._member_cache_bytes = {},0
        self._public_base = None
        self._public_table = self.prefix + "_public"
        self.maximum_batch_bytes = self.maximum_batch_members = self.enriched_bytes = self.batch_count = 0
        self.has_clinical = clinical_annotator is not None
        try:
            for suffix,definition in (
                ("members","id TEXT PRIMARY KEY,ordinal INTEGER UNIQUE,content BLOB"),
                ("reports","ordinal INTEGER PRIMARY KEY,population TEXT,clinical TEXT"),
                ("public","evidence_id TEXT,manifest TEXT,record_id TEXT,content TEXT,PRIMARY KEY(evidence_id,record_id)"),
                ("decisions","id TEXT PRIMARY KEY,subject TEXT,content TEXT"),
            ):
                table = self.prefix + "_" + suffix
                self.connection.execute(f"CREATE TABLE {table} ({definition})")
                self.tables.append(table)
            self.connection.execute(f"CREATE INDEX {self.prefix}_decision_subject ON {self.prefix}_decisions(subject,id)")
            if closed.public_owner is None:
                for row in closed.iter_public_records():
                    self._public(row)
            else:
                shared = getattr(writer, '_context_public_shared', None)
                if shared is None:
                    shared = dict(owner=closed.public_owner, table='produced_shared_public', root=None, digest=None)
                    self.connection.execute('CREATE TABLE produced_shared_public '
                        '(evidence_id TEXT,manifest TEXT,record_id TEXT,content TEXT,'
                        'PRIMARY KEY(evidence_id,record_id)) WITHOUT ROWID')
                    ordinary = self._public_table
                    self._public_table = shared['table']
                    for row in closed.iter_public_records():
                        self._public(row)
                    self._public_table = ordinary
                    shared['last'] = self.connection.execute('SELECT evidence_id,record_id '
                        'FROM produced_shared_public ORDER BY evidence_id DESC,record_id DESC LIMIT 1').fetchone()
                    writer._context_public_shared = shared
                elif shared['owner'] is not closed.public_owner:
                    raise EvidenceStoreError('shared context public owner changed')
                self._public_base = shared
            if selection is not None:
                for row in selection.iter_public_records():
                    self._public(row)
            limit = min(execution_config.public_queries.batch_rows,population_session.limits.max_alleles)

            def batches():
                members,size = [],0
                observe('member_read')
                for member, input_size in closed.iter_sized_members():
                    if selection is not None and not selection.retains(member):
                        continue
                    if input_size > MAX_BATCH_BYTES:
                        raise EvidenceStoreError("context member exceeds enrichment window",capacity=True)
                    if members and (size+input_size > MAX_BATCH_BYTES or len(members) == limit):
                        yield tuple(members),size
                        observe('member_read')
                        members,size = [],0
                    members.append(member)
                    size += input_size
                if members:
                    yield tuple(members),size

            def rows():
                count = 0
                for members,size in batches():
                    self.batch_count += 1
                    self.maximum_batch_bytes = max(self.maximum_batch_bytes,size)
                    self.maximum_batch_members = max(self.maximum_batch_members,len(members))
                    requests = []
                    for member in members:
                        allele = member.allele
                        if allele is None or allele.population or allele.clinical:
                            raise EvidenceStoreError("context member is already enriched")
                        normalized = tuple(i.identity_id for i in allele.identities if i.stage is IdentityStage.NORMALIZED)
                        if len(normalized) > 1:
                            raise EvidenceStoreError("context member normalized identity is ambiguous")
                        requests.append(PopulationRequest(allele,normalized[0] if normalized else None))
                    requests = tuple(requests)
                    observe('population_query')
                    population = (population_session.query(requests) if query_cache is None
                                  else query_cache.population_results(requests,session=population_session))
                    original = {m.allele.allele_id:m.allele for m in members}
                    observe('population_bind')
                    self._check_enrichment(original,population.alleles,"population")
                    clinical = None
                    enriched = population.alleles
                    if clinical_annotator is not None:
                        observe('clinical_query')
                        clinical = clinical_annotator(enriched)
                        observe('clinical_bind')
                        self._check_enrichment({a.allele_id:a for a in enriched},clinical.alleles,"clinical")
                        enriched = clinical.alleles
                        for source in clinical.public_sourced_evidence:
                            if not source.records:
                                self._public(ClosedPublicRecord(source.evidence_id,source.public_manifest_sha256,None))
                            for record in source.records.values():
                                self._public(ClosedPublicRecord(source.evidence_id,source.public_manifest_sha256,record))
                    observe('member_encode')
                    population_report = _canonical(population.report)
                    clinical_report = None if clinical is None else _canonical({"coverage":clinical.coverage,"bindings":clinical.bindings})
                    if len(population_report) + len(clinical_report or "") > MAX_RECORD_BYTES:
                        raise EvidenceStoreError("context annotation report exceeds capacity",capacity=True)
                    self.connection.execute(f"INSERT INTO {self.prefix}_reports VALUES (?,?,?)",
                                            (self.batch_count,population_report,clinical_report))
                    by_id = {a.allele_id:a for a in enriched}
                    for member in members:
                        value = replace(member,allele=by_id[member.allele.allele_id])
                        content = _member_bytes(value)
                        self.connection.execute(f"INSERT INTO {self.prefix}_members VALUES (?,?,?)",
                            (value.allele.allele_id,value.retained_ordinal,content))
                        self.enriched_bytes += len(content)
                        count += 1
                        yield {"admitted":member,"enriched":value,"batch":self.batch_count}
                if count != self.member_count:
                    raise EvidenceStoreError("enriched context member coverage changed")

            self.members_sha256 = writer.add_sequence(self.name+".members","context-enriched-members-v1",rows())
            observe('member_publication')
            self.reports_sha256 = writer.add_sequence(self.name+".reports","context-annotation-reports-v1",self.iter_reports())
            first = self.connection.execute(f'SELECT evidence_id,record_id FROM {self._public_table} '
                'ORDER BY evidence_id,record_id LIMIT 1').fetchone()
            unchanged_public = self._public_base is not None and first is None
            # Admission owners precede the context-specific clinical owners.
            # Share that whole nonempty prefix even when the clinical suffix
            # differs. Interleaved owners retain the ordinary exact merge.
            prefix = (self._public_base is not None and self._public_base['last'] is not None
                and (first is None or self._public_base['last'][0] < first[0]))
            if prefix:
                shared = self._public_base
                args = ({'prefix':shared['digest']} if shared['digest'] is not None else
                        {'prefix_values':self._public_rows(shared['table'])})
                self.public_sha256, parts = writer.add_sequence_prefix(self.name+'.public',
                    'context-public-records-v1', self._public_rows(self._public_table), **args)
                shared['digest'] = parts[0]
                self._public_segments = [(parts[0],shared['table'])]
                if len(parts) == 2:
                    self._public_segments.append((parts[1],self._public_table))
            elif unchanged_public and self._public_base['root'] is not None:
                self.public_sha256 = writer.add_alias(self.name+'.public', 'context-public-records-v1',
                    self._public_base['root'])
            else:
                self.public_sha256 = writer.add_sequence(self.name+'.public','context-public-records-v1',self.iter_public_records())
                if unchanged_public:
                    self._public_base['root'] = self.name+'.public'
            if not prefix:
                self._public_segments = [(self.public_sha256,None)]
            def decisions():
                count = 0
                for value in closed.iter_decisions():
                    self.connection.execute(f"INSERT INTO {self.prefix}_decisions VALUES (?,?,?)",
                        (value.decision_id,value.subject_id,_canonical(value)))
                    count += 1
                    yield value
                if count != len(closed.seed.decision_ids):
                    raise EvidenceStoreError("context admission proof closure is incomplete")
            self.decisions_sha256 = writer.add_sequence(self.name+".decisions","context-admission-decisions-v1",decisions())
            self.sha256 = writer.add(self.name,"context-evidence-v1",{
                "schema":"genosapiens.wgs-context-evidence","version":1 if selection is None else 2,"ordinal":self.ordinal,
                "seed":self.seed,"member_count":self.member_count,"input_bytes":self.input_bytes,
                "enriched_bytes":self.enriched_bytes,"batch_count":self.batch_count,"has_clinical":self.has_clinical,
                **({} if self.work_filter is None else {'work_filter': self.work_filter}),
                **{key:{"name":self.name+"."+key,"sha256":getattr(self,key+"_sha256")}
                   for key in ("members","reports","public","decisions")}})
        except BaseException as error:
            # A failed scientific projection must never yield a publishable
            # prefix, including failures after a member root was completed.
            writer.failed = True
            primary = error if failure_observer is None else failure_observer.capture(error)
            from genosapiens.wgs_failures import run_cleanup
            run_cleanup((('context_store', self.close),), primary=primary)
            raise primary from None

    @staticmethod
    def _check_enrichment(original, enriched, role):
        if (len(enriched) != len(original) or {a.allele_id for a in enriched} != set(original)
                or any(replace(a,**{role:getattr(original[a.allele_id],role)}) != original[a.allele_id] for a in enriched)):
            raise EvidenceStoreError("context annotation changed source evidence")

    def iter_decisions(self, *, subject_id=None):
        from genosapiens.candidate_admission import AdmissionDecision
        self._require_open()
        sql = f"SELECT content FROM {self.prefix}_decisions"
        args = ()
        if subject_id is not None:
            sql += " WHERE subject=?"
            args = (subject_id,)
        for (content,) in self.connection.execute(sql + " ORDER BY id", args):
            yield _typed(json.loads(content), AdmissionDecision, decode=True)

    def iter_work_decisions(self):
        self._require_open()
        if self.work_filter is None:
            return
        if getattr(self, '_work_selection', None) is not None:
            yield from self._work_selection.decisions()
        else:
            for content, in self.connection.execute(f'SELECT content FROM {self.prefix}_work_selection ORDER BY ordinal'):
                yield json.loads(content)

    def population_work_assessment(self, key, allele, policy):
        self._require_open()
        selection = getattr(self, '_work_selection', None)
        # Independent diagnostic replay still recomputes the existing kernel.
        return None if selection is None else selection.population_assessment(key, allele, policy)

    def _public(self, row):
        key = "" if row.record is None else row.record.source_record_id
        _typed(PublicSourcedEvidence(row.evidence_id,row.public_manifest_sha256,
            {} if row.record is None else {key:row.record}),PublicSourcedEvidence)
        content = _canonical(row.record)
        if len(content) > MAX_RECORD_BYTES:
            raise EvidenceStoreError("context public record exceeds capacity",capacity=True)
        base = getattr(self, '_public_base', None)
        if base is not None:
            owner = self.connection.execute(f"SELECT manifest FROM {base['table']} WHERE evidence_id=? LIMIT 1",
                                            (row.evidence_id,)).fetchone()
            if owner is not None and owner[0] != row.public_manifest_sha256:
                raise EvidenceStoreError('context public owner conflicts')
            prior = self.connection.execute(f"SELECT content FROM {base['table']} WHERE evidence_id=? AND record_id=?",
                                            (row.evidence_id,key)).fetchone()
            if prior is not None and key:
                if prior[0] != content:
                    raise EvidenceStoreError('context public record conflicts')
                return
            if not key and owner is not None:
                return
        owner = self.connection.execute(f"SELECT manifest FROM {self._public_table} WHERE evidence_id=? LIMIT 1",
                                        (row.evidence_id,)).fetchone()
        if owner is not None and owner[0] != row.public_manifest_sha256:
            raise EvidenceStoreError("context public owner conflicts")
        if key:
            self.connection.execute(f"DELETE FROM {self._public_table} WHERE evidence_id=? AND record_id=''",(row.evidence_id,))
        elif owner is not None:
            return
        prior = self.connection.execute(f"SELECT content FROM {self._public_table} WHERE evidence_id=? AND record_id=?",
                                        (row.evidence_id,key)).fetchone()
        if prior is not None and prior[0] != content:
            raise EvidenceStoreError("context public record conflicts")
        self.connection.execute(f"INSERT OR IGNORE INTO {self._public_table} VALUES (?,?,?,?)",
                                (row.evidence_id,row.public_manifest_sha256,key,content))

    def _require_open(self):
        if self.closed:
            raise EvidenceStoreError("enriched context owner is closed")

    def iter_members(self, *, allele_order=False):
        self._require_open()
        order = "id" if allele_order else "ordinal"
        for (content,) in self.connection.execute(f"SELECT content FROM {self.prefix}_members ORDER BY {order}"):
            yield self._decoded_member(content)

    def _decoded_member(self, content):
        self._require_open()
        cached = self._member_cache.get(content)
        if cached is not None:
            return cached[0]
        value = _member_from_bytes(content)
        # Consumers repeatedly scan the same context. Preserve the retained
        # subset when full; an LRU would thrash on each larger-than-cache scan.
        remaining = MAX_MEMBER_CACHE_BYTES - self._member_cache_bytes
        if len(self._member_cache) < MAX_MEMBER_CACHE_ENTRIES and 2048 + sys.getsizeof(content) <= remaining:
            size = _member_cache_size(content,value,limit=remaining)
            if size <= remaining:
                self._member_cache[content] = value,size
                self._member_cache_bytes += size
        return value

    def member(self, key):
        self._require_open()
        row = self.connection.execute(f"SELECT content FROM {self.prefix}_members WHERE id=?",(key,)).fetchone()
        if row is None:
            raise EvidenceStoreError("enriched context member is missing")
        return self._decoded_member(row[0])

    def iter_reports(self):
        self._require_open()
        for ordinal,population,clinical in self.connection.execute(f"SELECT * FROM {self.prefix}_reports ORDER BY ordinal"):
            yield {"ordinal":ordinal,"population":json.loads(population),"clinical":None if clinical is None else json.loads(clinical)}

    def iter_public_records(self):
        self._require_open()
        sql = f"SELECT * FROM {self._public_table}"
        base = getattr(self, '_public_base', None)
        if base is not None:
            sql += (f" UNION ALL SELECT * FROM {base['table']} b WHERE NOT "
                f"(b.record_id='' AND EXISTS(SELECT 1 FROM {self._public_table} d WHERE d.evidence_id=b.evidence_id))")
        for evidence_id,manifest,key,content in self.connection.execute(sql+' ORDER BY evidence_id,record_id'):
            yield ClosedPublicRecord(evidence_id,manifest,None if not key else _typed(json.loads(content),SourcedRecord,decode=True))

    def _public_rows(self, table):
        self._require_open()
        for evidence_id,manifest,key,content in self.connection.execute(
                f'SELECT * FROM {table} ORDER BY evidence_id,record_id'):
            yield ClosedPublicRecord(evidence_id,manifest,None if not key else _typed(json.loads(content),SourcedRecord,decode=True))

    def iter_unique_public_records(self, seen_table):
        """Consume every new authenticated component before marking it seen.

        Callers own the disk-backed seen table and still check every context's
        references. The table contains content hashes, never external trust flags.
        """
        self._require_open()
        root = self.connection.execute('SELECT kind,digest FROM roots WHERE name=?',
            (self.name+'.public',)).fetchone()
        if root != ('context-public-records-v1',self.public_sha256):
            raise EvidenceStoreError('scientific public collection binding changed')
        for digest,table in self._public_segments:
            if self.connection.execute(f'SELECT 1 FROM {seen_table} WHERE digest=?',(digest,)).fetchone():
                continue
            yield from self.iter_public_records() if table is None else self._public_rows(table)
            self.connection.execute(f'INSERT INTO {seen_table} VALUES (?)',(digest,))

    def close(self):
        self._member_cache.clear()
        self._member_cache_bytes = 0
        for table in reversed(self.tables):
            self.connection.execute("DROP TABLE "+table)
        self.tables.clear()
        self.closed = True

    def __enter__(self):
        return self

    def __exit__(self,*_):
        self.close()


class VerifiedContextEvidence(EnrichedContextEvidence):
    """Authenticate all typed member/report/proof closure before query access."""
    def __init__(self, reader, reference, *, resources=None):
        from genosapiens.candidate_admission import CandidateContextSeed, AdmissionDecision
        from genosapiens.clinical_annotation import ClinicalCoverage
        from genosapiens.wgs_evidence_store import EvidenceReader
        from genosapiens.wgs_pair_domain import MAX_CONTEXT_MEMBERS
        if type(reader) is not EvidenceReader or reference.kind != "context-evidence-v1":
            raise EvidenceStoreError("foreign context evidence owner")
        self.writer,self.connection,self.name,self.sha256 = reader,reader.connection,reference.name,reference.sha256
        self.prefix = "verified_enriched_"+reference.sha256
        self.tables,self.closed = [],False
        self._member_cache,self._member_cache_bytes = {},0
        self._public_base = None
        self._public_table = self.prefix + "_public"
        try:
            header = reader.reconstruct(reference)
            expected = {"schema","version","ordinal","seed","member_count","input_bytes","enriched_bytes",
                        "batch_count","has_clinical","members","reports","public","decisions"}
            if header.get('version') == 2:
                expected.add('work_filter')
            if (set(header) != expected or header["schema"] != "genosapiens.wgs-context-evidence"
                    or type(header["version"]) is not int or header["version"] not in (1, 2)
                    or any(type(header[k]) is not int or header[k] < 0
                           for k in ("ordinal","member_count","input_bytes","enriched_bytes","batch_count"))
                    or header["ordinal"] == 0 or header["member_count"] > MAX_CONTEXT_MEMBERS
                    or type(header["has_clinical"]) is not bool):
                raise EvidenceStoreError("invalid context evidence header")
            self.seed = _typed(header["seed"],CandidateContextSeed,decode=True)
            self.context_id,self.ordinal = self.seed.context.context_id,header["ordinal"]
            self.member_count,self.input_bytes,self.enriched_bytes = header["member_count"],header["input_bytes"],header["enriched_bytes"]
            self.batch_count,self.has_clinical = header["batch_count"],header["has_clinical"]
            self.work_filter = header.get('work_filter')
            if self.work_filter is not None:
                from genosapiens.wgs_population_work_filter import load_work_selection
                self.work_filter_header = load_work_selection(reader, self.work_filter, self, resources)
            for suffix,definition in (
                ("members","id TEXT PRIMARY KEY,ordinal INTEGER UNIQUE,content BLOB"),
                ("admitted","id TEXT PRIMARY KEY,ordinal INTEGER UNIQUE,content BLOB,batch INTEGER"),
                ("reports","ordinal INTEGER PRIMARY KEY,population TEXT,clinical TEXT"),
                ("public","evidence_id TEXT,manifest TEXT,record_id TEXT,content TEXT,PRIMARY KEY(evidence_id,record_id)"),
                ("decisions","id TEXT PRIMARY KEY,subject TEXT,content TEXT"),
            ):
                table = self.prefix+"_"+suffix
                self.connection.execute(f"CREATE TABLE {table} ({definition})")
                self.tables.append(table)
            self.connection.execute(f"CREATE INDEX {self.prefix}_decision_subject ON {self.prefix}_decisions(subject,id)")
            kinds = {"members":"context-enriched-members-v1","reports":"context-annotation-reports-v1",
                     "public":"context-public-records-v1","decisions":"context-admission-decisions-v1"}
            refs = {}
            for key,kind in kinds.items():
                value = header[key]
                if set(value) != {"name","sha256"} or value["name"] != self.name+"."+key:
                    raise EvidenceStoreError("context evidence reference is foreign")
                ref = reader.reference(value["name"],kind=kind)
                if ref.sha256 != value["sha256"]:
                    raise EvidenceStoreError("context evidence reference changed")
                refs[key] = ref
            self.public_sha256 = refs["public"].sha256
            count,raw_bytes,enriched_bytes,prior_ordinal,prior_batch = 0,0,0,0,0
            for row in reader.iter_array(refs["members"],max_item_bytes=MAX_MEMBER_BYTES+5*1024**2):
                if set(row) != {"admitted","enriched","batch"} or type(row["batch"]) is not int:
                    raise EvidenceStoreError("invalid context member row")
                admitted = _typed(row["admitted"],LedgerAdmissionInput,decode=True)
                enriched = _typed(row["enriched"],LedgerAdmissionInput,decode=True)
                if (admitted.allele is None or admitted.allele.population or admitted.allele.clinical
                        or enriched.allele is None or admitted.retained_ordinal <= prior_ordinal
                        or not 1 <= row["batch"] <= self.batch_count or row["batch"] not in (prior_batch,prior_batch+1)
                        or replace(enriched,allele=admitted.allele) != admitted):
                    raise EvidenceStoreError("context member identity or source order changed")
                original = admitted.allele
                if replace(enriched.allele,population=original.population,clinical=original.clinical) != original:
                    raise EvidenceStoreError("context enrichment changed original evidence")
                if not self.has_clinical and enriched.allele.clinical != original.clinical:
                    raise EvidenceStoreError("context clinical evidence has no report")
                raw,content = admitted.as_bytes(),_member_bytes(enriched)
                if self.work_filter is not None:
                    selected = self.connection.execute(f'SELECT ordinal,native,input_sha,size,retained '
                        f'FROM {self.prefix}_work_selection WHERE id=?', (admitted.allele.allele_id,)).fetchone()
                    if selected != (admitted.retained_ordinal, admitted.allele.native_identity_id,
                            hashlib.sha256(raw).hexdigest(), len(raw), 1):
                        raise EvidenceStoreError('retained member differs from authenticated work selection')
                self.connection.execute(f"INSERT INTO {self.prefix}_members VALUES (?,?,?)",
                    (enriched.allele.allele_id,enriched.retained_ordinal,content))
                self.connection.execute(f"INSERT INTO {self.prefix}_admitted VALUES (?,?,?,?)",
                    (admitted.allele.allele_id,admitted.retained_ordinal,raw,row["batch"]))
                count += 1
                raw_bytes += len(raw)
                enriched_bytes += len(content)
                prior_ordinal,prior_batch = admitted.retained_ordinal,row["batch"]
            if (count,raw_bytes,enriched_bytes,prior_batch) != (self.member_count,self.input_bytes,self.enriched_bytes,self.batch_count):
                raise EvidenceStoreError("context member coverage changed")
            if self.connection.execute(f"SELECT 1 FROM {self.prefix}_admitted GROUP BY batch HAVING sum(length(content))>? LIMIT 1",
                                       (MAX_BATCH_BYTES,)).fetchone():
                raise EvidenceStoreError("context enrichment batch exceeds capacity",capacity=True)
            from genosapiens.evidence_model import AlleleReference
            references = tuple(AlleleReference(m.allele.allele_id,m.allele.native_identity_id) for m in self.iter_members())
            if references != self.seed.allele_refs:
                raise EvidenceStoreError("context members disagree with seed closure")
            count = 0
            for row in reader.iter_array(refs["reports"],max_item_bytes=MAX_RECORD_BYTES+1024):
                count += 1
                if set(row) != {"ordinal","population","clinical"} or type(row["ordinal"]) is not int or row["ordinal"] != count:
                    raise EvidenceStoreError("context query report order changed")
                population = _typed(row["population"],SourcedRecord,decode=True)
                clinical = row["clinical"]
                if self.has_clinical != (clinical is not None):
                    raise EvidenceStoreError("context clinical report coverage changed")
                if clinical is not None:
                    if set(clinical) != {"coverage","bindings"}:
                        raise EvidenceStoreError("invalid context clinical report")
                    _typed(clinical["bindings"],SourcedRecord,decode=True)
                    coverage = _typed(clinical["coverage"],tuple[ClinicalCoverage,...],decode=True)
                    for item in coverage:
                        match = self.connection.execute(f"SELECT batch FROM {self.prefix}_admitted WHERE id=?",(item.allele_id,)).fetchone()
                        if match != (count,) or self.member(item.allele_id).allele.native_identity_id != item.native_identity_id:
                            raise EvidenceStoreError("context clinical report has a foreign member")
                self.connection.execute(f"INSERT INTO {self.prefix}_reports VALUES (?,?,?)",
                    (count,_canonical(population),None if clinical is None else _canonical(clinical)))
            if count != self.batch_count:
                raise EvidenceStoreError("context report coverage changed")
            # Cache one shared public array, including a shared prefix under
            # distinct complete roots. Selection uses authenticated DAG metadata,
            # and lifetime is bounded by the enclosing full-byte replay scope.
            shared = getattr(reader, '_context_public_shared', None)
            scoped = reader._read_scope_depth or reader._active_derivative is not None
            if scoped and shared is None:
                common = self.connection.execute("WITH candidates(digest) AS ("
                    "SELECT digest FROM roots WHERE kind='context-public-records-v1' UNION ALL "
                    "SELECT e.child FROM roots r JOIN edges e ON e.parent=r.digest "
                    "JOIN nodes n ON n.digest=e.child WHERE r.kind='context-public-records-v1' AND n.kind='array') "
                    "SELECT c.digest FROM candidates c JOIN nodes n ON n.digest=c.digest "
                    "WHERE n.kind='array' GROUP BY c.digest HAVING count(*)>1 "
                    "ORDER BY (count(*)-1)*n.size DESC,n.size DESC,c.digest LIMIT 1").fetchone()
                shared = dict(digest=None if common is None else common[0], ready=False,
                              table='verified_shared_public',last=None)
                reader._context_public_shared = shared
            reader._check_files()
            reader._require_reference(refs['public'])
            node = reader._node(refs['public'].sha256)
            reuse = bool(scoped and shared is not None and shared['digest'] is not None
                and (shared['digest'] == refs['public'].sha256 or
                     node['kind'] == 'arrays' and len(node['items']) == 2 and node['items'][0] == shared['digest']))
            def consume(digest, previous=None):
                for row in reader._iter_array_values(digest,max_item_bytes=MAX_RECORD_BYTES+1024):
                    value = _typed(row,ClosedPublicRecord,decode=True)
                    key = (value.evidence_id,'' if value.record is None else value.record.source_record_id)
                    if previous is not None and (key <= previous or key[0] == previous[0] and previous[1] == ''):
                        raise EvidenceStoreError('context public record order changed')
                    self._public(value)
                    previous = key
                return previous
            if reuse:
                if not shared['ready']:
                    ordinary = self._public_table
                    self.connection.execute('CREATE TABLE verified_shared_public '
                        '(evidence_id TEXT,manifest TEXT,record_id TEXT,content TEXT,'
                        'PRIMARY KEY(evidence_id,record_id)) WITHOUT ROWID')
                    self._public_table = shared['table']
                    shared['last'] = consume(shared['digest'])
                    self._public_table = ordinary
                    shared['ready'] = True
                self._public_base = shared
                self._public_segments = [(shared['digest'],shared['table'])]
                if refs['public'].sha256 != shared['digest']:
                    suffix = node['items'][1]
                    consume(suffix,shared['last'])
                    self._public_segments.append((suffix,self._public_table))
            else:
                consume(refs['public'].sha256)
                self._public_segments = [(refs['public'].sha256,self._public_table)]
            reader._check_files()
            count = 0
            # Admission proofs describe the original closed context, including
            # members whose later scientific work was conclusively excluded.
            decision_subjects = self.prefix + ('_members' if self.work_filter is None else '_work_selection')
            for row in reader.iter_array(refs["decisions"],max_item_bytes=MAX_RECORD_BYTES):
                value = _typed(row,AdmissionDecision,decode=True)
                if (value.decision_id not in self.seed.decision_ids or self.connection.execute(
                        f"SELECT 1 FROM {decision_subjects} WHERE id=?",(value.subject_id,)).fetchone() is None):
                    raise EvidenceStoreError("context admission proof has a foreign owner")
                self.connection.execute(f"INSERT INTO {self.prefix}_decisions VALUES (?,?,?)",
                    (value.decision_id,value.subject_id,_canonical(value)))
                count += 1
            if count != len(self.seed.decision_ids):
                raise EvidenceStoreError("context admission proof closure is incomplete")
            self.decision_count = count
        except BaseException as error:
            self.close()
            if isinstance(error,(ValueError,TypeError,KeyError,AttributeError,sqlite3.Error)):
                raise EvidenceStoreError("invalid typed context evidence closure") from None
            raise

    def iter_admitted(self):
        from genosapiens.wgs_admission import load_ledger_admission_input
        self._require_open()
        for (content,) in self.connection.execute(f"SELECT content FROM {self.prefix}_admitted ORDER BY ordinal"):
            yield load_ledger_admission_input(content)
