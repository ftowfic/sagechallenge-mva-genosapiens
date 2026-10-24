"""Disk-backed closure of every normalized scientific public reference.

Case resources are registered once; every complete context retains its own
public records and all source/member proofs. Typed references use the same
evidence ID, manifest and record semantics as EvidenceDocument. Complete pair
families inherit only the source templates that were actually checked here.
"""
import hashlib
from genosapiens.wgs_context_domain import next_context_ordinal, require_complete_domain
import json

from genosapiens.evidence_model import _validate_context_sources
from genosapiens.evidence_serialization import _canonical, _plain, _owned_canonical_bytes
from genosapiens.wgs_evidence_replay import EvidenceReplay
from genosapiens.wgs_evidence_store import EvidenceStoreError
from genosapiens.wgs_provenance import evidence_references


def _allele_digests(allele):
    """Full occurrence identity and invariant evidence across query batches.

    Only the known population adapter's batch envelope may vary. The exact
    variant-row hash, raw fields, resource bindings and every scientific field
    still participate in equality. Full occurrence bytes remain in each context.
    """
    from genosapiens.population_annotation import SCHEMA, JOINT_SOURCE_DATASET, SOURCE_RELEASE
    value = _plain(allele)
    digest = hashlib.sha256(_owned_canonical_bytes(value)).hexdigest()
    changed = False
    def sha(value):
        return type(value) is str and len(value) == 64 and all(c in '0123456789abcdef' for c in value)
    for row in value['population'] or ():
        source = row['source']
        fields = source['source_fields']
        if (source['source_dataset'] != JOINT_SOURCE_DATASET or source['source_release'] != SOURCE_RELEASE
                or type(fields) is not dict or fields.get('schema') != SCHEMA
                or not sha(fields.get('query_record_sha256'))):
            continue
        binding = fields.get('annotation_binding')
        if (type(binding) is not dict or type(binding.get('query_batch_id')) is not str
                or not binding['query_batch_id'] or not sha(binding.get('query_sha256'))):
            # A malformed row must not share the normalized domain even when
            # other rows in this allele have valid batch bindings.
            return digest, digest
        # Preserve the keys (and their presence) while normalizing only values.
        binding['query_batch_id'] = None
        binding['query_sha256'] = None
        changed = True
    identity = (hashlib.sha256(b'population-query-batch-independent-v1\n' +
                _owned_canonical_bytes(value)).hexdigest() if changed else digest)
    return digest, identity


def _reference_events(value):
    for ref in evidence_references(value):
        yield 'reference',ref.evidence_id,ref.record_id,ref.public_manifest_sha256


def _context_source_events(members, universe, scores, genes, observation_ids, *, seen_table):
    """Project live context sources; the ordered case owner checks conflicts.

    This is a private producing iterator, not a saved validation certificate.
    It leaves cross-context state and complete ordinal coverage with the case.
    Each event is scalar; public records are hashed in the producing process.
    """
    connection = members.connection
    if (members.closed or universe.closed or scores.closed
            or universe.members is not members or scores.universe is not universe
            or any(o.connection is not connection for o in (universe,scores,genes))):
        raise EvidenceStoreError("foreign scientific source projection")
    genes.require(members.writer,genes.prepared.evidence)
    root = connection.execute('SELECT kind,digest FROM roots WHERE name=?',
        (members.name+'.public',)).fetchone()
    if root != ('context-public-records-v1',members.public_sha256):
        raise EvidenceStoreError('scientific public collection binding changed')
    for row in members.iter_unique_public_records(seen_table):
        yield 'owner',row.evidence_id,row.public_manifest_sha256,'public'
        if row.record is not None:
            digest = hashlib.sha256((_canonical(row.record)+"\n").encode("ascii")).hexdigest()
            yield 'record',row.evidence_id,row.record.source_record_id,digest
    yield from _reference_events(members.seed)
    if members.work_filter is not None:
        for row in members.iter_work_decisions():
            yield 'original_member', row['allele_id'], row['original_input_sha256']
            if 'rescue_input' in row:
                from genosapiens.evidence_model import AlleleEvidence
                from genosapiens.evidence_serialization import _typed
                yield from _reference_events(_typed(row['rescue_input'], AlleleEvidence, decode=True))
    for member in members.iter_members():
        if member.allele.mechanism_ids:
            raise EvidenceStoreError("scientific member mechanism is unresolved")
        digest,identity_digest = _allele_digests(member.allele)
        yield 'allele',member.allele.allele_id,digest,identity_digest
        yield from _reference_events(member.allele)
    for decision in members.iter_decisions():
        yield from _reference_events(decision)
    if not set(universe.observation_ids) <= observation_ids:
        raise EvidenceStoreError("scientific hypothesis observation is unresolved")
    if universe.hypothesis_count:
        _validate_context_sources(members.seed.context,genes.prepared.evidence,
                                  prepared_gene_evidence=genes.prepared)
    for evidence,ledger in scores.singletons():
        yield from _reference_events((evidence,ledger))
    if universe.pairs is not None:
        yield from _reference_events(universe.pairs.context)
        for key,model,qualification in universe.pairs.member_templates():
            yield from _reference_events((model,qualification,scores.sources.charges(key)))
        yield from _reference_events(scores.sources.shared)


def _reference_query_events(queries):
    if queries.closed or queries.connection is None:
        raise EvidenceStoreError("closed scientific reference query projection")
    reference = hashlib.sha256(_canonical(queries.reference.binding).encode("ascii")).hexdigest()
    for row in queries.rows():
        digest = hashlib.sha256(_canonical(row["bases"]).encode("ascii")).hexdigest()
        yield 'query',reference,row["contig"],row["position"],row["length"],digest


def _source_context_header(members, universe, scores):
    from genosapiens.wgs_population_work_filter import work_domain
    return {"ordinal":members.ordinal,"context_id":members.context_id,
        **work_domain(members),
        "members":{"name":members.name,"sha256":members.sha256},
        "hypotheses":{"name":universe.name,"sha256":universe.sha256},
        "scores":{"name":scores.name,"sha256":scores.sha256}}


class CompleteSourceClosure:
    def __init__(self, writer, genes, phenotypes, *, name, context_domain=None, resource_prefix=None):
        self.context_domain = context_domain
        self.writer,self.connection,self.name = writer,writer.connection,name
        self.genes,self.phenotypes = genes,phenotypes
        self.resource_prefix = resource_prefix
        self.prefix = "source_closure_" + hashlib.sha256(name.encode()).hexdigest()
        self.tables,self.closed,self.published,self.context_count = [],False,False,0
        genes.require(writer,genes.prepared.evidence)
        if phenotypes is not None and phenotypes.connection is not self.connection:
            raise EvidenceStoreError("foreign phenotype source owner")
        self.observation_ids = frozenset(() if phenotypes is None else (o.observation_id for o in phenotypes.view.observations))
        try:
            for suffix,definition in (
                ("owners","id TEXT PRIMARY KEY,manifest TEXT,role TEXT"),
                ("public_sequences","digest TEXT PRIMARY KEY"),
                ("records","owner TEXT,id TEXT,digest TEXT,PRIMARY KEY(owner,id)"),
                ("alleles","id TEXT PRIMARY KEY,digest TEXT,identity_digest TEXT"),
                ("original_members","id TEXT PRIMARY KEY,digest TEXT"),
                ("queries","reference TEXT,contig TEXT,position INTEGER,length INTEGER,digest TEXT,PRIMARY KEY(reference,contig,position,length)"),
                ("contexts","ordinal INTEGER PRIMARY KEY,id TEXT UNIQUE,content TEXT"),
            ):
                table = self.prefix + "_" + suffix
                self.connection.execute(f"CREATE TABLE {table} ({definition})")
                self.tables.append(table)
            if resource_prefix is not None:
                from genosapiens.wgs_science_sources import _BoundSourcePrefix
                if type(resource_prefix) is not _BoundSourcePrefix:
                    raise EvidenceStoreError('foreign resource source prefix')
                resource_prefix.populate(self)
                return
            for public in genes.prepared.evidence:
                self._owner(public.evidence_id,public.public_manifest_sha256,"gene")
                for key,record in public.store.records.items():
                    self._record(public.evidence_id,key,record)
            if phenotypes is not None and phenotypes.profile is not None:
                public = phenotypes.profile
                self._owner(public.evidence_id,public.public_manifest_sha256,"phenotype")
                for annotation in public.store.annotations:
                    self._record(public.evidence_id,annotation.annotation_id,annotation)
            self._references(genes.prepared.evidence)
            if phenotypes is not None:
                self._references((phenotypes.profile,phenotypes.observations,phenotypes.view))
        except BaseException:
            writer.failed = True
            self.close()
            raise

    def _owner(self, key, manifest, role):
        row = self.connection.execute(f"SELECT manifest,role FROM {self.prefix}_owners WHERE id=?",(key,)).fetchone()
        if row is not None and row != (manifest,role):
            raise EvidenceStoreError("scientific public evidence owner conflicts")
        self.connection.execute(f"INSERT OR IGNORE INTO {self.prefix}_owners VALUES (?,?,?)",(key,manifest,role))

    def _record(self, owner, key, value):
        digest = hashlib.sha256((_canonical(value)+"\n").encode("ascii")).hexdigest()
        self._source_event(('record',owner,key,digest))

    def _references(self, value):
        for event in _reference_events(value):
            self._source_event(event)

    def _source_event(self, event):
        """Apply one live producer event with the original case conflict rules."""
        kind,*row = event
        if kind == 'owner' and len(row) == 3:
            self._owner(*row)
        elif kind == 'record' and len(row) == 3:
            owner,key,digest = row
            prior = self.connection.execute(f"SELECT digest FROM {self.prefix}_records WHERE owner=? AND id=?",(owner,key)).fetchone()
            if prior is not None and prior != (digest,):
                raise EvidenceStoreError("scientific public source record conflicts")
            self.connection.execute(f"INSERT OR IGNORE INTO {self.prefix}_records VALUES (?,?,?)",row)
        elif kind == 'reference' and len(row) == 3:
            owner,key,manifest = row
            prior = self.connection.execute(f"SELECT o.manifest FROM {self.prefix}_owners o JOIN {self.prefix}_records r "
                "ON o.id=r.owner WHERE o.id=? AND r.id=?",(owner,key)).fetchone()
            if prior != (manifest,):
                raise EvidenceStoreError("scientific evidence reference is unresolved")
        elif kind == 'allele' and len(row) == 3:
            key,digest,identity_digest = row
            prior = self.connection.execute(f"SELECT identity_digest FROM {self.prefix}_alleles WHERE id=?",(key,)).fetchone()
            if prior is not None and prior != (identity_digest,):
                raise EvidenceStoreError("scientific allele identity has conflicting context evidence")
            self.connection.execute(f"INSERT OR IGNORE INTO {self.prefix}_alleles VALUES (?,?,?)",row)
        elif kind == 'original_member' and len(row) == 2:
            key, digest = row
            prior = self.connection.execute(f"SELECT digest FROM {self.prefix}_original_members WHERE id=?", (key,)).fetchone()
            if prior is not None and prior != (digest,):
                raise EvidenceStoreError('original scientific member has conflicting context evidence')
            self.connection.execute(f"INSERT OR IGNORE INTO {self.prefix}_original_members VALUES (?,?)", row)
        elif kind == 'query' and len(row) == 5:
            key,digest = row[:-1],row[-1]
            prior = self.connection.execute(f"SELECT digest FROM {self.prefix}_queries "
                "WHERE reference=? AND contig=? AND position=? AND length=?",key).fetchone()
            if prior is not None and prior != (digest,):
                raise EvidenceStoreError("scientific reference query has conflicting context evidence")
            self.connection.execute(f"INSERT OR IGNORE INTO {self.prefix}_queries VALUES (?,?,?,?,?)",row)
        else:
            raise EvidenceStoreError("invalid producing source event")

    def add_context(self, members, universe, scores):
        if (self.closed or self.published or members.closed or universe.closed or scores.closed
                or universe.members is not members or scores.universe is not universe
                or any(owner.connection is not self.connection for owner in (members,universe,scores))
                or members.ordinal != next_context_ordinal(self)):
            raise EvidenceStoreError("foreign or incomplete scientific source context")
        try:
            self._append_context(_source_context_header(members,universe,scores),
                _context_source_events(members,universe,scores,self.genes,self.observation_ids,
                    seen_table=f'{self.prefix}_public_sequences'))
        except BaseException:
            self.writer.failed = True
            raise

    def _append_context(self, header, events):
        if self.closed or self.published or header["ordinal"] != next_context_ordinal(self):
            raise EvidenceStoreError("foreign or incomplete scientific source context")
        try:
            for event in events:
                self._source_event(event)
            self.context_count += 1
            self.connection.execute(f"INSERT INTO {self.prefix}_contexts VALUES (?,?,?)",
                (header["ordinal"],header["context_id"],_canonical(header)))
        except BaseException:
            self.writer.failed = True
            raise

    def add_reference_queries(self, queries):
        if (self.closed or self.published or queries.closed or queries.connection is not self.connection):
            raise EvidenceStoreError("foreign scientific reference query owner")
        try:
            for event in _reference_query_events(queries):
                self._source_event(event)
        except BaseException:
            self.writer.failed = True
            raise

    def publish(self):
        require_complete_domain(self)
        if self.closed or self.published:
            raise EvidenceStoreError("scientific source closure is unavailable")
        try:
            owners_sha = self.writer.add_sequence(self.name + ".owners","scientific-source-owners-v1",
                (tuple(row) for row in self.connection.execute(f"SELECT * FROM {self.prefix}_owners ORDER BY id")))
            records_sha = (self.writer.add_sequence(self.name + ".records","scientific-source-records-v1",
                (tuple(row) for row in self.connection.execute(f"SELECT * FROM {self.prefix}_records ORDER BY owner,id")))
                if self.resource_prefix is None else self.resource_prefix.publish_records(self, self.name + '.records'))
            contexts_sha = self.writer.add_sequence(self.name + ".contexts","scientific-source-contexts-v1",
                (json.loads(row[0]) for row in self.connection.execute(f"SELECT content FROM {self.prefix}_contexts ORDER BY ordinal")))
            alleles_sha = self.writer.add_sequence(self.name + ".alleles","scientific-source-alleles-v1",
                (tuple(row) for row in self.connection.execute(f"SELECT id,digest FROM {self.prefix}_alleles ORDER BY id")))
            queries_sha = self.writer.add_sequence(self.name + ".queries","scientific-source-reference-queries-v1",
                (tuple(row) for row in self.connection.execute(
                    f"SELECT * FROM {self.prefix}_queries ORDER BY reference,contig,position,length")))
            original = {}
            if any('population_work' in json.loads(content) for content, in self.connection.execute(
                    f'SELECT content FROM {self.prefix}_contexts ORDER BY ordinal')):
                digest = self.writer.add_sequence(self.name + '.original-members', 'scientific-original-members-v1',
                    (tuple(row) for row in self.connection.execute(f'SELECT * FROM {self.prefix}_original_members ORDER BY id')))
                original = {'original_members': {'name': self.name + '.original-members', 'sha256': digest}}
            self.sha256 = self.writer.add(self.name,"complete-scientific-source-closure-v1",{
                "schema":"genosapiens.wgs-scientific-source-closure","version":2 if original else 1,"context_count":self.context_count,
                **original,
                "genes":self.genes.reference,"phenotypes":None if self.phenotypes is None else self.phenotypes.reference,
                "owners":{"name":self.name + ".owners","sha256":owners_sha},
                "records":{"name":self.name + ".records","sha256":records_sha},
                "alleles":{"name":self.name + ".alleles","sha256":alleles_sha},
                "queries":{"name":self.name + ".queries","sha256":queries_sha},
                "contexts":{"name":self.name + ".contexts","sha256":contexts_sha}})
            self.published = True
            return self.sha256
        except BaseException:
            self.writer.failed = True
            raise

    def close(self):
        for table in reversed(self.tables):
            self.connection.execute("DROP TABLE " + table)
        self.tables.clear()
        self.closed = True

    def __enter__(self):
        return self

    def __exit__(self,*_):
        self.close()


class VerifiedSourceClosure(CompleteSourceClosure):
    def __init__(self, reader, reference, genes, phenotypes, contexts):
        if reference.kind != "complete-scientific-source-closure-v1":
            raise EvidenceStoreError("foreign scientific source closure root")
        reader._require_reference(reference)
        super().__init__(EvidenceReplay(reader),genes,phenotypes,name=reference.name)
        try:
            for members,universe,scores in contexts:
                self.add_context(members,universe,scores)
            self.publish()
        except BaseException:
            self.close()
            raise
