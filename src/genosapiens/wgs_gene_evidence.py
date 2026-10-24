"""One immutable case gene-resource owner, published as bounded record streams."""
from genosapiens.evidence_serialization import _owned_canonical_bytes
from genosapiens.evidence_model import DiseaseGeneEvidence, _typed_value, _validation_session
from genosapiens.genotype_evaluation import PreparedGeneEvidence, prepare_gene_evidence
from genosapiens.wgs_evidence_store import EvidenceStoreError, MAX_LOGICAL_BYTES

# Retained single-item reconstruction ceiling.
MAX_CASE_GENE_BYTES = 128 * 1024**2
_COLLECTIONS = ('records', 'assertions', 'views', 'relations', 'unassigned_record_ids')
_MAPPINGS = frozenset(('records', 'views'))
_SEQUENCE_KIND = 'case-gene-records-v2'


def _prepare_wire_gene_evidence(evidence):
    """Decode once within the existing strict immutable validation lifetime.

    Keep the scientific policy-bound preparation implementation unchanged.
    Its ordinary type check reuses only this session's freshly validated copy.
    """
    with _validation_session(isolated=True):
        owned = _typed_value(evidence, tuple[DiseaseGeneEvidence, ...], decode=True)
        return prepare_gene_evidence(owned)


def _bounded(values):
    for value in values:
        if len(_owned_canonical_bytes(value)) > MAX_CASE_GENE_BYTES:
            raise EvidenceStoreError('case gene evidence item exceeds resident capacity', capacity=True)
        yield value


def _write_sequence(writer, name, values, count):
    digest = writer.add_sequence(name, _SEQUENCE_KIND, _bounded(values))
    return dict(name=name, sha256=digest, count=count)


def _write_stores(writer, prepared, name):
    for ordinal, public in enumerate(prepared.evidence):
        prefix = f'{name}.store-{ordinal}'
        store = dict(source_binding=public.store.source_binding,
                     evidence_sha256=public.store.evidence_sha256)
        for field in _COLLECTIONS:
            values = getattr(public.store, field)
            rows = ((key, values[key]) for key in sorted(values)) if field in _MAPPINGS else iter(values)
            store[field] = _write_sequence(writer, prefix + '.' + field, rows, len(values))
        yield dict(evidence_id=public.evidence_id, public_manifest_sha256=public.public_manifest_sha256,
                   store=store, aggregations=_write_sequence(writer, prefix + '.aggregations',
                                                            iter(public.aggregations), len(public.aggregations)))


def _read_sequence(reader, descriptor, name):
    if (type(descriptor) is not dict or set(descriptor) != {'name', 'sha256', 'count'}
            or descriptor['name'] != name or type(descriptor['count']) is not int
            or not 0 <= descriptor['count'] <= MAX_LOGICAL_BYTES):
        raise EvidenceStoreError('invalid case gene sequence descriptor')
    reference = reader.reference(name, kind=_SEQUENCE_KIND)
    if reference.sha256 != descriptor['sha256']:
        raise EvidenceStoreError('case gene sequence root changed')
    count = 0
    for value in reader.iter_array(reference, max_item_bytes=MAX_CASE_GENE_BYTES):
        count += 1
        if count > descriptor['count']:
            raise EvidenceStoreError('case gene sequence count changed')
        yield value
    if count != descriptor['count']:
        raise EvidenceStoreError('case gene sequence count changed')


def _read_stores(reader, child, name, count):
    evidence = []
    for ordinal, public in enumerate(reader.iter_array(child, max_item_bytes=MAX_CASE_GENE_BYTES)):
        if (ordinal >= count or type(public) is not dict
                or set(public) != {'evidence_id', 'public_manifest_sha256', 'store', 'aggregations'}
                or type(public['store']) is not dict
                or set(public['store']) != {*_COLLECTIONS, 'source_binding', 'evidence_sha256'}):
            raise EvidenceStoreError('invalid case gene store descriptor')
        prefix = f'{name}.store-{ordinal}'
        store = public['store']
        for field in _COLLECTIONS:
            rows = _read_sequence(reader, store[field], prefix + '.' + field)
            if field in _MAPPINGS:
                values = {}
                previous = None
                for row in rows:
                    if (type(row) is not list or len(row) != 2 or type(row[0]) is not str
                            or (previous is not None and row[0] <= previous)):
                        raise EvidenceStoreError('case gene mapping keys changed')
                    previous = row[0]
                    values[row[0]] = row[1]
                store[field] = values
            else:
                store[field] = tuple(rows)
        public['aggregations'] = tuple(_read_sequence(reader, public['aggregations'], prefix + '.aggregations'))
        # The preparation factory strictly decodes this wire tree once, while
        # constructing its immutable owner and logical/resource indexes.
        evidence.append(public)
    if len(evidence) != count:
        raise EvidenceStoreError('case gene store count changed')
    return tuple(evidence)


class CaseGeneEvidence:
    def __init__(self, writer, prepared, *, name):
        if type(prepared) is not PreparedGeneEvidence:
            raise EvidenceStoreError('immutable case gene resource owner required')
        if prepared.canonical_size > MAX_LOGICAL_BYTES:
            raise EvidenceStoreError('case gene logical evidence exceeds capacity', capacity=True)
        self.writer,self.connection,self.prepared,self.name = writer,writer.connection,prepared,name
        prepared.require(prepared.evidence)
        try:
            evidence_sha = writer.add_sequence(name + '.stores', 'case-gene-stores-v2',
                                               _bounded(_write_stores(writer, prepared, name)))
            self.sha256 = writer.add(name, 'case-gene-evidence-v1', {
                'schema':'genosapiens.wgs-case-gene-evidence', 'version':2,
                'logical_sha256':prepared.sha256, 'canonical_size':prepared.canonical_size,
                'store_count':len(prepared.evidence),
                'stores':{'name':name + '.stores', 'sha256':evidence_sha}})
        except BaseException:
            writer.failed = True
            raise

    @property
    def reference(self):
        return {'name':self.name,'kind':'case-gene-evidence-v1','sha256':self.sha256,
                'logical_sha256':self.prepared.sha256}

    def require(self, writer, evidence):
        if writer.connection is not self.connection or self.connection is None:
            raise EvidenceStoreError('foreign case gene evidence owner')
        self.prepared.require(evidence)


class VerifiedCaseGeneEvidence(CaseGeneEvidence):
    def __init__(self, reader, reference, *, expected_logical_sha256):
        if reference.kind != 'case-gene-evidence-v1':
            raise EvidenceStoreError('foreign case gene evidence root')
        header = reader.reconstruct(reference, max_bytes=MAX_CASE_GENE_BYTES)
        if (set(header) != {'schema','version','logical_sha256','canonical_size','store_count','stores'}
                or header['schema'] != 'genosapiens.wgs-case-gene-evidence'
                or type(header['version']) is not int or header['version'] != 2
                or header['logical_sha256'] != expected_logical_sha256
                or type(header['canonical_size']) is not int
                or not 0 <= header['canonical_size'] <= MAX_LOGICAL_BYTES
                or type(header['store_count']) is not int or not 0 <= header['store_count'] <= MAX_LOGICAL_BYTES
                or type(header['stores']) is not dict or set(header['stores']) != {'name','sha256'}
                or header['stores']['name'] != reference.name + '.stores'):
            raise EvidenceStoreError('invalid case gene evidence header')
        child = reader.reference(header['stores']['name'], kind='case-gene-stores-v2')
        if child.sha256 != header['stores']['sha256']:
            raise EvidenceStoreError('case gene evidence store root changed')
        evidence = _read_stores(reader, child, reference.name, header['store_count'])
        prepared = _prepare_wire_gene_evidence(evidence)
        if (prepared.sha256 != expected_logical_sha256 or prepared.canonical_size != header['canonical_size']
                or len(prepared.evidence) != header['store_count']):
            raise EvidenceStoreError('case gene evidence logical closure changed')
        self.writer,self.connection,self.prepared,self.name = reader,reader.connection,prepared,reference.name
        self.sha256 = reference.sha256
