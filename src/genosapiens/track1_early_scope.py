"""Public, diagnosis-independent scope for original-source Track 1 preparation."""
from __future__ import annotations

from dataclasses import dataclass
from contextlib import ExitStack
import ast
import hashlib
import json
import os
from pathlib import Path
import re
import sqlite3
import tempfile
from array import array
from bisect import bisect_left

from genosapiens.bounded_json import load_bounded_json
from genosapiens.disease_gene_evidence import read_pinned_disease_gene_evidence
from genosapiens.track1_analysis_profiles import (
    build_profile_catalog_from_gene_evidence, load_analysis_profile_policy,
)
from genosapiens.wgs_ledger import canonical
from genosapiens.common_allele_public_indexes import RETENTION_EXACT_PROJECTION_VERSION

POLICY_ID = 'track1-early-scope/v1'
POLICY_PATH = 'config/track1/early-scope-policy.json'
MANIFEST_SCHEMA = 'genosapiens.track1-early-scope-manifest'
SOURCE_ROLES = ('catalog_source', 'coordinate_source', 'regulatory_source',
                'reference', 'exact_allele_source')
PRODUCER_MODULES = ('track1_early_scope', 'track1_analysis_profiles',
    'disease_gene_evidence', 'admission_links', 'retention_sources',
    'curated_alleles', 'variant_retention')
EXACT_PROJECTION_VERSION = RETENTION_EXACT_PROJECTION_VERSION
TABLES = ('binding', 'catalog_genes', 'source_intervals', 'merged_intervals',
          'known_allele_keys')
SQL_SCHEMA = '''
PRAGMA user_version=1;
CREATE TABLE binding(name TEXT PRIMARY KEY, value TEXT NOT NULL) WITHOUT ROWID;
CREATE TABLE catalog_genes(gene_id TEXT PRIMARY KEY, coordinate_status TEXT NOT NULL)
    WITHOUT ROWID;
CREATE TABLE source_intervals(kind TEXT NOT NULL, gene_id TEXT NOT NULL,
    contig TEXT NOT NULL, start INTEGER NOT NULL CHECK(start>=0),
    end INTEGER NOT NULL CHECK(end>start), source_record_id TEXT NOT NULL,
    source_manifest_sha256 TEXT NOT NULL,
    PRIMARY KEY(kind,gene_id,contig,start,end,source_record_id,source_manifest_sha256)) WITHOUT ROWID;
CREATE INDEX interval_sort ON source_intervals(contig,start,end);
CREATE TABLE merged_intervals(kind TEXT NOT NULL CHECK(kind IN ('region','known_sequence')),
    contig TEXT NOT NULL, start INTEGER NOT NULL CHECK(start>=0),
    end INTEGER NOT NULL CHECK(end>start), PRIMARY KEY(kind,contig,start)) WITHOUT ROWID;
CREATE TABLE known_allele_keys(key_kind TEXT NOT NULL CHECK(key_kind IN ('source','normalized')),
    contig TEXT NOT NULL, position INTEGER NOT NULL CHECK(position>0),
    reference TEXT NOT NULL, alternate TEXT NOT NULL, assertion_ref TEXT NOT NULL,
    PRIMARY KEY(key_kind,contig,position,reference,alternate,assertion_ref)) WITHOUT ROWID;
'''
MAX_ROW_BYTES = 8 * 1024**2
MAX_BATCH_BYTES = 64 * 1024**2
BATCH_ROWS = 1024
_POLICY = dict(
    schema='genosapiens.track1-early-scope-policy', version=1, policy_id=POLICY_ID,
    catalog='existing_track1_direct_catalog',
    gene_coordinates='open_targets_target_genomicLocation', gene_padding_bases=0,
    regulatory='all_explicit_stable_catalog_gene_links',
    known_alleles='global_exact_shared_known_disease_route',
    sequence_geometry='source_or_normalized_ref_span',
    structural_geometry='nominal_or_supplied_bounded_envelope',
    unknown_scope='defer_before_expensive_processing',
    outside_scope='defer_before_expensive_processing',
    phase_or_native_override=False, diagnosis_or_hpo_filter=False,
)


class EarlyScopeError(ValueError):
    """A static diagnostic for invalid scope resources or evidence."""


def _require(condition, message='Early scope resource or evidence is invalid.'):
    if not condition:
        raise EarlyScopeError(message)


def _digest(content):
    return hashlib.sha256(content).hexdigest()


def _sha(value):
    _require(type(value) is str and re.fullmatch('[0-9a-f]{64}', value) is not None)
    return value


def _relative(value):
    _require(type(value) is str and value and not Path(value).is_absolute()
             and '..' not in Path(value).parts and str(Path(value)) == value)
    return value


@dataclass(frozen=True, slots=True)
class EarlyScopePolicy:
    content: bytes

    def __post_init__(self):
        _require(type(self.content) is bytes and self.content == canonical(_POLICY),
                 'Early scope policy is invalid.')

    @property
    def sha256(self):
        return _digest(self.content)

    @property
    def document(self):
        return json.loads(self.content)


def load_early_scope_policy(repository_root):
    """Authenticate exact policy values and hash their canonical representation."""
    try:
        document, _ = load_bounded_json(Path(repository_root) / POLICY_PATH,
                                       maximum_bytes=16384)
        return EarlyScopePolicy(canonical(document))
    except (OSError, ValueError, TypeError):
        raise EarlyScopeError('Early scope policy is missing, stale or invalid.') from None


def load_scope_catalog(repository_root, reference_root):
    """Use the existing complete public catalog without diagnosis resolution."""
    from genosapiens import track1_public_validation as public
    repository = Path(repository_root)
    store = read_pinned_disease_gene_evidence(
        contract_path=repository / public.RARE_DISEASE_CONTRACT_PATH,
        reference_root=reference_root,
        expected_bundle_sha256=public.RARE_DISEASE_BUNDLE_SHA256)
    policy = load_analysis_profile_policy(repository / 'config/track1/analysis-profiles.json')
    catalog = build_profile_catalog_from_gene_evidence(store, policy=policy)
    _require(not catalog.coverage_issues, 'Early scope catalog coverage is unresolved.')
    return catalog


@dataclass(frozen=True, slots=True)
class ScopeRequest:
    content: bytes

    def __post_init__(self):
        document = json.loads(self.content)
        _require(canonical(document) == self.content
            and set(document) == {'schema', 'version', 'policy_sha256', 'catalog_policy_sha256',
                'retention_policy_sha256', 'source_bindings', 'producer_sha256s', 'exact_projection_version'}
            and document['schema'] == 'genosapiens.track1-early-scope-request'
            and type(document['version']) is int and document['version'] == 1
            and document['exact_projection_version'] == EXACT_PROJECTION_VERSION
            and set(document['source_bindings']) == set(SOURCE_ROLES)
            and set(document['producer_sha256s']) == set(PRODUCER_MODULES))
        for field in ('policy_sha256', 'catalog_policy_sha256', 'retention_policy_sha256'):
            _sha(document[field])
        for value in document['producer_sha256s'].values():
            _sha(value)

    @property
    def document(self):
        return json.loads(self.content)

    @property
    def sha256(self):
        return _digest(self.content)


def scope_request(repository_root, source_bindings):
    """Bind public inputs and producer code independently of runtime consumers."""
    from genosapiens.variant_retention import POLICY_SHA256
    repository = Path(repository_root)
    _require(type(source_bindings) is dict and set(source_bindings) == set(SOURCE_ROLES))
    policy = load_early_scope_policy(repository)
    catalog = load_analysis_profile_policy(repository / 'config/track1/analysis-profiles.json')
    return ScopeRequest(canonical(dict(schema='genosapiens.track1-early-scope-request', version=1,
        policy_sha256=policy.sha256, catalog_policy_sha256=catalog.sha256,
        retention_policy_sha256=POLICY_SHA256, source_bindings=source_bindings,
        producer_sha256s={name: _digest(_producer_bytes(repository, name))
                         for name in PRODUCER_MODULES},
        exact_projection_version=EXACT_PROJECTION_VERSION)))


def _producer_bytes(repository, name):
    content = (repository / 'src/genosapiens' / (name+'.py')).read_bytes()
    if name != 'track1_early_scope':
        return content
    # The colocated classifier is a runtime consumer of the immutable asset.
    # Preserve producer source text, including this selection implementation.
    lines = content.splitlines(keepends=True)
    nodes = ast.parse(content).body
    return b''.join(b''.join(lines[min([node.lineno, *(d.lineno for d in getattr(node, 'decorator_list', ()))])-1:node.end_lineno]) for node in nodes
        if not (isinstance(node, ast.ClassDef) and node.name in ('EarlyScopeDecision', 'EarlyScopeSession'))
        and not (isinstance(node, ast.ImportFrom) and node.module in ('array', 'bisect')))


def validate_scope_pin(value):
    _require(type(value) is dict and set(value) == {'manifest_path', 'manifest_sha256'})
    _sha(value['manifest_sha256'])
    _require(type(value['manifest_path']) is str and re.fullmatch(
        r'manifests/track1-early-scope/[0-9a-f]{64}\.json', value['manifest_path']) is not None)
    return dict(value)


@dataclass(frozen=True, slots=True)
class EarlyScopeManifest:
    content: bytes

    def __post_init__(self):
        value = json.loads(self.content)
        _require(canonical(value) == self.content and set(value) == {
            'schema', 'version', 'status', 'classification', 'request', 'request_sha256',
            'source_bindings', 'catalog_sha256', 'output', 'row_counts',
            'logical_sha256s', 'coverage'})
        _require(value['schema'] == MANIFEST_SCHEMA and type(value['version']) is int
            and value['version'] == 1 and value['status'] == 'complete'
            and value['classification'] == 'public')
        request = ScopeRequest(canonical(value['request']))
        _require(value['request_sha256'] == request.sha256
            and value['source_bindings'] == request.document['source_bindings'])
        _sha(value['catalog_sha256'])
        output = value['output']
        _require(set(output) == {'path', 'size_bytes', 'sha256', 'schema_version', 'identity_policy'}
            and output['path'] == f'derived/track1-early-scope/{request.sha256}/scope.sqlite'
            and type(output['size_bytes']) is int and output['size_bytes'] > 0
            and type(output['schema_version']) is int and output['schema_version'] == 1
            and output['identity_policy'] == 'trusted_completed_resource')
        _sha(output['sha256'])
        _require(set(value['row_counts']) == set(TABLES)
            and all(type(n) is int and n >= 0 for n in value['row_counts'].values())
            and set(value['logical_sha256s']) == set(TABLES)
            and type(value['coverage']) is dict)
        for digest in value['logical_sha256s'].values():
            _sha(digest)

    @property
    def document(self):
        return json.loads(self.content)

    @property
    def sha256(self):
        return _digest(self.content)


def _create_database(path):
    db = sqlite3.connect(path)
    db.execute('PRAGMA cache_size=-65536')
    db.execute('PRAGMA temp_store=FILE')
    db.execute('PRAGMA trusted_schema=OFF')
    db.executescript(SQL_SCHEMA)
    return db


def iter_known_keys(index, *, genes, catalog_sha256, reference_identity, reference_sha256):
    """Select exact assertions through the shared predicate, without diagnosis filters."""
    from dataclasses import asdict
    from genosapiens.retention_sources import clinical_from_index
    from genosapiens.curated_alleles import functional_from_index
    from genosapiens.variant_retention import LookupCoverage, RetentionFacts, select_variant
    bindings = index.retention_coverage()['source_bindings']
    _require(bindings['reference'] == reference_sha256, 'Early scope reference bindings differ.')
    coverage = tuple(LookupCoverage(dataset, bindings[key], True, 'complete')
        for dataset, key in (('ClinVar VCF', 'clinvar'), ('Open Targets evidence_eva', 'open_targets'),
                            ('gencc_submissions', 'curated'), ('evidence_clingen', 'curated')))
    genes = frozenset(genes)
    for row in index.iter_retention_exact_records():
        if row['source'] in {'clinvar', 'open_targets'}:
            fact = clinical_from_index(row, reference_identity=reference_identity,
                reference_sha256=reference_sha256, source_release_sha256=bindings[row['source']])
            clinical, functional = (fact,), ()
        else:
            _require(row['source'] == 'curated_functional')
            fact = functional_from_index(row)
            _require(fact.allele.reference_identity == reference_identity
                and fact.allele.reference_sha256 == reference_sha256
                and fact.source.release_sha256 == bindings['curated'])
            clinical, functional = (), (fact,)
        facts = RetentionFacts(_digest(canonical(asdict(fact.allele))), fact.allele,
            'assess', (fact.source.record_id,), (), clinical, functional, (), coverage)
        decision = select_variant(facts, established_genes=genes, catalog_sha256=catalog_sha256)
        if 'known_disease_allele' not in decision.routes:
            continue
        assertion_ref = canonical(asdict(fact.source)).decode('ascii')
        for kind, field in (('source', 'allele'), ('normalized', 'normalized_alias')):
            key = row[field]
            _require(type(key['chromosome']) is str and type(key['position']) is int
                and key['position'] > 0 and all(type(key[f]) is str and re.fullmatch('[ACGT]+', key[f])
                for f in ('reference', 'alternate')) and key['reference'] != key['alternate'])
            yield (kind, key['chromosome'], key['position'], key['reference'],
                   key['alternate'], assertion_ref)


def resolve_scope_sources(repository_root, index_pin):
    """Resolve source identities from the existing declarations, without source scans."""
    from genosapiens import track1_public_validation as public
    from genosapiens.installed_public_index import validate_index_pin
    repository = Path(repository_root)
    arguments, _ = load_bounded_json(repository / 'config/public-resources/track1-retention-index.json')
    regulatory, _ = load_bounded_json(repository / 'config/track1/regulatory-index-pin.json')
    pins = {}
    for role, name in (('coordinate_source', 'open_targets'), ('reference', 'reference')):
        pin = {key: arguments[name+'_'+key]
               for key in ('contract_path', 'contract_sha256', 'manifest_path', 'manifest_sha256')}
        _relative(pin['contract_path'])
        _relative(pin['manifest_path'])
        _sha(pin['manifest_sha256'])
        _, content = load_bounded_json(repository / pin['contract_path'])
        _require(_digest(content) == _sha(pin['contract_sha256']))
        pins[role] = pin
    _require(regulatory['reference_contract_sha256'] == pins['reference']['contract_sha256']
        and regulatory['reference_manifest_sha256'] == pins['reference']['manifest_sha256']
        and pins['reference']['manifest_sha256'] == public.REFERENCE_MANIFEST_SHA256)
    _relative(regulatory['index_path'])
    _sha(regulatory['manifest_sha256'])
    _, catalog_bytes = load_bounded_json(repository / public.RARE_DISEASE_CONTRACT_PATH)
    return dict(catalog_source=dict(contract_path=str(public.RARE_DISEASE_CONTRACT_PATH),
            contract_sha256=_digest(catalog_bytes), manifest_sha256=public.RARE_DISEASE_BUNDLE_SHA256),
        **pins, regulatory_source=regulatory, exact_allele_source=validate_index_pin(index_pin))


def _public_pin(repository, reference, value):
    from genosapiens.admission_links import PinnedPublicSource
    return PinnedPublicSource(repository / value['contract_path'], value['contract_sha256'],
        reference, reference / value['manifest_path'], value['manifest_sha256'])


def _batches(cursor):
    """Bound rows and serialized bytes for the compact coordinate projections."""
    # Coordinate projections contain only scalars. Each native database cursor
    # is consumed once. Clinical blobs use the stricter bounded iterator owner.
    while batch := cursor.fetchmany(min(BATCH_ROWS, MAX_BATCH_BYTES // MAX_ROW_BYTES)):
        total = 0
        for row in batch:
            size = len(canonical(row))
            _require(size <= MAX_ROW_BYTES, 'Early scope source row exceeds its bound.')
            total += size
            _require(total <= MAX_BATCH_BYTES, 'Early scope source batch exceeds its bound.')
            yield row


class _Rows:
    def __init__(self, db):
        self.db = db
        self.counts = dict.fromkeys(TABLES, 0)
        self.hashes = {table: hashlib.sha256() for table in TABLES}

    def emit(self, table, row):
        content = canonical(row)
        _require(len(content) <= MAX_ROW_BYTES)
        cursor = self.db.execute(f'INSERT INTO {table} VALUES ('+
            ','.join('?' for _ in row)+') ON CONFLICT DO NOTHING', row)
        if cursor.rowcount:
            self.counts[table] += 1
            self.hashes[table].update(content)
        return bool(cursor.rowcount)


def _contig(chromosome, lengths):
    _require(type(chromosome) is str and bool(chromosome))
    name = 'chrM' if chromosome in ('MT', 'M', 'chrM') else (
        chromosome if chromosome.startswith('chr') else 'chr'+chromosome)
    _require(name in lengths, 'Early scope coordinates use an unknown reference contig.')
    return name


def _interval(chrom, start, end, lengths):
    contig = _contig(chrom, lengths)
    _require(type(start) is int and type(end) is int and 0 <= start < end <= lengths[contig],
             'Early scope coordinates are malformed or outside the reference.')
    return contig, start, end


def _project_genes(rows, genes, lengths, source_sha256, output):
    found = {}
    for gene, chrom, start, end in rows:
        _require(type(gene) is str and re.fullmatch(r'ENSG[0-9]{11}(?:\.[0-9]+)?', gene))
        stable = gene.split('.')[0]
        if stable not in genes:
            continue
        _require(type(start) is int and type(end) is int)
        interval = _interval(chrom, start-1, end, lengths)
        _require(stable not in found or found[stable] == interval,
                 'Early scope gene coordinates conflict.')
        found[stable] = interval
        output.emit('source_intervals', ('gene', stable, *interval, gene, source_sha256))
    for gene in sorted(genes):
        output.emit('catalog_genes', (gene, 'present' if gene in found else 'missing_gene_coordinates'))
    return len(genes) - len(found)


def _project_regulatory(index, genes, lengths, output):
    index.require_stable()
    cursor = index.connection.execute('''SELECT
        CASE WHEN length(m.gene_id)<=256 THEN m.gene_id END,
        CASE WHEN length(e.chromosome)<=256 THEN e.chromosome END,e.start,e.end,
        CASE WHEN length(e.legacy_id)<=8192 THEN e.legacy_id END
        FROM memberships m JOIN elements e ON e.legacy_id=m.legacy_id
        WHERE m.target_kind='ensembl' ORDER BY m.gene_id,e.chromosome,e.start,e.end,e.legacy_id''')
    for gene, chrom, start, end, record in _batches(cursor):
        _require(type(gene) is str and re.fullmatch(r'ENSG[0-9]{11}(?:\.[0-9]+)?', gene))
        stable = gene.split('.')[0]
        if stable in genes:
            _require(type(record) is str and bool(record))
            interval = _interval(chrom, start, end, lengths)
            output.emit('source_intervals', ('regulatory', stable, *interval, record, index.manifest_sha256))
    index.require_stable()


def _project_regions(stack, repository, reference, bindings, catalog, output):
    import duckdb
    from genosapiens import admission_links as links
    from genosapiens.regulatory_index import open_regulatory_index
    files = links._Files(stack, reuse_published_outputs=True)
    pin = _public_pin(repository, reference, bindings['coordinate_source'])
    release, contract, manifest, artifacts, _ = links._open_targets_metadata(files, pin)
    targets = [a for a in artifacts if a.dataset == 'target']
    _require(bool(targets), 'Early scope gene-coordinate resource is absent.')
    outputs = {o.path: o for o in manifest.outputs}
    paths = []
    for artifact in targets:
        descriptor = outputs[f'{contract.output_prefix}/{artifact.relative_path}']
        path = release / artifact.relative_path
        files.open(path, descriptor.sha256, descriptor.size_bytes)
        paths.append(str(path))
    _, fasta, proof = links._reference(files, _public_pin(repository, reference, bindings['reference']))
    lengths = {name: item.length for name, item in fasta.index.items()}
    regulatory = bindings['regulatory_source']
    index = stack.enter_context(open_regulatory_index(reference / regulatory['index_path'],
        expected_manifest_sha256=regulatory['manifest_sha256']))
    _require(index.manifest['reference'] == proof, 'Early scope regulatory reference differs.')
    db = stack.enter_context(duckdb.connect(':memory:'))
    db.execute('SET threads=1')
    db.execute("SET memory_limit='256MB'")
    # The projection reads the target dataset once. DuckDB owns bounded sorting.
    cursor = db.execute('''SELECT CASE WHEN length(id)<=256 THEN id END,
        CASE WHEN length(genomicLocation.chromosome)<=256 THEN genomicLocation.chromosome END,genomicLocation.start,
        genomicLocation.end FROM read_parquet(?) ORDER BY id,genomicLocation.chromosome,
        genomicLocation.start,genomicLocation.end''', [paths])
    genes = frozenset(catalog.genes)
    missing = _project_genes(_batches(cursor), genes, lengths, pin.manifest_sha256, output)
    _project_regulatory(index, genes, lengths, output)
    return lengths, proof, missing


def _merge_intervals(output):
    # Disk-backed ordering is necessary production work, not an output audit.
    cursor = output.db.execute('''SELECT 'region',contig,start,end FROM source_intervals
        UNION ALL SELECT 'known_sequence',contig,start,end FROM known_spans
        ORDER BY 1,2,3,4''')
    previous = None
    for kind, contig, start, end in _batches(cursor):
        if previous is not None and previous[:2] == (kind, contig) and start <= previous[3]:
            previous = (*previous[:3], max(previous[3], end))
        else:
            if previous is not None:
                output.emit('merged_intervals', previous)
            previous = (kind, contig, start, end)
    if previous is not None:
        output.emit('merged_intervals', previous)


def _binding_rows(request, catalog_sha256):
    document = request.document
    return tuple(sorted((key, canonical(value).decode('ascii')) for key, value in dict(
        request_sha256=request.sha256, scope_policy_sha256=document['policy_sha256'],
        catalog_sha256=catalog_sha256,
        reference_manifest_sha256=document['source_bindings']['reference']['manifest_sha256'],
        source_bindings=document['source_bindings']).items()))


def build_early_scope(repository_root, reference_root, index_pin, destination):
    """Publish one immutable public scope, then atomically publish its small store pin."""
    from genosapiens.common_allele_public_indexes import (
        verify_common_allele_public_indexes, CommonAllelePublicIndexSession, _file_sha256, _publish_file,
    )
    from genosapiens.bounded_wgs_config import PublicQueryLimits
    from genosapiens.track1_setup import atomic_json
    repository, reference = Path(repository_root).resolve(), Path(reference_root).resolve(strict=True)
    bindings = resolve_scope_sources(repository, index_pin)
    request = scope_request(repository, bindings)
    manifest_path = f'manifests/track1-early-scope/{request.sha256}.json'
    manifest_file = reference / manifest_path
    if manifest_file.exists():
        _, content = load_bounded_json(manifest_file, reference_root=reference)
        pin = dict(manifest_path=manifest_path, manifest_sha256=_digest(content))
        completed = load_completed_early_scope(reference, pin)
        _require(completed.document['request_sha256'] == request.sha256)
        atomic_json(destination, dict(request_sha256=request.sha256, **pin))
        return pin
    catalog = load_scope_catalog(repository, reference)
    parent = reference / 'derived/track1-early-scope'
    parent.mkdir(parents=True, exist_ok=True)
    staging = Path(tempfile.mkdtemp(prefix='.incomplete-', dir=parent))
    # Keep interrupted output explicitly incomplete. It cannot supply a pin.
    with ExitStack() as stack:
        db = _create_database(staging / 'scope.sqlite')
        stack.callback(db.close)
        output = _Rows(db)
        db.execute('CREATE TEMP TABLE known_spans(contig TEXT,start INTEGER,end INTEGER)')
        for row in _binding_rows(request, catalog.sha256):
            output.emit('binding', row)
        lengths, proof, missing = _project_regions(stack, repository, reference, bindings, catalog, output)
        verified = verify_common_allele_public_indexes(reference_root=reference,
            manifest_path=index_pin['manifest_path'], expected_manifest_sha256=index_pin['manifest_sha256'],
            verify_sources=False, reuse_qualified_reconciliation=True)
        index = stack.enter_context(CommonAllelePublicIndexSession(verified,
            limits=PublicQueryLimits(1024, 10, MAX_BATCH_BYTES, 300)))
        for kind, chrom, pos, ref, alt, assertion in iter_known_keys(index, genes=catalog.genes,
                catalog_sha256=catalog.sha256, reference_identity=proof['source_reference_identity_id'],
                reference_sha256=proof['reference_manifest_sha256']):
            contig, start, end = _interval(chrom, pos-1, pos-1+len(ref), lengths)
            if output.emit('known_allele_keys', (kind, contig, pos, ref, alt, assertion)):
                db.execute('INSERT INTO known_spans VALUES (?,?,?)', (contig, start, end))
        _merge_intervals(output)
        db.commit()
    database = staging / 'scope.sqlite'
    output_descriptor = dict(path=f'derived/track1-early-scope/{request.sha256}/scope.sqlite',
        size_bytes=database.stat().st_size, sha256=_file_sha256(database), schema_version=1,
        identity_policy='trusted_completed_resource')
    manifest = EarlyScopeManifest(canonical(dict(schema=MANIFEST_SCHEMA, version=1,
        status='complete', classification='public', request=request.document, request_sha256=request.sha256,
        source_bindings=bindings, catalog_sha256=catalog.sha256, output=output_descriptor,
        row_counts=output.counts, logical_sha256s={key: h.hexdigest() for key, h in output.hashes.items()},
        coverage=dict(catalog_genes=len(catalog.genes), missing_gene_coordinates=missing,
            unmapped_or_unknown_locus_assertions='not_in_scope', reference=proof))))
    database.chmod(0o444)
    with database.open('rb') as stream:
        os.fsync(stream.fileno())
    final = parent / request.sha256
    # The ordinary setup lock owns this publication. Never replace old assets.
    if final.exists():
        # A crash can leave a closed database without the completion manifest.
        # Preserve that incomplete attempt before publishing a new attempt.
        _require(not final.is_symlink())
        orphan = Path(tempfile.mkdtemp(prefix='.incomplete-abandoned-', dir=parent))
        final.rename(orphan / 'asset')
    staging.rename(final)
    manifest_file.parent.mkdir(parents=True, exist_ok=True)
    manifest_stage = final / 'completed-manifest.json'
    with manifest_stage.open('xb') as stream:
        stream.write(manifest.content)
        stream.flush()
        os.fsync(stream.fileno())
    _publish_file(manifest_stage, manifest_file)
    pin = dict(manifest_path=manifest_path, manifest_sha256=manifest.sha256)
    atomic_json(destination, dict(request_sha256=request.sha256, **pin))
    return pin


def load_completed_early_scope(reference_root, pin):
    """Validate completed metadata and indexed bindings, without full database readback."""
    from genosapiens.common_allele_public_indexes import _regular
    try:
        pin = validate_scope_pin(pin)
        root = Path(reference_root).resolve(strict=True)
        _, content = load_bounded_json(root / pin['manifest_path'], reference_root=root)
        _require(_digest(content) == pin['manifest_sha256'])
        manifest = EarlyScopeManifest(content)
        value = manifest.document
        _require(pin['manifest_path'] == f"manifests/track1-early-scope/{value['request_sha256']}.json")
        descriptor = value['output']
        path = _regular(root, descriptor['path'], descriptor['sha256'], descriptor['size_bytes'],
                        trust_published=True)
        with path.open('rb') as stream:
            before = os.fstat(stream.fileno())
            db = sqlite3.connect(f'file:/proc/self/fd/{stream.fileno()}?mode=ro&immutable=1', uri=True)
            try:
                db.execute('PRAGMA trusted_schema=OFF')
                _require(db.execute('PRAGMA user_version').fetchone() == (1,))
                expected = _binding_rows(ScopeRequest(canonical(value['request'])), value['catalog_sha256'])
                _require(tuple(db.execute('SELECT name,value FROM binding ORDER BY name LIMIT 6')) == expected)
                for table in TABLES:
                    db.execute(f'SELECT 1 FROM {table} LIMIT 1').fetchone()
            finally:
                db.close()
            after = os.fstat(stream.fileno())
            _require((before.st_size, before.st_mtime_ns, before.st_ctime_ns) ==
                     (after.st_size, after.st_mtime_ns, after.st_ctime_ns))
        return manifest
    except (OSError, KeyError, TypeError, ValueError, sqlite3.Error):
        raise EarlyScopeError('Completed early scope is missing, stale or corrupt.') from None


@dataclass(frozen=True, slots=True)
class EarlyScopeDecision:
    source_allele_id: str
    normalized_allele_id: str | None
    status: str
    reason_codes: tuple[str, ...]
    matched_routes: tuple[str, ...]
    source_span: tuple[str, int, int] | None
    normalized_span: tuple[str, int, int] | None
    scope_policy_sha256: str
    scope_manifest_sha256: str
    reference_manifest_sha256: str
    matched_region: tuple[str, int, int] | None = None
    matched_exact_key: tuple[str, str, int, str, str, str] | None = None

    def __post_init__(self):
        _require(type(self.source_allele_id) is str and bool(self.source_allele_id)
            and (self.normalized_allele_id is None or type(self.normalized_allele_id) is str)
            and self.status in {'in_scope', 'outside_scope', 'unresolved_scope', 'not_applicable'}
            and type(self.reason_codes) is tuple and bool(self.reason_codes)
            and all(type(v) is str and bool(v) for v in self.reason_codes)
            and type(self.matched_routes) is tuple
            and self.matched_routes in ((), ('region',), ('known_disease_allele',))
            and (self.status == 'in_scope') == bool(self.matched_routes)
            and ('region' in self.matched_routes) == (self.matched_region is not None)
            and ('known_disease_allele' in self.matched_routes) == (self.matched_exact_key is not None))
        for sha in (self.scope_policy_sha256, self.scope_manifest_sha256, self.reference_manifest_sha256):
            _sha(sha)
        for span in (self.source_span, self.normalized_span, self.matched_region):
            if span is not None:
                _require(type(span) is tuple and len(span) == 3 and type(span[0]) is str
                    and type(span[1]) is int and type(span[2]) is int and 0 <= span[1] < span[2])
        if self.matched_exact_key is not None:
            key = self.matched_exact_key
            _require(type(key) is tuple and len(key) == 6 and key[0] in ('source','normalized')
                and type(key[2]) is int and key[2] > 0
                and all(type(key[i]) is str and key[i] for i in (1,3,4,5)))


class EarlyScopeSession:
    """One read-only scope owner with packed regional geometry and indexed exact keys."""
    def __init__(self, reference_root, pin, *, renderer, maximum_coordinate_bytes=64*1024**2):
        from genosapiens.submission_rendering import SubmissionContigRenderer
        _require(type(renderer) is SubmissionContigRenderer
            and type(maximum_coordinate_bytes) is int and 0 < maximum_coordinate_bytes <= 64*1024**2)
        self.manifest = load_completed_early_scope(reference_root, pin)
        value = self.manifest.document
        self.reference_sha256 = value['source_bindings']['reference']['manifest_sha256']
        _require(renderer.reference_manifest_sha256 == self.reference_sha256,
                 'Early scope and source reference bindings differ.')
        self.renderer = renderer
        self.names = dict(renderer.source_names)
        self.policy_sha256 = value['request']['policy_sha256']
        self.path = Path(reference_root).resolve() / value['output']['path']
        self._descriptor = os.open(self.path, os.O_RDONLY | os.O_NOFOLLOW)
        self._signature = self._stat_signature(os.fstat(self._descriptor))
        self.connection = None
        self.intervals = {}
        self.exact_queries = 0
        try:
            self.connection = sqlite3.connect(
                f'file:/proc/self/fd/{self._descriptor}?mode=ro&immutable=1', uri=True)
            self.connection.execute('PRAGMA cache_size=-65536')
            self.connection.execute('PRAGMA temp_store=FILE')
            self.connection.execute('PRAGMA query_only=ON')
            self.connection.execute('PRAGMA trusted_schema=OFF')
            self.connection.execute('PRAGMA mmap_size=0')
            size = 0
            cursor = self.connection.execute('SELECT kind,contig,start,end FROM merged_intervals ORDER BY kind,contig,start')
            for kind, contig, start, end in _batches(cursor):
                _require(kind in ('region','known_sequence') and type(contig) is str
                    and type(start) is int and type(end) is int and 0 <= start < end)
                if (kind, contig) not in self.intervals:
                    size += len(contig.encode('utf-8')) + len(kind)
                    self.intervals[kind, contig] = (array('Q'), array('Q'))
                starts, ends = self.intervals[kind, contig]
                _require(not ends or ends[-1] < start, 'Early scope intervals are not disjoint.')
                size += 16
                _require(size <= maximum_coordinate_bytes, 'Early scope coordinate capacity exceeded.')
                starts.append(start)
                ends.append(end)
            self.assert_stable()
        except BaseException:
            self.close()
            raise

    @staticmethod
    def _stat_signature(value):
        return value.st_dev, value.st_ino, value.st_size, value.st_mtime_ns, value.st_ctime_ns

    def assert_stable(self):
        _require(self._descriptor is not None and not self.path.is_symlink()
            and self._stat_signature(os.fstat(self._descriptor)) == self._signature
            and self._stat_signature(self.path.stat()) == self._signature,
            'Early scope resource changed during preparation.')

    def overlap(self, span, kind='region'):
        if span is None or (kind, span[0]) not in self.intervals:
            return None
        starts, ends = self.intervals[kind, span[0]]
        index = bisect_left(starts, span[2]) - 1
        if index >= 0 and ends[index] > span[1]:
            return span[0], starts[index], ends[index]
        return None

    def _span(self, contig, position, reference):
        entry = self.names.get(contig)
        _require(entry is not None, 'Early scope source contig lacks reference authentication.')
        _require(type(position) is int and position > 0
            and position+len(reference)-1 <= entry[1])
        canonical_contig = ('chrM' if entry[0] in ('M','MT') else
            'chr'+entry[0] if entry[0] in {*(str(n) for n in range(1,23)), 'X', 'Y'} else entry[0])
        return canonical_contig, position-1, position-1+len(reference)

    def classify(self, facts, route, structural_context=None):
        from genosapiens.wgs_source_routing import SourceOccurrenceFacts, SourceOccurrenceRoute
        _require(type(facts) is SourceOccurrenceFacts and type(route) is SourceOccurrenceRoute
            and facts.reference_manifest_sha256 == self.reference_sha256
            and facts.reference_identity_id == self.renderer.reference_binding.identity_id
            and facts.reference_role.value == self.renderer.reference_binding.role)
        self.assert_stable()
        normalized = facts.normalization.normalized_allele
        source_span = normalized_span = None
        def decision(status, reason, region=None, exact=None):
            return EarlyScopeDecision(facts.source.source_allele_id,
                normalized.normalized_allele_id if normalized else None, status, (reason,),
                ('region',) if region else ('known_disease_allele',) if exact else (),
                source_span, normalized_span, self.policy_sha256, self.manifest.sha256,
                self.reference_sha256, region, exact)
        if route.preparation_disposition in ('non_variant','quarantined'):
            return decision('not_applicable', 'source_route_not_applicable')
        if 'repeat' in route.native_class_tags or 'structural' in route.native_class_tags:
            return self._native_decision(facts, route, structural_context, decision)
        source_key = normalized_key = None
        if re.fullmatch('[ACGT]+', facts.record.reference) and re.fullmatch('[ACGT]+', facts.source.alternate):
            source_span = self._span(facts.record.contig, facts.record.position, facts.record.reference)
            source_key = ('source', source_span[0], facts.record.position, facts.record.reference, facts.source.alternate)
        if normalized is not None:
            normalized_span = self._span(normalized.contig, normalized.position, normalized.reference)
            normalized_key = ('normalized', normalized_span[0], normalized.position, normalized.reference, normalized.alternate)
        for span in (source_span, normalized_span):
            region = self.overlap(span)
            if region:
                return decision('in_scope', 'regional_overlap', region=region)
        for span, key in ((source_span, source_key), (normalized_span, normalized_key)):
            if key is not None and self.overlap(span, 'known_sequence'):
                self.exact_queries += 1
                try:
                    row = self.connection.execute('''SELECT assertion_ref FROM known_allele_keys
                        WHERE key_kind=? AND contig=? AND position=? AND reference=? AND alternate=?
                        ORDER BY assertion_ref LIMIT 1''', key).fetchone()
                except sqlite3.Error:
                    raise EarlyScopeError('Early scope exact-key query failed.') from None
                self.assert_stable()
                if row:
                    return decision('in_scope', 'exact_known_disease_allele', exact=(*key, row[0]))
        return decision('outside_scope', 'no_allowed_scope_match') if source_span or normalized_span else \
            decision('unresolved_scope', 'sequence_geometry_unavailable')

    def _native_decision(self, facts, route, context, decision):
        from dataclasses import replace
        from genosapiens.supplied_structural import parse_structural_record, structural_envelope
        from genosapiens.wgs_source_routing import REPEAT_IDENTITY_FORMAT
        from genosapiens import longtr_admission
        if context is None:
            return decision('unresolved_scope', 'native_geometry_unavailable')
        _require(type(context) is dict and set(context) == {
            'definitions', 'vcf_version', 'repeat_identity_format', 'header_lines'})
        lengths = {contig: value[1] for contig, value in self.names.items()}
        interval = envelope = None
        if 'repeat' in route.native_class_tags:
            try:
                if context['repeat_identity_format'] == REPEAT_IDENTITY_FORMAT:
                    from genosapiens.wgs_admission import _repeat_identity, WgsAdmissionError
                    try:
                        identity = _repeat_identity(facts, context['definitions'], lengths,
                            reference_role=facts.reference_role.value,
                            reference_identity_id=facts.reference_identity_id)
                    except WgsAdmissionError:
                        return decision('unresolved_scope', 'repeat_geometry_unavailable')
                elif context['repeat_identity_format'] == longtr_admission.FORMAT:
                    from types import SimpleNamespace
                    from genosapiens.vcf_qc import _parse_info
                    longtr_admission.validate_header(context['header_lines'], context['definitions'])
                    columns = facts.record.raw_line.rstrip('\r\n').split('\t')
                    source = SimpleNamespace(source_fields=dict(columns=columns,
                        info=_parse_info(columns[7]), format_keys=columns[8].split(':'),
                        sample_values={facts.selected_sample_name: columns[9+facts.selected_sample_index]}))
                    identity = longtr_admission.interpret(source, facts.sample,
                        role=facts.reference_role, reference_id=facts.reference_identity_id).parent
                else:
                    return decision('unresolved_scope', 'repeat_geometry_unavailable')
                interval = envelope = identity.interval
            except (ValueError, TypeError, KeyError, IndexError):
                return decision('unresolved_scope', 'repeat_geometry_unavailable')
        else:
            try:
                from types import SimpleNamespace
                # The source owner already checked record QC. This parser view
                # selects one sample without relabeling the complete QC record.
                selected = SimpleNamespace(raw_line=facts.record.raw_line,
                    alternates=facts.record.alternates, position=facts.record.position,
                    contig=facts.record.contig, sample_calls=(facts.sample,))
                calls = parse_structural_record(selected,
                    context['definitions'], lengths[facts.record.contig], role=facts.reference_role,
                    reference_id=facts.reference_identity_id, vcf_version=context['vcf_version'])
                identity = calls[0][1]
                interval, envelope = identity.interval, structural_envelope(identity)
            except (ValueError, TypeError, KeyError, OverflowError, IndexError):
                return decision('unresolved_scope', 'unsupported_structural_geometry')
        # Shared native parsers supply closed intervals. Preserve nominal
        # semantics, including when the possible envelope supplies the overlap.
        if interval is None:
            return decision('unresolved_scope', 'native_geometry_unavailable')
        source_span = self._span(interval.contig, interval.start, 'N')
        source_span = (*source_span[:2], interval.end)
        _require(interval.end <= lengths[interval.contig])
        spans = [source_span]
        if envelope is not None:
            possible = self._span(envelope.contig, envelope.start, 'N')
            _require(envelope.end <= lengths[envelope.contig])
            spans.append((*possible[:2], envelope.end))
        for number, span in enumerate(spans):
            region = self.overlap(span)
            if region:
                return replace(decision('in_scope', 'native_nominal_overlap' if number == 0 else
                    'native_bounded_envelope_overlap', region=region), source_span=source_span)
        return replace(decision('outside_scope', 'no_allowed_native_scope_match') if envelope is not None else
            decision('unresolved_scope', 'native_endpoint_uncertainty_unbounded'), source_span=source_span)

    def close(self):
        if self.connection is not None:
            self.connection.close()
            self.connection = None
        if self._descriptor is not None:
            os.close(self._descriptor)
            self._descriptor = None

    def __enter__(self):
        return self

    def __exit__(self, *_args):
        try:
            self.assert_stable()
        finally:
            self.close()
