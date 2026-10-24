"""Resolve legacy ClinVar unknowns from complete, ordered source loci.

The old projection emits exactly one row per ClinVar data record. Reconcile all
rows at a source position, in source record order, before resolving any unknown.
A nearby reference-only row alone never establishes the unknown's identity.
This selection-only correction does not change population filtering or indexes.
An N-containing sequence cannot match an allele with a different length change.
This rule requires known reference bases and excludes structural uncertainty.
"""
import json
import re

from genosapiens.admission_links import parse_clinvar_vcf_row
from genosapiens.tabix import TabixIndex


class ReferenceSiteError(ValueError):
    """The bounded lookup did not establish the source record and reference."""


def _require(value):
    if not value:
        raise ReferenceSiteError('ClinVar reference-site reconciliation failed')


def _signature(path):
    value = path.stat()
    _require(not path.is_symlink())
    return value.st_dev, value.st_ino, value.st_size, value.st_mtime_ns, value.st_ctime_ns


_STRUCTURAL_INFO = frozenset({'END', 'SVTYPE', 'SVLEN', 'IMPRECISE', 'CIPOS', 'CIEND',
    'CILEN', 'HOMLEN', 'HOMSEQ', 'BKPTID', 'MATEID', 'PARID', 'EVENT', 'CN', 'CICN',
    'MEINFO', 'METRANS', 'DUPHOMLEN', 'DUPHOMSEQ', 'DPADJ', 'CNADJ'})


def _sequence_length_nonmatch(source, query, *, chromosome, fetch):
    """Prove unequal haplotype lengths without guessing any N nucleotide."""
    if query is None:
        return False
    from genosapiens.common_allele_public_indexes import PublicIndexQuery
    _require(type(query) is PublicIndexQuery and query.chromosome == chromosome
        and type(query.position) is int and query.position > 0
        and all(re.fullmatch('[ACGT]{1,10000}', value) is not None
                for value in (query.reference, query.alternate))
        and query.reference != query.alternate
        and all(getattr(query, name) is None for name in
                ('normalized_alias', 'source_span', 'normalized_span', 'native_scope')))
    fields, info = source['columns'], source['info']
    ref, alt = fields[3:5]
    if (re.fullmatch('[ACGT]{1,10000}', ref) is None
            or re.fullmatch('[ACGTN]{1,10000}', alt) is None or 'N' not in alt
            or _STRUCTURAL_INFO.intersection(info)):
        return False
    delta = len(alt) - len(ref)
    source_class = info.get('CLNVC')
    if source_class is not None and not (
            source_class == 'single_nucleotide_variant' and len(ref) == len(alt) == 1
            or source_class == 'Deletion' and delta < 0
            or source_class == 'Insertion' and delta > 0
            or source_class == 'Indel'):
        return False
    if delta == len(query.alternate) - len(query.reference):
        return False
    _require(fetch(chromosome, int(fields[1]), len(ref)) == ref)
    _require(fetch(chromosome, query.position, len(query.reference)) == query.reference)
    return True


def reconcile_locus(projected, raw_rows, *, artifact, chromosome, position, fetch, query=None):
    """Return resolved source ordinals after complete, ordered locus agreement."""
    _require(0 < len(projected) == len(raw_rows) <= 64)
    ordinals, parsed = [], []
    for row, raw in zip(projected, raw_rows, strict=True):
        locator = row['source_locator']
        ordinal = locator['record_number']
        _require(type(ordinal) is int and ordinal > 0
                 and locator['path'] == artifact.path and locator['sha256'] == artifact.sha256
                 and row['source'] == 'clinvar')
        ordinals.append(ordinal)
        source = parse_clinvar_vcf_row(raw)
        fields = source['columns']
        _require(fields[0] == chromosome and int(fields[1]) == position)
        if row['record_type'] == 'exact':
            stored = row['source_records']['clinvar']
            _require(list(stored['columns']) == list(fields) and stored['info'] == source['info'])
            _require(row['allele'] == dict(chromosome=chromosome, position=position,
                                          reference=fields[3], alternate=fields[4]))
        else:
            _require(row['record_type'] == 'unknown' and row['reason'] == 'unassessable_source_record')
            _require(row['affected_locus'] == dict(chromosome=chromosome, start=position,
                                                  end=position + len(fields[3]) - 1))
        parsed.append(source)
    _require(ordinals == list(range(ordinals[0], ordinals[0] + len(ordinals))))
    resolved = {}
    for row, source, ordinal in zip(projected, parsed, ordinals, strict=True):
        fields = source['columns']
        if row['record_type'] == 'unknown' and fields[4] == '.':
            _require(re.fullmatch('[ACGT]{1,10000}', fields[3]) is not None)
            _require(fetch(chromosome, position, len(fields[3])) == fields[3])
            # Keep the original locator and source text in the bounded cache.
            resolved[ordinal] = dict(row, record_type='reference_site',
                reason='clinvar_no_alternate_allele', source_records=dict(clinvar=source))
        elif row['record_type'] == 'unknown' and _sequence_length_nonmatch(
                source, query, chromosome=chromosome, fetch=fetch):
            # This proof is query-specific. Retain both original source text and
            # the exact query in its cache entry. It is not a benign assertion.
            resolved[ordinal] = dict(row, record_type='nonmatching_sequence',
                reason='clinvar_ambiguous_sequence_length_differs', source_records=dict(clinvar=source),
                query_allele=dict(chromosome=query.chromosome, position=query.position,
                    reference=query.reference, alternate=query.alternate))
    return resolved


class ClinvarReferenceSites:
    """One clinical session's lazy, bounded reader for old immutable indexes."""
    def __init__(self, session):
        self.session = session
        self.tabix = None
        self.files = None

    def assert_stable(self):
        if self.files is not None:
            artifacts = (self.session._clinvar_artifact, self.session._clinvar_contract.artifact('tbi'))
            _require(tuple(_signature(self.session.clinvar.root / item.path) for item in artifacts) == self.files)

    def resolve(self, records, *, budget, query=None):
        pending = [row for row in records if row['record_type'] == 'unknown']
        if not pending or any(row.get('source') != 'clinvar'
                or row.get('reason') != 'unassessable_source_record' for row in pending):
            return None
        session = self.session
        artifact = session._clinvar_artifact
        tbi = session._clinvar_contract.artifact('tbi')
        paths = tuple(session.clinvar.root / value.path for value in (artifact, tbi))
        current = tuple(_signature(path) for path in paths)
        _require(current[0][2] == artifact.size_bytes and current[1][2] == tbi.size_bytes)
        if self.files is None:
            self.tabix = TabixIndex.load(paths[1], budget=budget)
            self.tabix.require_vcf_layout()
            self.files = current
        _require(current == self.files)
        resolved = {}
        loci = {(row['affected_locus']['chromosome'], row['affected_locus']['start']) for row in pending}
        _require(len(loci) <= 8)
        for chromosome, position in sorted(loci):
            _require(type(chromosome) is str and chromosome in {*(str(n) for n in range(1, 23)), 'X', 'Y', 'MT'}
                     and type(position) is int and 1 <= position < 2**29)
            connection = session._index.connection
            projected = []
            with budget.sqlite(connection):
                cursor = connection.execute('''SELECT length(r.canonical_row), r.canonical_row
                    FROM records r JOIN (
                      SELECT sequence FROM exact_records WHERE source='clinvar' AND chromosome=? AND position=?
                      UNION ALL
                      SELECT sequence FROM unknown_loci WHERE source='clinvar' AND chromosome=? AND start=?
                    ) p USING(sequence) LIMIT 65''', (chromosome, position, chromosome, position))
                for size, content in cursor:
                    budget.charge(size, 'reference_site_projection')
                    projected.append(json.loads(content))
            _require(0 < len(projected) <= 64)
            projected.sort(key=lambda row: row['source_locator']['record_number'])
            raw_rows = self.tabix.query(paths[0], chromosome, position, position,
                                        max_records=65, budget=budget, preserve_duplicates=True)
            def fetch(chrom, start, length):
                budget.charge(length, 'reference_site_reference')
                return session._fasta.fetch('chrM' if chrom == 'MT' else 'chr' + chrom, start, length)
            resolved.update(reconcile_locus(projected, raw_rows, artifact=artifact,
                chromosome=chromosome, position=position, fetch=fetch, query=query))
            for row in pending:
                if row['affected_locus']['chromosome'] == chromosome and row['affected_locus']['start'] == position:
                    _require(row in projected)
        _require(tuple(_signature(path) for path in paths) == self.files)
        session._index.assert_stable()
        if any(row['source_locator']['record_number'] not in resolved for row in pending):
            return None
        return tuple(resolved[row['source_locator']['record_number']]
                     if row['record_type'] == 'unknown' else row for row in records)
