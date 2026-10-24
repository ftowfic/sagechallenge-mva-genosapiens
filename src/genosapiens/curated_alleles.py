"""Bounded public narrative projection with explicit unresolved extraction gaps.

This adapter accepts a narrow, literal expression-loss statement next to one
genomic substitution. Other wording remains a curation gap. It never labels
every allele in a gene damaging and never requests case-derived evidence.
"""
from __future__ import annotations

from dataclasses import asdict, dataclass
import re

from genosapiens.common_allele_public_indexes import normalize_public_alias
from genosapiens.variant_retention import ExactAllele, FunctionalFact, SourceFact, RetentionError


EXTRACTION_POLICY = 'curated-exact-expression-substitution/1.0.0'
MAX_NARRATIVE_BYTES = 1024 * 1024
_HGVS = re.compile(r'(?<![\w.])(NC_[0-9]+\.[0-9]+):g\.([1-9][0-9]*)([ACGT])>([ACGT])(?![A-Za-z0-9])')
_NEGATION = re.compile(r'\b(?:not|no|neither|unconfirmed|uncertain|unresolved|unlikely|possibly|might|may)\b', re.I)


@dataclass(frozen=True, slots=True)
class AccessionReference:
    accession: str
    assembly: str
    reference_identity: str
    reference_sha256: str
    contig: str


@dataclass(frozen=True, slots=True)
class CuratedGap:
    source: SourceFact
    reason: str
    original_hgvs: str | None
    passage: str
    extraction_policy: str = EXTRACTION_POLICY


@dataclass(frozen=True, slots=True)
class CuratedProjection:
    assertions: tuple[FunctionalFact, ...]
    gaps: tuple[CuratedGap, ...]


def project_curated_narrative(row, *, source_sha256, accession_references, fetch_reference):
    """Project one frozen public row. The reference owner authenticates aliases.

    All source publication identifiers remain together. Their presence does not
    establish which publication supports each sentence. The extraction limit
    records this distinction for downstream scientific assessment.
    """
    if (type(row) is not dict or row.get('source_dataset') not in
            {'gencc_submissions', 'evidence_clingen'}):
        raise RetentionError('curated allele source is unsupported')
    source = SourceFact(row['source_dataset'], source_sha256, row['source_record_id'])
    notes = row.get('notes') or ''
    if type(notes) is not str or len(notes.encode('utf-8')) > MAX_NARRATIVE_BYTES:
        raise RetentionError('curated narrative is invalid or exceeds its bound')
    references = {}
    for reference in accession_references:
        if type(reference) is not AccessionReference or reference.accession in references:
            raise RetentionError('curated accession mapping is ambiguous')
        references[reference.accession] = reference
    hits = tuple(_HGVS.finditer(notes))
    if not hits:
        return CuratedProjection((), (CuratedGap(source, 'no_supported_exact_genomic_substitution', None, notes),))
    boundaries = (0, *(match.end() for match in re.finditer(r'(?<=[.!?])\s+(?=[A-Z])', notes)), len(notes))
    assertions, gaps = [], []
    for hit in hits:
        # Accessions, g. notation, and "et al. (year)" do not end sentences.
        index = next(i for i in range(len(boundaries)-1)
                     if boundaries[i] <= hit.start() < boundaries[i+1])
        start, end = boundaries[index:index+2]
        next_end = boundaries[min(index+2, len(boundaries)-1)]
        sentence = notes[start:end]
        effect_sentence = notes[end:next_end].strip()
        passage = notes[start:next_end]
        def gap(reason):
            gaps.append(CuratedGap(source, reason, hit.group(), passage))
        reference = references.get(hit[1])
        if reference is None:
            gap('accession_reference_unresolved')
            continue
        assemblies = re.findall(r'\bGRCh[0-9]+(?:\.p[0-9]+)?\b', sentence)
        if assemblies and set(assemblies) != {reference.assembly}:
            gap('source_assembly_conflicts_with_accession')
            continue
        if len(_HGVS.findall(passage)) != 1:
            gap('multiple_alleles_in_effect_passage')
            continue
        gene, symbol, disease = row.get('ensembl_gene_id'), row.get('gene_symbol'), row.get('source_disease_id')
        publications = row.get('publications')
        if (not gene or not symbol or not disease or row.get('is_ontology_propagated') is not False
                or row.get('association_scope') != 'source_direct'):
            gap('gene_or_disease_applicability_unresolved')
            continue
        if (not isinstance(publications, (list, tuple)) or not publications
                or any(type(p) is not str or not p.isdecimal() for p in publications)):
            gap('publication_identity_unresolved')
            continue
        # Match the complete source statement. An arbitrary substring about
        # reduced expression cannot establish an assertion or resolve negation.
        expected = (r'This mutation is a hypomorphic mutation that down regulates? the expression of '
                    + re.escape(symbol) + r' gene\.')
        if _NEGATION.search(passage) or re.fullmatch(expected, effect_sentence) is None:
            gap('effect_or_negation_unresolved')
            continue
        try:
            alias = normalize_public_alias(reference.contig, int(hit[2]), hit[3], hit[4], fetch_reference)
            allele = ExactAllele(reference.reference_identity, reference.reference_sha256,
                alias.chromosome, alias.position, alias.reference, alias.alternate)
        except Exception:
            # A failed reference read is operational. Do not convert it into a
            # scientific no-match or a successfully completed empty index.
            raise RetentionError('curated reference mapping failed') from None
        assertions.append(FunctionalFact(source, allele, 'exact', gene, (disease,), 'damaging',
            False, '|'.join('PMID:' + p for p in publications), passage, EXTRACTION_POLICY,
            ('literal_expression_loss_statement', 'source_publications_not_sentence_specific',
             'no_genotype_or_mechanism_acceptance')))
    return CuratedProjection(tuple(assertions), tuple(gaps))


def accession_references_from_report(content, *, reference_identity, reference_sha256):
    """Read exact primary RefSeq accessions from an authenticated assembly report."""
    headers = dict((m[1].strip(), m[2].strip()) for line in content.splitlines()
                   if (m := re.fullmatch(r'# ([^:]+):\s*(.*)', line)))
    if (headers.get('Assembly name') != 'GRCh38' or headers.get('Taxid') != '9606'
            or headers.get('GenBank assembly accession') != 'GCA_000001405.15'):
        raise RetentionError('curated accession source assembly is invalid')
    result, seen = [], set()
    for line in content.splitlines():
        if not line or line.startswith('#'):
            continue
        row = line.split('\t')
        if len(row) != 10:
            raise RetentionError('curated accession source row is malformed')
        if row[1] != 'assembled-molecule' or row[2] not in {*map(str, range(1, 23)), 'X', 'Y', 'MT'}:
            continue
        if row[5] != '=' or re.fullmatch(r'NC_[0-9]+\.[0-9]+', row[6]) is None or row[6] in seen:
            raise RetentionError('curated accession source mapping is ambiguous')
        seen.add(row[6])
        result.append(AccessionReference(row[6], 'GRCh38', reference_identity, reference_sha256, row[2]))
    if not result:
        raise RetentionError('curated accession source mapping is empty')
    return tuple(result)


def project_curated_index(writer, rows, *, source_sha256, accession_references, fetch_reference):
    """Extend the existing public projection once, independently of case alleles."""
    from genosapiens.common_allele_public_indexes import _exact, _unknown, SequenceAlias
    counts = dict(source_rows=0, assertions=0, gaps=0)
    for row in rows:
        if row.get('source_dataset') not in {'gencc_submissions', 'evidence_clingen'}:
            continue
        counts['source_rows'] += 1
        projection = project_curated_narrative(row, source_sha256=source_sha256,
            accession_references=accession_references, fetch_reference=fetch_reference)
        for fact in projection.assertions:
            key = fact.allele
            alias = SequenceAlias(key.contig, key.position, key.reference, key.alternate)
            writer.emit(_exact('curated_functional', key.contig, key.position, key.reference,
                key.alternate, alias, asdict(fact.source), {'functional': asdict(fact)}))
            counts['assertions'] += 1
        for gap in projection.gaps:
            # Global source-curation gaps are inventoried without assigning an
            # unknown assertion to every allele in the source gene.
            writer.emit({**_unknown('curated_curation_gap', gap.reason, asdict(gap.source)),
                         'curation_gap': asdict(gap)})
            counts['gaps'] += 1
    return counts


def functional_from_index(record):
    """Strict current-format decode of an authenticated exact functional row."""
    try:
        if record['source'] != 'curated_functional' or record['record_type'] != 'exact':
            raise ValueError
        value = dict(record['evidence']['functional'])
        value['source'] = SourceFact(**value['source'])
        value['allele'] = ExactAllele(**value['allele'])
        for name in ('disease_ids', 'extraction_limits', 'conflicts'):
            value[name] = tuple(value[name])
        fact = FunctionalFact(**value)
        key = record['normalized_alias']
        if (fact.allele.contig, fact.allele.position, fact.allele.reference, fact.allele.alternate) != (
                key['chromosome'], key['position'], key['reference'], key['alternate']):
            raise ValueError
        if fact.extraction_policy != EXTRACTION_POLICY:
            raise ValueError
        return fact
    except (KeyError, TypeError, ValueError):
        raise RetentionError('curated functional index row is malformed or stale') from None
