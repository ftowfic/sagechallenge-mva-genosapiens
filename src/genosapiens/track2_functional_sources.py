"""Frozen, gene-independent functional lookup with conservative protein matching."""
from dataclasses import replace
import gzip
import io
import json
from pathlib import Path
import re
import sqlite3

from genosapiens.evidence_model import PublicSourcedEvidence, SourcedRecord
from genosapiens.evidence_serialization import _plain
from genosapiens.track2_operator_inputs import require, small_json
from genosapiens.track2_published_genotype import canonical, sha256

SCHEMA = 'genosapiens.track2-uniprot-functional-catalog/v1'
DATASET = 'UniProt reviewed human natural variants'
QUERY_SCHEMA = 'genosapiens.track2-uniprot-gene-query/v1'
POLICY = 'track2-uniprot-protein-applicability/2.0.0'
PIN = 'config/track2/functional-sources.json'


def installed_catalog(repository, reference):
    path = Path(repository)/PIN
    if not path.exists():
        return None
    pin = small_json(path)
    require(set(pin) == {'schema','manifest_path','manifest_sha256'}
        and pin['schema'] == SCHEMA, 'The functional source declaration is invalid.')
    root = Path(reference).resolve(strict=True)
    manifest = root/pin['manifest_path']
    require(manifest.is_relative_to(root) and manifest.resolve() == manifest.absolute(),
        'The functional source location is invalid.')
    descriptor = dict(manifest_path=str(manifest), manifest_sha256=pin['manifest_sha256'])
    data = catalog_manifest(descriptor)
    database = manifest.parent/data['database']
    require(database.is_file() and database.resolve() == database.absolute(),
        'The declared functional source index is missing or invalid.')
    return descriptor


def catalog_manifest(descriptor):
    path = Path(descriptor['manifest_path'])
    with path.open('rb') as stream:
        raw = stream.read(16*1024*1024+1)
    require(len(raw) <= 16*1024*1024 and sha256(raw) == descriptor['manifest_sha256'],
        'The functional source manifest changed.')
    data = json.loads(raw)
    require(data['schema'] == SCHEMA and data['classification'] == 'public'
        and data['status'] == 'complete' and data['records'] > 0 and data['pages']
        and sum(p['records'] for p in data['pages']) == data['records']
        and 'cc_alternative_products' in data['fields'].split(',')
        and data['database'] == 'catalog.sqlite', 'The functional source is incomplete.')
    return data


def gene_sources(descriptor, gene_ids):
    """Verify only requested gene payloads against the complete source manifest."""
    manifest = catalog_manifest(descriptor)
    from genosapiens.track2_protein_supplement import installed_supplement, supplement_entry, PIN_SHA256
    supplement = installed_supplement()
    database = Path(descriptor['manifest_path']).parent/manifest['database']
    require(database.resolve() == database.absolute(), 'The functional index path is invalid.')
    sources = []
    connection = sqlite3.connect(database.as_uri()+'?mode=ro', uri=True)
    try:
        for gene in sorted(set(gene_ids)):
            require(re.fullmatch(r'ENSG[0-9]{11}', gene) is not None,
                'Functional lookup requires an exact gene identifier.')
            expected = manifest['genes'].get(gene)
            row = connection.execute('SELECT content FROM genes WHERE gene=?', (gene,)).fetchone()
            require((row is None) == (expected is None), 'Functional gene lookup coverage changed.')
            if row is None:
                entries = []
            else:
                with gzip.GzipFile(fileobj=io.BytesIO(row[0])) as stream:
                    raw = stream.read(32*1024*1024+1)
                require(len(raw) <= 32*1024*1024 and sha256(raw) == expected['sha256'],
                    'The functional gene record changed.')
                value = json.loads(raw)
                require(value['gene_id'] == gene and len(value['entries']) == expected['entries'],
                    'The functional gene record has a different scope.')
                entries = [supplement_entry(entry, manifest['release'], supplement) for entry in value['entries']]
            fields = dict(schema=QUERY_SCHEMA, gene_id=gene, entries=entries,
                coverage='gene_present' if entries else 'gene_not_in_snapshot',
                catalog_release=manifest['release'], catalog_records=manifest['records'],
                source_gene_sha256=None if expected is None else expected['sha256'],
                protein_representation_supplement_sha256=PIN_SHA256,
                limitations=manifest['limitations'])
            key = 'uniprot-gene:'+gene
            record = SourcedRecord(DATASET, manifest['release'], key, sha256(canonical(fields)),
                descriptor['manifest_sha256'], fields)
            sources.append(PublicSourcedEvidence(key, descriptor['manifest_sha256'], {key:record}))
    finally:
        connection.close()
    return tuple(sources)


def bind_gene_sources(policy, descriptor, gene_ids):
    if descriptor is None:
        return policy
    added = gene_sources(descriptor, gene_ids)
    requirements = replace(policy.requirements,
        public_evidence=policy.requirements.public_evidence+added)
    return replace(policy, requirements=requirements)


def scoped_policy(policy, gene):
    """Do not duplicate other genes' annotations in every scenario proof."""
    return replace(policy, public_evidence=tuple(public for public in policy.public_evidence
        if not public.evidence_id.startswith('uniprot-gene:') or public.evidence_id == 'uniprot-gene:'+str(gene)))


def description_effect(description):
    """Interpret literal effect labels only, without mining pathogenicity labels."""
    from genosapiens.track2_functional_observations import observations
    measured, issue = observations(description)
    if issue:
        return None, issue
    if measured:
        categories = {row['category'] for row in measured}
        return (next(iter(categories)) if len(categories) == 1 else 'compatible_functional_observations'), None
    text = description.lower().replace('–','-')
    labels = {'loss_of_function':r'loss[ -]of[ -]function',
        'gain_of_function':r'gain[ -]of[ -]function', 'dominant_negative':r'dominant[ -]negative',
        'hypomorphic':r'hypomorphic'}
    found = [name for name, pattern in labels.items() if re.search(r'\b'+pattern+r'\b', text)]
    if found:
        # A label with uncertain scope or negation cannot become an assertion.
        if (len(found) != 1 or re.search(r'\b(?:not|no|without|may|might|possibly|putative|potential|uncertain)\b', text)
                or re.search(r'\b(?:when|together|combination|same allele|cis|trans)\b', text)):
            return None, 'effect_scope_or_negation_unresolved'
        label = labels[found[0]]
        if not any(re.fullmatch(r'(?:(?:causes|results in) )?(?:a )?'+label+
                r'(?: mutation| variant| effect)?\.?', clause.strip()) for clause in text.split(';')):
            return None, 'functional_assertion_requires_interpretation'
        return found[0], None
    # Retain unparsed assay statements as potential evidence, never as no-match.
    if re.search(r'\b(?:function|activity|expression|binding|affinity|localization|localisation|stability|'
                 r'checkpoint|activation|inhibition|abundance|transport|signaling|signalling|folding|'
                 r'catalytic|splicing|rescue|complementation|channel|current|protein levels)\b', text):
        return None, 'functional_assertion_requires_interpretation'
    return None, None


def _check_gene(record, context, allele, transcripts):
    from genosapiens.model_requirements import _transcripts
    fields = _plain(record.source_fields)
    require(record.source_dataset == DATASET and fields['schema'] == QUERY_SCHEMA
        and sha256(canonical(fields)) == record.source_sha256
        and fields['catalog_release'] == record.source_release
        and fields['gene_id'] == context.gene_id,
        'The retained functional query changed.')
    source = dict(evidence_id=record.source_record_id, manifest_sha256=record.public_manifest_sha256,
        release=record.source_release, coverage=fields['coverage'], entries=len(fields['entries']))
    matches, gaps = [], []
    if not fields['entries']:
        return dict(state='complete', sources=[source], matches=[], gaps=[], contrary=False,
            rule=POLICY, coverage=fields['coverage'], limitations=fields['limitations'])
    selected, errors = _transcripts(context, allele, transcripts)
    if errors:
        gaps.append('relevant_transcript_unavailable')
    for tx, raw, transcript, terms in selected:
        linked = []
        for entry in fields['entries']:
            for crossref in entry.get('uniProtKBCrossReferences', []):
                if crossref['database'] != 'Ensembl':
                    continue
                props = {p['key']:p['value'] for p in crossref.get('properties', [])}
                if props.get('GeneId','').split('.')[0] != context.gene_id:
                    continue
                if crossref['id'] == transcript:
                    linked.append((entry, crossref, props))
        if not linked:
            gaps.append('versioned_transcript_mapping_unavailable')
            continue
        for entry, crossref, props in linked:
            displayed = {identifier for comment in entry.get('comments', [])
                if comment['commentType'] == 'ALTERNATIVE PRODUCTS'
                for isoform in comment.get('isoforms', []) if isoform.get('isoformSequenceStatus') == 'Displayed'
                for identifier in isoform['isoformIds']}
            isoform = crossref.get('isoformId')
            if isoform and isoform not in displayed:
                gaps.append('noncanonical_protein_sequence_unavailable')
                continue
            if (tx.ensembl_protein_id != props.get('ProteinId')
                    or entry['primaryAccession'] not in (tx.uniprot_accessions or ())):
                gaps.append('versioned_protein_mapping_unavailable')
                continue
            position = raw.get('Protein_position') or ''
            aa = (raw.get('Amino_acids') or '').split('/')
            if (not re.fullmatch(r'[1-9][0-9]*(?:-[1-9][0-9]*)?', position)
                    or len(aa) != 2 or not re.fullmatch('[ACDEFGHIKLMNPQRSTVWY]+', aa[0])):
                gaps.append('protein_alteration_unresolved')
                continue
            bounds = list(map(int, position.split('-')))
            start, end = bounds[0], bounds[-1]
            sequence=entry.get('sequence',{}).get('value')
            if not isinstance(sequence,str):
                gaps.append('protein_reference_sequence_unavailable')
                continue
            if sequence[start-1:end] != aa[0]:
                gaps.append('protein_reference_residue_differs')
                continue
            for feature in entry.get('features', []):
                if feature['type'] != 'Natural variant':
                    continue
                effect, issue = description_effect(feature.get('description',''))
                if effect is None and issue is None:
                    continue
                location = feature['location']
                if location.get('sequence') not in (None, isoform):
                    continue
                left, right = location['start'].get('value'), location['end'].get('value')
                if left is None or right is None:
                    gaps.append('functional_variant_position_unresolved')
                    continue
                if end < left or start > right:
                    continue
                alternative = feature.get('alternativeSequence', {})
                representation = entry.get('track2_representation_proof') or {}
                if representation.get('state') == 'unresolved':
                    gaps.append(representation['reason'])
                from genosapiens.track2_protein_loss import terminal_loss_record, compatible_stop
                alteration = terminal_loss_record(entry, feature)
                compatibility = compatible_stop(alteration, start=start, end=end, reference=aa[0],
                    alternate=aa[1], terms=terms, sequence=sequence)
                if alteration and alteration['state'] == 'unresolved':
                    gaps.append(alteration['reason'])
                if not compatibility and (left != start or right != end or location['start'].get('modifier') != 'EXACT'
                        or location['end'].get('modifier') != 'EXACT'
                        or not alternative.get('originalSequence') or not alternative.get('alternativeSequences')):
                    gaps.append('overlapping_functional_variant_identity_unresolved')
                    continue
                if not compatibility and (alternative['originalSequence'] != aa[0]
                        or aa[1] not in alternative['alternativeSequences']):
                    continue
                if set(terms or ()) & {'splice_acceptor_variant','splice_donor_variant','splice_region_variant','frameshift_variant'}:
                    gaps.append('nucleotide_or_frameshift_effect_not_assessed_by_protein_source')
                    continue
                citations = [e['id'] for e in feature.get('evidences', [])
                    if e.get('source') == 'PubMed' and e.get('evidenceCode') == 'ECO:0000269']
                if not citations:
                    issue = 'experimental_source_citation_unavailable'
                if issue:
                    gaps.append(issue)
                from genosapiens.track2_functional_observations import interpretation_record
                interpretation = interpretation_record(feature.get('description', ''), citations)
                matches.append(dict(accession=entry['primaryAccession'], record_id=feature.get('featureId'),
                    transcript_id=transcript, protein_id=props['ProteinId'], value=effect,
                    interpretation_state='unresolved' if issue else 'agent_interpretation' if interpretation else 'source_label',
                    interpretation=interpretation,
                    citations=citations, source_sha256=record.source_sha256,
                    description=feature.get('description',''), original_feature=feature,
                    applicability='protein_effect_compatibility' if compatibility else 'exact_protein_alteration_only',
                    identity_proof=compatibility, nucleotide_identity='not_established',
                    representation_proof=entry.get('track2_representation_proof'),
                    patient_assay_transfer='not_established', disease_applicability='not_transferred'))
    from genosapiens.track2_functional_observations import COMPATIBLE_WITH_LOSS
    return dict(state='incomplete' if gaps else 'complete', sources=[source], matches=matches,
        gaps=sorted(set(gaps)), contrary=any(m['interpretation_state'] in {'source_label', 'agent_interpretation'}
            and m['value'] not in COMPATIBLE_WITH_LOSS for m in matches), rule=POLICY,
        coverage=fields['coverage'], limitations=fields['limitations'])


def functional_check(context, allele, transcripts, policy):
    from genosapiens.model_requirements import frozen_functional_check
    native, other = [], []
    for public in policy.public_evidence:
        if public.evidence_id.startswith('uniprot-gene:'):
            if public.evidence_id == 'uniprot-gene:'+str(context.gene_id):
                require(len(public.records) == 1, 'The functional gene query is ambiguous.')
                native.append(_check_gene(next(iter(public.records.values())), context, allele, transcripts))
        else:
            other.append(public)
    legacy = frozen_functional_check(context, allele, transcripts, replace(policy, public_evidence=tuple(other)))
    results = native+([legacy] if legacy['state'] != 'unavailable' else [])
    if not results:
        return legacy
    return dict(state='complete' if all(r['state']=='complete' for r in results) else 'incomplete',
        sources=[s for r in results for s in r['sources']], matches=[m for r in results for m in r['matches']],
        contrary=any(r['contrary'] for r in results), checks=results, rule=POLICY)
