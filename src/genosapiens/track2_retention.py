"""Shared exact public evidence in the frozen diagnosis domain."""
from dataclasses import asdict
import re

from genosapiens.common_allele_public_indexes import retention_disease_ids
from genosapiens.curated_alleles import functional_from_index
from genosapiens.retention_sources import clinical_from_index
from genosapiens.track2_published_genotype import canonical, sha256
from genosapiens.variant_retention import LookupCoverage, RetentionFacts, select_variant


def scoped_disease_ids(values):
    """Normalize identifier notation only, without name or mapped-view inference."""
    return {('ORPHA' if namespace == 'Orphanet' else namespace) + ':' + number
        for value in values for namespace, number in re.findall(
            r'(?<![A-Za-z0-9_])(MONDO|OMIM|ORPHA|Orphanet)[:_]([0-9]+)(?![0-9A-Za-z])', value)}


def select_alleles(alleles, *, scope, annotation_sha256, clinical_session=None):
    """Select once before admission graphs and genotype expansion."""
    from genosapiens.variant_retention import SourceFact, TranscriptFact, POLICY_SHA256, POLICY_VERSION
    from genosapiens.track2_scoped_interpretation import select_scoped_retention
    concepts = set(scope['mapping_evidence']['concept_ids'])
    catalog_sha = scope['resource_bindings']['disease_gene']
    if alleles and annotation_sha256 is None:
        raise ValueError('Complete initial annotation is required before variant retention.')
    if alleles and clinical_session is None:
        raise ValueError('Complete indexed public lookups are required before variant retention.')
    public = []
    if clinical_session is not None:
        for offset in range(0, len(alleles), clinical_session.execution_config.public_queries.batch_rows):
            public.extend(clinical_session.retention_batch(alleles[
                offset:offset + clinical_session.execution_config.public_queries.batch_rows]))
    retained, decisions, deferred = [], [], []
    for allele, lookup in zip(alleles, public, strict=True):
        coverage = list(lookup.coverage)
        coverage.append(LookupCoverage('Ensembl VEP', annotation_sha256, True, 'complete'))
        transcripts = []
        for tx in allele.transcript_consequences or ():
            raw = tx.source.source_fields.get('csq', {})
            from genosapiens.annotation_admission import _gene
            gene = _gene(raw.get('Gene', ''))
            # Scope uses stable Ensembl IDs. Keep the versioned source in the
            # transcript evidence, but join its validated stable ID here.
            gene = gene.split('.')[0] if gene else None
            transcript = raw.get('Feature') if raw.get('Feature_type') == 'Transcript' else None
            transcripts.append(TranscriptFact(
                SourceFact('Ensembl VEP', annotation_sha256, tx.consequence_id), gene,
                transcript, tuple(v for v in tx.consequences or () if v),
                'relevant' if gene and transcript else 'unresolved',
                sift=raw.get('SIFT'), polyphen=raw.get('PolyPhen'), impact=raw.get('IMPACT')))
        occurrences = tuple(sorted(i.identity_id for i in allele.identities if i.stage.value == 'source'))
        facts = RetentionFacts(allele.allele_id, lookup.allele, 'assess', occurrences, tuple(transcripts),
            tuple(f for f in lookup.clinical if scoped_disease_ids(f.disease_ids) & concepts),
            tuple(f for f in lookup.functional if scoped_disease_ids(f.disease_ids) & concepts),
            (), tuple(coverage))
        decision = select_scoped_retention(facts, scope=scope, catalog_sha256=catalog_sha)
        if decision.retained:
            retained.append(allele)
            decisions.append(asdict(decision))
        else:
            native = next(i.identity for i in allele.identities if i.identity_id == allele.native_identity_id)
            from genosapiens.evidence_serialization import _plain
            deferred.append({'identity_id': allele.allele_id, 'identity': _plain(native),
                'source_occurrences': list(occurrences),
                'annotation_summary': sorted({c for tx in transcripts for c in tx.consequences}),
                'reason': decision.disposition, 'policy_version': POLICY_VERSION,
                'policy_sha256': POLICY_SHA256,
                'source_versions': [asdict(c) for c in coverage]})
    summary = {'policy_version': POLICY_VERSION, 'policy_sha256': POLICY_SHA256,
        'catalog_sha256': catalog_sha, 'retained': decisions,
        'counts': {'unit': 'normalized_sequence_alleles', 'scoped': len(alleles),
                   'retained': len(retained), 'deferred': len(deferred)},
        'occurrences': {'unit': 'source_ALT_occurrences',
            'retained': sum(len(d['facts']['occurrence_ids']) for d in decisions),
            'deferred': sum(len(d['source_occurrences']) for d in deferred)},
        'lookup_metrics': dict(clinical_session.retention_metrics) if clinical_session is not None else {},
        'limitations': ['Deferred alleles did not enter genotype, mechanism, or medication analysis.',
            'Deferral does not establish benignity or the absence of another disease allele.',
            'Retention does not establish pathogenicity, mechanism readiness, or medication evidence.']}
    from genosapiens.track2_codec import json_value
    return tuple(retained), json_value(summary), tuple(json_value(deferred))


def validate_summary(value, selection, scope):
    """Replay the retained domain before any package can enter assessment."""
    from genosapiens.variant_retention import RetentionDecision, POLICY_SHA256, POLICY_VERSION
    from genosapiens.track2_codec import decode
    from genosapiens.track2_scoped_interpretation import require
    require(type(value) is dict and set(value) == {'policy_version','policy_sha256','catalog_sha256',
        'retained','counts','occurrences','lookup_metrics','limitations'})
    require(value['policy_version'] == POLICY_VERSION and value['policy_sha256'] == POLICY_SHA256
        and value['catalog_sha256'] == scope['resource_bindings']['disease_gene'])
    counts, occurrences = value['counts'], value['occurrences']
    require(set(counts) == {'unit','scoped','retained','deferred'}
        and set(occurrences) == {'unit','retained','deferred'}
        and counts['unit'] == 'normalized_sequence_alleles'
        and occurrences['unit'] == 'source_ALT_occurrences'
        and all(type(counts[k]) is int and counts[k] >= 0 for k in ('scoped','retained','deferred'))
        and all(type(occurrences[k]) is int and occurrences[k] >= 0 for k in ('retained','deferred'))
        and counts['scoped'] == counts['retained'] + counts['deferred']
        and counts['retained'] == len(value['retained'])
        and occurrences['retained'] == len(selection['selected_alleles']))
    identities, source_ids = set(), set()
    for row in value['retained']:
        decision = decode(RetentionDecision, row)
        replay = select_variant(decision.facts, established_genes=frozenset(scope['genes']),
                                catalog_sha256=value['catalog_sha256'])
        require(decision == replay and decision.retained
            and decision.facts.identity_id not in identities
            and not source_ids.intersection(decision.facts.occurrence_ids))
        identities.add(decision.facts.identity_id)
        source_ids.update(decision.facts.occurrence_ids)
    require(source_ids == {row['source_allele_id'] for row in selection['selected_alleles']})


def known_scope(index, *, concepts, genes, reference_identity, reference_sha256, catalog_sha256):
    """Find exact disease alleles before restrictive sequence retrieval.

    The database returns source assertions. Only the shared predicate decides
    whether an assertion supplies a retention route. Gene-free clinical records
    remain gene-free and cannot create a genotype association.
    """
    coverage = index.retention_coverage()
    bindings = coverage['source_bindings']
    if bindings['reference'] != reference_sha256:
        raise ValueError('Known allele index uses a different reference.')
    checked = tuple(LookupCoverage(dataset, bindings[key], True, 'complete')
        for dataset, key in (('ClinVar VCF', 'clinvar'), ('Open Targets evidence_eva', 'open_targets'),
                            ('gencc_submissions', 'curated'), ('evidence_clingen', 'curated')))
    # Resolution retains all validated exact cross-references. The allele index
    # supports a narrower identifier vocabulary. Project only its exact keys,
    # without dropping mappings from diagnosis evidence or parsing ID substrings.
    query_concepts = {value for value in concepts
                      if re.fullmatch(r'(MONDO|OMIM|ORPHA):[0-9]+', value)}
    if not query_concepts:
        raise ValueError('Known allele scope has no supported disease identifier.')
    rows = index.retention_for_diseases(tuple(sorted(query_concepts)))
    assertions = []
    for row in rows:
        if not set(retention_disease_ids(row)) & query_concepts:
            raise ValueError('Known allele query returned an unrelated disease.')
        clinical, functional = (), ()
        if row['source'] in {'clinvar', 'open_targets'}:
            fact = clinical_from_index(row, reference_identity=reference_identity,
                reference_sha256=reference_sha256, source_release_sha256=bindings[row['source']])
            clinical = (fact,)
        elif row['source'] == 'curated_functional':
            fact = functional_from_index(row)
            if (fact.allele.reference_identity != reference_identity
                    or fact.allele.reference_sha256 != reference_sha256
                    or fact.source.release_sha256 != bindings['curated']):
                raise ValueError('Known functional allele has different source or reference bindings.')
            functional = (fact,)
        else:
            raise ValueError('Known allele query returned an unsupported source.')
        identity = sha256(canonical(asdict(fact.allele)))
        facts = RetentionFacts(identity, fact.allele, 'assess', (fact.source.record_id,),
            (), clinical, functional, (), checked)
        decision = select_variant(facts, established_genes=frozenset(genes), catalog_sha256=catalog_sha256)
        if decision.retained:
            assertions.append({'allele': asdict(fact.allele), 'assertion': asdict(fact),
                'routes': list(decision.routes), 'disease_ids': list(retention_disease_ids(row))})
    return {'index_sha256': index.verified.manifest_sha256,
        'disease_lookup_sha256': getattr(getattr(index, 'disease_lookup', None), 'manifest_sha256', None),
        'source_coverage': coverage, 'query_disease_ids': sorted(query_concepts),
        'unqueried_disease_ids': sorted(set(concepts) - query_concepts),
        'assertions': sorted(assertions, key=canonical),
        'limitation': 'Known allele lookup uses only explicit MONDO, OMIM, and ORPHA identifiers. '
            'Other mapped identifiers remain in diagnosis evidence. Unmapped disease names remain unresolved.'}
