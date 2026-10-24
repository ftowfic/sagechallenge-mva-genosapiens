"""Local interpretation of an already closed diagnosis-guided source selection.

Only the source adapter and graph restriction are Track 2 specific. Scientific
association, QC, annotation, population and exception semantics stay with their
existing owners. Regional membership alone never becomes an association proof.
"""
from collections import defaultdict
from dataclasses import dataclass, replace
from types import SimpleNamespace

from genosapiens import candidate_admission as admission, annotation_admission
from genosapiens import admission_links as links, population_annotation, population_policy
from genosapiens.allele_identity import ReferenceRole, SequenceAllele
from genosapiens.evidence_model import (
    AlleleEvidence, CallEvidence, DiseaseGeneEvidence, IdentityMapping,
    IdentityStage, ScoringView, SourcedRecord,
)
from genosapiens.evidence_serialization import _plain
from genosapiens.regulatory_index import indexed_regulatory_report
from genosapiens.track2_diagnosis_scope import resolve_diagnosis
from genosapiens.track2_published_genotype import canonical, sha256
from genosapiens.track2_vcf_selection import SourceSelection, sequence_quality
from genosapiens.track2_vcf_evidence import source_numbers
from genosapiens.track2_annotation_artifacts import FileAnnotationInput
from genosapiens.variant_normalization import VARIANT_NORMALIZATION_SCHEMA
from genosapiens.vcf_qc import FilterStatus, technical_flags_for


class ScopedInterpretationError(ValueError):
    """Static messages; no case-derived fields in operator diagnostics."""


def require(ok, message='scoped interpretation input is missing or inconsistent'):
    if not ok:
        raise ScopedInterpretationError(message)


def selected_alleles(selection):
    """Project actual QC/normalization records; merge equivalent observations.

    Original ALT indices and source calls survive. Duplicate records are not
    independent copies; the existing genotype kernel owns copy reconciliation.
    """
    require(type(selection) is SourceSelection)
    bound = selection.execution['binding']
    wanted = {r['source_allele_id']: r for r in selection.selected_sources}
    # Annotation describes supplied calls. Scope and retention changes must not
    # alter these calls' identities when the selected requests remain identical.
    # The interpretation and genotype bindings still bind the full selection.
    projection_sha256=sha256(canonical({key:bound[key] for key in (
        'source_sha256','sample','reference_sha256','reference_identity_id',
        'reference_role','concordance_sha256','vcf_policy_sha256','pid_contract_sha256')}))
    require(len(wanted) == len(selection.selection['selected_alleles']))
    groups, seen = defaultdict(list), set()
    for record, normalized in selection.records:
        sample = next((s for s in record.sample_calls if s.sample_name == bound['sample']), None)
        require(sample is not None and normalized.source_record_id == record.source_record_id)
        for item in normalized.alleles:
            if item.source_allele_id not in wanted:
                continue
            n = item.normalized_allele
            require(n is not None and n.normalized_allele_id == wanted[item.source_allele_id]['normalized_allele_id']
                    and item.source_allele_id not in seen)
            seen.add(item.source_allele_id)
            flags = technical_flags_for(record, sample, item.allele_number)
            source = SourcedRecord('genosapiens.track2-scoped-vcf', '1', item.source_allele_id,
                bound['source_sha256'], None, {
                    'source_record_id': record.source_record_id, 'source_allele_id': item.source_allele_id,
                    'source_input_sha256': bound['source_sha256'], 'record_number': record.record_number,
                    'allele_number': item.allele_number, 'raw_line_sha256': record.raw_line_sha256,
                    'filter_value': record.filter_value, 'technical_flags': flags,
                    'automatic_quality': sequence_quality(record, sample, item.allele_number),
                    'quality_alternate_count': len(record.alternates),
                    'read_validation_required': bool(flags), 'quality': record.quality,
                    'numeric_evidence':source_numbers(record),
                    'source_projection_sha256': projection_sha256,
                    'normalization': _plain(item), 'reference_concordance_sha256': bound['concordance_sha256']})
            original = SequenceAllele(ReferenceRole(n.reference_role), n.reference_identity_id,
                record.contig, record.position, record.reference, item.alternate)
            identity = IdentityMapping(item.source_allele_id, IdentityStage.SOURCE, original,
                bound['reference_sha256'], source)
            call = CallEvidence(source, sample,
                None if record.filter_status is FilterStatus.NOT_APPLIED else record.filters, record.quality)
            groups[n.normalized_allele_id].append((n, identity, call))
    require(seen == set(wanted), 'selected source allele accounting is incomplete')
    result = []
    for key, rows in sorted(groups.items()):
        n = rows[0][0]
        require(all(row[0] == n for row in rows))
        sequence = SequenceAllele(ReferenceRole(n.reference_role), n.reference_identity_id,
            n.contig, n.position, n.reference, n.alternate)
        originals = tuple(row[1] for row in rows)
        provenance = SourcedRecord('genosapiens.track2-scoped-normalization', '2', key,
            sha256(canonical(_plain(rows))), None, {'source_projection_sha256':projection_sha256})
        identities = (*originals, IdentityMapping(key, IdentityStage.NORMALIZED, sequence,
            bound['reference_sha256'], provenance, tuple(i.identity_id for i in originals), VARIANT_NORMALIZATION_SCHEMA),
            IdentityMapping('native-'+key, IdentityStage.NATIVE, sequence, bound['reference_sha256'],
                            provenance, (key,), 'genosapiens.track2-scoped-sequence/v1'))
        result.append(AlleleEvidence(key, 'native-'+key, identities, tuple(row[2] for row in rows)))
    return tuple(result)


@dataclass(frozen=True)
class ScopedSources:
    annotation: admission.AnnotationInput | FileAnnotationInput | None = None
    regulatory_index: object | None = None
    reference: links.PinnedPublicSource | None = None
    clinvar: links.PinnedPublicSource | None = None
    open_targets: links.PinnedPublicSource | None = None
    population: object | None = None
    clinical_session: object | None = None


@dataclass(frozen=True)
class ScopedInterpretation:
    alleles: tuple
    contexts: tuple
    decisions: tuple
    annotation_proofs: tuple
    link_proofs: tuple
    population_assessments: dict
    binding: dict
    limitations: tuple
    population_policy_json: str | None = None
    retention: dict | None = None
    deferred: tuple = ()


def _association(proof, adapter):
    return admission._Association(proof.allele_id, proof.rule.replace('_','-'), adapter,
        {'admitted':'match', 'evaluated_no_match':'no_match', 'unresolved':'unresolved'}[proof.outcome],
        proof.target_namespace, proof.target_id, getattr(proof, 'condition_id', None),
        proof.reasons, (*proof.source_refs, proof.provenance))


def _scope_graph(graph, concepts, genes):
    """Restrict existing context alternatives, retaining their original ledgers."""
    keys = {k for k,(c,namespace) in graph.candidates.items()
            if c.native_disease_id in concepts and namespace == 'Ensembl'
            and c.gene_id and c.gene_id.split('.')[0] in genes}
    changes = {name:{k:v for k,v in getattr(graph,name).items() if k in keys}
        for name in ('candidates','ledgers','supported','phenotype_state','seed_contexts')}
    changes.update({name:{k:tuple(v for v in values if v in keys)
                         for k,values in getattr(graph,name).items()}
                    for name in ('candidate_by_target','candidate_by_condition')})
    # HPO never supplies admission or ordering in this origin. It is handled as
    # separate optional review context after genotype construction.
    changes['pheno'] = ()
    return replace(graph, **changes)


def select_scoped_retention(facts, *, scope, catalog_sha256):
    """Select already scoped, source-bound facts before genotype construction.

    The source adapter owns diagnosis matching and complete exact-allele intake.
    This boundary applies no second consequence list or phenotype predicate.
    """
    from genosapiens.variant_retention import RetentionFacts, select_variant
    require(type(facts) is RetentionFacts and isinstance(scope, dict)
        and scope.get('schema') == 'genosapiens.track2-diagnosis-scope/v1'
        and isinstance(scope.get('genes'), (list, tuple)))
    return select_variant(facts, established_genes=frozenset(scope['genes']), catalog_sha256=catalog_sha256)


def interpret_selection(selection, *, scope, store, sources=ScopedSources()):
    """Consume local source owners, not caller-supplied scientific verdicts."""
    require(type(selection) is SourceSelection and type(sources) is ScopedSources)
    require(selection.selection['scope_sha256'] == sha256(canonical(scope))
        and scope['resource_bindings']['disease_gene'] == store.evidence_sha256
        and scope['resource_bindings']['disease_gene_sources'] == sha256(canonical(_plain(store.source_binding))))
    resolved, concepts, _ = resolve_diagnosis(store, scope['original_diagnosis'])
    require(resolved == scope['resolved_diagnosis'] and set(scope['mapping_evidence']['concept_ids']) == concepts)
    alleles = selected_alleles(selection)
    base = {a.allele_id:a for a in alleles}
    annotations, proofs, associations = (), [], []
    resource_bindings = {}
    if sources.reference is not None:
        require(sources.reference.manifest_sha256 == scope['reference_sha256']
            and sources.reference.contract_sha256 == scope['resource_bindings']['reference_contract'])
    if sources.annotation is not None:
        a = sources.annotation
        require(tuple(r.allele for r in a.requests) == alleles,
                'annotation requests differ from selected source alleles')
        validate = (annotation_admission.validate_annotation_admission_files
            if type(a) is FileAnnotationInput else annotation_admission.validate_annotation_admission)
        intake = validate(a.artifacts,
            expected_binding_sha256=a.expected_binding_sha256, requests=a.requests)
        annotations, alleles = intake.proofs(), intake.alleles
        require({p.allele_id for p in annotations} == set(base))
        associations.extend(_association(p, 'annotation') for p in annotations)
        resource_bindings['annotation'] = a.expected_binding_sha256
    if sources.regulatory_index is not None:
        require(sources.reference is not None
            and sources.regulatory_index.manifest_sha256 == scope['resource_bindings']['regulatory'])
        rows = tuple(indexed_regulatory_report(alleles, index=sources.regulatory_index,
            reference=sources.reference, selected_targets={a.allele_id:tuple(scope['genes']) for a in alleles}))
        require({p.allele_id for p in rows} == set(base))
        associations.extend(_association(p, 'encode') for p in rows); proofs.extend(rows)
        resource_bindings['regulatory'] = sources.regulatory_index.manifest_sha256
    if sources.clinical_session is not None:
        from genosapiens.track2_clinical_batches import clinical_batches
        session = sources.clinical_session
        require(session.reference == sources.reference)
        enriched = []
        for association_batch, annotation_batch in clinical_batches(session, alleles):
            for item in association_batch.proofs:
                associations.append(_association(item.proof, item.adapter))
                proofs.append(item.proof)
            enriched.extend(annotation_batch.alleles)
        alleles = tuple(enriched)
        resource_bindings['retention_index'] = session.indexes.manifest_sha256
    population = sources.population
    prepared_policy = None
    if population is not None:
        from genosapiens.genotype_hypotheses import HypothesisPopulationInputs
        require(type(population) is HypothesisPopulationInputs)
        prepared_policy = population_policy.prepare_population_policy(population.policy_bytes,
            expected_sha256=population.policy_sha256)
        population_policy.validate_population_curation(prepared_policy.policy, population.public_curation_bytes)
        requests = tuple(population_annotation.PopulationRequest(a,a.allele_id) for a in alleles)
        annotated = population_annotation.run_population_annotation(requests, joint=population.joint,
            all_sites=population.all_sites, reference=population.reference,
            runtime_image=population.runtime_image, limits=population.limits)
        alleles = annotated.alleles
        resource_bindings['population_policy'] = population.policy_sha256
        resource_bindings['population_query'] = sha256(canonical(_plain(annotated.report)))
    from genosapiens.track2_retention import select_alleles
    alleles, retention, deferred = select_alleles(alleles, scope=scope,
        annotation_sha256=resource_bindings.get('annotation'), clinical_session=sources.clinical_session)
    retained_ids = {a.allele_id for a in alleles}
    associations = [p for p in associations if p.subject in retained_ids]
    annotations = tuple(p for p in annotations if p.allele_id in retained_ids)
    proofs = [p for p in proofs if p.allele_id in retained_ids]
    # An exact curated assertion can link a distal allele to its source gene.
    # Clinical assertions without a gene remain unassigned.
    for decision in retention['retained']:
        for fact in decision['facts']['functional']:
            if fact['gene_id'] not in scope['genes']:
                continue
            source = SourcedRecord(fact['source']['dataset'], fact['extraction_policy'],
                fact['source']['record_id'], fact['source']['release_sha256'],
                fact['source']['release_sha256'], fact)
            from genosapiens.track2_retention import scoped_disease_ids
            for condition in sorted(scoped_disease_ids(fact['disease_ids']) & concepts):
                associations.append(admission._Association(decision['facts']['identity_id'],
                    'explicit-clinical-link', 'curated', 'match', 'Ensembl', fact['gene_id'],
                    condition, ('exact_curated_disease_allele',), (source,)))
    emitted = {(p.subject,p.adapter) for p in associations}
    for allele in alleles:
        for adapter,rule in (('annotation','sequence-transcript-overlap'),('encode','regulatory-link'),
                             ('open_targets','explicit-clinical-link'),('clinvar','explicit-clinical-link')):
            if (allele.allele_id,adapter) not in emitted:
                associations.append(admission._Association(allele.allele_id,rule,adapter,'unresolved',
                    None,None,None,('source-unavailable',)))
    from genosapiens.disease_gene_evidence import track2_requirement_projection, TRACK2_ALLELE_POLICY
    projected_store=track2_requirement_projection(store)
    genes = (DiseaseGeneEvidence('diagnosis-scope-genes', store.source_binding['bundle_manifest_sha256'], projected_store),)
    prepared = admission._prepare_admission_phenotype(ScoringView((),(),{}), None, genes, None)
    graph = admission._prepare_admission_graph(prepared,
        reached_targets=admission.association_targets(tuple(associations)))
    graph = _scope_graph(graph, concepts, set(scope['genes']))
    binding = {'scope_sha256':sha256(canonical(scope)), 'selection_sha256':sha256(canonical(selection.selection)),
        'allele_requirement_policy':TRACK2_ALLELE_POLICY,'requirement_projection_sha256':projected_store.evidence_sha256,
        'resources':resource_bindings, 'admission_policy_sha256':admission.POLICY_SHA256,
        'interpretation':'shared-kernels-fixed-diagnosis/v2', 'phenotype':'context_only',
        'retention_sha256':sha256(canonical(retention))}
    inventory = SimpleNamespace(alleles=alleles,source_alleles=())
    contexts, decisions, _, limits = admission._evaluate_graph(inventory,tuple(associations),prepared,
        graph,binding,include_phenotype=False)
    assessments = {}
    if prepared_policy is not None:
        by_id = {a.allele_id:a for a in alleles}
        for seed in contexts:
            assessments[seed.context.context_id] = population_policy.evaluate_population_members(
                population_policy.ContextKey.from_context(seed.context),
                tuple(by_id[r.allele_id] for r in seed.allele_refs), prepared_policy=prepared_policy)
    return ScopedInterpretation(alleles,contexts,decisions,annotations,tuple(proofs),assessments,binding,
        tuple(sorted(set((*limits, 'SV/CNV/repeat/symbolic interpretation unsupported',
            'Outside-scope partners unassessed', 'Regional inclusion is not causal evidence',
            *(('Population and exception assessment unavailable',) if population is None else ()))))),
        None if population is None else population.policy_bytes.decode('utf-8'), retention, deferred)
