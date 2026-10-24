"""Shared genotype kernels over the complete frozen diagnosis-scoped domain."""
from dataclasses import dataclass
from types import SimpleNamespace

from genosapiens import genotype_properties as properties, genotype_hypotheses as hypotheses
from genosapiens import population_policy as frequency
from genosapiens.evidence_model import DiseaseGeneEvidence, SourcedRecord
from genosapiens.evidence_serialization import _plain
from genosapiens.partner_relevance import ContextRelevanceAssessment, SourceContextRelevance
from genosapiens.track2_published_genotype import canonical, sha256
from genosapiens.track2_scoped_interpretation import ScopedInterpretation, require, selected_alleles
from genosapiens.variant_normalization import SourceAllele


@dataclass(frozen=True)
class ScopedGenotypes:
    properties: properties.SuppliedGenotypeProperties
    searches: tuple
    hypotheses: tuple
    binding: dict
    population_assessments: tuple = ()
    requirement_policy: object | None = None


def _inventory(selection, alleles):
    mapped={i.identity_id:a.allele_id for a in alleles for i in a.identities if i.stage.value=='source'}
    entries,records=[],[]
    for record,normal in selection.records:
        originals=[]
        for item in normal.alleles:
            if item.source_allele_id not in mapped:
                continue
            source=SourceAllele(item.source_allele_id,normal.source_record_id,normal.input_sha256,
                               normal.record,item.allele_number,item.alternate)
            entries.append(SimpleNamespace(source=source,derived_allele_ids=(mapped[item.source_allele_id],),status='derived'))
            originals.append(source.source_allele_id)
        if originals:
            source=SourcedRecord('genosapiens.track2-scoped-vcf','1',record.source_record_id,
                normal.input_sha256,None,{'raw_line_sha256':record.raw_line_sha256})
            records.append(SimpleNamespace(qc=record,source=source))
    require({entry.source.source_allele_id for entry in entries}==set(mapped))
    return SimpleNamespace(records=tuple(records),source_alleles=tuple(entries),alleles=alleles,
        source_set=None,input_sha256=selection.selection['source_sha256'],
        reference_manifest_sha256=selection.selection['reference_sha256'])


def _relevance(selection, interpreted, seed):
    context=seed.context
    binding={'schema':'genosapiens.track2-scoped-context-relevance/v1',
        'interpretation_sha256':sha256(canonical(_plain(interpreted))),
        'context_id':context.context_id,'scope':'complete_retained_domain; deferred_and_outside_scope_unknown'}
    # Existing shared gene relevance requires a positive target association;
    # absent/other target links cannot prove exclusion from a gene context.
    matching={d.subject_id for d in interpreted.decisions if d.association_state=='match'
        and d.target_namespace=='Ensembl' and d.target_id
        and d.target_id.split('.')[0]==context.gene_id.split('.')[0]
        and (d.condition_id is None or context.context_id in d.context_ids)}
    entries=[]
    for allele in interpreted.alleles:
        for call in allele.calls or ():
            fields=call.source.source_fields
            relevant=allele.allele_id in matching
            entries.append(SourceContextRelevance(fields['source_record_id'],fields['source_allele_id'],
                (allele.allele_id,),'derived',(),'relevant' if relevant else 'unresolved_relevance',(),
                (call.source,),context.evidence_refs,(),(),(),
                ('source_target_association_context_only_not_qualification' if relevant
                 else 'no_source_proved_context_target_or_exclusion',),sha256(canonical(binding))))
    return ContextRelevanceAssessment(context,tuple(sorted(entries,key=lambda e:e.source_allele_id)),binding),tuple(sorted(matching))


def construct_scoped_genotypes(selection, interpreted, *, store, policy, maximum_pairs=100000):
    """Close all partitions before pairing; no top-N or shortened success path."""
    require(type(interpreted) is ScopedInterpretation and type(policy) is hypotheses.HypothesisPolicy
        and policy.maximum_pairs_per_context is None, 'scoped genotype construction cannot truncate partners')
    require(interpreted.binding['selection_sha256']==sha256(canonical(selection.selection)))
    require(type(maximum_pairs) is int and 0<=maximum_pairs<=1000000,'invalid genotype capacity bound')
    original={a.allele_id:(a.identities,a.calls) for a in selected_alleles(selection)}
    require(interpreted.retention is not None
        and interpreted.binding['retention_sha256'] == sha256(canonical(interpreted.retention)),
        'genotype construction requires completed variant retention')
    retained = {r['facts']['identity_id'] for r in interpreted.retention['retained']}
    require(retained <= set(original)
        and {a.allele_id:(a.identities,a.calls) for a in interpreted.alleles}
            == {k:v for k,v in original.items() if k in retained},
            'genotype source domain is incomplete or changed')
    inventory=_inventory(selection,interpreted.alleles)
    digest=sha256(canonical({'selection':selection.selection,'alleles':_plain(interpreted.alleles)}))
    sample=selection.selection['sample'];by_id={a.allele_id:a for a in interpreted.alleles}
    contexts,all_pairs={},set()
    for seed in interpreted.contexts:
        relevance,keys=_relevance(selection,interpreted,seed)
        require(not hypotheses._requires_pairs(seed.context) or len(keys)*(len(keys)-1)//2<=maximum_pairs,
                'genotype construction incomplete: complete partner domain exceeds capacity')
        universe=hypotheses._pairs(keys,by_id,None) if hypotheses._requires_pairs(seed.context) else None
        if universe:all_pairs.update(universe.pairs)
        contexts[seed.context.context_id]=(relevance,keys,universe)
    genotypes=properties._evaluate(inventory,sample,tuple(sorted(all_pairs)),inventory_sha256=digest)
    from genosapiens.disease_gene_evidence import track2_requirement_projection
    projected_store=track2_requirement_projection(store)
    require(projected_store.evidence_sha256==interpreted.binding['requirement_projection_sha256'],
            'Allele requirement derivation differs from interpretation.')
    gene_evidence=(DiseaseGeneEvidence('diagnosis-scope-genes',store.source_binding['bundle_manifest_sha256'],projected_store),)
    policy_sha=sha256(hypotheses.hypothesis_policy_bytes(policy,
        admission_policy_sha256=interpreted.binding['admission_policy_sha256']))
    searches,results=[],[]
    for seed in interpreted.contexts:
        relevance,keys,universe=contexts[seed.context.context_id]
        search,generated=hypotheses.evaluate_complete_context(context=seed.context,relevance=relevance,keys=keys,
            alleles=by_id,genotypes=genotypes,policy=policy,disease_gene_evidence=gene_evidence,sample_name=sample,
            observation_ids=(),input_sha256=selection.selection['source_sha256'],inventory_sha256=digest,
            policy_sha256=policy_sha,pair_universe=universe,source_inventory_status='assessed-domain-complete',
            outside_context_partner_uncertainty=True)
        searches.append(search);results.extend(generated)
    results=tuple(sorted(results,key=lambda h:h.hypothesis.hypothesis_id))
    population=()
    if interpreted.population_policy_json is not None:
        prepared=frequency.prepare_population_policy(interpreted.population_policy_json.encode('utf-8'),
            expected_sha256=interpreted.binding['resources']['population_policy'])
        population=tuple(frequency.evaluate_population_policy(item.hypothesis,
            tuple(by_id[key] for key in item.hypothesis.evidence_subject_ids),prepared_policy=prepared)
            for item in results)
    return ScopedGenotypes(genotypes,tuple(searches),results,
        {'scope_sha256':selection.selection['scope_sha256'],'selection_sha256':sha256(canonical(selection.selection)),
         'interpretation_sha256':sha256(canonical(_plain(interpreted))),'policy_sha256':policy_sha,
         'genotype_policy_sha256':properties.POLICY_SHA256,'inventory_sha256':digest,
         'population_assessments_sha256':sha256(canonical(_plain(population))),
         'outside_scope_partners':'unknown','phenotype_ordering':'unused'},population,policy.requirements)
