"""Compose the minimal package from the actual scoped scientific producers."""
import copy
import json
from dataclasses import asdict

from genosapiens.evidence_serialization import _plain
from genosapiens.track2_diagnosis_package import scenario_identity, validate_components
from genosapiens.track2_mechanism import INPUT_FIELDS, build_mechanism_model
from genosapiens.track2_phenotype import OptionalPhenotype, ingest_optional_phenotype
from genosapiens.track2_published_genotype import canonical, sha256
from genosapiens.track2_scoped_interpretation import require
from genosapiens.track2_scoped_mechanisms import _evidence


MODEL_SCHEMA='genosapiens.track2-scoped-mechanism/v1'


def construct_components(*, scope, selection, interpreted, genotypes, mechanisms, phenotype=None):
    """Preserve every selected source finding and all constructed alternatives.

    A genotype uses one deterministic representative source ID per normalized
    allele. All duplicate observations remain in variants and original evidence;
    their explicit map never promotes extra genotype copies or independent study
    support. Scientific genotype identity remains the shared kernel's identity.
    """
    phenotype=ingest_optional_phenotype() if phenotype is None else phenotype
    require(type(phenotype) is OptionalPhenotype
        and selection.selection['scope_sha256']==sha256(canonical(scope))
        and genotypes.binding['interpretation_sha256']==sha256(canonical(_plain(interpreted))))
    by_id={a.allele_id:a for a in interpreted.alleles}
    originals={a.allele_id:tuple(sorted(i.identity_id for i in a.identities if i.stage.value=='source'))
               for a in interpreted.alleles}
    representative={key:values[0] for key,values in originals.items()}
    source_to_normal={source:key for key,values in originals.items() for source in values}
    require(set(source_to_normal)<={r['source_allele_id'] for r in selection.selection['selected_alleles']})
    require(interpreted.retention is not None
        and set(source_to_normal) == {occurrence for row in interpreted.retention['retained']
                                     for occurrence in row['facts']['occurrence_ids']})
    evidence={row['evidence_id']:row for row in (_evidence(a) for a in interpreted.alleles)}
    allele_evidence={key:_evidence(a)['evidence_id'] for key,a in by_id.items()}
    def retain(row):
        previous=evidence.setdefault(row['evidence_id'],row)
        require(previous==row,'evidence identity has conflicting source content')
    interpretation_raw=canonical({**{key:_plain(getattr(interpreted,key)) for key in (
        'contexts','decisions','annotation_proofs','link_proofs','population_assessments','binding','limitations')},
        'genotype_population_assessments':_plain(genotypes.population_assessments)})
    interpretation_id='scoped-interpretation:'+sha256(interpretation_raw)[:24]
    retain({'evidence_id':interpretation_id,'source_identity':sha256(interpretation_raw),
        'lineage_key':'interpretation:'+selection.selection['subset_sha256'],
        'assertion':interpretation_raw.decode('ascii'),'evidence_role':'source_observation'})
    variants=[]
    admitted={r.allele_id for seed in interpreted.contexts for r in seed.allele_refs}
    for row in selection.selection['selected_alleles']:
        if row['source_allele_id'] not in source_to_normal:
            continue
        normalized=source_to_normal[row['source_allele_id']]
        genes=sorted({d.target_id.split('.')[0] for d in interpreted.decisions if d.subject_id==normalized
            and d.association_state=='match' and d.target_namespace=='Ensembl' and d.target_id
            and d.target_id.split('.')[0] in scope['genes']})
        variants.append({**row,'gene_ids':genes,'assessment':'admitted' if normalized in admitted else 'uncertain',
                         'evidence_ids':[allele_evidence[normalized],interpretation_id]})
    hypotheses={h.hypothesis.hypothesis_id:h for h in genotypes.hypotheses}
    scenarios=[]
    for mechanism in mechanisms:
        require(mechanism.genotype_id in hypotheses)
        genotype=hypotheses[mechanism.genotype_id]
        require(set(mechanism.model.joint_model.allele_ids)==set(genotype.hypothesis.allele_ids))
        for row in mechanism.evidence:retain(row)
        statement={'schema':MODEL_SCHEMA,'mechanism':_plain(mechanism),'genotype':_plain(genotype),
            'normalized_to_package_allele_ids':{key:representative[key] for key in genotype.hypothesis.allele_ids},
            'construction_binding':genotypes.binding}
        raw=canonical(statement);proof_id='scoped-model:'+sha256(raw)[:24]
        retain({'evidence_id':proof_id,'source_identity':sha256(raw),'lineage_key':mechanism.genotype_id,
            'assertion':raw.decode('ascii'),'evidence_role':'source_observation'})
        facts=json.loads(canonical({k:v for k,v in asdict(mechanism.readiness).items() if k in INPUT_FIELDS}))
        for allele in facts['allele_evidence']:allele['allele_id']=representative[allele['allele_id']]
        scenario={'scenario_id':'pending','disease_id':genotype.hypothesis.context.native_disease_id,
            'gene_id':genotype.hypothesis.context.gene_id.split('.')[0],
            'genotype_id':mechanism.genotype_id,'allele_ids':[representative[key] for key in genotype.hypothesis.allele_ids],
            'evidence_ids':sorted({proof_id,interpretation_id,*(e['evidence_id'] for e in mechanism.evidence)}),
            'mechanism_inputs':facts,'limitations':list(mechanism.limitations)}
        scenario['scenario_id']=scenario_identity(scenario);facts['scenario_id']=scenario['scenario_id']
        scenarios.append(scenario)
    retained_selection = copy.deepcopy(selection.selection)
    retained_selection['selected_alleles'] = [row for row in retained_selection['selected_alleles']
        if row['source_allele_id'] in source_to_normal]
    retained_selection['subset_sha256'] = sha256(canonical(retained_selection['selected_alleles']))
    retained_selection['coverage']['selected_allele_count'] = len(retained_selection['selected_alleles'])
    retained_selection['coverage']['retention'] = copy.deepcopy(interpreted.retention)
    from genosapiens.track2_deferred import describe
    retained_selection['coverage']['deferred_inventory'] = describe(
        interpreted.deferred, retained_selection, interpreted.retention)
    parts={'scope':copy.deepcopy(scope),'selection':retained_selection,'variants':variants,
        'structural':copy.deepcopy(selection.structural),
        'scenarios':sorted(scenarios,key=lambda r:r['scenario_id']),'phenotype':copy.deepcopy(phenotype.component),
        'evidence':[evidence[k] for k in sorted(evidence)]}
    validate_components(parts,classification=selection.selection['classification'],
                        evidence_origin=selection.selection['evidence_origin'])
    return parts


def scenario_model(components, scenario):
    """Read the source model bound by this exact package scenario's evidence."""
    selected=set(scenario['evidence_ids']);proofs=[]
    for row in components['evidence']:
        if row['evidence_id'] not in selected or not row['evidence_id'].startswith('scoped-model:'):
            continue
        require(sha256(row['assertion'].encode())==row['source_identity'])
        proof=json.loads(row['assertion']);require(proof['schema']==MODEL_SCHEMA)
        require(proof['mechanism']['genotype_id']==scenario['genotype_id'])
        proofs.append(proof)
    require(len(proofs)==1,'scenario requires exactly one source-bound mechanism model')
    proof=proofs[0];model=copy.deepcopy(proof['mechanism']['model'])
    mapping=proof['normalized_to_package_allele_ids']
    require(set(mapping.values())==set(scenario['allele_ids']))
    model['scenario_id']=scenario['scenario_id']
    for allele in model['allele_models']:allele['allele_id']=mapping[allele['allele_id']]
    model['joint_model']['allele_ids']=[mapping[key] for key in model['joint_model']['allele_ids']]
    return build_mechanism_model(model)
