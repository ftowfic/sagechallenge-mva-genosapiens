"""Published-genotype admission to the shared native discovery engine."""
from dataclasses import asdict
from datetime import datetime
import json
from pathlib import Path

from genosapiens.track2_operator_inputs import require
from genosapiens.track2_published_genotype import canonical, sha256, load_published_genotype
from genosapiens.track2_scientific_storage import json_identity, publish_json

SCHEMA = 'genosapiens.track2-public-discovery-plan/v1'


def reviewed_mechanism(bundle, plan_text, review_text):
    """Reuse actual review of the same public genotype and mechanism only."""
    from genosapiens.track2_input import source_mechanism_inputs
    from genosapiens.track2_mechanism import assess_mechanism_readiness
    plan, review = json.loads(plan_text), json.loads(review_text)
    require(review['schema']=='genosapiens.track2-plan-review/v1'
        and review['decision']=='accepted' and review['reviewer_type']=='human'
        and review['independence'] in {'independent','self_review_disclosed'}
        and review['reviewer'] and review['rationale']
        and datetime.fromisoformat(review['reviewed_at'].replace('Z','+00:00')).tzinfo is not None
        and review['plan_sha256']==sha256(plan_text.encode()), 'The retained public mechanism review is invalid.')
    require(plan['classification']=='public' and plan['input_origin']=='published_patient_genotype'
        and all(review['source_hashes'].get(k)==v for k,v in bundle.source_hashes.items()),
        'The retained mechanism review belongs to different public sources.')
    origin=asdict(bundle.admitted_origin())
    records=plan['tables']['mechanism_records']
    require(len(records)==1 and records[0]['context']['origin']==origin,
        'The retained mechanism belongs to a different reviewed genotype.')
    source_mechanism_inputs(bundle,plan['mechanism_inputs'])
    facts=plan['mechanism_inputs'][origin['scenario_id']]
    decision=assess_mechanism_readiness(facts)
    current=json.loads(canonical(asdict(decision)))
    require('readiness_state' in records[0]['readiness']
        and all(k in current and current[k]==v for k,v in records[0]['readiness'].items()),
        'The retained mechanism decision differs from its scientific inputs.')
    evidence={e['evidence_id'] for e in bundle.data['evidence']}
    require(all(set(a['evidence_ids'])<=evidence for a in facts['allele_evidence']),
        'The public mechanism requires evidence from its reviewed genotype.')
    return facts


def public_assessment(bundle, facts, mapping, descriptor, policy):
    from genosapiens.track2_exploratory_planning import discovery_eligibility
    from genosapiens.track2_mechanism import assess_mechanism_readiness
    symbols={a['gene'] for a in bundle.data['genotype']['alleles']}
    require(symbols=={mapping['gene_symbol']} and mapping['gene_id'].startswith('ENSG'),
        'The discovery gene mapping differs from the published genotype.')
    scenario=dict(scenario_id=bundle.data['scenario_id'],gene_id=mapping['gene_id'],
        mechanism_inputs=facts,evidence_ids=sorted({e for a in facts['allele_evidence'] for e in a['evidence_ids']}))
    scenario['discovery_eligibility']=discovery_eligibility(scenario,assess_mechanism_readiness(facts),policy)
    return dict(schema='genosapiens.track2-public-mechanism-assessment/v1',classification='public',
        construction_sha256=bundle.input_sha256,construction_review_sha256=bundle.review_sha256,
        scenarios=[scenario],evidence=[],automatic_inputs=dict(public_packet=None,
            discovery=dict(resources_sha256=sha256(canonical(descriptor)),policy_sha256=sha256(canonical(policy)))))


def checked_public_provenance(plan, *, bundle=None, allow_synthetic=False):
    """Replay retained public provenance without claiming a new human review."""
    from genosapiens.track2_published_genotype import PublishedGenotype, validate_genotype
    from genosapiens.track2_exploratory_planning import discovery_eligibility
    require(set(plan)=={'schema','classification','input_origin','public_provenance',
        'mechanism_assessment','exploratory_discovery','methods'} and plan['schema']==SCHEMA
        and plan['classification']=='public' and plan['input_origin']=='published_patient_genotype',
        'The public discovery plan contract is invalid.')
    proof=plan['public_provenance']
    require(set(proof)=={'genotype_json','input_review_json','mechanism_plan_json','mechanism_review_json',
        'source_hashes','gene_mapping','policy'},'The public provenance contract is invalid.')
    data=json.loads(proof['genotype_json']);review=json.loads(proof['input_review_json'])
    input_hash=sha256(proof['genotype_json'].encode());review_hash=sha256(proof['input_review_json'].encode())
    validate_genotype(data,None,retained_source_hashes=proof['source_hashes'],allow_synthetic=allow_synthetic)
    require(data['classification']=='public' and review['schema']=='genosapiens.published-genotype-review/v1'
        and review['input_sha256']==input_hash and review['source_hashes']==proof['source_hashes']
        and review['decision']=='accepted' and review['reviewer_type']=='human'
        and review['independence'] in {'independent','self_review_disclosed'} and review['reviewer']
        and review['rationale'] and review['limitations']
        and set(review['evidence_ids'])=={e['evidence_id'] for e in data['evidence']}
        and datetime.fromisoformat(review['reviewed_at'].replace('Z','+00:00')).tzinfo is not None,
        'The retained public input review is invalid.')
    if bundle is not None:
        require(bundle.input_sha256==input_hash and bundle.review_sha256==review_hash
            and bundle.source_hashes==proof['source_hashes'],'The public source input changed.')
    else:
        bundle=PublishedGenotype(data,input_hash,proof['source_hashes'],review,review_hash,True)
    facts=reviewed_mechanism(bundle,proof['mechanism_plan_json'],proof['mechanism_review_json'])
    policy=proof['policy']
    from genosapiens.track2_exploratory_planning import installed_policy
    require(policy==installed_policy(sha256(canonical(policy))),'The public discovery policy changed.')
    # The nested plan retains the descriptor hash. The full resource checker
    # authenticates the descriptor and exact gene mapping before execution.
    assessment=plan['mechanism_assessment']
    mapping=proof['gene_mapping']
    scenario=dict(scenario_id=data['scenario_id'],gene_id=mapping['gene_id'],mechanism_inputs=facts,
        evidence_ids=sorted({e for a in facts['allele_evidence'] for e in a['evidence_ids']}))
    from genosapiens.track2_mechanism import assess_mechanism_readiness
    scenario['discovery_eligibility']=discovery_eligibility(scenario,assess_mechanism_readiness(facts),policy)
    expected=dict(schema='genosapiens.track2-public-mechanism-assessment/v1',classification='public',
        construction_sha256=input_hash,construction_review_sha256=review_hash,scenarios=[scenario],evidence=[],
        automatic_inputs=dict(public_packet=None,discovery=dict(
            resources_sha256=plan['exploratory_discovery']['resources']['sha256'],
            policy_sha256=sha256(canonical(policy)))))
    require({a['gene'] for a in data['genotype']['alleles']}=={mapping['gene_symbol']}
        and json_identity(assessment)==json_identity(expected),'The public mechanism provenance differs.')
    return bundle


def check_public_plan(plan, resources, *, bundle=None):
    from genosapiens.track2_exploratory_planning import validate_exploratory_plan
    bundle=checked_public_provenance(plan,bundle=bundle)
    mapping=plan['public_provenance']['gene_mapping']
    rows=resources['connection'].execute('SELECT gene_id FROM genes WHERE approved_symbol=?',
        (mapping['gene_symbol'],)).fetchall()
    require([r[0] for r in rows]==[mapping['gene_id']], 'The public gene mapping is ambiguous or changed.')
    validate_exploratory_plan(plan['exploratory_discovery'],plan['mechanism_assessment'],resources)
    return bundle


def prepare_public_plan(bundle_root, reviewed_plan, descriptor, reference_root, output, runtime_images):
    from genosapiens.track2_discovery_resources import load_discovery_resources
    from genosapiens.track2_exploratory_planning import installed_policy, _prepare
    bundle_root,reviewed_plan,output=Path(bundle_root),Path(reviewed_plan),Path(output)
    bundle=load_published_genotype(bundle_root)
    plan_text=reviewed_plan.read_text();review_text=reviewed_plan.with_name('plan-review.json').read_text()
    facts=reviewed_mechanism(bundle,plan_text,review_text)
    policy=installed_policy();connection=load_discovery_resources(descriptor,reference_root,matrix_required=False)
    resources=dict(connection=connection,descriptor=descriptor,policy=policy,
        descriptor_path='discovery-resources.json',policy_path='exploratory-policy.json',runtime_images=runtime_images)
    try:
        symbol=bundle.data['genotype']['alleles'][0]['gene']
        rows=connection.execute('SELECT gene_id FROM genes WHERE approved_symbol=?',(symbol,)).fetchall()
        require(len(rows)==1,'The public gene needs one exact catalog mapping.')
        mapping=dict(gene_symbol=symbol,gene_id=rows[0][0])
        assessment=public_assessment(bundle,facts,mapping,descriptor,policy)
        plan=dict(schema=SCHEMA,classification='public',input_origin='published_patient_genotype',
            public_provenance=dict(genotype_json=(bundle_root/'genotype.json').read_text(),
                input_review_json=(bundle_root/'review.json').read_text(),mechanism_plan_json=plan_text,
                mechanism_review_json=review_text,source_hashes=bundle.source_hashes,gene_mapping=mapping,policy=policy),
            mechanism_assessment=assessment,exploratory_discovery=_prepare(assessment,resources,policy),
            methods=dict(scope='Reviewed published genotype with current native discovery across all eligible cells.',
                review='Earlier input and mechanism reviews retain their original scope. This run has no new human output review.',
                limits='The publication lacks a versioned transcript mapping. No source VCF or new genomic coordinates are inferred.'))
        check_public_plan(plan,resources,bundle=bundle)
    finally:connection.close()
    output.mkdir(parents=True,mode=0o700)
    for name,data in [('plan.json',plan),('discovery-resources.json',descriptor),('exploratory-policy.json',policy)]:
        publish_json(output/name,data)
    return plan


def source_section(plan):
    proof=plan['public_provenance'];data=json.loads(proof['genotype_json'])
    lines=['## Published genotype and scope','',
        'This analysis uses reviewed published observations. Challenge-case matching remains unassessed.',
        'The input and mechanism reviews retain their original scope. This reconstruction has no new human output review.','',
        '| Gene | Reported coding change | Reported protein change | Copies |','|---|---|---|---|']
    for a in data['genotype']['alleles']:
        lines.append('| '+' | '.join(str(a[k]) for k in ('gene','reported_coding','reported_protein','copies'))+' |')
    lines+=['','The genotype remains publication-scoped. This run does not infer a genomic coordinate or a versioned transcript.',
        'The current catalog defines the medication and cell coverage. Earlier candidate counts do not determine this run.',
        'Functional evidence retains its source tissue and assay limits. Expression agreement does not establish protein restoration.','',
        'Source titles and allele identifiers retain their exact forms for traceability.','']
    lines += ['- ['+s['citation']+']('+s['public_url']+')' for s in data['sources']]
    return '\n'.join(lines)+'\n\n'
