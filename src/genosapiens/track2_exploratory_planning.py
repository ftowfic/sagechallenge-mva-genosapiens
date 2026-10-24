"""Study-independent hypotheses over unchanged mechanism acceptance."""
from dataclasses import asdict
import json
from pathlib import Path

from genosapiens.track2_mechanism import assess_mechanism_readiness
from genosapiens.track2_operator_inputs import require
from genosapiens.track2_published_genotype import canonical, sha256

SCHEMA = 'genosapiens.track2-exploratory-plan/v1'
ROUTES = ('qualified_target_drug_linkage', 'lincs_cell_specific_expression')


def installed_policy(policy_sha256=None):
    root = Path(__file__).resolve().parents[2]/'config/track2'
    for version in ('v2', 'v1'):
        policy=json.loads((root/('exploratory-policy-'+version+'.json')).read_bytes())
        if policy_sha256 is None or sha256(canonical(policy))==policy_sha256:
            return policy
    raise ValueError('The native discovery policy version is unsupported.')


def objective_id(scenario, policy):
    return 'exploratory-objective:'+sha256(canonical(dict(scenario_id=scenario['scenario_id'],
        gene_id=scenario['gene_id'],direction=policy['eligibility']['desired_modulation'],policy_version=policy['version'])))


def discovery_eligibility(scenario, mechanism_decision, policy):
    decision = asdict(mechanism_decision) if not isinstance(mechanism_decision,dict) else mechanism_decision
    accepted = decision['readiness_state'] in {'ready','prediction_accepted'}
    direction = decision['direction']==policy['eligibility']['mechanism']
    outcome = 'eligible' if accepted and direction else 'excluded' if accepted or decision['readiness_state']=='not_ready' else 'unresolved'
    reasons = (['accepted_loss_of_function_hypothesis'] if outcome=='eligible' else
        ['mechanism_outside_exploratory_policy'] if accepted else
        ['mechanism_'+decision['readiness_state'], *decision.get('conflicts',[])])
    return dict(outcome=outcome,reasons=sorted(set(reasons)),
        mechanism_decision_sha256=sha256(canonical(decision)),policy_sha256=sha256(canonical(policy)),
        objective_ids=[objective_id(scenario,policy)] if outcome=='eligible' else [])


def restoration_objective(scenario, gene, decision, policy):
    eligibility = discovery_eligibility(scenario,decision,policy)
    require(eligibility['outcome']=='eligible', 'Only accepted loss of function can create a native hypothesis.')
    data = asdict(decision) if not isinstance(decision,dict) else decision
    model = scenario.get('model') or {}
    alleles=model.get('allele_models',[])
    absent = bool(alleles) and all(row.get('residual_product')=='absent' for row in alleles)
    return dict(objective_id=objective_id(scenario,policy),scenario_id=scenario['scenario_id'],
        gene_id=scenario['gene_id'],gene_symbol=gene['gene_symbol'],kind='restoration_or_compensation_hypothesis',
        mechanism='loss_of_function',desired_modulation='increase',authority='machine',
        evidence_basis='policy_prediction' if data['readiness_state']=='prediction_accepted' else 'curated_mechanism',
        mechanism_decision_sha256=eligibility['mechanism_decision_sha256'],
        source_evidence_ids=sorted(set(scenario['evidence_ids'])),policy_sha256=eligibility['policy_sha256'],
        assumptions=['Restoration or compensation is a research hypothesis.'],
        limitations=['Patient tissue relevance remains unknown.', 'Rescue efficacy remains unknown.',
            'Residual protein availability remains unknown.' if not absent else 'Evidence reports absent gene product.',
            'Expression matching does not establish production of absent protein.'],
        known_conflicts=['target_product_absent'] if absent else [])


def plan_routes(objective, resources, policy):
    from genosapiens.track2_discovery_resources import gene_coverage
    coverage = gene_coverage(resources['connection'],objective['gene_id'],policy)
    executable = any(s['method'] and s['drug_signature_ids'] for s in coverage['selections'])
    routes = []
    for route in ROUTES:
        available = route==ROUTES[0] or executable
        routes.append(dict(objective_id=objective['objective_id'],route=route,
            planned_state='executable' if available else 'unavailable',
            reason_codes=[] if available else coverage['reason_codes'],
            request_id='discovery-request:'+sha256(canonical(dict(objective_id=objective['objective_id'],route=route,
                resources=resources['descriptor']['build_id'],policy_sha256=sha256(canonical(policy)))))))
    return routes, coverage


def _prepare(assessment, resources, policy):
    if policy['schema']=='genosapiens.track2-exploratory-policy/v2':
        from genosapiens.track2_tissue_planning import prepare_tissue_plan
        return prepare_tissue_plan(assessment, resources, policy)
    from genosapiens.track2_discovery_resources import resource_bindings
    from genosapiens.track2_lincs_selection import discovery_job_key
    descriptor = resources['descriptor']
    require(descriptor['policy_sha256']==sha256(canonical(policy)), 'The discovery policy differs from its catalog binding.')
    require(assessment['automatic_inputs']['discovery']==dict(resources_sha256=sha256(canonical(descriptor)),
        policy_sha256=sha256(canonical(policy))), 'The mechanism assessment binds different discovery inputs.')
    from genosapiens.track2_scientific_records import stored_json
    from genosapiens.track2_scientific_storage import Records
    def prospective_scope():
        for scenario in sorted(assessment['scenarios'],key=lambda r:r['scenario_id']):
            decision=assess_mechanism_readiness(scenario['mechanism_inputs'])
            eligibility=discovery_eligibility(scenario,decision,policy)
            require(scenario.get('discovery_eligibility')==eligibility,'Native eligibility differs from mechanism acceptance.')
            if eligibility['outcome']!='eligible':
                yield dict(excluded=dict(scenario_id=scenario['scenario_id'],gene_id=scenario['gene_id'],**eligibility))
                continue
            gene=resources['connection'].execute('SELECT approved_symbol FROM genes WHERE gene_id=?',(scenario['gene_id'],)).fetchone()
            objective=restoration_objective(scenario,dict(gene_symbol=gene[0] if gene else scenario['gene_id']),decision,policy)
            planned,observed=plan_routes(objective,resources,policy)
            yield dict(objective=objective,planned=planned,observed=observed)
    records=stored_json(Records(prospective_scope()))
    scope, routes, selections, scenario_jobs, excluded, coverage = [], [], {}, [], [], []
    for record in records:
        if 'excluded' in record:
            excluded.append(record['excluded']);continue
        objective,observed=record['objective'],record['observed']
        scope.append(objective);routes.extend(record['planned'])
        coverage.append(dict(objective_id=objective['objective_id'],**observed))
        for selection in observed['selections']:
            if not selection['method'] or not selection['drug_signature_ids']:continue
            job=discovery_job_key(selection,resource_bindings(descriptor),sha256(canonical(policy)))
            selections[job]=dict(job_id=job,**selection)
            scenario_jobs.append(dict(scenario_id=objective['scenario_id'],objective_id=objective['objective_id'],job_id=job))
    return dict(schema=SCHEMA,mode='native_exploratory',
        policy=dict(path=resources['policy_path'],sha256=sha256(canonical(policy))),
        resources=dict(path=resources['descriptor_path'],sha256=sha256(canonical(descriptor)),build_id=descriptor['build_id']),
        objectives=scope,routes=routes,genetic_selections=[selections[k] for k in sorted(selections)],
        scenario_job_bindings=scenario_jobs,excluded_scenarios=excluded,coverage=coverage,runtime_images=resources['runtime_images'])


def prepare_exploratory_plan(bundle, assessment, resources, policy):
    require(assessment['schema']=='genosapiens.track2-mechanism-assessment/v3'
        and assessment['construction_sha256']==bundle.input_sha256, 'Native discovery requires its bound mechanism assessment.')
    return _prepare(assessment,resources,policy)


def validate_exploratory_plan(plan, assessment, resources):
    from genosapiens.json_schema import validate_json_schema
    name='tissue-plan.schema.json' if plan.get('schema')=='genosapiens.track2-exploratory-plan/v2' else 'exploratory-plan.schema.json'
    validate_json_schema(plan,json.loads((Path(__file__).resolve().parents[2]/'schemas/track2'/name).read_bytes()))
    policy = resources.get('policy') or installed_policy()
    expected=_prepare(assessment,resources,policy)
    from genosapiens.track2_speedup_recovery import validate_prior_plan
    require(plan==expected or validate_prior_plan(plan,expected,resources), 'Native prospective scope differs from its bound inputs.')
    return plan


def prepare_native_plan(bundle, *, assessment, resources, methods, output, draft=False):
    from genosapiens.track2_diagnosis_plan import input_only_plan
    from genosapiens.track2_mechanism_review import scope_from_assessment
    from genosapiens.track2_operator_inputs import write_new
    from genosapiens.track2_scientific_storage import json_identity
    from genosapiens.track2_scientific_records import stored_json
    policy = resources.get('policy') or installed_policy()
    output = Path(output)
    output.mkdir(parents=True,exist_ok=True,mode=0o700)
    for name,data in [(resources['descriptor_path'],resources['descriptor']),(resources['policy_path'],policy)]:
        require(Path(name).name==name, 'Native plan inputs must use local file names.')
        write_new(output/name,canonical(data))
    nested=prepare_exploratory_plan(bundle,assessment,resources,policy)
    validate_exploratory_plan(nested,assessment,resources)
    plan=input_only_plan(bundle,assessment=assessment,public_scope_sha256=json_identity(nested)['sha256'],methods=methods,draft=draft,streamed=True)
    plan.update(schema='genosapiens.track2-synthesis-plan/v5',assessment_scope=scope_from_assessment(assessment),
        exploratory_discovery=nested)
    return stored_json(plan), nested['excluded_scenarios']


def load_plan_resources(plan_path, nested, reference_root):
    from genosapiens.track2_discovery_resources import load_discovery_resources
    from genosapiens.track2_published_genotype import _source_file
    parent=Path(plan_path).parent
    descriptor_raw=_source_file(parent,nested['resources']['path']).read_bytes()
    policy_raw=_source_file(parent,nested['policy']['path']).read_bytes()
    require(sha256(descriptor_raw)==nested['resources']['sha256'] and sha256(policy_raw)==nested['policy']['sha256'],
        'A native plan resource or policy binding changed.')
    descriptor=json.loads(descriptor_raw)
    require(descriptor['build_id']==nested['resources']['build_id'], 'The native catalog build differs.')
    return dict(descriptor=descriptor,policy=json.loads(policy_raw),descriptor_path=nested['resources']['path'],
        policy_path=nested['policy']['path'],runtime_images=nested['runtime_images'],
        connection=load_discovery_resources(descriptor,reference_root,matrix_required=False))


def validate_native_plan(plan, resources):
    require(plan['input_origin']=='diagnosis_guided_vcf' and plan['lincs'] is None and plan['frozen_assessment'] is None,
        'Native and legacy discovery inputs cannot be mixed.')
    require(not any(plan[k] for k in ['candidate_inputs','routes','requests','objective_candidates',
        'candidate_route_dependencies','profile_inputs','evidence_imports'])
        and not any(plan['tables'][k] for k in ['endpoint_records','target_records','candidate_records']),
        'A native hypothesis cannot fabricate study endpoints or candidate evidence.')
    validate_exploratory_plan(plan['exploratory_discovery'],plan['mechanism_assessment'],resources)
    from genosapiens.track2_scientific_storage import json_identity
    require(plan['scope_sha256']==json_identity(plan['exploratory_discovery'])['sha256'], 'Native scope identity changed.')
