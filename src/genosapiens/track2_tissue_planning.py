"""Shared prospective tissue scope for both native input origins."""
from genosapiens.track2_operator_inputs import require
from genosapiens.track2_tissue_records import identity

SCHEMA = 'genosapiens.track2-exploratory-plan/v2'


def prepare_tissue_plan(assessment, resources, policy):
    from genosapiens.track2_exploratory_planning import discovery_eligibility, restoration_objective, ROUTES
    from genosapiens.track2_mechanism import assess_mechanism_readiness
    from genosapiens.track2_tissue_admission import TissueAdmission
    from genosapiens.track2_discovery_resources import resource_bindings
    from genosapiens.track2_discovery import scientific_identity
    descriptor = resources['descriptor']; connection = resources['connection']
    require(descriptor['schema']=='genosapiens.track2-discovery-resources/v2' and descriptor['query_policy']==policy,
        'Native tissue discovery requires its prepared companion and compatible policy binding.')
    require(assessment['automatic_inputs']['discovery']==dict(resources_sha256=identity(descriptor),
        policy_sha256=identity(policy)), 'The mechanism assessment binds different tissue discovery inputs.')
    admission = TissueAdmission(descriptor, policy, connection, counters=resources.setdefault('execution_counters', {}))
    scope, routes, coverage, excluded, bindings = [], [], [], [], []
    queries, jobs = {}, {}
    for scenario in sorted(assessment['scenarios'], key=lambda r:r['scenario_id']):
        decision = assess_mechanism_readiness(scenario['mechanism_inputs'])
        eligibility = discovery_eligibility(scenario, decision, policy)
        require(scenario['discovery_eligibility']==eligibility, 'Native eligibility differs from mechanism acceptance.')
        if eligibility['outcome']!='eligible':
            excluded.append(dict(scenario_id=scenario['scenario_id'], gene_id=scenario['gene_id'], **eligibility))
            continue
        gene = connection.execute('SELECT approved_symbol FROM genes WHERE gene_id=?', (scenario['gene_id'],)).fetchone()
        objective = restoration_objective(scenario, dict(gene_symbol=gene[0] if gene else scenario['gene_id']), decision, policy)
        scope.append(objective)
        selected = admission.genetic(objective['gene_id'])
        tissue_inventory = []
        for selection in selected['tissues']:
            query_id = 'tissue-query-'+identity(dict(selection=selection, resources=resource_bindings(descriptor),
                policy_sha256=identity(policy), scientific_code_sha256=scientific_identity()))
            queries[query_id] = dict(selection, schema='genosapiens.track2-tissue-query-selection/v1', query_job_id=query_id)
            drug_jobs = []
            for member in selection['membership']:
                conditions = admission.drugs(member['cell_id'])
                if not conditions['conditions']: continue
                job = dict(schema='genosapiens.track2-tissue-drug-job/v1', gene_id=objective['gene_id'],
                    tissue_id=selection['tissue_id'], query_job_id=query_id, cell_id=member['cell_id'],
                    method=selection['method'], conditions=conditions['conditions'], drug_exclusions=conditions['exclusions'],
                    drug_signature_ids=[r['sig_id'] for r in conditions['conditions']], runtime=resources['runtime_images']['discovery'])
                job_id = 'tissue-drug-'+identity(job)
                jobs[job_id] = dict(job, job_id=job_id)
                bindings.append(dict(scenario_id=objective['scenario_id'], objective_id=objective['objective_id'], job_id=job_id))
                drug_jobs.append(job_id)
            tissue_inventory.append(dict(tissue_id=selection['tissue_id'], query_job_id=query_id,
                method=selection['method'], drug_job_ids=drug_jobs))
        coverage.append(dict(objective_id=objective['objective_id'], gene_id=objective['gene_id'],
            tissues=tissue_inventory, unresolved_cells=selected['unresolved_cells']))
        executable = any(r['method'] and r['drug_job_ids'] for r in tissue_inventory)
        for route in ROUTES:
            available = route==ROUTES[0] or executable
            routes.append(dict(objective_id=objective['objective_id'], route=route,
                planned_state='executable' if available else 'unavailable',
                reason_codes=[] if available else ['no_executable_tissue_comparison'],
                request_id='discovery-request:'+identity(dict(objective=objective['objective_id'], route=route,
                    resources=identity(descriptor), policy=identity(policy)))))
    return dict(schema=SCHEMA, mode='native_exploratory',
        policy=dict(path=resources['policy_path'], sha256=identity(policy)),
        resources=dict(path=resources['descriptor_path'], sha256=identity(descriptor), build_id=descriptor['build_id']),
        objectives=scope, routes=routes, tissue_queries=[queries[k] for k in sorted(queries)],
        genetic_selections=[jobs[k] for k in sorted(jobs)], scenario_job_bindings=bindings,
        excluded_scenarios=excluded, coverage=coverage, runtime_images=resources['runtime_images'])
