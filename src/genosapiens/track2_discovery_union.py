"""Canonical medication union, objective overlap, and separate method recurrence."""
from dataclasses import asdict
import json

from genosapiens.track2_candidate_union import CandidateRouteObservation,build_candidate_union
from genosapiens.track2_drug_identity import resolve_drug_identity
from genosapiens.track2_medication_catalog import identity_assertions


def _identity(resolution):
    assertions=identity_assertions(resolution)
    return resolve_drug_identity(assertions[0].identifier,assertions)


def mechanism_observations(results):
    observations=[]
    for result in results:
        for nomination in result['nominations']:
            for record in nomination['records']:
                for identity in nomination['identities']:
                    if identity['compound_id'] not in record['source']['chemblIds']:continue
                    observations.append(CandidateRouteObservation(_identity(identity),result['objective_id'],
                        'qualified_target_drug_linkage',record['mechanism_id'],
                        'Native restoration hypothesis: '+record['action'],
                        ('native_hypothesis',result['gene_id'],record['mechanism_id']),
                        (record['mechanism_id'],),('OpenTargets/ChEMBL',)))
    return observations


def lincs_observations(results,scenario_bindings):
    observations=[]
    for binding in scenario_bindings:
        result=results[binding['job_id']]
        for medication in result['medications']:
            if medication['state']!='pass':continue
            query=result['query_sha256']
            observations.append(CandidateRouteObservation(_identity(medication['identity']),binding['objective_id'],
                'lincs_cell_specific_expression','lincs-query:'+query,
                'Native expression hypothesis: '+result['method'],('native_hypothesis',result['gene_id'],'lincs-query:'+query),
                ('lincs-query:'+query,),('GSE92742',)))
    return observations


def build_exploratory_union(observations,negative_searches=()):
    return build_candidate_union(observations,negative_searches)


def objective_concordance(candidate,objective,annotations):
    observations=[r for r in candidate.rationales if r.objective_id==objective['objective_id']]
    routes={r.route for r in observations};cells=annotations['cells']
    passes={method:sorted({r['cell_id'] for r in cells if r['method']==method and r['state']=='pass'})
        for method in ['paired_oe_kd','kd_only_reversal','oe_only_matching']}
    failed=sorted({r['cell_id'] for r in cells if r['state']=='fail'})
    undefined=sorted({r['cell_id'] for r in cells if r['state']=='unassessed'})
    conflicts=sorted({r['cell_id'] for r in cells if r.get('conflict')})
    paired=passes['paired_oe_kd'];kd=passes['kd_only_reversal'];oe=passes['oe_only_matching']
    return dict(canonical_active_ingredient_id=candidate.active_ingredient_id,gene_id=objective['gene_id'],
        gene_symbol=objective['gene_symbol'],scenario_id=objective['scenario_id'],objective_id=objective['objective_id'],
        mechanism_nominated='qualified_target_drug_linkage' in routes,lincs_nominated='lincs_cell_specific_expression' in routes,
        both_routes=len(routes)==2,direction_agreement=annotations['direction_agreement'],
        lincs_methods=sorted({r['method'] for r in cells if r['state']=='pass'}),paired_passing_cells=paired,
        kd_only_passing_cells=kd,oe_only_passing_cells=oe,failed_cells=failed,undefined_cells=undefined,
        condition_conflicts=conflicts,target_group_warnings=annotations['target_group_warnings'],
        applicability_conflicts=annotations['conflicts'],source_record_ids=[r['mechanism_id'] for r in annotations['records']],
        paired_research_eligible=len(paired)>=2 and not failed and not conflicts,
        fallback_recurrence='mixed_single_direction' if kd and oe else 'repeated_kd_only' if len(kd)>=2 else
            'repeated_oe_only' if len(oe)>=2 else 'single_direction' if kd or oe else 'none',
        cells=cells,mechanisms=annotations['records'])


def reconcile_results(plan,mechanisms,jobs,connection,policy,outcomes=(),*,queries=None):
    if plan['schema']=='genosapiens.track2-exploratory-plan/v2':
        from genosapiens.track2_tissue_union import reconcile_tissue_results
        return reconcile_tissue_results(plan,mechanisms,jobs,queries,connection,policy)
    from genosapiens.track2_drug_mechanisms import annotate_candidate
    from genosapiens.track2_candidate_union import NegativeRouteSearch
    # Add FDA identity records to retained scientific results without modifying job files.
    local={}
    for key,result in jobs.items():
        medications=[]
        for row in result['medications']:
            found=connection.execute("SELECT data FROM compound_identity WHERE provider='lincs' AND ingredient_key=? ORDER BY compound_id LIMIT 1",
                (row['medication_entity_key'],)).fetchone()
            medications.append({**row,'identity':json.loads(found[0])})
        local[key]={**result,'medications':medications}
    observations=mechanism_observations(mechanisms)+lincs_observations(local,plan['scenario_job_bindings'])
    negatives=[NegativeRouteSearch(r['objective_id'],r['route'],'unavailable' if r['state']=='unavailable' else 'no_result',
        tuple(r['query_ids']),('GSE92742' if r['route']=='lincs_cell_specific_expression' else 'OpenTargets/ChEMBL',),
        ', '.join(r['reason_codes_result']) or 'Completed search found no candidates.')
        for r in outcomes if r['state'] in {'completed_no_candidates','unavailable'}]
    union=build_exploratory_union(observations,negatives)
    identities={_identity(row['identity']).active_ingredient_ids[0]:row['medication_entity_key']
        for result in local.values() for row in result['medications']}
    identities.update({_identity(identity).active_ingredient_ids[0]:identity['ingredient_key']
        for result in mechanisms for nomination in result['nominations'] for identity in nomination['identities']})
    from genosapiens.track2_scientific_records import stored_json,ordered_records
    from genosapiens.track2_scientific_storage import Records
    def candidate_rows():
        for candidate in union.candidates:
            key=identities[candidate.active_ingredient_id]
            medication=json.loads(connection.execute('SELECT data FROM medications WHERE ingredient_key=?',(key,)).fetchone()[0])
            for objective in plan['objectives']:
                if not any(r.objective_id==objective['objective_id'] for r in candidate.rationales):continue
                annotation=annotate_candidate(key,objective,connection,policy)
                annotation['cells']=[m for binding in plan['scenario_job_bindings'] if binding['objective_id']==objective['objective_id']
                    for m in jobs[binding['job_id']]['medications'] if m['medication_entity_key']==key]
                yield dict(ingredient_key=key,name=' / '.join(medication['names']),
                    **objective_concordance(candidate,objective,annotation))
    concordance=stored_json(Records(ordered_records(candidate_rows(),
        key=lambda r:(not r['both_routes'],not r['mechanism_nominated'],r['canonical_active_ingredient_id'],r['objective_id']))))
    return dict(candidate_union=asdict(union),concordance=concordance)
