"""Complete tissue evidence stays separate from positive nomination authority."""
from dataclasses import asdict
from contextlib import contextmanager
import json
from pathlib import Path
import sqlite3
import tempfile

from genosapiens.track2_candidate_union import CandidateRouteObservation, build_candidate_union
from genosapiens.track2_discovery_union import _identity, mechanism_observations
from genosapiens.track2_drug_mechanisms import annotate_candidate
from genosapiens.track2_scientific_records import stored_json
from genosapiens.track2_scientific_storage import Records, json_chunks
from genosapiens.track2_tissue_aggregation import tissue_result, objective_score


@contextmanager
def _indexed_evidence(jobs, queries):
    """Index each input once, then decode one query and medication at a time.

    Query validation remains in tissue_result. Loading its unchanged JSON as
    one record avoids repeated scalar reads from the scientific disk view.
    Duplicate summaries retain their original order and validation behavior.
    """
    with tempfile.TemporaryDirectory(prefix='track2-tissue-evidence-') as temporary:
        path=Path(temporary)/'evidence.sqlite'
        connection=sqlite3.connect(path)
        path.chmod(0o600)
        try:
            connection.execute('PRAGMA cache_size=-2048')
            connection.execute('PRAGMA temp_store=FILE')
            connection.execute('CREATE TABLE summaries (job TEXT, medication TEXT, ordinal INTEGER, data TEXT)')
            connection.execute('CREATE INDEX summary_lookup ON summaries(job,medication,ordinal)')
            connection.execute('CREATE TABLE queries (job TEXT PRIMARY KEY, data TEXT)')
            for job_id, result in jobs.items():
                for ordinal, row in enumerate(result.get('medications', [])):
                    connection.execute('INSERT INTO summaries VALUES (?,?,?,?)',
                        (job_id,row['medication_entity_key'],ordinal,b''.join(json_chunks(row)).decode()))
            for query_id, result in queries.items():
                query=result.get('query')
                if query is not None:
                    connection.execute('INSERT INTO queries VALUES (?,?)',
                        (query_id,b''.join(json_chunks(query)).decode()))
            connection.commit()
            def summaries(job_id, medication):
                return (json.loads(row[0]) for row in connection.execute(
                    'SELECT data FROM summaries WHERE job=? AND medication=? ORDER BY ordinal',
                    (job_id,medication)))
            def query(query_id):
                row=connection.execute('SELECT data FROM queries WHERE job=?',(query_id,)).fetchone()
                return json.loads(row[0]) if row else None
            yield summaries, query
        finally:
            connection.close()


def reconcile_tissue_results(plan,mechanisms,jobs,queries,connection,policy):
    queries=queries or {}
    planned={r['job_id']:r for r in plan['genetic_selections']}
    coverage={r['objective_id']:r for r in plan['coverage']}
    mechanism_by_objective={r['objective_id']:r for r in mechanisms}
    def all_scores():
        for objective in plan['objectives']:
            inventory=coverage[objective['objective_id']]['tissues']
            keys={c['medication_entity_key'] for b in plan['scenario_job_bindings']
                if b['objective_id']==objective['objective_id'] for c in planned[b['job_id']]['conditions']}
            keys.update(r['ingredient_key'] for n in mechanism_by_objective[objective['objective_id']]['nominations']
                for r in n['identities'])
            for key in sorted(keys):
                scope=dict(medication_entity_key=key,**{k:objective[k] for k in ('gene_id','objective_id','scenario_id')})
                results=[]
                for tissue in inventory:
                    query=lookup_query(tissue['query_job_id'])
                    incomplete=[tissue['query_job_id']] if query is None else []
                    summaries=[]
                    for job_id in tissue['drug_job_ids']:
                        result=jobs.get(job_id,{})
                        if result.get('state') not in {'completed','unavailable'}:
                            incomplete.append(job_id)
                        summaries.extend(lookup_summaries(job_id,key))
                    results.append(tissue_result(scope,tissue['tissue_id'],query,summaries,incomplete_jobs=incomplete))
                score=objective_score(scope,[r['tissue_id'] for r in inventory],results)
                medication=json.loads(connection.execute('SELECT data FROM medications WHERE ingredient_key=?',(key,)).fetchone()[0])
                resolved=json.loads(connection.execute('SELECT data FROM compound_identity WHERE ingredient_key=? ORDER BY provider,compound_id LIMIT 1',
                    (key,)).fetchone()[0])
                yield dict(score,name=' / '.join(medication['names']),canonical_active_ingredient_id=_identity(resolved).active_ingredient_ids[0])
    with _indexed_evidence(jobs,queries) as (lookup_summaries,lookup_query):
        evaluated=stored_json(Records(all_scores()))
    observations=mechanism_observations(mechanisms)
    identities={_identity(r).active_ingredient_ids[0]:r['ingredient_key'] for m in mechanisms
        for n in m['nominations'] for r in n['identities']}
    for score in evaluated:
        if not score['nominated']: continue
        found=connection.execute("SELECT data FROM compound_identity WHERE provider='lincs' AND ingredient_key=? ORDER BY compound_id LIMIT 1",
            (score['medication_entity_key'],)).fetchone()
        resolved=_identity(json.loads(found[0])); identities[resolved.active_ingredient_ids[0]]=score['medication_entity_key']
        for tissue in score['coverage']:
            if not tissue['support'] and not tissue['weaker_nomination']: continue
            ref='lincs-query:'+tissue['query_sha256']
            observations.append(CandidateRouteObservation(resolved,score['objective_id'],'lincs_cell_specific_expression',ref,
                'Tissue expression hypothesis: '+tissue['method'],('native_tissue_hypothesis',score['gene_id'],ref),
                (ref,),('GSE92742',)))
    union=build_candidate_union(observations)
    def concordance():
        for candidate in union.candidates:
            key=identities[candidate.active_ingredient_id]
            medication=json.loads(connection.execute('SELECT data FROM medications WHERE ingredient_key=?',(key,)).fetchone()[0])
            for objective in plan['objectives']:
                local=[r for r in candidate.rationales if r.objective_id==objective['objective_id']]
                if not local: continue
                score=next(r for r in evaluated if r['medication_entity_key']==key and r['objective_id']==objective['objective_id'])
                annotation=annotate_candidate(key,objective,connection,policy)
                routes={r.route for r in local}
                yield dict(ingredient_key=key,name=' / '.join(medication['names']),
                    canonical_active_ingredient_id=candidate.active_ingredient_id,
                    **{k:objective[k] for k in ('gene_id','gene_symbol','scenario_id','objective_id')},
                    mechanism_nominated='qualified_target_drug_linkage' in routes,
                    lincs_nominated=score['nominated'],both_routes=len(routes)==2,
                    direction_agreement=annotation['direction_agreement'],lincs_methods=score['method_labels'],
                    paired_passing_cells=[],kd_only_passing_cells=[],oe_only_passing_cells=[],failed_cells=[],undefined_cells=[],
                    condition_conflicts=score['conflicting_tissues'],target_group_warnings=annotation['target_group_warnings'],
                    applicability_conflicts=annotation['conflicts'],source_record_ids=[r['mechanism_id'] for r in annotation['records']],
                    paired_research_eligible=False,legacy_research_eligibility='not_applicable',
                    fallback_recurrence='tissue_single_direction' if any(r['weaker_nomination'] for r in score['coverage']) else 'none',
                    cells=[],mechanisms=annotation['records'],tissue_score=score)
    return dict(candidate_union=asdict(union),concordance=stored_json(Records(concordance())),evaluated_tissue_scores=evaluated)
