"""Replay retained tissue counts without the original catalog or matrix."""
from pathlib import Path
import sqlite3
import tempfile

from genosapiens.track2_operator_inputs import require
from genosapiens.track2_scientific_reader import read_scientific
from genosapiens.track2_scientific_storage import json_identity, json_chunks
from genosapiens.track2_tissue_aggregation import tissue_result, objective_score, SCOPE, CheckedTissueQuery
from genosapiens.track2_tissue_union import _indexed_evidence


def check_tissue_scores(root):
    root=Path(root)
    plan=read_scientific(root/'plan-source.json').get('exploratory_discovery')
    if not plan or plan['schema']!='genosapiens.track2-exploratory-plan/v2':return
    result=read_scientific(root/'discovery/results.json')
    jobs={}
    for item in result['jobs']:
        if item['completion_sha256']:
            paths=list((root/'lincs').glob('*/'+item['job_id']+'/science/result.json'))
            require(len(paths)==1, 'A tissue score has no retained scientific job.')
            jobs[item['job_id']]=read_scientific(paths[0])
        else:jobs[item['job_id']]={'state':item['state'],'medications':[]}
    queries={q['query_job_id']:q['result'] or {} for q in read_scientific(root/'discovery/tissue-queries.json')}
    planned={r['job_id']:r for r in plan['genetic_selections']}
    coverage={r['objective_id']:r for r in plan['coverage']}
    mechanisms={r['objective_id']:r for r in read_scientific(root/'discovery/mechanism-results.json')}
    with tempfile.TemporaryDirectory(prefix='track2-report-tissue-proof-') as temporary:
        db=sqlite3.connect(Path(temporary)/'scores.sqlite')
        try:
            db.execute('PRAGMA cache_size=-2048')
            db.execute('CREATE TABLE scores (objective TEXT, ingredient TEXT, record BLOB, PRIMARY KEY(objective,ingredient))')
            for score in read_scientific(root/'discovery/tissue-scores.json'):
                try:db.execute('INSERT INTO scores VALUES (?,?,?)',(score['objective_id'],score['medication_entity_key'],b''.join(json_chunks(score))))
                except sqlite3.IntegrityError:raise ValueError('A retained tissue score is duplicated.') from None
            db.commit();checked=set();checked_queries={}
            with _indexed_evidence(jobs,queries) as (summaries,query):
                for objective in plan['objectives']:
                    oid=objective['objective_id'];inventory=coverage[oid]['tissues']
                    ingredients={c['medication_entity_key'] for b in plan['scenario_job_bindings'] if b['objective_id']==oid
                        for c in planned[b['job_id']]['conditions']}
                    ingredients.update(r['ingredient_key'] for n in mechanisms[oid]['nominations'] for r in n['identities'])
                    for ingredient in sorted(ingredients):
                        scope=dict(medication_entity_key=ingredient,**{k:objective[k] for k in SCOPE if k!='medication_entity_key'})
                        tissues=[]
                        for tissue in inventory:
                            query_id=tissue['query_job_id']
                            if query_id not in checked_queries:
                                raw_query=query(query_id)
                                checked_queries[query_id]=CheckedTissueQuery(raw_query) if raw_query is not None else None
                            genetic=checked_queries[query_id];pending=[] if genetic is not None else [query_id]
                            cells=[]
                            for job in tissue['drug_job_ids']:
                                if jobs.get(job,{}).get('state') not in {'completed','unavailable'}:pending.append(job)
                                cells.extend(summaries(job,ingredient))
                            tissues.append(tissue_result(scope,tissue['tissue_id'],genetic,cells,incomplete_jobs=pending))
                        expected=objective_score(scope,[r['tissue_id'] for r in inventory],tissues)
                        raw=db.execute('SELECT record FROM scores WHERE objective=? AND ingredient=?',(oid,ingredient)).fetchone()
                        require(raw is not None, 'A complete evaluated tissue score is missing.')
                        import json
                        saved=json.loads(raw[0])
                        require(json_identity(expected)==json_identity({k:saved[k] for k in expected}),
                                'A tissue score differs from its retained directional evidence.')
                        checked.add((oid,ingredient))
            require(db.execute('SELECT count(*) FROM scores').fetchone()[0]==len(checked),
                    'The evaluated tissue inventory differs from the declared scientific scope.')
            # Presentation copies cannot invent another version of the historical counts.
            for path,rows in (
                ('discovery/results.json',result['evaluated_tissue_scores']),
                ('discovery/concordance.json',(r['tissue_score'] for r in read_scientific(root/'discovery/concordance.json'))),
                ('supplements/objective-scores.json',(r['tissue_score'] for r in read_scientific(root/'supplements/objective-scores.json')))):
                for row in rows:
                    raw=db.execute('SELECT record FROM scores WHERE objective=? AND ingredient=?',
                        (row['objective_id'],row['medication_entity_key'])).fetchone()
                    require(raw is not None, 'A report copy contains an unknown tissue objective.')
                    import json
                    saved=json.loads(raw[0])
                    require(all(row[k]==saved[k] for k in row), 'A report copy changed the retained tissue score.')
        finally:db.close()
