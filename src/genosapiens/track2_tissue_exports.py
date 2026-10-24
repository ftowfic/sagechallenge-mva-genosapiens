"""Complete condition, contributor, and objective exports, including nonnominees."""
from pathlib import Path

from genosapiens.track2_discovery_report import _csv
from genosapiens.track2_scientific_storage import publish_json, Records


def publish_tissue_exports(root,plan,result,discovery_path):
    from genosapiens.track2_discovery import _read
    from genosapiens.track2_figure_stream import condition_rows
    root=Path(root);root.mkdir(parents=True,exist_ok=True,mode=0o700)
    nested=plan.get('exploratory_discovery',plan)
    scores=result['evaluated_tissue_scores']
    publish_json(root/'tissue-scores.json',scores)
    _csv(root,'tissue-scores.csv',scores,['medication_entity_key','canonical_active_ingredient_id','name',
        'gene_id','scenario_id','objective_id','P','C','N','numerator','denominator','nominated','provisional',
        'declared_tissues','supporting_tissues','conflicting_tissues','evaluable_tissues','method_labels','coverage'])
    _csv(root,'tissue-results.csv',(r for score in scores for r in score['coverage']),
        ['medication_entity_key','gene_id','scenario_id','objective_id','tissue_id','query_sha256','method',
         'evaluable','support','conflict','contribution','state','witnesses','condition_ids','n_conditions','n_valid',
         'n_undefined','coverage_reasons','incomplete_jobs','weaker_nomination','legacy_research_eligibility'])
    def queries():
        selections={r['query_job_id']:r for r in nested['tissue_queries']}
        for row in result['queries']:
            yield dict(row,selection=selections[row['query_job_id']],
                result=_read(Path(discovery_path)/row['path']) if row['completion_sha256'] else None)
    publish_json(root/'tissue-queries.json',Records(queries()))
    def conditions():
        for job in result['jobs']:
            if not job['completion_sha256']: continue
            bindings=[b for b in nested['scenario_job_bindings'] if b['job_id']==job['job_id']]
            source=Path(discovery_path)/'jobs'/job['job_id']/'conditions.jsonl'
            for row in condition_rows(source):
                for binding in bindings:
                    yield dict(row,objective_id=binding['objective_id'],scenario_id=binding['scenario_id'],
                        job_id=job['job_id'],query_reference='discovery/tissue-queries.json',
                        condition_reference='lincs/'+row['gene_id']+'/'+job['job_id']+'/science/conditions.jsonl')
    _csv(root,'tissue-conditions.csv',conditions(),['medication_entity_key','gene_id','scenario_id','objective_id',
        'tissue_id','query_sha256','job_id','sig_id','pert_id','pert_iname','cell_id','pert_dose','pert_dose_unit',
        'pert_time','pert_time_unit','quality','gene_ids','values','r_oe','r_kd','r_desired','method','state','reason',
        'query_reference','condition_reference'])
