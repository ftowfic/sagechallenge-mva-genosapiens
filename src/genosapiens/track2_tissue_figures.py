"""Bounded tissue panels and complete condition pages from retained evidence."""
from itertools import groupby, islice

from genosapiens.track2_scientific_records import ordered_records
from genosapiens.track2_tissue_records import identity

TISSUE_CAPTION=(
    'Each cell summarizes expression evidence for one medication and LINCS primary_site tissue group.\n\n'
    '- Blue +1 indicates support.\n- Orange -1 indicates conflict.\n'
    '- Patterned zero indicates both.\n- Plain zero indicates neither.\n'
    '- Gray NA indicates unavailable evidence.\n- Pending marks unfinished scientific work.\n\n'
    'P, C, and N count supporting, conflicting, and evaluable tissues. S_L equals (P - C) / N. '
    'Zero evaluable tissues gives NA. Scores use exact fractions before three-decimal display.\n\n'
    'Expression agreement does not establish treatment benefit or a clinical tissue indication.')
CONDITION_CAPTION=(
    'Each row retains one medication condition with its cell, dose, unit, time, and signature identifier.\n\n'
    'Columns display rOE and -rKD. The second column negates raw KD correlation so positive values indicate KD reversal. '
    'Blue positive values indicate agreement with the restoration direction. Orange negative values indicate the opposite direction. '
    'Zero supplies neither sign. Gray NA indicates an undefined or unavailable comparison.\n\n'
    'Each valid paired condition requires both correlations. Different conditions can supply the two directional witnesses. '
    'The scale extends from -1 to +1. Exports retain raw rKD.\n\n'
    'The source query retains genetic contributors and feature exclusions. Expression agreement does not establish treatment benefit.')


def tissue_panels(scores, *, page_size=20, main=False):
    scope=lambda r:(r['gene_id'],r['scenario_id'],r['objective_id'])
    ordered=ordered_records(scores,key=lambda r:(*scope(r),r.get('name','').casefold(),r['medication_entity_key']))
    for key,group in groupby(ordered,key=scope):
        page=0
        while rows:=list(islice(group,page_size)):
            page+=1
            tissues=sorted({t for r in rows for t in r['declared_tissues']})
            panel=dict(gene_id=key[0],scenario_id=key[1],objective_id=key[2],tissues=tissues,rows=rows)
            yield dict(schema='genosapiens.track2-tissue-figure/v1',kind='tissues',
                id='tissues-'+identity(list(key))[:16]+'-p'+str(page),panels=[panel],rows=rows,
                gene_id=key[0],scenario_id=key[1],objective_id=key[2],tissue_ids=tissues,
                page=page,main=main,selection='main_representatives' if main else 'all_evaluated_medications',
                policy='normalized_tissue_score/v1',caption=TISSUE_CAPTION,
                alt='Medication and tissue evidence. Printed labels and patterns distinguish all evidence states.',
                source_links=['discovery/tissue-scores.json','discovery/tissue-results.csv','discovery/tissue-conditions.csv',
                    'discovery/tissue-queries.json'],missing='Gray NA. Pending means incomplete work.',clipping=False)


def tissue_condition_pages(query,local,path,*,page_size=15):
    from genosapiens.track2_figure_stream import condition_rows
    rows=ordered_records(condition_rows(path),key=lambda r:(r['medication_entity_key'],r['cell_id'],
        r['pert_time'],str(r['pert_dose']),r['pert_dose_unit'] or '',r['sig_id']))
    for medication,group in groupby(rows,key=lambda r:r['medication_entity_key']):
        page=0
        while selected:=list(islice(group,page_size)):
            page+=1
            values=[dict(r,kd_reversal=-r['r_kd'] if r['r_kd'] is not None else None,
                directional_witnesses=[name for name,passed in [('oe_positive',r['r_oe'] is not None and r['r_oe']>0),
                    ('kd_negative',r['r_kd'] is not None and r['r_kd']<0),('oe_negative',r['r_oe'] is not None and r['r_oe']<0),
                    ('kd_positive',r['r_kd'] is not None and r['r_kd']>0)] if passed and r['state']=='valid']) for r in selected]
            yield dict(schema='genosapiens.track2-tissue-figure/v1',kind='tissue_conditions',
                id='conditions-'+identity(medication)[:16]+'-p'+str(page),rows=values,
                gene_id=query['gene_id'],tissue_id=query['tissue_id'],cell_id=local['cell_id'],
                method=query['method'],query_sha256=query['query_sha256'],page=page,
                medication_entity_key=medication,selection='all_admitted_conditions',policy='tissue_condition_correlations/v1',
                caption=CONDITION_CAPTION,alt='Individual condition correlations with pooled OE and reversed KD.',
                source_links=['../science/query.json','../science/conditions.jsonl'],
                missing='Gray NA for undefined or unavailable correlation.',clipping=False)
