"""Frozen local refinement and correlation-only discovery; no scientific gate upgrades."""
from collections import Counter, defaultdict
from dataclasses import asdict, dataclass
import json
import math
from pathlib import Path
from statistics import mean, median

from genosapiens.track2_lincs import number, unique
from genosapiens.track2_published_genotype import canonical, sha256

POLICY = {
    'schema': 'genosapiens.track2-lincs-policy/v1',
    'minimum_absolute_modz': 1, 'minimum_features': 3,
    'genetic_qc': {'distil_nsample': 2, 'distil_cc_q75': .2, 'tas': .2},
    'drug_qc': {'distil_nsample': 2, 'distil_cc_q75': .2, 'tas': .2},
    'drug_times_hours': [6, 24], 'features': 'measured_landmarks',
    'correlation': 'spearman_average_tied_ranks',
    'joint': 'min(r_desired,-r_opposite)',
    'aggregation': 'median_all_conditions_per_time_then_equal_time_mean',
    'undefined_condition': 'retain_and_make_medication_unassessed',
    'rank_ties': 'competition_rank_stable_ingredient_display',
    'significance': 'not_assessed', 'top_n_display': 5,
}


def ranks(values):
    ordered = sorted(range(len(values)), key=values.__getitem__)
    result = [0.] * len(values)
    i = 0
    while i < len(ordered):
        j = i + 1
        while j < len(ordered) and values[ordered[j]] == values[ordered[i]]:
            j += 1
        for k in range(i, j):
            result[ordered[k]] = (i + 1 + j) / 2
        i = j
    return result


def spearman(a, b):
    if len(a) != len(b):
        raise ValueError('misaligned correlation vectors')
    if len(a) < 3 or any(v is None or not math.isfinite(v) for v in (*a, *b)):
        return None
    x, y = ranks(a), ranks(b)
    mx, my = mean(x), mean(y)
    xx, yy = sum((v-mx)**2 for v in x), sum((v-my)**2 for v in y)
    if not xx or not yy:
        return None
    return max(-1., min(1., sum((u-mx)*(v-my) for u,v in zip(x,y)) / math.sqrt(xx*yy)))


@dataclass(frozen=True)
class RankedVector:
    values: tuple
    mean: float
    squared_deviations: float


def prepare_ranks(values, *, counters=None, counter='rank_calculations'):
    """Retain the exact rank, mean, and reduction order used by spearman."""
    if len(values) < 3 or any(v is None or not math.isfinite(v) for v in values):
        return None
    x = tuple(ranks(values))
    mx = mean(x)
    xx = sum((v-mx)**2 for v in x)
    if counters is not None:
        counters[counter] = counters.get(counter, 0) + 1
    return RankedVector(x, mx, xx)


def correlate_ranks(x, y):
    if x is None or y is None:
        return None
    if len(x.values) != len(y.values):
        raise ValueError('misaligned correlation vectors')
    xx, yy = x.squared_deviations, y.squared_deviations
    if not xx or not yy:
        return None
    return max(-1., min(1., sum((u-x.mean)*(v-y.mean) for u,v in zip(x.values,y.values)) / math.sqrt(xx*yy)))


def refine_query(cell, cell_config, *, scenario_id, target, source_sha256):
    if not target.ready_for_drug_retrieval or target.required_modulation not in {'increase', 'decrease'}:
        raise ValueError('LINCS requires an already qualified signed target')
    if cell['cell_id'] != cell_config['cell_id']:
        raise ValueError('cross-cell query borrowing prohibited')
    metadata = unique(cell['genetic_metadata'], 'sig_id')
    oe, kd = cell_config['oe_signature'], cell_config['kd_signature']
    reason = None
    for sig, kind in [(oe, 'trt_oe'), (kd, 'trt_sh.cgs')]:
        m = metadata.get(sig)
        if m is None:
            reason = 'missing_local_genetic_direction'; break
        if m['cell_id'] != cell['cell_id'] or m['pert_iname'] != target.target_id or m['pert_type'] != kind:
            raise ValueError('genetic target/cell/perturbation identity mismatch')
        if any(number(m.get(k)) is None or number(m[k]) < v for k,v in POLICY['genetic_qc'].items()):
            reason = 'failed_or_missing_genetic_QC'; break
    rows = []
    desired_key = oe if target.required_modulation == 'increase' else kd
    opposite_key = kd if desired_key == oe else oe
    for source in cell['genetic']:
        o, k = number(source.get(oe)), number(source.get(kd))
        excluded = reason
        if not excluded:
            if o is None or k is None: excluded = 'missing_value'
            elif o*k >= 0: excluded = 'not_opposing'
            elif min(abs(o),abs(k)) < POLICY['minimum_absolute_modz']: excluded = 'below_magnitude'
        desired = o if desired_key == oe else k
        opposite = k if desired_key == oe else o
        rows.append({'gene_id':source['gene_id'], 'gene_symbol':source['gene_symbol'],
                     'feature_class':'measured_landmark', 'oe':o, 'kd':k,
                     'desired':desired, 'opposite':opposite, 'selected':not excluded,
                     'reason':excluded or 'opposing_and_magnitude',
                     'desired_direction':('up' if desired > 0 else 'down' if desired < 0 else 'zero')
                                          if desired is not None else 'missing'})
    selected = [r for r in rows if r['selected']]
    if not reason and (len(selected) < POLICY['minimum_features'] or
            spearman([r['desired'] for r in selected], [r['opposite'] for r in selected]) is None):
        reason = 'inadequate_query'
    query = {'schema':'genosapiens.track2-lincs-query/v1', 'cell_id':cell['cell_id'],
        'scenario_id':scenario_id, 'target_id':target.target_id, 'objective_id':target.objective_id,
        'required_modulation':target.required_modulation, 'target_sha256':sha256(canonical(asdict(target))),
        'source_sha256':source_sha256, 'policy':POLICY, 'oe_signature':oe, 'kd_signature':kd,
        'desired_signature':desired_key, 'opposite_signature':opposite_key,
        'genetic_context':[metadata[s] for s in [oe,kd] if s in metadata],
        'genes':rows, 'gene_count':len(selected), 'state':'unassessed' if reason else 'ready',
        'reason':reason, 'selection_counts':dict(sorted(Counter(r['reason'] for r in rows).items())),
        'full_oe_kd_spearman':spearman([r['oe'] for r in rows],[r['kd'] for r in rows]),
        'selected_oe_kd_spearman':spearman([r['oe'] for r in selected],[r['kd'] for r in selected]),
        'limitations':['Opposition was deliberately selected; exploratory, not independent validation.',
                      'Genetic timing and drug 6/24h timing differ; these are not patient-genotype rescue assays.',
                      'Response features are not drug targets; OE-like expression does not restore absent protein.']}
    query['query_sha256'] = sha256(canonical(query))
    return query


def validate_query(query, cell):
    from genosapiens.track2_scientific_storage import json_identity
    content = {k:v for k,v in query.items() if k != 'query_sha256'}
    expected=POLICY
    if query['schema']=='genosapiens.track2-lincs-query/v2':
        # Historical cell queries keep their original exact policy. Tissue
        # queries use a distinct schema and validator.
        expected=json.loads((Path(__file__).resolve().parents[2]/
            'config/track2/exploratory-policy-v1.json').read_bytes())['lincs']
    if json_identity(content)['sha256'] != query['query_sha256'] or query['policy'] != expected:
        raise ValueError('frozen query or policy changed')
    if query['cell_id'] != cell['cell_id'] or [r['gene_id'] for r in query['genes']] != cell['genes']:
        raise ValueError('query and drug gene/cell alignment mismatch')


def score_cell(query, cell):
    validate_query(query, cell)
    selected = [i for i,r in enumerate(query['genes']) if r['selected']]
    desired = [query['genes'][i]['desired'] for i in selected]
    opposite = [query['genes'][i]['opposite'] for i in selected]
    scores = []
    method=query.get('method','paired_oe_kd')
    paired=method=='paired_oe_kd'
    for condition in cell['conditions']:
        if condition['cell_id'] != query['cell_id'] or len(condition['values']) != len(cell['genes']):
            raise ValueError('condition gene/cell alignment mismatch')
        values = [condition['values'][i] for i in selected]
        rd = spearman(desired,values) if query['state'] == 'ready' else None
        ro = spearman(opposite,values) if query['state'] == 'ready' and paired else None
        valid = rd is not None and (not paired or ro is not None)
        scores.append({k:v for k,v in condition.items() if k != 'values'} | {
            'query_sha256':query['query_sha256'], 'gene_count':len(selected),
            'gene_ids':[query['genes'][i]['gene_id'] for i in selected], 'values':values,
            **({'method':method} if query['schema']=='genosapiens.track2-lincs-query/v2' else {}),
            'r_desired':rd, 'r_opposite':ro, 'joint':min(rd,-ro) if valid and paired else None,
            'state':'pass' if valid and rd > 0 and (not paired or ro < 0) else 'fail' if valid else 'unassessed',
            'reason':None if valid else query['reason'] or 'missing_or_constant_condition_vector'})
    return {'cell_id':cell['cell_id'], 'query_sha256':query['query_sha256'], 'gene_count':len(selected),
            **({'schema':'genosapiens.track2-lincs-result/v2','method':method} if query['schema']=='genosapiens.track2-lincs-query/v2' else {}),
            'conditions':scores, 'medications':summarize_conditions(scores), 'exclusions':cell['exclusions']}


def summarize_conditions(scores):
    unique(scores, 'sig_id')
    groups = defaultdict(list)
    for row in scores:
        groups[row['medication_entity_key']].append(row)
    result = []
    for entity, rows in sorted(groups.items()):
        if len({(r['cell_id'],r['query_sha256'],r.get('method','paired_oe_kd')) for r in rows}) != 1:
            raise ValueError('cannot pool cell/query scores')
        method=rows[0].get('method','paired_oe_kd');paired=method=='paired_oe_kd'
        times = sorted({number(r['pert_time']) for r in rows})
        valid = all(r['state'] != 'unassessed' for r in rows)
        time_summaries = []
        for time in times:
            local = [r for r in rows if number(r['pert_time']) == time]
            time_summaries.append({'time_hours':time, 'n_conditions':len(local),
                **{key:median(r[key] for r in local) if all(r[key] is not None for r in local) else None
                   for key in ['r_desired','r_opposite','joint']}})
        aggregate = {key:mean(t[key] for t in time_summaries) if valid and (paired or key=='r_desired') else None
                     for key in ['r_desired','r_opposite','joint']}
        agreements = sum(r['state'] == 'pass' for r in rows)
        passed = valid and agreements > 0 and aggregate['r_desired'] > 0 and (not paired or aggregate['r_opposite'] < 0 and aggregate['joint'] > 0)
        result.append({'medication_entity_key':entity,'name':rows[0]['pert_iname'],
            **({'method':method} if 'method' in rows[0] else {}),
            'cell_id':rows[0]['cell_id'],'query_sha256':rows[0]['query_sha256'],'gene_count':rows[0]['gene_count'],
            **aggregate,'state':'pass' if passed else 'fail' if valid else 'unassessed',
            'rank':None,'times_hours':times,'n_conditions':len(rows),'n_agree':agreements,
            'agreement_fraction':agreements/len(rows),'conflict':any(r['state']=='fail' for r in rows),
            'n_undefined':sum(r['state']=='unassessed' for r in rows),
            'joint_range':[min(r['joint'] for r in rows),max(r['joint'] for r in rows)] if valid and paired else None,
            'time_summaries':time_summaries,'signature_ids':sorted(r['sig_id'] for r in rows),
            'pert_ids':sorted({r['pert_id'] for r in rows})})
    score=lambda row:row['joint'] if row.get('method','paired_oe_kd')=='paired_oe_kd' else row['r_desired']
    if len({(r['cell_id'],r['query_sha256'],r.get('method','paired_oe_kd')) for r in result})>1:
        raise ValueError('cannot rank different cells, queries or methods together')
    passing = sorted((r for r in result if r['state']=='pass'), key=lambda r:(-score(r),r['medication_entity_key']))
    for i,row in enumerate(passing):
        row['rank'] = passing[i-1]['rank'] if i and score(row)==score(passing[i-1]) else i+1
    return sorted(result,key=lambda r:(r['rank'] is None,r['rank'] or 0,r['medication_entity_key']))


def refine_exploratory_query(selection, genetic, features, *, source_sha256, policy):
    """Freeze a same-cell query using genetic inputs only, before drug extraction."""
    from genosapiens.track2_lincs_selection import select_cell_method
    from genosapiens.track2_operator_inputs import require
    method=selection['method'];selected=selection['selected']
    require(method==select_cell_method(selected,policy) and method is not None,'Genetic method selection changed.')
    ids=[f['pr_gene_id'] for f in features]
    require(len(ids)==len(set(ids)) and all(f['measured_landmark']==1 for f in features),
        'Native queries require unique measured landmark features.')
    for direction,row in selected.items():
        require(row['cell_id']==selection['cell_id'] and row['gene_id']==selection['gene_id']
            and row['pert_type']==policy['lincs']['genetic_types'][direction], 'Genetic gene or cell identity differs.')
        require(row['sig_id'] in genetic and len(genetic[row['sig_id']])==len(features),'A selected matrix column is missing or misaligned.')
    oe=selected.get('oe',{}).get('sig_id');kd=selected.get('kd',{}).get('sig_id')
    rows=[];minimum=policy['lincs']['minimum_absolute_modz']
    for i,feature in enumerate(features):
        o=number(genetic[oe][i]) if oe else None;k=number(genetic[kd][i]) if kd else None
        values=[v for key,v in [('oe',o),('kd',k)] if key in selected]
        reason=('missing_value' if any(v is None for v in values) else
            'not_opposing' if method=='paired_oe_kd' and o*k>=0 else
            'below_magnitude' if any(abs(v)<minimum for v in values) else None)
        desired=-k if method=='kd_only_reversal' and k is not None else o
        rows.append(dict(gene_id=feature['pr_gene_id'],gene_symbol=feature['symbol'],feature_class='measured_landmark',
            oe=o,kd=k,desired=desired,opposite=k if method=='paired_oe_kd' else None,selected=reason is None,
            reason=reason or ('opposing_and_magnitude' if method=='paired_oe_kd' else 'magnitude'),
            desired_direction='missing' if desired is None else 'up' if desired>0 else 'down' if desired<0 else 'zero'))
    admitted=[r for r in rows if r['selected']]
    desired=[r['desired'] for r in admitted]
    reason=None
    if len(admitted)<policy['lincs']['minimum_features'] or spearman(desired,desired) is None or (
            method=='paired_oe_kd' and spearman(desired,[r['opposite'] for r in admitted]) is None):
        reason='inadequate_paired_query' if method=='paired_oe_kd' else 'inadequate_single_direction_query'
    query=dict(schema='genosapiens.track2-lincs-query/v2',gene_id=selection['gene_id'],cell_id=selection['cell_id'],
        method=method,required_modulation='increase',source_sha256=source_sha256,policy=policy['lincs'],
        oe_signature=oe,kd_signature=kd,desired_signature=kd if method=='kd_only_reversal' else oe,
        opposite_signature=kd if method=='paired_oe_kd' else None,
        desired_transform='negate' if method=='kd_only_reversal' else 'identity',
        absent_direction='oe' if method=='kd_only_reversal' else 'kd' if method=='oe_only_matching' else None,
        genetic_context=[selected[k] for k in sorted(selected)],genes=rows,gene_count=len(admitted),
        state='unassessed' if reason else 'ready',reason=reason,
        selection_counts=dict(sorted(Counter(r['reason'] for r in rows).items())),
        full_oe_kd_spearman=spearman([r['oe'] for r in rows],[r['kd'] for r in rows]) if method=='paired_oe_kd' else None,
        selected_oe_kd_spearman=spearman(desired,[r['opposite'] for r in admitted]) if method=='paired_oe_kd' else None,
        limitations=['Selected opposition is a query rule, not independent validation.' if method=='paired_oe_kd' else
            'Single-direction evidence lacks the opposing genetic perturbation.',
            'Expression matching does not establish patient rescue or protein production.'])
    from genosapiens.track2_scientific_storage import json_identity
    query['query_sha256']=json_identity(query)['sha256']
    return query


def cross_cell_summary(results):
    cells = sorted(results)
    entities = sorted({r['medication_entity_key'] for cell in results.values() for r in cell['medications']})
    indexes = {c:unique(results[c]['medications'],'medication_entity_key') for c in cells}
    rows = []
    for entity in entities:
        local = {c:indexes[c].get(entity, {'medication_entity_key':entity,'cell_id':c,
            'query_sha256':results[c]['query_sha256'],'gene_count':results[c]['gene_count'],
            'state':'unassessed','reason':'no_admitted_local_drug_condition','rank':None,
            'r_desired':None,'r_opposite':None,'joint':None,'n_conditions':0}) for c in cells}
        rows.append({'medication_entity_key':entity,'cells':local,
            'nominated_cells':[c for c in cells if local[c]['state']=='pass'],
            'tested_cells':[c for c in cells if local[c]['n_conditions'] > 0],
            'interpretation':'Descriptive local recurrence; no pooled score or independent replication claim.'})
    return rows
