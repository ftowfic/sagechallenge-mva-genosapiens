"""Unpooled drug observations on a frozen genetic tissue query."""
from collections import defaultdict
from dataclasses import dataclass
from statistics import mean, median

from genosapiens.track2_lincs import number
from genosapiens.track2_lincs_scoring import RankedVector, prepare_ranks, correlate_ranks
from genosapiens.track2_operator_inputs import require
from genosapiens.track2_tissue_records import schema_name, validate_tissue_record

WITNESSES = ('oe_positive', 'kd_negative', 'oe_negative', 'kd_positive')


@dataclass(frozen=True)
class PreparedTissueQuery:
    gene_id: str
    tissue_id: str
    method: str
    query_sha256: str
    gene_ids: tuple
    oe: tuple
    kd: tuple
    has_oe: bool
    has_kd: bool
    ready: bool
    reason: str | None
    oe_ranks: RankedVector | None
    kd_ranks: RankedVector | None


def prepare_tissue_query(query, *, counters=None):
    """Validate once and copy only immutable scoring inputs."""
    validate_tissue_record(query, 'query')
    genes = [r for r in query['genes'] if r['selected']]
    ids = tuple(r['gene_id'] for r in genes)
    require(len(ids) == len(set(ids)), 'A prepared query requires unique selected features.')
    oe, kd = tuple(r['oe'] for r in genes), tuple(r['kd'] for r in genes)
    ready = query['state'] == 'ready'
    oe_ranks = prepare_ranks(oe, counters=counters, counter='genetic_rank_calculations') if ready and query['contributors']['oe'] else None
    kd_ranks = prepare_ranks(kd, counters=counters, counter='genetic_rank_calculations') if ready and query['contributors']['kd'] else None
    context = PreparedTissueQuery(query['gene_id'], query['tissue_id'], query['method'],
        query['query_sha256'], ids, oe, kd,
        bool(query['contributors']['oe']), bool(query['contributors']['kd']),
        ready, query['reason'], oe_ranks, kd_ranks)
    if counters is not None:
        counters['query_preparations'] = counters.get('query_preparations', 0) + 1
    return context


def score_tissue_condition(query, condition, values):
    return score_prepared_tissue_condition(prepare_tissue_query(query), condition, values)


def score_prepared_tissue_condition(context, condition, values, *, counters=None):
    require(isinstance(context, PreparedTissueQuery), 'Scoring requires a prepared tissue query.')
    vector = [number(values.get(gene)) for gene in context.gene_ids]
    ready = context.ready
    drug_ranks = prepare_ranks(vector, counters=counters, counter='drug_rank_calculations') if ready else None
    ro = correlate_ranks(context.oe_ranks, drug_ranks)
    rk = correlate_ranks(context.kd_ranks, drug_ranks)
    desired = -rk if context.method=='kd_only_reversal' and rk is not None else ro
    valid = desired is not None and (context.method!='paired_oe_kd' or ro is not None and rk is not None)
    record = dict(schema=schema_name('condition'), gene_id=context.gene_id, tissue_id=context.tissue_id,
        method=context.method, query_sha256=context.query_sha256,
        **{k:condition[k] for k in ('medication_entity_key','sig_id','pert_id','pert_iname','cell_id',
            'pert_dose','pert_dose_unit','pert_time','pert_time_unit')},
        quality={k:condition[k] for k in ('tas','distil_cc_q75','distil_nsample')},
        gene_ids=list(context.gene_ids), values=vector, r_oe=ro, r_kd=rk, r_desired=desired,
        state='valid' if valid else 'undefined' if ready else 'unavailable',
        reason=None if valid else 'missing_or_constant_condition_vector' if ready else context.reason)
    return validate_tissue_record(record, 'condition')


def summarize_tissue_conditions(conditions):
    """Retain witnesses and counts without storing expression vectors in summaries."""
    groups = {}; seen = set(); scopes = set()
    for row in conditions:
        validate_tissue_record(row, 'condition')
        key = row['medication_entity_key']
        require(row['sig_id'] not in seen, 'A drug-cell job has duplicate conditions.')
        seen.add(row['sig_id'])
        scopes.add(tuple(row[k] for k in ('gene_id','tissue_id','query_sha256','method','cell_id')))
        require(len(scopes)==1, 'A drug-cell summary cannot combine query scopes or cells.')
        result = groups.setdefault(key, dict(medication_entity_key=key, name=row['pert_iname'],
            gene_id=row['gene_id'], tissue_id=row['tissue_id'], query_sha256=row['query_sha256'],
            cell_id=row['cell_id'], method=row['method'], condition_ids=[], n_valid=0, n_undefined=0,
            witnesses={k:[] for k in WITNESSES}, undefined_reasons={}, fallback_times=defaultdict(list)))
        result['condition_ids'].append(row['sig_id'])
        if row['state']!='valid':
            result['n_undefined'] += 1
            result['undefined_reasons'][row['sig_id']] = row['reason']
            continue
        result['n_valid'] += 1
        if row['method']=='paired_oe_kd':
            for name, passed in [('oe_positive',row['r_oe']>0), ('kd_negative',row['r_kd']<0),
                    ('oe_negative',row['r_oe']<0), ('kd_positive',row['r_kd']>0)]:
                if passed: result['witnesses'][name].append(row['sig_id'])
        else:
            result['fallback_times'][row['pert_time']].append(row['r_desired'])
    for row in groups.values():
        times = row.pop('fallback_times')
        row['weaker_nomination'] = bool(times) and row['n_undefined']==0 and mean(median(v) for v in times.values())>0
        row['n_conditions'] = len(row['condition_ids'])
        row['condition_ids'].sort()
        for values in row['witnesses'].values(): values.sort()
    return [groups[k] for k in sorted(groups)]
