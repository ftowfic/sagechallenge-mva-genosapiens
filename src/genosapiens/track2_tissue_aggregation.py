"""Independent directional flags, exact tissue counts, and explicit coverage."""
from dataclasses import dataclass, field, InitVar
from types import MappingProxyType

from genosapiens.track2_operator_inputs import require
from genosapiens.track2_tissue_records import identity, schema_name, validate_tissue_record
from genosapiens.track2_tissue_scoring import WITNESSES

SCOPE = ('medication_entity_key','gene_id','objective_id','scenario_id')


@dataclass(frozen=True)
class CheckedTissueQuery:
    """Immutable metadata copied only after the complete query passes validation."""
    source: InitVar[object]
    _scope: object = field(init=False,repr=False)

    def __post_init__(self, source):
        validate_tissue_record(source,'query')
        object.__setattr__(self,'_scope',MappingProxyType({key:source[key] for key in
            ('gene_id','tissue_id','method','query_sha256','state','reason')}))


def tissue_result(scope, tissue_id, query, summaries, *, incomplete_jobs=()):
    require(set(scope)==set(SCOPE), 'A tissue result requires its complete objective scope.')
    if isinstance(query,CheckedTissueQuery):query=query._scope
    elif query is not None:validate_tissue_record(query,'query')
    if query is not None:
        require(query['gene_id']==scope['gene_id'] and query['tissue_id']==tissue_id,
            'The genetic query belongs to another tissue or gene.')
    witnesses={k:set() for k in WITNESSES}; seen={}; signatures=set(); valid=0; undefined=0; weaker=False
    for row in summaries:
        require(query is not None and all(row[k]==v for k,v in dict(gene_id=scope['gene_id'],
            tissue_id=tissue_id, medication_entity_key=scope['medication_entity_key'],
            query_sha256=query['query_sha256'],method=query['method']).items()),
            'Directional evidence crosses a medication, gene, tissue, or query boundary.')
        key=row['cell_id']
        if key in seen:
            require(seen[key]==identity(row), 'Repeated cell evidence has inconsistent values.')
            continue
        seen[key]=identity(row)
        require(not signatures.intersection(row['condition_ids']), 'Different cells repeat a condition identity.')
        signatures.update(row['condition_ids']); valid+=row['n_valid']; undefined+=row['n_undefined']
        for name in WITNESSES: witnesses[name].update(row['witnesses'][name])
        weaker = weaker or row['weaker_nomination']
    pending=sorted(set(incomplete_jobs)); method=query['method'] if query else None
    evaluable=bool(not pending and query and query['state']=='ready' and method=='paired_oe_kd' and valid)
    support=evaluable and bool(witnesses['oe_positive'] and witnesses['kd_negative'])
    conflict=evaluable and bool(witnesses['oe_negative'] and witnesses['kd_positive'])
    reasons=[]
    if pending: reasons.append('incomplete_scientific_jobs')
    if query is None: reasons.append('query_incomplete')
    elif query['state']!='ready': reasons.append(query['reason'])
    elif method!='paired_oe_kd': reasons.append('missing_genetic_direction')
    if not signatures: reasons.append('absent_drug_coverage')
    elif not valid: reasons.append('wholly_undefined_comparisons')
    if undefined: reasons.append('undefined_conditions_retained')
    state=('pending' if pending else 'unavailable' if not evaluable else 'mixed' if support and conflict else
        'support' if support else 'conflict' if conflict else 'no_two_direction_evidence')
    result=dict(schema=schema_name('result'),**scope,tissue_id=tissue_id,
        query_sha256=query['query_sha256'] if query else None,method=method,evaluable=evaluable,
        support=support,conflict=conflict,contribution=int(support)-int(conflict) if evaluable else None,state=state,
        witnesses={k:sorted(v) for k,v in witnesses.items()},condition_ids=sorted(signatures),
        n_conditions=len(signatures),n_valid=valid,n_undefined=undefined,coverage_reasons=sorted(set(reasons)),
        incomplete_jobs=pending,weaker_nomination=bool(weaker and not pending),legacy_research_eligibility='not_applicable')
    return validate_tissue_record(result,'result')


def objective_score(scope, declared_tissues, results):
    require(set(scope)==set(SCOPE), 'An objective score requires its complete scope.')
    unique={}
    for row in results:
        validate_tissue_record(row,'result')
        require(all(row[k]==scope[k] for k in SCOPE), 'Tissue scores cannot combine objective scopes.')
        tissue=row['tissue_id']
        require(tissue not in unique or unique[tissue]==row, 'Repeated tissue evidence has inconsistent values.')
        unique[tissue]=row
    declared=sorted(set(declared_tissues))
    require(set(unique)==set(declared), 'A tissue score must retain the complete declared inventory.')
    supporting=sorted(t for t,r in unique.items() if r['support'])
    conflicting=sorted(t for t,r in unique.items() if r['conflict'])
    evaluated=sorted(t for t,r in unique.items() if r['evaluable'])
    p,c,n=len(supporting),len(conflicting),len(evaluated)
    record=dict(schema=schema_name('score'),**scope,declared_tissues=declared,supporting_tissues=supporting,
        conflicting_tissues=conflicting,evaluable_tissues=evaluated,P=p,C=c,N=n,numerator=p-c,denominator=n,
        nominated=bool(p or any(r['weaker_nomination'] for r in unique.values())),
        provisional=any(r['incomplete_jobs'] for r in unique.values()),coverage=[unique[t] for t in declared],
        method_labels=sorted({r['method'] for r in unique.values() if r['method']}),
        legacy_research_eligibility='not_applicable')
    return validate_tissue_record(record,'score')
