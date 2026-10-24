"""Explicit selected-gene calibration. Report readers perform lookup only."""
from collections import defaultdict
import hashlib
import json
import math
from statistics import mean, median

from genosapiens.track2_derived_strength import METHODS, policy
from genosapiens.track2_lincs_scoring import RankedVector, prepare_ranks, correlate_ranks


def canonical(value):
    return json.dumps(value,sort_keys=True,separators=(',', ':'),ensure_ascii=True,allow_nan=False).encode()


def identity(value):
    return hashlib.sha256(canonical(value)).hexdigest()


def query_content(query):
    """Exclude run provenance and objective labels from the randomization stream."""
    return dict(gene_id=query['gene_id'],tissue_id=query['tissue_id'],method=query['method'],
        source_sha256=query['source_sha256'],
        genes=[{k:g[k] for k in ('gene_id','oe','kd')} for g in query['genes'] if g['selected']])


class Permutations:
    """SHA-256 counter words with unbiased rejection and Fisher-Yates swaps."""
    def __init__(self, query_key, *, seed=0):
        self.key=canonical(dict(query_key=query_key,seed=seed,generator=policy()['calibration']['generator']))
        self.counter=0

    def below(self, limit):
        ceiling=(1<<256)-((1<<256)%limit)
        while True:
            value=int.from_bytes(hashlib.sha256(self.key+self.counter.to_bytes(16,'big')).digest(),'big')
            self.counter+=1
            if value<ceiling:return value%limit

    def next(self, size):
        indices=list(range(size))
        for index in range(size-1,0,-1):
            other=self.below(index+1)
            indices[index],indices[other]=indices[other],indices[index]
        return indices


def permuted_ranks(ranks, indices):
    if ranks is None:return None
    return RankedVector(tuple(ranks.values[i] for i in indices),ranks.mean,ranks.squared_deviations)


def prepare_distributions(queries, conditions, *, direction, eligible_ingredients, counters=None):
    """Construct 10,000 complete aggregate draws for every evaluable ingredient.

    Callers must admit the FDA scope and retain its binding before publication.
    This function never runs as a report lookup fallback.
    """
    if direction!='increase':raise ValueError('Calibration supports increase-function objectives only.')
    counters={} if counters is None else counters
    prepared={}; aliases={}; scopes={}
    for query in queries:
        content=query_content(query);key=identity(content)
        if query['method'] not in METHODS or query['state']!='ready':
            raise ValueError('Calibration requires supported ready tissue queries.')
        genes=content['genes'];ids=tuple(g['gene_id'] for g in genes)
        if len(ids)!=len(set(ids)):raise ValueError('A query contains repeated gene identifiers.')
        oe=prepare_ranks([g['oe'] for g in genes]) if query['method']!='kd_only_reversal' else None
        kd=prepare_ranks([g['kd'] for g in genes]) if query['method']!='oe_only_matching' else None
        required=(oe,kd) if query['method']=='paired_oe_kd' else (oe,) if oe is not None else (kd,)
        if not required or any(r is None or not r.squared_deviations for r in required):
            raise ValueError('Calibration requires finite nonconstant query ranks.')
        prepared[key]=(ids,oe,kd,Permutations(key,seed=policy()['calibration']['seed']))
        scopes[key]=(query['gene_id'],query['tissue_id'],query['method'])
        previous=aliases.setdefault(query['query_sha256'],key)
        if previous!=key:raise ValueError('A query identity has conflicting scientific content.')
    eligible=set(eligible_ingredients)
    seen={};rows=[];groups=defaultdict(list)
    for row in conditions:
        if row['medication_entity_key'] not in eligible:continue
        query_key=aliases[row['query_sha256']]
        if (row['gene_id'],row['tissue_id'],row['method'])!=scopes[query_key]:
            raise ValueError('A retained condition differs from its query scope.')
        key=(query_key,row['sig_id'])
        observation={k:v for k,v in row.items() if k not in {'scenario_id','objective_id','query_sha256'}}
        binding=identity(observation)
        if key in seen:
            if seen[key]!=binding:raise ValueError('Repeated conditions disagree.')
            continue
        seen[key]=binding
        if row['state']!='valid':continue
        ids,oe,kd,_=prepared[query_key]
        if tuple(row['gene_ids'])!=ids:raise ValueError('Retained condition gene order differs from its query.')
        ranks=prepare_ranks(row['values'])
        if ranks is None or not ranks.squared_deviations:raise ValueError('A valid condition lacks finite nonconstant ranks.')
        ro,rk=correlate_ranks(oe,ranks),correlate_ranks(kd,ranks)
        if (ro,rk)!=(row['r_oe'],row['r_kd']):raise ValueError('Retained correlations differ from their ranks.')
        if row['pert_time_unit']!='h':raise ValueError('Calibration requires retained time in hours.')
        index=len(rows);rows.append((query_key,row['method'],ranks))
        groups[(row['medication_entity_key'],row['method'],row['tissue_id'],row['cell_id'],row['pert_time'])].append(index)
    # Freeze group order once. Every draw uses the observed aggregation hierarchy.
    dose_groups=sorted(groups.items())
    cell_groups=defaultdict(list)
    for index,(key,_) in enumerate(dose_groups):cell_groups[key[:4]].append(index)
    cell_groups=sorted(cell_groups.items());tissue_groups=defaultdict(list)
    for index,(key,_) in enumerate(cell_groups):tissue_groups[key[:3]].append(index)
    tissue_groups=sorted(tissue_groups.items());ingredients=defaultdict(list)
    for index,(key,_) in enumerate(tissue_groups):ingredients[key[:2]].append(index)
    result={key:[] for key in sorted(ingredients)}
    for _ in range(policy()['calibration']['draws']):
        permuted={}
        for key,(ids,oe,kd,stream) in sorted(prepared.items()):
            indices=stream.next(len(ids))
            permuted[key]=(permuted_ranks(oe,indices),permuted_ranks(kd,indices))
            counters['query_permutations']=counters.get('query_permutations',0)+1
        values=[]
        for key,method,drug in rows:
            oe,kd=permuted[key]
            ro,rk=correlate_ranks(oe,drug),correlate_ranks(kd,drug)
            value=(ro-rk)/2 if method=='paired_oe_kd' else ro if method=='oe_only_matching' else -rk
            if not math.isfinite(value):raise ValueError('A required background draw is invalid.')
            values.append(value)
        doses=[median(values[i] for i in indices) for _,indices in dose_groups]
        cells=[mean(doses[i] for i in indices) for _,indices in cell_groups]
        tissues=[mean(cells[i] for i in indices) for _,indices in tissue_groups]
        for key,indices in ingredients.items():result[key].append(mean(tissues[i] for i in indices))
        counters['aggregate_iterations']=counters.get('aggregate_iterations',0)+1
    return {key:sorted(draws) for key,draws in result.items()}
