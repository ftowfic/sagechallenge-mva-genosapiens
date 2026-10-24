"""Versioned derived strength. Original scientific decisions remain external."""
from collections import defaultdict
from functools import lru_cache
import json
import math
from pathlib import Path
from statistics import mean, median

POLICY_SCHEMA = 'genosapiens.track2-derived-scoring-policy/v1'
SCORE_SCHEMA = 'genosapiens.track2-derived-strength/v1'
BACKGROUND_SCHEMA = 'genosapiens.track2-frozen-background/v1'
METHODS = ('paired_oe_kd', 'oe_only_matching', 'kd_only_reversal')


@lru_cache(maxsize=1)
def _policy_bytes():
    return (Path(__file__).resolve().parents[2]/'config/track2/derived-scoring-v1.json').read_bytes()


def policy():
    """Return a fresh policy so a consumer cannot mutate the shared contract."""
    value = json.loads(_policy_bytes())
    if value['schema'] != POLICY_SCHEMA or value['calibration']['schema'] != BACKGROUND_SCHEMA:
        raise ValueError('Unsupported derived scoring contract.')
    if tuple(value['methods']) != METHODS or value['supported_direction'] != 'increase':
        raise ValueError('Unsupported derived scoring method or direction.')
    return value


def validate_contract(record, kind):
    from genosapiens.json_schema import validate_json_schema
    if kind not in {'derived-strength', 'frozen-background'}:
        raise ValueError('Unknown derived contract.')
    schema = json.loads((Path(__file__).resolve().parents[2]/'schemas/track2'/(kind+'.schema.json')).read_bytes())
    validate_json_schema(record, schema)
    if kind == 'frozen-background':
        draws = record['draws']
        if any(not math.isfinite(x) for x in draws) or any(a>b for a,b in zip(draws,draws[1:])):
            raise ValueError('A frozen background requires sorted finite draws.')
    else:
        strength = record['strength']
        if strength is not None and not math.isfinite(strength):
            raise ValueError('Derived strength must be finite.')
        if record['n_conditions'] != record['n_valid']+record['n_undefined']:
            raise ValueError('Derived condition counts disagree.')
        if (record['direction'] != 'increase') != (record['status'] == 'unsupported_direction'):
            raise ValueError('Derived direction and status disagree.')
        if record['status'] in {'unavailable','unsupported_direction'} and strength is not None:
            raise ValueError('Unavailable evidence cannot have a strength.')
        if record['status'] == 'complete' and strength is None:
            raise ValueError('Complete evidence requires a strength.')
        if set(record['contributing_tissues']) | set(record['missing_tissues']) != set(record['declared_tissues']):
            raise ValueError('Derived coverage differs from the declared scope.')
    return record


def coefficient(value):
    if value is None:
        return None
    if isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value) or not -1 <= value <= 1:
        raise ValueError('A correlation must be finite and between -1 and 1.')
    return value


def condition_strength(condition, *, direction):
    method = condition['method']
    if method not in METHODS:
        raise ValueError('Unsupported expression method.')
    if direction != 'increase':
        return dict(status='unsupported_direction', strength=None)
    oe, kd = coefficient(condition.get('r_oe')), coefficient(condition.get('r_kd'))
    if condition['state'] != 'valid':
        return dict(status='unavailable', strength=None)
    required = (oe, kd) if method == 'paired_oe_kd' else (oe,) if method == 'oe_only_matching' else (kd,)
    if any(value is None for value in required):
        return dict(status='unavailable', strength=None)
    value = (oe-kd)/2 if method == 'paired_oe_kd' else oe if method == 'oe_only_matching' else -kd
    return dict(status='complete', strength=value)


def aggregate_strength(conditions, *, method, direction, declared_tissues, incomplete=False):
    """Deduplicate retained observations and balance dose, time, cell, and tissue."""
    if method not in METHODS:
        raise ValueError('Unsupported expression method.')
    declared = sorted(set(declared_tissues))
    groups = defaultdict(list)
    seen = {}
    unavailable = []
    scope = set()
    for source in conditions:
        if source['method'] != method or source['tissue_id'] not in declared:
            raise ValueError('A condition differs from the declared method or tissue scope.')
        scope.add((source['gene_id'], source['medication_entity_key']))
        if len(scope) != 1:
            raise ValueError('An aggregate cannot combine genes or ingredients.')
        key = (source['query_sha256'], source['sig_id'])
        # Objective labels do not define a new experimental observation.
        observation = {k:v for k,v in source.items() if k not in {'scenario_id','objective_id'}}
        serialized = json.dumps(observation, sort_keys=True, separators=(',', ':'), allow_nan=False)
        if key in seen:
            if seen[key] != serialized:
                raise ValueError('Repeated condition signatures contain conflicting evidence.')
            continue
        seen[key] = serialized
        result = condition_strength(source, direction=direction)
        if result['strength'] is None:
            unavailable.append(dict(query_sha256=key[0], sig_id=key[1], status=result['status']))
        else:
            if source['pert_time_unit'] != 'h':
                raise ValueError('A retained condition requires time in hours.')
            groups[(source['tissue_id'], source['cell_id'], source['pert_time'])].append(result['strength'])
    cells = defaultdict(list)
    for (tissue, cell, time), values in sorted(groups.items()):
        cells[(tissue, cell)].append(median(values))
    tissues = defaultdict(list)
    for (tissue, cell), times in sorted(cells.items()):
        tissues[tissue].append(mean(times))
    values = {tissue:mean(cells) for tissue,cells in sorted(tissues.items())}
    strength = mean(values.values()) if values else None
    status = ('unsupported_direction' if direction != 'increase' else 'provisional' if incomplete
              else 'complete' if strength is not None else 'unavailable')
    record = dict(schema=SCORE_SCHEMA, policy_id=policy()['policy_id'], method=method, direction=direction,
        status=status, strength=strength, declared_tissues=declared,
        contributing_tissues=sorted(tissues), missing_tissues=sorted(set(declared)-set(tissues)),
        tissue_strengths=values, n_conditions=len(seen), n_valid=sum(map(len,groups.values())),
        n_undefined=len(unavailable), n_cells=len(cells), n_times=len(groups),
        unavailable_conditions=sorted(unavailable,key=lambda r:(r['query_sha256'],r['sig_id'])))
    return validate_contract(record, 'derived-strength')
