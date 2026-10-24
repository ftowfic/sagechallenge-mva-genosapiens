"""Explicit tissue contracts. Legacy cell records require their original readers."""
from functools import lru_cache
from collections.abc import Mapping
import json
import math
from pathlib import Path

from genosapiens.json_schema import validate_json_schema
from genosapiens.track2_operator_inputs import require
from genosapiens.track2_published_genotype import canonical, sha256

KINDS = ('membership', 'query', 'condition', 'result', 'score')
NORMALIZATION_POLICY = 'lincs-primary-site/case-whitespace-v1'


def schema_name(kind):
    require(kind in KINDS, 'Unknown tissue record kind.')
    return 'genosapiens.track2-tissue-'+kind+'/v1'


@lru_cache(maxsize=None)
def _schema(kind):
    schema_name(kind)
    path = Path(__file__).resolve().parents[2]/'schemas/track2'/('tissue-'+kind+'.schema.json')
    return json.loads(path.read_bytes())


def identity(record):
    from genosapiens.track2_scientific_storage import json_identity
    return json_identity(record)['sha256']


def validate_tissue_record(record, kind):
    validate_json_schema(record, _schema(kind))
    def finite(value):
        if isinstance(value, float):
            require(math.isfinite(value), 'Tissue records cannot contain nonfinite numbers.')
        elif isinstance(value, Mapping):
            for item in value.values(): finite(item)
        elif isinstance(value, list):
            for item in value: finite(item)
    finite(record)
    if kind == 'membership':
        require((record['tissue_id'] is None) == (record['unresolved_reason'] is not None),
            'Unresolved tissue membership requires a reason.')
    elif kind == 'query':
        require(record['query_sha256'] == identity({k:v for k,v in record.items() if k!='query_sha256'}),
            'The pooled query identity changed.')
        require(record['gene_count'] == sum(bool(r['selected']) for r in record['genes']),
            'The pooled query feature count changed.')
    elif kind == 'condition':
        require(len(record['gene_ids']) == len(record['values']), 'The condition features are misaligned.')
        if record['state'] == 'valid':
            require(record['r_desired'] is not None and record['reason'] is None,
                'A valid condition requires a defined comparison.')
            require(record['method'] != 'paired_oe_kd' or
                record['r_oe'] is not None and record['r_kd'] is not None,
                'A paired condition requires both correlations.')
    elif kind == 'result':
        require(record['n_conditions'] == len(record['condition_ids']) and
            record['n_valid'] + record['n_undefined'] == record['n_conditions'],
            'The tissue condition inventory differs from its counts.')
        require(record['contribution'] == (int(record['support'])-int(record['conflict'])
            if record['evaluable'] else None), 'The tissue contribution differs from its flags.')
        require(record['evaluable'] or not (record['support'] or record['conflict']),
            'Unavailable tissue evidence cannot establish paired flags.')
        require(not record['incomplete_jobs'] or record['state']=='pending' and not record['evaluable'],
            'Incomplete tissue jobs cannot establish finalized flags.')
    elif kind == 'score':
        declared = set(record['declared_tissues'])
        evaluated = set(record['evaluable_tissues'])
        require(evaluated <= declared, 'Evaluable tissues exceed the declared inventory.')
        for name, field in [('P','supporting_tissues'), ('C','conflicting_tissues'), ('N','evaluable_tissues')]:
            require(record[name] == len(record[field]) and set(record[field]) <= evaluated,
                'Tissue counts differ from their exact sets.')
        require(record['numerator'] == record['P']-record['C'] and record['denominator'] == record['N'],
            'The normalized score differs from its counts.')
    return record
