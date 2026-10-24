"""Separate percentile-first ranks without changing source nominations."""
from fractions import Fraction
import math

from genosapiens.track2_derived_strength import METHODS, policy


def exact_pair(row):
    numerator,denominator=row['numerator'],row['denominator']
    if (type(numerator) is not int or type(denominator) is not int or denominator!=20000
            or not 0<=numerator<=denominator):
        raise ValueError('A calibrated rank requires exact 10,000-draw percentile counts.')
    strength=row['strength']
    if isinstance(strength,bool) or not isinstance(strength,(int,float)) or not math.isfinite(strength) or not -1<=strength<=1:
        raise ValueError('A calibrated rank requires finite raw strength.')
    if not row.get('entry_key') or not row.get('distribution_sha256'):
        raise ValueError('A calibrated rank requires its frozen background identity.')
    return Fraction(numerator,denominator),strength


def rank_methods(rows, *, source_selected, nominated):
    """Keep original ingredient/objective selections and all unranked evidence.

    source_selected maps (ingredient, method) to the source-selected objective.
    nominated contains original (ingredient, objective, method) nominations.
    Complete source and alternative objective rows remain in the full export.
    """
    from genosapiens.track2_compact_ranking import CompactRankStore
    from genosapiens.track2_scientific_records import stored_json
    from genosapiens.track2_scientific_storage import Records
    nominated=set(nominated)
    with CompactRankStore() as store:
        for source in rows:
            method=source['method']
            if method not in METHODS:raise ValueError('Unknown derived ranking method.')
            key=(source['medication_entity_key'],source['objective_id'],method)
            row=dict(source,derived_ranking_policy=policy()['policy_id'],method_rank=None,
                source_nominated=key in nominated,
                source_representative=source_selected.get((key[0],method))==key[1])
            record=store.add(row,unique=key)
            if row['status']!='complete' or row.get('numerator') is None:continue
            pair=exact_pair(row)
            if row['source_nominated'] and row['source_representative']:
                score=(-pair[0],-pair[1])
                store.key(method,key[0],record,score,(*score,row.get('name','').casefold(),key[0],key[1]))
        store.finish(representatives=False)
        return stored_json(dict(schema='genosapiens.track2-derived-rankings/v1',
            methods={method:Records(store.rows(method,derived=True)) for method in METHODS},
            unranked=Records(row for row in store.complete() if row['status']!='complete' or row.get('numerator') is None),
            complete_export=Records(store.complete())))
