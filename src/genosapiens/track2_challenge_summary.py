"""Compact source projections that preserve absent and unfinished evidence."""
from collections import Counter

from genosapiens.track2_operator_inputs import require
from genosapiens.track2_report_ranking import ROUTES

COMPLETE = {'completed_with_candidates', 'completed_no_candidates', 'not_applicable'}


def route_category(state):
    return {'completed_with_candidates': 'positive', 'completed_no_candidates': 'negative',
            'not_applicable': 'not_applicable', 'unavailable': 'unavailable',
            'failed': 'failed'}.get(state, 'incomplete')


def _ties(rows, rank_field):
    counts = Counter(row[rank_field] for row in rows if row.get(rank_field) is not None)
    return [dict(rank=rank, count=count) for rank, count in sorted(counts.items()) if count > 1]


def summarize(context, coverage, historical, *, derived=None, route_names=None, required_routes=None):
    """Project admitted records without loading condition vectors or copying scores.

    References are JSON pointers into the artifacts bound by the presentation.
    A missing route is unavailable. No candidate count implies a completed search.
    """
    route_names=tuple(ROUTES.values()) if route_names is None else tuple(route_names)
    required_routes=set(route_names if required_routes is None else required_routes)
    require(required_routes <= set(route_names), 'The summary selected routes exceed the original policy.')
    hypotheses = []; objectives = {}; routes = {}; totals = {r:Counter() for r in route_names}
    for index, hypothesis in enumerate(context['hypotheses']):
        local = []
        for j, objective in enumerate(hypothesis['objectives']):
            oid = objective['objective_id']
            require(oid not in objectives, 'A summary objective occurs more than once.')
            require(objective['scenario_id'] == hypothesis['scenario_id'], 'A summary objective has another hypothesis.')
            objectives[oid] = dict(objective_id=oid, scenario_id=hypothesis['scenario_id'],
                gene_id=objective['gene_id'], desired_modulation=objective.get('desired_modulation'),
                source_pointer='/hypotheses/'+str(index)+'/objectives/'+str(j), routes=[])
            local.append(oid)
        hypotheses.append(dict(scenario_id=hypothesis['scenario_id'], label=hypothesis['label'],
            genes=[dict(gene_id=g['gene_id'], symbol=g.get('symbol'), metadata_state=g['metadata_state'])
                   for g in hypothesis['genes']], objective_ids=local,
            blocker_count=len(hypothesis.get('mechanism_blockers', [])),
            source_pointer='/hypotheses/'+str(index)))
    require(len({h['scenario_id'] for h in hypotheses}) == len(hypotheses), 'A summary hypothesis is duplicated.')
    for index, row in enumerate(coverage):
        key = row['objective_id'], row['route']
        require(key[0] in objectives and key[1] in route_names, 'A summary route has an unknown objective or method.')
        require(key not in routes, 'A summary route occurs more than once.')
        state = row['state']
        routes[key] = dict(route=key[1], state=state, category=route_category(state) if key[1] in required_routes else 'outside_scope',
                           source_pointer='/'+str(index))
    for oid, objective in objectives.items():
        for route in route_names:
            value = routes.get((oid, route), dict(route=route, state='missing_record',
                category='unavailable' if route in required_routes else 'outside_scope', source_pointer=None))
            value['required']=route in required_routes
            objective['routes'].append(value)
            totals[route][value['category']] += 1
        objective['processing_complete'] = bool(required_routes) and all(r['state'] in COMPLETE for r in objective['routes'] if r['required'])
    historical_ties = {method:_ties(selection['ranking'], 'rank') for method, selection in historical.items()}
    method_summaries = []
    if derived is not None:
        from genosapiens.track2_derived_strength import METHODS
        for method in METHODS:
            ranked = derived['methods'][method]
            unranked = Counter()
            coverage_counts = Counter()
            for row in derived['complete_export']:
                if row['method'] != method or not row['source_nominated'] or not row['source_representative']:
                    continue
                state = row['status']
                if state != 'complete' or row.get('numerator') is None:
                    unranked['uncalibrated' if state == 'complete' else state] += 1
                # Counts describe objective/ingredient summaries, never independent experiments.
                if row.get('missing_tissues'):coverage_counts['with_missing_tissues'] += 1
                if row.get('n_undefined', 0):coverage_counts['with_undefined_conditions'] += 1
            method_summaries.append(dict(method=method, ranked_candidates=len(ranked),
                unranked_candidates=dict(sorted(unranked.items())), ties=_ties(ranked, 'method_rank'),
                coverage_limits=dict(sorted(coverage_counts.items()))))
    return dict(schema='genosapiens.track2-challenge-summary/v1', hypotheses=hypotheses,
        objectives=list(objectives.values()), hypothesis_count=len(hypotheses),
        unresolved_hypothesis_count=sum(not h['objective_ids'] for h in hypotheses),
        objective_count=len(objectives), completed_objective_count=sum(o['processing_complete'] for o in objectives.values()),
        route_totals={route:dict(sorted(counts.items())) for route,counts in totals.items()},
        historical_ties=historical_ties, derived_methods=method_summaries)
