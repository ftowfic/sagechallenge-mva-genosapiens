"""Medication assessment completion, independent of candidate eligibility.

Scope records precede results. Route records come from validated portfolios,
including source-bound frozen searches. Review acceptance remains separate.
"""
from dataclasses import asdict, dataclass
from typing import Mapping


@dataclass(frozen=True, slots=True)
class AssessmentReason:
    code: str
    scenario_id: str | None = None
    objective_id: str | None = None
    route: str | None = None


@dataclass(frozen=True, slots=True)
class AssessmentStatus:
    status: str
    reasons: tuple[AssessmentReason, ...]
    scenarios: tuple[dict, ...]
    objectives: tuple[dict, ...]

    @property
    def complete(self):
        return self.status == 'complete'

    def as_dict(self):
        return asdict(self)


def route_evidence(records, selected_routes, *, native=False):
    """Retain completed searches, failed attempts and exclusions separately."""
    if native:
        return tuple({**row,'required':row['route'] in selected_routes} for row in records)
    return tuple({
        'route': row.route, 'required': row.route in selected_routes,
        'state': row.attempt.state if row.attempt else 'not_applicable',
        'executed': bool(row.attempt and row.attempt.executed),
        'query_ids': list(row.attempt.query_ids) if row.attempt else [],
        'evidence_ids': list(row.evidence_ids),
        'review_record_id': row.review_record_id,
        'rationale': row.biological_rationale,
    } for row in sorted(records, key=lambda item: item.route))


def assess_objective(objective_id, scenario_id, supported, routes, *, native=False):
    if native:
        required=[r for r in routes if r['required']]
        complete=[r for r in required if r['state'] in {'completed_with_candidates','completed_no_candidates'} and r['executed'] and r['query_ids']]
        failures=[r for r in required if r['state'] not in {'completed_with_candidates','completed_no_candidates','unavailable'}
            or r['state']=='unavailable' and not r.get('reason_codes_result')]
        reasons=[]
        if not supported:reasons.append(AssessmentReason('unsupported_native_hypothesis',scenario_id,objective_id))
        if not complete:reasons.append(AssessmentReason('no_completed_applicable_search',scenario_id,objective_id))
        reasons.extend(AssessmentReason('required_discovery_incomplete',scenario_id,objective_id,r['route']) for r in failures)
        status='complete' if not reasons else 'incomplete' if failures or any(r['executed'] for r in required) else 'not_performed'
        return dict(objective_id=objective_id,scenario_id=scenario_id,supported=supported,status=status,
            reasons=[asdict(r) for r in reasons],routes=list(routes))
    reasons = []
    required = [r for r in routes if r['required']]
    completed = [r for r in required if r['executed'] and r['query_ids']
                 and r['state'] in {'attempted_complete', 'no_result'}]
    started = supported and any(r['executed'] for r in required)
    if not supported:
        reasons.append(AssessmentReason('unsupported_objective', scenario_id, objective_id))
    if not completed:
        reasons.append(AssessmentReason('no_completed_applicable_search', scenario_id, objective_id))
    for row in required:
        if row['state'] not in {'attempted_complete', 'no_result', 'not_applicable'}:
            reasons.append(AssessmentReason('required_route_' + row['state'], scenario_id, objective_id, row['route']))
    status = 'complete' if not reasons else 'incomplete' if started else 'not_performed'
    return {'objective_id': objective_id, 'scenario_id': scenario_id, 'supported': supported,
            'status': status, 'reasons': [asdict(r) for r in reasons], 'routes': list(routes)}


def assess_scope(scenarios, objectives, *, native=False):
    """Derive status from required scenarios and supported objective searches.

    An excluded scenario requires a prospective scope rationale. Optional routes
    stay visible but do not impose a universal discovery or cell-line gate.
    """
    scenarios = tuple(sorted(scenarios, key=lambda r: r['scenario_id']))
    objectives = tuple(sorted(objectives, key=lambda r: r['objective_id']))
    by_scenario = {s['scenario_id']: s for s in scenarios}
    if len(by_scenario) != len(scenarios) or len({o['objective_id'] for o in objectives}) != len(objectives):
        raise ValueError('Assessment scope contains duplicate identities.')
    if any(o['scenario_id'] not in by_scenario for o in objectives):
        raise ValueError('Assessment objective has no source scenario.')
    if any(not s['required'] and not s['rationale'] for s in scenarios):
        raise ValueError('Excluded scenarios require a prospective scope rationale.')
    if any(not by_scenario[o['scenario_id']]['required'] for o in objectives):
        raise ValueError('Excluded scenarios cannot acquire assessment objectives.')
    reasons = []
    required = [s for s in scenarios if s['required']]
    if not required:
        reasons.append(AssessmentReason('no_eligible_exploratory_objectives' if native else 'no_required_scenarios'))
    if not objectives:
        reasons.append(AssessmentReason('no_supported_objectives'))
    for scenario in required:
        sid = scenario['scenario_id']
        if not scenario['supported']:
            reasons.append(AssessmentReason('unresolved_required_scenario', sid))
        if not any(o['scenario_id'] == sid for o in objectives):
            reasons.append(AssessmentReason('scenario_without_objective', sid))
    replayed = []
    for objective in objectives:
        row = assess_objective(objective['objective_id'], objective['scenario_id'],
            by_scenario[objective['scenario_id']]['supported'] and objective['supported'], objective['routes'],native=native)
        replayed.append(row)
        reasons.extend(AssessmentReason(**r) for r in row['reasons'])
    started = any(o['status'] != 'not_performed' for o in replayed)
    status = 'complete' if not reasons else 'incomplete' if started else 'not_performed'
    return AssessmentStatus(status, tuple(reasons), scenarios, tuple(replayed))


def assess_plan(mechanisms: Mapping, objective_scenarios: Mapping, portfolios: Mapping, scope=None, *, native_plan=None,discovery=None):
    if native_plan is not None:
        return assess_native_plan(native_plan,discovery)
    if scope is not None and set(scope)!=set(mechanisms):
        raise ValueError('Assessment scenario scope differs from retained mechanisms.')
    if scope is not None and any(set(row)!={'required','rationale'} or type(row['required']) is not bool
                                 or not isinstance(row['rationale'],str) or not row['rationale'].strip()
                                 for row in scope.values()):
        raise ValueError('Assessment scope requires explicit inclusion and a rationale.')
    scenarios = [{'scenario_id': sid, 'required': True, 'supported': ready.permits_planning,
                  'rationale': 'All retained mechanism scenarios are in the assessment scope.',
                  **({} if scope is None else scope[sid])}
                 for sid, ready in mechanisms.items()]
    objectives = [assess_objective(oid, objective_scenarios[oid],
        mechanisms[objective_scenarios[oid]].permits_planning,
        route_evidence(portfolio.assessments, portfolio.selected_routes))
        for oid, portfolio in portfolios.items()]
    return assess_scope(scenarios, objectives)


def assess_native_plan(plan,discovery):
    nested=plan.get('exploratory_discovery',plan)
    scenarios={o['scenario_id']:dict(scenario_id=o['scenario_id'],required=True,supported=True,
        rationale='Accepted native restoration hypothesis.') for o in nested['objectives']}
    scenarios.update({s['scenario_id']:dict(scenario_id=s['scenario_id'],required=False,supported=False,
        rationale=', '.join(s['reasons'])) for s in nested['excluded_scenarios']})
    outcomes={(r['objective_id'],r['route']):r for r in (discovery or {}).get('route_outcomes',[])}
    objectives=[]
    for objective in nested['objectives']:
        routes=[]
        for planned in nested['routes']:
            if planned['objective_id']!=objective['objective_id']:continue
            row=outcomes.get((objective['objective_id'],planned['route']))
            if row is None:
                row={**planned,'state':'failed' if planned['planned_state']=='executable' else 'unavailable',
                    'executed':False,'query_ids':[],'reason_codes_result':planned['reason_codes'] or ['completion_missing']}
            routes.append({**row,'required':True})
        objectives.append(assess_objective(objective['objective_id'],objective['scenario_id'],True,routes,native=True))
    return assess_scope(scenarios.values(),objectives,native=True)
