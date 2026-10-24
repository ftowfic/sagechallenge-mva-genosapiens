"""Descriptive review scores. These functions do not change scientific decisions."""
from itertools import groupby

from genosapiens.track2_operator_inputs import require
from genosapiens.track2_report_model import presentation_policy
from genosapiens.track2_scientific_storage import Records, json_identity
from genosapiens.track2_scientific_records import ordered_records, stored_json

METHODS = ('paired_oe_kd', 'kd_only_reversal', 'oe_only_matching')
ROUTES = {'lincs': 'lincs_cell_specific_expression', 'mechanism': 'qualified_target_drug_linkage'}


def score_key(row):
    return tuple(row[k] for k in ('canonical_active_ingredient_id', 'gene_id', 'scenario_id', 'objective_id'))


def _normalized(value):
    return ' '.join(value.split()).upper()


def score_objective(rows, objective, coverage, selected_cells=()):
    """Count only retained observations within one ingredient and objective.

    Source pointers identify every row and contribution. Duplicate cells and
    mechanism records remain available without multiplying their points.
    """
    rows = list(rows)
    if rows and 'tissue_score' in rows[0]:
        from genosapiens.track2_tissue_ranking import tissue_report_score
        return tissue_report_score(rows,objective,coverage)
    require(bool(rows) and len({score_key(r) for r in rows}) == 1, 'Report scoring cannot pool different hypotheses.')
    first = rows[0]
    require(all(first[k] == objective[k] for k in ('objective_id', 'gene_id', 'scenario_id')),
            'Report contributions belong to another objective.')
    paired, single, failed, conflicts, cell_records, links, excluded = set(), set(), set(), set(), [], {}, []
    for row in rows:
        source = row['source_reference']
        for index, cell in enumerate(row['cells']):
            pointer = source+'/cells/'+str(index)
            state, method, cid = cell['state'], cell['method'], cell['cell_id']
            require(state in {'pass','fail','unassessed'} and method in METHODS, 'A report cell has an unsupported result state or method.')
            components = []
            if state == 'pass':
                if method == 'paired_oe_kd':
                    paired.add(cid); components.append('B')
                else:
                    single.add(cid); components.append('U')
            if state == 'fail':
                failed.add(cid); components.append('C')
            if cell.get('conflict'):
                conflicts.add(cid)
                if 'C' not in components: components.append('C')
            cell_records.append(dict(cell=cell, source_reference=pointer, components=components,
                                     reason='observed' if components else 'unassessed'))
        for index, record in enumerate(row['mechanisms']):
            pointer = source+'/mechanisms/'+str(index)
            targets = sorted({_normalized(t) for t in record.get('targets', [])})
            applicable_targets = {_normalized(t) for t in objective.get('target_ids', [objective['gene_id']])}
            desired = objective.get('desired_modulation', '').removeprefix('partial_')
            reported = {'reported_positive': 'increase', 'reported_negative': 'decrease'}.get(record['direction_class'])
            reason = None
            if reported is None: reason = record['direction_class']
            elif not applicable_targets.intersection(targets): reason = 'unrelated_target'
            elif record.get('group_warning') or len(targets) != 1: reason = 'unresolved_target_group'
            elif reported != desired: reason = 'opposed_or_unresolved_objective_direction'
            if reason:
                excluded.append(dict(record=record, source_reference=pointer, reason=reason))
                continue
            action = _normalized(record['action'])
            require(bool(action), 'A counted mechanism link has no action.')
            key = (first['canonical_active_ingredient_id'], objective['objective_id'], targets[0], action)
            link = links.setdefault(key, dict(canonical_active_ingredient_id=key[0], objective_id=key[1],
                target_id=key[2], action=key[3], reported_direction=reported, desired_direction=desired,
                source_record_ids=[], source_references=[], references=[]))
            if record['mechanism_id'] not in link['source_record_ids']: link['source_record_ids'].append(record['mechanism_id'])
            link['source_references'].append(pointer)
            for reference in record.get('references', []):
                if reference not in link['references']: link['references'].append(reference)
    observed = {r['cell']['cell_id'] for r in cell_records}
    for cell in selected_cells:
        if cell['cell_id'] not in observed:
            cell_records.append(dict(cell=dict(cell, state='unassessed', conflict=False),
                source_reference=cell['source_reference'], components=[], reason='local_result_unavailable'))
    counted = dict(B=sorted(paired), U=sorted(single), C=sorted(failed | conflicts))
    components = dict(I_L=int(bool(paired | single)), I_M=int(bool(links)), B=len(paired), U=len(single),
                      C=len(failed | conflicts), M=len(links))
    weights = presentation_policy()['weights']
    lincs = sum(weights[k]*components[k] for k in ('I_L','B','U','C'))
    mechanism = sum(weights[k]*components[k] for k in ('I_M','M'))
    route_states = [r for r in coverage if r['objective_id'] == objective['objective_id']]
    complete_states = {'completed_with_candidates', 'completed_no_candidates', 'not_applicable'}
    provisional = any(not any(r['route'] == route and r['state'] in complete_states for r in route_states)
                      for route in ROUTES.values())
    names = sorted({r['name'] for r in rows}, key=lambda value:(value.casefold(), value))
    return dict(**{k:first[k] for k in ('canonical_active_ingredient_id','gene_id','gene_symbol','scenario_id','objective_id')},
        name=names[0], names=names, ingredient_keys=sorted({r['ingredient_key'] for r in rows}),
        mechanism_nominated=any(r['mechanism_nominated'] for r in rows),
        lincs_nominated=any(r['lincs_nominated'] for r in rows),
        components=components, E=components['I_L']+components['I_M'],
        half_points=dict(S_L=lincs, S_M=mechanism, S=lincs+mechanism),
        scores=dict(S_L=lincs/2, S_M=mechanism/2, S=(lincs+mechanism)/2),
        counted_cells=counted, cells=cell_records,
        counted_links=[links[k] for k in sorted(links)], excluded_links=excluded,
        conflicts=sorted({value for row in rows for key in ('applicability_conflicts','target_group_warnings','condition_conflicts')
                          for value in row.get(key, [])}),
        source_references=[r['source_reference'] for r in rows],
        scientific_decisions=[dict(source_reference=r['source_reference'], **{k:r[k] for k in
            ('paired_research_eligible','gates','clinical_tier','eligibility','direction_agreement') if k in r}) for r in rows],
        coverage=route_states, provisional=provisional)


def score_rows(concordance, objectives, coverage, selected_cells=None):
    """Retain all objective scores using bounded external grouping."""
    objectives = {o['objective_id']:o for o in objectives}
    selected_cells = selected_cells or {}
    indexed = (dict(row, source_reference=row.get('source_reference', 'discovery/concordance.json#/'+str(i)))
               for i,row in enumerate(concordance))
    def scores():
        for key, rows in groupby(ordered_records(indexed, key=score_key), key=score_key):
            require(key[3] in objectives, 'A report nomination has no retained objective.')
            yield score_objective(rows, objectives[key[3]], coverage, selected_cells.get(key[3], ()))
    return stored_json(Records(scores()))


def selected_native_cells(plan):
    if plan['schema']=='genosapiens.track2-exploratory-plan/v2':
        return {}
    selections = {row['job_id']:row for row in plan['genetic_selections']}
    result = {}
    for index, binding in enumerate(plan['scenario_job_bindings']):
        selected = selections[binding['job_id']]
        result.setdefault(binding['objective_id'], []).append(dict(cell_id=selected['cell_id'], method=selected['method'],
            job_id=binding['job_id'], source_reference='plan-source.json#/exploratory_discovery/scenario_job_bindings/'+str(index)))
    return result


def select_tables(scores, *, tissue=False):
    """Choose each metric's best objective without pooling across objectives."""
    from itertools import chain
    iterator=iter(scores);first=next(iterator,None)
    rows=chain((first,),iterator) if first is not None else ()
    from genosapiens.track2_compact_ranking import select_compact_tables
    return select_compact_tables(rows,tissue=tissue or first is not None and 'tissue_score' in first)


def publish_rankings(root, scores, selections):
    """Export every objective score and every table's complete representatives."""
    from pathlib import Path
    from genosapiens.track2_scientific_storage import publish_json
    root = Path(root)
    publish_json(root/'objective-scores.json', scores)
    for table, selection in selections.items():
        publish_json(root/(table+'-ranking.json'), selection['ranking'])
