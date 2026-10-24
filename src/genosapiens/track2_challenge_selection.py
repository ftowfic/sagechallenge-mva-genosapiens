"""One bounded table, with independent expression methods and retained source ranks."""
from genosapiens.track2_derived_ranking import exact_pair
from genosapiens.track2_derived_strength import METHODS
from genosapiens.track2_operator_inputs import require


def _key(row):
    return row['canonical_active_ingredient_id'], row['objective_id']


def _compact(row, pointer, *, group, derived=None):
    value = dict(group=group, canonical_active_ingredient_id=row['canonical_active_ingredient_id'],
        name=row['name'], objective_id=row['objective_id'], scenario_id=row['scenario_id'],
        gene_id=row['gene_id'], gene_symbol=row.get('gene_symbol'), source_reference=pointer,
        historical_rank=row['rank'], historical_metric=row['metric'],
        historical_scores=dict(row['scores']), historical_components=dict(row['components']),
        provisional=row['provisional'], conflict_count=len(row['conflicts']),
        conflicts=list(row['conflicts']),
        excluded_link_count=len(row['excluded_links']), mechanism=None, expression=None)
    value['mechanism_actions'] = [dict(target=link['target_id'], action=link['action'])
        for link in row.get('counted_links', [])]
    value['excluded_mechanism_reasons'] = sorted({link['reason'] for link in row['excluded_links'] if isinstance(link,dict)})
    value['blocking_findings'] = []
    for decision in row.get('scientific_decisions', []):
        for gate in decision.get('gates', []):
            if gate['gate_id']=='G6' and gate['outcome']=='contraindicating':
                value['blocking_findings'].append(dict(kind='blocking_hazard',reasons=list(gate['reason_codes'])))
            elif gate['gate_id']=='G3' and gate['outcome']=='fail':
                value['blocking_findings'].append(dict(kind='action_conflict',reasons=list(gate['reason_codes'])))
    if row['components']['M']:
        value['mechanism'] = dict(component=row['scores']['S_M'], historical_rank=None,
            counted_link_count=row['components']['M'], source_reference=pointer)
    if derived is not None:
        value['expression'] = {k:derived[k] for k in ('method','strength','numerator','denominator',
            'entry_key','distribution_sha256','method_rank','tie_count','status')}
        value['expression']['degenerate_background'] = derived.get('degenerate_background', False)
        value['expression']['coverage'] = {k:derived[k] for k in ('declared_tissues','contributing_tissues',
            'missing_tissues','n_conditions','n_valid','n_undefined','n_cells','n_times')}
    return value


def select_evidence(historical, *, derived=None):
    """Select five per method and mechanism, preserving source representatives.

    The caller retains complete exports. Only compact selected rows enter the model.
    Legacy cell reports use their original expression selection and rank labels.
    """
    selected = []
    if derived is None:
        for index, row in enumerate(historical['lincs']['rows']):
            require('tissue_score' not in row, 'Tissue candidates require explicit derived ranking status.')
            selected.append(_compact(row, 'supplements/lincs-ranking.json#/'+str(index), group='legacy_cell'))
    else:
        wanted = {}
        for method in METHODS:
            group = derived['methods'][method]
            for index, row in enumerate(group[:5]):
                require(row['method'] == method and row['status'] == 'complete'
                        and row['source_nominated'] and row['source_representative'],
                        'A challenge expression row is not a complete source nomination.')
                exact_pair(row)
                key = row['medication_entity_key'], row['objective_id']
                wanted.setdefault(key, []).append((method, index, row))
        matched = {}
        # Keep at most fifteen source representatives while scanning the retained ranking.
        for index, row in enumerate(historical['lincs']['ranking']):
            tissue = row.get('tissue_score')
            if not tissue:continue
            key = tissue['medication_entity_key'], row['objective_id']
            if key not in wanted:continue
            require(key not in matched, 'A source representative is duplicated.')
            matched[key] = row, index
        require(set(matched) == set(wanted), 'A derived row changed the source-selected objective.')
        for method in METHODS:
            for derived_row in derived['methods'][method][:5]:
                key = derived_row['medication_entity_key'], derived_row['objective_id']
                row, index = matched[key]
                selected.append(_compact(row, 'supplements/lincs-ranking.json#/'+str(index),
                    group=method, derived=derived_row))
    require(len(selected) <= 15, 'The expression table exceeds its group limits.')
    for index, row in enumerate(historical['mechanism']['rows']):
        require(index < 5, 'The source mechanism selection exceeds five rows.')
        mechanism = dict(component=row['scores']['S_M'], historical_rank=row['rank'],
            counted_link_count=row['components']['M'], source_reference='supplements/mechanism-ranking.json#/'+str(index))
        matches = [r for r in selected if _key(r) == _key(row)]
        if matches:
            # Retain each expression score. A shared mechanism does not join methods.
            for match in matches:match['mechanism'] = dict(mechanism)
        else:
            value = _compact(row, mechanism['source_reference'], group='mechanism')
            value['mechanism'] = mechanism
            selected.append(value)
    require(len({_key(row)+(row['group'],) for row in selected}) == len(selected),
            'A challenge evidence group repeats an ingredient and objective.')
    return selected


def check_selection(rows, *, has_derived_policy):
    """Check bounded display records without accepting them as source science."""
    from collections import Counter
    import math
    from genosapiens.track2_report_model import relative_destination
    counts = Counter(); seen = set(); order = []
    group_order = {name:i for i,name in enumerate((*METHODS, 'legacy_cell', 'mechanism'))}
    for row in rows:
        group = row['group']
        require(group in group_order, 'A challenge evidence group is unsupported.')
        counts[group] += 1; order.append(group_order[group])
        key = _key(row)+(group,)
        require(key not in seen and all(row[k] for k in ('objective_id','scenario_id','gene_id')),
                'A challenge row has a duplicated or missing hypothesis identity.')
        seen.add(key)
        relative_destination(row['source_reference'])
        require(row['conflict_count'] == len(row['conflicts']), 'A challenge conflict count differs.')
        for value in row['historical_scores'].values():
            require(value is None or type(value) in (int,float) and math.isfinite(value),
                    'A historical score is invalid.')
        expression = row['expression']
        if group in METHODS:
            require(has_derived_policy and expression is not None and expression['method'] == group
                    and expression['status'] == 'complete', 'An expression row lacks its derived policy or method.')
            exact_pair(expression)
            require(all(type(expression[k]) is int and expression[k] >= 1 for k in ('method_rank','tie_count')),
                    'A challenge method rank or tie count is invalid.')
            coverage = expression['coverage']
            require(all(type(coverage[k]) is int and coverage[k] >= 0 for k in
                        ('n_conditions','n_valid','n_undefined','n_cells','n_times'))
                    and coverage['n_conditions'] == coverage['n_valid']+coverage['n_undefined'],
                    'A challenge condition count is invalid.')
        else:
            require(expression is None, 'A legacy or mechanism row cannot invent a derived score.')
        if row['mechanism'] is not None:
            relative_destination(row['mechanism']['source_reference'])
    require(order == sorted(order) and all(n <= 5 for n in counts.values()),
            'Challenge groups exceed their limits or change navigation order.')
    require(not (counts['legacy_cell'] and any(counts[m] for m in METHODS)),
            'A challenge table cannot relabel legacy cell evidence as tissue evidence.')
    return rows
