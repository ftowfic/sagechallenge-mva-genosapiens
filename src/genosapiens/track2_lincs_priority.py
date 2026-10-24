"""Descriptive review ordering over unchanged, independently scored cell queries."""
import math


CAVEAT = (
    '“Repeated agreement with no observed conflict” qualifies a medication for research under the selected LINCS policy. '
    'Expression agreement alone permits this eligibility. Separate assessments cover mechanism, drug action, achievable exposure, safety and validation. '
    'This eligibility does not establish clinical benefit or a case match.'
)
ELIGIBILITY_POLICY = {
    'version': 'lincs-research-eligibility/v1',
    'basis': 'lincs_repeated_consistent',
    'agreement_group': 'repeated_consistent',
    'minimum_passing_cells': 2,
    'failed_cells_allowed': 0,
    'condition_conflicts_allowed': 0,
    'unassessed_cells': 'disclose_without_imputing_failure',
    'scientific_gate_or_tier_override': False,
    'drug_specific_review_required_for_this_basis': False,
}
PRIORITY_POLICY = {
    'version': 'lincs-review-priority/v1',
    'groups': ['repeated_consistent', 'repeated_mixed_conditions', 'mixed_cells',
               'single_cell', 'no_local_support', 'unassessed'],
    'within_group': ['passing_cells_desc', 'joint_priority_desc', 'desired_direction_cells_desc',
                     'desired_priority_desc', 'conflicting_cells_asc', 'medication_entity_key_asc'],
    'joint_priority': 'number_of_passing_cells * sum(finite_local_joint_scores_including_negative)',
    'desired_priority': 'number_of_cells_with_r_desired_gt_0 * sum(finite_local_r_desired)',
    'missing': 'omit from sums; retain missing and tested denominators',
    'ties': 'same review rank for identical sorting components; stable ingredient display',
    'interpretation': 'Review ordering only; not a pooled correlation, efficacy estimate or scientific tier.',
}
GROUPS = {
    'repeated_consistent': ('🟢', 'Repeated agreement with no observed conflict'),
    'repeated_mixed_conditions': ('🟡', 'Repeated agreement with some conflicting conditions'),
    'mixed_cells': ('🟠', 'Cell lines disagree'),
    'single_cell': ('🔵', 'Support in one cell line only'),
    'no_local_support': ('🔴', 'No cell line passes both checks'),
    'unassessed': ('⚪', 'Not assessed'),
}


def finite(value):
    return value is not None and math.isfinite(value)


def review_priority(result):
    """Include every tested identity; failed/absent evidence cannot disappear."""
    output = []
    for source in result['cross_cell']:
        cells = source['cells']
        passing = [c for c, r in cells.items() if r['state'] == 'pass']
        failed = [c for c, r in cells.items() if r['state'] == 'fail']
        conflicts = [c for c, r in cells.items() if r.get('conflict', False)]
        finite_cells = [c for c, r in cells.items() if finite(r['joint'])]
        desired = [c for c, r in cells.items() if finite(r['r_desired']) and r['r_desired'] > 0]
        opposed = [c for c, r in cells.items() if finite(r['r_desired']) and r['r_desired'] < 0]
        if len(passing) >= 2 and not failed and not conflicts:
            group = 'repeated_consistent'
        elif len(passing) >= 2 and not failed:
            group = 'repeated_mixed_conditions'
        elif passing and failed:
            group = 'mixed_cells'
        elif passing:
            group = 'single_cell'
        elif finite_cells:
            group = 'no_local_support'
        else:
            group = 'unassessed'
        joint_sum = math.fsum(r['joint'] for r in cells.values() if finite(r['joint']))
        desired_sum = math.fsum(r['r_desired'] for r in cells.values() if finite(r['r_desired']))
        output.append({
            **source, 'name': next((r['name'] for r in cells.values() if r.get('name')), source['medication_entity_key']),
            'agreement_group': group, 'passing_cells': len(passing), 'desired_direction_cells': len(desired),
            'opposed_direction_cells': len(opposed), 'finite_score_cells': len(finite_cells),
            'tested_cell_count': sum(r['n_conditions'] > 0 for r in cells.values()),
            'selected_cell_count': len(cells), 'missing_cells': [c for c, r in cells.items() if r['state'] == 'unassessed'],
            'failed_cells': failed, 'conflicting_cells': conflicts,
            'joint_score_sum': joint_sum if finite_cells else None,
            'joint_priority': len(passing) * joint_sum if finite_cells else None,
            'desired_priority': len(desired) * desired_sum if any(finite(r['r_desired']) for r in cells.values()) else None,
            'scientific_interpretation': PRIORITY_POLICY['interpretation'],
        })

    def key(row):
        # With zero passes the proposed product is zero. Use desired-direction
        # count/score to distinguish partial alignment, neutral and opposing rows.
        return (PRIORITY_POLICY['groups'].index(row['agreement_group']), -row['passing_cells'],
                -(row['joint_priority'] or 0), -row['desired_direction_cells'],
                -(row['desired_priority'] or 0), len(row['conflicting_cells']))
    output.sort(key=lambda r: (*key(r), r['medication_entity_key']))
    for i, row in enumerate(output):
        row['review_rank'] = output[i-1]['review_rank'] if i and key(row) == key(output[i-1]) else i+1
    return output


def md(value):
    return str(value).replace('|', '\\|').replace('\n', ' ')


def score(value):
    return '—' if value is None else f'{value:+.3f}'


def priority_table(rows, cells):
    labels = [f"{cell} ({rows[0]['cells'][cell]['gene_count']} genes)" if rows else cell for cell in cells]
    lines = ['| Rank | Medication | Agreement | Research eligibility | Pass / tested | Direction + | Review score | '
             + ' | '.join(labels) + ' |',
             '|---:|---|---|---|---|---:|---:|' + '---|' * len(cells)]
    for row in rows:
        color, label = GROUPS[row['agreement_group']]
        local = []
        for cell in cells:
            r = row['cells'][cell]
            if r['state'] == 'unassessed':
                local.append('⚪ Not tested' if not r['n_conditions'] else '⚪ Undefined')
            else:
                symbol = '✓' if r['state'] == 'pass' else '×'
                suffix = ' (mixed conditions)' if r.get('conflict') else ''
                local.append(f"{symbol} {score(r['joint'])}{suffix}")
        eligibility = ('Eligible — LINCS agreement' if row.get('category_eligibility', {}).get('eligible')
                       else 'Eligible — reviewed evidence' if row.get('research_eligible')
                       else 'Not eligible on current evidence')
        lines.append(f"| {row['review_rank']} | {md(row['name'])} | {color} {label} | {eligibility} | "
                     f"{row['passing_cells']} / {row['tested_cell_count']} of {row['selected_cell_count']} | "
                     f"{row['desired_direction_cells']} | {score(row['joint_priority'])} | " + ' | '.join(local) + ' |')
    return '\n'.join(lines)


def priority_explanation():
    return ('The table orders agreement groups as follows:\n\n'
            '1. Repeated agreement with no observed conflict.\n'
            '2. Repeated agreement with some conflicting conditions.\n'
            '3. Disagreement between cell lines.\n'
            '4. Support in one cell line only.\n'
            '5. No cell line passes both checks.\n'
            '6. Not assessed.\n\n'
            'Within each group, the table uses these measures in order:\n\n'
            '1. More cell lines that pass both checks.\n'
            '2. Higher review score.\n'
            '3. More cell lines with positive correlation to the desired pattern.\n'
            '4. Higher score for the desired direction.\n'
            '5. Fewer cell lines with conflicting conditions.\n\n'
            'Ties share a rank.\n\n'
            '**Review score = passing cell lines × sum of all finite local joint scores**, including negative scores. '
            '“Direction +” counts positive correlations to the desired pattern. This alone cannot pass the rule that requires both checks. '
            'The score for the desired direction equals this count × sum of finite correlations to the desired pattern.\n\n'
            'The sums exclude missing cells. Coverage counts retain missing cells, which do not mean failed results. '
            'Cell columns show the local joint score. ✓ means both checks pass. × means a failure.\n\n'
            'Colors describe observed agreement, not safety or benefit. Stronger correlations are not stronger clinical evidence. '
            'HA1E and VCAP use different gene sets. This order does not represent pooled correlation, independent replication, statistical significance or therapeutic rank.')
