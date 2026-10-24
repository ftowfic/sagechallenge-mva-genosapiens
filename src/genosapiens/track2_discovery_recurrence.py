"""Report views over retained objective evidence, without scientific recomputation."""
from genosapiens.track2_lincs_priority import PRIORITY_POLICY,finite,review_priority


METHODS = ('paired_oe_kd','kd_only_reversal','oe_only_matching')
VIEWS = {
    'paired_oe_kd': 'Paired recurrence',
    'kd_only_reversal': 'KD-only recurrence',
    'oe_only_matching': 'OE-only recurrence',
    'mixed_single_direction': 'Mixed single-direction recurrence',
}
GROUP_LABELS = {
    'repeated_consistent': 'Repeated agreement without observed conflicts in this view',
    'repeated_mixed_conditions': 'Repeated agreement with condition conflicts',
    'mixed_cells': 'Passing and failed cells',
    'single_cell': 'One passing cell',
    'no_local_support': 'No passing cells',
    'unassessed': 'Unassessed',
}
FALLBACK_ORDER = ['agreement_group','passing_cells_desc','ingredient_key_asc']


def _fallback_row(source):
    cells=source['cells']
    passing=[c for c,r in cells.items() if r['state']=='pass']
    failed=[c for c,r in cells.items() if r['state']=='fail']
    conflicts=[c for c,r in cells.items() if r.get('conflict')]
    if len(passing)>=2 and not failed and not conflicts:group='repeated_consistent'
    elif len(passing)>=2 and not failed:group='repeated_mixed_conditions'
    elif passing and failed:group='mixed_cells'
    elif passing:group='single_cell'
    elif any(finite(r['r_desired']) for r in cells.values()):group='no_local_support'
    else:group='unassessed'
    return dict(source,agreement_group=group,passing_cells=len(passing),failed_cells=failed,
        conflicting_cells=conflicts,selected_cell_count=len(cells),
        tested_cell_count=sum(r['n_conditions']>0 for r in cells.values()),
        missing_cells=[c for c,r in cells.items() if r['state']=='unassessed'])


def recurrence_views(plan,concordance):
    """Keep every nomination in its own gene, scenario, objective, and method."""
    if plan['schema']=='genosapiens.track2-exploratory-plan/v2':
        from genosapiens.track2_scientific_records import stored_json
        from genosapiens.track2_scientific_storage import Records
        return stored_json(dict(schema='genosapiens.track2-tissue-recurrence/v1',
            legacy_research_eligibility='not_applicable',scopes=Records(dict(
                **{k:o[k] for k in ('gene_id','gene_symbol','scenario_id','objective_id')},
                views={view:[] for view in VIEWS},
                tissue_scores=Records(r['tissue_score'] for r in concordance if r['objective_id']==o['objective_id']))
                for o in plan['objectives'])))
    selections={r['job_id']:r for r in plan['genetic_selections']}
    selected={}
    for binding in plan['scenario_job_bindings']:
        row=selections[binding['job_id']]
        selected.setdefault(binding['objective_id'],{})[row['cell_id']]=row['method']
    from genosapiens.track2_scientific_records import stored_json
    from genosapiens.track2_scientific_storage import Records
    def scopes():
        for objective in plan['objectives']:
            oid=objective['objective_id'];views={view:[] for view in VIEWS}
            for candidate in concordance:
                if candidate['objective_id']!=oid:continue
                observed={r['cell_id']:r for r in candidate['cells']}
                cells={cell:dict(observed[cell]) if cell in observed else dict(cell_id=cell,method=method,
                    state='unassessed',n_conditions=0,r_desired=None,r_opposite=None,joint=None,conflict=False,
                    reason='local_result_unavailable') for cell,method in sorted(selected.get(oid,{}).items())}
                common={k:candidate[k] for k in ['name','ingredient_key','canonical_active_ingredient_id',
                    'gene_id','gene_symbol','scenario_id','objective_id','fallback_recurrence']}
                common.update(medication_entity_key=candidate['ingredient_key'],
                    objective_failed_cells=candidate['failed_cells'],
                    objective_condition_conflicts=candidate['condition_conflicts'])
                for view in VIEWS:
                    methods=METHODS[1:] if view=='mixed_single_direction' else (view,)
                    local={c:r for c,r in cells.items() if r['method'] in methods}
                    if not local:continue
                    if view=='mixed_single_direction' and not all(
                            any(r['method']==method and r['state']=='pass' for r in local.values()) for method in methods):continue
                    views[view].append(dict(common,cells=local,
                        paired_research_eligible=candidate['paired_research_eligible'] if view=='paired_oe_kd' else False))
            for view,rows in views.items():
                if view=='paired_oe_kd':
                    names={r['ingredient_key']:r['name'] for r in rows}
                    views[view]=[dict(r,name=names[r['ingredient_key']]) for r in review_priority(dict(cross_cell=rows))]
                else:
                    rows=[_fallback_row(row) for row in rows]
                    rows.sort(key=lambda r:(PRIORITY_POLICY['groups'].index(r['agreement_group']),
                        -r['passing_cells'],r['ingredient_key']))
                    for index,row in enumerate(rows,1):row['review_rank']=index
                    views[view]=rows
            yield dict(**{k:objective[k] for k in ['objective_id','gene_id','gene_symbol','scenario_id']},views=views)
    return stored_json(dict(schema='genosapiens.track2-lincs-recurrence/v1',paired_order=PRIORITY_POLICY,
        fallback_order=FALLBACK_ORDER,agreement_groups=PRIORITY_POLICY['groups'],scopes=Records(scopes())))


def recurrence_chunks(data):
    from genosapiens.track2_discovery_report import _text
    from genosapiens.track2_lincs_priority import score
    lines=['## LINCS recurrence',
        'Each view retains one gene, scenario, and objective. Tables include nominations from either route.',
        'Paired ordering uses the existing LINCS review policy. Fallback tables use this order:\n\n'
        '1. Agreement group.\n2. More passing cells.\n3. Ingredient key.',
        'Agreement groups use this order:\n\n'+'\n'.join(
            str(i)+'. '+GROUP_LABELS[group]+'.' for i,group in enumerate(data['agreement_groups'],1)),
        'Paired ordering uses these measures within each agreement group:\n\n'
        '1. More passing cells.\n2. Higher joint priority.\n3. More cells with positive desired correlation.\n'
        '4. Higher desired priority.\n5. Fewer cells with condition conflicts.\n6. Ingredient key.',
        'Joint priority equals passing cell count multiplied by the sum of finite local joint scores, including negative scores.',
        'Desired priority equals positive desired correlation cell count multiplied by the sum of finite local desired correlations.',
        'Paired ties share a rank. Fallback tables use ingredient keys to resolve ties.',
        'Scores remain local to each method. The tables do not calculate a pooled cross-method score.',
        'Agreement groups describe cells within each view. Objective conflicts include failed cells and condition conflicts from every method.',
        'Fallback recurrence never establishes paired research eligibility. Review ordering does not establish therapeutic rank.',
        '[Complete recurrence export](discovery/lincs-recurrence.json).']
    yield '\n\n'.join(lines)+'\n\n'
    for view,title in VIEWS.items():
        yield '### '+title+'\n\n[Complete table](discovery/'+view+'-recurrence.csv).\n\n'
        for scope in data['scopes']:
            rows=scope['views'][view]
            if not rows:continue
            yield '#### '+_text(scope['gene_symbol'])+' / '+_text(scope['scenario_id'])+' / '+_text(scope['objective_id'])+'\n\n'
            headers=['Rank','Medication','Agreement group','Pass / tested / selected','Paired research eligibility',
                'Objective failed cells','Objective condition conflicts','Local results']
            yield '| '+' | '.join(headers)+' |\n'+'|'+'---|'*len(headers)+'\n'
            for row in rows:
                local=[]
                for cell,r in row['cells'].items():
                    local.append(cell+' ('+r['method']+'): '+r['state']+', desired '+score(r['r_desired'])+
                        ', opposite '+score(r['r_opposite'])+', joint '+score(r['joint'])+
                        (', condition conflict' if r.get('conflict') else ''))
                values=[row['review_rank'],row['name'],GROUP_LABELS[row['agreement_group']],
                    f"{row['passing_cells']} / {row['tested_cell_count']} / {row['selected_cell_count']}",
                    row['paired_research_eligible'],row['objective_failed_cells'],row['objective_condition_conflicts'],
                    ' / '.join(local)]
                yield '| '+' | '.join(_text(v) for v in values)+' |\n'
            yield '\n'
        if not any(scope['views'][view] for scope in data['scopes']):yield 'No nominations meet the scope of this view.\n\n'


def recurrence_section(data):
    return ''.join(recurrence_chunks(data))
