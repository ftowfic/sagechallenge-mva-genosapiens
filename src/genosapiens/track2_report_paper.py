"""Deterministic short papers over retained scientific evidence."""
from pathlib import Path

from genosapiens.track2_report_context import anchor, context_markdown, display, project_context
from genosapiens.track2_report_model import build_model, publish_model, table_caption
from genosapiens.track2_report_ranking import score_rows, select_tables, publish_rankings
from genosapiens.track2_scientific_storage import publish_chunks, publish_json, file_identity

COLUMNS = ['Rank','Medication','Hypothesis','S','S_L','S_M','E','B','U','C','M','Coverage']
DEFINITIONS = {
    'Rank':'Competition rank for this table. Equal scores share rank.',
    'Medication':'Preferred retained name and canonical active ingredient identifier.',
    'Hypothesis':'Gene and objective representative, with a link to the opening context.',
    'S':'Combined descriptive score in points. S = S_L + S_M.',
    'S_L':'LINCS score in points. S_L = I_L + 2B + 0.5U - 0.5C.',
    'S_M':'Mechanism score in points. S_M = I_M + 2M.',
    'E':'Number of supporting evidence types. E = I_L + I_M, from zero to two.',
    'B':'Distinct passing cells with paired overexpression and knockdown results.',
    'U':'Distinct passing cells with a knockdown-only or overexpression-only result.',
    'C':'Distinct cells with failures or condition conflicts. Each cell contributes one penalty.',
    'M':'Distinct reported medication-target-action links that support this objective.',
    'Coverage':'Provisional means at least one required evidence route lacks a complete result.',
}
RULES = [
    'I_L equals one when at least one cell passes its actual method. Otherwise it equals zero.',
    'I_M equals one when M exceeds zero. Otherwise it equals zero.',
    'Reported inhibition supports decreased activity. Reported activation supports increased activity.',
    'Assumed actions and unresolved target groups earn no mechanism points.',
    'Missing cells earn no support and no penalty. Passing cells with conflicts retain support and one penalty.',
    'Each table selects at most five distinct ingredients across the report.',
    'Each ingredient uses its highest scoring objective for this metric. Different tables can select different objectives.',
    'Equal scores use name, ingredient identifier, and objective identifier for display order.',
    'Display order does not establish stronger evidence. Scores do not change scientific eligibility.',
]


def caption_markdown(label, caption):
    yield '**'+label+'. '+caption['title']+'**\n\n'
    yield caption['scope']+'\n\n'
    for column in caption['columns']:
        yield '**'+column+':** '+caption['definitions'][column]+'\n\n'
    for rule in caption['rules']:
        yield rule+'\n\n'
    yield '[Complete data and column definitions]('+caption['data']+').\n\n'


def medication_anchor(cid):
    return 'medication-'+anchor(cid).removeprefix('context-')


def table_rows(selection, context):
    labels = {h['scenario_id']:h for h in context['hypotheses']}
    for row in selection['rows']:
        hypothesis = labels[row['scenario_id']]
        if 'tissue_score' in row:
            tissue=row['tissue_score']
            values=[str(row['rank']),'['+display(row['name'])+'](#'+medication_anchor(row['canonical_active_ingredient_id'])+')',
                '['+display(row['gene_symbol'] or row['gene_id'])+' / '+hypothesis['label']+'](#'+hypothesis['anchor']+')',
                *['NA' if row['scores'][k] is None else format(row['scores'][k],'.3f') for k in ('S','S_L','S_M')],str(row['E']),
                *[str(row['components'][k]) for k in ('P','C','N','M')],display(tissue['method_labels']),
                ('Provisional' if row['provisional'] else 'Complete within scope')+' / '+str(tissue['N'])+' of '+str(len(tissue['declared_tissues']))+' tissues']
            yield '| '+' | '.join(values)+' |\n'
            continue
        values = [str(row['rank']), '['+display(row['name'])+'](#'+medication_anchor(row['canonical_active_ingredient_id'])+')',
            '['+display(row['gene_symbol'] or row['gene_id'])+' / '+hypothesis['label']+'](#'+hypothesis['anchor']+')',
            *[str(row['scores'][k]) for k in ('S','S_L','S_M')], str(row['E']),
            *[str(row['components'][k]) for k in ('B','U','C','M')],
            'Provisional' if row['provisional'] else 'Complete within scope']
        yield '| '+' | '.join(values)+' |\n'


def paired_research_summary(model):
    """Distinguish paired cell support from the saved scientific decision."""
    total = supported = eligible = 0
    if model.get('schema')=='genosapiens.track2-report-presentation/v2':
        return dict(nominations=0,paired_support=0,eligible=0)
    for row in model['scores']:
        decisions = {d['paired_research_eligible'] for d in row['scientific_decisions']
                     if 'paired_research_eligible' in d}
        if not decisions:
            continue
        from genosapiens.track2_operator_inputs import require
        require(len(decisions) == 1, 'A report objective has conflicting paired research decisions.')
        total += 1
        cells = {c['cell']['cell_id'] for c in row['cells']
                 if c['cell']['method'] == 'paired_oe_kd' and c['cell']['state'] == 'pass'}
        supported += len(cells) >= 2
        eligible += decisions.pop()
    return dict(nominations=total, paired_support=supported, eligible=eligible)


def paper_markdown(model, title, readiness, *, hypothesis_context=None, integrated_workflow=False, calibration_preparation=None):
    if model['schema'] == 'genosapiens.track2-report-presentation/v3':
        from genosapiens.track2_challenge_paper import challenge_markdown
        yield from challenge_markdown(model, title, readiness, hypothesis_context=hypothesis_context,
            integrated_workflow=integrated_workflow,calibration_preparation=calibration_preparation)
        return
    tissue=model['schema']=='genosapiens.track2-report-presentation/v2'
    yield title+'\n\n'
    yield from context_markdown(model['context'])
    yield 'Assessment status: **'+display(readiness['assessment_status'])+'**.\n\n'
    yield 'Research report ready: **'+str(readiness['research_report_ready']).lower()+'**. Submission ready: **false**.\n\n'
    yield '## Abstract\n\n'
    if model['disclosure'].get('evidence_origin')=='independently_synthetic':
        yield 'This report uses independently synthetic evidence for software validation. It does not represent a patient analysis.\n\n'
    yield 'This report reviews medication hypotheses from retained genetic and scientific evidence. It does not establish a treatment.\n\n'
    yield ('The two evidence routes assess tissue expression agreement and reported medication actions. ' if tissue else
        'The two evidence routes assess local expression agreement and reported medication actions. ')
    yield 'Each score belongs to one ingredient and one objective.\n\n'
    yield 'The main tables contain up to five ingredients each. Complete results retain other objectives, failures, conflicts, and unavailable evidence.\n\n'
    paired = paired_research_summary(model)
    if paired['nominations']:
        yield str(paired['paired_support'])+' nominations pass in at least two paired cell lines. '
        yield str(paired['eligible'])+' meet the complete paired research rule.\n\n'
        yield 'The complete rule excludes any failed cell or condition conflict across all assessed methods. '
        yield 'Paired support and paired research eligibility describe different results.\n\n'
    yield '## Biological rationale\n\n'
    yield 'The opening context identifies each retained variant interpretation and intervention objective. Alternative hypotheses remain separate.\n\n'
    yield 'An expression match does not establish functional rescue. Direction agreement alone does not establish exposure, safety, or clinical benefit.\n\n'
    yield from main_figure(model, 1)
    for number, (kind, heading) in enumerate((('combined','Combined'),('lincs','LINCS'),('mechanism','Mechanism')), 1):
        selection = model['selections'][kind]
        yield '## '+heading+' medications\n\n'
        yield from caption_markdown('Table '+str(number), model['captions'][kind])
        if selection['rows']:
            columns=model['captions'][kind]['columns']
            yield '| '+' | '.join(columns)+' |\n| '+' | '.join('---' for _ in columns)+' |\n'
            yield from table_rows(selection, model['context'])
            yield '\n'
        else:
            yield 'This route has no retained nominated medication. [Actual route outcomes](supplements/coverage.json) distinguish incomplete work from completed negative searches.\n\n'
        if len(selection['cutoff_ties']):
            yield str(len(selection['cutoff_ties']))+' additional ingredients share the cutoff score. The complete ranking retains these ties.\n\n'
        if kind in {'combined','lincs'}:
            yield from main_figure(model, 2 if kind=='combined' else 3)
    yield '## Discussion\n\n'
    yield 'These scores organize research review. The scientific records retain their original gates, eligibility decisions, conflicts, and uncertainty.\n\n'
    seen = set()
    for selection in model['selections'].values():
        for row in selection['rows']:
            cid = row['canonical_active_ingredient_id']
            if cid in seen: continue
            seen.add(cid)
            yield '<a id="'+medication_anchor(cid)+'"></a>\n\n### '+display(row['name'])+'\n\n'
            yield 'Canonical ingredient: `'+display(cid)+'`.\n\n'
            for kind, other in model['selections'].items():
                chosen = next((r for r in other['rows'] if r['canonical_active_ingredient_id'] == cid), None)
                if chosen is None: continue
                yield kind.capitalize()+' representative objective: `'+display(chosen['objective_id'])+'`.\n\n'
                comp = chosen['components']
                if tissue:
                    yield f"Tissue counts are P = {comp['P']}, C = {comp['C']}, and N = {comp['N']}.\n\n"
                    yield f"The mechanism component counts {comp['M']} applicable reported links.\n\n"
                    reasons=sorted({reason for r in chosen['tissue_score']['coverage'] for reason in r['coverage_reasons']})
                    if reasons: yield 'Coverage reasons: '+display(reasons)+'.\n\n'
                    yield 'Tissue support does not satisfy the legacy same-cell eligibility rule.\n\n'
                else:
                    yield f"Observed support includes {comp['B']} paired cells, {comp['U']} single-direction cells, and {comp['M']} applicable reported mechanism links.\n\n"
                    yield f"The penalty includes {comp['C']} failed or conflicting cells.\n\n"
                if chosen['excluded_links']:
                    yield 'Some source actions do not earn points. [Excluded contributions and original records](supplements/objective-scores.json) retain their reasons.\n\n'
                if chosen['conflicts']:
                    yield 'Retained conflicts: '+display(chosen['conflicts'])+'.\n\n'
                if chosen['provisional']:
                    yield 'This result is provisional because required route evidence remains incomplete or unavailable.\n\n'
            yield '[All objective-specific contributions and scientific decisions](supplements/objective-scores.json) remain available for this ingredient.\n\n'
    yield 'Establish the variant effect in a relevant model before interpreting a proposed intervention as rescue.\n\n'
    yield 'Test reported target action and functional rescue independently. Assess exposure and hazards before any treatment interpretation.\n\n'
    yield '## Methods\n\n'
    yield 'The report uses retained input interpretations and selected transcripts. It creates no new variant mapping or transcript selection.\n\n'
    if tissue:
        yield 'The query pools OE and KD separately through per-feature medians. Each selected cell contributes once per direction.\n\n'
        yield 'The query excludes features with missing contributor values. It applies magnitude and opposition rules after pooling.\n\n'
        yield 'Drug comparisons use the frozen feature set. A failed paired query does not switch to a single-direction method.\n\n'
        yield 'Scores use exact fractions before three-decimal display. Complete exports retain every evaluated medication and admitted condition.\n\n'
        yield 'Single-direction results remain weaker and stay outside the paired tissue denominator. The mechanism formula remains I_M + 2M.\n\n'
    else:
        yield 'The report scores the actual LINCS method in each cell. It does not convert a failed paired result into single-direction support.\n\n'
        yield 'The complete local plots retain their original selection and pagination. Main table selection does not change local plot selection.\n\n'
        yield 'Scores use integer half-point arithmetic. Repeated cells and equivalent medication-target-action links do not multiply their contribution.\n\n'
    yield 'Deterministic templates generate this report without a runtime language-model call.\n\n'
    yield 'Codex assisted implementation of the presentation code. This statement does not describe earlier research or establish human review.\n\n'
    for key in ('description','llm_disclosure','effort','limitations','validation_followup'):
        yield '**'+{'description':'Retained method','llm_disclosure':'Actual AI disclosure','effort':'Recorded effort',
                    'limitations':'Retained limitations','validation_followup':'Recorded follow-up'}[key]+' (source text):**\n\n'
        methods=model['disclosure']['methods']
        yield methods.get(key, methods.get('scope','Unavailable.') if key=='description' else 'Unavailable.')+'\n\n'
    if not model['disclosure']['methods'].get('effort'):
        yield 'This report does not estimate undocumented effort.\n\n'
    for key,label in (('input_review','Input review'),('plan_review','Plan review')):
        review=model['disclosure'].get(key)
        if review:
            yield label+' decision: '+display(review.get('decision'))+'.\n\n'
            yield 'Recorded reviewer: '+display(review.get('reviewer'))+'. Reviewer type: '+display(review.get('reviewer_type'))+'.\n\n'
        else:
            yield label+': no manual review record accompanies this input.\n\n'
    yield 'The [detailed scientific record](scientific-details.md) retains actual review attribution and source versions.\n\n'
    yield 'Report generation does not record human acceptance. Input review, plan review, and final report review remain separate.\n\n'
    yield '## References\n\n'
    yield 'The [source catalog](scientific-details.md) retains exact citations and provenance.\n\n'
    yield 'The [configuration](plan-source.json), [coverage](supplements/coverage.json), and [presentation model](presentation.json) retain reproducibility records.\n\n'
    from genosapiens.track2_report_language import report_terms, report_acknowledgement
    yield report_terms(lincs=True, diagnosis=model['input_origin']=='diagnosis_guided_vcf', exploratory=True,tissue=tissue).replace(
        '## Terms used in this report', '## Glossary\n\nTable 4 defines report terms. Term identifies the expression. Meaning gives its report-specific definition.')+'\n\n'
    if tissue:
        from genosapiens.track2_tissue_report_text import glossary
        yield '| Term | Meaning |\n|---|---|'+glossary()+'\n'
    yield report_acknowledgement()+'\n\n'
    yield '## Supplementary index\n\n'
    for item in model['supplements']:
        yield '- ['+item['title']+']('+item['path']+').\n'


def compose_paper(root, *, bundle, plan, objectives, nominations, coverage, selected_cells,
                  readiness, title, gene_metadata=None):
    """Compose a revision without reading matrices or changing scientific results."""
    model = export_scientific_supplements(root, bundle=bundle, plan=plan,
        objectives=objectives, nominations=nominations, coverage=coverage,
        selected_cells=selected_cells, gene_metadata=gene_metadata)
    publish_model(Path(root)/'presentation.json', model)
    return ''.join(paper_markdown(model, title, readiness)).encode()


def export_scientific_supplements(root, *, bundle, plan, objectives, nominations,
                                  coverage, selected_cells, gene_metadata=None):
    """Export full scientific records without composing a report presentation."""
    root = Path(root)
    folder = root/'supplements'
    folder.mkdir(mode=0o700)
    if 'exploratory_discovery' not in plan:
        publish_json(folder/'nomination-inputs.json', nominations)
    context = project_context(bundle, plan, objectives, gene_metadata)
    tissue=plan.get('exploratory_discovery',{}).get('schema')=='genosapiens.track2-exploratory-plan/v2'
    scores = score_rows(nominations, objectives, coverage, selected_cells)
    selections = select_tables(scores,tissue=tissue)
    publish_rankings(folder, scores, selections)
    publish_json(folder/'context.json', context)
    publish_json(folder/'coverage.json', coverage)
    columns,definitions,rules=COLUMNS,DEFINITIONS,RULES
    if tissue:
        from genosapiens.track2_tissue_report_text import COLUMNS as columns, DEFINITIONS as definitions, RULES as rules
    captions = {kind:table_caption(kind.capitalize()+' medication review order', columns, definitions,
        scope='Scores describe retained evidence within each representative objective.', rules=rules,
        data='supplements/'+kind+'-ranking.json') for kind in selections}
    from genosapiens.track2_report_supplements import publish_figure_index
    figures, figure_caption = publish_figure_index(root, plan)
    captions['figures'] = figure_caption
    from genosapiens.track2_report_captions import annotate_supplements
    captions.update(annotate_supplements(root))
    captions['glossary'] = table_caption('Report terms', ['Term','Meaning'],
        {'Term':'Expression defined in this report.','Meaning':'Definition within this report.'},
        scope='Terms used in the main paper.',rules=[],data='supplements/column-dictionaries.json')
    from genosapiens.track2_report_captions import export_index
    captions=export_index(root,captions)
    publish_json(folder/'column-dictionaries.json', captions)
    supplements = [dict(title='Complete '+kind+' ranking', path='supplements/'+kind+'-ranking.json') for kind in selections]
    supplements += [dict(title=title, path=path) for title, path in (
        ('All objective scores, counted links, exclusions, and conflicts','supplements/objective-scores.json'),
        ('Full genotype, transcript, quality, and alternative hypotheses','supplements/context.json'),
        ('Route coverage and source records','supplements/coverage.json'),
        ('Column dictionaries and counting rules','supplements/column-dictionaries.json'),
        ('Complete exports and shared column dictionaries','supplements/export-index.json'),
        ('Complete scientific details and local figures','scientific-details.md'),
        ('Every local figure, caption, source record, and image','supplements/figures.md'),
        ('Original plan and configuration','plan-source.json'))]
    model = build_model(input_origin=bundle.data['input_origin'], classification=plan['classification'],
        sources=dict(input=bundle.input_sha256, plan=file_identity(root/'plan-source.json')['sha256']),
        context=context, scores=scores, selections=selections, captions=captions, figures=figures,
        supplements=supplements, coverage=coverage, tissue=tissue, disclosure=dict(methods=plan['methods'],overview_state='pending',
            input_review=bundle.review, evidence_origin=bundle.data.get('evidence_origin'),
            plan_review=__import__('json').loads((root/'plan-review-source.json').read_bytes())
                if (root/'plan-review-source.json').is_file() else None))
    return model


def main_figure(model, number):
    for figure in model['figures']:
        if not figure['main'] or figure['number'] != number:
            continue
        yield figure['label']+' summarizes the corresponding retained evidence.\n\n'
        yield '!['+figure['label']+']('+figure['png']+')\n\n'
        yield figure['caption']+'\n\n'
        yield '[SVG]('+figure['svg']+') and [source data]('+figure['data']+') retain the full figure.\n\n'
        yield '[Complete local figures](supplements/figures.md) retain the separate local display scope.\n\n'
