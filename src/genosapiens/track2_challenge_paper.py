"""Findings-first prose from compact, source-bound challenge records."""
from genosapiens.track2_report_context import display
from genosapiens.track2_report_language import report_acknowledgement

METHOD_LABELS = {'paired_oe_kd':'Paired OE/KD', 'oe_only_matching':'OE-only',
                 'kd_only_reversal':'KD-only', 'legacy_cell':'Legacy cell evidence', 'mechanism':'Mechanism'}
COLUMNS = ['Medication and approval source','Gene and hypothesis','Expression evidence',
           'Reported mechanism','Conflict or evidence gap','Research interpretation']


def kd_label(row):
    expression = row.get('expression')
    if not expression or expression['method'] != 'kd_only_reversal' or expression['status'] != 'complete' or expression['strength'] <= 0:
        return None
    gene = row.get('gene_symbol') or row['gene_id']
    return 'KD-reversal hypothesis; '+str(gene)+' increase unconfirmed'


def kd_explanation(gene):
    gene = display(gene)
    return ('"Knockdown" means reducing a gene\'s activity in laboratory cells.\n\n'
        'We compared separate experiments that reduced '+gene+' activity or exposed cells to the drug.\n\n'
        'The drug\'s pattern opposed the pattern linked to reduced '+gene+' activity.\n\n'
        'This suggests a possible way to counter some effects of reduced '+gene+' activity.\n\n'
        'It does not prove that the drug increases '+gene+' or restores its function.\n\n'
        'This comparison has no OE evidence.\n\n')


def _link(model, name, label):
    return '['+label+']('+model['artifacts'][name]['path']+')'


def displayed_ties(model):
    """Count displayed ingredients separately for each method and exact rank."""
    from collections import defaultdict
    groups = {}
    displayed = defaultdict(set)
    historical = model['summary']['historical_ties']
    for row in model['selections']:
        expression = row['expression']
        if expression:
            key = row['group'], expression['method_rank']
            groups[key] = expression['tie_count']
            displayed[key].add(row['canonical_active_ingredient_id'])
        elif row['group'] == 'legacy_cell':
            key = 'legacy_cell', row['historical_rank']
            groups[key] = next((g['count'] for g in historical.get('lincs', []) if g['rank'] == key[1]), 1)
            displayed[key].add(row['canonical_active_ingredient_id'])
        mechanism = row['mechanism']
        if mechanism and mechanism['historical_rank'] is not None:
            key = 'mechanism', mechanism['historical_rank']
            groups[key] = next((g['count'] for g in historical.get('mechanism', []) if g['rank'] == key[1]), 1)
            displayed[key].add(row['canonical_active_ingredient_id'])
    return {key:dict(total=total, displayed=len(displayed[key]), omitted=total-len(displayed[key]))
            for key,total in groups.items() if total > 1}


def hypothesis_markdown(context, destination):
    """Stream retained context with shared facts and package-local references."""
    from genosapiens.track2_report_context import context_markdown
    from genosapiens.track2_scientific_storage import json_identity
    shared = {}; genes = {}; hypotheses = []
    for hypothesis in context['hypotheses']:
        unique = []; links = []
        for gene in hypothesis['genes']:
            key = json_identity(gene)['sha256']
            if key in genes:
                links.append(genes[key])
            else:
                genes[key] = hypothesis['anchor']
                unique.append(gene)
        hypotheses.append(dict(hypothesis, genes=unique))
        shared['<a id="'+hypothesis['anchor']+'"></a>'] = sorted(set(links))
    projection = dict(observations=context['observations'], hypotheses=hypotheses)
    for chunk in context_markdown(projection):
        yield chunk.replace('(supplements/context.json)', '('+destination+')').replace('](#', '](report.md#')
        if chunk.startswith('<a id="'):
            for target in shared.get(chunk.splitlines()[0], []):
                yield '- [Shared gene and location](report.md#'+target+'). This hypothesis retains its own interpretation.\n\n'


def evidence_rows(model):
    revised = model.get('narrative_version', 1) == 2
    ties = displayed_ties(model) if revised else {}
    labels={oid:h['label']+(', objective '+str(i) if len(h['objective_ids'])>1 else '')
        for h in model['summary']['hypotheses'] for i,oid in enumerate(h['objective_ids'],1)}
    for row in model['selections']:
        source = '['+display(row['name'])+']('+row['source_reference']+')'
        approval = row.get('approval_source')
        if approval:
            source += '. [Retained '+display(approval['authority'])+' snapshot]('+approval['source_reference']+').'
            source += ' Approval date, jurisdiction, and indication: unavailable in this projection.'
        else:source += '. Approval authority, date, jurisdiction, and indication: unavailable in this projection.'
        gene = display(row.get('gene_symbol') or row['gene_id'])
        hypothesis = '['+gene+' / '+display(labels.get(row['objective_id'],'Objective'))+']('+row['source_reference']+')'
        expression = row['expression']
        if expression:
            coverage = expression['coverage']
            text = METHOD_LABELS[row['group']]+'. Strength '+format(expression['strength'],'.6g')
            text += '. Percentile '+format(100*expression['numerator']/expression['denominator'],'.6g')
            text += '. Method rank '+str(expression['method_rank'])+'. Tie count '+str(expression['tie_count'])+'.'
            if expression.get('degenerate_background'):text += ' Degenerate background: all frozen strengths are equal.'
            text += ' Tissues '+str(len(coverage['contributing_tissues']))+'/'+str(len(coverage['declared_tissues']))
            text += ', cells '+str(coverage['n_cells'])+', valid conditions '+str(coverage['n_valid'])+'.'
            text += ' Undefined conditions '+str(coverage['n_undefined'])+'.'
            text += ' Contributing tissues: '+display(coverage['contributing_tissues'])+'.'
            if coverage['missing_tissues']:text += ' Missing tissues: '+display(coverage['missing_tissues'])+'.'
            if row['group']=='paired_oe_kd' and not revised:
                text += ' Historical S_L '+display(row['historical_scores']['S_L'])+'.'
                text += ' '+', '.join(k+' = '+str(row['historical_components'][k]) for k in ('P','C','N'))+'.'
        elif row['group']=='legacy_cell':
            text='Legacy cell evidence. Historical rank '+str(row['historical_rank'])+'.'
            if not revised:text+=' S_L '+display(row['historical_scores']['S_L'])+'.'
            text+=' '+', '.join(k+' = '+str(row['historical_components'][k]) for k in ('B','U','C') if k in row['historical_components'])+'.'
        else:text='No selected calibrated expression score for this objective.'
        if revised:
            historical = row['historical_scores']['S_L']
            text += ' Historical S_L = '+('NA' if historical is None else display(historical))+'.'
            components = row['historical_components']
            if all(k in components for k in ('P', 'C', 'N')):
                text += ' '+', '.join(k+' = '+str(components[k]) for k in ('P', 'C', 'N'))+'.'
            rank = expression['method_rank'] if expression else row['historical_rank']
            if row['group'] != 'mechanism' and (row['group'], rank) in ties:
                text += ' Display example.'
        mechanism=row['mechanism']
        mechanism_text = ('S_M = '+display(mechanism['component'])+'. '+str(mechanism['counted_link_count'])+
            ' counted reported links. [Sources]('+mechanism['source_reference']+').' if mechanism else
            'No retained applicable reported action supports this objective.')
        if row.get('mechanism_actions'):
            mechanism_text += ' '+'. '.join(display(r['action'])+' at '+display(r['target']) for r in row['mechanism_actions'])+'.'
        if revised:
            if mechanism is None:
                mechanism_text += ' Historical S_M = '+display(row['historical_scores']['S_M'])+'.'
            elif ('mechanism', mechanism['historical_rank']) in ties:
                mechanism_text += ' Display example.'
        gaps = str(row['conflict_count'])+' retained conflicts. '+str(row['excluded_link_count'])+' excluded mechanism links.'
        if row.get('conflicts'):gaps+=' Conflict records: '+display(row['conflicts'])+'.'
        if row.get('excluded_mechanism_reasons'):gaps+=' Excluded actions: '+display(row['excluded_mechanism_reasons'])+'.'
        for finding in row.get('blocking_findings', []):
            gaps += (' Retained blocking hazard. ' if finding['kind']=='blocking_hazard' else ' Retained action conflict. ')
            gaps += display(finding['reasons'])+'.'
        if row['provisional']:gaps+=' The source assessment remains provisional.'
        interpretation=kd_label(row) or 'Research hypothesis. Functional rescue remains unconfirmed.'
        interpretation=display(interpretation)+' [Complete evidence]('+row['source_reference']+').'
        yield '| '+' | '.join((source,hypothesis,text,mechanism_text,gaps,interpretation))+' |\n'


def challenge_markdown(model, title, readiness, *, hypothesis_context=None, integrated_workflow=False, calibration_preparation=None):
    from genosapiens.track2_operator_inputs import require
    version = model.get('narrative_version', 1)
    require(type(version) is int and version in (1, 2), 'The challenge narrative version is unsupported.')
    if version == 2:
        require(hypothesis_context is not None, 'The challenge narrative requires retained hypothesis context.')
        require([h['scenario_id'] for h in hypothesis_context['hypotheses']] ==
                [h['scenario_id'] for h in model['summary']['hypotheses']],
                'The narrative context differs from its hypothesis scope.')
    summary=model['summary']
    yield title+'\n\n## Research question and findings\n\n'
    origin={'diagnosis_guided_vcf':'diagnosis-guided variant analysis',
            'published_patient_genotype':'a published patient genotype'}[model['input_origin']]
    yield 'This report evaluates medication hypotheses from '+origin+'.\n\n'
    if summary.get('evidence_origin')=='independently_synthetic':
        yield 'This report uses independently synthetic evidence for software checks. It does not describe a patient analysis.\n\n'
    yield 'The research question asks which retained evidence supports a testable intervention hypothesis.\n\n'
    objectives=summary['objective_count'];hypotheses=summary['hypothesis_count'];unresolved=summary['unresolved_hypothesis_count']
    yield 'The report assesses '+str(objectives)+(' objective' if objectives==1 else ' objectives')
    yield ' across '+str(hypotheses)+(' hypothesis. ' if hypotheses==1 else ' hypotheses. ')
    yield str(unresolved)+(' hypothesis lacks' if unresolved==1 else ' hypotheses lack')+' an intervention objective.\n\n'
    yield 'Assessment status: **'+display(readiness['assessment_status'])+'**. '
    yield 'Research report ready: **'+str(readiness['research_report_ready']).lower()+'**. Submission ready: **false**.\n\n'
    for route,totals in summary['route_totals'].items():
        label={'lincs_cell_specific_expression':'Expression search','qualified_target_drug_linkage':'Mechanism search'}.get(route,route.replace('_',' ').capitalize())
        outcomes=', '.join(str(totals.get(key,0))+' '+label for key,label in
            (('positive','positive'),('negative','completed negative'),('unavailable','unavailable'),
             ('failed','failed'),('incomplete','incomplete'),('not_applicable','not applicable'),
             ('outside_scope','outside declared scope')) if totals.get(key,0))
        yield label+': '+(outcomes or 'no objectives in scope')+'.\n\n'
    yield 'A completed negative search found no candidate within its retained scope. Missing evidence does not establish a negative result.\n\n'
    for method in summary['derived_methods']:
        unranked=sum(method['unranked_candidates'].values())
        if unranked:
            yield METHOD_LABELS[method['method']]+': '+str(unranked)+' nominated representatives lack complete calibrated ranks. '
            yield _link(model,'derived_rankings','Unranked results')+' retain their raw scores and status.\n\n'
    for method,groups in summary['historical_ties'].items():
        if groups:
            yield 'Historical '+display(method)+' ties: '+', '.join(str(g['count'])+' ingredients at rank '+str(g['rank']) for g in groups)+'.\n\n'
    yield 'The method groups do not form one shared efficacy ranking. Exact ties remain ties despite alphabetical display order.\n\n'
    yield from figure_markdown(model,'coverage')
    if version == 2:
        yield from hypothesis_markdown(hypothesis_context, model['artifacts']['context']['path'])
    yield '## Mechanism rationale\n\n'
    yield 'Each hypothesis retains its original interpretation and intervention direction. '
    yield _link(model,'context','Complete hypothesis context')+' preserves observations, predictions, and unresolved steps.\n\n'
    genes={g['gene_id']:g for h in summary['hypotheses'] for g in h['genes']}
    if genes:
        yield 'Retained genes: '+display([g.get('symbol') or g['gene_id'] for g in genes.values()])+'.\n\n'
    for hypothesis in summary['hypotheses']:
        yield '- '+display(hypothesis['label'])+': '
        count=len(hypothesis['objective_ids'])
        yield str(count)+(' objective. ' if count==1 else ' objectives. ')+str(hypothesis['blocker_count'])+' retained mechanism blockers.\n'
    yield '\nExpression agreement does not establish a drug mechanism or functional rescue.\n\n'
    yield from figure_markdown(model,'mechanism')
    yield '## Medication hypotheses\n\n'
    yield 'Table 1 separates expression methods and reported mechanism evidence. '
    yield 'Each expression group displays at most five source-nominated representatives.\n\n'
    yield 'A retained approval record does not establish current or disease-specific approval.\n\n'
    if version == 2:
        for (method, rank), tie in displayed_ties(model).items():
            yield METHOD_LABELS[method]+' at rank '+str(rank)+': '+str(tie['total'])+' ingredients share the score. '
            yield str(tie['displayed'])+' displayed ingredients are display examples. '
            yield 'This group omits '+str(tie['omitted'])+' tied ingredients.\n\n'
    if model['selections']:
        yield '| '+' | '.join(COLUMNS)+' |\n| '+' | '.join('---' for _ in COLUMNS)+' |\n'
        yield from evidence_rows(model)
        yield '\n'
    else:
        yield 'No medication has a selected row under this presentation policy. '
        yield _link(model,'coverage','Route outcomes')+' distinguish completed searches from missing or unfinished work.\n\n'
    yield 'S_M is the historical mechanism component. It equals I_M + 2M. '
    yield _link(model,'objective_scores','Historical components and decisions')+' retain exact arithmetic, eligibility, and exclusions.\n\n'
    if version == 2:
        yield 'NA means the historical expression component is unavailable. Missing calibration does not erase a measured historical component.\n\n'
    explained=set()
    for row in model['selections']:
        gene=row.get('gene_symbol') or row['gene_id']
        if kd_label(row) and gene not in explained:
            explained.add(gene);yield kd_explanation(gene)
    yield 'OE matching means resemblance to an overexpression pattern. It does not prove that the drug increases the gene.\n\n'
    yield 'Paired matching compares separate OE and KD expression experiments. It does not demonstrate restored protein function.\n\n'
    yield 'When available, expression evidence compares drug conditions with retained genetic patterns. '
    yield ('The tissue summary' if model['derived_policy'] else 'The cell summary')+' does not describe one jointly measured experiment.\n\n'
    yield from figure_markdown(model,'evidence')
    yield '## Research value and next experiment\n\n'
    yield 'Potential Impact carries 25 percent of the challenge score. Any proposed benefit remains conditional on independent validation.\n\n'
    yield 'A relevant experiment would compare a disease model with controls before and after drug exposure. '
    yield 'Failure to improve the predefined functional outcome would reject that rescue hypothesis within the tested conditions.\n\n'
    yield '## Methods and contribution\n\n'
    yield 'Scientific Rigor carries 35 percent of the challenge score. '
    yield 'The retained source records define variant interpretation, search scope, nomination, and original scientific decisions.\n\n'
    yield 'Innovation carries 25 percent of the challenge score. This report makes no unsupported claim of methodological novelty.\n\n'
    if model['derived_policy']:
        yield 'Paired strength equals (r_oe - r_kd) / 2. OE-only strength equals r_oe. KD-only strength equals -r_kd.\n\n'
        yield 'The calculation uses a median across doses, then equal means across times, cells, and tissues. Missing directions remain missing.\n\n'
        yield 'Percentiles compare observed strength with 10,000 frozen shuffled-query strengths. Equal strengths receive half credit.\n\n'
        yield 'The percentile depends on the selected genes and shuffle rule. It is not Broad tau, statistical significance, or a probability of benefit.\n\n'
        yield 'Gene dependence limits the shuffled-query comparison. The draws do not represent independent biological experiments.\n\n'
    yield _link(model,'disclosures','Retained methods and disclosures')+' records actual source use, manual review, AI use, and effort where available.\n\n'
    methods=summary.get('methods',{})
    for key,label in (('description','Method description'),('llm_disclosure','AI provider, tier, and settings'),
                      ('effort','Recorded effort'),('limitations','Source limitations'),('validation_followup','Source validation plan')):
        value=methods.get(key,methods.get('scope') if key=='description' else None)
        if value:
            yield '**'+label+' (source quotation):**\n\n'
            yield '\n'.join('> '+line for line in str(value).splitlines())+'\n\n'
        else:yield label+': unavailable in the retained methods fields.\n\n'
    yield 'Deterministic templates generate this presentation without a runtime language-model call.\n\n'
    yield '## Reuse and execution\n\n'
    yield 'Scalability carries 15 percent of the challenge score. '
    yield 'The revision retains source evidence and local assets for independent inspection.\n\n'
    if integrated_workflow:
        yield 'The report workflow prepares missing applicable calibration automatically before rendering. Compatible retries reuse completed backgrounds without new permutations.\n\n'
        if calibration_preparation is not None:
            value=calibration_preparation
            yield 'This composition used '+str(value['entries_new'])+' new calibration entries and '+str(value['entries_reused'])+' reused entries. '
            yield 'Preparation performed '+str(value['operation_counts']['query_permutations'])+' query permutations. '
            yield 'It imported '+str(value['entries_imported'])+' entries from a retained report.\n\n'
            yield '[Preparation records](calibration-preparation.json) retain resource use and entry counts. '
            yield '[Applicability records](calibration-applicability.json) retain scope limits and exclusion reasons.\n\n'
    else:
        yield 'Ordinary reporting performs no background permutations. Missing calibration requires a separate preparation operation.\n\n'
    yield 'This presentation does not estimate undocumented runtime, cost, or effort.\n\n'
    yield '## References\n\n'+_link(model,'scientific_details','Original scientific details and references')+' retains source attribution and limitations.\n\n'
    yield '## Glossary\n\n| Term | Meaning |\n| --- | --- |\n'
    for term,meaning in (
        ('Gene expression','A measure of how active genes are in the experiment.'),
        ('OE','Overexpression, an experiment that increases gene expression.'),
        ('KD','Knockdown, an experiment that reduces gene activity.'),
        ('Paired evidence','Expression comparisons with both OE and KD experiments.'),
        ('Strength','A method-specific summary of expression correlations.'),
        ('r_oe','The Spearman correlation with the retained overexpression pattern.'),
        ('r_kd','The Spearman correlation with the retained knockdown pattern.'),
        ('Random-query percentile','The observed strength\'s position among shuffled-query strengths, with half credit for ties.'),
        ('Frozen background','A retained distribution bound to its scientific inputs and scoring policy.'),
        ('Method rank','A competition rank within one expression method. Equal score pairs share rank.'),
        ('Historical rank','The original rank under the source scoring policy.'),
        ('S_L','The historical expression component under the original tissue or cell policy.'),
        ('S_M','The historical mechanism component under the original policy.'),
        ('M','The number of distinct applicable reported medication-target-action links.'),
        ('I_M','One when M exceeds zero. Otherwise zero.'),
        ('P','The number of supporting paired tissues under the original policy.'),
        ('C','The number of conflicting paired tissues, or failed or conflicting legacy cells, under the original policy.'),
        ('N','The number of evaluable paired tissues under the original policy.'),
        ('B','The number of passing paired cells under the legacy policy.'),
        ('U','The number of passing single-direction cells under the legacy policy.'),
        ('Provisional','Required scientific work remains incomplete.')):
        yield '| '+term+' | '+meaning+' |\n'
    yield '\n## Disclosures\n\n'
    yield 'Original review records apply to their original artifacts. They do not establish acceptance of this revision.\n\n'
    for key,label in (('input_review','Input review'),('plan_review','Plan review')):
        review=summary.get('reviews',{}).get(key)
        if review:
            yield label+': '+display(review.get('decision'))+'. Recorded reviewer: '+display(review.get('reviewer'))+'. '
            yield 'Reviewer type: '+display(review.get('reviewer_type'))+'.\n\n'
        else:yield label+': no review record accompanies the source presentation.\n\n'
    yield _link(model,'disclosures','Source disclosures')+' identifies recorded review and methods. Missing provider settings, effort, and costs remain unavailable.\n\n'
    yield report_acknowledgement()+'\n\n'
    yield '## Supplementary index\n\n'
    for supplement in model['supplements']:
        yield '- ['+display(supplement['title'])+']('+supplement['artifact']['path']+').\n'


def figure_markdown(model, kind):
    for figure in model['figures']:
        if figure['provenance']['kind']!=kind:continue
        yield '!['+figure['label']+']('+figure['png']['path']+')\n\n'
        yield figure['caption']+'\n\n'
        yield '['+figure['label']+' vector image]('+figure['svg']['path']+'). '
        yield '['+figure['label']+' source data]('+figure['data']['path']+').\n\n'
