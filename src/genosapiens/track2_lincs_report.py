"""Deterministic linked scientific figures and their mandatory report bindings."""
import csv
import io
import json
from pathlib import Path
import subprocess
import shutil

from genosapiens.track2_lincs import ROOT
from genosapiens.track2_lincs_scoring import POLICY
from genosapiens.track2_lincs_priority import (CAVEAT, PRIORITY_POLICY, ELIGIBILITY_POLICY, priority_table, priority_explanation)
from genosapiens.track2_lincs_integration import research_priority
from genosapiens.track2_published_genotype import canonical, sha256


def csv_bytes(rows, fields):
    stream=io.StringIO(newline='')
    writer=csv.DictWriter(stream,fieldnames=fields,extrasaction='ignore');writer.writeheader()
    writer.writerows(rows)
    return stream.getvalue().encode()


def native_render_binding(discovery):
    from genosapiens.track2_discovery_resources import file_digest
    from genosapiens.track2_construction_dependencies import python_dependencies
    from genosapiens.track2_challenge_overviews import pinned_runtime
    return dict(job_contract='genosapiens.track2-figure-job/v1',runtime=pinned_runtime(),
        discovery_sha256=file_digest(Path(discovery)/'discovery-results.json'),renderer=file_digest(ROOT/'scripts/render_track2_lincs.R'),
        tissue_renderer=file_digest(ROOT/'scripts/render_track2_overviews.R'),
        code=python_dependencies(ROOT,[('genosapiens.track2_figure_stream','figure_pages'),
            ('genosapiens.track2_tissue_figures','tissue_condition_pages'),('genosapiens.track2_tissue_figures','tissue_panels'),
            ('genosapiens.track2_figure_workers','render_jobs')]))


def render_native_jobs(discovery,output, *, figure_workers=2):
    """Render retained vectors in the LINCS image. No matrix access is needed."""
    from genosapiens.track2_discovery import discovery_lock
    from genosapiens.track2_discovery_outcome import record_failure
    from genosapiens.track2_operator_inputs import OperatorError
    with discovery_lock(output):
        try:
            return _render_native_jobs(discovery,output,figure_workers=figure_workers)
        except BaseException as error:
            record_failure(output,stage='figures',role='figure-artifacts',error=error)
            if isinstance(error,(KeyboardInterrupt,SystemExit)):raise
            raise OperatorError('Native figure publication failed. Read the current diagnostics.') from None


def _render_native_jobs(discovery,output, *, figure_workers=2):
    from genosapiens.track2_checked_pages import CheckedPageStreams
    context=CheckedPageStreams()
    try:return _render_native_jobs_checked(discovery,output,context,figure_workers=figure_workers)
    finally:context.close()


def _render_native_jobs_checked(discovery,output,page_context, *, figure_workers=2):
    from genosapiens.track2_discovery import _read,validate_job,load_completed_discovery
    from genosapiens.track2_run_state import RunState
    from genosapiens.track2_scientific_storage import publish_chunks,json_chunks,json_identity,file_identity
    from genosapiens.track2_figure_stream import check_figure_pages
    from genosapiens.track2_figure_workers import render_jobs,effective_workers
    from genosapiens.track2_challenge_overviews import pinned_runtime
    effective_workers(figure_workers)
    discovery=Path(discovery);output=Path(output);state=RunState(output,direct_stage='figures')
    result=_read(discovery/'discovery-results.json')
    load_completed_discovery(discovery,result['bindings'],allow_partial=True)
    binding=native_render_binding(discovery)
    completed,reason=state.completed('figures',binding)
    if completed:return completed
    state.restart('figures',reason)
    attempts=output.parent/(output.name+'-attempts')
    def jobs():
        for job in result['jobs']:
            if not job['completion_sha256']:continue
            source=discovery/'jobs'/job['job_id'];local=validate_job(source);query=_read(source/'query.json')
            target=output/local['gene_id']/job['job_id']
            def prepare(folder,query=query,local=local,source=source):
                page_context.publish(folder,query,local,source/'conditions.jsonl')
                return page_context._key(query,local,source/'conditions.jsonl')
            def check(folder,query=query,local=local,source=source):
                check_figure_pages(folder,query,local,source/'conditions.jsonl',context=page_context)
            script='render_track2_overviews.R' if query['schema']=='genosapiens.track2-tissue-query/v1' else 'render_track2_lincs.R'
            yield dict(target=target,renderer=ROOT/'scripts'/script,runtime=pinned_runtime(),prepare=prepare,check=check)
        if 'evaluated_tissue_scores' in result:
            from genosapiens.track2_tissue_figures import tissue_panels
            expected={}
            def prepare(folder):
                publish_chunks(folder/'figure-data.jsonl',(chunk for figure in tissue_panels(result['evaluated_tissue_scores']) for chunk in json_chunks(figure)))
                expected.update(file_identity(folder/'figure-data.jsonl'))
                return binding
            def check(folder):
                require(file_identity(folder/'figure-data.jsonl')==expected,'Tissue overview page data changed.')
                with (folder/'figure-data.jsonl').open() as stream:
                    for line in stream:
                        figure=json.loads(line)
                        for suffix in ('.png','.svg'):
                            require(file_identity(folder/(figure['id']+suffix))['size']>0,'A tissue overview is missing.')
            from genosapiens.track2_operator_inputs import require
            yield dict(target=output/'tissues',renderer=ROOT/'scripts/render_track2_overviews.R',runtime=pinned_runtime(),prepare=prepare,check=check)
    outcome=render_jobs(jobs(),attempts,workers=figure_workers)
    targets={Path(j['target']).absolute() for j in outcome['jobs']}
    # Preserve stale targets outside the new complete artifact inventory.
    import uuid
    for path in list(output.rglob('figure-data.jsonl')):
        if path.parent.absolute() not in targets:
            path.parent.rename(attempts/('obsolete-'+uuid.uuid4().hex))
    count=0
    def index():
        nonlocal count
        for target in sorted(targets):
            with (target/'figure-data.jsonl').open() as stream:
                for line in stream:
                    figure=json.loads(line)
                    yield from json_chunks(dict(path=str(target.relative_to(output.absolute())/figure['id']),caption=figure['caption'],alt=figure['alt']))
                    count+=1
    publish_chunks(output/'index.jsonl',index())
    return state.publish('figures',binding=binding,result=dict(figures=count))


def figure_records(result):
    figures=[]
    for cell_id,cell in sorted(result['cells'].items()):
        q=result['queries'][cell_id]
        genes=[g for g in q['genes'] if g['selected']]
        displayed=[r for r in cell['medications'] if r['state']=='pass'][:POLICY['top_n_display']]
        ids={r['medication_entity_key'] for r in displayed}
        conditions=[c for c in cell['conditions'] if c['medication_entity_key'] in ids]
        # Stable display order; all admitted conditions behind displayed medications.
        order={r['medication_entity_key']:i for i,r in enumerate(displayed)}
        conditions.sort(key=lambda c:(order[c['medication_entity_key']],float(c['pert_time']),c['sig_id']))
        measured=[g[k] for g in genes for k in ['oe','kd'] if g[k] is not None]
        measured += [v for c in conditions for v in c['values'] if v is not None]
        limit=max([1.]+[abs(v) for v in measured])
        common={'cell_id':cell_id,'query_sha256':q['query_sha256'],
                'gene_ids':[g['gene_id'] for g in genes],'top_n_display':POLICY['top_n_display'],
                'palette':['#2166AC','#F7F7F7','#D97706'],'missing':'gray labeled NA',
                'expression_limit':limit,'clipping':False}
        timing=', '.join(f"{m['pert_type']} {m['pert_time']} h" for m in q['genetic_context'])
        caption=(f"{cell_id}: {timing}. Drug profiles are 6/24 h. KD means reduced target expression. "
            'OE means extra target expression.\n\nThe query deliberately selects opposite responses. This selection does not provide independent validation. '
            'An expression match is an unproven research hypothesis, not functional or disease rescue.')
        figures.append({**common,'id':cell_id+'-selection','kind':'selection','query':q,
            'caption':caption,'alt':f"{cell_id} selection: {len(q['genes'])} measured genes, {len(genes)} retained. Status: {q['state']}."})
        if genes:
            for gi in range(0,len(genes),24):
                for ci in range(0,max(1,len(conditions)),8):
                    subset=conditions[ci:ci+8]
                    figures.append({**common,'id':f'{cell_id}-heatmap-{gi//24+1}-{ci//8+1}',
                        'kind':'heatmap','genes':genes[gi:gi+24], 'gene_offset':gi,
                        'conditions':subset, 'condition_ids':[c['sig_id'] for c in subset],
                        'caption':caption+'\n\nBlue/orange show lower/higher standardized MODZ relative to controls, '
                            'not fold change or benefit. Desired arrows are categorical. The figure does not clip values. '
                            f"Top {POLICY['top_n_display']} local matches at most. The pages include all their admitted conditions.",
                        'alt':f"{cell_id} gene-by-condition map, gene page {gi//24+1}, condition page {ci//8+1}. Local KD/OE and desired arrows."})
        figures.append({**common,'id':cell_id+'-ranking','kind':'ranking','medications':displayed,
            'n_tested':len(cell['medications']), 'n_pass':sum(r['state']=='pass' for r in cell['medications']),
            'n_fail':sum(r['state']=='fail' for r in cell['medications']),
            'n_unassessed':sum(r['state']=='unassessed' for r in cell['medications']),
            'caption':caption+'\n\nSpearman measures similarity of response ordering. Both checks must agree on one condition. '
                'The joint score is the weaker check per condition. Each component and joint score uses the median across all conditions at each time. '
                'The next step takes the mean of the time medians with equal weights.\n\n'
                'The label identifies the negated opposite correlation. Conflicting conditions remain in the downloads. '
                'The scores do not imply statistical significance or a probability of efficacy.',
            'alt':f"{cell_id} local expression ranks and signed component/joint scores. The figure shows {len(displayed)} matches."})
    return figures


def render_lincs(result, root, *, retained_figures_root=None):
    root=Path(root);folder=root/'lincs';folder.mkdir(exist_ok=False)
    figures=figure_records(result)
    (folder/'results.json').write_bytes(canonical(result))
    (folder/'figure-data.json').write_bytes(canonical(figures))
    values=[]; selection=[]; rankings=[]; cross=[]
    for cell_id,cell in sorted(result['cells'].items()):
        query=result['queries'][cell_id]
        for g in query['genes']: selection.append({'cell_id':cell_id,'query_sha256':query['query_sha256'],**g})
        for c in cell['conditions']:
            for i,gene in enumerate(c['gene_ids']):
                values.append({'cell_id':cell_id,'query_sha256':query['query_sha256'],'gene_id':gene,
                    'sig_id':c['sig_id'],'medication_entity_key':c['medication_entity_key'],
                    'drug':c['pert_iname'],'dose':c['pert_dose'],'dose_unit':c['pert_dose_unit'],
                    'time_hours':c['pert_time'],'modz':c['values'][i], 'state':c['state'],
                    **{k:c[k] for k in ['r_desired','r_opposite','joint']}})
        rankings.extend(cell['medications'])
    for row in result['cross_cell']:
        for c,local in row['cells'].items():cross.append({'medication_entity_key':row['medication_entity_key'],**local})
    fields=['cell_id','query_sha256','medication_entity_key','name','gene_count','rank','r_desired','r_opposite','joint',
            'state','n_conditions','n_agree','agreement_fraction','n_undefined','conflict','times_hours']
    exports={'gene-selection.csv':(selection,['cell_id','query_sha256','gene_id','gene_symbol','feature_class','oe','kd','desired','opposite','desired_direction','selected','reason']),
             'gene-condition-values.csv':(values,['cell_id','query_sha256','gene_id','sig_id','medication_entity_key','drug','dose','dose_unit','time_hours','modz','state','r_desired','r_opposite','joint']),
             'medication-rankings.csv':(rankings,fields), 'cross-cell.csv':(cross,fields)}
    priority = research_priority(result)
    priority_fields = ['review_rank','medication_entity_key','name','agreement_group','passing_cells',
        'tested_cell_count','selected_cell_count','finite_score_cells','desired_direction_cells',
        'opposed_direction_cells','joint_score_sum','joint_priority','desired_priority','failed_cells',
        'conflicting_cells','missing_cells','research_eligible','eligibility_bases']
    exports['review-priority.csv'] = (priority, priority_fields)
    (folder/'review-priority.json').write_bytes(canonical({'policy':PRIORITY_POLICY,
        'research_eligibility_policy':ELIGIBILITY_POLICY,'rows':priority}))
    (folder/'review-priority.md').write_text('# LINCS medication review order\n\n' + CAVEAT + '\n\n'
        + priority_explanation() + '\n\n' + priority_table(priority, sorted(result['cells'])) + '\n')
    for name,(rows,columns) in exports.items():(folder/name).write_bytes(csv_bytes(rows,columns))
    retained_graphics=None
    if retained_figures_root is None:
        subprocess.run(['Rscript',str(ROOT/'scripts/render_track2_lincs.R'),str(folder/'figure-data.json'),str(folder)],
                       check=True,capture_output=True,text=True)
    else:
        source=Path(retained_figures_root)
        original=json.loads((source/'lincs/results.json').read_bytes())
        validate_figures(source,original)
        if figure_records(original)!=figures:
            raise ValueError('Frozen graphics do not represent the current scientific figure records')
        for figure in figures:
            for suffix in ('.png','.svg'):
                shutil.copyfile(source/'lincs'/(figure['id']+suffix),folder/(figure['id']+suffix))
        retained_graphics={'source_manifest_sha256':sha256((source/'lincs/manifest.json').read_bytes()),
                           'source_result_sha256':sha256(canonical(original))}
    names=['lincs/'+p.name for p in sorted(folder.iterdir()) if p.is_file()]
    manifest={'schema':'genosapiens.track2-lincs-figures/v1','result_sha256':sha256(canonical(result)),
        'figure_data_sha256':sha256(canonical(figures)), 'figures':[{k:f[k] for k in
            ['id','kind','cell_id','query_sha256','gene_ids','top_n_display','palette','missing','expression_limit','clipping','caption','alt']} for f in figures],
        'outputs':{name:sha256((root/name).read_bytes()) for name in names}}
    if retained_graphics is not None:
        manifest['retained_graphics']=retained_graphics
    (folder/'manifest.json').write_bytes(canonical(manifest));names.append('lincs/manifest.json')
    return lincs_section(result, figures), names


def lincs_section(result, figures):
    """Render retained LINCS results without executing a scientific query."""
    lines=['','## Cell-specific LINCS expression nominations','',
        CAVEAT,'',
        'Earlier cross-cell transfer experiments and reagent-sensitivity findings remain limitations. '
        'These independently derived local queries are a different experiment. Small queries have coarse '
        'correlation resolution. The gene count is not a replicate count.','',
        '[All structured results and review decisions](lincs/results.json) · [All selected and excluded genes](lincs/gene-selection.csv) · '
        '[All scored condition values](lincs/gene-condition-values.csv) · [Local rankings](lincs/medication-rankings.csv) · '
        '[Cross-cell coverage](lincs/cross-cell.csv)','',
        '| Cell | Query genes | Tested medications | Local nominations | Failed | Unassessed |',
        '|---|---:|---:|---:|---:|---:|']
    for c,cell in sorted(result['cells'].items()):
        ms=cell['medications'];lines.append(f"| {c} | {cell['gene_count']} | {len(ms)} | {sum(m['state']=='pass' for m in ms)} | {sum(m['state']=='fail' for m in ms)} | {sum(m['state']=='unassessed' for m in ms)} |")
    lines += ['', '### How to read the medication order', '', priority_explanation(), '',
              '[Complete ordered medication table](lincs/review-priority.md) · [Review-priority values and policy](lincs/review-priority.json)']
    for f in figures:
        cell = f['cell_id']
        timing = ', '.join(('extra target expression' if m['pert_type']=='trt_oe' else 'reduced target expression')
                          + f" at {m['pert_time']} h" for m in result['queries'][cell]['genetic_context'])
        captions = {
            'selection': f"{cell}: {timing}. The query deliberately selects opposite responses. This selection does not provide independent validation.",
            'heatmap': f"{cell}: blue shows lower expression. Orange shows higher expression. Arrows show the desired direction. "
                       'The genetic and drug columns share a zero-centred scale. The pages include all admitted conditions for the top five local matches.',
            'ranking': f"{cell}: each drug must match the desired pattern and oppose the opposite pattern. "
                       'The weaker check determines its joint score. These local ranks differ from the review order across cells. They do not establish functional rescue.',
        }
        lines += ['',f"![{f['alt']}](lincs/{f['id']}.png)",'',captions[f['kind']],
                  f"[Scalable figure](lincs/{f['id']}.svg) · [Full caption and figure data](lincs/manifest.json)"]
    lines += ['','### Scientific review queue','',
        f"Local queries nominated {len(result.get('candidate_decisions',[]))} resolved medications. "
        f"{sum(r['research_eligible'] for r in result.get('candidate_decisions',[]))} qualify for research under the recorded eligibility policy. "
        'The report identifies gaps in independent action, mechanism, current warnings, exposure, hazards and validation plans. '
        'Synthesis does not issue live queries for each medication.','']
    return '\n'.join(lines)


def validate_figures(root, result):
    root=Path(root);folder=root/'lincs'
    manifest=json.loads((folder/'manifest.json').read_bytes())
    expected=figure_records(result)
    if (manifest['result_sha256']!=sha256(canonical(result))
            or manifest['figure_data_sha256']!=sha256(canonical(expected))
            or (folder/'figure-data.json').read_bytes()!=canonical(expected)
            or (folder/'results.json').read_bytes()!=canonical(result)):
        raise ValueError('LINCS figure records are stale or inconsistent')
    required={'lincs/'+f['id']+suffix for f in expected for suffix in ['.png','.svg']}
    required |= {'lincs/'+name for name in ['results.json','figure-data.json','gene-selection.csv',
                 'gene-condition-values.csv','medication-rankings.csv','cross-cell.csv',
                 'review-priority.csv','review-priority.json','review-priority.md']}
    if set(manifest['outputs'])!=required:
        raise ValueError('required LINCS figure assets missing')
    expected_priority = {'policy':PRIORITY_POLICY,'research_eligibility_policy':ELIGIBILITY_POLICY,
                         'rows':research_priority(result)}
    if (folder/'review-priority.json').read_bytes() != canonical(expected_priority):
        raise ValueError('LINCS review priority is stale or inconsistent')
    for name,digest in manifest['outputs'].items():
        if not (root/name).is_file() or sha256((root/name).read_bytes())!=digest:
            raise ValueError('LINCS figure asset is missing or changed')
