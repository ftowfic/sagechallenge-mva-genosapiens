"""Overview requests depend on presentation selections, never local plot ranks."""
from pathlib import Path
import subprocess

from genosapiens.track2_report_context import display
from genosapiens.track2_scientific_storage import publish_json, file_identity, json_identity

SCRIPT = Path(__file__).resolve().parents[2]/'scripts/render_track2_overviews.R'


def scope_caption(rows):
    sentences=[]
    for row in rows:
        sentences.extend(['Medication: '+display(row['name'])+'.',
            'Gene: '+display(row.get('gene_symbol') or row['gene_id'])+'.',
            'Hypothesis: '+display(row['scenario_id'])+'.', 'Objective: '+display(row['objective_id'])+'.'])
    return '\n\n'.join(sentences)


def overview_records(model):
    selected = {r['scenario_id'] for s in model['selections'].values() for r in s['rows']}
    observations = {o['observation_id']:o for o in model['context']['observations']}
    context_rows = []
    for hypothesis in model['context']['hypotheses']:
        if selected and hypothesis['scenario_id'] not in selected: continue
        variants = [observations[key] for key in hypothesis['observation_ids']]
        if not variants: continue
        context_rows.append(dict(hypothesis=hypothesis['label'], scenario_id=hypothesis['scenario_id'],
            genes=hypothesis['genes'],
            variants=[dict(label=v['label'], identity=v['observation_id'],
                description=v.get('reported_coding') or v['observation_id']) for v in variants],
            effect=hypothesis['mechanism_inputs'].get('direction', 'unknown'),
            evidence_status=hypothesis['mechanism_inputs'].get('joint_evidence', 'unknown'),
            process='Unresolved in this projection',
            objectives=[dict(objective_id=o['objective_id'], desired_modulation=o.get('desired_modulation','unknown'))
                        for o in hypothesis['objectives']]))
    records = []
    if context_rows:
        records.append(dict(id='overview-context',number=1,kind='context',rows=context_rows,
            caption='Figure 1 connects retained variants, interpreted effects, and intervention objectives. '
                'Boxes identify source status and unresolved process mapping. Arrows indicate the retained hypothesis, not proven rescue. '
                'The data identify each gene, hypothesis, allele, and objective. No quantitative axis applies.'))
    combined = model['selections']['combined']['rows']
    if combined:
        records.append(dict(id='overview-contributions',number=2,kind='contributions',rows=combined,
            caption='Figure 2 displays the combined-table representatives. The horizontal axis gives descriptive score points. '
                'Positive segments show I_L, 2B, 0.5U, I_M, and 2M. The negative segment shows -0.5C. '
                'I_L indicates any passing cell. B counts passing paired cells. U counts passing single-direction cells. '
                'C counts distinct failed or conflicting cells. I_M indicates any supporting reported action link. '
                'M counts applicable reported action links. Assumed actions earn no points. '
                'Missing evidence contributes no support or penalty. '
                'Labels identify every component and total. Each row retains one gene, hypothesis, and objective. '
                'Scores do not establish clinical benefit. Source data retain conflicts and excluded links.\n\n'+scope_caption(combined)))
    lincs = model['selections']['lincs']['rows']
    if lincs:
        cells = {}
        for row in lincs:
            for contribution in row['cells']:
                cell = contribution['cell']
                key=(row['canonical_active_ingredient_id'],row['objective_id'],cell['cell_id'],cell['method'])
                item=cells.setdefault(key,dict(ingredient=row['canonical_active_ingredient_id'],name=row['name'],
                    gene_id=row['gene_id'],scenario_id=row['scenario_id'],objective_id=row['objective_id'],
                    cell_id=cell['cell_id'],method=cell['method'],states=[],conflict=False,
                    query_gene_count=cell.get('query_gene_count',cell.get('gene_count')), source_references=[]))
                if cell['state'] not in item['states']:item['states'].append(cell['state'])
                item['conflict'] |= bool(cell.get('conflict'))
                item['source_references'].append(contribution['source_reference'])
        for item in cells.values():
            item['states']=sorted(item['states'])
            item['state']=item['states'][0] if len(item['states'])==1 else 'mixed'
            item['symbols']=[('P' if item['method']=='paired_oe_kd' else 'U') if state=='pass'
                else 'F' if state=='fail' else 'NA' for state in item['states']]
        records.append(dict(id='overview-cells',number=3,kind='cells',rows=list(cells.values()),
            caption='Figure 3 displays retained cells for the LINCS-table representatives. Rows identify medication, cell, actual method, and available query size. '
                'P means paired support. U means single-direction support. F means failure. NA means unavailable evidence. '
                'An exclamation mark identifies condition conflict, including conflict in a passing cell. '
                'The label n gives the selected query-gene count. Rows summarize all admitted conditions for each medication and cell. '
                'Color repeats these text labels. Query sizes do not imply a shared gene set. '
                'Source data retain gene, hypothesis, objective, and contribution references.\n\n'+scope_caption(lincs)))
    if model.get('schema')=='genosapiens.track2-report-presentation/v2':
        from genosapiens.track2_tissue_figures import tissue_panels, TISSUE_CAPTION
        records=[r for r in records if r['number']==1]
        if combined:
            records.append(dict(schema='genosapiens.track2-tissue-figure/v1',id='overview-contributions',number=2,
                kind='tissue_contributions',rows=combined,caption='Figure 2 displays S_L and S_M for the combined-table representatives. '
                    'S_L uses the normalized tissue fraction. S_M equals I_M + 2M. S equals their sum. '
                    'Unavailable S_L contributes zero to S and still displays NA. Signed labels retain negative contributions. '
                    'Scores do not establish treatment benefit.'))
        if lincs:
            pages=list(tissue_panels((dict(r['tissue_score'],name=r['name']) for r in lincs),main=True,page_size=5))
            records.append(dict(schema='genosapiens.track2-tissue-figure/v1',id='overview-tissues',number=3,kind='tissues',
                rows=[r for page in pages for r in page['rows']],panels=[p for page in pages for p in page['panels']],
                caption=TISSUE_CAPTION))
    if context_rows:
        records[0]['caption']+='\n\n'+'\n\n'.join('Hypothesis: '+r['hypothesis']+'. Gene: '+
            display([g['symbol'] or g['gene_id'] for g in r['genes']])+'. Evidence status: '+display(r['evidence_status'])+'.'
            for r in context_rows)
    return records


def overview_binding(model):
    records = overview_records(model)
    return dict(policy_sha256=model['policy_sha256'], data_sha256=json_identity(records)['sha256'],
                renderer_sha256=file_identity(SCRIPT)['sha256'])


def render_overviews(root, model):
    """Run only in a renderer runtime with R and jsonlite available."""
    root = Path(root)
    root.mkdir(mode=0o700,parents=True,exist_ok=True)
    records = overview_records(model)
    if (root/'manifest.json').is_file():
        from genosapiens.track2_scientific_reader import read_scientific
        previous = read_scientific(root/'manifest.json')
        if previous['binding'] == overview_binding(model) and all(
                (root/name).is_file() and file_identity(root/name)['sha256'] == digest
                for name,digest in previous['outputs'].items()):
            return previous
    publish_json(root/'overview-data.json', records)
    if records:
        subprocess.run(['Rscript',str(SCRIPT),str(root/'overview-data.json'),str(root)],
                       check=True,capture_output=True,text=True)
    manifest = dict(schema='genosapiens.track2-report-overviews/v1', binding=overview_binding(model),
        figures=records, outputs={name:file_identity(root/name)['sha256'] for name in
            ['overview-data.json', *(r['id']+suffix for r in records for suffix in ('.png','.svg'))]})
    for record in records:
        for suffix in ('.png','.svg'):
            if record['id']+suffix not in manifest['outputs']:
                raise ValueError('An overview image is missing.')
    publish_json(root/'manifest.json', manifest)
    return manifest
