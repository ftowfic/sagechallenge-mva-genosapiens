"""Complete figure navigation, without changing local scientific plot selection."""
from pathlib import Path

from genosapiens.track2_report_context import display
from genosapiens.track2_report_model import table_caption
from genosapiens.track2_scientific_reader import read_scientific
from genosapiens.track2_scientific_records import stored_json
from genosapiens.track2_scientific_storage import Records, publish_json, publish_chunks, file_identity

FIGURE_COLUMNS = ['Figure','Gene','Hypothesis','Cell','Method','Scope','Type','Page','Assets','Data','Provenance']
FIGURE_DEFINITIONS = dict(zip(FIGURE_COLUMNS, [
    'Stable supplementary figure label.', 'Retained target gene identifier.',
    'All objectives that use this local experiment.', 'Cell line used in this experiment.',
    'Actual genetic comparison method. Missing directions remain missing.',
    'Original local medication display selection. Main table selection does not change this scope.',
    'Selection, heatmap, or local medication ranking.', 'Original figure identifier, including pagination where applicable.',
    'Complete PNG and SVG figure files.', 'Retained figure data, including the original caption.',
    'Source query hash and scientific output directory.'
]))


def local_figures(root, plan):
    root = Path(root)
    if 'exploratory_discovery' in plan:
        bindings = {}
        for row in plan['exploratory_discovery']['scenario_job_bindings']:
            bindings.setdefault(row['job_id'], []).append(row['objective_id'])
        from genosapiens.track2_figure_stream import condition_rows
        result = read_scientific(root/'discovery/results.json')
        for job in result['jobs']:
            if not job['completion_sha256']: continue
            path = next((root/'lincs').glob('*/'+job['job_id']+'/figures/figure-data.jsonl'))
            digest = file_identity(path)['sha256']
            for index, record in enumerate(condition_rows(path)):
                yield _figure(root, path, record, job['job_id']+'-'+record['id'], bindings[job['job_id']], index, digest)
        tissue_path=root/'lincs/tissues/figure-data.jsonl'
        if tissue_path.is_file():
            digest=file_identity(tissue_path)['sha256']
            for index,record in enumerate(condition_rows(tissue_path)):
                yield _figure(root,tissue_path,record,record['id'],[record['objective_id']],index,digest)
    elif (root/'lincs/figure-data.json').is_file():
        path = root/'lincs/figure-data.json'
        result = read_scientific(root/'lincs/results.json')
        digest = file_identity(path)['sha256']
        for index, record in enumerate(read_scientific(path)):
            query = result['queries'][record['cell_id']]
            objectives = sorted({new for b in result.get('case_objective_bindings', [])
                for old,new in b['objectives'].items() if old == query['objective_id']}) or [query['objective_id']]
            yield _figure(root, path, dict(record, gene_id=query['target_id'], method=query.get('method','paired_oe_kd')),
                record['id'], objectives, index, digest)


def _figure(root, path, record, identity, objectives, index, digest):
    prefix = str(path.parent.relative_to(root))+'/'+record['id']
    if record.get('schema')=='genosapiens.track2-tissue-figure/v1':
        return dict(label='Figure S-'+identity,main=False,gene_id=record['gene_id'],objectives=objectives,
            tissue_ids=record.get('tissue_ids',[record.get('tissue_id')]),cell_id=record.get('cell_id'),
            method=record.get('method'),kind=record['kind'],page=record['id'],
            scope='All evaluated medications.' if record['kind']=='tissues' else 'All admitted conditions for this medication and drug cell.',
            caption=record['caption'],alt=record['alt'],png=prefix+'.png',svg=prefix+'.svg',
            data=str(path.relative_to(root))+'#/'+str(index),
            provenance=dict(query_sha256=record.get('query_sha256'),data_sha256=digest,
                directory=str(path.parent.relative_to(root)),source_links=record['source_links']),
            units='Signed tissue contributions.' if record['kind']=='tissues' else 'Spearman correlations from -1 to +1.',
            missing=record['missing'],clipping=False)
    return dict(label='Figure S-'+identity, main=False, gene_id=record['gene_id'], objectives=objectives,
        cell_id=record['cell_id'], method=record['method'], kind=record['kind'], page=record['id'],
        scope='Original local top '+str(record['top_n_display'])+' display and all its admitted conditions.',
        caption=record['caption'], alt=record['alt'], png=prefix+'.png', svg=prefix+'.svg',
        data=str(path.relative_to(root))+'#/'+str(index), provenance=dict(query_sha256=record['query_sha256'],
            data_sha256=digest, directory=str(path.parent.relative_to(root))),
        units='Retained source units. Standardized expression values are not fold changes.',
        missing=record['missing'], clipping=record['clipping'])


def publish_figure_index(root, plan):
    root = Path(root)
    figures = stored_json(Records(local_figures(root, plan)))
    tissue=plan.get('exploratory_discovery',{}).get('schema')=='genosapiens.track2-exploratory-plan/v2'
    columns=['Cell or tissue' if c=='Cell' else c for c in FIGURE_COLUMNS]
    definitions=dict(FIGURE_DEFINITIONS)
    if tissue:
        definitions['Cell or tissue']='Drug cell for condition pages, or all declared tissue groups for tissue panels.'
        del definitions['Cell']
        definitions['Scope']='Complete evaluated medication or admitted condition inventory for this page.'
        definitions['Type']='Tissue evidence panel or individual condition correlations.'
    else: columns=FIGURE_COLUMNS
    caption = table_caption('Complete local figure index', columns, definitions,
        scope='Every generated local figure and page appears here.',
        rules=['All evaluated medications appear in the tissue pages.' if tissue else 'Original local selection and pagination remain unchanged.',
               'Each figure retains its complete caption and source data.'], data='supplements/figures.json')
    publish_json(root/'supplements/figures.json', figures)
    def gallery():
        from genosapiens.track2_report_paper import caption_markdown
        yield '# Complete local figures\n\n'
        yield from caption_markdown('Table S1', dict(caption, data='figures.json'))
        yield '| '+' | '.join(columns)+' |\n| '+' | '.join('---' for _ in columns)+' |\n'
        for figure in figures:
            location=figure.get('cell_id') or ', '.join(figure.get('tissue_ids',[]))
            values = [figure['label'],figure['gene_id'],', '.join(figure['objectives']),location,figure['method'],
                figure['scope'],figure['kind'],figure['page'],
                '[PNG](../'+figure['png']+') / [SVG](../'+figure['svg']+')',
                '[Data](../'+figure['data']+')',figure['provenance'].get('query_sha256') or figure['provenance']['data_sha256']]
            yield '| '+' | '.join(display(v) for v in values)+' |\n'
        yield '\n'
        for figure in figures:
            yield '## '+figure['label']+'\n\n'
            yield '!['+display(figure['alt'])+'](../'+figure['png']+')\n\n'
            yield figure['caption']+'\n\n'
            if 'tissue_ids' in figure:
                yield 'Gene: '+figure['gene_id']+'. Tissue: '+', '.join(figure['tissue_ids'])+'.\n\n'
                for link in figure['provenance']['source_links']:
                    if link.startswith('../science/'):
                        link=str(Path(figure['provenance']['directory']).parent/'science'/Path(link).name)
                    yield '[Source record](../'+link+').\n\n'
            else:
                yield 'Gene: '+figure['gene_id']+'. Cell: '+figure['cell_id']+'. Method: '+figure['method']+'.\n\n'
            yield 'Objectives: '+', '.join(figure['objectives'])+'.\n\n'
            yield figure['units']+'\n\n'
            yield 'Missing values: '+figure['missing']+'.\n\n'
            yield '[SVG](../'+figure['svg']+') and [source data](../'+figure['data']+') retain this page.\n\n'
    publish_chunks(root/'supplements/figures.md', (v.encode() for v in gallery()))
    return figures, caption
