"""Presentation consistency checks over retained report products only."""
from pathlib import Path
import re

from genosapiens.track2_operator_inputs import require
from genosapiens.track2_report_model import check_model, relative_destination
from genosapiens.track2_report_ranking import select_tables
from genosapiens.track2_scientific_reader import read_scientific
from genosapiens.track2_scientific_storage import file_identity, json_identity


def check_paper(root):
    root = Path(root)
    model = check_model(read_scientific(root/'presentation.json'))
    import json
    from genosapiens.json_schema import validate_json_schema
    schema_path=Path(__file__).resolve().parents[2]/'schemas/track2/report-presentation.schema.json'
    validate_json_schema(model,json.loads(schema_path.read_bytes()))
    check_scientific_supplements(root, model)
    if model['disclosure']['overview_state']=='complete':
        from genosapiens.track2_report_overviews import overview_binding
        overview=read_scientific(root/'overviews/manifest.json')
        require(overview['binding']==overview_binding(model), 'An overview has stale presentation dependencies.')
        require(all(file_identity(root/'overviews'/name)['sha256']==digest for name,digest in overview['outputs'].items()),
                'An overview asset differs from its renderer receipt.')
    from genosapiens.track2_report_paper import paper_markdown
    actual = (root/'report.md').read_text()
    expected_text = ''.join(paper_markdown(model, actual.splitlines()[0], read_scientific(root/'readiness.json')))
    require(actual == expected_text, 'Report prose differs from its presentation model.')
    positions = [actual.index('## '+heading+'\n') for heading in model['policy']['sections']]
    require(positions == sorted(positions), 'The report section order differs.')
    for link in re.findall(r'\]\(([^)]+)\)', actual):
        relative_destination(link)
        path, _, fragment = link.partition('#')
        require(not path or (root/path).is_file(), 'A report destination is missing.')
        if not path and fragment:
            require('id="'+fragment+'"' in actual, 'A report context anchor is missing.')
    return model


def check_scientific_supplements(root, model):
    """Replay export arithmetic and references without requiring report prose."""
    root = Path(root)
    check_model(model)
    require(model['sources']['plan'] == file_identity(root/'plan-source.json')['sha256'], 'Presentation source plan differs.')
    expected = select_tables(model['scores'],tissue=model['schema']=='genosapiens.track2-report-presentation/v2')
    require(json_identity(expected) == json_identity(model['selections']), 'Presentation selections differ from objective scores.')
    for row in model['scores']:
        c = row['components']
        if 'tissue_score' in row:
            from fractions import Fraction
            from genosapiens.track2_tissue_records import validate_tissue_record
            from genosapiens.track2_tissue_ranking import exact_score
            tissue=validate_tissue_record(row['tissue_score'],'score')
            require(all(c[k]==tissue[k] for k in ('P','C','N')) and c['M']==len(row['counted_links']),
                'Presentation tissue counts differ from retained evidence.')
            expected_lincs=Fraction(tissue['numerator'],tissue['N']) if tissue['N'] else None
            require(exact_score(row,'S_L')==expected_lincs and exact_score(row,'S_M')==int(c['M']>0)+2*c['M']
                and exact_score(row,'S')==(expected_lincs or Fraction(0))+exact_score(row,'S_M'),
                'Presentation tissue score arithmetic differs.')
            require(all(row['scores'][k]==(float(exact_score(row,k)) if exact_score(row,k) is not None else None)
                for k in ('S_L','S_M','S')), 'Presentation tissue decimal values differ.')
            continue
        require(c['B'] == len(row['counted_cells']['B']) and c['U'] == len(row['counted_cells']['U'])
                and c['C'] == len(row['counted_cells']['C']) and c['M'] == len(row['counted_links']),
                'Presentation score counts differ from contributions.')
        require(c['I_L'] == int(c['B']+c['U'] > 0) and c['I_M'] == int(c['M'] > 0)
                and row['E'] == c['I_L']+c['I_M'], 'Presentation evidence-type count differs.')
        half = dict(S_L=2*c['I_L']+4*c['B']+c['U']-c['C'], S_M=2*c['I_M']+4*c['M'])
        half['S'] = half['S_L']+half['S_M']
        require(half == row['half_points'] and {k:v/2 for k,v in half.items()} == row['scores'],
                'Presentation score arithmetic differs.')
    for kind in expected:
        require(json_identity(read_scientific(root/'supplements'/(kind+'-ranking.json'))) ==
                json_identity(expected[kind]['ranking']), 'The complete ranking export differs.')
    require(json_identity(read_scientific(root/'supplements/objective-scores.json')) == json_identity(model['scores']),
            'The objective score export differs.')
    require(json_identity(read_scientific(root/'supplements/column-dictionaries.json')) == json_identity(model['captions']),
            'A presentation column dictionary differs.')
    from genosapiens.track2_report_captions import annotate_tables
    cell_ids=list(read_scientific(root/'lincs/results.json')['cells']) if (root/'lincs/results.json').is_file() else []
    for path in [*root.glob('*.md'),*(root/'lincs').glob('*.md')]:
        if path.name=='report.md':continue
        text=path.read_text()
        annotated,captions=annotate_tables(text,str(path.relative_to(root)),cell_ids)
        require(text==annotated and all(model['captions'].get(key)==value for key,value in captions.items()),
                'Supplementary table captions differ from their column dictionaries.')
    import csv
    for entry in read_scientific(root/'supplements/export-index.json'):
        path=root/relative_destination(entry['path'])
        require(path.is_file() and all(key in model['captions'] for key in entry['captions']),
                'A complete export or its column dictionary is missing.')
        if path.suffix=='.csv':
            with path.open() as stream:columns=next(csv.reader(stream))
            require(any(model['captions'][key]['columns']==columns for key in entry['captions']),
                    'An export header differs from its column dictionary.')
    from genosapiens.track2_report_supplements import local_figures
    from genosapiens.track2_scientific_storage import Records
    plan = read_scientific(root/'plan-source.json')
    local = (f for f in model['figures'] if not f['main'])
    require(json_identity(Records(local_figures(root, plan))) == json_identity(Records(local)),
            'The complete local figure index differs.')
    for figure in model['figures']:
        for key in ('png','svg','data'):
            path = relative_destination(figure[key]).split('#',1)[0]
            require((root/path).is_file(), 'A presentation figure asset is missing.')
    return model
