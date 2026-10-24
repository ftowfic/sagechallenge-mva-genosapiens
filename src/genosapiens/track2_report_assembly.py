"""Attach checked overview assets before operator publication, without discovery."""
from pathlib import Path
import shutil
from itertools import chain

from genosapiens.track2_operator_inputs import require
from genosapiens.track2_scientific_reader import read_scientific
from genosapiens.track2_scientific_storage import publish_json, publish_chunks, file_identity, Records
from genosapiens.track2_report_overviews import overview_binding


def presentation_ready(root):
    path = Path(root)/'presentation.json'
    if not path.is_file():return True
    model=read_scientific(path)
    if model.get('schema')=='genosapiens.track2-report-presentation/v3':
        from genosapiens.track2_challenge_model import check_challenge_model
        check_challenge_model(model)
        return bool(model['figures']) and (Path(root)/'figures/manifest.json').is_file()
    return model['disclosure'].get('overview_state') == 'complete'


def attach_overviews(root, graphics):
    root, graphics = Path(root), Path(graphics)
    from genosapiens.track2_report_model import publish_model
    from genosapiens.track2_report_paper import paper_markdown
    from genosapiens.track2_scientific_records import stored_json
    model = dict(read_scientific(root/'presentation.json'))
    manifest = read_scientific(graphics/'manifest.json')
    require(manifest['schema'] == 'genosapiens.track2-report-overviews/v1'
        and manifest['binding'] == overview_binding(model), 'Overview graphics differ from the report selections.')
    for name, digest in manifest['outputs'].items():
        require(Path(name).name == name and file_identity(graphics/name)['sha256'] == digest,
                'An overview asset is missing or changed.')
    destination = root/'overviews'
    destination.mkdir(mode=0o700,exist_ok=True)
    for name in [*manifest['outputs'], 'manifest.json']:
        shutil.copyfile(graphics/name, destination/name)
    figures = []
    for i, record in enumerate(manifest['figures']):
        for suffix in ('.png','.svg'):
            require(record['id']+suffix in manifest['outputs'], 'An overview image is missing.')
        figures.append(dict(label='Figure '+str(record['number']),number=record['number'], main=True,
            caption=record['caption'],alt=record['caption'],png='overviews/'+record['id']+'.png',
            svg='overviews/'+record['id']+'.svg',data='overviews/overview-data.json#/'+str(i),
            provenance=manifest['binding']))
    old = model['figures']
    model['figures'] = stored_json(Records(chain(figures, (f for f in old if not f['main']))))
    model['disclosure'] = dict(model['disclosure'],overview_state='complete')
    supplement = dict(title='Main overview figures and source data',path='overviews/manifest.json')
    model['supplements'] = [*model['supplements']] if supplement in model['supplements'] else [*model['supplements'],supplement]
    with (root/'report.md').open() as source:
        title = source.readline().rstrip('\n')
    paper = ''.join(paper_markdown(model,title,read_scientific(root/'readiness.json'))).encode()
    publish_model(root/'presentation.json', model)
    (root/'report.md').chmod(0o600)
    publish_chunks(root/'report.md',[paper])
    (root/'manifest.json').chmod(0o600)
    report_manifest = dict(read_scientific(root/'manifest.json'))
    if report_manifest['schema'] == 'genosapiens.track2-scientific-completion/v1':
        from genosapiens.track2_scientific_completion import publish_completion
        publish_completion(root,stage='native-report',binding=report_manifest['binding'],filename='manifest.json',result=report_manifest['result'])
    else:
        names=['report.md','presentation.json', *('overviews/'+name for name in [*manifest['outputs'],'manifest.json'])]
        report_manifest['outputs'] = dict(report_manifest['outputs'])
        report_manifest['outputs'].update({name:file_identity(root/name)['sha256'] for name in names})
        publish_json(root/'manifest.json',report_manifest)
    from genosapiens.track2_report_checks import check_paper
    check_paper(root)
    return file_identity(root/'manifest.json')['sha256']


def complete_operator_presentation(root, recovery_binding=None):
    from genosapiens.track2_report_finalization import save_production, checked_synthesis, production_matches
    from genosapiens.track2_run_state import atomic_record
    root = Path(root)
    report = root/'synthesis/controlled-work/report'
    journal_path = root/'synthesis/presentation-assembly.json'
    binding = recovery_binding or read_scientific(root/'current-routing.json')['report_binding']
    if journal_path.is_file():
        journal = read_scientific(journal_path)
        require(production_matches(journal['binding'],binding), 'Presentation recovery dependencies changed.')
        require(file_identity(root/'synthesis/overviews/manifest.json')['sha256'] == journal['overview_sha256'],
                'Presentation recovery graphics changed.')
        result = dict(journal['result'])
    else:
        checked_synthesis(root,binding)
        result = dict(read_scientific(root/'synthesis/controlled-work/result.json'))
        journal = dict(binding=binding,result=result,
            overview_sha256=file_identity(root/'synthesis/overviews/manifest.json')['sha256'])
        atomic_record(journal_path,journal)
    result['manifest_sha256'] = attach_overviews(report,root/'synthesis/overviews')
    from genosapiens.track2_report_metadata import checked_report_metadata
    checked_report_metadata(report)
    atomic_record(root/'synthesis/controlled-work/result.json',result)
    save_production(root/'synthesis',binding)
    journal_path.unlink()
    return 0
