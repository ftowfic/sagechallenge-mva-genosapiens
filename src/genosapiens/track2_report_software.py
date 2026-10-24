"""A reporting software freeze independent of the original scientific package."""
from pathlib import Path
import shutil
import uuid

from genosapiens.track2_operator_inputs import require
from genosapiens.track2_package_files import runtime_files,SCRIPTS
from genosapiens.track2_challenge_revision import inventory
from genosapiens.track2_scientific_reader import read_scientific
from genosapiens.track2_scientific_storage import publish_json,file_identity


def freeze_reporter(repository,destination, *, retained_files=None):
    root=Path(repository).absolute();destination=Path(destination).absolute()
    require(destination.resolve()==destination and not destination.exists(), 'Choose a new reporting software directory.')
    destination.mkdir(mode=0o700,parents=True)
    for source in runtime_files(root):
        target=destination/source.relative_to(root);target.parent.mkdir(mode=0o700,parents=True,exist_ok=True)
        shutil.copyfile(source,target);target.chmod(0o700 if source.relative_to(root).as_posix() in SCRIPTS else 0o600)
    for name,source in (retained_files or {}).items():
        require(Path(name).name==name and not (destination/name).exists(),
            'A retained reporting input cannot replace runtime software.')
        shutil.copyfile(source,destination/name);(destination/name).chmod(0o600)
    runtime=read_scientific(destination/'config/track2/discovery-runtime.json')['image']
    value=dict(schema='genosapiens.track2-reporting-software/v1',runtime_image=runtime,files=inventory(destination))
    publish_json(destination/'reporter.json',value)
    expected=file_identity(destination/'reporter.json')['sha256']
    check_reporter(destination,expected)
    return dict(sha256=expected,runtime_image=runtime)


def check_reporter(root,expected):
    root=Path(root)
    require(file_identity(root/'reporter.json')['sha256']==expected, 'The reporting software identity changed.')
    value=read_scientific(root/'reporter.json')
    require(value['schema']=='genosapiens.track2-reporting-software/v1'
        and value['files']==inventory(root,omit={'reporter.json'})
        and value['runtime_image']==read_scientific(root/'config/track2/discovery-runtime.json')['image'],
        'The reporting software files or runtime changed.')
    return value


def select_reporter(repository, operation, *, refresh=False, published=False):
    """Retain old software when refreshing an operation that failed admission."""
    operation=Path(operation)
    selection=operation/'reporter-selection.json'
    original=operation/'reporter-identity.json'
    if selection.exists():
        selected=read_scientific(selection)
        relative=Path(selected['software'])
        require(len(relative.parts)==3 and relative.parts[0]=='reporters'
            and relative.parts[2]=='software' and '..' not in relative.parts,
            'The reporting software selection has an invalid path.')
        software=operation/relative
        reporter=selected['reporter']
    elif original.exists():
        software=operation/'software';reporter=read_scientific(original)
    else:
        require(not refresh, 'Reporter refresh requires an existing failed report operation.')
        software=operation/'software';reporter=freeze_reporter(repository,software)
        publish_json(original,reporter)
        return software,reporter
    check_reporter(software,reporter['sha256'])
    if not refresh:return software,reporter
    require(not published and not any((operation/name).exists() for name in (
        'source-admission.json','stages','attempts','workflow','package','result.json','current-routing.json')),
        'Reporter refresh requires failure before source admission. Existing work remains preserved.')
    replacement=operation/'reporters'/uuid.uuid4().hex/'software'
    current=freeze_reporter(repository,replacement)
    record=dict(software=replacement.relative_to(operation).as_posix(),reporter=current)
    publish_json(replacement.parent/'selection.json',dict(selected=record,
        previous=dict(software=software.relative_to(operation).as_posix(),reporter=reporter)))
    publish_json(selection,record)
    return replacement,current


def public_reporting_software(repository, source, operation, run):
    """Select a current reporter for the exact retained public source baseline.

    Scientific phases keep the original software. Report export still checks
    discovery dependencies and semantics before consuming retained results.
    """
    from genosapiens.track2_construction_dependencies import python_dependencies
    from genosapiens.track2_report_producer import code_identity,PRECHANGE_NATIVE
    from genosapiens.track2_report_finalization import VALIDATORS
    from genosapiens.track2_report_gene_metadata import FILENAME
    repository=Path(repository);source=Path(source);operation=Path(operation)
    require(inventory(source)==run['software_files']
        and code_identity(run['software_files'])==run['software_sha256'],
        'The frozen public software differs from its original identity.')
    if 'src/genosapiens/track2_report_fresh.py' in run['software_files']:
        return source
    code=python_dependencies(source,[('genosapiens.track2_discovery_report','synthesize_native'),
        ('genosapiens.track2_report_assembly','attach_overviews')],ignored=VALIDATORS,explicit_empty_ast=True)
    require(code_identity(code)==PRECHANGE_NATIVE,
        'Unsupported unfinished public producer. Preserve the original software and checkpoints.')
    # Reject changed policies and runtimes before any scientific worker starts.
    for name in ('source-domain-matrix-v1.json','decision-policy-g0-g7-v1.json',
                 'public-expansion-policy-v1.json','exploratory-policy-v2.json',
                 'diagnosis-runtime.json','discovery-runtime.json'):
        relative=Path('config/track2')/name
        require(file_identity(source/relative)==file_identity(repository/relative),
            'The public reporting transition changed a scientific policy or runtime.')
    require(all(run['images'][kind]==read_scientific(source/('config/track2/'+kind+'-runtime.json'))['image']
        for kind in ('diagnosis','discovery')), 'The public run runtime differs from its frozen software.')
    retained={name:source/name for name in ('discovery-resources.json',FILENAME)}
    require(all(path.is_file() for path in retained.values()), 'The public reporter lacks retained resource metadata.')
    binding=dict(schema='genosapiens.track2-public-report-transition/v1',
        scientific_software_sha256=run['software_sha256'],producer_code_sha256=PRECHANGE_NATIVE,
        retained_files={name:file_identity(path)['sha256'] for name,path in retained.items()})
    from genosapiens.track2_report_stages import stage
    def freeze(output):
        freeze_reporter(repository,output,retained_files=retained)
    software,_=stage(operation,'software',binding,freeze)
    check_reporter(software,file_identity(software/'reporter.json')['sha256'])
    return software
