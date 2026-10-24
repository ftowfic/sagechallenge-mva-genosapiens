"""Keep compatible study LINCS products outside replaceable presentation output."""
from pathlib import Path
import shutil

from genosapiens.track2_operator_inputs import require
from genosapiens.track2_scientific_storage import json_identity, file_identity, publish_json
from genosapiens.track2_scientific_reader import read_scientific


def science_dependencies(package):
    from genosapiens.track2_construction_dependencies import python_dependencies
    return python_dependencies(package,[('genosapiens.track2_lincs_integration','execute')])


def science_binding(package, plan, input_sha256):
    return dict(input_sha256=input_sha256,config=plan['lincs'],targets=plan['tables']['target_records'],
                code=science_dependencies(package))


def retain_study_science(root, package):
    root,package=Path(root),Path(package)
    report=root/'synthesis/controlled-work/report'
    proof_path=root/'synthesis/production.json'
    if not proof_path.is_file() or not (report/'lincs/results.json').is_file():return
    manifest=read_scientific(report/'manifest.json')
    if manifest['schema']!='genosapiens.track2-research-manifest/v1':return
    plan=read_scientific(report/'plan-source.json')
    if plan['lincs'] is None:return
    proof=read_scientific(proof_path)
    if proof['manifest_sha256']!=file_identity(report/'manifest.json')['sha256']:return
    dependencies=science_dependencies(package)
    if any(proof['binding']['code'].get(key)!=value for key,value in dependencies.items()):return
    from genosapiens.track2_synthesis import validate_study_report
    validate_study_report(report,presentation=False)
    binding=science_binding(package,plan,manifest['bindings']['input'])
    target=root/'study-science'/json_identity(binding)['sha256']
    if not (target/'science.json').is_file():
        # An interrupted copy has no completion record and no accepted output.
        if target.exists():
            shutil.rmtree(target)
        target.mkdir(parents=True,mode=0o700)
        shutil.copytree(report/'lincs',target/'lincs')
        publish_json(target/'science.json',dict(schema='genosapiens.track2-retained-study-science/v1',
            binding=binding,result_sha256=file_identity(target/'lincs/results.json')['sha256'],
            figure_manifest_sha256=file_identity(target/'lincs/manifest.json')['sha256'],
            source_bindings={k:v for k,v in manifest['bindings'].items() if k.startswith('lincs:')}))


def retained_study_science(cache, package, plan, input_sha256):
    if cache is None:return None
    binding=science_binding(package,plan,input_sha256)
    root=Path(cache)/json_identity(binding)['sha256']
    if not (root/'science.json').is_file():return None
    record=read_scientific(root/'science.json')
    require(record['schema']=='genosapiens.track2-retained-study-science/v1' and record['binding']==binding,
            'Retained study science has incompatible dependencies.')
    require(record['result_sha256']==file_identity(root/'lincs/results.json')['sha256']
            and record['figure_manifest_sha256']==file_identity(root/'lincs/manifest.json')['sha256'],
            'Retained study science changed.')
    result=read_scientific(root/'lincs/results.json')
    from genosapiens.track2_lincs_report import validate_figures
    validate_figures(root,result)
    # The composer updates top-level decisions. Scientific arrays stay lazy.
    return result.copy(),root,record['source_bindings']
