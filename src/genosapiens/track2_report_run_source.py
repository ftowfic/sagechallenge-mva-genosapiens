"""Read-only admission of completed runs without opening original inputs."""
from pathlib import Path

from genosapiens.track2_operator_inputs import require,small_json
from genosapiens.track2_scientific_storage import file_identity
from genosapiens.track2_report_metadata import checked_report_metadata


def admit_completed_run(root, *, recover_unfinished=False):
    from genosapiens.track2_report_producer import admit_producer
    root=Path(root).absolute()
    require(root.resolve()==root, 'The source run contains a path link.')
    report=root/'synthesis/controlled-work/report'
    proof_path=root/'synthesis/production.json'
    if not proof_path.is_file():proof_path=root/'retained-report-routing.json'
    require(proof_path.is_file(), 'The source run lacks retained production proof.')
    proof=small_json(proof_path)
    require(proof.get('schema')=='genosapiens.track2-report-production/v1', 'The source production schema is unsupported.')
    metadata=checked_report_metadata(report,expected_manifest_sha256=proof['manifest_sha256'],for_regeneration=True)
    fields={'snapshot_sha256','code','gene_metadata','overview_renderer'}
    if metadata.kind=='native':fields|={'discovery_sha256','figures_completion_sha256','discovery_completion_sha256','renderer'}
    require(set(proof['binding'])==fields, 'The production dependency record is incomplete or unsupported.')
    producer=admit_producer(proof_path,expected_production_sha256=file_identity(proof_path)['sha256'],
        expected_manifest_sha256=metadata.manifest_sha256,native=metadata.kind=='native')
    result_path=root/'synthesis/controlled-work/result.json';result=small_json(result_path)
    require(set(result)==set(metadata.readiness)|{'manifest_sha256'}|({'output'} if metadata.kind=='study' else set())
        and result['manifest_sha256']==metadata.manifest_sha256
        and all(result[k]==v for k,v in metadata.readiness.items())
        and proof['result_sha256']==file_identity(result_path)['sha256'],
        'The saved scientific result differs from its production proof.')
    route=small_json(root/'current-routing.json')
    require(route['report_binding']==proof['binding'] and route['native']==(metadata.kind=='native'),
        'The source routing differs from original production.')
    completed=route.get('report_complete') is True
    require(completed or recover_unfinished and route.get('report_complete') is False,
        'The source run has no matching completed report routing.')
    completion_names=('completions/report.json','synthesis/completion.json')
    completions=[small_json(root/name) for name in completion_names if (root/name).exists()]
    for completion in completions:
        require(completion['binding']['production']==proof['binding']
            and completion['result']==dict(readiness=metadata.readiness,manifest_sha256=metadata.manifest_sha256),
            'The source completion differs from production or scientific readiness.')
    require(not completions or all(value==completions[0] for value in completions),
        'The source run completion differs from the retained report completion.')
    from genosapiens.track2_terminal import assessment_outcome
    if completed:
        require(len(completions)==2 and route.get('synthesis_complete') is True,
            'The source run lacks complete report evidence.')
        require(route.get('finalization_exit_code')==assessment_outcome(metadata.readiness).exit_code,
            'The source terminal outcome differs from its scientific result.')
    else:
        require('finalization_exit_code' not in route,
            'The unfinished source claims a terminal report outcome.')
        require(producer['rule']!='current-report-producer-v1',
            'Unfinished report recovery requires an explicitly supported historical producer.')
    snapshot=small_json(root/'preparation/snapshot.json')
    require(file_identity(root/'preparation/snapshot.json')['sha256']==proof['binding']['snapshot_sha256'],
        'The source preparation snapshot differs from production.')
    if metadata.input_snapshot_sha256 is not None:
        require(metadata.input_snapshot_sha256==proof['binding']['snapshot_sha256']
            and metadata.bindings==snapshot['bindings'] and metadata.classification==snapshot['classification'],
            'The source evidence differs from its preparation snapshot.')
    else:
        from genosapiens.track2_freeze import load_frozen_public_envelope
        frozen=load_frozen_public_envelope(root/'preparation/snapshot.json')
        require(metadata.classification=='public' and metadata.bindings==frozen['manifest_checksums'],
            'The source study differs from its frozen envelope.')
    require(small_json(root/'construction/run.json')['input_sha256']==metadata.bindings['input'],
        'The retained construction identity differs from the report.')
    scientific_names=[]
    if metadata.kind=='native':
        for name,key in (('discovery/discovery-results.json','discovery_sha256'),
                         ('discovery-figures/completion.json','figures_completion_sha256'),
                         ('discovery/completion.json','discovery_completion_sha256')):
            if proof['binding'][key] is None:
                require(not (root/name).exists(), 'The native completion is missing or unexpected.')
            else:
                require(file_identity(root/name)['sha256']==proof['binding'][key],
                        'Retained native science differs from original production.')
                scientific_names.append(name)
    names=('synthesis/controlled-work/result.json',
        'current-routing.json','preparation/snapshot.json','construction/run.json',proof_path.relative_to(root).as_posix())
    admission=dict(schema='genosapiens.track2-report-run-admission/v1',report=str(report),
        source_validation='scientific_dependencies',
        report_manifest_sha256=metadata.manifest_sha256,classification=metadata.classification,
        source_kind=metadata.kind,readiness=metadata.readiness,producer_admission=producer,
        metadata={name:file_identity(root/name)['sha256'] for name in (*names,*scientific_names)})
    if completed:
        admission['metadata'].update({name:file_identity(root/name)['sha256'] for name in completion_names})
    else:
        admission.update(source_report_completion='unfinished',
            completion_records={name:file_identity(root/name)['sha256'] if (root/name).exists() else None
                for name in completion_names})
    return admission


def check_source_metadata(root, admission):
    root=Path(root)
    require(all(file_identity(root/name)['sha256']==digest for name,digest in admission['metadata'].items()),
            'Source run metadata changed during reporting.')
    require(file_identity(Path(admission['report'])/'manifest.json')['sha256']==admission['report_manifest_sha256'],
            'The source report changed during reporting.')
    require(all((file_identity(root/name)['sha256'] if (root/name).exists() else None)==digest
        for name,digest in admission.get('completion_records',{}).items()),
        'Source completion evidence changed during reporting.')
