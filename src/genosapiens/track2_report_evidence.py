"""Scientific exports and admission independent of report prose."""
from pathlib import Path

from genosapiens.track2_operator_inputs import require
from genosapiens.track2_report_contracts import EVIDENCE_SCHEMA, check_evidence_contract
from genosapiens.track2_scientific_storage import publish_json, file_identity, json_identity
from genosapiens.track2_scientific_reader import read_scientific


def _reference(root, path):
    return dict(path=path,sha256=file_identity(Path(root)/path)['sha256'])


def publish_native_evidence(root, *, model, bindings, readiness, snapshot_sha256,
                            discovery_bindings, source_versions):
    """Keep producer decisions and complete supplements without a legacy paper."""
    from genosapiens.track2_report_sources import source_versions as structured_versions
    publish_json(Path(root)/'scientific-bindings.json',dict(discovery_bindings=discovery_bindings))
    return _publish_evidence(root,model=model,bindings=bindings,readiness=readiness,
        snapshot_sha256=snapshot_sha256,
        source_versions=structured_versions(source_versions,origin='producer_structured',kind='native'),kind='native')


def publish_study_evidence(root, *, model, scientific_metadata):
    """Preserve the study producer's original review and readiness semantics."""
    metadata={k:v for k,v in scientific_metadata.items() if k not in {'schema','outputs'}}
    publish_json(Path(root)/'scientific-bindings.json',dict(
        schema='genosapiens.track2-study-evidence-source/v1',**metadata))
    from genosapiens.track2_report_sources import source_versions
    versions=source_versions([dict(source='study',bindings=metadata['bindings'],runtime=metadata['runtime'])],
        origin='retained_study',kind='study')
    return _publish_evidence(root,model=model,bindings=metadata['bindings'],readiness=metadata['readiness'],
        snapshot_sha256=metadata.get('input_snapshot_sha256'),source_versions=versions,kind='study')


def _publish_evidence(root, *, model, bindings, readiness, snapshot_sha256, source_versions, kind):
    from genosapiens.track2_challenge_revision import inventory
    root=Path(root)
    require(not (root/'report.md').exists() and not (root/'presentation.json').exists(),
            'Fresh evidence cannot contain a historical report presentation.')
    from genosapiens.track2_report_sources import frozen_eligibility
    payloads={
        'historical-selections.json':model['selections'],
        'report-methods.json':dict(methods=model['disclosure']['methods'],
            evidence_origin=model['disclosure'].get('evidence_origin')),
        'reviews.json':{key:model['disclosure'].get(key) for key in ('input_review','plan_review')},
        'source-versions.json':source_versions,
        'eligibility.json':frozen_eligibility(root,source_versions,kind=kind),
        'decisions.json':{key:read_scientific(root/(key+'.json')) for key in ('assessment','tables','outcomes')},
        'supplement-layout.json':dict(schema='genosapiens.track2-scientific-supplements/v1',
            tissue=model['schema']=='genosapiens.track2-report-presentation/v2',policy=model['policy'],
            policy_sha256=model['policy_sha256'],sources=model['sources'],
            captions=model['captions'],figures=model['figures'],supplements=model['supplements'])}
    for name,value in payloads.items():publish_json(root/name,value)
    paths=dict(plan='plan-source.json',context='supplements/context.json',
        objective_scores='supplements/objective-scores.json',historical_selections='historical-selections.json',
        coverage='supplements/coverage.json',decisions='decisions.json',methods='report-methods.json',reviews='reviews.json',
        source_versions='source-versions.json',eligibility='eligibility.json',figure_index='supplements/figures.json',
        scientific_details='scientific-details.md')
    manifest=dict(schema=EVIDENCE_SCHEMA,classification=model['classification'],source_kind=kind,
        input_origin=model['input_origin'],bindings=bindings,input_snapshot_sha256=snapshot_sha256,
        readiness=readiness,records={key:_reference(root,path) for key,path in paths.items()},
        historical_source=None,outputs=inventory(root))
    check_evidence(root,manifest)
    publish_json(root/'manifest.json',manifest)
    return dict(**readiness,manifest_sha256=file_identity(root/'manifest.json')['sha256'])


def scientific_model(root, manifest, *, _read=read_scientific):
    """Reconstruct validation inputs without publishing an old presentation."""
    root=Path(root)
    if manifest['historical_source'] is not None:
        original=_read(root/'historical/manifest.json')
        if original.get('schema')=='genosapiens.track2-challenge-package/v2':
            return scientific_model(root/'historical/evidence',_read(root/'historical/evidence/manifest.json'),_read=_read)
        return _read(scientific_root(root,manifest)/'presentation.json')
    def record(key):return _read(root/manifest['records'][key]['path'])
    layout=_read(root/'supplement-layout.json')
    require(layout['schema']=='genosapiens.track2-scientific-supplements/v1',
            'The scientific supplement contract is unsupported.')
    methods=record('methods');reviews=record('reviews')
    return dict(schema='genosapiens.track2-report-presentation/v2' if layout['tissue'] else
        'genosapiens.track2-report-presentation/v1',policy=layout['policy'],policy_sha256=layout['policy_sha256'],
        input_origin=manifest['input_origin'],classification=manifest['classification'],sources=layout['sources'],
        context=record('context'),scores=record('objective_scores'),selections=record('historical_selections'),
        captions=layout['captions'],figures=layout['figures'],supplements=layout['supplements'],coverage=record('coverage'),
        disclosure=dict(methods=methods['methods'],evidence_origin=methods['evidence_origin'],
            **reviews,overview_state='pending'))


def check_native_evidence(root, manifest, *, page_context=None):
    from genosapiens.track2_discovery_report import validate_native_science
    from genosapiens.track2_report_checks import check_scientific_supplements
    root=Path(root)
    check_evidence_contract(manifest,root=root)
    require(manifest['source_kind']=='native', 'The native evidence adapter received another source kind.')
    binding=read_scientific(root/'scientific-bindings.json')
    readiness=validate_native_science(root,dict(bindings=manifest['bindings'],
        input_snapshot_sha256=manifest['input_snapshot_sha256'],
        discovery_bindings=binding['discovery_bindings'],readiness=manifest['readiness']),page_context=page_context)
    decisions=read_scientific(root/manifest['records']['decisions']['path'])
    require(set(decisions)=={'assessment','tables','outcomes'} and all(
        json_identity(decisions[key])==json_identity(read_scientific(root/(key+'.json'))) for key in decisions),
        'The evidence decisions differ from the scientific exports.')
    model=scientific_model(root,manifest)
    check_scientific_supplements(root,model)
    return readiness


def check_study_evidence(root, manifest):
    from genosapiens.track2_synthesis import validate_study_science
    from genosapiens.track2_report_checks import check_scientific_supplements
    root=Path(root)
    check_evidence_contract(manifest,root=root)
    require(manifest['source_kind']=='study', 'The study evidence adapter received another source kind.')
    metadata=read_scientific(root/'scientific-bindings.json')
    require(metadata['schema']=='genosapiens.track2-study-evidence-source/v1'
        and metadata['bindings']==manifest['bindings'] and metadata['classification']==manifest['classification']
        and metadata['readiness']==manifest['readiness']
        and metadata.get('input_snapshot_sha256')==manifest['input_snapshot_sha256'],
        'The study evidence changed its original bindings or outcome.')
    readiness=validate_study_science(root,dict(metadata,outputs=manifest['outputs']))
    decisions=read_scientific(root/manifest['records']['decisions']['path'])
    require(set(decisions)=={'assessment','tables','outcomes'} and all(
        json_identity(decisions[key])==json_identity(read_scientific(root/(key+'.json'))) for key in decisions),
        'The evidence decisions differ from the scientific exports.')
    check_scientific_supplements(root,scientific_model(root,manifest))
    return readiness


def check_evidence(root, manifest, *, page_context=None):
    check_evidence_contract(manifest)
    from genosapiens.track2_report_sources import check_source_records
    check_source_records(root,manifest,scientific_root=scientific_root(root,manifest))
    if manifest['historical_source'] is not None:
        return check_historical_evidence(root,manifest)
    if manifest['source_kind']=='native':
        return check_native_evidence(root,manifest,page_context=page_context)
    return check_study_evidence(root,manifest)


def scientific_root(root, manifest):
    """Resolve the original scientific files inside a historical wrapper."""
    root=Path(root)
    if manifest['historical_source'] is None:return root
    source=root/'historical'
    original=read_scientific(source/'manifest.json')
    if original.get('schema')=='genosapiens.track2-challenge-package/v2':
        return scientific_root(source/'evidence',read_scientific(source/'evidence/manifest.json'))
    if original.get('schema')=='genosapiens.track2-challenge-package/v1':source=source/'source'
    return source


def _historical_contract(root, metadata, *, source_validation=None):
    from genosapiens.track2_challenge_revision import inventory
    root=Path(root)
    parent=read_scientific(root/'historical/manifest.json')
    source=root/'historical'
    if parent.get('schema')=='genosapiens.track2-challenge-package/v1':source=source/'source'
    modern=parent.get('schema')=='genosapiens.track2-challenge-package/v2'
    inner=read_scientific(source/'evidence/manifest.json') if modern else None
    model=scientific_model(source/'evidence',inner) if modern else read_scientific(source/'presentation.json')
    require(model['schema'] in {'genosapiens.track2-report-presentation/v1','genosapiens.track2-report-presentation/v2'},
            'The historical evidence lacks a supported scientific model.')
    prefix=source.relative_to(root).as_posix()+'/'
    paths=dict(plan='plan-source.json',context='supplements/context.json',
        objective_scores='supplements/objective-scores.json',historical_selections='presentation.json',
        coverage='supplements/coverage.json',decisions='outcomes.json',methods='presentation.json',reviews='presentation.json',
        source_versions='scientific-details.md',eligibility='plan-source.json',
        figure_index='supplements/figures.json',scientific_details='scientific-details.md')
    records=({key:_reference(root,prefix+'evidence/'+ref['path']) for key,ref in inner['records'].items()}
        if modern else {key:_reference(root,prefix+path) for key,path in paths.items()})
    records.update(source_versions=_reference(root,'source-versions.json'),eligibility=_reference(root,'eligibility.json'))
    historical=dict(manifest=_reference(root,'historical/manifest.json'),report=_reference(root,'historical/report.md'))
    if source_validation is not None:historical['validation']=source_validation
    return dict(schema=EVIDENCE_SCHEMA,classification=metadata.classification,source_kind=metadata.kind,
        input_origin=model['input_origin'],bindings=metadata.bindings,input_snapshot_sha256=metadata.input_snapshot_sha256,
        readiness=metadata.readiness,records=records,
        historical_source=historical,
        outputs=inventory(root,omit={'manifest.json'}))


def adapt_report_evidence(source, destination, *, expected_manifest_sha256, for_regeneration=False):
    """Admit historical science and preserve the supplied source files.

    Regeneration checks scientific dependencies instead of the old archive.
    The new wrapper records this choice and binds its own retained files.
    """
    import shutil
    from genosapiens.track2_challenge_revision import inventory
    from genosapiens.track2_report_metadata import checked_report_metadata
    source=Path(source).absolute();destination=Path(destination).absolute()
    require(source.resolve()==source and destination.resolve()==destination and not destination.exists()
        and not destination.is_relative_to(source) and not source.is_relative_to(destination),
        'Choose a new independent evidence directory.')
    before=inventory(source)
    metadata=checked_report_metadata(source,expected_manifest_sha256=expected_manifest_sha256,for_regeneration=for_regeneration)
    require(inventory(source)==before, 'The historical report changed during admission.')
    destination.mkdir(mode=0o700,parents=True)
    target=destination/'historical';target.mkdir(mode=0o700)
    for name in before:
        path=target/name;path.parent.mkdir(mode=0o700,parents=True,exist_ok=True)
        shutil.copyfile(source/name,path);path.chmod(0o600)
    require(inventory(target)==before and inventory(source)==before,
            'Historical evidence changed during copying.')
    from genosapiens.track2_report_sources import historical_source_versions,frozen_eligibility
    original=read_scientific(target/'manifest.json')
    if original.get('schema')=='genosapiens.track2-challenge-package/v2':
        inner=read_scientific(target/'evidence/manifest.json')
        science=scientific_root(target/'evidence',inner)
        versions=read_scientific(target/'evidence'/inner['records']['source_versions']['path'])
    else:
        science=target/'source' if original.get('schema')=='genosapiens.track2-challenge-package/v1' else target
        versions=historical_source_versions(science,kind=metadata.kind)
    publish_json(destination/'source-versions.json',versions)
    publish_json(destination/'eligibility.json',frozen_eligibility(science,versions,kind=metadata.kind))
    manifest=_historical_contract(destination,metadata,
        source_validation='scientific_dependencies' if for_regeneration else None)
    check_evidence_contract(manifest,root=destination)
    publish_json(destination/'manifest.json',manifest)
    return _reference(destination,'manifest.json')


def check_historical_evidence(root, manifest):
    from genosapiens.track2_report_metadata import checked_report_metadata
    root=Path(root)
    check_evidence_contract(manifest,root=root)
    require(manifest['historical_source']['manifest']['path']=='historical/manifest.json'
        and manifest['historical_source']['report']['path']=='historical/report.md',
        'The historical source location differs from its evidence contract.')
    metadata=checked_report_metadata(root/'historical',
        expected_manifest_sha256=manifest['historical_source']['manifest']['sha256'],
        for_regeneration=manifest['historical_source'].get('validation')=='scientific_dependencies')
    from genosapiens.track2_report_sources import historical_source_versions
    parent=read_scientific(root/'historical/manifest.json')
    if parent.get('schema')=='genosapiens.track2-challenge-package/v2':
        inner=read_scientific(root/'historical/evidence/manifest.json')
        expected_versions=read_scientific(root/'historical/evidence'/inner['records']['source_versions']['path'])
    else:expected_versions=historical_source_versions(scientific_root(root,manifest),kind=metadata.kind)
    require(json_identity(read_scientific(root/'source-versions.json'))==json_identity(expected_versions),
        'The adapted source versions differ from the historical report.')
    require(json_identity(manifest)==json_identity(_historical_contract(root,metadata,
        source_validation=manifest['historical_source'].get('validation'))),
            'The evidence wrapper changed its original scientific records.')
    return metadata.readiness
