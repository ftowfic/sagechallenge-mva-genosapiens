"""Separate challenge revisions over supported, self-contained source reports."""
import json
from pathlib import Path
import shutil
import stat

from genosapiens.track2_operator_inputs import require
from genosapiens.track2_scientific_reader import read_scientific
from genosapiens.track2_scientific_storage import file_identity, json_identity, publish_json, publish_chunks
from genosapiens.track2_challenge_model import build_challenge_model, check_reference
from genosapiens.track2_challenge_summary import summarize
from genosapiens.track2_challenge_selection import select_evidence

SCHEMA = 'genosapiens.track2-challenge-package/v1'


def inventory(root, *, omit=()):
    root=Path(root).absolute()
    require(root.resolve()==root, 'A challenge package path contains a symbolic link.')
    result={}
    for path in sorted(root.rglob('*')):
        require(not path.is_symlink(), 'A challenge package contains a symbolic link.')
        mode=path.stat().st_mode
        if stat.S_ISDIR(mode):continue
        require(stat.S_ISREG(mode), 'A challenge package contains an unsupported file type.')
        name=path.relative_to(root).as_posix()
        if name not in omit:result[name]=file_identity(path)['sha256']
    return result


def reference(root, name):
    return check_reference(dict(path=name,sha256=file_identity(Path(root)/name)['sha256']))


def _source_versions(source):
    """Read the supported native producer's exact retained source-version records."""
    text=(Path(source)/'scientific-details.md').read_text()
    require('## Source versions\n' in text, 'The source report lacks retained calibration source identities.')
    records={}
    for line in text.split('## Source versions\n',1)[1].split('\n## ',1)[0].splitlines():
        if line.startswith('`{') and line.endswith('}`'):
            row=json.loads(line[1:-1])
            if row.get('source') in {'fda','lincs'}:
                require(row['source'] not in records, 'A retained source identity is duplicated.')
                records[row['source']]=row
    require(set(records)=={'fda','lincs'}, 'Calibration requires the original FDA and LINCS source identities.')
    return dict(expected_fda_snapshot_sha256=records['fda']['sha256'], expected_source_release=records['lincs']['release'])


def _derived(source, legacy, *, calibration=None):
    if legacy['schema'] != 'genosapiens.track2-report-presentation/v2':
        require(calibration is None, 'Legacy cell reports cannot consume tissue calibration.')
        return None
    from genosapiens.track2_report_derived import derive_rankings, report_conditions, calibrate_rankings
    source=Path(source)
    queries=[read_scientific(path) for path in sorted((source/'lincs').glob('*/*/science/query.json'))]
    result=derive_rankings(conditions=report_conditions(source), queries=queries,
        evaluated_scores=read_scientific(source/'discovery/tissue-scores.json'),
        objectives=(o for h in legacy['context']['hypotheses'] for o in h['objectives']),
        historical_lincs=legacy['selections']['lincs']['ranking'])
    if calibration is not None:
        result=calibrate_rankings(result,**calibration,**_source_versions(source))
    return result


def _model(root, legacy, metadata, derived, calibration_reference, *, narrative_version=2,
           source_prefix='source', source_versions=None, scientific_manifest=None):
    root=Path(root)
    artifacts={key:reference(root,path.replace('source/',source_prefix+'/',1) if path.startswith('source/') else path) for key,path in (
        ('context','source/supplements/context.json'),('objective_scores','source/supplements/objective-scores.json'),
        ('historical_rankings','historical-rankings.json'),('coverage','source/supplements/coverage.json'),
        ('figure_index','source/supplements/figures.json'),('disclosures','disclosures.json'),
        ('scientific_details','source/scientific-details.md'))}
    if derived is not None:artifacts['derived_rankings']=reference(root,'derived-rankings.json')
    selected=select_evidence(legacy['selections'],derived=derived)
    versions=(source_versions if source_versions is not None else _source_versions(root/source_prefix)) if metadata.kind=='native' and selected else None
    for row in selected:
        row['source_reference']=source_prefix+'/'+row['source_reference']
        if row['mechanism']:row['mechanism']['source_reference']=source_prefix+'/'+row['mechanism']['source_reference']
        if metadata.kind=='native':
            row['approval_source']=dict(authority='FDA',snapshot_sha256=versions['expected_fda_snapshot_sha256'],
                source_reference=source_prefix+'/scientific-details.md#source-versions',
                jurisdiction=None,approval_date=None,approved_indication=None)
        else:row['approval_source']=None
    scope={}
    if metadata.kind=='study':
        scope=dict(route_names=read_scientific(root/'scientific-policy.json')['route_bounds'],
            required_routes=read_scientific(root/source_prefix/'coverage.json')['scope']['selected_routes'])
    summary=summarize(legacy['context'],legacy['coverage'],legacy['selections'],derived=derived,**scope)
    summary['evidence_origin']=legacy['disclosure'].get('evidence_origin')
    summary['methods']=legacy['disclosure']['methods']
    summary['reviews']={key:legacy['disclosure'].get(key) for key in ('input_review','plan_review')}
    supplements=[dict(title=key.replace('_',' ').capitalize(),artifact=value) for key,value in artifacts.items()]
    if (root/'producer-admission.json').is_file():
        supplements.append(dict(title='Original report producer',artifact=reference(root,'producer-admission.json')))
    if (root/'revision-parent.json').is_file():
        supplements.append(dict(title='Source revision identity',artifact=reference(root,'revision-parent.json')))
    from genosapiens.track2_derived_strength import policy
    return build_challenge_model(input_origin=legacy['input_origin'],classification=metadata.classification,
        source_kind=metadata.kind,scientific_policy=reference(root,'scientific-policy.json'),
        sources={key:reference(root,path) for key,path in (
            ('input',source_prefix+'/supplements/context.json'),('plan',source_prefix+'/plan-source.json'),
            ('scientific_manifest',scientific_manifest or source_prefix+'/manifest.json'))},
        artifacts=artifacts,summary=summary,
        selections=selected,figures=[],supplements=supplements,derived_policy=policy() if derived is not None else None,
        calibration=calibration_reference,narrative_version=narrative_version)


def publish_revision(source, destination, *, expected_source_manifest_sha256, calibration=None, context=None, figure_renderer=None, producer_admission=None):
    from genosapiens.track2_checked_evidence import CheckedEvidenceContext
    if context is None:
        with CheckedEvidenceContext() as owned:
            return _publish_revision(source,destination,expected_source_manifest_sha256=expected_source_manifest_sha256,
                calibration=calibration,context=owned,figure_renderer=figure_renderer,producer_admission=producer_admission)
    require(isinstance(context,CheckedEvidenceContext), 'Revision publication requires an invocation-owned context.')
    return _publish_revision(source,destination,expected_source_manifest_sha256=expected_source_manifest_sha256,
        calibration=calibration,context=context,figure_renderer=figure_renderer,producer_admission=producer_admission)


def _checked_derived(source, legacy, calibration, context):
    result=context.derived_report(source,legacy)
    if calibration is not None:
        require(result is not None, 'Legacy cell evidence cannot consume tissue calibration.')
        from genosapiens.track2_report_derived import calibrate_rankings
        result=calibrate_rankings(result,**calibration,**_source_versions(source))
    return result


def _publish_revision(source, destination, *, expected_source_manifest_sha256, calibration, context, figure_renderer, parent=None, producer_admission=None):
    """Preserve the source directory. Publish completion only after local checks."""
    from genosapiens.track2_report_metadata import checked_report_metadata
    from genosapiens.track2_report_model import check_model, publish_model
    from genosapiens.track2_report_paper import paper_markdown
    from genosapiens.track2_synthesis import policies
    from genosapiens.track2_published_genotype import canonical, sha256
    source=Path(source).absolute();destination=Path(destination).absolute()
    require(source.resolve()==source and destination.resolve()==destination and not destination.exists()
            and not destination.is_relative_to(source) and not source.is_relative_to(destination),
            'Choose a new independent challenge revision directory.')
    require(file_identity(source/'manifest.json')['sha256']==expected_source_manifest_sha256,
            'The source report differs from its independently expected identity.')
    if read_scientific(source/'manifest.json').get('schema')==SCHEMA:
        retained=_admit_revision_source(source,context)
        if calibration is None:calibration=retained['calibration']
        return _publish_revision(source/'source',destination,
            expected_source_manifest_sha256=retained['source_manifest_sha256'],calibration=calibration,
            context=context,figure_renderer=figure_renderer,
            parent=dict(path=source,identity=expected_source_manifest_sha256,inventory=inventory(source)),
            producer_admission=read_scientific(source/'producer-admission.json') if (source/'producer-admission.json').is_file() else producer_admission)
    before=inventory(source)
    metadata=context.checked_report(source)
    require(metadata.manifest_sha256==expected_source_manifest_sha256, 'The source report identity changed during admission.')
    legacy=check_model(read_scientific(source/'presentation.json'))
    require(legacy['schema'] in {'genosapiens.track2-report-presentation/v1','genosapiens.track2-report-presentation/v2'},
            'This revision transition requires a supported original presentation.')
    scientific_policy=policies()[0]
    require(sha256(canonical(scientific_policy))==metadata.bindings['policy'],
            'The source scientific policy lacks a supported transition.')
    require(inventory(source)==before, 'The source report changed during admission.')
    if producer_admission is not None:
        from genosapiens.track2_report_producer import check_retained_producer
        check_retained_producer(producer_admission,expected_source_manifest_sha256)
    destination.mkdir(mode=0o700,parents=True,exist_ok=False)
    if producer_admission is not None:publish_json(destination/'producer-admission.json',producer_admission)
    publish_json(destination/'started.json',dict(status='incomplete',source_manifest_sha256=expected_source_manifest_sha256))
    if parent is not None:
        shutil.copyfile(parent['path']/'manifest.json',destination/'revision-parent.json')
        (destination/'revision-parent.json').chmod(0o600)
    target=destination/'source';target.mkdir(mode=0o700)
    for name in before:
        path=target/name;path.parent.mkdir(mode=0o700,parents=True,exist_ok=True)
        shutil.copyfile(source/name,path);path.chmod(0o600)
    require(inventory(target)==before and inventory(source)==before, 'Source evidence changed during revision copying.')
    (destination/'scientific-policy.json').write_bytes(canonical(scientific_policy))
    publish_json(destination/'disclosures.json',legacy['disclosure'])
    publish_json(destination/'historical-rankings.json',{k:reference(destination,'source/supplements/'+k+'-ranking.json')
        for k in ('combined','lincs','mechanism')})
    derived=_checked_derived(target,legacy,calibration,context)
    calibration_reference=None
    if derived is not None:
        publish_json(destination/'derived-rankings.json',derived)
    if calibration is not None:
        from genosapiens.track2_report_calibration import publish_calibration
        observations=(r for r in derived['complete_export'] if r.get('numerator') is not None)
        calibration_reference=publish_calibration(destination/'calibration',observations=observations,
            classification=metadata.classification,context=context,**calibration)
    model=_model(destination,legacy,metadata,derived,calibration_reference)
    from genosapiens.track2_challenge_overviews import render_overviews, figure_references, pinned_runtime
    reuse_figures=False
    if parent is not None:
        from genosapiens.track2_challenge_overviews import check_overviews
        try:
            check_overviews(destination,model,parent['path']/'figures')
            reuse_figures=True
        except (ValueError,OSError,KeyError,TypeError):pass
    if reuse_figures:shutil.copytree(parent['path']/'figures',destination/'figures')
    else:render_overviews(destination,model,destination/'figures',runtime_image=pinned_runtime(),run_renderer=figure_renderer)
    model['figures']=figure_references(destination,model,destination/'figures')
    publish_model(destination/'presentation.json',model)
    publish_json(destination/'readiness.json',metadata.readiness)
    publish_chunks(destination/'report.md',(chunk.encode() for chunk in
        paper_markdown(model,'# Challenge research report',metadata.readiness,hypothesis_context=legacy['context'])))
    outputs=inventory(destination)
    manifest=dict(schema=SCHEMA,classification=metadata.classification,source_kind=metadata.kind,
        source_manifest_sha256=expected_source_manifest_sha256,source_bindings=metadata.bindings,
        input_snapshot_sha256=metadata.input_snapshot_sha256,readiness=metadata.readiness,
        scientific=reference(destination,'source/manifest.json'),calibration=calibration_reference,
        figures=reference(destination,'figures/manifest.json'),
        revision_parent=reference(destination,'revision-parent.json') if parent is not None else None,
        presentation=reference(destination,'presentation.json'),outputs=outputs)
    # Admission is separate from successful publication. An exception leaves the
    # new directory and started record without a successful manifest.
    _check_contents(destination,manifest,metadata=metadata,context=context)
    if parent is not None:
        require(reference(destination,'revision-parent.json')['sha256']==parent['identity']
            and inventory(parent['path'])==parent['inventory'],'The source revision changed during publication.')
    publish_json(destination/'manifest.json',manifest)
    return reference(destination,'manifest.json')


def _check_contents(root, manifest, *, metadata, context):
    from genosapiens.track2_report_model import check_model
    from genosapiens.track2_report_paper import paper_markdown
    from genosapiens.track2_published_genotype import canonical, sha256
    root=Path(root)
    require(set(manifest)-{'revision_parent'}=={'schema','classification','source_kind','source_manifest_sha256',
        'source_bindings','input_snapshot_sha256','readiness','scientific','calibration','figures','presentation','outputs'}
        and manifest['schema']==SCHEMA, 'The challenge package contract is unsupported.')
    require(reference(root,'source/manifest.json')['sha256']==metadata.manifest_sha256,
            'The copied source manifest differs from its admitted identity.')
    require(manifest['outputs']==inventory(root,omit={'manifest.json'}), 'The challenge package inventory differs.')
    require(manifest['source_manifest_sha256']==metadata.manifest_sha256
            and manifest['source_bindings']==metadata.bindings and manifest['source_kind']==metadata.kind
            and manifest['classification']==metadata.classification and manifest['readiness']==metadata.readiness
            and manifest['input_snapshot_sha256']==metadata.input_snapshot_sha256,
            'The challenge package changed source identity, classification, or readiness.')
    require(sha256(canonical(read_scientific(root/'scientific-policy.json')))==metadata.bindings['policy'],
            'The retained original scientific policy differs.')
    _check_parent(root,manifest.get('revision_parent'),metadata)
    legacy=read_scientific(root/'source/presentation.json')
    require(legacy.get('schema') in {'genosapiens.track2-report-presentation/v1','genosapiens.track2-report-presentation/v2'},
            'The challenge package has an unsupported source presentation transition.')
    require(json_identity(read_scientific(root/'disclosures.json'))==json_identity(legacy['disclosure']),
            'The revision changed source disclosures.')
    historical={k:reference(root,'source/supplements/'+k+'-ranking.json') for k in ('combined','lincs','mechanism')}
    require(read_scientific(root/'historical-rankings.json')==historical, 'The revision changed historical ranking references.')
    # Establish original-condition aggregates before checking their retained
    # calibration copies. Exact scientific fingerprints connect these proofs.
    context.derived_report(root/'source',legacy)
    calibration=None
    if manifest['calibration'] is not None:
        from genosapiens.track2_report_calibration import check_calibration
        require(manifest['calibration']['path']=='calibration/manifest.json', 'The calibration component path differs.')
        component=context.checked_calibration(root/'calibration',expected_manifest_sha256=manifest['calibration']['sha256'])
        require(component['classification']==metadata.classification, 'The calibration component classification differs.')
        calibration=dict(bank=root/'calibration',resource_contract=root/'calibration/resource-contract.json',
            expected_resource_sha256=component['resource_contract_sha256'])
    derived=_checked_derived(root/'source',legacy,calibration,context)
    if derived is not None:
        require(json_identity(read_scientific(root/'derived-rankings.json'))==json_identity(derived),
                'The revision derived scores differ from retained scientific evidence.')
    model=check_model(read_scientific(root/'presentation.json'))
    expected=_model(root,legacy,metadata,derived,manifest['calibration'],
        narrative_version=model.get('narrative_version',1))
    from genosapiens.track2_challenge_overviews import figure_references
    require(manifest['figures']==reference(root,'figures/manifest.json'), 'The figure component binding differs.')
    expected['figures']=figure_references(root,expected,root/'figures')
    require(json_identity(model)==json_identity(expected), 'The challenge presentation differs from retained evidence.')
    require(manifest['presentation']==reference(root,'presentation.json')
            and manifest['scientific']==reference(root,'source/manifest.json'), 'A challenge component binding differs.')
    require(read_scientific(root/'readiness.json')==metadata.readiness, 'The revision changed scientific readiness.')
    require((root/'report.md').read_text()==''.join(paper_markdown(model,'# Challenge research report',metadata.readiness,
        hypothesis_context=legacy['context'])),
            'Challenge prose differs from its presentation contract.')
    from genosapiens.track2_report_links import check_links
    check_links(root,(root/'report.md').read_text())
    return manifest


def check_revision(root, *, expected_manifest_sha256, context=None):
    from genosapiens.track2_checked_evidence import CheckedEvidenceContext
    if context is None:
        with CheckedEvidenceContext() as owned:
            return check_revision(root,expected_manifest_sha256=expected_manifest_sha256,context=owned)
    require(isinstance(context,CheckedEvidenceContext), 'Revision admission requires an invocation-owned context.')
    root=Path(root)
    require(file_identity(root/'manifest.json')['sha256']==expected_manifest_sha256,
            'The challenge package differs from its independently expected identity.')
    manifest=read_scientific(root/'manifest.json')
    require(manifest.get('schema')==SCHEMA, 'The challenge package version is unsupported.')
    metadata=context.checked_report(root/'source')
    if (root/'producer-admission.json').is_file():
        from genosapiens.track2_report_producer import check_retained_producer
        check_retained_producer(read_scientific(root/'producer-admission.json'),metadata.manifest_sha256)
    return _check_contents(root,manifest,metadata=metadata,context=context)


def _check_parent(root, parent, metadata):
    root=Path(root)
    require((parent is not None)==(root/'revision-parent.json').is_file(),'The source revision identity is missing or unexpected.')
    if parent is None:return
    require(parent==reference(root,'revision-parent.json'),'The source revision identity differs.')
    value=read_scientific(root/'revision-parent.json')
    require(value.get('schema')==SCHEMA and value['source_kind']==metadata.kind
        and value['source_manifest_sha256']==metadata.manifest_sha256
        and value['source_bindings']==metadata.bindings and value['classification']==metadata.classification
        and value['readiness']==metadata.readiness,'The source revision has incompatible scientific evidence.')


def _admit_revision_source(root, context):
    """Admit retained science without requiring the predecessor's current layout.

    The caller supplies the independently expected parent manifest identity.
    Old figure and presentation bytes remain authenticated, but new versions
    need not reproduce their old layout to reuse compatible scientific records.
    """
    root=Path(root);manifest=read_scientific(root/'manifest.json');before=inventory(root)
    require(manifest['outputs']=={k:v for k,v in before.items() if k!='manifest.json'},
        'The source revision inventory differs.')
    for key,path in (('scientific','source/manifest.json'),('presentation','presentation.json'),('figures','figures/manifest.json')):
        require(manifest[key]==reference(root,path),'A source revision component differs.')
    metadata=context.checked_report(root/'source')
    if (root/'producer-admission.json').is_file():
        from genosapiens.track2_report_producer import check_retained_producer
        check_retained_producer(read_scientific(root/'producer-admission.json'),metadata.manifest_sha256)
    require(manifest['source_manifest_sha256']==metadata.manifest_sha256
        and manifest['source_bindings']==metadata.bindings and manifest['source_kind']==metadata.kind
        and manifest['classification']==metadata.classification and manifest['readiness']==metadata.readiness
        and manifest['input_snapshot_sha256']==metadata.input_snapshot_sha256,
        'The source revision changed scientific identity or readiness.')
    _check_parent(root,manifest.get('revision_parent'),metadata)
    legacy=read_scientific(root/'source/presentation.json');context.derived_report(root/'source',legacy)
    calibration=None
    if manifest['calibration'] is not None:
        require(manifest['calibration']==reference(root,'calibration/manifest.json'),'The source calibration binding differs.')
        component=context.checked_calibration(root/'calibration',expected_manifest_sha256=manifest['calibration']['sha256'])
        require(component['classification']==metadata.classification,'The source calibration classification differs.')
        calibration=dict(bank=root/'calibration',resource_contract=root/'calibration/resource-contract.json',expected_resource_sha256=component['resource_contract_sha256'])
    derived=_checked_derived(root/'source',legacy,calibration,context)
    model=read_scientific(root/'presentation.json')
    require(model.get('schema')=='genosapiens.track2-report-presentation/v3','The source presentation transition is unsupported.')
    require(model['scientific_policy']==reference(root,'scientific-policy.json')
        and model['scientific_policy']['sha256']==metadata.bindings['policy'],
        'The source scientific policy lacks a compatible transition.')
    from genosapiens.track2_derived_strength import policy
    require(model['derived_policy']==(policy() if derived is not None else None),'The source derived policy lacks a compatible transition.')
    if derived is not None:
        require(json_identity(read_scientific(root/'derived-rankings.json'))==json_identity(derived),
            'The source derived scores differ from retained scientific evidence.')
    require(inventory(root)==before,'The source revision changed during scientific admission.')
    return dict(source_manifest_sha256=metadata.manifest_sha256,calibration=calibration)
