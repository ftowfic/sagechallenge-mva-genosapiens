"""Challenge packages composed directly from independently admitted evidence."""
from contextlib import ExitStack
from pathlib import Path
from genosapiens.track2_report_local_model import _LocalModel
import shutil

from genosapiens.track2_operator_inputs import require
from genosapiens.track2_scientific_reader import read_scientific
from genosapiens.track2_scientific_storage import file_identity,json_identity,publish_json,publish_chunks
from genosapiens.track2_challenge_revision import inventory,reference,_model
from genosapiens.track2_report_contracts import PACKAGE_SCHEMA,checked_contract_file,check_package_contract,_complete_package_contract
from genosapiens.track2_report_evidence import check_evidence,scientific_root,scientific_model


def _versions(root, evidence):
    versions=read_scientific(Path(root)/evidence['records']['source_versions']['path'])['calibration']
    return dict(expected_fda_snapshot_sha256=versions['fda_snapshot_sha256'],
        expected_source_release=versions['source_release'])


def derive_evidence(root, evidence, calibration):
    from genosapiens.track2_report_derived import derive_rankings,report_conditions,calibrate_rankings
    root=Path(root);source=scientific_root(root,evidence);model=scientific_model(root,evidence)
    if model['schema']!='genosapiens.track2-report-presentation/v2':return None
    derived=derive_rankings(conditions=report_conditions(source),
        queries=[read_scientific(p) for p in sorted((source/'lincs').glob('*/*/science/query.json'))],
        evaluated_scores=read_scientific(source/'discovery/tissue-scores.json'),
        objectives=(o for h in model['context']['hypotheses'] for o in h['objectives']),
        historical_lincs=model['selections']['lincs']['ranking'])
    if calibration is not None:
        derived=calibrate_rankings(derived,**calibration,**_versions(root,evidence))
    for row in derived['complete_export']:
        require(row['status']!='complete' or row['direction']!='increase' or row['numerator'] is not None,
                'An applicable derived score lacks completed calibration.')
    return derived


def evidence_model(root, evidence, derived, calibration_reference, *, _scientific=None):
    """Borrow the caller's model within its owned helper chain.

    Without an owner, load records locally and return an independent disk view.
    Borrowed views remain open until their caller finishes or fails.
    """
    if _scientific is None:
        from genosapiens.track2_scientific_records import stored_json
        with _LocalModel(Path(root)/'evidence',evidence) as owned:
            return stored_json(evidence_model(root,evidence,derived,calibration_reference,_scientific=owned))
    from genosapiens.track2_report_metadata import ReportMetadata
    root=Path(root);source=root/'evidence';legacy=_scientific.borrow(source,evidence)
    metadata=ReportMetadata(evidence['source_kind'],evidence['bindings'],evidence['input_snapshot_sha256'],
        evidence['classification'],evidence['readiness'],reference(root,'evidence/manifest.json')['sha256'])
    prefix=scientific_root(source,evidence).relative_to(root).as_posix()
    model=_model(root,legacy,metadata,derived,calibration_reference,source_prefix=prefix,
        source_versions=_versions(source,evidence),scientific_manifest='evidence/manifest.json')
    model['supplements'].append(dict(title='Calibration preparation',artifact=reference(root,'calibration-preparation.json')))
    model['supplements'].append(dict(title='Calibration applicability',artifact=reference(root,'calibration-applicability.json')))
    _scientific.borrow(source,evidence)
    return model


_UNPREPARED=object()


def check_producer_lineage(evidence, manifest, proof):
    """A later revision must retain the proof already checked in its parent."""
    from genosapiens.track2_report_producer import check_retained_producer
    original=manifest['historical_source']
    require(original is not None, 'Producer transition proof requires historical evidence.')
    if proof['source_manifest_sha256']==original['manifest']['sha256']:
        check_retained_producer(proof,original['manifest']['sha256'])
        return
    parent=Path(evidence)/'historical'
    value=read_scientific(parent/'manifest.json')
    require(value.get('schema') in {'genosapiens.track2-challenge-package/v1',PACKAGE_SCHEMA}
        and 'producer-admission.json' in value['outputs']
        and read_scientific(parent/'producer-admission.json')==proof,
        'The revision changed its original producer evidence.')
    # check_evidence has independently admitted the entire parent package.
    expected=value['source_manifest_sha256'] if value['schema']=='genosapiens.track2-challenge-package/v1' else proof['source_manifest_sha256']
    check_retained_producer(proof,expected)


def compose_content(evidence, destination, *, manifest, calibration, applicability, producer_admission=None,
                    derived_result=_UNPREPARED, _checked_result=None):
    """Create content with one owned scientific model for adjacent consumers.

    The owner closes source views on success or failure. The result owns a
    separate disk view of composition.json, so it survives source cleanup.
    """
    from genosapiens.track2_synthesis import policies
    from genosapiens.track2_published_genotype import canonical,sha256
    from genosapiens.track2_report_calibration import publish_calibration
    evidence=Path(evidence).absolute();destination=Path(destination).absolute()
    require(destination.resolve()==destination and evidence.resolve()==evidence and not destination.exists()
        and not destination.is_relative_to(evidence) and not evidence.is_relative_to(destination),
        'Choose an independent new report content directory.')
    check_evidence(evidence,manifest)
    before=inventory(evidence)
    policy=policies()[0]
    require(sha256(canonical(policy))==manifest['bindings']['policy'],
            'The source scientific policy lacks a compatible transition.')
    destination.mkdir(mode=0o700,parents=True)
    shutil.copytree(evidence,destination/'evidence')
    for path in (destination/'evidence').rglob('*'):path.chmod(0o700 if path.is_dir() else 0o600)
    require(inventory(destination/'evidence')==before and inventory(evidence)==before,
            'Evidence changed during report copying.')
    (destination/'scientific-policy.json').write_bytes(canonical(policy))
    source=scientific_root(destination/'evidence',manifest)
    prefix=source.relative_to(destination).as_posix()
    with ExitStack() as stack:
        owned=stack.enter_context(_LocalModel(destination/'evidence',manifest))
        legacy=owned.borrow(destination/'evidence',manifest)
        publish_json(destination/'disclosures.json',legacy['disclosure'])
        publish_json(destination/'historical-rankings.json',{k:reference(destination,prefix+'/supplements/'+k+'-ranking.json')
            for k in ('combined','lincs','mechanism')})
        if producer_admission is not None:
            check_producer_lineage(evidence,manifest,producer_admission)
            publish_json(destination/'producer-admission.json',producer_admission)
        calibration_args={k:calibration[k] for k in ('bank','resource_contract','expected_resource_sha256')}
        if derived_result is _UNPREPARED:
            from genosapiens.track2_report_local_derivation import _checked_derivation
            derived,_checked_result=stack.enter_context(_checked_derivation(destination/'evidence',manifest,calibration_args))
        else:derived=derived_result
        calibration_ref=None
        if derived is not None:
            publish_json(destination/'derived-rankings.json',derived)
            from genosapiens.track2_checked_evidence import CheckedEvidenceContext
            from genosapiens.track2_report_local_derivation import _CheckedDerivation
            calibration_context=stack.enter_context(CheckedEvidenceContext())
            if _checked_result is not None:
                require(isinstance(_checked_result,_CheckedDerivation), 'Calibration requires a current checked derivation handoff.')
                _checked_result.seed(calibration_context,evidence=manifest,calibration=calibration_args,derived=derived)
            calibration_ref=publish_calibration(destination/'calibration',**calibration_args,context=calibration_context,
                observations=(r for r in derived['complete_export'] if r.get('numerator') is not None),
                classification=manifest['classification'])
        publish_json(destination/'calibration-applicability.json',applicability)
        publish_json(destination/'calibration-preparation.json',calibration['operations'])
        publish_json(destination/'readiness.json',manifest['readiness'])
        publish_json(destination/'report-review.json',dict(schema='genosapiens.track2-report-review/v1',
            report_human_review='not_performed',original_reviews=manifest['records']['reviews'],
            original_reviews_root='evidence',scientific_decisions='unchanged'))
        model=evidence_model(destination,manifest,derived,calibration_ref,_scientific=owned)
        publish_json(destination/'composition.json',model)
        owned.borrow(destination/'evidence',manifest)
    return read_scientific(destination/'composition.json')


def render_content(root, *, figure_renderer=None):
    from genosapiens.track2_challenge_overviews import render_overviews,figure_references,pinned_runtime
    from genosapiens.track2_report_model import publish_model
    from genosapiens.track2_report_paper import paper_markdown
    root=Path(root);evidence=read_scientific(root/'evidence/manifest.json')
    model=dict(read_scientific(root/'composition.json'))
    if not (root/'figures').exists():
        render_overviews(root,model,root/'figures',runtime_image=pinned_runtime(),run_renderer=figure_renderer)
    model['figures']=figure_references(root,model,root/'figures')
    publish_model(root/'presentation.json',model)
    legacy=scientific_model(root/'evidence',evidence)
    publish_chunks(root/'report.md',(s.encode() for s in paper_markdown(model,'# Challenge research report',
        evidence['readiness'],hypothesis_context=legacy['context'],integrated_workflow=True,
        calibration_preparation=read_scientific(root/'calibration-preparation.json'))))
    return model


def package_manifest(root):
    """Construct and check an envelope with one immediate inventory scan."""
    root=Path(root);evidence=read_scientific(root/'evidence/manifest.json')
    value=dict(schema=PACKAGE_SCHEMA,classification=evidence['classification'],source_kind=evidence['source_kind'],
        input_origin=evidence['input_origin'],source_bindings=evidence['bindings'],
        input_snapshot_sha256=evidence['input_snapshot_sha256'],readiness=evidence['readiness'],
        evidence=reference(root,'evidence/manifest.json'),scientific_policy=reference(root,'scientific-policy.json'),
        presentation=reference(root,'presentation.json'),figures=reference(root,'figures/manifest.json'),
        calibration=reference(root,'calibration/manifest.json') if (root/'calibration').exists() else None,
        calibration_applicability=reference(root,'calibration-applicability.json'),
        derived_rankings=reference(root,'derived-rankings.json') if (root/'derived-rankings.json').exists() else None,
        producer_admission=reference(root,'producer-admission.json') if (root/'producer-admission.json').exists() else None,
        report_review=reference(root,'report-review.json'))
    return _complete_package_contract(value,root=root)


def publish_package(evidence,destination,*,manifest,calibration,applicability,figure_renderer=None,producer_admission=None):
    compose_content(evidence,destination,manifest=manifest,calibration=calibration,
        applicability=applicability,producer_admission=producer_admission)
    render_content(destination,figure_renderer=figure_renderer)
    value=package_manifest(destination)
    check_contents(destination,value)
    publish_json(Path(destination)/'manifest.json',value)
    return reference(destination,'manifest.json')


def check_contents(root,value, *, _admitted=None):
    """Check current files and borrow one model until semantic checks finish.

    Independent invocations open fresh documents. Closing integrity checks
    still run after semantic processing and before publication.
    """
    from genosapiens.track2_report_calibration_inputs import extract_calibration_inputs
    from genosapiens.track2_report_calibration import check_calibration
    from genosapiens.track2_challenge_overviews import figure_references
    from genosapiens.track2_report_paper import paper_markdown
    from genosapiens.track2_report_links import check_links
    from genosapiens.track2_synthesis import policies
    from genosapiens.track2_published_genotype import canonical,sha256
    root=Path(root);check_package_contract(value,root=root)
    require(value['evidence']==reference(root,'evidence/manifest.json'), 'The report evidence binding differs.')
    evidence=checked_contract_file(root/'evidence',expected_manifest_sha256=value['evidence']['sha256'])
    extracted=extract_calibration_inputs(root/'evidence',evidence) if _admitted is None else _admitted[0]
    for key,source_key in (('classification','classification'),('source_kind','source_kind'),('input_origin','input_origin'),
                          ('source_bindings','bindings'),('input_snapshot_sha256','input_snapshot_sha256'),('readiness','readiness')):
        require(value[key]==evidence[source_key], 'The report changed original scientific identity or readiness.')
    policy=read_scientific(root/'scientific-policy.json')
    require(sha256(canonical(policy))==evidence['bindings']['policy'] and policy==policies()[0],
            'The report changed its original scientific policy.')
    require(read_scientific(root/'calibration-applicability.json')==extracted['applicability'],
            'The calibration applicability differs from retained evidence.')
    preparation=read_scientific(root/'calibration-preparation.json')
    counts=('entries_new','entries_reused','entries_imported','entry_count')
    require(preparation['status']=='complete' and all(type(preparation[k]) is int and preparation[k]>=0 for k in counts)
        and preparation['entries_new']+preparation['entries_reused']==preparation['entry_count']
        and preparation['entries_imported']<=preparation['entries_reused']
        and all(type(preparation['operation_counts'][k]) is int and preparation['operation_counts'][k]>=0
                for k in ('query_permutations','aggregate_iterations')),
        'The retained calibration preparation record is incomplete or inconsistent.')
    calibration=None
    if value['calibration'] is not None:
        component=check_calibration(root/'calibration',expected_manifest_sha256=value['calibration']['sha256'])
        require(component['classification']==evidence['classification'], 'Calibration classification differs.')
        calibration=dict(bank=root/'calibration',resource_contract=root/'calibration/resource-contract.json',
            expected_resource_sha256=component['resource_contract_sha256'])
    derived=derive_evidence(root/'evidence',evidence,calibration) if _admitted is None else _admitted[1]
    require((derived is not None)==(value['derived_rankings'] is not None), 'The derived ranking component is missing or unexpected.')
    if derived is not None:
        require(json_identity(derived)==json_identity(read_scientific(root/'derived-rankings.json')),
                'The report derived scores differ from retained science.')
    with _LocalModel(root/'evidence',evidence) as owned:
        legacy=owned.borrow(root/'evidence',evidence)
        require(json_identity(read_scientific(root/'disclosures.json'))==json_identity(legacy['disclosure']),
                'The report changed source disclosures.')
        prefix=scientific_root(root/'evidence',evidence).relative_to(root).as_posix()
        require(read_scientific(root/'historical-rankings.json')=={k:reference(root,prefix+'/supplements/'+k+'-ranking.json')
            for k in ('combined','lincs','mechanism')}, 'The report changed historical ranking references.')
        model=evidence_model(root,evidence,derived,value['calibration'],_scientific=owned)
        require(json_identity(model)==json_identity(read_scientific(root/'composition.json')), 'The report composition differs.')
        model['figures']=figure_references(root,model,root/'figures')
        require(json_identity(model)==json_identity(read_scientific(root/'presentation.json')), 'The report presentation differs.')
        require(read_scientific(root/'readiness.json')==evidence['readiness'], 'The report readiness differs.')
        require(read_scientific(root/'report-review.json')==dict(schema='genosapiens.track2-report-review/v1',
            report_human_review='not_performed',original_reviews=evidence['records']['reviews'],
            original_reviews_root='evidence',scientific_decisions='unchanged'), 'Report review attribution differs.')
        if value['producer_admission'] is not None:
            check_producer_lineage(root/'evidence',evidence,read_scientific(root/'producer-admission.json'))
        require((root/'report.md').read_text()==''.join(paper_markdown(model,'# Challenge research report',
            evidence['readiness'],hypothesis_context=legacy['context'],integrated_workflow=True,
            calibration_preparation=preparation)), 'The report prose differs from its evidence.')
        check_links(root,(root/'report.md').read_text())
        require(value==package_manifest(root), 'A report component binding differs.')
        owned.borrow(root/'evidence',evidence)
    return value


def check_package(root, *, expected_manifest_sha256):
    value=checked_contract_file(root,expected_manifest_sha256=expected_manifest_sha256,package=True)
    check_contents(root,value)
    require(file_identity(Path(root)/'manifest.json')['sha256']==expected_manifest_sha256,
            'The report manifest changed during admission.')
    return value
