"""Shared reporting sequence over independently admitted retained evidence."""
from contextlib import ExitStack
import time
import uuid
import shutil
from pathlib import Path

from genosapiens.track2_operator_inputs import require
from genosapiens.track2_scientific_storage import file_identity,publish_json
from genosapiens.track2_scientific_reader import read_scientific


def revise_report(source, destination, *, expected_source_manifest_sha256, production=None,
                  operation=None, resume=False, figure_renderer=None, require_review=False, controls=None):
    """Share the automatic workflow with the retained-report command."""
    from genosapiens.track2_report_metadata import checked_report_metadata
    from genosapiens.track2_report_producer import admit_producer
    from genosapiens.track2_report_evidence import adapt_report_evidence
    from genosapiens.track2_report_stages import destination_lock,stage,code_binding,atomic_publish
    source=Path(source).absolute();destination=Path(destination).absolute()
    require(not source.is_relative_to(destination) and not destination.is_relative_to(source),
            'Choose an independent report destination.')
    metadata=checked_report_metadata(source,expected_manifest_sha256=expected_source_manifest_sha256)
    value=read_scientific(source/'manifest.json')
    if value.get('schema') in {'genosapiens.track2-challenge-package/v1','genosapiens.track2-challenge-package/v2'}:
        require(production is None, 'A later revision must retain its original producer evidence.')
        proof=read_scientific(source/'producer-admission.json')
    else:
        require(production is not None, 'A first revision requires an independently bound original production record.')
        proof=admit_producer(production['record'],expected_production_sha256=production['expected_sha256'],
            expected_manifest_sha256=expected_source_manifest_sha256,native=metadata.kind=='native')
    evidence=destination.with_name(destination.name+'.evidence')
    evidence_operation=destination.with_name(destination.name+'.evidence-operation')
    require(not evidence_operation.is_relative_to(source) and not source.is_relative_to(evidence_operation)
        and not evidence.is_relative_to(source) and not source.is_relative_to(evidence),
        'Report evidence work must remain outside the source.')
    identity=dict(source=str(source),manifest_sha256=expected_source_manifest_sha256)
    with destination_lock(evidence,evidence_operation,identity=identity,resume=resume):
        if not evidence.exists():
            copied,_=stage(evidence_operation,'adapt',dict(identity,code=code_binding('track2_report_evidence')),
                lambda output:adapt_report_evidence(source,output,expected_manifest_sha256=expected_source_manifest_sha256))
            atomic_publish(copied,evidence)
        from genosapiens.track2_challenge_revision import inventory
        require(inventory(evidence)==read_scientific(evidence_operation/'stages/adapt.json')['outputs'],
                'Adapted report evidence changed.')
    return report_workflow(evidence,destination,expected_evidence_sha256=file_identity(evidence/'manifest.json')['sha256'],
        operation=operation,resume=resume,figure_renderer=figure_renderer,producer_admission=proof,
        require_review=require_review,controls=controls)


def prepare_and_compose(evidence, destination, operation, *, expected_evidence_sha256, compose,
                        prepare=None):
    """Applicable calibration must finish before composition can start."""
    from genosapiens.track2_report_calibration_inputs import extract_calibration_inputs
    from genosapiens.track2_report_worker import prepare_bounded_calibration
    evidence=Path(evidence);operation=Path(operation);destination=Path(destination)
    require(file_identity(evidence/'manifest.json')['sha256']==expected_evidence_sha256,
            'The report evidence differs from its independently admitted identity.')
    manifest=read_scientific(evidence/'manifest.json')
    operation.mkdir(mode=0o700,parents=True,exist_ok=True)
    attempts=operation/'workflow-attempts';attempts.mkdir(mode=0o700,exist_ok=True)
    attempt=attempts/uuid.uuid4().hex;attempt.mkdir(mode=0o700)
    publish_json(attempt/'started.json',dict(status='incomplete',classification=manifest['classification'],
        source_manifest_sha256=expected_evidence_sha256))
    start=time.monotonic();stage='evidence_admission'
    try:
        extracted=extract_calibration_inputs(evidence,manifest)
        require(file_identity(evidence/'manifest.json')['sha256']==expected_evidence_sha256,
                'The report evidence changed during admission.')
        stage='calibration'
        prepare=prepare_bounded_calibration if prepare is None else prepare
        calibration=prepare(operation/'calibration',extracted['inputs'],classification=manifest['classification'],
            source_manifest_sha256=expected_evidence_sha256)
        require(calibration['operations']['status']=='complete', 'Calibration preparation did not complete.')
        if extracted['applicability']['applicable_count']:
            require(calibration['operations']['entry_count']>0, 'Applicable evidence has no completed calibration entries.')
        stage='composition'
        report=compose(evidence,destination,manifest=manifest,calibration=calibration,
            applicability=extracted['applicability'])
        result=dict(status='composed',report=report,calibration_operations=calibration['operations'],
            scientific_outcome=manifest['readiness']['assessment_status'])
    except BaseException as error:
        publish_json(attempt/'outcome.json',dict(status='incomplete',stage=stage,error_type=type(error).__name__,
            elapsed_seconds=time.monotonic()-start))
        raise
    publish_json(attempt/'outcome.json',dict(status='complete',stage=stage,elapsed_seconds=time.monotonic()-start))
    return result


def report_workflow(evidence, destination, *, expected_evidence_sha256, operation=None, resume=False,
                    controls=None, cancelled=None, figure_renderer=None, producer_admission=None, on_stage=None,
                    require_review=False):
    """Run and recover report stages without changing the original evidence."""
    from genosapiens.track2_report_stages import destination_lock,stage,code_binding,atomic_publish
    from genosapiens.track2_report_contracts import checked_contract_file
    from genosapiens.track2_report_evidence import check_evidence
    from genosapiens.track2_report_calibration_inputs import extract_calibration_inputs
    from genosapiens.track2_report_worker import prepare_bounded_calibration,bounded_renderer,limits
    from genosapiens.track2_report_package import derive_evidence,compose_content,render_content,package_manifest,check_contents
    from genosapiens.track2_challenge_revision import inventory,reference
    from genosapiens.track2_challenge_overviews import render_overviews,pinned_runtime,check_overviews
    from genosapiens.track2_scientific_storage import json_identity
    evidence=Path(evidence).absolute();destination=Path(destination).absolute()
    operation=Path(operation).absolute() if operation is not None else destination.with_name(destination.name+'.operation')
    require(not destination.is_relative_to(evidence) and not evidence.is_relative_to(destination)
        and not operation.is_relative_to(evidence) and not evidence.is_relative_to(operation),
        'Report work must be independent of its evidence source.')
    controls=limits(**(controls or {}));started=time.monotonic()
    def remaining():
        require(cancelled is None or not cancelled(), 'The report invocation was cancelled.')
        effective=dict(controls)
        if effective['deadline_seconds'] is not None:
            effective['deadline_seconds']-=time.monotonic()-started
            require(effective['deadline_seconds']>0, 'The report invocation deadline expired.')
        return effective
    identity=dict(evidence=str(evidence),evidence_sha256=expected_evidence_sha256,
        producer_admission=producer_admission)
    reused=[]
    with ExitStack() as owned, destination_lock(destination,operation,identity=identity,resume=resume):
        manifest=checked_contract_file(evidence,expected_manifest_sha256=expected_evidence_sha256)
        check_evidence(evidence,manifest)
        attempts=operation/'workflow-attempts';attempts.mkdir(mode=0o700,exist_ok=True)
        attempt=attempts/uuid.uuid4().hex;attempt.mkdir(mode=0o700)
        publish_json(attempt/'started.json',dict(status='incomplete',requested_limits=controls,identity=identity))
        active='evidence'
        def perform(name,binding,produce):
            nonlocal active
            active=name;remaining()
            output,reuse=stage(operation,name,binding,produce,on_stage=on_stage)
            if reuse:reused.append(name)
            return output
        try:
            def extract(output):
                output.mkdir(mode=0o700)
                publish_json(output/'inputs.json',extract_calibration_inputs(evidence,manifest))
            admitted=perform('evidence',dict(source=expected_evidence_sha256,code=code_binding(
                'track2_report_calibration_inputs','track2_report_sources','track2_report_evidence')),extract)
            extracted=read_scientific(admitted/'inputs.json')
            active='calibration'
            retained=[]
            if manifest['historical_source'] is not None:
                parent=evidence/'historical';prior=read_scientific(parent/'manifest.json')
                if prior.get('schema') in {'genosapiens.track2-challenge-package/v1','genosapiens.track2-challenge-package/v2'} and prior.get('calibration'):
                    retained=[dict(root=str(parent/'calibration'),expected_manifest_sha256=prior['calibration']['sha256'])]
            calibration=prepare_bounded_calibration(operation/'calibration',extracted['inputs'],
                classification=manifest['classification'],source_manifest_sha256=expected_evidence_sha256,
                controls=remaining(),cancelled=cancelled,retained_calibration=retained)
            require(calibration['operations']['status']=='complete' and
                (not extracted['applicability']['applicable_count'] or calibration['operations']['entry_count']>0),
                'Applicable report calibration is incomplete.')
            if on_stage is not None:on_stage('calibration')
            calibration_args={k:calibration[k] for k in ('bank','resource_contract','expected_resource_sha256')}
            derived_binding=dict(evidence=expected_evidence_sha256,calibration=calibration['expected_resource_sha256'],
                code=code_binding('track2_report_package','track2_report_derived','track2_derived_strength',
                    'track2_derived_ranking','track2_calibration_bank','track2_report_local_derivation'))
            checked_result=None
            def derive(output):
                nonlocal checked_result
                from genosapiens.track2_report_local_derivation import _checked_derivation
                output.mkdir(mode=0o700)
                current,checked_result=owned.enter_context(_checked_derivation(evidence,manifest,calibration_args))
                publish_json(output/'rankings.json',current)
            derived_path=perform('derived',derived_binding,derive)
            derived=read_scientific(derived_path/'rankings.json')
            composition_binding=dict(derived=file_identity(derived_path/'rankings.json')['sha256'],
                evidence=expected_evidence_sha256,calibration=calibration['expected_resource_sha256'],
                producer_admission=producer_admission,code=code_binding('track2_report_package',
                    'track2_challenge_revision','track2_challenge_model','track2_challenge_selection','track2_challenge_summary',
                    'track2_report_local_model','track2_report_local_derivation','track2_report_calibration'))
            content=perform('composition',composition_binding,lambda output:compose_content(evidence,output,
                manifest=manifest,calibration=calibration,applicability=extracted['applicability'],
                producer_admission=producer_admission,derived_result=derived,_checked_result=checked_result))
            model=read_scientific(content/'composition.json')
            figure_binding=dict(composition=reference(content,'composition.json'),runtime=pinned_runtime(),
                code=code_binding('track2_challenge_overviews'),renderer=file_identity(
                    Path(__file__).resolve().parents[2]/'scripts/render_track2_challenge_overviews.R')['sha256'])
            renderer=figure_renderer or bounded_renderer(operation/'figure-workers',controls=remaining(),cancelled=cancelled)
            figures=perform('figures',figure_binding,lambda output:render_overviews(content,model,output,
                runtime_image=pinned_runtime(),run_renderer=renderer))
            check_overviews(content,model,figures)
            final_binding=dict(composition=json_identity(inventory(content))['sha256'],
                figures=reference(figures,'manifest.json'),code=code_binding('track2_report_package',
                    'track2_report_paper','track2_challenge_paper','track2_report_links','track2_report_contracts','track2_report_local_model'))
            active='finalization'
            if destination.exists():
                record=read_scientific(operation/'stages/finalization.json')
                require(record['binding']==final_binding and record['outputs']==inventory(destination),
                        'Published report differs from its independent finalization record.')
                reused.append('finalization')
            else:
                def finalize(output):
                    shutil.copytree(content,output)
                    shutil.copytree(figures,output/'figures')
                    render_content(output)
                    value=package_manifest(output)
                    check_contents(output,value,_admitted=(extracted,derived))
                    publish_json(output/'manifest.json',value)
                final=perform('finalization',final_binding,finalize)
                remaining();atomic_publish(final,destination)
                if on_stage is not None:on_stage('publication')
            from genosapiens.track2_terminal import assessment_outcome
            result=dict(status='complete',report=reference(destination,'manifest.json'),
                readiness=manifest['readiness'],outcome=assessment_outcome(manifest['readiness'],require_review=require_review).as_dict(),
                reused_stages=reused,calibration_operations=calibration['operations'])
        except BaseException as error:
            publish_json(attempt/'outcome.json',dict(status='incomplete',stage=active,error_type=type(error).__name__,
                elapsed_seconds=time.monotonic()-started,reused_stages=reused))
            raise
        publish_json(attempt/'outcome.json',dict(result,elapsed_seconds=time.monotonic()-started))
        return result
