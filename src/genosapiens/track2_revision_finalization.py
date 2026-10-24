"""Finalize separate challenge revisions without changing original run state."""
from pathlib import Path
import uuid

from genosapiens.track2_operator_inputs import require
from genosapiens.track2_run_state import atomic_record
from genosapiens.track2_scientific_reader import read_scientific
from genosapiens.track2_scientific_storage import file_identity


def operation_directory(report, operation=None, *, protected=()):
    report=Path(report).absolute()
    operation=Path(operation).absolute() if operation is not None else report.with_name(report.name+'-operation')
    require(operation.resolve()==operation and not operation.is_relative_to(report)
        and not report.is_relative_to(operation),'Revision operations require an independent directory without symbolic links.')
    for source in protected:
        source=Path(source).absolute()
        require(not operation.is_relative_to(source) and not source.is_relative_to(operation),
            'Revision operations cannot write into the source report.')
    marker=dict(schema='genosapiens.track2-revision-operation/v1',report_name=report.name)
    if operation.exists():
        require((operation/'operation.json').is_file()
            and read_scientific(operation/'operation.json')==marker,
            'Choose a new operation directory or the matching revision operation directory.')
    else:
        operation.mkdir(mode=0o700,parents=True,exist_ok=False)
        atomic_record(operation/'operation.json',marker)
    return operation


def finalize_revision(report, *, expected_manifest_sha256, context=None, operation=None, require_review=False):
    from genosapiens.track2_checked_evidence import CheckedEvidenceContext
    if context is None:
        with CheckedEvidenceContext() as owned:
            return finalize_revision(report,expected_manifest_sha256=expected_manifest_sha256,
                context=owned,operation=operation,require_review=require_review)
    require(isinstance(context,CheckedEvidenceContext),'Revision finalization requires an invocation-owned context.')
    report=Path(report);operation=operation_directory(report,operation)
    attempt=operation/('finalization-'+uuid.uuid4().hex+'.json')
    initial=dict(schema='genosapiens.track2-revision-finalization/v1',state='incomplete',
        expected_manifest_sha256=expected_manifest_sha256)
    atomic_record(attempt,initial);atomic_record(operation/'current.json',initial)
    try:
        from genosapiens.track2_report_metadata import checked_report_metadata
        from genosapiens.track2_challenge_revision import SCHEMA
        require(read_scientific(report/'manifest.json').get('schema') in {SCHEMA,'genosapiens.track2-challenge-package/v2'},
            'Separate revision finalization requires the supported challenge package contract.')
        metadata=checked_report_metadata(report,expected_manifest_sha256=expected_manifest_sha256,context=context)
        from genosapiens.track2_report_assembly import presentation_ready
        require(presentation_ready(report),'The challenge presentation is incomplete.')
        from genosapiens.track2_terminal import assessment_outcome
        outcome=assessment_outcome(metadata.readiness,require_review=require_review)
        record=dict(initial,state='complete',classification=metadata.classification,
            source_kind=metadata.kind,scientific_readiness=metadata.readiness,outcome=outcome.as_dict(),
            validator_sha256=context._validators)
        atomic_record(attempt,record);atomic_record(operation/'current.json',record)
        return outcome
    except BaseException as error:
        from genosapiens.track2_discovery_outcome import safe_diagnostic
        record=dict(initial,state='interrupted' if isinstance(error,(KeyboardInterrupt,TimeoutError)) else 'failed',
            diagnostic=safe_diagnostic('revision-finalization','report-completion',error))
        atomic_record(attempt,record);atomic_record(operation/'current.json',record)
        raise


def publish_and_finalize_revision(source, destination, *, expected_source_manifest_sha256,
        calibration=None, figure_renderer=None, require_review=False, operation=None, on_publish=None, production=None):
    """Share one admitted context across assembly and finalization.

    The operational journal is not an authenticity anchor. Standalone recovery
    still requires the independently retained publication digest.
    """
    from genosapiens.track2_checked_evidence import CheckedEvidenceContext
    from genosapiens.track2_challenge_revision import publish_revision
    source=Path(source).absolute();destination=Path(destination).absolute()
    require(source.resolve()==source and destination.resolve()==destination and not destination.exists()
        and not destination.is_relative_to(source) and not source.is_relative_to(destination),
        'Choose a new independent challenge revision directory.')
    from genosapiens.track2_challenge_revision import SCHEMA
    from genosapiens.track2_report_producer import admit_producer,check_retained_producer
    source_manifest=read_scientific(source/'manifest.json')
    producer_admission=None
    if source_manifest.get('schema')==SCHEMA:
        require(production is None,'A later revision must retain its original producer evidence.')
        require((source/'producer-admission.json').is_file(),'This source revision lacks explicit producer evidence.')
        check_retained_producer(read_scientific(source/'producer-admission.json'),source_manifest['source_manifest_sha256'])
    else:
        require(production is not None,'A first revision requires an independently bound original production record.')
        producer_admission=admit_producer(production['record'],
            expected_production_sha256=production['expected_sha256'],
            expected_manifest_sha256=expected_source_manifest_sha256,
            native=source_manifest.get('schema')=='genosapiens.track2-scientific-completion/v1')
    operation=operation_directory(destination,operation,protected=(source,))
    journal=operation/('assembly-'+uuid.uuid4().hex+'.json')
    record=dict(schema='genosapiens.track2-revision-assembly/v1',state='incomplete',
        expected_source_manifest_sha256=expected_source_manifest_sha256,
        expected_calibration_resource_sha256=calibration['expected_resource_sha256'] if calibration is not None else None)
    atomic_record(journal,record);atomic_record(operation/'current.json',record)
    with CheckedEvidenceContext() as context:
        try:
            binding=publish_revision(source,destination,expected_source_manifest_sha256=expected_source_manifest_sha256,
                calibration=calibration,context=context,figure_renderer=figure_renderer,producer_admission=producer_admission)
            # Emit the publication identity before finalization so a caller can
            # retain an independent recovery anchor if the final stage fails.
            if on_publish is not None:on_publish(binding)
            record=dict(record,state='assembled',report_binding=binding)
            atomic_record(journal,record);atomic_record(operation/'current.json',record)
            outcome=finalize_revision(destination,expected_manifest_sha256=binding['sha256'],context=context,
                operation=operation,require_review=require_review)
            atomic_record(journal,dict(record,state='complete',outcome=outcome.as_dict()))
            return dict(report_binding=binding,outcome=outcome.as_dict(),semantic_passes=dict(context.counters))
        except BaseException as error:
            from genosapiens.track2_discovery_outcome import safe_diagnostic
            failure=dict(record,state='interrupted' if isinstance(error,(KeyboardInterrupt,TimeoutError)) else 'failed',
                diagnostic=safe_diagnostic('revision-assembly','report-completion',error))
            atomic_record(journal,failure);atomic_record(operation/'current.json',failure)
            raise
