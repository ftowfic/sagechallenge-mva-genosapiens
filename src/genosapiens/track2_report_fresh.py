"""The shared fresh-run boundary from scientific exports to the final report."""
from pathlib import Path

from genosapiens.track2_scientific_storage import file_identity
from genosapiens.track2_scientific_reader import read_scientific


def run_fresh_report(output, *, evidence_binding, export, require_review=False, controls=None, figure_renderer=None):
    from genosapiens.track2_report_stages import stage,code_binding,private_directory
    from genosapiens.track2_report_workflow import report_workflow
    output=Path(output).absolute();parent=output.parent
    private_directory(parent)
    evidence_operation=parent/(output.name+'.evidence-operation')
    evidence,_=stage(evidence_operation,'export',dict(science=evidence_binding,
        code=code_binding('track2_discovery_report','track2_synthesis','track2_report_evidence',
            'track2_report_sources','track2_report_paper')),export)
    operation=parent/(output.name+'.operation')
    completed=report_workflow(evidence,output,expected_evidence_sha256=file_identity(evidence/'manifest.json')['sha256'],
        operation=operation,resume=operation.exists(),controls=controls,require_review=require_review,figure_renderer=figure_renderer)
    result=dict(completed['readiness'],manifest_sha256=completed['report']['sha256'])
    if read_scientific(evidence/'manifest.json')['source_kind']=='study':result['output']=str(output)
    return result


def coordinator_metadata(output, *, expected_manifest_sha256):
    """Use the independently retained finalization output expectation internally.

    The caller must first admit original production dependencies and saved result.
    External package admission continues to replay complete scientific semantics.
    """
    from genosapiens.track2_operator_inputs import require
    from genosapiens.track2_report_contracts import checked_contract_file
    from genosapiens.track2_report_metadata import ReportMetadata
    from genosapiens.track2_challenge_revision import inventory
    from genosapiens.track2_report_stages import code_binding
    output=Path(output)
    record=read_scientific(output.with_name(output.name+'.operation')/'stages/finalization.json')
    require(record['status']=='complete' and record['stage']=='finalization'
        and record['outputs']==inventory(output)
        and record['binding']['code']==code_binding('track2_report_package','track2_report_paper',
            'track2_challenge_paper','track2_report_links','track2_report_contracts'),
        'The report finalization expectation or its dependencies changed.')
    value=checked_contract_file(output,expected_manifest_sha256=expected_manifest_sha256,package=True)
    return ReportMetadata(value['source_kind'],dict(value['source_bindings']),value['input_snapshot_sha256'],
        value['classification'],value['readiness'],expected_manifest_sha256)
