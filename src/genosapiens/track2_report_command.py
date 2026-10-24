"""Offline report regeneration with source mounts owned by the host launcher."""
import sys
from pathlib import Path

from genosapiens.track2_operator_inputs import require
from genosapiens.track2_scientific_reader import read_scientific
from genosapiens.track2_scientific_storage import publish_json,file_identity


def regenerate(source, operation, *, expected_classification, reporter_sha256, resume=False,
               memory_gib=4, deadline_seconds=None, require_review=False, recover_unfinished=False):
    from genosapiens.track2_report_software import check_reporter
    from genosapiens.track2_report_run_source import admit_completed_run,check_source_metadata
    from genosapiens.track2_report_evidence import adapt_report_evidence,check_evidence
    from genosapiens.track2_report_workflow import report_workflow
    source=Path(source);operation=Path(operation)
    check_reporter(Path(__file__).resolve().parents[2],reporter_sha256)
    admitted=admit_completed_run(source,recover_unfinished=recover_unfinished)
    require(admitted['classification']==expected_classification, 'The source classification differs from launcher admission.')
    record=operation/'source-admission.json'
    if record.exists():
        require(read_scientific(record)==admitted, 'The resumed source admission changed.')
    else:publish_json(record,admitted)
    try:
        from genosapiens.track2_report_stages import stage,code_binding
        evidence,_=stage(operation,'source-evidence',dict(source=admitted['report_manifest_sha256'],
            producer=admitted['producer_admission'],code=code_binding('track2_report_evidence','track2_report_sources')),
            lambda output:adapt_report_evidence(admitted['report'],output,
                expected_manifest_sha256=admitted['report_manifest_sha256'],for_regeneration=True))
        check_evidence(evidence,read_scientific(evidence/'manifest.json'))
        result=report_workflow(evidence,operation/'package',expected_evidence_sha256=file_identity(evidence/'manifest.json')['sha256'],
            operation=operation/'workflow',resume=(operation/'workflow').exists() and resume,
            controls=dict(memory_bytes=memory_gib*1024**3,deadline_seconds=deadline_seconds),
            producer_admission=admitted['producer_admission'],require_review=require_review)
        publish_json(operation/'result.json',result)
        return result['outcome']['exit_code']
    finally:check_source_metadata(source,admitted)


def main():
    import argparse
    parser=argparse.ArgumentParser()
    parser.add_argument('--source',type=Path,required=True);parser.add_argument('--operation',type=Path,required=True)
    parser.add_argument('--classification',choices=('public','controlled'),required=True)
    parser.add_argument('--reporter-sha256',required=True);parser.add_argument('--resume',action='store_true')
    parser.add_argument('--recover-unfinished',action='store_true')
    parser.add_argument('--memory-gib',type=int,default=4);parser.add_argument('--deadline-seconds',type=float)
    parser.add_argument('--require-review',action='store_true');args=parser.parse_args()
    return regenerate(args.source,args.operation,expected_classification=args.classification,
        reporter_sha256=args.reporter_sha256,resume=args.resume,memory_gib=args.memory_gib,
        deadline_seconds=args.deadline_seconds,require_review=args.require_review,recover_unfinished=args.recover_unfinished)


if __name__=='__main__':sys.exit(main())
