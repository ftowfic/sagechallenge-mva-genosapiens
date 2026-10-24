#!/usr/bin/env python3
"""Create or finalize a separate Track 2 report from retained evidence."""
import argparse
import json
from pathlib import Path
import signal
import sys
import uuid

ROOT=Path(__file__).resolve().parents[1]
sys.path[:0]=[str(ROOT),str(ROOT/'src')]


def parser():
    value=argparse.ArgumentParser(description=__doc__)
    commands=value.add_subparsers(dest='command',required=True)
    create=commands.add_parser('create',help='Create a new revision and finalize its current scientific outcome.')
    create.add_argument('--source',type=Path,required=True)
    create.add_argument('--output',type=Path,required=True)
    create.add_argument('--expected-source-manifest-sha256',required=True)
    create.add_argument('--production-record',type=Path)
    create.add_argument('--expected-production-sha256')
    create.add_argument('--calibration-bank',type=Path)
    create.add_argument('--calibration-resource',type=Path)
    create.add_argument('--expected-calibration-resource-sha256')
    create.add_argument('--resume',action='store_true')
    finish=commands.add_parser('finalize',help='Check an existing revision with an independently retained publication identity.')
    finish.add_argument('--report',type=Path,required=True)
    finish.add_argument('--expected-manifest-sha256',required=True)
    for command in (create,finish):
        command.add_argument('--operation',type=Path)
        command.add_argument('--require-review',action='store_true')
        command.add_argument('--timeout-seconds',type=int,default=900)
    return value


def docker_renderer(operation):
    from scripts.run_track2 import Docker,mount
    from genosapiens.track2_challenge_overviews import pinned_runtime
    def render(script,data,destination):
        docker=Docker(pinned_runtime(),4)
        options=['--entrypoint','Rscript','--mount',mount(script,'/renderer.R'),
            '--mount',mount(destination,'/figures',False)]
        log=operation/('renderer-'+uuid.uuid4().hex+'.log')
        try:
            code=docker.run(options,['/renderer.R','/figures/'+data.name,'/figures'],log=log)
            if code:raise RuntimeError('The pinned overview renderer failed. Its operational log remains available.')
        finally:docker.cleanup()
    return render


def main(argv=None):
    p=parser();args=p.parse_args(argv)
    if not 1<=args.timeout_seconds<=86400:p.error('--timeout-seconds requires an integer from 1 through 86400.')
    from genosapiens.track2_calibration_bank import _digest
    if args.command=='create':
        _digest(args.expected_source_manifest_sha256)
        if bool(args.production_record)!=bool(args.expected_production_sha256):p.error('Supply both production record options together.')
        if args.expected_production_sha256:_digest(args.expected_production_sha256)
        supplied=[args.calibration_bank,args.calibration_resource,args.expected_calibration_resource_sha256]
        if any(supplied) and not all(supplied):p.error('Supply all three calibration options together.')
        if args.expected_calibration_resource_sha256:_digest(args.expected_calibration_resource_sha256)
    else:_digest(args.expected_manifest_sha256)
    def deadline(signum,frame):raise TimeoutError('The revision invocation reached its declared deadline.')
    def cancelled(signum,frame):raise KeyboardInterrupt('The revision invocation was cancelled.')
    previous_term=signal.signal(signal.SIGTERM,cancelled)
    previous=signal.signal(signal.SIGALRM,deadline);signal.setitimer(signal.ITIMER_REAL,args.timeout_seconds)
    try:
        from genosapiens.track2_revision_finalization import publish_and_finalize_revision,finalize_revision
        if args.command=='finalize':
            outcome=finalize_revision(args.report,expected_manifest_sha256=args.expected_manifest_sha256,
                operation=args.operation,require_review=args.require_review)
            print(json.dumps(dict(status='finalized',outcome=outcome.as_dict())),flush=True)
            return outcome.exit_code
        output=args.output.absolute()
        operation=args.operation.absolute() if args.operation else output.with_name(output.name+'.operation')
        calibration=None if args.calibration_bank is None else dict(bank=args.calibration_bank,
            resource_contract=args.calibration_resource,expected_resource_sha256=args.expected_calibration_resource_sha256)
        if calibration is None:
            from genosapiens.track2_report_workflow import revise_report
            result=revise_report(args.source,output,expected_source_manifest_sha256=args.expected_source_manifest_sha256,
                production=None if args.production_record is None else dict(record=args.production_record,
                    expected_sha256=args.expected_production_sha256),operation=operation,resume=args.resume,
                figure_renderer=docker_renderer(operation),require_review=args.require_review,
                controls=dict(deadline_seconds=args.timeout_seconds))
            print(json.dumps(dict(status='finalized',**result)),flush=True)
            return result['outcome']['exit_code']
        if args.resume:p.error('Legacy explicit-bank assembly requires a new destination.')
        def published(binding):print(json.dumps(dict(status='assembled',report_binding=binding)),flush=True)
        result=publish_and_finalize_revision(args.source,output,
            expected_source_manifest_sha256=args.expected_source_manifest_sha256,calibration=calibration,
            figure_renderer=docker_renderer(operation),require_review=args.require_review,
            operation=operation,on_publish=published,
            production=None if args.production_record is None else dict(record=args.production_record,expected_sha256=args.expected_production_sha256))
        print(json.dumps(dict(status='finalized',**result)),flush=True)
        return result['outcome']['exit_code']
    finally:
        signal.setitimer(signal.ITIMER_REAL,0);signal.signal(signal.SIGALRM,previous);signal.signal(signal.SIGTERM,previous_term)


if __name__=='__main__':
    try:raise SystemExit(main())
    except (ValueError,OSError,RuntimeError,TimeoutError) as error:
        print(str(error),file=sys.stderr);raise SystemExit(3)
