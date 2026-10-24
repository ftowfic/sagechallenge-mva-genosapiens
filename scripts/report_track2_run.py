#!/usr/bin/env python3
"""Regenerate a completed Track 2 report without reopening scientific inputs."""
import argparse
import fcntl
import json
import os
from pathlib import Path
import signal
import sys
import time
import uuid

ROOT=Path(__file__).resolve().parents[1]
sys.dont_write_bytecode=True
sys.path[:0]=[str(ROOT),str(ROOT/'src')]


def parser():
    from scripts.run_track2 import positive_hours
    value=argparse.ArgumentParser(description=__doc__)
    value.add_argument('--run',type=Path,required=True)
    value.add_argument('--output',type=Path,required=True)
    value.add_argument('--resume',action='store_true')
    value.add_argument('--recover-unfinished',action='store_true',
        help='Recover supported historical production without claiming original report completion.')
    value.add_argument('--refresh-reporter',action='store_true',
        help='Preserve the old reporter and select current code after failure before source admission. Requires --resume.')
    value.add_argument('--settings',type=Path)
    value.add_argument('--require-review',action='store_true')
    value.add_argument('--memory-gib',type=int,default=4)
    value.add_argument('--max-hours',type=positive_hours)
    return value


def worker_options(source, operation, software, destination):
    from scripts.run_track2 import mount
    options=['--mount',mount(source,'/source'), '--mount',mount(software,'/software'),
        '--mount',mount(operation,'/work',False),'--env','PYTHONPATH=/software/src',
        '--env','TMPDIR=/work/scratch','--env','SQLITE_TMPDIR=/work/scratch']
    # Publication is a host rename. The worker sees only its private workspace,
    # never a writable mount of the parent containing the original source run.
    if destination.exists():options+=['--mount',mount(destination,'/work/package')]
    return options


def historical_resume(args, *, public=False):
    """Route retained historical production before scientific input resolution."""
    from scripts.run_track2 import read_json
    if not args.resume or args.output is None:return None
    scientific=('bundle','reviewed_plan') if public else ('vcf','hpo','concordance','diagnosis','sample','review','mechanism_review')
    if any(getattr(args,key,None) is not None for key in scientific):return None
    root=args.output.absolute();manifest=root/'synthesis/controlled-work/report/manifest.json'
    if not manifest.is_file():return None
    from genosapiens.track2_scientific_reader import read_scientific
    schema=read_scientific(manifest).get('schema')
    if schema not in {'genosapiens.track2-research-manifest/v1','genosapiens.track2-scientific-completion/v1',
                       'genosapiens.track2-challenge-package/v1'}:return None
    route=read_json(root/'current-routing.json')
    output=root.with_name(root.name+'-report-revision')
    command=['--run',str(root),'--output',str(output),'--memory-gib',str(args.memory_gib)]
    if not route.get('report_complete'):command.append('--recover-unfinished')
    if output.with_name(output.name+'.operation').exists():command.append('--resume')
    if getattr(args,'settings',None) is not None:command+=['--settings',str(args.settings)]
    if getattr(args,'require_review',False):command.append('--require-review')
    if getattr(args,'max_hours',None) is not None:command+=['--max-hours',str(args.max_hours)]
    print('The historical report remains preserved. Creating or resuming its separate challenge revision.',flush=True)
    return main(command)


def main(argv=None):
    from scripts.run_track2 import Docker,private,read_json,atomic
    from genosapiens.track2_operator_inputs import require
    from genosapiens.track2_report_software import select_reporter
    from genosapiens.track2_report_stages import atomic_publish
    from genosapiens.track2_scientific_storage import file_identity
    os.umask(0o077);args=parser().parse_args(argv)
    require(os.getuid()!=0 and os.getgid()!=0, 'Run as your ordinary non-root account.')
    require(1<=args.memory_gib<=65536, '--memory-gib requires an integer from 1 through 65536.')
    require(not args.refresh_reporter or args.resume, '--refresh-reporter requires --resume.')
    source=private(args.run,directory=True);destination=args.output.absolute()
    operation=destination.with_name(destination.name+'.operation')
    for path in (destination,operation):
        require(path.resolve()==path and not path.is_relative_to(source) and not source.is_relative_to(path),
                'Choose report output independent of the complete source run, without links.')
    snapshot=read_json(source/'preparation/snapshot.json')
    classification=snapshot.get('classification','public' if (source/'run.json').is_file() else 'controlled')
    require(classification in {'public','controlled'}, 'Unsupported source classification.')
    settings=None
    if classification=='controlled':
        if args.settings is None:
            current=read_json(source/'current-run.json')
            reference=Path(current['reference_root'])
            args.settings=reference/'controlled-work'/('.track2-operator-'+str(os.getuid()))/'command/track2-settings.json'
        settings=read_json(private(args.settings))
        from genosapiens.track2_operator import authorize
        authorize(settings)
        inventory=Path(settings['inventory']).absolute()
        require(not inventory.is_relative_to(source) and not inventory.is_relative_to(destination)
            and not inventory.is_relative_to(operation), 'Controlled inventory must remain outside source and report outputs.')
    require(args.resume or not operation.exists() and not destination.exists(), 'Existing report output requires --resume.')
    operation.mkdir(mode=0o700,parents=True,exist_ok=True);private(operation,directory=True)
    fd=os.open(operation/'.launcher.lock',os.O_CREAT|os.O_RDWR|os.O_NOFOLLOW,0o600);inventory_lock=None
    try:
        try:fcntl.flock(fd,fcntl.LOCK_EX|fcntl.LOCK_NB)
        except BlockingIOError:raise ValueError('Another invocation owns this report destination.') from None
        marker=dict(schema='genosapiens.track2-report-launch/v1',source=str(source),destination=str(destination),
            classification=classification,settings=None if settings is None else str(args.settings.absolute()))
        if args.recover_unfinished:marker['recover_unfinished']=True
        marker_path=operation/'launch.json'
        if marker_path.exists():
            require(args.resume and read_json(marker_path)==marker, 'Report resume has another source or destination.')
        else:
            require(not args.resume and not destination.exists(), 'The report lacks matching operation ownership.')
            atomic(marker_path,marker)
        if settings is not None:
            inventory_parent=private(Path(settings['inventory']).parent,directory=True)
            inventory_lock=os.open(inventory_parent/'.track2.lock',os.O_CREAT|os.O_RDWR|os.O_NOFOLLOW,0o600)
            try:fcntl.flock(inventory_lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
            except BlockingIOError:raise ValueError('Another invocation owns the controlled inventory.') from None
        software,reporter=select_reporter(ROOT,operation,refresh=args.refresh_reporter,published=destination.exists())
        docker=Docker(reporter['runtime_image'],args.memory_gib)
        previous={};started=time.monotonic()
        def stopped(signum,frame):
            docker.cleanup()
            raise TimeoutError('The report invocation deadline expired.') if signum==signal.SIGALRM else KeyboardInterrupt()
        for sig in (signal.SIGINT,signal.SIGTERM,signal.SIGALRM):previous[sig]=signal.signal(sig,stopped)
        if args.max_hours is not None:signal.setitimer(signal.ITIMER_REAL,args.max_hours*3600)
        attempt=operation/('invocation-'+uuid.uuid4().hex);attempt.mkdir(mode=0o700)
        (operation/'scratch').mkdir(mode=0o700,exist_ok=True)
        atomic(attempt/'started.json',dict(status='incomplete',reporter=reporter,classification=classification))
        try:
            command=['-m','genosapiens.track2_report_command','--source','/source','--operation','/work',
                '--classification',classification,'--reporter-sha256',reporter['sha256'],'--memory-gib',str(args.memory_gib)]
            if args.resume:command+=['--resume']
            if args.recover_unfinished:command+=['--recover-unfinished']
            if args.require_review:command+=['--require-review']
            if args.max_hours is not None:command+=['--deadline-seconds',str(args.max_hours*3600)]
            code=docker.run(worker_options(source,operation,software,destination),command,log=attempt/'worker.log')
            require(code in (0,2), 'Report regeneration failed. The operation log and partial work remain available.')
            result=read_json(operation/'result.json')
            require(result['status']=='complete' and result['outcome']['exit_code']==code,
                    'The report worker lacks matching finalization proof.')
            staged=destination if destination.exists() else operation/'package'
            require(file_identity(staged/'manifest.json')['sha256']==result['report']['sha256'],
                    'The staged report differs from its publication identity.')
            if not destination.exists():atomic_publish(staged,destination)
            atomic(operation/'current-routing.json',dict(schema='genosapiens.track2-report-routing/v1',
                source_run=str(source),active_report=str(destination/'report.md'),
                manifest_sha256=result['report']['sha256'],report_complete=True,finalization_exit_code=code))
            atomic(attempt/'outcome.json',dict(status='complete',exit_code=code,report=result['report'],
                elapsed_seconds=time.monotonic()-started))
            print('Report: '+str(destination/'report.md'),flush=True)
            return code
        except BaseException as error:
            atomic(attempt/'outcome.json',dict(status='incomplete',error_type=type(error).__name__,
                elapsed_seconds=time.monotonic()-started))
            raise
        finally:
            docker.cleanup();signal.setitimer(signal.ITIMER_REAL,0)
            for sig,handler in previous.items():signal.signal(sig,handler)
            if settings is not None:
                from genosapiens.track2_storage import ControlledStore
                from genosapiens.track2_operator import record_outputs
                store=ControlledStore(settings['reference_root'],settings['inventory'],settings['case_id'])
                record_outputs(store,operation)
                if destination.exists():record_outputs(store,destination)
    finally:
        if inventory_lock is not None:os.close(inventory_lock)
        os.close(fd)


if __name__=='__main__':sys.exit(main())
