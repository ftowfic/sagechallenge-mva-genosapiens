#!/usr/bin/env python3
"""Run reviewed public genotypes through current native discovery and reporting."""
import argparse
import fcntl
import hashlib
import json
import os
from pathlib import Path
import shutil
import sys
import tempfile

ROOT=Path(__file__).resolve().parents[1]
sys.path[:0]=[str(ROOT),str(ROOT/'src')]


def prepare_report_metadata(repository, reference, software, *, memory_gib=4):
    """Copy prepared Open Targets mappings, or parse the installed source once."""
    from genosapiens.track2_report_gene_metadata import (
        FILENAME, copy_prepared_gene_metadata, installed_gene_source)
    if copy_prepared_gene_metadata(repository, reference, Path(software)/FILENAME):
        return
    from scripts.run_track2 import Docker, mount
    pin = installed_gene_source(repository, reference)
    contract = json.loads(pin.contract_path.read_bytes())
    release = Path(reference)/contract['output_prefix']
    docker = Docker(json.loads((Path(repository)/'config/track2/diagnosis-runtime.json').read_bytes())['image'], memory_gib)
    # Only existing public target inputs enter this preparation worker.
    options = ['--env','PYTHONPATH=/software/src',
        '--mount',mount(software,'/software'), '--mount',mount(software,'/output',False),
        '--mount',mount(pin.manifest_path), '--mount',mount(release)]
    try:
        code = docker.run(options, ['-c',
            'from genosapiens.track2_report_gene_metadata import installed_gene_source, prepare_gene_metadata, FILENAME; '
            'from pathlib import Path; '
            'prepare_gene_metadata(installed_gene_source(Path("/software"), Path(__import__("sys").argv[1])), Path("/output")/FILENAME)',
            str(reference)], log=Path(software).parent/'gene-metadata-preparation.log')
        if code:
            raise RuntimeError('The installed gene metadata projection failed. Its output remains incomplete.')
    finally:
        docker.cleanup()


def worker(phase, *, root=Path('/output'), software=Path('/software'),
           reference=Path('/discovery-resources'), published=Path('/published'),lincs_workers=2,figure_workers=2):
    from genosapiens.track2_run_state import atomic_record
    from genosapiens.track2_scientific_reader import read_scientific
    from genosapiens.track2_scientific_storage import file_identity, publish_json
    from genosapiens.track2_public_discovery import prepare_public_plan, checked_public_provenance
    from genosapiens.track2_discovery_sources import native_snapshot
    from genosapiens.track2_synthesis import policies
    from genosapiens.track2_published_genotype import canonical, sha256
    run=json.loads((root/'run.json').read_text())
    descriptor=json.loads((software/'discovery-resources.json').read_text())
    if phase=='prepare':
        plan=prepare_public_plan(published,published/run['reviewed_plan_relative'],
            descriptor,reference,root/'plan',run['images'])
        bundle=checked_public_provenance(plan)
        bindings=dict(input=bundle.input_sha256,input_review=bundle.review_sha256,
            plan=file_identity(root/'plan/plan.json')['sha256'],plan_review=None,
            code=run['software_sha256'],policy=sha256(canonical(policies()[0])),**bundle.source_hashes)
        snapshot=native_snapshot(plan,bindings,dict(images=run['images'],code_revision=run['software_sha256']),descriptor)
        publish_json(root/'preparation/snapshot.json',snapshot)
        atomic_record(root/'construction/run.json',dict(schema='genosapiens.track2-public-input-admission/v1',
            input_origin='published_patient_genotype',input_sha256=bundle.input_sha256,
            input_review_sha256=bundle.review_sha256,source_hashes=bundle.source_hashes))
        nested=plan['exploratory_discovery']
        summary=dict(objectives=len(nested['objectives']),jobs=len(nested['genetic_selections']),
            cells=[dict(cell=s['cell_id'],method=s['method'],drug_signatures=len(s['drug_signature_ids']))
                for s in nested['genetic_selections']])
        atomic_record(root/'scope.json',summary);print(json.dumps(summary),flush=True)
        return 0
    if phase=='discover':
        from genosapiens.track2_discovery import discover_file
        result=discover_file(root/'plan/plan.json',reference,root/'discovery',lincs_workers=lincs_workers)
        return 0 if result['completion_state']=='complete' else 2
    if phase=='figures':
        from genosapiens.track2_lincs_report import render_native_jobs
        render_native_jobs(root/'discovery',root/'discovery-figures',figure_workers=figure_workers);return 0
    if phase=='report':
        from genosapiens.track2_report_gene_metadata import FILENAME
        if not (software/FILENAME).is_file():
            raise ValueError('The frozen public software lacks installed gene metadata. Create a report revision with the prepared mappings.')
        from genosapiens.track2_report_finalization import prepare_synthesis, save_production, finalize_report
        from genosapiens.track2_discovery_report import export_native_evidence
        from genosapiens.track2_report_fresh import run_fresh_report
        prepare_synthesis(root,software,native=True)
        route=json.loads((root/'current-routing.json').read_text())
        if route['report_complete']:return route['finalization_exit_code']
        if route.get('synthesis_complete'):return 0
        plan=read_scientific(root/'plan/plan.json');bundle=checked_public_provenance(plan)
        snapshot=json.loads((root/'preparation/snapshot.json').read_text())
        def export(evidence):
            return export_native_evidence(bundle,plan,None,snapshot['bindings'],root/'plan/plan.json',
            evidence,input_snapshot_path=root/'preparation/snapshot.json',
            expected_snapshot_sha256=file_identity(root/'preparation/snapshot.json')['sha256'],
            controlled_store=None,discovery_path=root/'discovery',reference_root=reference,
            figures_root=root/'discovery-figures')
        result=run_fresh_report(root/'synthesis/controlled-work/report',evidence_binding=route['report_binding'],export=export,
            controls=dict(memory_bytes=int(os.environ.get('GENOSAPIENS_TRACK2_REPORT_MEMORY_BYTES',4*1024**3))))
        atomic_record(root/'synthesis/controlled-work/result.json',result)
        save_production(root/'synthesis',route['report_binding'])
        return 0
    if phase=='finalize':
        from genosapiens.track2_report_finalization import finalize_report
        return finalize_report(root,software)
    raise ValueError('Unknown public analysis phase.')


def parser():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--bundle',type=Path)
    parser.add_argument('--reviewed-plan',type=Path)
    parser.add_argument('--output',type=Path)
    parser.add_argument('--reference-root',type=Path,
        help='Host reference directory. Fresh default: /mnt/reference. Resume restores the saved root.')
    parser.add_argument('--memory-gib',type=int,default=4)
    parser.add_argument('--lincs-workers',type=int,choices=(1,2),default=2)
    parser.add_argument('--figure-workers',type=int,choices=(1,2),default=2)
    parser.add_argument('--resume',action='store_true')
    parser.add_argument('--prepare-only',action='store_true')
    parser.add_argument('--worker',choices=['prepare','discover','figures','report','finalize'])
    return parser


def resume_bindings(args, root):
    from scripts.portable_paths import bound_json, reference_selection
    run=bound_json(root/'run.json')
    if run.get('schema') != 'genosapiens.track2-public-run/v1':
        raise ValueError('The public run binding is incompatible.')
    args.reference_root=reference_selection(args.reference_root, run.get('reference_root'), resume=True)
    if args.bundle is not None and str(args.bundle.resolve(strict=True)) != run['bundle']:
        raise ValueError('The explicit public bundle conflicts with the saved run.')
    if args.reviewed_plan is not None:
        expected=Path(run['bundle'])/run['reviewed_plan_relative']
        if args.reviewed_plan.resolve(strict=True) != expected:
            raise ValueError('The explicit reviewed plan conflicts with the saved run.')
    return run


def main(argv=None):
    p=parser()
    args=p.parse_args(argv)
    if args.worker:return worker(args.worker,lincs_workers=args.lincs_workers,figure_workers=args.figure_workers)
    from scripts.report_track2_run import historical_resume
    revised=historical_resume(args,public=True)
    if revised is not None:return revised
    if not 1 <= args.memory_gib <= 65536:p.error('--memory-gib requires an integer from 1 through 65536.')
    from scripts.run_track2 import Docker,mount,discovery_mounts,atomic
    from scripts.prepare_track2_discovery import ensure_discovery_resources
    os.umask(0o077)
    if args.output is None:p.error('--output is required')
    root=args.output.absolute()
    if root.resolve()!=root:p.error('Use an output path without symbolic links.')
    if args.resume:
        resume_bindings(args, root)
    else:
        from scripts.portable_paths import reference_selection
        args.reference_root=reference_selection(args.reference_root)
    root.mkdir(mode=0o700,parents=True,exist_ok=args.resume)
    with (root/'.public-run.lock').open('a') as lock:
        fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
        software=root/'software'
        if not args.resume:
            if args.bundle is None or args.reviewed_plan is None:p.error('Supply the public bundle and reviewed plan.')
            bundle=args.bundle.resolve();source_plan=args.reviewed_plan.resolve()
            relative=source_plan.relative_to(bundle)
            # Check the public declaration before mounting any source bundle.
            if json.loads((bundle/'genotype.json').read_text()).get('classification')!='public':
                raise ValueError('This command accepts only public genotype evidence.')
            installed=ensure_discovery_resources(ROOT,args.reference_root,memory_gib=args.memory_gib)
            software.mkdir(mode=0o700)
            from genosapiens.track2_package_files import runtime_files
            for source in runtime_files(ROOT):
                destination=software/source.relative_to(ROOT)
                destination.parent.mkdir(parents=True,exist_ok=True)
                shutil.copyfile(source,destination)
            shutil.copyfile(args.reference_root/installed['path'],software/'discovery-resources.json')
            prepare_report_metadata(ROOT,args.reference_root,software,memory_gib=args.memory_gib)
            hashes={str(p.relative_to(software)):hashlib.sha256(p.read_bytes()).hexdigest()
                for p in sorted(software.rglob('*')) if p.is_file()}
            images={kind:json.loads((ROOT/'config/track2'/(kind+'-runtime.json')).read_text())['image']
                for kind in ('diagnosis','discovery')}
            run=dict(schema='genosapiens.track2-public-run/v1',bundle=str(bundle),
                reviewed_plan_relative=str(relative),reference_root=str(args.reference_root.absolute()),
                images=images,software_files=hashes,
                software_sha256=hashlib.sha256(json.dumps(hashes,sort_keys=True,separators=(',',':')).encode()).hexdigest())
            atomic(root/'run.json',run)
        else:
            run=resume_bindings(args, root)
            for name,digest in run['software_files'].items():
                if hashlib.sha256((software/name).read_bytes()).hexdigest()!=digest:
                    raise ValueError('The frozen public analysis software changed.')
            if args.bundle is not None and str(args.bundle.resolve())!=run['bundle']:
                raise ValueError('The resumed public input differs.')
        reporting_software=software
        if args.resume:
            from genosapiens.track2_report_software import public_reporting_software
            reporting_software=public_reporting_software(ROOT,software,root/'reporting-software',run)
        logs=root/'logs';logs.mkdir(exist_ok=True)
        phases=(['prepare'] if not (root/'preparation/snapshot.json').is_file() else [])
        if not args.prepare_only:phases+=['discover','figures','report','finalize']
        for phase in phases:
            selected_software=reporting_software if phase in {'report','finalize'} else software
            number=1+len(list(logs.glob(phase+'-*.log')))
            options=['--env','PYTHONDONTWRITEBYTECODE=1','--env','PYTHONPATH=/software/src',
                '--env','GENOSAPIENS_TRACK2_REPORT_MEMORY_BYTES='+str(args.memory_gib*1024**3),
                '--mount',mount(selected_software,'/software'), '--mount',mount(root,'/output',False),
                *discovery_mounts(selected_software,Path(run['reference_root']),matrix=phase=='discover')]
            if phase=='prepare':options+=['--mount',mount(run['bundle'],'/published')]
            image=run['images']['discovery' if phase in {'discover','figures','report'} else 'diagnosis']
            docker=Docker(image,args.memory_gib)
            print('Public analysis stage: '+phase+'.',flush=True)
            try:
                command=['/software/scripts/run_track2_public.py','--lincs-workers',str(args.lincs_workers),'--figure-workers',str(args.figure_workers),'--worker',phase]
                log=logs/(phase+'-'+str(number)+'.log')
                if phase in {'report','finalize'}:
                    # Scientific JSON readers keep temporary SQLite indexes.
                    # Use the run filesystem instead of the bounded /tmp tmpfs.
                    with tempfile.TemporaryDirectory(prefix='.report-tmp-',dir=root) as scratch:
                        temporary_options=[*options,
                            '--mount',mount(scratch,'/track2-scratch',False),
                            '--env','TMPDIR=/track2-scratch',
                            '--env','SQLITE_TMPDIR=/track2-scratch']
                        code=docker.run(temporary_options,command,log=log)
                else:
                    code=docker.run(options,command,log=log)
            except BaseException as error:
                try:
                    atomic(root/'current-public-outcome.json',dict(stage=phase,
                        state='interrupted' if isinstance(error,(KeyboardInterrupt,TimeoutError)) else 'failed'))
                except BaseException:
                    print('The public failure record is unavailable.',file=sys.stderr)
                if phase in {'overviews','finalize'}:
                    from genosapiens.track2_report_finalization import record_finalization_failure
                    record_finalization_failure(root,error)
                raise
            finally:docker.cleanup()
            accepted=code==0 or code==2 and phase in {'discover','report','finalize'}
            atomic(root/'current-public-outcome.json',dict(stage=phase,exit_code=code,
                state='complete' if code==0 else 'partial' if accepted else 'failed'))
            if not accepted:
                if phase in {'overviews','finalize'}:
                    from scripts.run_track2 import record_incomplete_finalization
                    record_incomplete_finalization(software,root,code)
                    report=root/'synthesis/controlled-work/report/report.md'
                    if report.is_file():print('Unfinalized report: '+str(report),flush=True)
                return code
            if phase=='finalize':
                route=json.loads((root/'current-routing.json').read_bytes())
                if not route.get('report_complete') or route.get('finalization_exit_code')!=code:
                    atomic(root/'current-public-outcome.json',dict(stage=phase,exit_code=1,state='failed'))
                    raise ValueError('Public report finalization has no matching completion.')
                print('Public research report: '+str(root/'synthesis/controlled-work/report/report.md'),flush=True)
                return code
        return 0


if __name__=='__main__':
    raise SystemExit(main())
