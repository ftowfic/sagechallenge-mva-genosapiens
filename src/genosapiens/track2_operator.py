"""Operator-owned orchestration of the existing diagnosis-guided scientific APIs.

The host owns Docker. These fixed phases never start Docker or contact a model.
Only the construction phase receives the explicitly selected raw inputs.
"""

from datetime import datetime, timezone
import os
from pathlib import Path
import socket
import time
import uuid

from genosapiens.track2_operator_inputs import (OperatorError, require, parser, private_path,
    small_json, write_new, intake, review_worksheet, read_review, input_review_record, plan_review_record)
from genosapiens.track2_published_genotype import canonical, sha256
from genosapiens.track2_terminal import assessment_outcome, publish_outcome, pending_review


def authorize(settings):
    from genosapiens.participant_governance import load_signed_participant_governance, authorize_controlled_access
    from genosapiens.incident_response import assert_processing_allowed
    require(os.getuid()!=0, 'Run as your ordinary non-root account.')
    governance,_ = load_signed_participant_governance(settings['governance'], settings['verification_key'],
        settings['verification_key_sha256'], datetime.now(timezone.utc), settings['rules_sha256'])
    authorize_controlled_access(governance, settings['participant_id'], os.getuid())
    assert_processing_allowed(settings['incident_root'], settings['environment_id'])
    require(sha256(Path(settings['rules']).read_bytes())==settings['rules_sha256'],
            'The installed rules binding changed. Refresh the authorized installation.')
    return governance


def storage(settings, args, *, local=False):
    from genosapiens.track2_storage import ControlledStore
    bindings=[]
    for name in ('output','vcf','hpo','concordance','review','mechanism_review'):
        host=os.environ.get('GENOSAPIENS_TRACK2_HOST_'+name.upper())
        if host is not None and getattr(args,name) is not None:
            bindings.append((getattr(args,name),Path(host)))
    return ControlledStore(settings['reference_root'],
        args.output/'current-inventories'/((local if isinstance(local,str) else 'construct')+'.json')
            if local else settings['inventory'],settings['case_id'],
        output_root=args.output, input_paths=tuple(p for p in
            (args.vcf,args.hpo,args.concordance,args.review,args.mechanism_review) if p is not None),path_bindings=tuple(bindings))


def record_source(store, path, digest, signature):
    """Inventory an already authenticated source without hashing the VCF again."""
    from genosapiens.governance import CaseInventory, write_inventory
    from genosapiens.track2_vcf_selection import signature as current_signature
    path=store.path(path)
    require(current_signature(path)==signature, 'The selected source changed during construction.')
    relative=store.relative_path(path)
    record=store.record_for(path,kind='track2-operator-input',digest=digest,size=path.stat().st_size)
    records={r.relative_path:r for r in store.inventory().artifacts}
    require(relative not in records or records[relative]==record, 'An inventoried source identity changed.')
    records[relative]=record
    write_inventory(CaseInventory(store.case_id,tuple(records[k] for k in sorted(records))),store.inventory_path)


def record_outputs(store, root):
    """Only this selected run's outputs; never enumerate the shared storage root."""
    files=[]
    for path in sorted(Path(root).rglob('*')):
        require(not path.is_symlink(), 'Run outputs must not contain symbolic links.')
        if path.is_file() and path!=store.inventory_path:
            path.chmod(0o600);files.append(path)
        elif path.is_dir():
            path.chmod(0o700)
    from genosapiens.governance import CaseInventory, write_inventory
    known={row.relative_path:row for row in store.inventory().artifacts}
    present=set()
    for path in files:
        relative=store.relative_path(path)
        present.add(relative)
        known[relative]=store.record_for(path,kind='track2-operator-output')
    host_root=store.host_path(root)
    for relative, previous in tuple(known.items()):
        host=Path(previous.operator_path) if previous.operator_path else store.reference_root/relative
        if host.is_relative_to(host_root) and relative not in present and previous.status=='present':
            known[relative]=previous.mark_deleted()
    write_inventory(CaseInventory(store.case_id,tuple(known[k] for k in sorted(known))),store.inventory_path)


def construct(args, settings, package_manifest):
    from genosapiens.reference_concordance import (ReferenceConcordanceReport, load_reference_identity,
        open_installed_reference, evaluate_vcf_reference_concordance, concordance_report_as_json)
    from genosapiens.vcf_qc import _file_sha256
    from genosapiens.track2_vcf_selection import signature, select_vcf
    from genosapiens.track2_operator_resources import open_resources, interpret, phenotype
    from genosapiens.track2_diagnosis_scope import build_diagnosis_scope
    from genosapiens.track2_scoped_genotypes import construct_scoped_genotypes
    from genosapiens.track2_scoped_mechanisms import assemble_scoped_mechanisms
    from genosapiens.track2_diagnosis_builder import construct_components
    from genosapiens.track2_diagnosis_package import write_diagnosis_package, load_diagnosis_package, render_diagnosis_dossier
    # Scientific workers own a local inventory. Only the separate control phase
    # receives the host's shared inventory directory.
    store=storage(settings,args,local=True)
    output=store.path(args.output);package=Path(settings['package'])
    private_path(args.vcf)
    if args.hpo is not None: private_path(args.hpo)
    if args.concordance is not None: private_path(args.concordance)
    selected,hpo_bytes=intake(args.vcf,args.sample,args.hpo)
    bound={'source_path':str(store.host_path(args.vcf)),'source_signature':signature(args.vcf), 'sample':selected,
        'diagnosis':args.diagnosis,'hpo_path':None if args.hpo is None else str(store.host_path(args.hpo)),
        'hpo_sha256':None if hpo_bytes is None else sha256(hpo_bytes),
        'explicit_concordance':None if args.concordance is None else sha256(Path(args.concordance).read_bytes()),
        'package_sha256':settings['package_sha256'],'case_id':settings['case_id']}
    from genosapiens.track2_run_state import RunState, atomic_record
    from genosapiens.track2_construction_dependencies import construction_dependencies
    state=RunState(output)
    inputs={k:v for k,v in bound.items() if k!='package_sha256'}
    dependency={'inputs':inputs,'image':settings['image'],
        'dependencies':construction_dependencies(package)}
    completed,reason=state.completed('construct',dependency)
    if completed is not None:
        saved=small_json(output/'construction/run.json')
        bundle=load_diagnosis_package(output/'construction/package',require_review=False,
            expected_input_sha256=saved['input_sha256'])
        print('Construction is complete. Reuse the bound package.',flush=True)
        return
    print('Restart construction: '+reason+'.',flush=True)
    target=state.restart('construct',reason)
    store.inventory_path.unlink(missing_ok=True)
    # Save the selected sample once intake resolves a single-sample VCF.
    current=small_json(output/'current-run.json')
    current['inputs']['sample']=selected
    current['source_versions']={'image':settings['image'],'package_sha256':settings['package_sha256']}
    atomic_record(output/'current-run.json',current)
    runtime={'image':settings['image'],'code_revision':'sha256:'+settings['package_sha256']}
    started=time.monotonic()
    print('Open the installed public resources.',flush=True)
    with open_resources(package/'resources.json',package_manifest['files']['resources.json'],diagnosis=args.diagnosis) as resources:
        require(resources.retention_index is not None,
            'Prepare the shared public selection index before case construction.')
        print('Resolve the diagnosis and freeze its gene and regulatory scope.',flush=True)
        scope=build_diagnosis_scope(args.diagnosis,store=resources.store,gene_source=resources.genes,
            regulatory_index=resources.index,reference=resources.reference,retention_index=resources.retention_index)
        reference=resources.reference;contract=load_reference_identity(reference.contract_path)
        with open_installed_reference(contract,identity_path=reference.contract_path,
                reference_manifest_path=reference.manifest_path,reference_root=reference.root,
                reuse_published_validation=True) as (manifest_sha,fetch,contigs):
            from genosapiens.track2_construction_dependencies import python_dependencies
            concordance_binding={'source_signature':bound['source_signature'],
                'reference_manifest_sha256':manifest_sha,'explicit':bound['explicit_concordance'],
                'code':python_dependencies(package,(
                    ('genosapiens.reference_concordance','evaluate_vcf_reference_concordance'),))}
            completed,reason=state.completed('concordance',concordance_binding)
            if completed is not None:
                concordance=ReferenceConcordanceReport.from_mapping(small_json(output/'concordance/report.json'))
                print('Reuse completed source/reference concordance.',flush=True)
            elif args.concordance is None:
                print('Check the source REF bases against the installed reference.',flush=True)
                concordance=evaluate_vcf_reference_concordance(args.vcf,fetch,contigs,
                    contract=contract,reference_manifest_sha256=manifest_sha)
            else:
                concordance=ReferenceConcordanceReport.from_mapping(small_json(args.concordance))
                digest,_=_file_sha256(args.vcf)
                require(digest==concordance.validated_input_sha256,'The concordance report does not match this source.')
            require(concordance.accepted and concordance.reference_manifest_sha256==manifest_sha
                and concordance.input_sha256==concordance.validated_input_sha256
                and concordance.reference_role=='challenge_alignment'
                and concordance.checked_records==concordance.matching_records==concordance.total_records
                and concordance.unresolved_records==concordance.mismatching_records==0,
                'Complete accepted source/reference concordance is required.')
            if completed is None:
                location=state.restart('concordance',reason)
                write_new(location/'report.json',concordance_report_as_json(concordance).encode())
                state.publish('concordance',binding=concordance_binding,result={})
            write_new(target/'concordance.json',concordance_report_as_json(concordance).encode())
            print('Select small variants within the declared scope.',flush=True)
            selection=select_vcf(args.vcf,scope=scope,concordance=concordance,sample=selected,
                fetch_reference=fetch,reference_mapping=resources.data['reference_mapping'],
                exclusion_root=target/'selection-exclusions',
                classification='controlled',evidence_origin='controlled_case',deadline_seconds=None)
        print('Interpret scoped variants and construct genotype alternatives.',flush=True)
        interpreted=interpret(selection,scope,resources,target,runtime,
            annotation_root=output/'annotation-checkpoint')
        print(f'Retained alleles: {len(interpreted.alleles)}. Deferred alleles: {len(interpreted.deferred)}.',flush=True)
        from genosapiens.track2_functional_sources import bind_gene_sources
        functional_policy=bind_gene_sources(resources.policy,resources.data['functional_catalog'],
            {seed.context.gene_id for seed in interpreted.contexts if seed.context.gene_id is not None})
        genotypes=construct_scoped_genotypes(selection,interpreted,store=resources.store,policy=functional_policy)
        mechanisms=assemble_scoped_mechanisms(genotypes,interpreted,evidence_origin='controlled_case')
        parts=construct_components(scope=scope,selection=selection,interpreted=interpreted,
            genotypes=genotypes,mechanisms=mechanisms,
            phenotype=phenotype(args.hpo,resources,bound['hpo_sha256']))
        from genosapiens.track2_public_matching import load_public_evidence, attach_matching
        print('Match the constructed genotypes to reviewed public evidence.',flush=True)
        parts=attach_matching(parts,load_public_evidence(package/'public-evidence'))
        write_diagnosis_package(target/'package',components=parts,classification='controlled',
            exclusion_source=selection.exclusion_root,
            deferred_rows=interpreted.deferred,
            evidence_origin='controlled_case',runtime=runtime,builder_code_sha256=settings['package_sha256'],
            policy_sha256=scope['policy_sha256'],controlled_store=store)
    bundle=load_diagnosis_package(target/'package',require_review=False)
    write_new(target/'input-dossier.md',render_diagnosis_dossier(bundle).encode())
    (output/'input-review.tsv').unlink(missing_ok=True)
    write_new(output/'input-review.tsv',review_worksheet('input',bundle.input_sha256))
    record_source(store,args.vcf,concordance.validated_input_sha256,bound['source_signature'])
    if args.hpo is not None:
        store.record([args.hpo],kind='track2-operator-hpo')
    write_new(target/'run.json',{'binding':bound,'input_sha256':bundle.input_sha256,
        'author':settings['participant_id'],'construction_seconds':time.monotonic()-started,
        'runtime':runtime})
    state.publish('construct',binding=dependency,result={'input_sha256':bundle.input_sha256})


METHODS={
    'description':'The offline runner constructs a diagnosis-guided package and generates a source-bound report.',
    'llm_disclosure':'No language model receives case inputs or controlled reports.',
    'limitations':'The scope covers supplied small variants. Unknown mechanisms and missing therapeutic evidence remain unassessed.',
    'validation_followup':'Review the findings locally. Final scientific review and release approval remain separate.'}


def review(args,settings,governance):
    from genosapiens.track2_report_finalization import retain_production_routing
    retain_production_routing(args.output)
    from genosapiens.track2_diagnosis_package import load_diagnosis_package, render_diagnosis_dossier
    from genosapiens.track2_scientific_planning import prepare_scientific_plan
    from genosapiens.track2_storage import stage_diagnosis_input
    from genosapiens.track2_synthesis import freeze_input_snapshot
    from genosapiens.track2_freeze import load_frozen_public_envelope
    store=storage(settings,args,local='review');root=store.path(args.output);run=small_json(root/'construction/run.json')
    print('Check the completed construction for assessment.',flush=True)
    package=root/'construction/package';bundle=load_diagnosis_package(package,require_review=False,
        expected_input_sha256=run['input_sha256'])
    dossier=root/'construction/input-dossier.md'
    if not dossier.exists():
        write_new(dossier,render_diagnosis_dossier(bundle).encode())
    plan_dir=root/'plan';plan_path=plan_dir/'plan.json'
    plan_bindings=None
    kind='input' if not bundle.human_review_ready else 'plan'
    if args.review is not None:
        import csv, io
        from genosapiens.bounded_json import read_bounded_bytes
        rows=list(csv.DictReader(io.StringIO(read_bounded_bytes(private_path(args.review), maximum_bytes=65536).decode()),delimiter='\t'))
        require(len(rows)==1 and rows[0].get('object_id') in {'input','plan'},'Use the generated review worksheet.')
        kind=rows[0]['object_id']
    reviewed=package/'review.json' if kind=='input' else plan_dir/'plan-review.json'
    if args.review is not None and not reviewed.exists():
        require(args.resume,'Review intake requires explicit --resume.')
        digest=bundle.input_sha256 if kind=='input' else sha256(plan_path.read_bytes())
        row=read_review(args.review,kind=kind,object_sha256=digest,governance=governance,host_uid=os.getuid())
        require(row['independence']!='independent' or row['reviewer']!=run['author'],
                'The run author must disclose self review.')
        audit=root/'reviews'/(kind+'-'+uuid.uuid4().hex)
        audit.mkdir(parents=True,mode=0o700)
        write_new(audit/'worksheet.tsv',Path(args.review).read_bytes())
        write_new(audit/'receipt.json',row)
        if row['decision']!='accepted':
            print('Review is pending or rejected. Dependent analysis remains stopped.',flush=True)
            return pending_review('review_pending_or_rejected')
        record=input_review_record(bundle,row) if kind=='input' else plan_review_record(plan_path,bundle,row)
        write_new(reviewed,record)
        bundle=load_diagnosis_package(package)
    elif args.review is not None:
        # Completed decisions are immutable. Repeating the same worksheet is
        # harmless, but a new review cannot silently replace an accepted one.
        row=read_review(args.review,kind=kind,object_sha256=bundle.input_sha256 if kind=='input' else sha256(plan_path.read_bytes()),
            governance=governance,host_uid=os.getuid())
        submitted=input_review_record(bundle,row) if kind=='input' else plan_review_record(plan_path,bundle,row)
        existing=small_json(reviewed)
        require(all(submitted[k]==existing[k] for k in submitted if k!='reviewed_at'),
            'This object already has a different accepted review. Its decision remains preserved.')
    if not bundle.human_review_ready and args.require_review:
        print('Review construction/input-dossier.md in this output directory. Save a completed copy of input-review.tsv.',flush=True)
        return pending_review('input_review_required')
    if (Path(settings['package'])/'discovery-resources.json').is_file():
        from genosapiens.track2_discovery_operator import review_native
        return review_native(args,settings,governance,bundle,store,run)
    from genosapiens.track2_run_state import RunState, atomic_record
    from genosapiens.track2_construction_dependencies import python_dependencies
    state=RunState(root)
    mechanism_record=root/'mechanism'/(bundle.input_sha256+'.json')
    if args.mechanism_review is not None:
        from genosapiens.track2_public_matching import load_public_evidence
        from genosapiens.track2_mechanism_review import read_mechanism_review
        submitted=read_mechanism_review(args.mechanism_review,bundle,
            load_public_evidence(Path(settings['package'])/'public-evidence'),
            governance=governance,host_uid=os.getuid())
        if mechanism_record.exists():
            require(submitted==small_json(mechanism_record),
                'This input already has a different accepted mechanism review.')
        else:
            audit=root/'reviews'/('mechanism-'+uuid.uuid4().hex)
            write_new(audit/'worksheet.tsv',Path(args.mechanism_review).read_bytes())
            write_new(audit/'receipt.json',submitted)
            if not submitted['all_decisions_accepted']:
                return pending_review('mechanism_review_pending_or_rejected')
            write_new(mechanism_record,submitted)
    assessment_binding={'input_sha256':bundle.input_sha256,'input_review_sha256':bundle.review_sha256,
        'require_review':args.require_review,'image':settings['image'],
        'public_evidence':sha256((Path(settings['package'])/'public-evidence/evidence.json').read_bytes()),
        'envelope':sha256((Path(settings['package'])/'envelope.json').read_bytes()),
        'mechanism_review':None if not mechanism_record.exists() else sha256(mechanism_record.read_bytes()),
        'code':python_dependencies(settings['package'],(
            ('genosapiens.track2_mechanism_assessment','build_assessment'),
            ('genosapiens.track2_scientific_planning','prepare_scientific_plan'),
            ('genosapiens.track2_synthesis','plan_acceptance')))}
    completed,reason=state.completed('assessment',assessment_binding)
    if completed is None:
        print('Restart assessment: '+reason+'.',flush=True)
        state.restart('assessment',reason)
        (root/'plan-review.tsv').unlink(missing_ok=True)
        from genosapiens.track2_public_matching import load_public_evidence
        from genosapiens.track2_mechanism_review import worksheet, read_mechanism_review
        from genosapiens.track2_mechanism_assessment import build_assessment
        public_root=Path(settings['package'])/'public-evidence'
        print('Open the retained public evidence for assessment.',flush=True)
        public=load_public_evidence(public_root)
        generated=worksheet(bundle,public) if bundle.input_report_ready else None
        if generated is not None:
            (root/'mechanism-review.tsv').unlink(missing_ok=True)
            write_new(root/'mechanism-review.tsv',generated)
        if generated is not None and not mechanism_record.exists() and args.require_review:
            print('Review mechanism-review.tsv with retained case and public evidence.',flush=True)
            return pending_review('mechanism_and_objective_review_required')
        from genosapiens.track2_scientific_planning import frozen_packet
        from genosapiens.track2_mechanism_review import scope_from_assessment
        print('Assess mechanism evidence and objective applicability.',flush=True)
        mechanism=small_json(mechanism_record) if mechanism_record.exists() else {
            'assessment':build_assessment(bundle,automatic_inputs={
                'execution':{'code_sha256':settings['package_sha256'],'image':settings['image']},
                'timestamp':datetime.now(timezone.utc).isoformat(),
                'public_packet':frozen_packet(public_root) if public['summary'] is not None else None}),
            'endorsements':{}}
        mechanism.setdefault('assessment_scope',scope_from_assessment(mechanism['assessment']))
        # prepare_scientific_plan validates the assessment through report_origins.
        frozen=load_frozen_public_envelope(Path(settings['package'])/'envelope.json')
        # Publish completion only after replay validates the current plan.
        prepared_plan=plan_dir
        print('Prepare the scientific plan and record unresolved evidence.',flush=True)
        plan,gaps=prepare_scientific_plan(bundle,assessment=mechanism['assessment'],public=public,
            public_root=public_root,frozen=frozen,endorsements=mechanism['endorsements'],output=prepared_plan,
            assessment_scope=mechanism.get('assessment_scope'),
            draft=not args.require_review,methods={
            **METHODS,'effort':f"Construction took {run['construction_seconds']:.3f} seconds. "
                'Review records contain submission times only after a reviewer submits a worksheet.'})
        digest=write_new(prepared_plan/'plan.json',plan)
        write_new(prepared_plan/'planning-gaps.json',gaps)
        write_new(prepared_plan/'review-notes.md',(
            '# Plan review\n\nReview the input dossier and the structured plan.\n\n'
            'This plan preserves genotype alternatives and unresolved mechanisms.\n\n'
            'Qualified objectives use accepted mechanism evidence and compatible frozen public evidence. '
            'The planning gaps identify unresolved scenario work.\n\n'
            'Acceptance confirms this stated scope. It does not establish a treatment or case applicability.\n').encode())
        from genosapiens.track2_synthesis import plan_acceptance
        write_new(prepared_plan/'plan-acceptance.json',plan_acceptance(plan,bundle.input_sha256,
            timestamp=datetime.now(timezone.utc).isoformat(),image=settings['image']))
        # Release preparation objects before the independent disk replay. The
        # stage remains unpublished until both the plan and acceptance pass.
        del plan, gaps, mechanism, public, frozen, generated
        from genosapiens.track2_synthesis import load_plan
        print('Check the generated plan against its evidence.',flush=True)
        plan_bindings=load_plan(package,prepared_plan/'plan.json',require_review=False)[3]
        assessment_binding['mechanism_review']=None if not mechanism_record.exists() else sha256(mechanism_record.read_bytes())
        state.publish('assessment',binding=assessment_binding,result={'plan_sha256':digest})
        write_new(root/'plan-review.tsv',review_worksheet('plan',digest))
        if args.require_review:
            print('Review the private plan. Save a completed copy of plan-review.tsv.',flush=True)
            return pending_review('plan_review_required')
    if not (plan_dir/'plan-review.json').exists() and args.require_review:
        if not (root/'plan-review.tsv').exists():
            write_new(root/'plan-review.tsv',review_worksheet('plan',sha256(plan_path.read_bytes())))
        print('Plan review is required. Resume with its completed worksheet.',flush=True)
        return pending_review('plan_review_required')
    from genosapiens.track2_synthesis import load_plan
    if plan_bindings is None:
        plan_bindings=load_plan(package,plan_path,require_review=args.require_review)[3]
    binding={'image':settings['image'],
        'code':plan_bindings['code'],
        'input_sha256':bundle.input_sha256,'input_review_sha256':bundle.review_sha256,
        'automatic_input_sha256':bundle.automatic_acceptance_sha256,
        'plan_sha256':plan_bindings['plan'],
        'automatic_plan_sha256':plan_bindings['automatic_plan'],
        'plan_review_sha256':plan_bindings['plan_review']}
    completed,reason=state.completed('review',binding)
    if completed is None:
        prepared=state.restart('review',reason)
        store.inventory_path.unlink(missing_ok=True)
        print('Stage the checked construction and freeze the current runtime.',flush=True)
        stage_diagnosis_input(package,prepared/'staged',expected_input_sha256=binding['input_sha256'],
            expected_review_sha256=binding['input_review_sha256'],controlled_store=store,draft=not args.require_review)
        freeze_input_snapshot(prepared/'staged/package',plan_path,Path(settings['package'])/'envelope.json',
            prepared/'snapshot.json',runtime={'image':settings['image'],
                'code_revision':'sha256:'+settings['package_sha256']},
            controlled_store=store,draft=not args.require_review)
        state.publish('review',binding=binding,result={})
    else:
        prepared=root/'preparation'
    from genosapiens.track2_report_finalization import prepare_synthesis
    prepare_synthesis(root,settings['package'],native=False,require_review=args.require_review)
    print('The package and plan are ready for offline report synthesis.',flush=True)
    return 0


def synthesize_stage(settings, *, require_review):
    from genosapiens.track2_storage import ControlledStore
    from genosapiens.track2_synthesis import export_study_evidence
    from genosapiens.track2_report_fresh import run_fresh_report
    from genosapiens.track2_report_finalization import save_production
    # The host mounts ONLY these selected inputs and the new stage output.
    store=ControlledStore('/output','/output/controlled-work/inventory.json',settings['case_id'])
    output=Path('/output/controlled-work/report')
    def export(evidence):
        return export_study_evidence('/package','/plan/plan.json',Path(settings['package'])/'envelope.json',evidence,
        input_snapshot_path='/inputs/snapshot.json',
        expected_snapshot_sha256=sha256(Path('/inputs/snapshot.json').read_bytes()),controlled_store=store,
        public_evidence_path=Path(settings['package'])/'public-evidence',draft=not require_review,
        discovery_path='/discovery-results' if (Path(settings['package'])/'discovery-resources.json').is_file() else None,
            retained_science='/study-science' if Path('/study-science').is_dir() else None)
    result=run_fresh_report(output,evidence_binding=small_json('/inputs/production.json')['report_binding'],
        export=export,require_review=require_review,
        controls=dict(memory_bytes=int(os.environ.get('GENOSAPIENS_TRACK2_REPORT_MEMORY_BYTES',4*1024**3))))
    result_path=Path('/output/controlled-work/result.json')
    if result_path.exists():require(small_json(result_path)==result, 'The saved report result changed.')
    else:write_new(result_path,result)
    save_production('/output',small_json('/inputs/production.json')['report_binding'])
    return publish_outcome(assessment_outcome(result,require_review=require_review))


def run_operator_command(argv=None):
    phase=os.environ.get('GENOSAPIENS_TRACK2_PHASE')
    if phase is None:
        import sys
        launcher=Path(__file__).resolve().parents[2]/'track2'
        require(launcher.is_file(),'Use the delivered ./track2 launcher.')
        os.execv(str(launcher),[str(launcher),'run',*(sys.argv[1:] if argv is None else argv)])
    args=parser().parse_args(argv)
    from genosapiens.track2_operator_inputs import validate_command
    validate_command(args)
    args.require_review=args.require_review or args.review is not None or args.mechanism_review is not None
    require(args.review is None or args.resume,'Use --resume when supplying a review worksheet.')
    from genosapiens.track2_operator_install import load_installation
    settings,manifest=load_installation(os.environ['GENOSAPIENS_TRACK2_SETTINGS'])
    require(all(name=='lo' for _,name in socket.if_nameindex()),'Case execution requires an offline network namespace.')
    if phase=='synthesize':
        require(not Path(args.vcf).exists(),'Synthesis must not mount raw genomic inputs.')
        return synthesize_stage(settings,require_review=args.require_review)
    governance=authorize(settings)
    if phase=='construct':
        try:
            construct(args,settings,manifest)
            return 0
        finally:
            local=storage(settings,args,local=True)
            record_outputs(local,args.output/'construction')
            record_outputs(local,args.output/'annotation-checkpoint')
    store=storage(settings,args);output=store.path(args.output)
    private_path(output,directory=True)
    if phase=='authorize':
        for path in (args.vcf,args.hpo,args.concordance):
            if path is not None:
                # Host intake already checks these paths. This phase has no
                # raw mounts and must not resolve paths in its own filesystem.
                require(not path.is_relative_to(output),'Keep source inputs outside the run output directory.')
        return 0
    if phase=='review':
        return review(args,settings,governance)
    if phase=='prepare-report':
        from genosapiens.track2_discovery_operator import prepare_report
        return prepare_report(output,settings['package'],require_review=args.require_review)
    if phase=='inventory':
        from genosapiens.governance import load_inventory, CaseInventory, write_inventory
        records={r.relative_path:r for r in store.inventory().artifacts}
        for local in (output/'current-inventories').glob('*.json'):
            for row in load_inventory(local).artifacts:
                host=Path(row.operator_path) if row.operator_path else store.reference_root/row.relative_path
                if host.is_relative_to(output):
                    continue  # Reconcile current outputs once, below.
                previous=records.get(row.relative_path)
                require(previous is None or previous==row,'An inventoried source input changed.')
                records[row.relative_path]=row
        write_inventory(CaseInventory(store.case_id,tuple(records[k] for k in sorted(records))),store.inventory_path)
        record_outputs(store,output)
        return 0
    if phase=='finalize':
        from genosapiens.track2_report_finalization import finalize_report
        return finalize_report(output,settings['package'],require_review=args.require_review)
    if phase=='complete-report':
        from genosapiens.track2_report_assembly import complete_operator_presentation
        return complete_operator_presentation(output)
    raise OperatorError('The launcher selected an unsupported execution stage.')
