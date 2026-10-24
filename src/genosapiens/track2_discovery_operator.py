"""Native assessment and report preparation inside isolated diagnosis workers."""
from datetime import datetime,timezone
import json
from pathlib import Path

from genosapiens.track2_operator_inputs import require,small_json,write_new
from genosapiens.track2_published_genotype import canonical,sha256
from genosapiens.track2_run_state import RunState,atomic_record


def review_native(args,settings,governance,bundle,store,run):
    from genosapiens.track2_report_finalization import retain_production_routing
    retain_production_routing(args.output)
    from genosapiens.track2_construction_dependencies import python_dependencies
    from genosapiens.track2_discovery import scientific_identity
    from genosapiens.track2_discovery_resources import load_discovery_resources
    from genosapiens.track2_exploratory_planning import installed_policy,prepare_native_plan
    from genosapiens.track2_mechanism_assessment import build_assessment
    from genosapiens.track2_public_matching import load_public_evidence
    from genosapiens.track2_scientific_planning import frozen_packet
    from genosapiens.track2_synthesis import load_plan,plan_acceptance,freeze_input_snapshot
    from genosapiens.track2_operator import METHODS,review_worksheet
    from genosapiens.track2_terminal import pending_review
    root=Path(args.output);package=Path(settings['package']);state=RunState(root)
    descriptor=small_json(package/'discovery-resources.json');policy=installed_policy()
    public_root=package/'public-evidence';public=load_public_evidence(public_root)
    mechanism_record=root/'mechanism'/(bundle.input_sha256+'.json')
    if args.mechanism_review is not None:
        from genosapiens.track2_mechanism_review import read_mechanism_review
        submitted=read_mechanism_review(args.mechanism_review,bundle,public,governance=governance,host_uid=__import__('os').getuid())
        if mechanism_record.exists():
            prior=small_json(mechanism_record)
            require(prior['worksheet_sha256']==submitted['worksheet_sha256'],'An accepted mechanism review cannot change.')
        else:
            write_new(root/'reviews'/('mechanism-'+submitted['worksheet_sha256'])/'receipt.json',submitted)
            if not submitted['all_decisions_accepted']:return pending_review('mechanism_review_pending_or_rejected')
            write_new(mechanism_record,submitted)
    binding=dict(input_sha256=bundle.input_sha256,input_review_sha256=bundle.review_sha256,
        require_review=args.require_review,image=settings['image'],images=settings['images'],
        scientific_code_sha256=scientific_identity(),resources=sha256(canonical(descriptor)),
        policy=sha256(canonical(policy)),public_evidence=public['sha256'],
        mechanism_review=sha256(mechanism_record.read_bytes()) if mechanism_record.exists() else None,
        code=python_dependencies(package,[('genosapiens.track2_mechanism_assessment','build_assessment'),
            ('genosapiens.track2_exploratory_planning','prepare_exploratory_plan'),
            ('genosapiens.track2_synthesis','plan_acceptance')]))
    completed,reason=state.completed('assessment',binding)
    plan_path=root/'plan/plan.json'
    if completed is None:
        target=state.restart('assessment',reason)
        inputs=dict(execution=dict(code_sha256=sha256(canonical(binding['code'])),image=settings['image']),
            timestamp=datetime.now(timezone.utc).isoformat(),public_packet=frozen_packet(public_root) if public['summary'] else None,
            discovery=dict(resources_sha256=binding['resources'],policy_sha256=binding['policy']))
        effects=small_json(mechanism_record)['assessment']['effects'] if mechanism_record.exists() else ()
        assessment=build_assessment(bundle,effects,automatic_inputs=inputs,streamed=True)
        connection=load_discovery_resources(descriptor,'/discovery-resources',matrix_required=False)
        try:
            resources=dict(descriptor=descriptor,connection=connection,policy=policy,descriptor_path='discovery-resources.json',
                policy_path='exploratory-policy.json',runtime_images=settings['images'])
            plan,gaps=prepare_native_plan(bundle,assessment=assessment,resources=resources,methods=METHODS,
                output=target,draft=not args.require_review)
        finally:connection.close()
        from genosapiens.track2_scientific_storage import publish_json
        digest=publish_json(plan_path,plan)['sha256'];publish_json(target/'planning-gaps.json',gaps)
        write_new(target/'plan-acceptance.json',plan_acceptance(plan,bundle.input_sha256,
            timestamp=inputs['timestamp'],image=settings['image']))
        load_plan(bundle.root,plan_path,require_review=False)
        state.publish('assessment',binding=binding,result=dict(plan_sha256=digest))
        (root/'plan-review.tsv').unlink(missing_ok=True)
        write_new(root/'plan-review.tsv',review_worksheet('plan',digest))
    if args.require_review and not (root/'plan/plan-review.json').exists():
        return pending_review('plan_review_required')
    bindings=load_plan(bundle.root,plan_path,require_review=args.require_review)[3]
    review_binding=dict(bindings=bindings,images=settings['images'])
    complete,reason=state.completed('review',review_binding)
    if complete is None:
        from genosapiens.track2_storage import stage_diagnosis_input
        target=state.restart('review',reason);store.inventory_path.unlink(missing_ok=True)
        stage_diagnosis_input(bundle.root,target/'staged',expected_input_sha256=bundle.input_sha256,
            expected_review_sha256=bundle.review_sha256,controlled_store=store,draft=not args.require_review)
        freeze_input_snapshot(target/'staged/package',plan_path,None,target/'snapshot.json',
            runtime=dict(images=settings['images'],code_revision='sha256:'+settings['package_sha256']),
            controlled_store=store,draft=not args.require_review)
        state.publish('review',binding=review_binding,result={})
    atomic_record(root/'current-routing.json',dict(native=True,report_complete=False))
    return 0


def prepare_report(root,package,*,require_review=False):
    from genosapiens.track2_report_finalization import prepare_synthesis
    return prepare_synthesis(root,package,native=True,require_review=require_review)
