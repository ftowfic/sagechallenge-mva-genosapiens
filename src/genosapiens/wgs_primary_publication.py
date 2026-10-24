"""Assemble the ordinary bounded primary output from authenticated products."""

from dataclasses import replace
import hashlib
import json

from genosapiens.track1_epcr import assign_sharded_epcr
from genosapiens.track1_methods import (
    BoundedWgsMethodsDisclosure, Track1MethodsBindings, Track1MethodsFacts,
    build_track1_operator_packet,
)
from genosapiens.track1_preflight import (
    CompetitiveRunEvidence, prepare_bounded_release_submission,
)
from genosapiens.track1_publication import FinalArtifact, FinalRole, PrimaryPublicationAuthority
from genosapiens.track1_installation import load_run_versions as load_base_release
from genosapiens.wgs_reporting import profile_source_accounting


def _sha(content):
    return hashlib.sha256(content).hexdigest()


def assemble_primary_publication(context, products, *, builder, scientific_reader,
        release_manifest_bytes, public_qualification_authority=None):
    """Reconstruct run authority after checkpointing. never serialize a live owner."""
    if scientific_reader is None:
        raise ValueError('primary assembly requires the live authenticated scientific reader')
    release = load_base_release(release_manifest_bytes,
        expected_sha256=_sha(release_manifest_bytes))
    prepared = context.prepared
    policies = builder.policies
    assignment = assign_sharded_epcr(products.selection,
        policy_bytes=policies.epcr_policy_bytes, policy_sha256=policies.epcr_policy_sha256,
        derivation_bytes=policies.epcr_derivation_bytes,
        expected_derivation_sha256=policies.epcr_policy.derivation_sha256)
    intake = None if context.read_intake is None else context.read_intake.manifest
    rescued = intake is not None and bool(intake.deltas)
    source = context.source_manifest.source
    run = CompetitiveRunEvidence(
        run_id='track1-run:bounded-' + prepared.manifest.sha256[:20],
        input_manifest_sha256=prepared.manifest.sha256,
        early_snapshot_sha256=prepared.early_snapshot.sha256,
        primary_snapshot_sha256=products.primary.sha256,
        selection_sha256=products.selection.sha256,
        epcr_assignment_sha256=_sha(assignment.as_bytes()),
        # Preparation authenticates an accepted, complete concordance and binds
        # it to the original source. These are reference states, not diagnosis.
        accepted_reference_status='accepted',
        accepted_reference_identity_id=context.source_routing_context.contract.reference_identity_id,
        accepted_reference_manifest_sha256=source.reference_manifest_sha256,
        concordance_status=prepared.diagnostic.concordance_status,
        concordance_report_sha256=prepared.manifest.concordance_sha256,
        assessment_scope_sha256=release.assessment_scope_sha256,
        # Selection is derived only after the declared scientific domain closes;
        # preflight authenticates that closure with this same live reader.
        completely_enumerated_hypothesis_ids=tuple(sorted({
            group.primary_hypothesis_id for group in products.selection.selected})),
        technical_assessments=products.reports.diagnostic.selected_technical_assessments,
        additional_class_status='not_assessed', rescue_recommendations=(),
        rescue_execution_status='performed' if rescued else 'not_performed',
        rescue_authorization_status='authorized' if rescued else 'not_applicable',
        rescue_authorization_sha256=intake.sha256 if rescued else None,
        rescue_plan_status=('complete' if not intake.mandatory_stage_blockers else 'incomplete')
            if rescued else 'not_applicable',
        release_manifest_sha256=release.sha256, contradiction_warning_sha256s=(),
    )
    submission = prepare_bounded_release_submission(products.selection, assignment,
        run, products.scientific, products.order,
        scientific_directory=products.scientific_directory,
        order_directory=products.order_directory,
        selection_directory=products.selection_directory,
        release_manifest_bytes=release_manifest_bytes,
        expected_release_manifest_sha256=release.sha256,
        evaluator_source_path=builder.resources.evaluator_source_path,
        public_reference=builder.resources.submission_reference,
        read_intake_manifest=intake, eligibility_root=products.eligibility_root,
        profile_root=products.profile_root, scientific_reader=scientific_reader,
        public_qualification_authority=public_qualification_authority)
    policy = builder.eligibility_common_policy
    qualification = public_qualification_authority is not None
    disclosure = BoundedWgsMethodsDisclosure(
        scope_policy_sha256=context.filter_manifest.scope_policy_sha256,
        scope_manifest_sha256=context.filter_manifest.scope_manifest_sha256,
        early_policy_id=policy.policy_id, early_policy_version=str(policy.version),
        early_policy_sha256=prepared.manifest.prefilter_policy_sha256,
        early_population_source=policy.frequency.dataset, early_population_release='4.1',
        early_population_metric=policy.frequency.metric,
        early_threshold=str(policy.threshold) + ' (strictly greater, precision bound required)',
        exception_scope='Exact clinical protections and curated allele/gene intervals. authenticated policy bound above',
        wgs_coverage='Supplied calls only. source-wide accounting is in the audit. no whole-WGS runtime claim',
        unknown_retention='Unknown evidence retains its status. no benignity or no-partner inference',
        phased_bypass='In-scope phase and protected composite representations retain their population-policy bypass',
        native_bypass='In-scope native classes retain their population-policy bypass',
        downstream_population_limitation='Synthetic development separator. not empirically calibrated',
        storage_qualification='Bounded durable scientific shards with authenticated selected details',
        resource_qualification=('Public qualification only. candidate release remains unaccepted. no operator authorization'
            if qualification else 'Installed versions: ' + release.sha256),
        predecessor_effect='Early filtering changes predecessor candidates and rankings. original source accounting retained',
    )
    facts = Track1MethodsFacts(
        team_name=None, approach_name=None, submitter_registration_confirmed=None,
        repository_url=None, code_revision=builder.executable_stack_sha256,
        automation_statement='Automated bounded primary analysis, EPCR, preflight and publication',
        manual_review_statement=None, proprietary_source_statement=None,
        ai_provider_plan_data_use=None,
        runtime_and_cost='Invocation runtime phase receipts record measured costs. final publication is still in progress',
        required_acknowledgement=None, dataset_citation=None, retention_coordinator=None,
        video_status='unconfirmed',
    )
    packet = build_track1_operator_packet(facts, Track1MethodsBindings(
        run_id=run.run_id,
        csv_sha256=None if submission.csv_bytes is None else _sha(submission.csv_bytes),
        audit_json_sha256=_sha(products.reports.audit_json_index),
        audit_markdown_sha256=_sha(products.reports.audit_markdown_index),
        preflight_sha256=_sha(submission.preflight.as_bytes()),
        selected_row_count=0 if submission.csv_bytes is None else len(products.selection.selected),
        software_ready=submission.preflight.ready and not qualification,
        public_sources=tuple(sorted('Manifest SHA-256: ' + value
            for value in set(release.resource_manifest_sha256s))), license_notices=(),
    ), bounded_wgs=disclosure, analysis_profile=context.profile_selection.profile,
        profile_source_scope=profile_source_accounting(
            builder.resources.regulatory_index, builder.resources.profile_catalog))
    methods = packet.methods_markdown
    if hasattr(scientific_reader, 'scientific_scope'):
        from genosapiens.wgs_science_reporting import scope_text
        methods += '\n## Scientific coverage\n\n' + scope_text(scientific_reader.scientific_scope.header) + '\n'
    methods += '\n## Bound primary run evidence\n\n```json\n' + run.as_bytes().decode('ascii') + '\n```\n'
    if builder.focused_work_contract is not None:
        methods += '\nFocused work bindings: ' + json.dumps(builder.focused_work_contract.bindings(), sort_keys=True) + '\n'
    artifacts = (*products.final_artifacts,
        FinalArtifact(FinalRole.EPCR, assignment.as_bytes()),
        FinalArtifact(FinalRole.PREFLIGHT, submission.preflight.as_bytes()),
        FinalArtifact(FinalRole.METHODS_MARKDOWN, methods.encode('utf-8')),
        *((FinalArtifact(FinalRole.CSV, submission.csv_bytes),) if submission.csv_bytes is not None else ()))
    authority = None if not submission.preflight.ready else PrimaryPublicationAuthority(
        context.source_manifest, context.filter_manifest, context.working_manifest,
        prepared.manifest.case_manifest_sha256, prepared.manifest.concordance_sha256,
        run, builder.resources.submission_reference, builder.resources.evaluator_source_path)
    return replace(products, final_artifacts=artifacts, primary_authority=authority)
