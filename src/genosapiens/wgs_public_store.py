"""Content-addressed persistence and read-only replay of bounded public products."""

import hashlib
import json
from contextlib import contextmanager
from pathlib import Path
from tempfile import TemporaryDirectory
from types import SimpleNamespace

from genosapiens.bounded_wgs_config import bounded_wgs_execution_config_as_json, load_bounded_wgs_execution_config
from genosapiens.evidence_model import _typed, _sha
from genosapiens.track1_epcr import EpcrAssignment
from genosapiens.track1_preflight import ValidationOnlyPreflight, prepare_bounded_validation_only_submission
from genosapiens.wgs_annotation import load_fixed_annotation_source, _load_partition_closure
from genosapiens.wgs_failures import WgsFailure
from genosapiens.wgs_ledger import Artifact, JournalEntry, load_model, read_journal
from genosapiens.wgs_native_annotation import native_index_rows
from genosapiens.wgs_public_case import load_public_case
from genosapiens.wgs_public_execution import load_public_execution
from genosapiens.wgs_scientific_artifact import scientific_artifact_from_bytes
from genosapiens.wgs_scientific_products import (
    ScientificOrderProduct, ScientificPrimaryProduct, ScientificSelectionProduct,
    ScientificProducts, scientific_product_from_bytes, open_scientific_products, product_descriptors,
)
from genosapiens.wgs_shards import read_bounded_metadata
from genosapiens.wgs_validation_audit import load_validation_audit
from genosapiens.wgs_reporting import (
    BoundedDiagnosticReport, PreparedBoundedReports, audit_manifest_from_json, verify_bounded_audit,
)


def persist_public_execution(execution, *, objects):
    """Publish complete sharded dependencies before the execution completion leaf."""
    from genosapiens.track1_public_validation_coordinator import publish_validation_artifact, publish_validation_file
    objects = Path(objects)

    def metadata(content):
        publish_validation_artifact(objects / hashlib.sha256(content).hexdigest(), content)

    def artifact(root, row):
        path = Path(root) / (row.relative_path if hasattr(row, "relative_path") else row.name)
        size = row.size_bytes
        publish_validation_file(objects / row.sha256, path, expected_sha256=row.sha256, expected_size_bytes=size)

    products = execution.products
    metadata(execution.public_case.manifest.as_bytes())
    metadata(execution.public_case.authority.as_bytes())
    for source in execution.public_case.sources:
        metadata(source.prepared.completed_source.descriptor.as_bytes())
        metadata(source.prepared.prepared_index.as_bytes())
        metadata(source.prepared.inventory.as_bytes())
    metadata(bounded_wgs_execution_config_as_json(products.execution_config))
    evidence_content = read_bounded_metadata(products.scientific_directory / "scientific-evidence-manifest.json")
    evidence = products.scientific.evidence
    if evidence_content != evidence.as_bytes():
        raise WgsFailure("bounded_products", "malformed_evidence", operation="bounded_product_validate", exception_kind="project_failure")
    metadata(evidence_content)
    from genosapiens.wgs_analysis_eligibility import scientific_eligibility_component,eligibility_retained_files
    from genosapiens.wgs_association_projection import assessment_retained_files
    from genosapiens.wgs_profile_projection import scientific_profile_component, profile_retained_files
    eligibility=scientific_eligibility_component(products.scientific,products.scientific_directory,eligibility_root=getattr(products,"eligibility_root",None))
    profile=scientific_profile_component(products.scientific,products.scientific_directory,profile_root=products.profile_root)
    if eligibility is None or profile is None or dict(execution.manifest.output_bindings)['analysis_profile'] != profile['manifest_sha256']:
        raise WgsFailure('bounded_products','malformed_evidence',operation='bounded_product_validate',exception_kind='project_failure')
    from itertools import chain
    for path,relative,digest,size in chain(eligibility_retained_files(products.scientific_directory.parent,eligibility,eligibility_root=getattr(products,"eligibility_root",None)),
            assessment_retained_files(products.scientific_directory.parent,dict(products.scientific.parent_bindings)),
            profile_retained_files(products.scientific_directory.parent,profile,profile_root=products.profile_root)):
        publish_validation_file(objects/digest,path,expected_sha256=digest,expected_size_bytes=size)
    groups = ((products.scientific, products.scientific_directory, (*evidence.node_shards, *evidence.root_shards)),
        (products.order, products.order_directory, ()),
        (products.selection, products.selection_directory, ()),
        (products.reports.audit, products.reports.output_directory,
            (*products.reports.audit.json_shards, *products.reports.audit.markdown_shards)),
        (execution.audit, execution.audit_directory, execution.audit.shards))
    for manifest, root, shards in groups:
        for row in shards:
            artifact(root, row)
        metadata(manifest.as_bytes())
    for value in (products.primary, products.reports.diagnostic, execution.epcr_assignment, execution.draft.preflight):
        metadata(value.as_bytes())
    if execution.draft.csv_bytes is not None:
        metadata(execution.draft.csv_bytes)
    annotation = products.annotation_source
    if annotation.native_index is not None:
        artifact(annotation.output_root, annotation.native_index)
        for row in native_index_rows(annotation.output_root, annotation.native_index):
            for item in row["artifacts"]:
                artifact(annotation.output_root / row["directory"], SimpleNamespace(**item))
    for partition in annotation.plan.partitions:
        artifact(annotation.plan_root, partition.transport)
        artifact(annotation.plan_root, partition.expected_requests)
    metadata(annotation.plan.as_bytes())
    for row in annotation.result.manifest.partition_closures:
        content = read_bounded_metadata(annotation.output_root / row.name)
        if hashlib.sha256(content).hexdigest() != row.sha256:
            raise WgsFailure("bounded_products", "malformed_evidence", operation="bounded_product_validate", exception_kind="project_failure")
        closure = _load_partition_closure(content)
        artifact(annotation.output_root, closure.output)
        if closure.stderr is not None:
            artifact(annotation.output_root, closure.stderr)
        parsed = read_bounded_metadata(annotation.output_root / "parsed" / closure.parsed_manifest.name)
        if hashlib.sha256(parsed).hexdigest() != closure.parsed_manifest.sha256:
            raise WgsFailure("bounded_products", "malformed_evidence", operation="bounded_product_validate", exception_kind="project_failure")
        from genosapiens.wgs_annotation import _load_parsed_partition
        for shard in _load_parsed_partition(parsed).shards:
            artifact(annotation.output_root / "parsed", shard)
        metadata(parsed)
        metadata(content)
    artifact(annotation.output_root.parent / "state", annotation.result.manifest.state_artifact)
    artifact(annotation.output_root, annotation.result.journal_head.artifacts[0])
    for entry in read_journal(annotation.output_root, stage="annotation", parent_sha256=annotation.plan.sha256,
            expected_head_sha256=annotation.result.journal_head.sha256):
        metadata(entry.as_bytes())
    metadata(annotation.result.manifest.as_bytes())
    metadata(annotation.as_bytes())
    metadata(execution.as_bytes())


@contextmanager
def open_public_execution(objects, *, expected_sha256, resources, policies, scratch_root, population_selector=None):
    """Replay saved scientific integrity and diagnostics without source/VEP work."""
    from genosapiens.track1_public_validation_coordinator import publish_validation_file
    from genosapiens.wgs_public_analysis import (
        VEP_BUNDLE_SHA256, VEP_CONTRACT_PATH, _resource_session_sha256,
    )
    from genosapiens.wgs_validation_metrics import open_validation_projection
    from genosapiens.wgs_annotation import VepTransportPlan, _load_stage_manifest
    objects = Path(objects)

    def content(digest, *, max_bytes=8*1024**2):
        _sha(digest)
        value = read_bounded_metadata(objects / digest,max_bytes=max_bytes)
        if hashlib.sha256(value).hexdigest() != digest:
            raise WgsFailure("bounded_products", "malformed_evidence", operation="bounded_product_validate", exception_kind="project_failure")
        return value

    execution = load_public_execution(content(expected_sha256), expected_sha256=expected_sha256)
    bindings = dict(execution.output_bindings)
    case_content = content(execution.public_case_sha256)
    config_sha = json.loads(case_content)["manifest"]["execution_config_sha256"]
    config = load_bounded_wgs_execution_config(content(config_sha), config_sha)
    scientific = scientific_artifact_from_bytes(content(bindings["scientific"]), expected_sha256=bindings["scientific"])
    evidence = scientific.evidence
    if content(evidence.sha256) != evidence.as_bytes():
        raise WgsFailure("bounded_products", "malformed_evidence", operation="bounded_product_validate", exception_kind="project_failure")
    def product(role, kind):
        return scientific_product_from_bytes(content(bindings[role]), expected_sha256=bindings[role], product_type=kind)
    order = product("order", ScientificOrderProduct)
    primary = product("primary", ScientificPrimaryProduct)
    selection = product("selection", ScientificSelectionProduct)
    report_audit = audit_manifest_from_json(content(bindings["audit_report"]), expected_sha256=bindings["audit_report"])
    audit = load_validation_audit(content(bindings["validation_audit"]), expected_sha256=bindings["validation_audit"])
    assignment = _typed(json.loads(content(bindings["epcr"]))["result"], EpcrAssignment, decode=True)
    preflight = _typed(json.loads(content(bindings["preflight"]))["preflight"], ValidationOnlyPreflight, decode=True)
    if assignment.as_bytes() != content(bindings["epcr"]) or preflight.as_bytes() != content(bindings["preflight"]):
        raise WgsFailure("bounded_products", "malformed_evidence", operation="bounded_product_validate", exception_kind="project_failure")
    diagnostic_content = content(bindings["diagnostic"])
    diagnostic = _typed(json.loads(diagnostic_content)["report"], BoundedDiagnosticReport, decode=True)
    if diagnostic.as_bytes() != diagnostic_content:
        raise WgsFailure("bounded_products", "malformed_evidence", operation="bounded_product_validate", exception_kind="project_failure")
    scratch = Path(scratch_root)
    scratch.mkdir(mode=0o700, parents=True, exist_ok=True)
    with TemporaryDirectory(prefix="bounded-readback-", dir=scratch) as temporary:
        root = Path(temporary)
        case = load_public_case(case_content, expected_sha256=execution.public_case_sha256,
            authority_bytes=content(execution.public_authority_sha256), execution_config=config,
            scratch_path=root / "prepared.sqlite", resources=resources)
        authority = case.authority
        runtime = next((row for row in policies.runtimes if row.candidate.bindings.sha256 == execution.candidate_stack_sha256), None)
        if (runtime is None or runtime.candidate.bindings not in authority.candidate_bindings
                or authority.case_id != execution.case_id or authority.runner_view_sha256 != execution.ranker_case_view_sha256
                or tuple(row.sha256 for row in authority.inputs) != execution.runner_input_sha256s
                or authority.runner_resource_sha256s != execution.runner_resource_descriptor_sha256s
                or authority.policy_family_sha256 != policies.policy_sha256):
            raise WgsFailure("prepared_authority", "malformed_evidence", operation="prepared_authority_validate", exception_kind="project_failure")
        from genosapiens.track1_public_validation_coordinator import executable_stack_bytes
        stack = resources.executable_stack_content or executable_stack_bytes(resources.repository_root)
        if authority.executable_stack_sha256 != hashlib.sha256(stack).hexdigest():
            raise WgsFailure("prepared_authority", "malformed_evidence", operation="prepared_authority_validate", exception_kind="project_failure")
        source = case.sources[0]
        for retained in case.sources:
            if content(retained.prepared.manifest.completed_source_input.sha256) != retained.prepared.completed_source.descriptor.as_bytes():
                raise WgsFailure("prepared_authority", "malformed_evidence",
                    operation="prepared_authority_validate", exception_kind="project_failure")
        expected_parents = dict(prepared_case=source.prepared.manifest.sha256, source_manifest=source.source.sha256,
            prefilter_manifest=source.filtered.sha256, working_set=source.working.sha256,
            execution_config=config_sha, policy_bindings=execution.candidate_stack_sha256, read_intake=case.intake.manifest.sha256)
        def restore(directory, name, digest, size=None):
            target = root / directory / name
            if size is None:
                value = content(digest)
                size = len(value)
            publish_validation_file(target, objects / digest, expected_sha256=digest, expected_size_bytes=size)

        from genosapiens.wgs_analysis_eligibility import (
            load_analysis_eligibility_manifest,scientific_eligibility_component,eligibility_component_parents)
        from genosapiens.wgs_association_projection import (
            AssessmentProjection,assessment_parents,load_association_projection,MAX_DISCOVERY_MANIFEST_BYTES)
        from genosapiens.wgs_scientific import _load_manifest
        from genosapiens.wgs_profile_projection import (
            _json, _profile_component, profile_component_parents, scientific_profile_component,
            verify_analysis_profile_component, verify_profile_assessment_binding)
        from genosapiens.wgs_shards import ScientificShard
        parents=dict(scientific.parent_bindings)
        required={'analysis_eligibility_policy','analysis_eligibility_manifest','analysis_eligibility_rescue_manifest',
            'association_projection','assessment_projection','analysis_profile_manifest',
            'analysis_profile_policy','analysis_profile_catalog','analysis_profile_regulatory_index',
            'analysis_profile_adapter'}
        if parents.keys()!=expected_parents.keys()|required:
            raise WgsFailure('bounded_products','malformed_evidence',operation='bounded_product_validate',exception_kind='project_failure')
        if bindings['analysis_profile'] != parents['analysis_profile_manifest']:
            raise WgsFailure('bounded_products','malformed_evidence',operation='bounded_product_validate',exception_kind='project_failure')
        profile_content=content(bindings['analysis_profile'])
        profile_manifest=_json(profile_content)
        profile=_profile_component(profile_manifest,bindings['analysis_profile'])
        if profile_component_parents(profile) != {key:value for key,value in parents.items() if key.startswith('analysis_profile_')}:
            raise WgsFailure('bounded_products','malformed_evidence',operation='bounded_product_validate',exception_kind='project_failure')
        restore('analysis/profile','manifest.json',bindings['analysis_profile'],len(profile_content))
        for item in profile_manifest['shards']:
            shard=ScientificShard(**item)
            restore('analysis/profile',shard.relative_path,shard.sha256,shard.size_bytes)
        profile=scientific_profile_component(scientific,root/'scientific')
        verify_analysis_profile_component(root,profile,full=False,scratch_directory=root)
        expected_parents.update(profile_component_parents(profile))
        eligibility=load_analysis_eligibility_manifest(content(parents['analysis_eligibility_manifest']),parents['analysis_eligibility_manifest'])
        restore('eligibility','manifest.json',eligibility.sha256)
        restore('eligibility','rescue-manifest.json',parents['analysis_eligibility_rescue_manifest'])
        for digest,size in sorted({(digest,size) for _,_,digest,size in eligibility.dependencies}):
            restore('eligibility',f'dependency-{digest}.json',digest,size)
        for shard in (*eligibility.shards,*eligibility.witness_shards):
            restore('eligibility',shard.relative_path,shard.sha256,shard.size_bytes)
        expected_parents.update(eligibility_component_parents(scientific_eligibility_component(scientific,root/'scientific')))
        association_content=content(parents['association_projection'],max_bytes=MAX_DISCOVERY_MANIFEST_BYTES)
        associations=load_association_projection(association_content,expected_sha256=parents['association_projection'])
        assessments=_load_manifest(content(parents['assessment_projection']),expected_sha256=parents['assessment_projection'],
            result_type=AssessmentProjection,result_key='manifest')
        restore('associations','manifest.json',associations.sha256,len(association_content))
        restore('associations/assessments','manifest.json',assessments.sha256)
        for directory,shards in (('associations',(*associations.shards,*associations.subject_shards)),
                ('associations/assessments',(*assessments.evidence.node_shards,*assessments.evidence.root_shards))):
            for shard in shards: restore(directory,shard.relative_path,shard.sha256,shard.size_bytes)
        expected_parents.update(assessment_parents(associations,assessments))
        verify_profile_assessment_binding((associations,assessments),profile)
        if parents!=expected_parents:
            raise WgsFailure('bounded_products','malformed_evidence',operation='bounded_product_validate',exception_kind='project_failure')

        groups = (("scientific", (*evidence.node_shards, *evidence.root_shards)),
            ("order", ()), ("selection", ()), ("audit", audit.shards),
            ("reports", (*report_audit.json_shards, *report_audit.markdown_shards)))
        for name, shards in groups:
            (root / name).mkdir(mode=0o700, exist_ok=True)
            for shard in shards:
                restore(name, shard.relative_path, shard.sha256, shard.size_bytes)
        restore("scientific", "scientific-evidence-manifest.json", evidence.sha256)
        restore("reports", "diagnostics.json", diagnostic.sha256)
        restore("reports", "audit-manifest.json", report_audit.sha256)
        report_index = root / "reports/audit-index.md"
        report_index.write_bytes(report_audit.markdown_index())
        report_index.chmod(0o600)
        reports = PreparedBoundedReports(diagnostic, report_audit, report_audit.as_bytes(),
            report_audit.markdown_index(), root / "reports")
        verify_bounded_audit(reports)
        with open_scientific_products(scientific, ScientificProducts(order, primary, selection),
                directory=root / "scientific", expected_parent_bindings=expected_parents, scratch_directory=root) as reader:
            from genosapiens.wgs_scientific_products import scientific_report_roots
            if reports.audit.scientific_roots != scientific_report_roots(reader,
                    ScientificProducts(order, primary, selection)):
                raise WgsFailure("bounded_products", "malformed_evidence", operation="bounded_product_validate", exception_kind="project_failure")
            draft = prepare_bounded_validation_only_submission(selection, assignment, scientific, order,
                scientific_directory=root / "scientific", order_directory=root / "order", selection_directory=root / "selection",
                evaluator_source_path=resources.evaluator_source_path, expected_preflight=preflight, scientific_reader=reader)
            if (None if draft.csv_bytes is None else hashlib.sha256(draft.csv_bytes).hexdigest()) != execution.draft_csv_sha256:
                raise WgsFailure("bounded_products", "malformed_evidence", operation="bounded_product_validate", exception_kind="project_failure")
            if execution.draft_csv_sha256 is not None and content(execution.draft_csv_sha256) != draft.csv_bytes:
                raise WgsFailure("bounded_products", "malformed_evidence", operation="bounded_product_validate", exception_kind="project_failure")

            fixed_content = content(bindings["annotation_source"])
            fixed = json.loads(fixed_content)
            plan_content = content(fixed["plan_sha256"])
            plan = _typed(json.loads(plan_content)["data"], VepTransportPlan, decode=True)
            restore("annotation-plan", "manifest.json", fixed["plan_sha256"])
            for partition in plan.partitions:
                for row in (partition.transport, partition.expected_requests):
                    restore("annotation-plan", row.name, row.sha256, row.size_bytes)
            stage_content = content(fixed["stage_manifest_sha256"])
            stage = _load_stage_manifest(stage_content)
            restore("annotation", "annotation-manifest.json", fixed["stage_manifest_sha256"])
            completion = load_model(content(fixed["journal_head_sha256"]), JournalEntry)
            if (completion.event_type != "stage_completed"
                    or len(completion.artifacts) != 1
                    or completion.artifacts[0].name != "annotation-stage-closure.json"):
                raise WgsFailure("bounded_products", "malformed_evidence",
                                 operation="bounded_product_validate",
                                 exception_kind="project_failure")
            seal = completion.artifacts[0]
            restore("annotation", seal.name, seal.sha256, seal.size_bytes)
            for row in stage.partition_closures:
                restore("annotation", row.name, row.sha256, row.size_bytes)
                closure = _load_partition_closure(content(row.sha256))
                for artifact in (closure.output, closure.stderr):
                    if artifact is not None:
                        restore("annotation", artifact.name, artifact.sha256, artifact.size_bytes)
                restore("annotation/parsed", closure.parsed_manifest.name,
                        closure.parsed_manifest.sha256, closure.parsed_manifest.size_bytes)
                from genosapiens.wgs_annotation import _load_parsed_partition
                parsed = _load_parsed_partition(content(closure.parsed_manifest.sha256))
                for shard in parsed.shards:
                    restore("annotation/parsed", shard.name, shard.sha256, shard.size_bytes)
            restore("state", stage.state_artifact.name, stage.state_artifact.sha256,
                    stage.state_artifact.size_bytes)
            head = fixed["journal_head_sha256"]
            previous_ordinal = None
            while head is not None:
                entry = load_model(content(head), JournalEntry)
                if entry.stage != "annotation" or previous_ordinal is not None and entry.ordinal != previous_ordinal - 1:
                    raise WgsFailure("bounded_products", "malformed_evidence", operation="bounded_product_validate", exception_kind="project_failure")
                restore("annotation", f"annotation-{entry.ordinal:06d}.journal.json", head)
                previous_ordinal, head = entry.ordinal, entry.previous_entry_sha256
            if fixed.get("native_index") is not None:
                native_index = Artifact(**fixed["native_index"])
                if native_index.name != "native-index.jsonl":
                    raise ValueError("invalid native annotation index name")
                restore("annotation", native_index.name, native_index.sha256, native_index.size_bytes)
                for row in native_index_rows(root / "annotation", native_index):
                    for item in row["artifacts"]:
                        restore("annotation/" + row["directory"], item["name"], item["sha256"], item["size_bytes"])
            annotation = load_fixed_annotation_source(fixed_content, expected_sha256=bindings["annotation_source"],
                plan_root=root / "annotation-plan", output_root=root / "annotation")
            annotation.verify(expected_plan=plan, source_manifest=source.source, working_manifest=source.working,
                resource_session_sha256=_resource_session_sha256(resources, executable_stack_sha256=authority.executable_stack_sha256,
                    runtime_image_digest=authority.runtime_image_digest,
                    vep_contract_bytes=read_bounded_metadata(
                        resources.repository_root / VEP_CONTRACT_PATH
                    )), bundle_manifest_sha256=VEP_BUNDLE_SHA256,
                executable_stack_sha256=authority.executable_stack_sha256, runtime_image_digest=authority.runtime_image_digest,
                execution_config=config, intake_binding=case.intake.manifest.annotation_binding,
                prepared_index=source.prepared.prepared_index,
                prepared_state_root=source.prepared.root.parent / "state")
            products = SimpleNamespace(execution_config=config, scientific=scientific, scientific_directory=root / "scientific",
                order=order, order_directory=root / "order", selection=selection, selection_directory=root / "selection",
                primary=primary, reports=reports,
                annotation_source=annotation)
            with open_validation_projection(products, audit, audit_directory=root / "audit", scratch_path=root / "metrics.sqlite",
                    epcr_assignment=assignment, population_selector=population_selector, scientific_reader=reader) as recovery:
                yield SimpleNamespace(manifest=execution, public_case=case, products=products, audit=audit,
                    audit_directory=root / "audit", epcr_assignment=assignment, draft=draft,
                    metrics=recovery.metrics(), scratch_root=root, scientific_reader=reader, recovery=recovery)


def read_public_execution(objects, **arguments):
    """Return only finite verified metrics; temporary shard projections are removed."""
    with open_public_execution(objects, **arguments) as result:
        return result.manifest, result.metrics


def public_secondary_resources(execution, *, resources):
    """Reconstruct the exact primary public-resource subset from its saved audit."""
    from genosapiens.wgs_validation_audit import iter_validation_audit
    from genosapiens.wgs_public_analysis import _disease_subset
    from genosapiens.wgs_scientific import _hash
    with open_scientific_products(execution.products.scientific, product_descriptors(execution.products),
            directory=execution.products.scientific_directory,
            expected_parent_bindings=dict(execution.products.scientific.parent_bindings),
            scratch_directory=execution.products.scientific_directory.parent,
            eligibility_root=execution.products.eligibility_root,
            profile_root=execution.products.profile_root) as reader:
        for row in iter_validation_audit(execution.audit, directory=execution.audit_directory,
                execution_config=execution.products.execution_config,
                expected_parents=dict(execution.audit.parent_bindings), scientific_reader=reader):
            pass
        if reader.assessment_projections is None:
            raise WgsFailure("bounded_products", "malformed_evidence", operation="bounded_product_validate", exception_kind="project_failure")
        # The primary subset includes selected regulatory and clinical targets,
        # including targets absent from molecular transcript annotations.
        profile, sources = _disease_subset(resources, reader.assessment_projections[0].targets)
    primary = execution.products.primary
    if _hash(profile) != primary.disease_profile_sha256 or _hash(sources) != primary.disease_gene_evidence_sha256:
        raise WgsFailure("bounded_products", "malformed_evidence", operation="bounded_product_validate", exception_kind="project_failure")
    return profile, sources
