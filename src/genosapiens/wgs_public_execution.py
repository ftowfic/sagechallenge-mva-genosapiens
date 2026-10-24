"""One public scientific invocation over saved bounded preparation and annotation."""

from dataclasses import dataclass, replace
import hashlib
from pathlib import Path
from typing import Literal

from genosapiens.evidence_model import _sha
from genosapiens.track1_epcr import EpcrAssignment, assign_sharded_epcr
from genosapiens.track1_preflight import PreparedValidationOnlySubmission, prepare_bounded_validation_only_submission
from genosapiens.wgs_analysis import BoundedAnalysisProducts, analyze_bounded_case
from genosapiens.wgs_failures import WgsFailure
from genosapiens.wgs_ledger import JournalEntry, load_model, read_journal
from genosapiens.wgs_prepared_case import CompletedJournalBinding
from genosapiens.wgs_public_case import BoundedPublicCase, load_public_case
from genosapiens.wgs_scientific import _bytes, _load_manifest
from genosapiens.wgs_shards import read_bounded_metadata
from genosapiens.wgs_validation_audit import ValidationAuditManifest, ValidationAuditWriter


OUTPUT_ROLES = frozenset(("scientific", "order", "primary",
    "selection", "diagnostic", "audit_report", "validation_audit", "annotation_source",
    "epcr", "preflight", "analysis_profile"))


@dataclass(frozen=True, slots=True, repr=False)
class PublicExecutionManifest:
    case_id: str
    candidate_stack_sha256: str
    ranker_case_view_sha256: str
    runner_input_sha256s: tuple[str, ...]
    runner_resource_descriptor_sha256s: tuple[str, ...]
    public_case_sha256: str
    public_authority_sha256: str
    output_bindings: tuple[tuple[str, str], ...]
    draft_csv_sha256: str | None
    schema: Literal["genosapiens.wgs-public-development-execution"] = "genosapiens.wgs-public-development-execution"
    version: Literal[4] = 4
    classification: Literal["controlled"] = "controlled"
    llm_policy: Literal["no-llm"] = "no-llm"
    analysis_mode: Literal["validation_only"] = "validation_only"
    controlled_use_ready: Literal[False] = False
    competitive_output_allowed: Literal[False] = False

    def __post_init__(self):
        for value in (self.candidate_stack_sha256, self.ranker_case_view_sha256,
                self.public_case_sha256, self.public_authority_sha256, *self.runner_input_sha256s,
                *self.runner_resource_descriptor_sha256s, *(value for _, value in self.output_bindings)):
            _sha(value)
        if self.draft_csv_sha256 is not None:
            _sha(self.draft_csv_sha256)
        if (type(self.case_id) is not str or not self.case_id or not self.runner_input_sha256s
                or not self.runner_resource_descriptor_sha256s
                or self.output_bindings != tuple(sorted(self.output_bindings))
                or len(self.output_bindings) != len(OUTPUT_ROLES)
                or {key for key, _ in self.output_bindings} != OUTPUT_ROLES
                or self.schema != "genosapiens.wgs-public-development-execution"
                or type(self.version) is not int or self.version != 4
                or self.classification != "controlled" or self.llm_policy != "no-llm"
                or self.analysis_mode != "validation_only" or self.controlled_use_ready is not False
                or self.competitive_output_allowed is not False):
            raise WgsFailure("bounded_products", "malformed_evidence", operation="bounded_product_validate", exception_kind="project_failure")

    def as_bytes(self):
        content = _bytes({"classification": "controlled", "llm_policy": "no-llm", "manifest": self})
        if len(content) > 8 * 1024 * 1024:
            raise WgsFailure("bounded_products", "response_capacity", operation="bounded_product_validate", exception_kind="project_failure")
        return content

    @property
    def sha256(self):
        return hashlib.sha256(self.as_bytes()).hexdigest()


def load_public_execution(content, *, expected_sha256):
    return _load_manifest(content, expected_sha256=expected_sha256,
        result_type=PublicExecutionManifest, result_key="manifest")


@dataclass(frozen=True, slots=True, repr=False)
class BoundedPublicExecution:
    manifest: PublicExecutionManifest
    public_case: BoundedPublicCase
    products: BoundedAnalysisProducts
    audit: ValidationAuditManifest
    audit_directory: Path
    epcr_assignment: EpcrAssignment
    draft: PreparedValidationOnlySubmission

    def __post_init__(self):
        if (self.manifest.public_case_sha256 != self.public_case.manifest.sha256
                or self.manifest.public_authority_sha256 != self.public_case.authority.sha256
                or self.manifest.output_bindings != _outputs(self.products, self.audit, self.epcr_assignment, self.draft)
                or self.manifest.draft_csv_sha256 != (None if self.draft.csv_bytes is None
                    else hashlib.sha256(self.draft.csv_bytes).hexdigest())):
            raise WgsFailure("bounded_products", "malformed_evidence", operation="bounded_product_validate", exception_kind="project_failure")

    def as_bytes(self):
        return self.manifest.as_bytes()


def _outputs(products, audit, epcr, draft):
    if products.annotation_source is None:
        raise WgsFailure("bounded_products", "malformed_evidence", operation="bounded_product_validate", exception_kind="project_failure")
    objects = dict(scientific=products.scientific, order=products.order, primary=products.primary,
        selection=products.selection, diagnostic=products.reports.diagnostic, audit_report=products.reports.audit, validation_audit=audit,
        annotation_source=products.annotation_source, epcr=epcr, preflight=draft.preflight)
    from genosapiens.wgs_profile_projection import scientific_profile_component
    profile = scientific_profile_component(products.scientific, products.scientific_directory,
        profile_root=products.profile_root)
    if profile is None:
        raise WgsFailure("bounded_products", "malformed_evidence", operation="bounded_product_validate", exception_kind="project_failure")
    return tuple(sorted((*((key, hashlib.sha256(value.as_bytes()).hexdigest()) for key, value in objects.items()),
        ("analysis_profile", profile['manifest_sha256']))))


def execute_public_case(public_case, *, builder, output_root):
    """Reauthenticate saved sources, then invoke the bounded builder exactly once."""
    root = Path(output_root)
    config = public_case.sources[0].execution_config
    public_case.authority.require_execution(original_input_path=public_case.sources[0].original_input_path,
        hpo_input_sha256=public_case.sources[0].prepared.manifest.hpo_input_sha256, output_path=root,
        executable_stack_sha256=builder.executable_stack_sha256, runtime_image_digest=builder.runtime_image_digest)
    root.mkdir(mode=0o700, parents=True, exist_ok=True)
    # In-memory handles never replace independent source/intake readback.
    current = load_public_case(public_case.manifest.as_bytes(), expected_sha256=public_case.manifest.sha256,
        authority_bytes=public_case.authority.as_bytes(), execution_config=config,
        scratch_path=root / "public-case-readback.sqlite", resources=builder.resources)
    authority = current.authority
    if (builder.runtime.candidate.bindings not in authority.candidate_bindings
            or builder.policies.policy_sha256 != authority.policy_family_sha256
            or builder.executable_stack_sha256 != authority.executable_stack_sha256
            or builder.runtime_image_digest != authority.runtime_image_digest
            or builder.sample_name != current.intake.manifest.sample_id):
        raise WgsFailure("prepared_authority", "malformed_evidence", operation="prepared_authority_validate", exception_kind="project_failure")
    source = current.sources[0]
    parent = source.prepared.manifest
    root.mkdir(mode=0o700, parents=True, exist_ok=True)
    audit_root = root / "validation-audit"
    audit_root.mkdir(mode=0o700)
    writer = ValidationAuditWriter(audit_root, execution_config=config,
        prepared_sha256=parent.sha256, policy_sha256=builder.runtime.candidate.bindings.sha256)
    products = None

    def build(context):
        nonlocal products
        products = replace(builder, validation_audit=writer)(context)
        return products

    journals = []
    for head in parent.journal_heads:
        stage = head.role.removesuffix("_journal")
        directory = source.prepared.root.parent / "journals" / stage
        first = load_model(read_bounded_metadata(directory / (stage + "-000001.journal.json")), JournalEntry)
        read_journal(directory, stage=stage, parent_sha256=first.stage_parent_sha256, expected_head_sha256=head.sha256)
        journals.append(CompletedJournalBinding(stage, directory, first.stage_parent_sha256, head.sha256))
    resource_bindings = {row.role: row.sha256 for row in parent.resources}
    try:
        result = analyze_bounded_case(output_root=root / "analysis", publication_root=source.prepared.root,
            expected_manifest_sha256=parent.sha256, expected_early_snapshot_sha256=source.prepared.early_snapshot.sha256,
            original_input_path=source.original_input_path, source_root=source.source_root,
            source_manifest=source.source, source_accounting=source.source_accounting,
            source_routing_context=source.source_routing_context,
            filter_root=source.filter_root, filter_manifest=source.filtered,
            working_root=source.working_root, working_manifest=source.working,
            policy_sha256=parent.prefilter_policy_sha256,
            **{name + "_sha256": resource_bindings[name] for name in (
                "population_manifest", "clinical_manifest", "exception_manifest", "precision_manifest")},
            case_manifest_sha256=parent.case_manifest_sha256, selected_sample_sha256=parent.selected_sample_sha256,
            hpo_input_sha256=parent.hpo_input_sha256, phenotype_ingestion_sha256=source.prepared.diagnostic.phenotype_ingestion_sha256,
            ontology_manifest_sha256=parent.ontology_manifest_sha256, execution_config=config,
            base_release_bytes=None, expected_base_release_sha256=None,
            executable_stack_sha256=authority.executable_stack_sha256, runtime_image_digest=authority.runtime_image_digest,
            journals=tuple(journals), product_builder=build, public_validation_authority_bytes=authority.as_bytes(),
            expected_public_validation_authority_sha256=authority.sha256, read_intake=current.intake)
        if result.receipt.status != "passed" or products is None:
            raise WgsFailure(
                result.receipt.terminal_stage, "malformed_evidence",
                operation="analysis_stage_action", exception_kind="project_failure",
            )
        audit = writer.finish(products.scientific.sha256)
    finally:
        writer.close()
    policies = builder.policies
    epcr = assign_sharded_epcr(products.selection, policy_bytes=policies.epcr_policy_bytes,
        policy_sha256=policies.epcr_policy_sha256, derivation_bytes=policies.epcr_derivation_bytes,
        expected_derivation_sha256=policies.epcr_policy.derivation_sha256)
    draft = prepare_bounded_validation_only_submission(products.selection, epcr, products.scientific,
        products.order, scientific_directory=products.scientific_directory, order_directory=products.order_directory,
        selection_directory=products.selection_directory, evaluator_source_path=builder.resources.evaluator_source_path,
        eligibility_root=products.eligibility_root, profile_root=products.profile_root)
    manifest = PublicExecutionManifest(authority.case_id, builder.runtime.candidate.bindings.sha256,
        authority.runner_view_sha256, tuple(row.sha256 for row in authority.inputs),
        authority.runner_resource_sha256s, current.manifest.sha256, authority.sha256,
        _outputs(products, audit, epcr, draft), None if draft.csv_bytes is None else hashlib.sha256(draft.csv_bytes).hexdigest())
    return BoundedPublicExecution(manifest, current, products, audit, audit_root, epcr, draft)
