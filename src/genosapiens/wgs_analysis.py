"""Prepared-authority orchestration for bounded Track 1 analysis."""

from __future__ import annotations

from dataclasses import asdict, dataclass, field, replace
from contextlib import ExitStack
import hashlib
import json
import os
from pathlib import Path
import stat
from typing import Callable, Literal, Protocol, TYPE_CHECKING

if TYPE_CHECKING:
    from genosapiens.wgs_analysis_inputs import ResolvedCompletedAnalysisInput

from genosapiens.bounded_wgs_config import BoundedWgsExecutionConfig
from genosapiens.track1_analysis_profiles import (
    AnalysisProfileSelection, POLICY_SHA256 as ANALYSIS_PROFILE_POLICY_SHA256,
    authenticate_analysis_profile,
)
from genosapiens.track1_publication import (
    BoundedPublicationInputs,
    FinalArtifact,
    FinalPublication,
    PrimaryPublicationAuthority,
    FinalRole,
    publish_bounded_final_analysis,
)
from genosapiens.wgs_ledger import (
    Artifact,
    FilterManifest,
    JournalEntry,
    LedgerError,
    LedgerManifest,
    ParentIdentity,
    StageCommitState,
    StageRunBinding,
    WorkingSetManifest,
    canonical,
    load_model,
    publish_bytes,
    read_artifact,
    stage_lock,
)
from genosapiens.wgs_preparation import SourceAccountingClosure
from genosapiens.wgs_failures import publish_failure_evidence, static_failure
from genosapiens.wgs_handoff import HandoffPublicationAuthority
from genosapiens.wgs_read_rescue import MergedReadIntake
from genosapiens.wgs_annotation import BoundedAnnotationSource
from genosapiens.wgs_prepared_case import (
    CompletedJournalBinding,
    PreparedCasePublication,
    load_prepared_case,
)
from genosapiens.wgs_source_routing import SourceRoutingContext
from genosapiens.wgs_reporting import PreparedBoundedReports
from genosapiens.wgs_scientific_artifact import ScientificArtifact, scientific_artifact_complete
from genosapiens.wgs_scientific_products import (
    ScientificOrderProduct, ScientificPrimaryProduct,
    ScientificSelectionProduct, open_scientific_products, product_descriptors,
)

from genosapiens.wgs_runtime import RuntimePhaseRecorder, AnalysisProgressRecorder, analysis_progress, runtime_phase


ANALYSIS_RECEIPT_SCHEMA = "genosapiens.wgs-bounded-analysis-result"
ANALYSIS_RECEIPT_NAME = "controlled-analysis-result.json"
PHENOTYPE_NAME = "phenotype-ingestion.json"
PRODUCT_CHECKPOINT_NAME = "products-descriptor.json"
_MAX_PHENOTYPE_BYTES = 2 * 1024 * 1024


class WgsAnalysisError(RuntimeError):
    """Static bounded-analysis failure without controlled values or paths."""


def _digest(value: object) -> bool:
    return (
        type(value) is str
        and len(value) == 64
        and all(character in "0123456789abcdef" for character in value)
    )


def _private_directory(path: Path) -> None:
    try:
        if not path.exists():
            path.mkdir(mode=0o700, parents=True)
        metadata = path.lstat()
        if (
            not path.is_absolute()
            or path.is_symlink()
            or not stat.S_ISDIR(metadata.st_mode)
            or stat.S_IMODE(metadata.st_mode) != 0o700
            or metadata.st_uid != os.getuid()
        ):
            raise ValueError
    except (OSError, ValueError):
        raise WgsAnalysisError("private analysis directory is unsafe") from None


def _verify_phenotype(parent: PreparedCasePublication) -> tuple[Path, str]:
    try:
        entries = {value.role: value for value in parent.inventory.entries}
        expected = entries["phenotype_ingestion"].sha256
        path = parent.root.parent / "phenotype" / PHENOTYPE_NAME
        metadata = path.lstat()
        if (
            path.is_symlink()
            or not stat.S_ISREG(metadata.st_mode)
            or stat.S_IMODE(metadata.st_mode) != 0o600
            or metadata.st_uid != os.getuid()
            or metadata.st_size <= 0
            or metadata.st_size > _MAX_PHENOTYPE_BYTES
        ):
            raise ValueError
        digest = hashlib.sha256()
        with path.open("rb") as source:
            for block in iter(lambda: source.read(65536), b""):
                digest.update(block)
        if digest.hexdigest() != expected:
            raise ValueError
        return path, expected
    except (KeyError, OSError, ValueError):
        raise WgsAnalysisError("prepared phenotype artifact is invalid") from None


@dataclass(frozen=True, slots=True, repr=False)
class PreparedAnalysisContext:
    prepared: PreparedCasePublication
    source_root: Path
    source_manifest: LedgerManifest
    source_accounting: SourceAccountingClosure
    source_routing_context: SourceRoutingContext
    filter_root: Path
    filter_manifest: FilterManifest
    working_root: Path
    working_manifest: WorkingSetManifest
    phenotype_ingestion_path: Path
    phenotype_ingestion_sha256: str
    work_root: Path
    execution_config: BoundedWgsExecutionConfig
    read_intake: MergedReadIntake | None = None
    phase_recorder: RuntimePhaseRecorder | None = None
    resume: bool = False
    completed_analysis_input: ResolvedCompletedAnalysisInput | None = None
    completed_admission_input: object | None = None
    profile_selection: AnalysisProfileSelection = field(default_factory=authenticate_analysis_profile)
    profile_manifest_sha256: str | None = None
    analysis_progress: AnalysisProgressRecorder | None = None

    def with_profile_projection(self, reader):
        """Bind the authenticated gate before selected discovery or science."""
        from genosapiens.wgs_profile_projection import ProfileProjectionReader
        if (type(reader) is not ProfileProjectionReader
                or reader.root!=self.work_root/'analysis/profile'
                or reader.identity['profile']!=self.profile_selection.profile
                or reader.identity['policy_sha256']!=self.profile_selection.policy.sha256
                or reader.identity['parent']['source_manifest_sha256']!=self.source_manifest.sha256
                or self.profile_manifest_sha256 is not None
                    and self.profile_manifest_sha256!=reader.manifest_sha256):
            raise WgsAnalysisError('selected profile context differs')
        reader.require_stable()
        return replace(self,profile_manifest_sha256=reader.manifest_sha256)

    @property
    def retained_projection_path(self):
        """Original read-only projection for imported inputs; local state otherwise."""
        if self.completed_analysis_input is not None:
            return Path(self.completed_analysis_input.descriptor.retained_projection.path)
        return self.work_root / "state" / "retained.sqlite"

    @property
    def prepared_index(self):
        return self.prepared.prepared_index

    @property
    def completed_source(self):
        return self.prepared.completed_source

    def __post_init__(self) -> None:
        if type(self.profile_selection) is not AnalysisProfileSelection:
            raise WgsAnalysisError('analysis profile selection is unauthenticated')
        if self.profile_manifest_sha256 is not None and not _digest(self.profile_manifest_sha256):
            raise WgsAnalysisError('selected profile context is invalid')
        if self.completed_admission_input is not None:
            from genosapiens.wgs_admission_inputs import ResolvedCompletedAdmissionInput
            if (type(self.completed_admission_input) is not ResolvedCompletedAdmissionInput
                    or self.completed_analysis_input is None or self.read_intake is not None or self.resume
                    or self.completed_admission_input.descriptor.file('completed_analysis_input').sha256
                       != self.completed_analysis_input.descriptor.sha256):
                raise WgsAnalysisError("completed admission context is invalid")
            self.completed_admission_input.require_stable()
        if self.completed_analysis_input is not None:
            from genosapiens.wgs_analysis_inputs import ResolvedCompletedAnalysisInput
            if (type(self.completed_analysis_input) is not ResolvedCompletedAnalysisInput
                    or self.completed_analysis_input.prepared != self.prepared
                    or self.read_intake is not None):
                raise WgsAnalysisError("completed analysis context is invalid")
        if (
            type(self.prepared) is not PreparedCasePublication
            or self.source_root != self.prepared.completed_source.root
            or type(self.source_manifest) is not LedgerManifest
            or type(self.source_accounting) is not SourceAccountingClosure
            or type(self.source_routing_context) is not SourceRoutingContext
            or self.source_accounting.source_manifest_sha256
            != self.source_manifest.sha256
            or self.source_routing_context.contract.parent_source_manifest_sha256
            != self.source_manifest.sha256
            or type(self.filter_manifest) is not FilterManifest
            or type(self.working_manifest) is not WorkingSetManifest
            or self.prepared.manifest.source_manifest.sha256
            != self.source_manifest.sha256
            or self.prepared.manifest.prefilter_manifest.sha256
            != self.filter_manifest.sha256
            or self.prepared.manifest.working_sets[0].sha256
            != self.working_manifest.sha256
            or self.prepared.manifest.prepared_index.sha256
            != self.prepared_index.sha256
            or self.read_intake is not None and (
                type(self.read_intake) is not MergedReadIntake
                or self.prepared.manifest.sha256 not in {
                    row.prepared_sha256 for row in self.read_intake.manifest.source_bindings})
            or self.phase_recorder is not None
            and type(self.phase_recorder) is not RuntimePhaseRecorder
            or self.analysis_progress is not None and type(self.analysis_progress) is not AnalysisProgressRecorder
            or type(self.resume) is not bool
            or not _digest(self.phenotype_ingestion_sha256)
            or any(not value.is_absolute() for value in (
                self.source_root,
                self.filter_root,
                self.working_root,
                self.phenotype_ingestion_path,
                self.work_root,
            ))
        ):
            raise WgsAnalysisError("prepared analysis context is invalid")


class AnalysisProductBuilder(Protocol):
    def __call__(self, context: PreparedAnalysisContext) -> "BoundedAnalysisProducts": ...


@dataclass(frozen=True, slots=True, repr=False)
class BoundedAnalysisProducts:
    scientific: ScientificArtifact
    scientific_directory: Path
    order: ScientificOrderProduct
    order_directory: Path
    primary: ScientificPrimaryProduct
    selection: ScientificSelectionProduct
    selection_directory: Path
    reports: PreparedBoundedReports
    final_artifacts: tuple[FinalArtifact, ...]
    execution_config: BoundedWgsExecutionConfig
    handoff_authority: HandoffPublicationAuthority | None = None
    primary_authority: PrimaryPublicationAuthority | None = None
    read_intake: MergedReadIntake | None = None
    annotation_source: BoundedAnnotationSource | None = None
    eligibility_root: Path | None = None
    profile_root: Path | None = None

    def __post_init__(self) -> None:
        if (
            not scientific_artifact_complete(self.scientific)
            or type(self.execution_config) is not BoundedWgsExecutionConfig
            or self.annotation_source is not None and type(self.annotation_source) is not BoundedAnnotationSource
            or self.eligibility_root is not None and (not isinstance(self.eligibility_root,Path)
                or not self.eligibility_root.is_absolute() or self.eligibility_root.name!='eligibility')
            or self.profile_root is not None and (not isinstance(self.profile_root,Path)
                or not self.profile_root.is_absolute() or self.profile_root.name!='profile')
            or self.handoff_authority is not None
            and type(self.handoff_authority) is not HandoffPublicationAuthority
            or self.primary_authority is not None and type(self.primary_authority) is not PrimaryPublicationAuthority
            or type(self.order) is not ScientificOrderProduct
            or type(self.primary) is not ScientificPrimaryProduct
            or type(self.selection) is not ScientificSelectionProduct
            or type(self.reports) is not PreparedBoundedReports
            or type(self.final_artifacts) is not tuple
            or any(type(value) is not FinalArtifact for value in self.final_artifacts)
            or self.order.scientific_manifest_sha256 != self.scientific.sha256
            or self.primary.scientific_manifest_sha256 != self.scientific.sha256
            or self.primary.order_manifest_sha256 != self.order.sha256
            or self.selection.scientific_manifest_sha256 != self.scientific.sha256
            or self.selection.order_manifest_sha256 != self.order.sha256
            or any(value.evidence_manifest_sha256 != self.scientific.evidence.sha256 for value in (
                self.order, self.primary, self.selection))
            or any(value.hypothesis_count != self.scientific.hypothesis_count for value in (
                self.order, self.primary))
            or self.selection.ranked_hypothesis_count != self.scientific.hypothesis_count
            or any(not value.is_absolute() for value in (
                self.scientific_directory,
                self.order_directory,
                self.selection_directory,
            ))
        ):
            raise WgsAnalysisError("bounded analysis products are invalid")

    @property
    def publication_inputs(self) -> BoundedPublicationInputs:
        return BoundedPublicationInputs(
            self.scientific,
            self.scientific_directory,
            self.order,
            self.order_directory,
            self.primary,
            self.selection,
            self.selection_directory,
            self.reports,
            self.execution_config,
            self.read_intake,
            self.eligibility_root,
            self.profile_root,
        )


def _product_stage_binding(
    context: PreparedAnalysisContext,
    *,
    builder_identity: str,
    executable_stack_sha256: str,
    runtime_image_digest: str,
    policy_sha256: str,
    resources: tuple[str, ...],
    run_options_sha256: str,
) -> StageRunBinding:
    from genosapiens.bounded_wgs_config import bounded_wgs_execution_config_as_json

    limits = hashlib.sha256(
        bounded_wgs_execution_config_as_json(context.execution_config)
    ).hexdigest()
    parent = context.prepared.manifest
    plan = hashlib.sha256(canonical({
        "schema": "genosapiens.wgs-bounded-product-plan", "version": 2,
        "prepared_case": parent.sha256,
        "builder_identity": builder_identity,
        "analysis_profile": context.profile_selection.profile,
        "analysis_profile_policy_sha256": context.profile_selection.policy.sha256,
        "selected_profile_manifest_sha256":context.profile_manifest_sha256,
        "completed_analysis_input": (None if context.completed_analysis_input is None
            else context.completed_analysis_input.descriptor.sha256),
        "completed_admission_input": (None if context.completed_admission_input is None
            else context.completed_admission_input.binding()),
    })).hexdigest()
    return StageRunBinding(
        "bounded_products", parent.sha256, executable_stack_sha256,
        runtime_image_digest, limits, run_options_sha256,
        context.source_manifest.source.input_sha256,
        context.source_manifest.source.reference_manifest_sha256,
        policy_sha256, resources + (builder_identity,),
        (
            ParentIdentity("stage_parent", parent.schema, parent.version, parent.sha256),
            ParentIdentity("partition_plan", "genosapiens.wgs-bounded-product-plan", 2, plan),
            ParentIdentity("execution_limits", context.execution_config.schema,
                           context.execution_config.version, limits),
            ParentIdentity("executable_stack", "genosapiens.executable-stack", 1,
                           executable_stack_sha256),
            ParentIdentity("runtime_image", "genosapiens.runtime-image", 1,
                           runtime_image_digest.removeprefix("sha256:")),
        ),
    )


def _product_checkpoint_paths(context: PreparedAnalysisContext) -> tuple[Path, Path]:
    return (context.work_root / "product-checkpoint",
            context.work_root / "state" / "bounded-products.sqlite")


def _checkpoint_products(
    context: PreparedAnalysisContext, products: BoundedAnalysisProducts,
    binding: StageRunBinding, *, resume: bool = False,
) -> None:
    if products.handoff_authority is not None or products.primary_authority is not None:
        raise WgsAnalysisError("product checkpoint cannot retain handoff authority")
    root, state_path = _product_checkpoint_paths(context)
    _private_directory(root)
    _private_directory(state_path.parent)
    exists = state_path.exists() or state_path.is_symlink()
    if exists and not resume:
        raise WgsAnalysisError("product checkpoint already exists")
    with stage_lock(root, "bounded_products") as lock, StageCommitState(
            root, state_path, binding, lock=lock, create=not exists,
    ) as state:
        entries = [("primary.json", products.primary.as_bytes())]
        entries.extend((f"final-{item.role.value}.json", item.content)
                       for item in products.final_artifacts)
        if len({name for name, _ in entries}) != len(entries):
            raise WgsAnalysisError("duplicate product checkpoint role")
        if exists:
            if state.committed("products") is not None:
                raise WgsAnalysisError("product checkpoint was already committed")
            incomplete = state.open_attempt("products")
            if incomplete is None:
                raise WgsAnalysisError("product checkpoint has no recoverable attempt")
            state.archive_orphans(
                names=tuple(name for name, _ in entries) + (PRODUCT_CHECKPOINT_NAME,),
                attempt_id=incomplete.attempt_id,
            )
            state.fail(incomplete, failure_code="process_lost", interrupted=True)
        start = state.begin(partition_id="products", scientific=False,
                            max_attempts=2 if resume else 1)
        outputs = tuple(publish_bytes(root, name, content) for name, content in entries)
        descriptor = publish_bytes(root, PRODUCT_CHECKPOINT_NAME, canonical({
            "schema": "genosapiens.wgs-bounded-product-checkpoint", "version": 2,
            "artifacts": [{"name": item.name, "sha256": item.sha256,
                           "size_bytes": item.size_bytes} for item in outputs],
        }))
        state.commit_partition(start, cursor=(1,), outputs=outputs + (descriptor,),
                               counts={"artifacts": len(outputs)})
        state.finish(expected_partitions=("products",), output=descriptor,
                     process_id=start.process_id, closure_check=lambda _: True)


def _load_checkpointed_products(
    context: PreparedAnalysisContext, binding: StageRunBinding,
) -> BoundedAnalysisProducts | None:
    root, state_path = _product_checkpoint_paths(context)
    if not state_path.exists() and not state_path.is_symlink():
        if root.exists() or root.is_symlink():
            raise WgsAnalysisError("product checkpoint state is missing")
        return None
    _private_directory(root)
    try:
        with stage_lock(root, "bounded_products") as lock, StageCommitState(
                root, state_path, binding, lock=lock,
        ) as state:
            committed = state.committed("products")
            if committed is None:
                if state.open_attempt("products") is None:
                    raise WgsAnalysisError("product checkpoint is incomplete")
                return None
            tail = state.connection.execute(
                "SELECT event_bytes FROM events ORDER BY ordinal DESC LIMIT 1"
            ).fetchone()
            event = load_model(tail[0], JournalEntry)
            if event.event_type != "stage_completed":
                event = state.finish(expected_partitions=("products",),
                    output=committed.outputs[-1],
                    process_id=committed.completion.process_id,
                    closure_check=lambda _: True)
            if event.artifacts != (committed.outputs[-1],):
                raise WgsAnalysisError("product checkpoint closure differs")
            descriptor = json.loads(read_artifact(root, committed.outputs[-1]))
            if descriptor != {
                "schema": "genosapiens.wgs-bounded-product-checkpoint", "version": 2,
                "artifacts": [{"name": item.name, "sha256": item.sha256,
                               "size_bytes": item.size_bytes}
                              for item in committed.outputs[:-1]],
            }:
                raise WgsAnalysisError("product checkpoint descriptor differs")
            content = {item.name: read_artifact(root, item)
                       for item in committed.outputs[:-1]}
    except (LedgerError, OSError, ValueError, KeyError, TypeError):
        raise WgsAnalysisError("product checkpoint is inconsistent") from None

    return _products_from_checkpoint_content(context, content)


def _products_from_checkpoint_content(
    context: PreparedAnalysisContext, content: dict[str, bytes],
) -> BoundedAnalysisProducts:
    if set(content)-{'primary.json',*(f'final-{role.value}.json' for role in FinalRole)}:
        raise WgsAnalysisError('product checkpoint contains a foreign product')
    from genosapiens.wgs_reporting import (
        PreparedBoundedReports, audit_manifest_from_json, diagnostic_from_json,
    )
    from genosapiens.wgs_scientific_artifact import scientific_artifact_from_bytes, scientific_artifact_complete
    from genosapiens.wgs_scientific_products import scientific_product_from_bytes
    def final(role: FinalRole) -> bytes:
        return content[f"final-{role.value}.json"]
    def decode(value: bytes, loader):
        return loader(value, expected_sha256=hashlib.sha256(value).hexdigest())
    def product(value, kind):
        return scientific_product_from_bytes(value, expected_sha256=hashlib.sha256(value).hexdigest(), product_type=kind)
    scientific = decode(final(FinalRole.HYPOTHESES), scientific_artifact_from_bytes)
    order = product(final(FinalRole.RANKING), ScientificOrderProduct)
    primary = product(content["primary.json"], ScientificPrimaryProduct)
    selection = product(final(FinalRole.SELECTION), ScientificSelectionProduct)
    diagnostic = decode(final(FinalRole.EVIDENCE), diagnostic_from_json)
    audit = decode(final(FinalRole.AUDIT_JSON), audit_manifest_from_json)
    reports = PreparedBoundedReports(diagnostic, audit, audit.as_bytes(),
                                     audit.markdown_index(), context.work_root / "reports")
    artifacts = tuple(FinalArtifact(
        FinalRole(name.removeprefix("final-").removesuffix(".json")), value,
    ) for name, value in content.items() if name.startswith("final-"))
    if context.read_intake is None and any(item.role is FinalRole.READ_INTAKE for item in artifacts):
        raise WgsAnalysisError("product checkpoint read intake differs")
    return BoundedAnalysisProducts(
        scientific, context.work_root / "scientific",
        order, context.work_root / "order", primary,
        selection, context.work_root / "selection", reports, artifacts,
        context.execution_config, read_intake=context.read_intake,
        profile_root=context.work_root/'analysis/profile',
    )


@dataclass(frozen=True, slots=True, repr=False)
class AnalysisRunReceipt:
    status: Literal["passed", "failed", "interrupted"]
    terminal_stage: str
    base_release_sha256: str | None
    prepared_case_sha256: str
    output_manifest_sha256: str | None
    failure_code: str | None
    public_validation_authority_sha256: str | None = None
    schema: Literal["genosapiens.wgs-bounded-analysis-result"] = (
        ANALYSIS_RECEIPT_SCHEMA
    )
    version: Literal[1] = 1
    classification: Literal["controlled"] = "controlled"
    llm_policy: Literal["no-llm"] = "no-llm"

    def __post_init__(self) -> None:
        passed = self.status == "passed"
        if (
            self.schema != ANALYSIS_RECEIPT_SCHEMA
            or self.version != 1
            or self.classification != "controlled"
            or self.llm_policy != "no-llm"
            or self.status not in {"passed", "failed", "interrupted"}
            or not self.terminal_stage
            or (self.base_release_sha256 is None) == (self.public_validation_authority_sha256 is None)
            or not _digest(self.base_release_sha256 or self.public_validation_authority_sha256)
            or not _digest(self.prepared_case_sha256)
            or passed != (self.output_manifest_sha256 is not None)
            or passed == (self.failure_code is not None)
            or self.output_manifest_sha256 is not None
            and not _digest(self.output_manifest_sha256)
            or self.failure_code is not None
            and self.failure_code not in {"analysis_failed", "analysis_interrupted"}
        ):
            raise WgsAnalysisError("bounded analysis receipt is invalid")

    def as_bytes(self) -> bytes:
        return canonical(asdict(self))


@dataclass(frozen=True, slots=True, repr=False)
class BoundedAnalysisResult:
    root: Path
    receipt: AnalysisRunReceipt
    receipt_artifact: Artifact
    publication: FinalPublication | None

    def __post_init__(self) -> None:
        if (
            not self.root.is_absolute()
            or self.receipt_artifact.name != ANALYSIS_RECEIPT_NAME
            or self.receipt_artifact.sha256
            != hashlib.sha256(self.receipt.as_bytes()).hexdigest()
            or (self.receipt.status == "passed")
            != (self.publication is not None)
        ):
            raise WgsAnalysisError("bounded analysis result is invalid")


def _terminal(
    root: Path,
    *,
    status: Literal["passed", "failed", "interrupted"],
    stage: str,
    base_release_sha256: str | None,
    prepared_sha256: str,
    publication: FinalPublication | None = None,
    failure=None,
    public_validation_authority_sha256: str | None = None,
    phase_recorder: RuntimePhaseRecorder | None = None,
    progress: AnalysisProgressRecorder | None = None,
) -> BoundedAnalysisResult:
    receipt = AnalysisRunReceipt(
        status,
        stage,
        base_release_sha256,
        prepared_sha256,
        publication.output_manifest.sha256 if publication is not None else None,
        None if status == "passed" else (
            "analysis_interrupted" if status == "interrupted" else "analysis_failed"
        ),
        public_validation_authority_sha256=(None if base_release_sha256 is not None
                                            else public_validation_authority_sha256),
    )
    with runtime_phase(phase_recorder, "analysis_terminal_publication"):
        try:
            artifact = publish_bytes(root, ANALYSIS_RECEIPT_NAME, receipt.as_bytes())
        except (Exception, KeyboardInterrupt) as error:
            if failure is not None:
                failure.add_secondary(static_failure(
                    error, "final_publication", operation="receipt_publish",
                ))
                raise failure from None
            raise
        if failure is not None:
            from genosapiens.wgs_admission_diagnostics import current, digest
            owner = current()
            try:
                saved_failure = publish_failure_evidence(
                    root, stage="analysis", receipt_content=receipt.as_bytes(), failure=failure,
                )
                if saved_failure is not None and owner is not None and owner.value['primary'] is not None:
                    from genosapiens.bounded_json import read_bounded_bytes
                    content = read_bounded_bytes(root/'failure-analysis.json', maximum_bytes=16384)
                    owner.terminal_binding = (digest(receipt.as_bytes()), digest(content))
            except (Exception, KeyboardInterrupt) as error:
                failure.add_secondary(static_failure(error, 'final_publication', operation='diagnostic_publish'))
                if owner is not None:
                    owner.secondary(error, 'diagnostic_publish')
                raise failure from None
    if progress is not None and (status != 'passed' or not getattr(phase_recorder, 'analysis_inventory_pending', False)):
        try:
            progress.finish('complete' if status == 'passed' else status)
        except (Exception, KeyboardInterrupt) as error:
            if failure is None:
                raise
            failure.add_secondary(static_failure(error, 'final_publication', operation='diagnostic_publish'))
            from genosapiens.wgs_admission_diagnostics import current
            owner = current()
            if owner is not None:
                owner.secondary(error, 'progress_publish')
    return BoundedAnalysisResult(root, receipt, artifact, publication)


def analyze_bounded_case(
    *,
    output_root: str | Path,
    publication_root: str | Path,
    expected_manifest_sha256: str,
    expected_early_snapshot_sha256: str,
    original_input_path: str | Path,
    source_root: str | Path,
    source_manifest: LedgerManifest,
    source_accounting: SourceAccountingClosure,
    source_routing_context: SourceRoutingContext,
    filter_root: str | Path,
    filter_manifest: FilterManifest,
    working_root: str | Path,
    working_manifest: WorkingSetManifest,
    policy_sha256: str,
    population_manifest_sha256: str,
    clinical_manifest_sha256: str,
    exception_manifest_sha256: str,
    precision_manifest_sha256: str,
    case_manifest_sha256: str,
    selected_sample_sha256: str,
    hpo_input_sha256: str,
    phenotype_ingestion_sha256: str,
    ontology_manifest_sha256: str,
    execution_config: BoundedWgsExecutionConfig,
    base_release_bytes: bytes | None,
    expected_base_release_sha256: str | None,
    executable_stack_sha256: str,
    runtime_image_digest: str,
    journals: tuple[CompletedJournalBinding, ...],
    product_builder: AnalysisProductBuilder,
    resume: bool = False,
    public_qualification_authority_bytes: bytes | None = None,
    expected_public_qualification_authority_sha256: str | None = None,
    public_validation_authority_bytes: bytes | None = None,
    expected_public_validation_authority_sha256: str | None = None,
    read_intake: MergedReadIntake | None = None,
    phase_recorder: RuntimePhaseRecorder | None = None,
    completed_analysis_input: ResolvedCompletedAnalysisInput | None = None,
    completed_annotation_prefix=None,
    completed_admission_input=None,
    analysis_profile: str = 'focused',
    analysis_profile_policy_sha256: str = ANALYSIS_PROFILE_POLICY_SHA256,
    preparation_versions_bytes: bytes | None = None,
    report_revision: str | None = None,
) -> BoundedAnalysisResult:
    """Authenticate preparation, run bounded products, and publish final leaves.

    Public completed inputs retain producer preparation/plan/state read-only;
    the independently checked consumer owns new annotation and final publication.
    Ordinary preparation and controlled handoff keep their existing authority.
    """

    selection = authenticate_analysis_profile(analysis_profile,
        expected_policy_sha256=analysis_profile_policy_sha256)
    if base_release_bytes is not None:
        from genosapiens.track1_installation import load_run_versions as load_base_release
        release = load_base_release(base_release_bytes, expected_sha256=expected_base_release_sha256)
        release.require_analysis_profile(selection)
        work = getattr(product_builder, 'focused_work_contract', None)
        if work is not None:
            release.require_focused_work(work)
        else:
            from genosapiens.population_filter_candidate import WORK_POLICY_SHA256
            if WORK_POLICY_SHA256 in release.resource_manifest_sha256s:
                raise WgsAnalysisError('release requires active focused work construction')
    root = Path(output_root)
    if (
        not root.is_absolute()
        or not callable(product_builder)
        or type(resume) is not bool
        or phase_recorder is not None
        and type(phase_recorder) is not RuntimePhaseRecorder
    ):
        raise WgsAnalysisError("bounded analysis invocation is invalid")
    _private_directory(root)
    stage = "prepared_authority"
    prepared_sha256 = expected_manifest_sha256
    progress = None
    try:
        validation = None
        if public_validation_authority_bytes is not None:
            from genosapiens.wgs_validation_authority import load_public_validation_authority
            if (base_release_bytes is not None or expected_base_release_sha256 is not None
                    or public_qualification_authority_bytes is not None or resume):
                raise WgsAnalysisError("mixed or resumed public validation authority")
            validation = load_public_validation_authority(public_validation_authority_bytes,
                expected_sha256=expected_public_validation_authority_sha256)
            validation.require_execution(original_input_path=original_input_path,
                hpo_input_sha256=hpo_input_sha256, output_path=root,
                executable_stack_sha256=executable_stack_sha256,
                runtime_image_digest=runtime_image_digest,
                verify_input_bytes=False)
        elif expected_public_validation_authority_sha256 is not None:
            raise WgsAnalysisError("public validation authority is incomplete")
        if (
            (public_qualification_authority_bytes is None)
            != (expected_public_qualification_authority_sha256 is None)
        ):
            raise WgsAnalysisError(
                "public qualification authority binding is incomplete"
            )
        qualification = None
        if public_qualification_authority_bytes is not None:
            from genosapiens.wgs_public_qualification import (
                load_public_qualification_authority,
                require_public_qualification_execution,
            )

            qualification = load_public_qualification_authority(
                public_qualification_authority_bytes,
                expected_sha256=expected_public_qualification_authority_sha256,
            )
            require_public_qualification_execution(
                qualification,
                original_input_path=original_input_path,
                hpo_input_path=None,
                output_path=root,
                input_sha256=source_manifest.source.input_sha256,
                input_size_bytes=source_manifest.source.input_size_bytes,
                hpo_input_sha256=hpo_input_sha256,
                execution_config=execution_config,
                executable_stack_sha256=executable_stack_sha256,
                runtime_image_digest=runtime_image_digest,
                policy_sha256=policy_sha256,
                population_manifest_sha256=population_manifest_sha256,
                clinical_manifest_sha256=clinical_manifest_sha256,
                exception_manifest_sha256=exception_manifest_sha256,
                precision_manifest_sha256=precision_manifest_sha256,
            )
        with runtime_phase(phase_recorder, "prepared_authority_validation"):
            if completed_analysis_input is not None:
                from genosapiens.wgs_analysis_inputs import ResolvedCompletedAnalysisInput
                from genosapiens.track1_release import load_base_release
                if (type(completed_analysis_input) is not ResolvedCompletedAnalysisInput
                        or qualification is None or validation is not None or read_intake is not None):
                    raise WgsAnalysisError("completed input requires public qualification authority")
                current_release = load_base_release(base_release_bytes,
                    expected_sha256=expected_base_release_sha256)
                completed_analysis_input.require_consumer(qualification, current_release)
                bound_prefix = any(r.role == "completed_annotation_prefix" for r in qualification.resource_bindings)
                if bound_prefix != (completed_annotation_prefix is not None):
                    raise WgsAnalysisError("completed annotation prefix is missing or unbound")
                if completed_annotation_prefix is not None:
                    completed_annotation_prefix.require_consumer(qualification, current_release)
                bound_admission = any(r.role == "completed_admission_input" for r in qualification.resource_bindings)
                if bound_admission != (completed_admission_input is not None):
                    raise WgsAnalysisError("completed admission input is missing or unbound")
                if completed_admission_input is not None:
                    from genosapiens.wgs_admission_inputs import ResolvedCompletedAdmissionInput
                    if (type(completed_admission_input) is not ResolvedCompletedAdmissionInput
                            or completed_annotation_prefix is None or resume):
                        raise WgsAnalysisError("completed admission requires fresh public discovery")
                    completed_admission_input.require_consumer(qualification, current_release,
                        completed_analysis_input=completed_analysis_input, annotation_prefix=completed_annotation_prefix)
                prepared = completed_analysis_input.prepared
                if (prepared.manifest.sha256 != expected_manifest_sha256
                        or prepared.early_snapshot.sha256 != expected_early_snapshot_sha256
                        or prepared.root != Path(publication_root)
                        or prepared.manifest.hpo_input_sha256 != hpo_input_sha256
                        or prepared.manifest.ontology_manifest_sha256 != ontology_manifest_sha256
                        or source_manifest != prepared.completed_source.manifest
                        or source_accounting != prepared.completed_source.accounting
                        or filter_manifest != completed_analysis_input.filter_manifest
                        or working_manifest != completed_analysis_input.working_manifest
                        or journals != completed_analysis_input.journals
                        or prepared.completed_source.manifest.source.case_manifest_sha256 != case_manifest_sha256
                        or prepared.completed_source.manifest.source.selected_sample_sha256 != selected_sample_sha256
                        or Path(filter_root) != Path(completed_analysis_input.descriptor.preparation_root) / "prefilter"
                        or Path(working_root) != Path(completed_analysis_input.descriptor.preparation_root) / "working-set"):
                    raise WgsAnalysisError("completed input differs from analysis parents")
                content = completed_analysis_input.binding_bytes(annotation_prefix=completed_annotation_prefix,
                    completed_admission_input=completed_admission_input)
                path = root / "analysis-input-binding.json"
                if resume:
                    if read_artifact(root, Artifact(path.name, hashlib.sha256(content).hexdigest(), len(content))) != content:
                        raise WgsAnalysisError("completed input checkpoint binding changed")
                else:
                    publish_bytes(root, path.name, content)
            else:
                if completed_annotation_prefix is not None or completed_admission_input is not None:
                    raise WgsAnalysisError("completed annotation/admission requires completed public analysis input")
                if qualification is not None and any(row.role == "completed_analysis_input"
                        for row in qualification.resource_bindings):
                    raise WgsAnalysisError("completed analysis input is missing")
                preparation_versions = None
                if preparation_versions_bytes is not None:
                    from genosapiens.track1_installation import load_run_versions
                    preparation_versions = load_run_versions(preparation_versions_bytes,
                        expected_sha256=hashlib.sha256(preparation_versions_bytes).hexdigest())
                prepared = load_prepared_case(
                    publication_root=publication_root,
                    expected_manifest_sha256=expected_manifest_sha256,
                    expected_early_snapshot_sha256=expected_early_snapshot_sha256,
                    original_input_path=original_input_path,
                    source_root=source_root,
                    source_manifest=source_manifest,
                    source_accounting=source_accounting,
                    filter_root=filter_root,
                    filter_manifest=filter_manifest,
                    working_root=working_root,
                    working_manifest=working_manifest,
                    policy_sha256=policy_sha256,
                    scope_policy_sha256=filter_manifest.scope_policy_sha256,
                    scope_manifest_sha256=filter_manifest.scope_manifest_sha256,
                    population_manifest_sha256=population_manifest_sha256,
                    clinical_manifest_sha256=clinical_manifest_sha256,
                    exception_manifest_sha256=exception_manifest_sha256,
                    precision_manifest_sha256=precision_manifest_sha256,
                    case_manifest_sha256=case_manifest_sha256,
                    selected_sample_sha256=selected_sample_sha256,
                    hpo_input_sha256=hpo_input_sha256,
                    phenotype_ingestion_sha256=phenotype_ingestion_sha256,
                    ontology_manifest_sha256=ontology_manifest_sha256,
                    execution_config=execution_config,
                    base_release_bytes=base_release_bytes if preparation_versions is None else preparation_versions_bytes,
                    expected_base_release_sha256=expected_base_release_sha256 if preparation_versions is None else preparation_versions.sha256,
                    executable_stack_sha256=executable_stack_sha256 if preparation_versions is None else preparation_versions.executable_stack_sha256,
                    runtime_image_digest=runtime_image_digest,
                    journals=journals,
                    public_qualification_authority_bytes=(
                        public_qualification_authority_bytes
                    ),
                    expected_public_qualification_authority_sha256=(
                        expected_public_qualification_authority_sha256
                    ),
                    public_validation_authority_bytes=public_validation_authority_bytes,
                    expected_public_validation_authority_sha256=expected_public_validation_authority_sha256,
                )
            run_options_sha256 = prepared.manifest.preparation_run_options.sha256
            if (validation is not None and
                    validation.preparation_run_options_sha256 != run_options_sha256
                    or qualification is not None and
                    qualification.preparation_run_options_sha256 != run_options_sha256):
                raise WgsAnalysisError("prepared audit mode differs from public authority")
        prepared_sha256 = prepared.manifest.sha256
        with runtime_phase(phase_recorder, "analysis_resource_initialization"):
            phenotype_path, phenotype_sha256 = _verify_phenotype(prepared)
            if completed_analysis_input is not None and phenotype_sha256 != phenotype_ingestion_sha256:
                raise WgsAnalysisError("completed phenotype binding differs from analysis")
            work = root / "analysis-work"
            _private_directory(work)
        if phase_recorder is not None:
            pilot = getattr(product_builder, 'public_pilot_inputs', None)
            progress = AnalysisProgressRecorder(path=work/'state/analysis-progress.json',
                attempt_id=phase_recorder.attempt_id,
                candidate_sha256=None if qualification is None else qualification.engineering_candidate_sha256,
                restart_stage=('scientific_publication' if getattr(pilot, 'completed_scientific_inputs', None) is not None
                    else 'assessment_publication' if getattr(pilot, 'completed_assessment_inputs', None) is not None
                    else 'association_discovery' if getattr(getattr(pilot, 'completed_context_inputs', None), 'profile_only', False)
                    else 'context_assembly' if getattr(pilot, 'completed_context_inputs', None) is not None else 'analysis_entry'))
            phase_recorder.analysis_progress = progress
        context = PreparedAnalysisContext(
            prepared,
            Path(source_root),
            source_manifest,
            source_accounting,
            source_routing_context,
            Path(filter_root),
            filter_manifest,
            Path(working_root),
            working_manifest,
            phenotype_path,
            phenotype_sha256,
            work,
            execution_config,
            read_intake,
            phase_recorder,
            resume,
            completed_analysis_input,
            completed_admission_input,
            selection,
            analysis_progress=progress,
        )
        stage = "bounded_products"
        with ExitStack() as product_lifetime:
            publication_scope = None
            publication_operation = getattr(product_builder, 'publication_operation', None)
            if callable(publication_operation):
                product_builder, publication_scope = product_lifetime.enter_context(publication_operation())
            checkpoint_identity = getattr(product_builder, "checkpoint_identity", None)
            binding = None
            if callable(checkpoint_identity):
                builder_identity = checkpoint_identity(context)
                if not _digest(builder_identity):
                    raise WgsAnalysisError("analysis builder checkpoint identity is invalid")
                binding = _product_stage_binding(
                    context, builder_identity=builder_identity,
                    executable_stack_sha256=executable_stack_sha256,
                    runtime_image_digest=runtime_image_digest,
                    policy_sha256=policy_sha256,
                    resources=(population_manifest_sha256, clinical_manifest_sha256,
                               exception_manifest_sha256, precision_manifest_sha256,
                               ontology_manifest_sha256, phenotype_ingestion_sha256,
                               hpo_input_sha256, case_manifest_sha256,
                               selected_sample_sha256),
                    run_options_sha256=hashlib.sha256(canonical({
                        "preparation_run_options":
                            prepared.manifest.preparation_run_options.sha256,
                        "preparation_versions": (expected_base_release_sha256 if preparation_versions_bytes is None
                            else hashlib.sha256(preparation_versions_bytes).hexdigest()),
                        "report_revision": report_revision,
                        "qualification": expected_public_qualification_authority_sha256,
                        "validation": expected_public_validation_authority_sha256,
                        "read_intake": None if read_intake is None else read_intake.manifest.sha256,
                    })).hexdigest(),
                )
            if resume and binding is None:
                raise WgsAnalysisError("analysis builder has no checkpointed resume")
            current_run = getattr(product_builder, 'current_run', None)
            if resume and current_run is not None and current_run.reusable('science'):
                products = product_builder.reuse_current_scientific(context)
                _checkpoint_products(context, products, binding)
            else:
                products = _load_checkpointed_products(context, binding) if resume and binding else None
            if products is None:
                checkpoint_state = _product_checkpoint_paths(context)[1]
                if resume and (checkpoint_state.exists() or checkpoint_state.is_symlink()):
                    recover = getattr(product_builder, "load_completed_products", None)
                    if not callable(recover):
                        raise WgsAnalysisError("product checkpoint cannot be reconstructed")
                    products = recover(context)
                    if type(products) is not BoundedAnalysisProducts:
                        raise WgsAnalysisError("reconstructed products are invalid")
                    _checkpoint_products(context, products, binding, resume=True)
                else:
                    if resume:
                        preflight = getattr(product_builder, "resume_preflight", None)
                        if not callable(preflight):
                            raise WgsAnalysisError("analysis builder has no checkpointed resume")
                        preflight(context)
                    products = product_builder(context)
                    if type(products) is BoundedAnalysisProducts and binding is not None:
                        _checkpoint_products(context, products, binding)
            if type(products) is not BoundedAnalysisProducts:
                raise WgsAnalysisError("analysis builder returned invalid products")
            if resume and publication_scope is not None and publication_scope.reader is None:
                bind_completed = getattr(product_builder, 'bind_checkpointed_products', None)
                if not callable(bind_completed) or binding is None:
                    raise WgsAnalysisError('completed scientific publication owner is unavailable')
                bind_completed(context, products, binding)
            stage = "final_publication"
            assemble_publication = getattr(product_builder, 'assemble_primary_publication', None)
            if validation is None and callable(assemble_publication):
                # Persist scientific products before deriving invocation-specific
                # publication authority. Resume reconstructs it from current parents.
                with runtime_phase(phase_recorder, "final_publication"):
                    products = assemble_publication(context, products,
                        release_manifest_bytes=base_release_bytes,
                        public_qualification_authority=qualification)
            if validation is not None:
                reference_root = Path(validation.output_root)
                run_root = reference_root / "runs" / products.scientific.sha256
            elif qualification is None:
                run_root = (
                    root / "controlled-work" / "track1" / "runs"
                    / ("bounded-" + prepared.manifest.sha256[:20]
                       + ("-" + report_revision[:20] if report_revision is not None else ""))
                )
                reference_root = root
            else:
                reference_root = Path(qualification.output_root)
                run_root = reference_root / "qualification-final"
            analysis_progress(progress, 'publication', 'final_artifacts')
            with runtime_phase(phase_recorder, "final_publication"):
                publication = publish_bounded_final_analysis(
                    prepared,
                    products.final_artifacts,
                    bounded_inputs=products.publication_inputs,
                    reference_root=reference_root,
                    run_root=run_root,
                    base_release_bytes=base_release_bytes,
                    handoff_authority=products.handoff_authority,
                    primary_authority=products.primary_authority,
                    public_qualification_authority=qualification,
                    public_validation_authority_bytes=public_validation_authority_bytes,
                    completed_analysis_input=completed_analysis_input,
                    expected_base_release_sha256=expected_base_release_sha256,
                    progress=progress,
                    **({} if publication_scope is None or publication_scope.reader is None
                        else {'_scientific_reader': publication_scope.reader}),
                )
        return _terminal(
            root,
            status="passed",
            stage=stage,
            base_release_sha256=expected_base_release_sha256,
            prepared_sha256=prepared_sha256,
            publication=publication,
            public_validation_authority_sha256=expected_public_validation_authority_sha256,
            phase_recorder=phase_recorder, progress=progress,
        )
    except KeyboardInterrupt as error:
        return _terminal(
            root,
            status="interrupted",
            stage=stage,
            base_release_sha256=expected_base_release_sha256,
            prepared_sha256=prepared_sha256,
            failure=static_failure(
                error, stage, operation={
                    "prepared_authority": "runner_authority_authentication",
                    "bounded_products": "analysis_stage_action",
                    "final_publication": "final_publish",
                }[stage],
            ),
            public_validation_authority_sha256=expected_public_validation_authority_sha256,
            phase_recorder=phase_recorder, progress=progress,
        )
    except Exception as error:
        failure = static_failure(
            error, stage, operation={
                "prepared_authority": "runner_authority_authentication",
                "bounded_products": "analysis_stage_action",
                "final_publication": "final_publish",
            }[stage],
        )
        return _terminal(
            root,
            status=(
                "interrupted" if failure.category == "interrupted" else "failed"
            ),
            stage=stage,
            base_release_sha256=expected_base_release_sha256,
            prepared_sha256=prepared_sha256,
            failure=failure,
            public_validation_authority_sha256=expected_public_validation_authority_sha256,
            phase_recorder=phase_recorder, progress=progress,
        )
