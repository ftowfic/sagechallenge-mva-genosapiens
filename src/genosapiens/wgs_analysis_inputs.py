"""Current-format, read-only preparation inputs for public analysis successors.

The producer reader authenticates preparation under its original release. The
consumer remains independently bound to its own authority and executable. This
module never constructs source, filter, retained state, transport or annotation.
"""

from __future__ import annotations

from contextlib import closing
from dataclasses import asdict, dataclass, fields
import hashlib
import json
import os
from pathlib import Path
import stat
from typing import Literal

from genosapiens.bounded_wgs_config import bounded_wgs_execution_config_as_json
from genosapiens.track1_release import BaseRelease, load_base_release
from genosapiens.wgs_annotation import (
    VepTransportPublication, load_completed_vep_plan, open_retained_projection,
    _projection_binding,
)
from genosapiens.wgs_ledger import (
    MAX_MANIFEST_BYTES, Artifact, FilterManifest, JournalEntry, PreparedCaseManifest,
    WorkingSetManifest, _json, canonical, load_model, read_artifact, read_journal,
)
from genosapiens.wgs_preparation import (
    COMPLETED_SOURCE_SCHEMA, CompletedSourceArtifact, PreparationRunReceipt,
    completed_source_readonly_paths, load_completed_source_input, resolve_completed_source,
)
from genosapiens.wgs_prepared_case import (
    CompletedJournalBinding, PreparedCasePublication, load_prepared_case,
)
from genosapiens.wgs_public_qualification import (
    PublicQualificationAuthority, load_public_qualification_authority,
    require_public_qualification_release,
)
from genosapiens.wgs_shards import read_bounded_metadata
from genosapiens.wgs_runtime import (
    PreparationScratchException, load_stage_runtime_evidence, load_stage_runtime_sample,
    require_capacity_qualified_runtime,
)


SCHEMA = "genosapiens.wgs-completed-analysis-input"
_RESOLVED = object()
_PREPARATION_PHASES = ("input_validation", "hpo_initialization", "source_ingestion",
    "source_accounting", "source_audit", "filter_publication", "filter_index_commit", "filter_closure",
    "working_set_publication", "working_set_closure", "preparation_terminal_publication")


class WgsAnalysisInputError(RuntimeError):
    """Static failure of the completed public input boundary."""


def _require(condition):
    if not condition:
        raise WgsAnalysisInputError("completed analysis input binding is invalid")


def _sha(content):
    return hashlib.sha256(content).hexdigest()


def _digest(value):
    return type(value) is str and len(value) == 64 and all(c in "0123456789abcdef" for c in value)


def _path(value):
    _require(type(value) is str and bool(value) and "\x00" not in value)
    path = Path(value)
    _require(path.is_absolute() and str(path) == value and ".." not in path.parts
             and not {"controlled", "controlled-work", "data"}.intersection(path.parts))
    return path


def _safe_file(path, size):
    path = _path(str(path))
    for parent in reversed(path.parents):
        _require(not parent.is_symlink() and parent.is_dir())
    metadata = path.lstat()
    _require(stat.S_ISREG(metadata.st_mode) and metadata.st_nlink == 1
        and metadata.st_uid == os.getuid() and stat.S_IMODE(metadata.st_mode) == 0o600
        and metadata.st_size == size)
    return path


@dataclass(frozen=True, slots=True, repr=False)
class AnalysisInputFile:
    path: str
    sha256: str
    size_bytes: int

    def __post_init__(self):
        _path(self.path)
        _require(_digest(self.sha256) and type(self.size_bytes) is int
                 and 0 < self.size_bytes <= MAX_MANIFEST_BYTES)

    def read(self):
        path = _safe_file(self.path, self.size_bytes)
        content = read_bounded_metadata(path)
        _require(len(content) == self.size_bytes and _sha(content) == self.sha256)
        return content


@dataclass(frozen=True, slots=True, repr=False)
class CompletedRetainedProjection:
    path: str
    user_version: int
    binding_sha256: str
    size_bytes: int
    identities: int
    occurrences: int

    def __post_init__(self):
        _path(self.path)
        _require(type(self.user_version) is int and self.user_version == 1
            and _digest(self.binding_sha256) and type(self.size_bytes) is int and self.size_bytes > 0
            and type(self.identities) is int and self.identities >= 0
            and type(self.occurrences) is int and self.occurrences >= self.identities)


@dataclass(frozen=True, slots=True, repr=False)
class CompletedAnalysisInput:
    preparation_root: str
    producer_authority: AnalysisInputFile
    producer_release: AnalysisInputFile
    preparation_receipt: AnalysisInputFile
    prepared_manifest: AnalysisInputFile
    early_snapshot: AnalysisInputFile
    annotation_plan_root: str
    annotation_plan_manifest: AnalysisInputFile
    retained_projection: CompletedRetainedProjection
    preparation_runtime: AnalysisInputFile
    preparation_runtime_sample: AnalysisInputFile
    preparation_gate_exception: AnalysisInputFile
    schema: Literal["genosapiens.wgs-completed-analysis-input"] = SCHEMA
    version: Literal[1] = 1
    classification: Literal["public"] = "public"

    def __post_init__(self):
        prep, plan = _path(self.preparation_root), _path(self.annotation_plan_root)
        from genosapiens import wgs_public_qualification as public
        public._public_path(str(prep), roots=(*public._PUBLIC_INPUT_ROOTS, *public._public_output_roots()))
        _require(self.schema == SCHEMA and type(self.version) is int and self.version == 1
                 and self.classification == "public" and prep.name == "preparation"
                 and plan == prep.parent / "analysis/analysis-work/annotation-plan"
                 and type(self.retained_projection) is CompletedRetainedProjection
                 and self.retained_projection.path == str(plan.parent / "state/retained.sqlite"))
        for name in _FILE_FIELDS:
            _require(type(getattr(self, name)) is AnalysisInputFile)
        for name, path in (
            ("preparation_receipt", prep / "controlled-preparation-result.json"),
            ("prepared_manifest", prep / "prepared/prepared-case-manifest.json"),
            ("early_snapshot", prep / "prepared/ledger-bound-early-snapshot.json"),
            ("annotation_plan_manifest", plan / "manifest.json"),
            ("preparation_runtime", prep.parent / "runtime-preparation.json"),
            ("preparation_runtime_sample", prep.parent / ".runtime-sample-preparation.json"),
        ):
            _require(getattr(self, name).path == str(path))

    def as_bytes(self):
        return canonical(asdict(self))

    @property
    def sha256(self):
        return _sha(self.as_bytes())


_FILE_FIELDS = ("producer_authority", "producer_release", "preparation_receipt",
    "prepared_manifest", "early_snapshot", "annotation_plan_manifest", "preparation_runtime",
    "preparation_runtime_sample", "preparation_gate_exception")


def load_completed_analysis_input(content, *, expected_sha256):
    """Decode bounded canonical metadata without reading any declared payload."""
    try:
        _require(_digest(expected_sha256) and _sha(content) == expected_sha256)
        value = _json(content, 64 * 1024)
        _require(type(value) is dict and set(value) == {f.name for f in fields(CompletedAnalysisInput)})
        for name in _FILE_FIELDS:
            value[name] = AnalysisInputFile(**value[name])
        value["retained_projection"] = CompletedRetainedProjection(**value["retained_projection"])
        result = CompletedAnalysisInput(**value)
        _require(result.as_bytes() == content)
        return result
    except (KeyError, ValueError, TypeError, RuntimeError):
        raise WgsAnalysisInputError("completed analysis descriptor is invalid") from None


def _source_binding(descriptor):
    manifest = load_model(descriptor.prepared_manifest.read(), PreparedCaseManifest)
    content = read_bounded_metadata(Path(manifest.completed_source_input_path))
    _require(_sha(content) == manifest.completed_source_input.sha256)
    return CompletedSourceArtifact(manifest.completed_source_input_path,
        COMPLETED_SOURCE_SCHEMA, 1, manifest.completed_source_input.sha256,
        len(content))


def completed_analysis_readonly_paths(descriptor):
    """Exact original input/mount dependencies; no annotation payload discovery."""
    _require(type(descriptor) is CompletedAnalysisInput)
    return tuple(sorted({descriptor.preparation_root, descriptor.annotation_plan_root,
        descriptor.retained_projection.path, *(getattr(descriptor, n).path for n in _FILE_FIELDS),
        *completed_source_readonly_paths(_source_binding(descriptor))}))


def _journal(root, stage, head):
    """Reuse the externally bound completed head without replaying its payloads."""
    from genosapiens.wgs_ledger import read_journal_head
    terminal=read_journal_head(root,stage=stage,expected_head_sha256=head)
    _require(terminal.event_type == "stage_completed")
    return CompletedJournalBinding(stage, root, terminal.stage_parent_sha256, head)


# Dependencies of retained preparation/transport meaning, not downstream scores.
# The plan also authenticates exact source/working-set parents, and the retained
# index binds its source routing contract and selected sample through the source.
UPSTREAM_AUTHORITY_FIELDS = (
    "input_path", "input_sha256", "input_size_bytes", "hpo_input_path",
    "hpo_input_sha256", "hpo_input_size_bytes", "runtime_image_digest",
    "execution_config_sha256", "preparation_run_options_sha256", "completed_source_input",
)
UPSTREAM_RELEASE_FIELDS = (
    "reference_identity_id", "reference_manifest_sha256", "reference_bindings_sha256",
    "runtime_image_digest",
)
PREPARATION_RESOURCE_ROLES = (
    "common_allele_policy", "common_allele_exception_manifest", "common_allele_public_index_manifest",
    "gnomad_projection_precision_report", "gnomad_projection_precision_contract",
    "gnomad_joint_manifest", "gnomad_all_sites_manifest", "clinvar_manifest", "hpo_manifest", "hpo_store",
    "reference_identity_contract", "reference_manifest",
)
ANNOTATION_RESOURCE_ROLES = (*PREPARATION_RESOURCE_ROLES, "completed_analysis_input", "vep_bundle")


def upstream_dependencies(authority, release, *, annotation=False):
    """Exact current stage dependency projection of an authenticated producer."""
    require_public_qualification_release(authority, release)
    return dict(authority={name:asdict(authority)[name] for name in UPSTREAM_AUTHORITY_FIELDS},
        reference={name:getattr(release,name) for name in UPSTREAM_RELEASE_FIELDS},
        resources={role:authority.binding(role) for role in
            (ANNOTATION_RESOURCE_ROLES if annotation else PREPARATION_RESOURCE_ROLES)})


def require_upstream_dependencies(producer, producer_release, consumer, consumer_release, *, annotation=False):
    """Each release authenticates its own code/policies; reused meaning is equal."""
    _require(upstream_dependencies(producer,producer_release,annotation=annotation)
        == upstream_dependencies(consumer,consumer_release,annotation=annotation))


def _same_upstream(producer, producer_release, consumer, consumer_release, descriptor, config):
    _require(type(consumer) is PublicQualificationAuthority and type(consumer_release) is BaseRelease
        and producer.invocation_id == "complete-hg002-deepvariant-wgs"
        and (consumer.invocation_id == "complete-hg002-deepvariant-wgs"
             or consumer.is_finite_profile_pilot)
        and producer.preparation_attempts == 1 and consumer.preparation_attempts == 0
        and consumer.binding("completed_analysis_input") == descriptor.sha256
        and not any(row.role == "completed_analysis_input" for row in producer.resource_bindings)
        and Path(producer.output_root) / "preparation" == Path(descriptor.preparation_root)
        and producer.output_root != consumer.output_root
        and not Path(consumer.output_root).is_relative_to(Path(producer.output_root))
        and not Path(producer.output_root).is_relative_to(Path(consumer.output_root)))
    require_upstream_dependencies(producer,producer_release,consumer,consumer_release)
    _require(config.version == 5 and _sha(bounded_wgs_execution_config_as_json(config))
             == consumer.execution_config_sha256)


def validate_completed_preparation_runtime(descriptor, producer, source):
    """Validate unchanged terminal telemetry and only its bound scratch exception.

    No evidence is rewritten and no preparation capacity pass is produced. The
    old analysis attempt remains historical and cannot authorize a new start.
    """
    from genosapiens.wgs_ledger import _pairs, _constant
    exception = json.loads(descriptor.preparation_gate_exception.read(),
        object_pairs_hook=_pairs, parse_constant=_constant)
    expected_fields = {"schema", "version", "classification", "status", "recorded_at_utc",
        "operator_instruction", "scope", "continuation_review_sha256", "candidate_sha256",
        "preparation_runtime_sha256", "analysis_attempt_id", "analysis_attempt_budget",
        "automatic_retries", "preparation_invocations", "preparation_capacity_pass",
        "controlled_use_authorized", "qualified_release_authorized", "historical_telemetry_rewrite_authorized"}
    _require(type(exception) is dict and set(exception) == expected_fields
        and exception["schema"] == "genosapiens.public-postprep-disk-exception"
        and type(exception["version"]) is int and exception["version"] == 1
        and exception["classification"] == "public" and exception["status"] == "authorized"
        and exception["operator_instruction"] == "exception granted to the disk-telemetry gate"
        and exception["scope"] == "preparation-scratch-coverage-only"
        and exception["candidate_sha256"] == producer.engineering_candidate_sha256
        and exception["preparation_runtime_sha256"] == descriptor.preparation_runtime.sha256
        and _digest(exception["continuation_review_sha256"])
        and type(exception["recorded_at_utc"]) is str and bool(exception["recorded_at_utc"])
        and all(type(exception[name]) is int and exception[name] == count for name, count in (
            ("analysis_attempt_budget", 1), ("automatic_retries", 0), ("preparation_invocations", 0)))
        and all(exception[name] is False for name in ("preparation_capacity_pass",
            "controlled_use_authorized", "qualified_release_authorized", "historical_telemetry_rewrite_authorized")))
    runtime_bytes = descriptor.preparation_runtime.read()
    value = _json(runtime_bytes, 65536)
    runtime = load_stage_runtime_evidence(runtime_bytes, expected_stage="preparation",
        expected_attempt_id=value["attempt_id"], expected_execution_config_sha256=producer.execution_config_sha256)
    sample = load_stage_runtime_sample(descriptor.preparation_runtime_sample.read(),
        expected_stage="preparation", expected_attempt_id=runtime.attempt_id,
        expected_execution_config_sha256=producer.execution_config_sha256, expected_exit_code=0)
    _require(runtime.sample_sha256 == descriptor.preparation_runtime_sample.sha256
             and sample.metrics == runtime.metrics)
    allowance = PreparationScratchException(descriptor.preparation_gate_exception.sha256,
        descriptor.preparation_runtime.sha256, producer.engineering_candidate_sha256,
        exception["analysis_attempt_id"])
    require_capacity_qualified_runtime(runtime, required_phases=_PREPARATION_PHASES,
        completed_source=source if load_model(descriptor.prepared_manifest.read(),
            PreparedCaseManifest).source_origin == "completed_input" else None,
        preparation_scratch_exception=allowance)
    return runtime, sample, allowance


@dataclass(frozen=True, slots=True, repr=False)
class ResolvedCompletedAnalysisInput:
    descriptor: CompletedAnalysisInput
    prepared: PreparedCasePublication
    transport: VepTransportPublication
    filter_manifest: FilterManifest
    working_manifest: WorkingSetManifest
    journals: tuple[CompletedJournalBinding, ...]
    producer_authority: PublicQualificationAuthority
    producer_release: BaseRelease
    consumer_authority: PublicQualificationAuthority
    consumer_release: BaseRelease
    _token: object

    def __post_init__(self):
        _require(self._token is _RESOLVED)

    def require_consumer(self, authority, release, *, prepared=None):
        _require(type(authority) is PublicQualificationAuthority and type(release) is BaseRelease
            and authority == self.consumer_authority and release == self.consumer_release
            and authority.binding("completed_analysis_input") == self.descriptor.sha256
            and (prepared is None or prepared == self.prepared)
            and self.prepared.manifest.base_release_sha256 == self.descriptor.producer_release.sha256
            and self.prepared.manifest.sha256 == self.descriptor.prepared_manifest.sha256)
        require_upstream_dependencies(self.producer_authority,self.producer_release,authority,release)

    def binding_bytes(self, *, annotation_prefix=None, completed_admission_input=None):
        if annotation_prefix is not None:
            annotation_prefix.require_consumer(self.consumer_authority, self.consumer_release)
        if completed_admission_input is not None:
            from genosapiens.wgs_admission_inputs import ResolvedCompletedAdmissionInput
            _require(type(completed_admission_input) is ResolvedCompletedAdmissionInput
                     and annotation_prefix is not None)
            completed_admission_input.require_consumer(self.consumer_authority, self.consumer_release,
                completed_analysis_input=self, annotation_prefix=annotation_prefix)
        return canonical({
            "schema": "genosapiens.wgs-analysis-input-binding", "version": 5,
            "producer_policy_bindings_sha256": self.producer_release.policy_bindings.sha256,
            "consumer_policy_bindings_sha256": self.consumer_release.policy_bindings.sha256,
            "upstream_dependencies": upstream_dependencies(self.producer_authority,self.producer_release),
            "classification": "public", "status": "consumed_completed_input",
            "descriptor_sha256": self.descriptor.sha256,
            "producer_authority": asdict(self.descriptor.producer_authority),
            "producer_release": asdict(self.descriptor.producer_release),
            "prepared_manifest": asdict(self.descriptor.prepared_manifest),
            "annotation_plan": asdict(self.descriptor.annotation_plan_manifest),
            "retained_projection": asdict(self.descriptor.retained_projection),
            "preparation_runtime": asdict(self.descriptor.preparation_runtime),
            "preparation_runtime_sample": asdict(self.descriptor.preparation_runtime_sample),
            "preparation_gate_exception": asdict(self.descriptor.preparation_gate_exception),
            "preparation_capacity_pass": False,
            "source_producer": asdict(self.prepared.completed_source.descriptor.producer),
            "consumer_authority_sha256": self.consumer_authority.sha256,
            "consumer_release_sha256": _sha(self.consumer_release.as_bytes()),
            "consumer_executable_stack_sha256": self.consumer_authority.executable_stack_sha256,
            "execution_config_sha256": self.consumer_authority.execution_config_sha256,
            "runtime_image_digest": self.consumer_authority.runtime_image_digest,
            "planning_operation": "reuse_completed_plan", "preparation_invocations": 0,
            "source_audit_invocations": 0,
            "completed_admission_input": None if completed_admission_input is None else completed_admission_input.binding(),
            "admission_operation": "produce_fresh" if completed_admission_input is None else "reuse_completed_admission",
            "admission_projection_invocations": int(completed_admission_input is None),
            "annotation_start_partition": 1 if annotation_prefix is None else None,
            **({"annotation_operation": "reuse_completed_annotation", "annotation_child_invocations": 0}
                if annotation_prefix is not None else {}),
            **({} if annotation_prefix is None else {
                "completed_annotation_prefix_sha256": annotation_prefix.descriptor.sha256,
                "annotation_prefix_producer_authority_sha256": annotation_prefix.producer.sha256,
                "annotation_prefix_producer_release_sha256": annotation_prefix.producer_release.sha256,
                "reused_annotation_partitions": annotation_prefix.count}),
            "controlled_use_authorized": False,
        })


def resolve_completed_analysis_input(descriptor, *, consumer_authority, consumer_release,
                                     execution_config, source_routing_context):
    """Authenticate original producer inputs and exact reused stage dependencies.

    Uses the real producer prepared/source readers, completed plan metadata and
    retained singleton state. Database pages are never copied, hashed or rebuilt.
    """
    try:
        from genosapiens.vep_annotation import _readonly
        _require(type(descriptor) is CompletedAnalysisInput)
        for path in completed_analysis_readonly_paths(descriptor):
            _readonly(Path(path))
        authority_bytes = descriptor.producer_authority.read()
        producer = load_public_qualification_authority(authority_bytes,
            expected_sha256=descriptor.producer_authority.sha256)
        release_bytes = descriptor.producer_release.read()
        release = load_base_release(release_bytes, expected_sha256=descriptor.producer_release.sha256)
        _same_upstream(producer, release, consumer_authority, consumer_release, descriptor, execution_config)
        prep = Path(descriptor.preparation_root)
        manifest = load_model(descriptor.prepared_manifest.read(), PreparedCaseManifest)
        receipt_bytes = descriptor.preparation_receipt.read()
        value = _json(receipt_bytes, 65536)
        value["completed_journal_heads"] = tuple(value["completed_journal_heads"])
        receipt = PreparationRunReceipt(**value)
        _require(receipt.as_bytes() == receipt_bytes and receipt.status == "passed"
            and receipt.audit_requested is True and receipt.public_validation_authority_sha256 is None
            and receipt.base_release_sha256 == descriptor.producer_release.sha256
            and receipt.prepared_case_sha256 == descriptor.prepared_manifest.sha256
            and receipt.early_snapshot_sha256 == descriptor.early_snapshot.sha256
            and receipt.source_origin == manifest.source_origin
            and receipt.completed_source_input_sha256 == manifest.completed_source_input.sha256
            and receipt.preparation_run_options_sha256 == producer.preparation_run_options_sha256
            and receipt.source_accounting_sha256 == manifest.source_accounting.sha256
            and receipt.source_audit_sha256 == manifest.source_audit.sha256
            and receipt.prepared_index_sha256 == manifest.prepared_index.sha256)
        source_binding = _source_binding(descriptor)
        source_input = load_completed_source_input(source_binding.read())
        source = resolve_completed_source(source_input, expected_source=source_input.producer.source,
                                           audit_requested=True)
        _require(source_input.audit_requested and source_input.audit_status == "passed"
            and (producer.completed_source_input == source_binding
                 or producer.completed_source_input is None and manifest.source_origin == "produced_this_run")
            and source.manifest.source.input_sha256 == producer.input_sha256
            and source.manifest.source.input_size_bytes == producer.input_size_bytes
            and source_routing_context.contract.parent_source_manifest_sha256 == source.manifest.sha256
            and manifest.hpo_input_sha256 == producer.hpo_input_sha256
            and manifest.ontology_manifest_sha256 == producer.binding("hpo_manifest"))
        filtered = load_model(read_bounded_metadata(prep / "prefilter/manifest.json"), FilterManifest)
        working = load_model(read_bounded_metadata(prep / "working-set/manifest.json"), WorkingSetManifest)
        journals = tuple(_journal(
            Path(source_input.source_journal_root) if stage == "source" else
            Path(source_input.source_audit_journal_root) if stage == "source_audit" else
            prep / "journals" / stage, stage, parent.sha256)
            for parent in manifest.journal_heads for stage in [parent.role.removesuffix("_journal")])
        _require(tuple(sorted(receipt.completed_journal_heads))
                 == tuple(sorted(row.sha256 for row in manifest.journal_heads)))
        phenotype_sha = _sha(read_bounded_metadata(prep / "phenotype/phenotype-ingestion.json"))
        prepared = load_prepared_case(publication_root=prep / "prepared",
            expected_manifest_sha256=descriptor.prepared_manifest.sha256,
            expected_early_snapshot_sha256=descriptor.early_snapshot.sha256,
            original_input_path=Path(producer.input_path), source_root=source.root,
            source_manifest=source.manifest, source_accounting=source.accounting,
            filter_root=prep / "prefilter", filter_manifest=filtered,
            working_root=prep / "working-set", working_manifest=working,
            policy_sha256=producer.binding("common_allele_policy"),
            scope_policy_sha256=filtered.scope_policy_sha256,
            scope_manifest_sha256=filtered.scope_manifest_sha256,
            population_manifest_sha256=producer.binding("gnomad_joint_manifest"),
            clinical_manifest_sha256=producer.binding("clinvar_manifest"),
            exception_manifest_sha256=producer.binding("common_allele_exception_manifest"),
            precision_manifest_sha256=producer.binding("gnomad_projection_precision_report"),
            case_manifest_sha256=source.manifest.source.case_manifest_sha256,
            selected_sample_sha256=source.manifest.source.selected_sample_sha256,
            hpo_input_sha256=producer.hpo_input_sha256, phenotype_ingestion_sha256=phenotype_sha,
            ontology_manifest_sha256=manifest.ontology_manifest_sha256,
            execution_config=execution_config, base_release_bytes=release_bytes,
            expected_base_release_sha256=descriptor.producer_release.sha256,
            executable_stack_sha256=producer.executable_stack_sha256,
            runtime_image_digest=producer.runtime_image_digest, journals=journals,
            public_qualification_authority_bytes=authority_bytes,
            expected_public_qualification_authority_sha256=descriptor.producer_authority.sha256)
        transport = load_completed_vep_plan(output_root=descriptor.annotation_plan_root,
            expected_manifest_sha256=descriptor.annotation_plan_manifest.sha256,
            source_manifest=source.manifest, working_manifest=working, execution_config=execution_config,
            prepared_index=prepared.prepared_index, prepared_state_root=prep / "state")
        descriptor.annotation_plan_manifest.read()
        projection = descriptor.retained_projection
        _safe_file(projection.path, projection.size_bytes)
        _require(projection.binding_sha256 == _projection_binding(prepared.prepared_index,
            source_routing_context, execution_config)
            and projection.identities == working.summary.retained_unique_identities
            and projection.occurrences == working.summary.retained_occurrences
            and transport.plan.total_retained_identities == projection.identities)
        with closing(open_retained_projection(Path(projection.path), prepared_index=prepared.prepared_index,
                routing=source_routing_context, execution_config=execution_config,
                expected_identities=projection.identities, expected_occurrences=projection.occurrences)):
            pass
        # These terminal files and the historical narrow exception are mandatory
        # inputs. Their runtime semantics are checked independently of science.
        validate_completed_preparation_runtime(descriptor, producer, source)
        return ResolvedCompletedAnalysisInput(descriptor, prepared, transport, filtered, working,
            journals, producer, release, consumer_authority, consumer_release, _RESOLVED)
    except (OSError, ValueError, TypeError, KeyError, AttributeError, RuntimeError):
        raise WgsAnalysisInputError("completed analysis inputs cannot be resolved") from None
