"""Read-only reuse of a complete original annotation stage.

Original nested producer descriptors retain their unchanged schemas. The current
consumer accepts only complete annotation and cannot execute a historical suffix.
"""
from __future__ import annotations

from dataclasses import asdict, dataclass
from contextlib import closing
import hashlib
import json
from pathlib import Path

from genosapiens.bounded_wgs_config import bounded_wgs_execution_config_as_json
from genosapiens.evidence_model import _typed
from genosapiens.wgs_analysis_inputs import (
    AnalysisInputFile, _safe_file, _path, require_upstream_dependencies, WgsAnalysisInputError,
)
from genosapiens.wgs_ledger import (MAX_MANIFEST_BYTES, Artifact, JournalEntry,
    canonical, iter_artifact_lines, load_model, read_journal_head)
from genosapiens.wgs_public_qualification import (load_public_qualification_authority,
    require_public_qualification_release)
from genosapiens.track1_release import load_base_release
from genosapiens.wgs_runtime import load_stage_runtime_evidence

SCHEMA = "genosapiens.wgs-completed-annotation-prefix"
_RESOLVED = object()


class AnnotationPrefixError(RuntimeError):
    """Static rejection of an unqualified public prefix."""


def _require(ok):
    if not ok:
        raise AnnotationPrefixError("completed annotation prefix is invalid")


def _sha(content):
    return hashlib.sha256(content).hexdigest()


@dataclass(frozen=True, slots=True, repr=False)
class CompletedAnnotationPrefix:
    annotation_root: str
    plan_sha256: str
    completed_partitions: int
    completed_requests: int
    state_size_bytes: int
    producer_authority: AnalysisInputFile
    producer_release: AnalysisInputFile
    producer_start: AnalysisInputFile
    producer_failure: AnalysisInputFile
    producer_runtime: AnalysisInputFile
    terminal_journal: AnalysisInputFile
    schema: str = SCHEMA
    version: int = 1
    classification: str = "public"
    stage_manifest: AnalysisInputFile | None = None

    def __post_init__(self):
        _path(self.annotation_root)
        _require(self.schema == SCHEMA and type(self.version) is int and self.version in (1, 2)
            and ((self.version == 1 and self.stage_manifest is None)
                 or (self.version == 2 and type(self.stage_manifest) is AnalysisInputFile
                     and Path(self.stage_manifest.path) == Path(self.annotation_root) / "annotation-manifest.json"))
            and self.classification == "public" and type(self.plan_sha256) is str
            and len(self.plan_sha256) == 64 and all(c in "0123456789abcdef" for c in self.plan_sha256)
            and all(type(v) is int and v > 0 for v in (self.completed_partitions,
                self.completed_requests, self.state_size_bytes))
            and self.completed_partitions < 65536
            and all(type(v) is AnalysisInputFile for v in (self.producer_authority,
                self.producer_release, self.producer_start, self.producer_failure,
                self.producer_runtime, self.terminal_journal)))

    def as_bytes(self):
        data = asdict(self)
        if self.version == 1:
            data.pop("stage_manifest")
        return canonical(data)

    @property
    def sha256(self):
        return _sha(self.as_bytes())

    def readonly_paths(self):
        return (self.annotation_root, str(Path(self.annotation_root).parent / "state/annotation.sqlite"),
            *(v.path for v in (self.producer_authority, self.producer_release, self.producer_start,
                self.producer_failure, self.producer_runtime, self.terminal_journal)),
            *((self.stage_manifest.path,) if self.stage_manifest is not None else ()))


def load_completed_annotation_prefix(content: bytes, *, expected_sha256: str):
    """Sole current external contract: a fully completed annotation stage."""
    value = _load_producer_annotation(content, expected_sha256=expected_sha256)
    _require(value.version == 2)
    return value


def _load_producer_annotation(content: bytes, *, expected_sha256: str):
    """Authenticate unchanged upstream schemas inside an original completed stage."""
    try:
        _require(type(content) is bytes and len(content) <= MAX_MANIFEST_BYTES and _sha(content) == expected_sha256)
        data = json.loads(content)
        if data.get("version") == 1 and "stage_manifest" not in data:
            data["stage_manifest"] = None
        value = _typed(data, CompletedAnnotationPrefix, decode=True)
        _require(value.as_bytes() == content)
        return value
    except (KeyError, TypeError, ValueError):
        raise AnnotationPrefixError("completed annotation prefix is invalid") from None


@dataclass(frozen=True, slots=True, repr=False)
class ResolvedAnnotationPrefix:
    descriptor: CompletedAnnotationPrefix
    producer: object
    producer_release: object
    resource_session_sha256: str
    closures: tuple[Artifact, ...]
    parsed_manifests: tuple[Artifact, ...]
    partition_counts: tuple[tuple[int, int, int], ...]
    _proof: object
    inherited: ResolvedAnnotationPrefix | None = None

    def __post_init__(self):
        _require(self._proof is _RESOLVED)

    @property
    def count(self):
        return self.descriptor.completed_partitions

    @property
    def root(self):
        return Path(self.descriptor.annotation_root)

    def partition_source(self, partition):
        if self.inherited is not None and partition.ordinal <= self.inherited.count:
            return self.inherited.partition_source(partition)
        return self.root, self.resource_session_sha256, self.partition_counts[partition.ordinal - 1]

    def require_consumer(self, authority, release):
        _require(self.descriptor.version == 2)
        self._require_stage_dependencies(authority, release)

    def _require_stage_dependencies(self, authority, release):
        old = self.producer
        try:
            require_upstream_dependencies(old,self.producer_release,authority,release,annotation=True)
        except (WgsAnalysisInputError,ValueError):
            raise AnnotationPrefixError("completed annotation stage dependencies differ") from None
        _require(authority.binding("completed_annotation_prefix") == self.descriptor.sha256
            and old.invocation_id == "complete-hg002-deepvariant-wgs"
            and (authority.invocation_id == "complete-hg002-deepvariant-wgs"
                 or authority.is_finite_profile_pilot)
            and authority.preparation_attempts == old.preparation_attempts == 0
            and authority.output_root != old.output_root
            and not Path(authority.output_root).is_relative_to(old.output_root)
            and not Path(old.output_root).is_relative_to(authority.output_root))
        # Replay the original resource-session identity with the current contract
        # bytes. Only the actual producer's executable participates in that old
        # identity; consumer code remains independently bound by its authority.
        from genosapiens.track1_public_validation import VEP_CONTRACT_PATH
        from genosapiens.wgs_shards import read_bounded_metadata
        contract = read_bounded_metadata(Path(__file__).resolve().parents[2] / VEP_CONTRACT_PATH)
        _require(self.resource_session_sha256 == annotation_resource_session_sha256(
            reference_manifest_sha256=self.producer_release.reference_manifest_sha256,
            vep_bundle_manifest_sha256=old.binding("vep_bundle"),
            vep_contract_sha256=_sha(contract),executable_stack_sha256=old.executable_stack_sha256,
            runtime_image_digest=old.runtime_image_digest))


def annotation_resource_session_sha256(*,reference_manifest_sha256,vep_bundle_manifest_sha256,
        vep_contract_sha256,executable_stack_sha256,runtime_image_digest):
    """Unchanged annotation resource identity shared by producer and consumer."""
    return _sha(canonical(dict(schema="genosapiens.wgs-public-annotation-resource-session",version=1,
        reference_manifest_sha256=reference_manifest_sha256,vep_bundle_manifest_sha256=vep_bundle_manifest_sha256,
        vep_contract_sha256=vep_contract_sha256,executable_stack_sha256=executable_stack_sha256,
        runtime_image_digest=runtime_image_digest)))


def resolve_completed_annotation_prefix(descriptor, *, plan, execution_config):
    _require(type(descriptor) is CompletedAnnotationPrefix and descriptor.version == 2)
    return _resolve_producer_annotation(descriptor, plan=plan, execution_config=execution_config)


def _resolve_producer_annotation(descriptor, *, plan, execution_config, _depth=0):
    """Reuse committed metadata and validate requested evidence at actual use.

    Shard payloads retain the ordinary strict ordered reader at actual use. The
    private read-only SQLite singleton/journal/index is the existing trusted-state
    boundary; this does not manufacture a completed old stage or runtime pass.
    """
    from genosapiens.vep_annotation import _readonly
    from genosapiens.wgs_annotation import (_annotation_parents, read_completed_partition_closure,
        _open_annotation_state, _annotation_state_binding)
    _require(_depth < 8 and type(descriptor) is CompletedAnnotationPrefix and descriptor.plan_sha256 == plan.sha256
        and (descriptor.completed_partitions < len(plan.partitions) if descriptor.version == 1
             else descriptor.completed_partitions == len(plan.partitions))
        and descriptor.completed_requests == sum(p.requests for p in plan.partitions[:descriptor.completed_partitions]))
    for path in descriptor.readonly_paths():
        _readonly(Path(path))
    root = Path(descriptor.annotation_root)
    _safe_file(root.parent / "state/annotation.sqlite", descriptor.state_size_bytes)
    producer = load_public_qualification_authority(descriptor.producer_authority.read(),
        expected_sha256=descriptor.producer_authority.sha256)
    release = load_base_release(descriptor.producer_release.read(), expected_sha256=descriptor.producer_release.sha256)
    require_public_qualification_release(producer, release)
    config_sha = _sha(bounded_wgs_execution_config_as_json(execution_config))
    _require(producer.invocation_id == "complete-hg002-deepvariant-wgs" and producer.preparation_attempts == 0
        and root == Path(producer.output_root) / "analysis/analysis-work/annotation"
        and producer.execution_config_sha256 == config_sha == plan.execution_limits_sha256)
    start = json.loads(descriptor.producer_start.read())
    _require(start["authority_sha256"] == producer.sha256 and start["candidate_release_sha256"] == release.sha256
        and start["candidate_sha256"] == producer.engineering_candidate_sha256
        and start["executable_stack_sha256"] == producer.executable_stack_sha256
        and start["execution_config_sha256"] == config_sha
        and start["stage"] == "analysis" and start["automatic_retries"] == 0)
    runtime = load_stage_runtime_evidence(descriptor.producer_runtime.read(), expected_stage="analysis",
        expected_attempt_id=start["runtime_attempt_id"], expected_execution_config_sha256=config_sha)
    _require(runtime.exit_code == 1 and runtime.oom_killed is False)
    failure = json.loads(descriptor.producer_failure.read())
    if descriptor.version == 2:
        _require(failure["failure_code"] == "bounded_products_malformed_evidence"
            and failure["operation"] == "annotation_replay" and failure["stage"] == "analysis")
        return _resolve_completed_stage(descriptor, plan, execution_config, producer, release, start, _depth)
    _require(failure["failure_code"] in {"annotation_child_malformed_evidence", "annotation_child_missing_response"}
        and failure["operation"] == "annotation_child_execute" and failure["stage"] == "analysis")
    terminal = load_model(descriptor.terminal_journal.read(), JournalEntry)
    _require(Path(descriptor.terminal_journal.path).parent == root
        and terminal.event_type == "attempt_failed"
        and terminal.partition_id == f"partition-{descriptor.completed_partitions + 1:06d}"
        and terminal.failure_code == failure["failure_code"])
    _require(read_journal_head(root, stage="annotation", parent_sha256=plan.sha256,
        expected_head_sha256=terminal.sha256) == terminal
        and terminal.ordinal == 2 * (descriptor.completed_partitions + 1))
    resource = next(p.sha256 for p in terminal.parents if p.role == "vep_resource_session")
    bundle = producer.binding("vep_bundle")
    parents = _annotation_parents(plan, execution_limits_sha256=config_sha,
        execution_limits_version=execution_config.version, resource_session_sha256=resource,
        bundle_manifest_sha256=bundle, executable_stack_sha256=producer.executable_stack_sha256,
        runtime_image_digest=producer.runtime_image_digest)
    _require(terminal.parents == parents and terminal.attempt_number == 1)
    # Plan location comes from the completed-analysis descriptor bound by the
    # original authority; caller supplies the already resolved exact plan root.
    from genosapiens.wgs_analysis_inputs import load_completed_analysis_input
    analysis_path = Path(start["completed_analysis_input"]["path"])
    _readonly(analysis_path)
    from genosapiens.wgs_shards import read_bounded_metadata
    analysis = load_completed_analysis_input(read_bounded_metadata(analysis_path),
        expected_sha256=producer.binding("completed_analysis_input"))
    plan_root = Path(analysis.annotation_plan_root)
    _readonly(plan_root)
    _require(analysis.annotation_plan_manifest.sha256 == plan.sha256)
    # Existing partition completions supply committed counts. No aggregate scan
    # of locations/expected/responses is needed to reopen completed partitions.
    state_binding = _annotation_state_binding(plan.sha256, resource, bundle, config_sha)
    partition_counts = tuple((p.requests,) * 3 for p in plan.partitions[:descriptor.completed_partitions])
    closures, parsed = [], []
    for partition in plan.partitions[:descriptor.completed_partitions]:
        with closing(_open_annotation_state(root, state_binding, execution_config, writable=False)) as state:
            committed = state.execute('SELECT e.event_bytes FROM commits c JOIN events e ON e.ordinal=c.ordinal WHERE c.partition_id=?',
                (f'partition-{partition.ordinal:06d}',)).fetchone()
        _require(committed is not None)
        ended = load_model(committed[0], JournalEntry)
        _require(ended.event_type == 'partition_completed' and ended.parents == parents
            and ended.partition_id == f'partition-{partition.ordinal:06d}'
            and ended.attempt_number == 1 and len(ended.artifacts) == 1)
        artifact = ended.artifacts[0]
        closure = read_completed_partition_closure(root, artifact, partition, plan=plan, plan_root=plan_root,
            resource_session_sha256=resource, bundle_manifest_sha256=bundle, execution_config=execution_config)
        closures.append(artifact)
        parsed.append(closure.parsed_manifest)
    return ResolvedAnnotationPrefix(descriptor, producer, release, resource, tuple(closures), tuple(parsed), partition_counts, _RESOLVED)


def _resolve_completed_stage(descriptor, plan, config, producer, release, start, depth):
    """Reuse a completed sequence stage, including its immutable prior prefix.

    Version two is an explicit post-annotation handoff, not an old checkpoint
    resumed under different code. Validate each producer before following its
    partition references; consumers still read and authenticate every shard.
    """
    from genosapiens.vep_annotation import _readonly
    from genosapiens.wgs_analysis_inputs import load_completed_analysis_input
    from genosapiens.wgs_annotation import (_load_stage_manifest, _annotation_parents,
        _verify_state_descriptor, _open_annotation_state, _annotation_state_binding,
        read_completed_partition_closure, read_artifact)
    from genosapiens.wgs_shards import read_bounded_metadata
    root = Path(descriptor.annotation_root)
    manifest = _load_stage_manifest(descriptor.stage_manifest.read())
    terminal = load_model(descriptor.terminal_journal.read(), JournalEntry)
    _require(Path(descriptor.terminal_journal.path).parent == root
        and terminal.event_type == "stage_completed" and len(terminal.artifacts) == 1
        and terminal.artifacts[0].name == "annotation-stage-closure.json")
    analysis_path = Path(start["completed_analysis_input"]["path"])
    _readonly(analysis_path)
    analysis = load_completed_analysis_input(read_bounded_metadata(analysis_path),
        expected_sha256=producer.binding("completed_analysis_input"))
    plan_root = Path(analysis.annotation_plan_root)
    _readonly(plan_root)
    _require(analysis.annotation_plan_manifest.sha256 == plan.sha256)
    inherited = None
    bound_prefix = next((r.sha256 for r in producer.resource_bindings
                         if r.role == "completed_annotation_prefix"), None)
    if manifest.prefix_source is not None:
        old = _load_producer_annotation(read_artifact(root, manifest.prefix_source),
            expected_sha256=manifest.prefix_source.sha256)
        _require(old.sha256 == bound_prefix and old.completed_partitions < descriptor.completed_partitions)
        inherited = _resolve_producer_annotation(old, plan=plan, execution_config=config, _depth=depth + 1)
        inherited._require_stage_dependencies(producer, release)
    else:
        _require(bound_prefix is None)
    inherited_count = 0 if inherited is None else inherited.count
    resource = manifest.resource_session_sha256
    bundle, limits = producer.binding("vep_bundle"), producer.execution_config_sha256
    parents = _annotation_parents(plan, execution_limits_sha256=limits,
        execution_limits_version=config.version, resource_session_sha256=resource,
        bundle_manifest_sha256=bundle, executable_stack_sha256=producer.executable_stack_sha256,
        runtime_image_digest=producer.runtime_image_digest, prefix_sha256=bound_prefix,
        prefix_version=1 if inherited is None else inherited.descriptor.version)
    _require(manifest.parents == parents and manifest.parent_plan_sha256 == plan.sha256
        and manifest.bundle_manifest_sha256 == bundle and manifest.execution_limits_sha256 == limits
        and manifest.sequence_requests == descriptor.completed_requests
        and manifest.sequence_partitions == descriptor.completed_partitions
        and manifest.state_artifact.size_bytes == descriptor.state_size_bytes)
    _verify_state_descriptor(root, manifest.state_artifact)
    _require(read_journal_head(root, stage="annotation", parent_sha256=plan.sha256,
        expected_head_sha256=terminal.sha256) == terminal
        and terminal.ordinal == 2 * (descriptor.completed_partitions - inherited_count) + 1
        and terminal.attempt_number is None and terminal.parents == parents)
    seal_bytes = read_artifact(root, terminal.artifacts[0])
    seal = json.loads(seal_bytes)
    _require(canonical(seal) == seal_bytes
        and seal.get("schema") == "genosapiens.wgs-annotation-stage-closure" and seal.get("version") == 1
        and seal.get("plan_sha256") == plan.sha256
        and seal.get("partition_closures") == [a.sha256 for a in manifest.partition_closures])
    partition_counts = (() if inherited is None else inherited.partition_counts) + tuple(
        (p.requests,) * 3 for p in plan.partitions[inherited_count:])
    for partition, artifact, parsed in zip(plan.partitions, manifest.partition_closures, manifest.parsed_manifests, strict=True):
        if partition.ordinal <= inherited_count:
            _require(artifact == inherited.closures[partition.ordinal - 1]
                and parsed == inherited.parsed_manifests[partition.ordinal - 1])
            continue
        closure = read_completed_partition_closure(root, artifact, partition, plan=plan, plan_root=plan_root,
            resource_session_sha256=resource, bundle_manifest_sha256=bundle, execution_config=config)
        _require(closure.parsed_manifest == parsed)
    return ResolvedAnnotationPrefix(descriptor, producer, release, resource,
        manifest.partition_closures, manifest.parsed_manifests, partition_counts, _RESOLVED, inherited)
