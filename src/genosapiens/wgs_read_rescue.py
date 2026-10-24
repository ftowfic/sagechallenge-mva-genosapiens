"""Bounded read-evidence deltas tied to complete prepared source ledgers.

The producer remains separately invoked. This module authenticates its persisted
results and exposes merged intake only after replay; it never invokes a caller.
"""

from __future__ import annotations

from collections import Counter
from collections.abc import Mapping
from dataclasses import dataclass, replace
import hashlib
import json
import re
from pathlib import Path
import sqlite3
from typing import Literal

from genosapiens.evidence_model import _sha, _typed
from genosapiens.governance import ArtifactRecord, ArtifactClassification
from genosapiens.track1_authorized_evidence import (
    AddedReadInventory, AuthorizedEvidenceExpectations, BatchAllele,
    ClassAssessmentState, ClassAssessmentStatus, EvidenceUseBinding,
    ProducerStack, ReadAssessment, ReadVcfIdentityReport, ReferenceAuthorizationBinding,
    SelectedPlan, StageState, StageStatus,
    ArtifactReplay, ReferenceAuthorizationReplay, IdentityStatus, EvidenceUse,
    AssessmentOutcome, _verify_artifacts, _identity_reference, _assessment_target,
)
from genosapiens.allele_identity import ReferenceRole, SequenceAllele, MitochondrialAllele, MosaicAllele, StructuralAllele
from genosapiens.case_manifest import CaseFileClass, case_manifest_as_bytes, load_case_manifest
from genosapiens.evidence_model import AlleleEvidence, IdentityMapping, IdentityStage, SourcedRecord
from genosapiens.reference_concordance import require_rescue_reference_concordance
from genosapiens.variant_inventory import WgsClass, ScreenOutcome
from genosapiens.wgs_admission import _open_index
from genosapiens.wgs_admission import (
    NativeSourceContract, LedgerAdmissionInput, iter_ledger_admission_inputs,
    load_ledger_admission_input, _occurrence, _normalized_sequence,
    SourceTechnicalObligation, RULE_IDS,
)
from genosapiens.wgs_annotation import AnnotationIntakeBinding, RetainedAllele
from genosapiens.bounded_wgs_config import BoundedWgsExecutionConfig
from genosapiens.wgs_ledger import (
    LedgerManifest, FilterManifest, WorkingSetManifest, PreparedCaseManifest,
    JournalEntry, load_model, read_journal,
    iter_ledger, iter_filter_shard,
)
from genosapiens.wgs_prepared_case import (
    CompletedJournalBinding, PreparedCasePublication, load_prepared_case,
)
from genosapiens.wgs_preparation import SourceAccountingClosure, load_source_accounting
from genosapiens.wgs_source_routing import SourceRoutingContext
from genosapiens.wgs_failures import WgsFailure
from genosapiens.wgs_provenance import execution_sha256
from genosapiens.wgs_scientific import _bytes, _read_rows, _load_manifest, _ShardWriter, _write_manifest
from genosapiens.wgs_shards import ScientificShard, read_bounded_metadata
from genosapiens.vcf_qc import VcfQcStreamHeader


class WgsReadRescueError(WgsFailure):
    def __init__(self, message, *, category="malformed_evidence"):
        super().__init__(
            "bounded_products", category,
            operation="bounded_product_validate", exception_kind="project_failure",
        )
        RuntimeError.__init__(self, message)


@dataclass(frozen=True, slots=True, repr=False)
class IntakeSourceBinding:
    prepared_sha256: str
    early_snapshot_sha256: str
    source_manifest_sha256: str
    prefilter_manifest_sha256: str
    working_set_sha256: str
    input_sha256: str
    native_contract_sha256: str
    original_occurrences: int
    selected_sample_name: str
    source_header: VcfQcStreamHeader

    def __post_init__(self):
        for value in (self.prepared_sha256, self.early_snapshot_sha256,
                self.source_manifest_sha256, self.prefilter_manifest_sha256,
                self.working_set_sha256, self.input_sha256, self.native_contract_sha256):
            _sha(value)
        if (type(self.original_occurrences) is not int or self.original_occurrences < 0
                or type(self.selected_sample_name) is not str or not self.selected_sample_name):
            raise WgsReadRescueError("rescue original occurrence count is invalid")
        if (type(self.source_header) is not VcfQcStreamHeader
                or self.source_header.input_sha256 != self.input_sha256
                or self.selected_sample_name not in self.source_header.sample_names):
            raise WgsReadRescueError("rescue source header binding is invalid")
        self.source_header.__post_init__()


@dataclass(frozen=True, slots=True, repr=False)
class PreparedIntakeSource:
    prepared: PreparedCasePublication
    original_input_path: Path
    source_root: Path
    source: LedgerManifest
    source_accounting: SourceAccountingClosure
    filter_root: Path
    filtered: FilterManifest
    working_root: Path
    working: WorkingSetManifest
    source_routing_context: SourceRoutingContext
    execution_config: BoundedWgsExecutionConfig

    def __post_init__(self):
        if (self.prepared.manifest.source_manifest.sha256 != self.source.sha256
                or self.source_root != self.prepared.completed_source.root
                or self.prepared.manifest.prefilter_manifest.sha256 != self.filtered.sha256
                or self.prepared.manifest.working_sets[0].sha256 != self.working.sha256
                or self.source_accounting.source_manifest_sha256 != self.source.sha256
                or self.prepared.manifest.source_accounting.sha256 != self.source_accounting.sha256
                or self.source_routing_context.contract.parent_source_manifest_sha256 != self.source.sha256
                or self.source_routing_context.contract.reference_manifest_sha256
                != self.source.source.reference_manifest_sha256
                or self.prepared.manifest.execution_limits_sha256 != execution_sha256(self.execution_config)):
            raise WgsReadRescueError("prepared rescue source binding is invalid")

    @property
    def native_contract(self):
        return self.source_routing_context.contract

    @property
    def binding(self):
        from genosapiens.wgs_ledger import read_artifact
        envelope = json.loads(read_artifact(self.source_root, self.source.header))
        sample = envelope["selected_sample_name"]
        return IntakeSourceBinding(self.prepared.manifest.sha256, self.prepared.early_snapshot.sha256,
            self.source.sha256, self.filtered.sha256, self.working.sha256,
            self.source.source.input_sha256, self.native_contract.sha256,
            self.source.summary.total_physical_alternate_slots, sample,
            _typed(envelope['header'], VcfQcStreamHeader, decode=True))


def load_prepared_intake_source(*, preparation_root, original_input_path,
        expected_prepared_sha256, expected_early_snapshot_sha256, source_routing_context,
        execution_config, base_release_bytes, public_validation_authority_bytes=None):
    """Replay the same prepared authority used by bounded scientific execution."""
    root = Path(preparation_root)
    from genosapiens.wgs_preparation import load_completed_source_input
    filter_root, working_root = (root / name for name in ("prefilter", "working-set"))
    filtered = load_model(read_bounded_metadata(filter_root / "manifest.json"), FilterManifest)
    working = load_model(read_bounded_metadata(working_root / "manifest.json"), WorkingSetManifest)
    prepared = load_model(read_bounded_metadata(root / "prepared/prepared-case-manifest.json"), PreparedCaseManifest)
    if prepared.sha256 != expected_prepared_sha256:
        raise WgsReadRescueError("prepared rescue parent changed")
    source_input = load_completed_source_input(read_bounded_metadata(Path(prepared.completed_source_input_path)))
    if source_input.sha256 != prepared.completed_source_input.sha256:
        raise WgsReadRescueError("prepared source input changed")
    source_root = Path(source_input.source_root)
    source = load_model(source_input.source_manifest.read(), LedgerManifest)
    closure = load_source_accounting(read_bounded_metadata(source_root / "accounting.json"))
    if closure.sha256 != prepared.source_accounting.sha256:
        raise WgsReadRescueError("prepared rescue accounting binding changed")
    journals = []
    for head in prepared.journal_heads:
        stage = head.role.removesuffix("_journal")
        directory = (Path(source_input.source_journal_root) if stage == "source" else
                     Path(source_input.source_audit_journal_root) if stage == "source_audit" else
                     root / "journals" / stage)
        paths = sorted(directory.glob(stage + "-*.journal.json"))
        if not paths:
            raise WgsReadRescueError("prepared rescue journal is missing")
        first = load_model(read_bounded_metadata(paths[0]), JournalEntry)
        entries = read_journal(directory, stage=stage,
            parent_sha256=first.stage_parent_sha256, expected_head_sha256=head.sha256)
        if not entries or entries[-1].event_type != "stage_completed":
            raise WgsReadRescueError("prepared rescue journal is incomplete")
        journals.append(CompletedJournalBinding(stage, directory, first.stage_parent_sha256, head.sha256))
    resources = {row.role: row.sha256 for row in prepared.resources}
    phenotype = root / "phenotype/phenotype-ingestion.json"
    publication = load_prepared_case(publication_root=root / "prepared",
        expected_manifest_sha256=expected_prepared_sha256,
        expected_early_snapshot_sha256=expected_early_snapshot_sha256,
        original_input_path=original_input_path, source_root=source_root, source_manifest=source,
        source_accounting=closure,
        filter_root=filter_root, filter_manifest=filtered,
        working_root=working_root, working_manifest=working,
        policy_sha256=prepared.prefilter_policy_sha256,
        scope_policy_sha256=filtered.scope_policy_sha256,
        scope_manifest_sha256=filtered.scope_manifest_sha256,
        **{name + "_sha256": resources[name] for name in (
            "population_manifest", "clinical_manifest", "exception_manifest", "precision_manifest")},
        case_manifest_sha256=prepared.case_manifest_sha256,
        selected_sample_sha256=prepared.selected_sample_sha256,
        hpo_input_sha256=prepared.hpo_input_sha256,
        phenotype_ingestion_sha256=hashlib.sha256(read_bounded_metadata(phenotype)).hexdigest(),
        ontology_manifest_sha256=prepared.ontology_manifest_sha256,
        execution_config=execution_config, base_release_bytes=base_release_bytes,
        expected_base_release_sha256=prepared.base_release_sha256,
        executable_stack_sha256=prepared.executable_stack_sha256,
        runtime_image_digest=prepared.runtime_image_digest,
        public_validation_authority_bytes=public_validation_authority_bytes,
        expected_public_validation_authority_sha256=prepared.public_validation_authority_sha256,
        journals=tuple(sorted(journals, key=lambda row: row.stage)))
    return PreparedIntakeSource(publication, Path(original_input_path), source_root, source,
        closure, filter_root, filtered, working_root, working,
        source_routing_context, execution_config)


@dataclass(frozen=True, slots=True, repr=False)
class ReadDeltaRecord:
    stage_id: str
    addition: BatchAllele | None = None
    assessment: ReadAssessment | None = None

    def __post_init__(self):
        if (not self.stage_id or (self.addition is None) == (self.assessment is None)
                or self.addition is not None and type(self.addition) is not BatchAllele
                or self.assessment is not None and type(self.assessment) is not ReadAssessment):
            raise WgsReadRescueError("read delta record is invalid")

    @property
    def key(self):
        return ("addition:" + self.addition.evidence.allele_id if self.addition is not None
                else "assessment:" + self.assessment.assessment_id)

    @property
    def variant_class(self):
        return (self.addition or self.assessment).variant_class


@dataclass(frozen=True, slots=True, repr=False)
class ReadDeltaStage:
    stage_id: str
    output_artifact: ArtifactRecord
    shards: tuple[ScientificShard, ...]
    addition_count: int
    assessment_count: int
    class_counts: tuple[tuple[WgsClass, int], ...]
    stream_sha256: str

    def __post_init__(self):
        _sha(self.stream_sha256)
        if (not self.stage_id
                or self.output_artifact.classification is not ArtifactClassification.CONTROLLED
                or self.output_artifact.status != "present"
                or self.output_artifact.deleted_on is not None
                or self.output_artifact.retention != "delete"
                or any(type(count) is not int or count < 0
                       for count in (self.addition_count, self.assessment_count))
                or any(row.kind != "read_evidence" for row in self.shards)
                or len({row.relative_path for row in self.shards}) != len(self.shards)
                or tuple(kind for kind, _ in self.class_counts) != tuple(sorted(WgsClass))
                or any(type(count) is not int or count < 0 for _, count in self.class_counts)
                or sum(count for _, count in self.class_counts) != self.addition_count + self.assessment_count
                or sum(row.row_count for row in self.shards) != self.addition_count + self.assessment_count):
            raise WgsReadRescueError("read delta stage is invalid")

    @property
    def sha256(self):
        return hashlib.sha256(_bytes(self)).hexdigest()


def write_read_delta_stage(records, *, directory, stage_id, output_artifact, execution_config):
    """Externally order producer records and reject duplicate typed identities."""
    root = Path(directory)
    if not root.is_absolute() or any(part.is_symlink() for part in (root, *root.parents)):
        raise WgsReadRescueError("read delta directory is invalid")
    root.mkdir(parents=True, mode=0o700, exist_ok=True)
    scratch = root / "delta-sort.sqlite"
    connection = _open_index(scratch, execution_config)
    writer = _ShardWriter(root, "read_evidence",
        max_rows=execution_config.source_shards.max_uncompressed_rows,
        max_bytes=execution_config.source_shards.max_uncompressed_bytes)
    counts = Counter()
    additions = assessments = 0
    try:
        connection.execute("CREATE TABLE records(key TEXT PRIMARY KEY,payload BLOB) WITHOUT ROWID")
        for record in records:
            if type(record) is not ReadDeltaRecord or record.stage_id != stage_id:
                raise WgsReadRescueError("read delta stage membership changed")
            content = _bytes(record)
            if len(content) > execution_config.source_shards.max_uncompressed_bytes:
                raise WgsReadRescueError("read delta record exceeds capacity", category="response_capacity")
            connection.execute("INSERT INTO records VALUES(?,?)", (record.key, content))
            counts[record.variant_class] += 1
            additions += record.addition is not None
            assessments += record.assessment is not None
        connection.commit()
        for payload, in connection.execute("SELECT payload FROM records ORDER BY key"):
            record = _typed(json.loads(payload), ReadDeltaRecord, decode=True)
            writer.write(record, key=record.key)
        stage = ReadDeltaStage(stage_id, output_artifact, writer.close(), additions, assessments,
            tuple((kind, counts[kind]) for kind in sorted(WgsClass)), writer.stream.hexdigest())
        _write_manifest(root, "read-delta-stage.json", _bytes(stage))
        return stage
    except sqlite3.Error:
        raise WgsReadRescueError("read delta typed identity is duplicated or storage failed") from None
    finally:
        writer._close()
        connection.close()
        scratch.unlink(missing_ok=True)


@dataclass(frozen=True, slots=True, repr=False)
class AuthorizedReadDelta:
    parent_prepared_sha256: str
    parent_early_snapshot_sha256: str
    case_manifest_sha256: str
    case_id: str
    sample_id: str
    source_bindings: tuple[IntakeSourceBinding, ...]
    execution_config_sha256: str
    use_binding: EvidenceUseBinding
    added_read_inventory: AddedReadInventory
    reference_authorizations: tuple[ReferenceAuthorizationBinding, ...]
    read_identity: ReadVcfIdentityReport
    producer_stack: ProducerStack
    selected_plan: SelectedPlan
    stage_states: tuple[StageState, ...]
    stage_results: tuple[ReadDeltaStage, ...]
    class_assessments: tuple[ClassAssessmentState, ...]
    analysis_mode: Literal["diagnosis_blind_primary"] = "diagnosis_blind_primary"
    schema: Literal["genosapiens.wgs-authorized-read-delta"] = "genosapiens.wgs-authorized-read-delta"
    version: Literal[1] = 1

    def __post_init__(self):
        for digest in (self.parent_prepared_sha256, self.parent_early_snapshot_sha256,
                       self.case_manifest_sha256, self.execution_config_sha256):
            _sha(digest)
        states = {row.stage_id: row for row in self.stage_states}
        results = {row.stage_id: row for row in self.stage_results}
        sources = tuple(row.input_sha256 for row in self.source_bindings)
        if (self.schema != "genosapiens.wgs-authorized-read-delta" or self.version != 1
                or self.analysis_mode != "diagnosis_blind_primary" or not self.case_id or not self.sample_id
                or not sources or sources != tuple(sorted(set(sources)))
                or self.parent_prepared_sha256 not in {row.prepared_sha256 for row in self.source_bindings}
                or tuple(row.stage_id for row in self.stage_states)
                != tuple(row.stage_id for row in self.selected_plan.stages)
                or tuple(results) != tuple(row.stage_id for row in self.stage_states
                                          if row.status is StageStatus.SUCCEEDED)
                or len(results) != len(self.stage_results)
                or any(states[row.stage_id].result_sha256 != row.sha256 for row in self.stage_results)
                or tuple(row.variant_class for row in self.class_assessments) != tuple(sorted(WgsClass))):
            raise WgsReadRescueError("authorized read delta parent or plan is invalid")
        for assessment in self.class_assessments:
            stage_ids = tuple(row.stage_id for row in self.selected_plan.stages
                              if assessment.variant_class in row.variant_classes)
            statuses = {states[key].status for key in stage_ids}
            if not stage_ids:
                status, reason = ClassAssessmentStatus.NOT_ASSESSED, "class-not-selected"
            elif StageStatus.FAILED in statuses:
                status, reason = ClassAssessmentStatus.FAILED, "stage-failed"
            elif StageStatus.INCONCLUSIVE in statuses:
                status, reason = ClassAssessmentStatus.INCONCLUSIVE, "stage-inconclusive"
            elif StageStatus.NOT_RUN in statuses:
                status, reason = ClassAssessmentStatus.NOT_ASSESSED, "stage-not-run"
            else:
                activity = any(dict(results[key].class_counts)[assessment.variant_class] for key in stage_ids)
                status, reason = (ClassAssessmentStatus.ASSESSED if activity else ClassAssessmentStatus.ABSENT), None
            if (assessment.stage_ids != stage_ids or assessment.status is not status
                    or assessment.reason_code != reason or assessment.scope != "selected-plan"):
                raise WgsReadRescueError("authorized read delta class accounting is invalid")

    def as_bytes(self):
        content = _bytes({"classification": "controlled", "llm_policy": "no-llm", "manifest": self})
        if len(content) > 8 * 1024 * 1024:
            raise WgsReadRescueError("authorized read delta exceeds capacity", category="response_capacity")
        return content

    @property
    def sha256(self):
        return hashlib.sha256(self.as_bytes()).hexdigest()


def load_authorized_read_delta(content, *, expected_sha256):
    return _load_manifest(content, expected_sha256=expected_sha256,
                          result_type=AuthorizedReadDelta, result_key="manifest")


def iter_read_delta_stage(stage, *, directory, execution_config):
    if any(row.size_bytes > execution_config.source_shards.max_uncompressed_bytes
            or row.row_count > execution_config.source_shards.max_uncompressed_rows for row in stage.shards):
        raise WgsReadRescueError("read delta shard exceeds capacity", category="response_capacity")
    def records():
        for _, _, row in _read_rows(Path(directory), stage.shards, "read_evidence", ReadDeltaRecord,
                expected_stream_sha256=stage.stream_sha256, row_key=lambda row: row.key):
            yield row
    counts = Counter()
    additions = assessments = 0
    previous = None
    for row in records():
        if row.stage_id != stage.stage_id or previous is not None and row.key <= previous:
            raise WgsReadRescueError("read delta record order or membership changed")
        previous = row.key
        counts[row.variant_class] += 1
        additions += row.addition is not None
        assessments += row.assessment is not None
    if ((additions, assessments) != (stage.addition_count, stage.assessment_count)
            or tuple((kind, counts[kind]) for kind in sorted(WgsClass)) != stage.class_counts):
        raise WgsReadRescueError("read delta record accounting changed")
    yield from records()


@dataclass(frozen=True, slots=True, repr=False)
class ReadDeltaReplay:
    delta: AuthorizedReadDelta
    stage_directories: tuple[tuple[str, Path], ...]
    artifacts: tuple[ArtifactReplay, ...]
    reference_authorizations: tuple[ReferenceAuthorizationReplay, ...]
    expectations: AuthorizedEvidenceExpectations

    def __post_init__(self):
        if (tuple(key for key, _ in self.stage_directories)
                != tuple(row.stage_id for row in self.delta.stage_results)
                or any(not path.is_absolute() for _, path in self.stage_directories)
                or tuple(row.artifact_id for row in self.artifacts)
                != tuple(sorted({row.artifact_id for row in self.artifacts}))
                or tuple(row.source_input_sha256 for row in self.reference_authorizations)
                != tuple(row.input_sha256 for row in self.delta.source_bindings)):
            raise WgsReadRescueError("read delta replay roster is invalid")

    def records(self, config):
        for stage, (_, directory) in zip(self.delta.stage_results, self.stage_directories, strict=True):
            yield from iter_read_delta_stage(stage, directory=directory, execution_config=config)


def _signature(identity, reference_sha256):
    return hashlib.sha256(_bytes({"identity": identity,
        "reference_manifest_sha256": reference_sha256})).hexdigest()


def _target(allele, call):
    fields = call.source.source_fields or {}
    return (fields.get("source_input_sha256"),
        None if call.sample is None else call.sample.sample_name,
        allele.allele_id, allele.native_identity_id, fields.get("source_record_id"),
        fields.get("source_allele_id"),
        tuple(fields.get("technical_flags", ())) or ("source-allele-read-validation-required",),
        "supplied-exact-event")


@dataclass(frozen=True, slots=True, repr=False)
class IntakeOccurrence:
    source_input_sha256: str
    source_manifest_sha256: str
    source_record_id: str
    source_allele_id: str
    payload_sha256: str
    preparation_disposition: str
    early_filter_disposition: str | None

    def __post_init__(self):
        for digest in (self.source_input_sha256, self.source_manifest_sha256, self.payload_sha256):
            _sha(digest)
        if not self.source_record_id or not self.source_allele_id or not self.preparation_disposition:
            raise WgsReadRescueError("intake source occurrence is invalid")


@dataclass(frozen=True, slots=True, repr=False)
class IntakeRecord:
    key: str
    occurrence: IntakeOccurrence | None = None
    admission: LedgerAdmissionInput | None = None
    read_evidence: ReadDeltaRecord | None = None
    delta_sha256: str | None = None

    def __post_init__(self):
        if sum(value is not None for value in (self.occurrence, self.admission, self.read_evidence)) != 1:
            raise WgsReadRescueError("intake record is invalid")
        if self.occurrence is not None:
            expected = "occurrence:" + self.occurrence.source_allele_id
        elif self.admission is not None:
            expected = "admission:" + self.admission.subject_id
        else:
            _sha(self.delta_sha256)
            expected = "read:" + self.delta_sha256 + ":" + self.read_evidence.key
        if self.key != expected or (self.read_evidence is None and self.delta_sha256 is not None):
            raise WgsReadRescueError("intake record identity is invalid")


@dataclass(frozen=True, slots=True, repr=False)
class MergedIntakeManifest:
    source_bindings: tuple[IntakeSourceBinding, ...]
    deltas: tuple[AuthorizedReadDelta, ...]
    case_manifest_sha256: str
    execution_config_sha256: str
    shards: tuple[ScientificShard, ...]
    record_counts: tuple[tuple[str, int], ...]
    admitted_source_occurrences: int
    remaining_technical_obligations: int
    mandatory_stage_blockers: tuple[str, ...]
    competitive_eligibility: Literal["not-established", "ineligible"]
    stream_sha256: str
    sample_id: str
    schema: Literal["genosapiens.wgs-merged-read-intake"] = "genosapiens.wgs-merged-read-intake"
    version: Literal[1] = 1

    def __post_init__(self):
        for digest in (self.case_manifest_sha256, self.execution_config_sha256, self.stream_sha256):
            _sha(digest)
        counts = dict(self.record_counts)
        blockers = tuple(sorted({"mandatory-stage:" + state.stage_id + ":" + state.status.value
            for delta in self.deltas for definition, state in zip(delta.selected_plan.stages, delta.stage_states, strict=True)
            if definition.mandatory and state.status is not StageStatus.SUCCEEDED}))
        if (self.schema != "genosapiens.wgs-merged-read-intake" or type(self.version) is not int or self.version != 1
                or type(self.sample_id) is not str or not self.sample_id
                or not self.source_bindings
                or tuple(row.input_sha256 for row in self.source_bindings)
                    != tuple(sorted({row.input_sha256 for row in self.source_bindings}))
                or tuple(row.sha256 for row in self.deltas) != tuple(sorted({row.sha256 for row in self.deltas}))
                or any(row.source_bindings != self.source_bindings
                       or row.sample_id != self.sample_id
                       or row.case_manifest_sha256 != self.case_manifest_sha256
                       or row.execution_config_sha256 != self.execution_config_sha256 for row in self.deltas)
                or tuple(counts) != ("admission", "assessment", "addition", "occurrence")
                or len(counts) != len(self.record_counts)
                or any(type(count) is not int or count < 0 for count in counts.values())
                or type(self.admitted_source_occurrences) is not int
                or self.admitted_source_occurrences < counts.get("admission", 0)
                or type(self.remaining_technical_obligations) is not int
                or self.remaining_technical_obligations < 0
                or self.mandatory_stage_blockers != blockers
                or counts.get("occurrence") != sum(source.original_occurrences for source in self.source_bindings)
                or any(row.kind != "intake" for row in self.shards)
                or len({row.relative_path for row in self.shards}) != len(self.shards)
                or sum(row.row_count for row in self.shards) != sum(counts.values())):
            raise WgsReadRescueError("merged intake manifest is invalid")
        eligible = (not self.remaining_technical_obligations and not self.mandatory_stage_blockers
            and all(row.use_binding.intended_use is EvidenceUse.CONTROLLED_PRIMARY for row in self.deltas))
        if self.competitive_eligibility != ("not-established" if eligible else "ineligible"):
            raise WgsReadRescueError("merged intake eligibility is invalid")

    def as_bytes(self):
        content = _bytes({"classification": "controlled", "llm_policy": "no-llm", "manifest": self})
        if len(content) > 8 * 1024 * 1024:
            raise WgsReadRescueError("merged intake metadata exceeds capacity", category="response_capacity")
        return content

    @property
    def sha256(self):
        return hashlib.sha256(self.as_bytes()).hexdigest()

    @property
    def annotation_binding(self):
        return AnnotationIntakeBinding(self.sha256, dict(self.record_counts)["admission"],
                                      self.admitted_source_occurrences)


def load_merged_intake(content, *, expected_sha256):
    return _load_manifest(content, expected_sha256=expected_sha256,
                          result_type=MergedIntakeManifest, result_key="manifest")


def validate_read_intake_run(manifest, *, run_evidence, release):
    """Check current rescue authority independently of the final CSV readiness."""
    run = run_evidence
    if (type(manifest) is not MergedIntakeManifest
            or manifest.competitive_eligibility != "not-established"
            or run.input_manifest_sha256 not in {row.prepared_sha256 for row in manifest.source_bindings}
            or manifest.execution_config_sha256 != release.policy_bindings.bounded_wgs_execution_sha256):
        raise WgsReadRescueError("competitive read intake authority is incomplete")
    expected = (("performed", "authorized", "complete", manifest.sha256) if manifest.deltas
                else ("not_performed", "not_applicable", "not_applicable", None))
    if (run.rescue_execution_status, run.rescue_authorization_status,
            run.rescue_plan_status, run.rescue_authorization_sha256) != expected:
        raise WgsReadRescueError("competitive intake rescue declaration changed")
    sources = {row.input_sha256: row for row in manifest.source_bindings}
    for delta in manifest.deltas:
        parent = next(row for row in manifest.source_bindings if row.prepared_sha256 == delta.parent_prepared_sha256)
        identity, inventory = delta.read_identity, delta.added_read_inventory
        if (delta.use_binding.intended_use is not EvidenceUse.CONTROLLED_PRIMARY
                or delta.use_binding.candidate_stack_sha256 != release.executable_stack_sha256
                or delta.use_binding.release_freeze_sha256 != release.sha256
                or delta.parent_early_snapshot_sha256 != parent.early_snapshot_sha256
                or identity.status is not IdentityStatus.PASSED
                or identity.added_read_inventory_sha256 != inventory.sha256
                or identity.source_input_sha256s != tuple(sources)
                or tuple(row.source_input_sha256 for row in delta.reference_authorizations) != tuple(sources)
                or identity.validated_input_sha256s != tuple(sorted(row.authorization.validated_input_sha256
                                                                  for row in delta.reference_authorizations))
                or any((value.parent_early_snapshot_sha256, value.case_id, value.sample_id)
                       != (parent.early_snapshot_sha256, delta.case_id, delta.sample_id)
                       for value in (identity, inventory))
                or any(row.authorization.reference_manifest_sha256 != run.accepted_reference_manifest_sha256
                       or row.authorization.reference_identity_id != run.accepted_reference_identity_id
                       for row in delta.reference_authorizations)):
            raise WgsReadRescueError("competitive read delta release or source authority changed")
    parent = next(row for row in manifest.source_bindings if row.prepared_sha256 == run.input_manifest_sha256)
    if parent.early_snapshot_sha256 != run.early_snapshot_sha256:
        raise WgsReadRescueError("competitive read intake early parent changed")


class _SourceIndex:
    """Index the full source universe, then its admitted subset, on bounded disk."""

    def __init__(self, sources, path, config):
        self.path, self.config = path, config
        self.connection = _open_index(path, config)
        self.connection.executescript("""
            CREATE TABLE occurrences(id TEXT PRIMARY KEY, input_sha TEXT, record_id TEXT,
                classes BLOB, payload BLOB NOT NULL);
            CREATE TABLE identities(id TEXT PRIMARY KEY);
            CREATE TABLE geometries(signature TEXT PRIMARY KEY);
            CREATE TABLE admissions(id TEXT PRIMARY KEY, payload BLOB NOT NULL);
            CREATE TABLE obligations(target BLOB PRIMARY KEY, classes BLOB NOT NULL);
            CREATE TABLE assessments(target BLOB PRIMARY KEY, id TEXT UNIQUE, payload BLOB NOT NULL);
            CREATE TABLE additions(id TEXT PRIMARY KEY, native_id TEXT UNIQUE,
                signature TEXT UNIQUE, payload BLOB NOT NULL);
            CREATE TABLE output(key TEXT PRIMARY KEY, payload BLOB NOT NULL);
        """)
        try:
            for source in sources:
                self._source(source)
        except BaseException:
            self.close()
            raise

    def close(self):
        self.connection.close()
        self.path.unlink(missing_ok=True)

    def _bounded(self, content):
        if len(content) > min(self.config.annotation.max_response_row_bytes,
                              self.config.source_shards.max_uncompressed_bytes):
            raise WgsReadRescueError("merged intake row exceeds capacity", category="response_capacity")
        return content

    def _source(self, source):
        db = self.connection
        for row in iter_ledger(source.source_root, source.source):
            for source_id in row.source_allele_ids:
                occurrence = _occurrence(row.payload_json, source_allele_id=source_id,
                    source_manifest=source.source, contract=source.native_contract)
                normalized = occurrence.normalization.normalized_allele
                if normalized is not None:
                    db.execute("INSERT OR IGNORE INTO geometries VALUES(?)", (
                        _signature(_normalized_sequence(occurrence), source.native_contract.reference_manifest_sha256),))
                    for identity_id in (normalized.normalized_allele_id, "native-" + normalized.normalized_allele_id):
                        db.execute("INSERT OR IGNORE INTO identities VALUES(?)", (identity_id,))
                db.execute("INSERT INTO identities VALUES(?)", (source_id,))
                classes = set(occurrence.inventory.represented_classes)
                # Mosaic is a selected-call screen outcome, separate from the
                # record's static representation classes. The decoder already
                # authenticates the original screen and selected sample.
                if occurrence.screen is not None and occurrence.screen.outcome is ScreenOutcome.CANDIDATE:
                    classes.add(WgsClass.MOSAIC)
                db.execute("INSERT INTO occurrences VALUES(?,?,?,?,?)", (source_id,
                    source.source.source.input_sha256, row.source_record_id,
                    _bytes(tuple(sorted(classes))), b""))
            db.commit()
        for shard in source.filtered.shards:
            for decision in iter_filter_shard(source.filter_root, shard):
                row = db.execute("SELECT input_sha,record_id FROM occurrences WHERE id=?", (decision.source_allele_id,)).fetchone()
                if row is None:
                    raise WgsReadRescueError("source disposition has no original occurrence")
                # The original payload hash is read below; disposition is never an absence claim.
                db.execute("UPDATE occurrences SET payload=? WHERE id=?", (
                    _bytes((decision.preparation_disposition, decision.early_filter_disposition)), decision.source_allele_id))
        for row in iter_ledger(source.source_root, source.source):
            for source_id in row.source_allele_ids:
                disposition = db.execute("SELECT payload FROM occurrences WHERE id=?", (source_id,)).fetchone()[0]
                if not disposition:
                    raise WgsReadRescueError("original source disposition is missing")
                preparation, early = json.loads(disposition)
                occurrence = IntakeOccurrence(source.source.source.input_sha256, source.source.sha256,
                    row.source_record_id, source_id, hashlib.sha256(row.payload_json.encode("ascii")).hexdigest(),
                    preparation, early)
                self.output(IntakeRecord("occurrence:" + source_id, occurrence=occurrence))
        for value in iter_ledger_admission_inputs(source_root=source.source_root,
                source_manifest=source.source, working_root=source.working_root,
                working_manifest=source.working,
                source_routing_context=source.source_routing_context,
                scratch_path=self.path.with_name("rescue-admission.sqlite"), execution_config=self.config):
            previous = db.execute("SELECT payload FROM admissions WHERE id=?", (value.subject_id,)).fetchone()
            combined = value if previous is None else _combine_admissions(load_ledger_admission_input(previous[0]), value)
            db.execute("INSERT OR REPLACE INTO admissions VALUES(?,?)", (value.subject_id, self._bounded(combined.as_bytes())))
            allele = value.allele
            if allele is None:
                continue
            for identity_id in (allele.allele_id, *(row.identity_id for row in allele.identities)):
                db.execute("INSERT OR IGNORE INTO identities VALUES(?)", (identity_id,))
            native = next(row for row in allele.identities if row.identity_id == allele.native_identity_id)
            db.execute("INSERT OR IGNORE INTO geometries VALUES(?)", (
                _signature(native.identity, native.reference_manifest_sha256),))
            for call in allele.calls or ():
                target = _target(allele, call)
                source_row = db.execute("SELECT input_sha,record_id,classes FROM occurrences WHERE id=?", (target[5],)).fetchone()
                if source_row is None or source_row[:2] != (target[0], target[4]):
                    raise WgsReadRescueError("source call occurrence binding is invalid")
                if (call.source.source_fields or {}).get("read_validation_required") is True:
                    db.execute("INSERT INTO obligations VALUES(?,?)", (_bytes(target), source_row[2]))
        db.commit()

    def output(self, row):
        self.connection.execute("INSERT INTO output VALUES(?,?)", (row.key, self._bounded(_bytes(row))))


def _validate_delta_authority(replay, sources, case, *, reference_root):
    if any(source.prepared.manifest.public_validation_authority_sha256 is not None for source in sources):
        raise WgsReadRescueError("public validation authority cannot authorize read rescue")
    delta, expected = replay.delta, replay.expectations
    config = sources[0].execution_config
    bindings = tuple(source.binding for source in sources)
    parents = {source.binding.prepared_sha256: source for source in sources}
    parent = parents.get(delta.parent_prepared_sha256)
    if (parent is None or delta.source_bindings != bindings
            or delta.parent_early_snapshot_sha256 != parent.binding.early_snapshot_sha256
            or delta.case_manifest_sha256 != hashlib.sha256(case_manifest_as_bytes(case)).hexdigest()
            or (delta.case_id, delta.sample_id) != (case.case_id, case.sample_id)
            or delta.execution_config_sha256 != execution_sha256(config)
            or delta.producer_stack != expected.producer_stack
            or delta.selected_plan != expected.selected_plan
            or delta.use_binding != expected.use_binding
            or delta.read_identity.identity_policy_sha256 != expected.identity_policy_sha256):
        raise WgsReadRescueError("authorized read delta parent or policy changed")
    inventory, identity = delta.added_read_inventory, delta.read_identity
    if any((item.parent_early_snapshot_sha256, item.case_id, item.sample_id)
            != (delta.parent_early_snapshot_sha256, case.case_id, case.sample_id)
            for item in (inventory, identity)):
        raise WgsReadRescueError("read inventory parent binding changed")
    expected_sources = tuple(row.input_sha256 for row in bindings)
    if tuple(row.source_input_sha256 for row in delta.reference_authorizations) != expected_sources:
        raise WgsReadRescueError("reference authorization omits an original source")
    for source, claim, item in zip(sources, delta.reference_authorizations, replay.reference_authorizations, strict=True):
        if item.input_path != source.original_input_path:
            raise WgsReadRescueError("reference authorization original source changed")
        authorization = require_rescue_reference_concordance(identity_path=item.identity_path,
            report_path=item.report_path, input_path=item.input_path,
            validated_input_path=item.validated_input_path, reference_manifest_path=item.reference_manifest_path,
            reference_root=reference_root, order_adjudication_path=item.order_adjudication_path,
            source_report_path=item.source_report_path, source_reference_manifest_path=item.source_reference_manifest_path,
            source_reference_root=item.source_reference_root)
        if (authorization != claim.authorization
                or authorization.reference_manifest_sha256 != case.reference_bundle.manifest_sha256
                or authorization.reference_manifest_sha256 != source.native_contract.reference_manifest_sha256
                or authorization.reference_identity_id != source.native_contract.reference_identity_id):
            raise WgsReadRescueError("reference authorization replay changed")
    if (identity.status is not IdentityStatus.PASSED
            or identity.source_input_sha256s != expected_sources
            or identity.validated_input_sha256s != tuple(sorted(row.authorization.validated_input_sha256
                                                              for row in delta.reference_authorizations))
            or identity.added_read_inventory_sha256 != inventory.sha256):
        raise WgsReadRescueError("read and original call identity did not pass")
    artifacts = tuple(row.artifact for row in inventory.reads) + identity.supporting_artifacts + tuple(
        row.output_artifact for row in delta.stage_results)
    _verify_artifacts(artifacts, replay.artifacts, reference_root=reference_root)


def _validate_addition(record, delta, stage, index):
    evidence = record.addition.evidence
    db = index.connection
    native = next(row for row in evidence.identities if row.identity_id == evidence.native_identity_id)
    signature = _signature(native.identity, native.reference_manifest_sha256)
    references = {(ReferenceRole.CHALLENGE_ALIGNMENT, row.authorization.reference_identity_id,
                   row.authorization.reference_manifest_sha256) for row in delta.reference_authorizations}
    provenance = tuple(row.source for row in evidence.identities) + tuple(row.source for row in evidence.calls or ())
    if (not evidence.calls or evidence.population is not None or evidence.clinical is not None
            or evidence.transcript_consequences is not None or evidence.mechanism_ids or evidence.evidence_refs
            or any((*_identity_reference(row.identity), row.reference_manifest_sha256) not in references
                   for row in evidence.identities)
            or any(row.source_sha256 != stage.output_artifact.sha256 or row.public_manifest_sha256 is not None
                   or row.source_dataset != delta.producer_stack.producer_id
                   or row.source_release != delta.producer_stack.release for row in provenance)
            or any(call.sample is None or call.sample.sample_name != delta.sample_id for call in evidence.calls)
            or any("authorized_read_evidence" in (row.source_fields or {})
                   or "authorized_read_assessment" in (row.source_fields or {}) for row in provenance)
            or any(db.execute("SELECT 1 FROM identities WHERE id=?", (identity_id,)).fetchone()
                   for identity_id in (evidence.allele_id, *(row.identity_id for row in evidence.identities)))
            or db.execute("SELECT 1 FROM geometries WHERE signature=?", (signature,)).fetchone()):
        raise WgsReadRescueError("read addition collides or has invalid provenance")
    db.execute("INSERT INTO additions VALUES(?,?,?,?)", (evidence.allele_id, evidence.native_identity_id,
        signature, index._bounded(_bytes(record))))
    for identity_id in {evidence.allele_id, *(row.identity_id for row in evidence.identities)}:
        db.execute("INSERT INTO identities VALUES(?)", (identity_id,))
    db.execute("INSERT INTO geometries VALUES(?)", (signature,))
    admission = _addition_input(evidence, delta=delta, stage=stage)
    for source_id in admission.source_allele_ids:
        if db.execute("SELECT 1 FROM occurrences WHERE id=?", (source_id,)).fetchone():
            raise WgsReadRescueError("read addition source identity collides")
    db.execute("INSERT INTO admissions VALUES(?,?)", (admission.subject_id, index._bounded(admission.as_bytes())))


def _combine_admissions(first, second):
    """Preserve calls from multiple original VCFs sharing one normalized identity."""
    if (first.allele is None or second.allele is None or first.transport != second.transport
            or first.identity_id != second.identity_id or first.allele.native_identity_id != second.allele.native_identity_id
            or set(first.source_allele_ids).intersection(second.source_allele_ids)):
        raise WgsReadRescueError("original source admission identities conflict")
    identities = {row.identity_id: row for row in first.allele.identities}
    for row in second.allele.identities:
        previous = identities.get(row.identity_id)
        if previous is not None:
            if (previous.identity != row.identity or previous.reference_manifest_sha256 != row.reference_manifest_sha256
                    or previous.stage is not row.stage or previous.transformation_id != row.transformation_id):
                raise WgsReadRescueError("original source identity geometry conflicts")
            parents = tuple(sorted(set((*previous.derived_from_identity_ids, *row.derived_from_identity_ids))))
            source = SourcedRecord("genosapiens.wgs-merged-normalization", "1", row.identity_id,
                hashlib.sha256(_bytes((previous.source, row.source))).hexdigest(), None,
                {"contributing_source_records": json.loads(_bytes((previous.source, row.source)))})
            row = replace(row, derived_from_identity_ids=parents, source=source)
        identities[row.identity_id] = row
    calls = tuple(sorted((*first.allele.calls, *second.allele.calls),
                         key=lambda call: call.source.source_fields["source_allele_id"]))
    obligations = tuple(sorted((*first.technical_obligations, *second.technical_obligations),
                               key=lambda row: row.source_allele_id))
    return replace(first, source_allele_ids=tuple(row.source_allele_id for row in obligations),
        technical_obligations=obligations,
        allele=replace(first.allele, identities=tuple(identities[key] for key in sorted(identities)), calls=calls))


def _addition_input(evidence, *, delta, stage):
    native = next(row for row in evidence.identities if row.identity_id == evidence.native_identity_id)
    identities = {row.identity_id: row for row in evidence.identities}
    calls, obligations = [], []
    authority = json.loads(_bytes({
        "parent_prepared_sha256": delta.parent_prepared_sha256,
        "parent_early_snapshot_sha256": delta.parent_early_snapshot_sha256,
        "added_read_inventory_sha256": delta.added_read_inventory.sha256,
        "read_identity_report_sha256": delta.read_identity.sha256,
        "producer_stack_sha256": delta.producer_stack.sha256,
        "selected_plan_sha256": delta.selected_plan.sha256,
        "stage_output_sha256": stage.output_artifact.sha256,
        "read_delta_sha256": delta.sha256,
        "reference_authorizations": delta.reference_authorizations,
    }))
    for call in evidence.calls:
        fields = dict(call.source.source_fields or {})
        source_id = fields.get("source_allele_id")
        mapping = identities.get(source_id)
        if mapping is None:
            source_id = "source-allele-" + hashlib.sha256(_bytes((stage.output_artifact.sha256,
                call.source.source_record_id, evidence.allele_id, call.sample.sample_name))).hexdigest()[:24]
            mapping = IdentityMapping(source_id, IdentityStage.SOURCE, native.identity,
                                      native.reference_manifest_sha256, call.source)
            identities[source_id] = mapping
        if mapping.stage is not IdentityStage.SOURCE or mapping.source != call.source:
            raise WgsReadRescueError("read addition source call lineage is invalid")
        flags = tuple(fields.get("technical_flags", ()))
        required = fields.get("read_validation_required", bool(flags))
        if (type(required) is not bool or any(type(flag) is not str or not flag for flag in flags)
                or type(fields.get("allele_number")) is not int or fields["allele_number"] < 1):
            raise WgsReadRescueError("read addition technical facts are invalid")
        record_id = "read-record-" + hashlib.sha256(_bytes((stage.output_artifact.sha256,
                                                          call.source.source_record_id))).hexdigest()[:24]
        fields.update(source_allele_id=source_id, source_record_id=record_id,
            producer_source_record_id=call.source.source_record_id,
            source_input_sha256=stage.output_artifact.sha256, authorized_read_evidence=authority,
            read_validation_required=required, technical_flags=flags, storage_envelope="bounded-read-delta-v1")
        source = replace(call.source, source_fields=fields)
        identities[source_id] = replace(mapping, source=source)
        calls.append(replace(call, source=source))
        obligations.append(SourceTechnicalObligation(source_id, call.sample.sample_name, flags, required,
                                                     bool(call.filters and call.filters != ("PASS",))))
    if len({row.source_allele_id for row in obligations}) != len(obligations):
        raise WgsReadRescueError("read addition contains duplicate calls")
    if not native.derived_from_identity_ids:
        identities[native.identity_id] = replace(native,
            derived_from_identity_ids=tuple(sorted(row.source_allele_id for row in obligations)))
    # Native-only evidence remains native-only: normalization is never invented.
    sequence = native.identity.sequence if isinstance(native.identity, (MitochondrialAllele, MosaicAllele)) else native.identity
    sequence_route = type(sequence) is SequenceAllele
    transport = "sequence_vep" if sequence_route else "native_vep" if type(native.identity) is StructuralAllele else "repeat_catalog"
    rules = (RULE_IDS[0], RULE_IDS[2], RULE_IDS[4]) if sequence_route else (
        (RULE_IDS[1], RULE_IDS[2], RULE_IDS[4]) if transport == "native_vep" else (RULE_IDS[2], RULE_IDS[3], RULE_IDS[4]))
    ordered = sorted(zip(calls, obligations), key=lambda pair: pair[1].source_allele_id)
    return LedgerAdmissionInput(1, evidence.allele_id, tuple(row.source_allele_id for _, row in ordered),
        transport, replace(evidence, identities=tuple(identities[key] for key in sorted(identities)),
                           calls=tuple(call for call, _ in ordered)), None, rules,
        tuple(row for _, row in ordered), delta.sha256, delta.sha256,
        hashlib.sha256(_bytes(delta.source_bindings)).hexdigest(), "pending", ())


def validate_screening_read_lineage(call, obligation, manifest):
    """Bind a read-source call to the intake's already replayed delta and stage.

    MergedIntakeIndex remains the complete payload/measurement validation owner.
    This check preserves those exact parents across the admission worker boundary.
    """
    fields = call.source.source_fields or {}
    authority = fields.get('authorized_read_evidence')
    if not isinstance(authority, Mapping):
        raise WgsReadRescueError('screening read lineage is absent')
    matches = [(delta, stage) for delta in manifest.deltas for stage in delta.stage_results
        if stage.output_artifact.sha256 == fields.get('source_input_sha256')
        and delta.sha256 == authority.get('read_delta_sha256')]
    if len(matches) != 1:
        raise WgsReadRescueError('screening read lineage is undeclared')
    delta, stage = matches[0]
    expected = dict(parent_prepared_sha256=delta.parent_prepared_sha256,
        parent_early_snapshot_sha256=delta.parent_early_snapshot_sha256,
        added_read_inventory_sha256=delta.added_read_inventory.sha256,
        read_identity_report_sha256=delta.read_identity.sha256,
        producer_stack_sha256=delta.producer_stack.sha256,
        selected_plan_sha256=delta.selected_plan.sha256,
        stage_output_sha256=stage.output_artifact.sha256,
        read_delta_sha256=delta.sha256, reference_authorizations=delta.reference_authorizations,
        merged_intake_sha256=manifest.sha256)
    if (_bytes(authority) != _bytes(expected)
            or call.sample is None or call.sample.sample_name != delta.sample_id
            or obligation.sample_name != delta.sample_id
            or fields.get('source_allele_id') != obligation.source_allele_id
            or fields.get('producer_source_record_id') != call.source.source_record_id
            or fields.get('technical_flags') != obligation.technical_flags
            or fields.get('read_validation_required') is not obligation.read_validation_required
            or type(fields.get('allele_number')) is not int or fields['allele_number'] < 1):
        raise WgsReadRescueError('screening read lineage differs')


def iter_deferred_intake_targets(eligibility_root, *, manifest, source):
    """Address original unresolved targets; this view supplies no read evidence."""
    from genosapiens.wgs_analysis_eligibility import AnalysisEligibilityManifest,iter_refcall_rescue_candidates,analysis_eligibility_component
    if (type(source) is not PreparedIntakeSource or type(manifest) is not AnalysisEligibilityManifest
            or manifest.source_manifest_sha256!=source.binding.source_manifest_sha256
            or manifest.reference_manifest_sha256!=source.native_contract.reference_manifest_sha256):
        raise WgsReadRescueError('foreign deferred rescue source')
    analysis_eligibility_component(eligibility_root,manifest)
    yield from iter_refcall_rescue_candidates(eligibility_root,manifest=manifest)


def merge_deferred_eligibility_intake(replays, *, eligibility_root, eligibility_manifest,
        source, sources, case_manifest_path, reference_root, output_directory):
    """Bind a distinct rescue run to the old domain and existing delta validator.

    The targets are audit scope, never a gene/variant whitelist or a carriage
    certificate. Missing or failed assessments leave original GT and technical
    obligations unresolved through the ordinary merge and admission owners.
    """
    sources=tuple(sources)
    if source not in sources:
        raise WgsReadRescueError('deferred source missing from rescue intake')
    for target in iter_deferred_intake_targets(eligibility_root,manifest=eligibility_manifest,source=source):
        pass
    return merge_prepared_intake(replays,sources=sources,case_manifest_path=case_manifest_path,
        reference_root=reference_root,output_directory=output_directory)


def evaluate_deferred_intake_target(target,value,*,source,index):
    """Reassess a verified distinct intake while preserving original GT/history."""
    from genosapiens.wgs_analysis_eligibility import AnalysisEligibilityDecision,evaluate_refcall_eligibility
    from genosapiens.wgs_admission import _source_header
    if (type(target) is not AnalysisEligibilityDecision or target.disposition!='defer_refcall'
            or type(source) is not PreparedIntakeSource or type(index) is not MergedIntakeIndex
            or index.connection is None or source.binding not in index.intake.manifest.source_bindings
            or target.source_manifest_sha256!=source.source.sha256
            or value.parent_source_manifest_sha256!=index.intake.manifest.sha256
            or value.identity_id!=target.identity_id or value.source_allele_ids!=target.source_allele_ids
            or value.allele is None):
        raise WgsReadRescueError('foreign or unvalidated deferred reassessment')
    index.validate_allele(value.allele)
    # The validated intake supplies the binding between the distinct run and
    # this original source. The temporary predicate view changes no bytes in
    # either run and confers no alternate genotype/carriage.
    original_view=replace(value,parent_source_manifest_sha256=source.source.sha256)
    return evaluate_refcall_eligibility(original_view,
        header=_source_header(source.source_root,source.source),selected_sample=source.binding.selected_sample_name)


def merge_prepared_intake(replays, *, sources, case_manifest_path, reference_root,
                                 output_directory):
    """Authenticate every source/delta and persist an intake without scientific claims."""
    root = Path(output_directory)
    if (not sources or not root.is_absolute()
            or any(path.is_symlink() for path in (root, *root.parents))):
        raise WgsReadRescueError("merged intake inputs are invalid")
    sources = tuple(sources)
    source_ids = tuple(row.binding.input_sha256 for row in sources)
    case = load_case_manifest(case_manifest_path)
    case_sha = hashlib.sha256(case_manifest_as_bytes(case)).hexdigest()
    config = sources[0].execution_config
    if (source_ids != tuple(sorted(set(source_ids)))
            or source_ids != tuple(sorted(row.sha256 for row in case.inputs if row.file_class is CaseFileClass.VCF))
            or any(row.prepared.manifest.case_manifest_sha256 != case_sha
                   or row.execution_config != config
                   or row.prepared.manifest.base_release_sha256 != sources[0].prepared.manifest.base_release_sha256
                   or row.prepared.manifest.public_validation_authority_sha256 != sources[0].prepared.manifest.public_validation_authority_sha256
                   or row.prepared.manifest.executable_stack_sha256 != sources[0].prepared.manifest.executable_stack_sha256
                   or row.prepared.manifest.hpo_input_sha256 != sources[0].prepared.manifest.hpo_input_sha256
                   for row in sources)):
        raise WgsReadRescueError("merged intake omits or changes an original source")
    unique = {}
    for replay in replays:
        # Replay even duplicate submissions before treating identical intake as idempotent.
        _validate_delta_authority(replay, sources, case, reference_root=reference_root)
        if replay.delta.sha256 in unique:
            for _ in replay.records(config):
                pass
        unique.setdefault(replay.delta.sha256, replay)
    ordered = tuple(unique[key] for key in sorted(unique))
    if len({row.delta.use_binding for row in ordered}) > 1:
        raise WgsReadRescueError("read deltas have incompatible release authority")
    root.mkdir(mode=0o700, parents=True, exist_ok=True)
    index = _SourceIndex(sources, root / "merge.sqlite", config)
    db = index.connection
    writer = _ShardWriter(root, "intake", max_rows=config.source_shards.max_uncompressed_rows,
                         max_bytes=config.source_shards.max_uncompressed_bytes)
    counts = Counter()
    blockers = set()
    remaining = admitted_occurrences = 0
    try:
        for replay in ordered:
            delta = replay.delta
            stages = {row.stage_id: row for row in delta.stage_results}
            definitions = {row.stage_id: row for row in delta.selected_plan.stages}
            blockers.update("mandatory-stage:" + row.stage_id + ":" + row.status.value
                for row in delta.stage_states if definitions[row.stage_id].mandatory
                and row.status is not StageStatus.SUCCEEDED)
            for record in replay.records(config):
                stage = stages[record.stage_id]
                if record.variant_class not in definitions[record.stage_id].variant_classes:
                    raise WgsReadRescueError("read delta stage class changed")
                if record.addition is not None:
                    _validate_addition(record, delta, stage, index)
                else:
                    item = record.assessment
                    target, payload = _bytes(_assessment_target(item)), _bytes(item)
                    obligation = db.execute("SELECT classes FROM obligations WHERE target=?", (target,)).fetchone()
                    if (obligation is None or item.variant_class.value not in json.loads(obligation[0])
                            or item.measurement_artifact_sha256 != stage.output_artifact.sha256):
                        raise WgsReadRescueError("read assessment does not match an exact call obligation")
                    prior = db.execute("SELECT payload FROM assessments WHERE target=? OR id=?", (target, item.assessment_id)).fetchall()
                    if prior:
                        if any(value[0] != payload for value in prior):
                            raise WgsReadRescueError("read assessments conflict")
                    else:
                        db.execute("INSERT INTO assessments VALUES(?,?,?)", (target, item.assessment_id, payload))
                index.output(IntakeRecord("read:" + delta.sha256 + ":" + record.key, read_evidence=record, delta_sha256=delta.sha256))
            db.commit()
        for payload, in db.execute("SELECT payload FROM admissions ORDER BY id"):
            value = load_ledger_admission_input(payload)
            admitted_occurrences += len(value.source_allele_ids)
            if value.allele is not None:
                calls, obligations = [], []
                for call, obligation in zip(value.allele.calls or (), value.technical_obligations, strict=True):
                    match = db.execute("SELECT payload FROM assessments WHERE target=?", (_bytes(_target(value.allele, call)),)).fetchone()
                    fields = dict(call.source.source_fields or {})
                    discharged = False
                    if match is not None:
                        assessment = _typed(json.loads(match[0]), ReadAssessment, decode=True)
                        fields["authorized_read_assessment"] = json.loads(match[0])
                        discharged = assessment.outcome is AssessmentOutcome.VALIDATED
                        if discharged:
                            fields["read_validation_required"] = False
                    remaining += obligation.read_validation_required and not discharged
                    calls.append(replace(call, source=replace(call.source, source_fields=fields)))
                    obligations.append(replace(obligation,
                        read_validation_required=obligation.read_validation_required and not discharged))
                value = replace(value, allele=replace(value.allele, calls=tuple(calls)), technical_obligations=tuple(obligations))
            else:
                remaining += sum(row.read_validation_required for row in value.technical_obligations)
            index.output(IntakeRecord("admission:" + value.subject_id, admission=value))
        db.commit()
        for payload, in db.execute("SELECT payload FROM output ORDER BY key"):
            row = _typed(json.loads(payload), IntakeRecord, decode=True)
            kind = ("occurrence" if row.occurrence else "admission" if row.admission
                    else "assessment" if row.read_evidence.assessment else "addition")
            counts[kind] += 1
            writer.write(row, key=row.key)
        eligible = (not remaining and not blockers
            and all(row.delta.use_binding.intended_use is EvidenceUse.CONTROLLED_PRIMARY for row in ordered))
        result = MergedIntakeManifest(tuple(source.binding for source in sources),
            tuple(row.delta for row in ordered), case_sha, execution_sha256(config), writer.close(),
            tuple((kind, counts[kind]) for kind in ("admission", "assessment", "addition", "occurrence")),
            admitted_occurrences,
            remaining, tuple(sorted(blockers)), "not-established" if eligible else "ineligible",
            writer.stream.hexdigest(), case.sample_id)
        _write_manifest(root, "merged-intake-manifest.json", result.as_bytes())
        return result
    except sqlite3.Error:
        raise WgsReadRescueError("read intake collision or disk join failed") from None
    finally:
        writer._close()
        index.close()


def iter_merged_intake(manifest, *, directory, execution_config):
    if (manifest.execution_config_sha256 != execution_sha256(execution_config)
            or any(row.size_bytes > execution_config.source_shards.max_uncompressed_bytes
                   or row.row_count > execution_config.source_shards.max_uncompressed_rows for row in manifest.shards)):
        raise WgsReadRescueError("merged intake execution binding changed")
    def rows():
        for _, _, row in _read_rows(Path(directory), manifest.shards, "intake", IntakeRecord,
                expected_stream_sha256=manifest.stream_sha256, row_key=lambda row: row.key):
            yield row
    counts, previous = Counter(), None
    deltas = {row.sha256 for row in manifest.deltas}
    for row in rows():
        if previous is not None and row.key <= previous or row.read_evidence and row.delta_sha256 not in deltas:
            raise WgsReadRescueError("merged intake record closure changed")
        previous = row.key
        counts["occurrence" if row.occurrence else "admission" if row.admission
               else "assessment" if row.read_evidence.assessment else "addition"] += 1
    if tuple((kind, counts[kind]) for kind, _ in manifest.record_counts) != manifest.record_counts:
        raise WgsReadRescueError("merged intake record accounting changed")
    yield from rows()


class MergedIntakeIndex:
    """Authenticate persisted stage and merge projections before downstream use."""

    def __init__(self, intake, *, execution_config, scratch_path):
        self.intake, self.config, self.path = intake, execution_config, Path(scratch_path)
        self.connection = None

    def __enter__(self):
        db = self.connection = _open_index(self.path, self.config)
        try:
            db.executescript("""
                CREATE TABLE admissions(id TEXT PRIMARY KEY,payload BLOB);
                CREATE TABLE occurrences(id TEXT PRIMARY KEY,input_sha TEXT,source_sha TEXT,payload BLOB,used INTEGER DEFAULT 0);
                CREATE TABLE records(delta_sha TEXT,stage TEXT,key TEXT,payload BLOB,PRIMARY KEY(delta_sha,stage,key));
                CREATE TABLE assessments(target BLOB PRIMARY KEY,id TEXT UNIQUE,payload BLOB,used INTEGER DEFAULT 0);
                CREATE TABLE additions(id TEXT PRIMARY KEY,payload BLOB,used INTEGER DEFAULT 0);
            """)
            counts = Counter()
            sources = {row.input_sha256: row for row in self.intake.manifest.source_bindings}
            deltas = {row.sha256: row for row in self.intake.manifest.deltas}
            for row in iter_merged_intake(self.intake.manifest, directory=self.intake.directory, execution_config=self.config):
                if row.admission is not None:
                    db.execute("INSERT INTO admissions VALUES(?,?)", (row.admission.subject_id, row.admission.as_bytes()))
                elif row.occurrence is not None:
                    value = row.occurrence
                    if value.source_input_sha256 not in sources or value.source_manifest_sha256 != sources[value.source_input_sha256].source_manifest_sha256:
                        raise WgsReadRescueError("merged occurrence source changed")
                    db.execute("INSERT INTO occurrences(id,input_sha,source_sha,payload) VALUES(?,?,?,?)",
                        (value.source_allele_id, value.source_input_sha256, value.source_manifest_sha256, _bytes(value)))
                    counts[value.source_input_sha256] += 1
                else:
                    record, delta = row.read_evidence, deltas[row.delta_sha256]
                    stage = next((stage for stage in delta.stage_results if stage.stage_id == record.stage_id), None)
                    if stage is None:
                        raise WgsReadRescueError("merged delta stage is undeclared")
                    db.execute("INSERT INTO records VALUES(?,?,?,?)", (row.delta_sha256, record.stage_id, record.key, _bytes(record)))
                    if record.addition is not None:
                        expected = _addition_input(record.addition.evidence, delta=delta, stage=stage)
                        db.execute("INSERT INTO additions(id,payload) VALUES(?,?)", (expected.subject_id, expected.as_bytes()))
                    else:
                        assessment = record.assessment
                        if assessment.measurement_artifact_sha256 != stage.output_artifact.sha256:
                            raise WgsReadRescueError("merged assessment measurement binding changed")
                        target, payload = _bytes(_assessment_target(assessment)), _bytes(assessment)
                        prior = db.execute("SELECT payload FROM assessments WHERE target=? OR id=?", (target, assessment.assessment_id)).fetchall()
                        if any(value[0] != payload for value in prior):
                            raise WgsReadRescueError("merged assessments conflict")
                        if not prior:
                            db.execute("INSERT INTO assessments(target,id,payload) VALUES(?,?,?)", (target, assessment.assessment_id, payload))
            if any(counts[source.input_sha256] != source.original_occurrences for source in sources.values()):
                raise WgsReadRescueError("merged intake omitted original occurrences")
            for delta in deltas.values():
                for stage in delta.stage_results:
                    self._stage(delta, stage)
            remaining = admitted = 0
            for id_, payload in db.execute("SELECT id,payload FROM admissions ORDER BY id"):
                value = load_ledger_admission_input(payload)
                admitted += len(value.source_allele_ids)
                addition = db.execute("SELECT payload FROM additions WHERE id=?", (id_,)).fetchone()
                if addition is not None:
                    expected = load_ledger_admission_input(addition[0])
                    if value != expected:
                        raise WgsReadRescueError("merged read addition projection changed")
                    db.execute("UPDATE additions SET used=1 WHERE id=?", (id_,))
                    remaining += sum(row.read_validation_required for row in value.technical_obligations)
                    continue
                for source_id in value.source_allele_ids:
                    source = db.execute("SELECT input_sha,payload FROM occurrences WHERE id=?", (source_id,)).fetchone()
                    if source is None:
                        raise WgsReadRescueError("merged original admission has no source occurrence")
                    original = _typed(json.loads(source[1]), IntakeOccurrence, decode=True)
                    if original.early_filter_disposition == "exclude_common":
                        raise WgsReadRescueError("excluded original was restored without authority")
                    db.execute("UPDATE occurrences SET used=1 WHERE id=?", (source_id,))
                if value.allele is None:
                    remaining += sum(row.read_validation_required for row in value.technical_obligations)
                    continue
                identities = {row.identity_id: row for row in value.allele.identities}
                for call, obligation in zip(value.allele.calls or (), value.technical_obligations, strict=True):
                    fields = call.source.source_fields or {}
                    source_id = fields.get("source_allele_id")
                    original = identities.get(source_id)
                    occurrence = db.execute("SELECT input_sha,payload FROM occurrences WHERE id=?", (source_id,)).fetchone()
                    if (original is None or occurrence is None or fields.get("source_input_sha256") != occurrence[0]
                            or source_id != obligation.source_allele_id or call.sample is None
                            or call.sample.sample_name != obligation.sample_name
                            or call.sample.sample_name != sources[occurrence[0]].selected_sample_name
                            or call.source.source_sha256 != original.source.source_sha256):
                        raise WgsReadRescueError("merged original call binding changed")
                    occurrence_value = _typed(json.loads(occurrence[1]), IntakeOccurrence, decode=True)
                    payload_sha = fields.get("source_payload_sha256", call.source.source_sha256)
                    if (fields.get("source_record_id") != occurrence_value.source_record_id
                            or payload_sha != occurrence_value.payload_sha256
                            or "source_payload_sha256" in fields and (
                                call.source.source_sha256 != occurrence[0]
                                or call.source.source_release != "VCFv4.1")):
                        raise WgsReadRescueError("merged call differs from original source occurrence")
                    before = original.source.source_fields or {}
                    target = _bytes(_target(value.allele, call))
                    found = db.execute("SELECT payload FROM assessments WHERE target=?", (target,)).fetchone()
                    assessment = None if found is None else _typed(json.loads(found[0]), ReadAssessment, decode=True)
                    expected_fields = dict(before)
                    if assessment is not None:
                        if before.get("read_validation_required") is not True:
                            raise WgsReadRescueError("assessment does not address an original obligation")
                        expected_fields["authorized_read_assessment"] = json.loads(found[0])
                        if assessment.outcome is AssessmentOutcome.VALIDATED:
                            expected_fields["read_validation_required"] = False
                        db.execute("UPDATE assessments SET used=1 WHERE target=?", (target,))
                    if (_bytes(fields) != _bytes(expected_fields)
                            or obligation.read_validation_required != fields.get("read_validation_required")
                            or obligation.technical_flags != tuple(fields.get("technical_flags", ()))):
                        raise WgsReadRescueError("merged technical discharge changed")
                    remaining += obligation.read_validation_required
            if (remaining != self.intake.manifest.remaining_technical_obligations
                    or admitted != self.intake.manifest.admitted_source_occurrences
                    or db.execute("SELECT count(*) FROM assessments WHERE used=0").fetchone()[0]
                    or db.execute("SELECT count(*) FROM additions WHERE used=0").fetchone()[0]):
                raise WgsReadRescueError("merged intake projection is incomplete")
            for payload, in db.execute("SELECT payload FROM occurrences WHERE used=0"):
                occurrence = _typed(json.loads(payload), IntakeOccurrence, decode=True)
                if occurrence.early_filter_disposition not in (None, "exclude_common"):
                    raise WgsReadRescueError("merged intake omitted retained original evidence")
            db.commit()
            return self
        except BaseException:
            self.__exit__(None, None, None)
            raise

    def _stage(self, delta, stage):
        if any(row.row_count > self.config.source_shards.max_uncompressed_rows
               or row.size_bytes > self.config.source_shards.max_uncompressed_bytes for row in stage.shards):
            raise WgsReadRescueError("merged producer shard exceeds execution capacity", category="response_capacity")
        cursor = iter(self.connection.execute("SELECT key,payload FROM records WHERE delta_sha=? AND stage=? ORDER BY key",
                                            (delta.sha256, stage.stage_id)))
        stream, classes = hashlib.sha256(), Counter()
        additions = assessments = 0
        for shard in stage.shards:
            digest, size, first, last = hashlib.sha256(), 0, None, None
            for _ in range(shard.row_count):
                key, payload = next(cursor, (None, None))
                if payload is None:
                    raise WgsReadRescueError("merged intake omitted producer evidence")
                record = _typed(json.loads(payload), ReadDeltaRecord, decode=True)
                if record.variant_class not in next(row.variant_classes for row in delta.selected_plan.stages if row.stage_id == stage.stage_id):
                    raise WgsReadRescueError("merged producer class changed")
                first, last = key if first is None else first, key
                digest.update(payload)
                stream.update(payload)
                size += len(payload)
                classes[record.variant_class] += 1
                additions += record.addition is not None
                assessments += record.assessment is not None
            if (digest.hexdigest(), size, first, last) != (shard.sha256, shard.size_bytes, shard.first_key, shard.last_key):
                raise WgsReadRescueError("merged producer shard projection changed")
        if (next(cursor, None) is not None or stream.hexdigest() != stage.stream_sha256
                or (additions, assessments) != (stage.addition_count, stage.assessment_count)
                or tuple((kind, classes[kind]) for kind in sorted(WgsClass)) != stage.class_counts):
            raise WgsReadRescueError("merged producer record closure changed")

    def validate_allele(self, allele):
        row = self.connection.execute("SELECT payload FROM admissions WHERE id=?", (allele.allele_id,)).fetchone()
        if row is None:
            raise WgsReadRescueError("scientific allele is outside merged intake")
        value = _bind_intake_admission(load_ledger_admission_input(row[0]), self.intake.manifest)
        expected = value.allele
        if (expected is None or allele.native_identity_id != expected.native_identity_id
                or allele.identities != expected.identities or allele.calls != expected.calls):
            raise WgsReadRescueError("scientific allele differs from merged intake")

    def __exit__(self, *_):
        if self.connection is not None:
            self.connection.close()
            self.connection = None
        self.path.unlink(missing_ok=True)


def _bind_intake_admission(value, manifest):
    if value.allele is not None:
        samples = {row.input_sha256: row.selected_sample_name for row in manifest.source_bindings}
        samples.update((stage.output_artifact.sha256, delta.sample_id)
                       for delta in manifest.deltas for stage in delta.stage_results)
        def bind(source):
            fields = dict(source.source_fields or {})
            source_sample = samples.get(fields.get("source_input_sha256"))
            if source_sample is not None:
                fields["merged_intake_sample"] = {"manifest_sha256": manifest.sha256,
                    "sample_id": manifest.sample_id, "source_sample_name": source_sample}
            if "authorized_read_evidence" in fields:
                fields["authorized_read_evidence"] = {
                    **fields["authorized_read_evidence"], "merged_intake_sha256": manifest.sha256}
            return replace(source, source_fields=fields)
        value = replace(value, allele=replace(value.allele,
            identities=tuple(replace(item, source=bind(item.source)) for item in value.allele.identities),
            calls=tuple(replace(call, source=bind(call.source)) for call in value.allele.calls or ())))
    return replace(value, parent_source_manifest_sha256=manifest.sha256,
        parent_working_set_manifest_sha256=manifest.sha256,
        native_source_contract_sha256=hashlib.sha256(_bytes(manifest.source_bindings)).hexdigest())


@dataclass(frozen=True, slots=True, repr=False)
class MergedReadIntake:
    manifest: MergedIntakeManifest
    directory: Path

    # A merged intake remains a distinct run. Its validated original calls and
    # read assessment provenance are replayed by admissions(), never inserted
    # into a previously completed admission database.

    def admissions(self, config):
        ordinal = 0
        for row in iter_merged_intake(self.manifest, directory=self.directory, execution_config=config):
            if row.admission is not None:
                ordinal += 1
                yield replace(_bind_intake_admission(row.admission, self.manifest), retained_ordinal=ordinal)

    def retained(self, config, *, source_manifest, working_manifest):
        for value in self.admissions(config):
            sequence = None
            if value.allele is not None:
                native = next(row.identity for row in value.allele.identities
                              if row.identity_id == value.allele.native_identity_id)
                sequence = native.sequence if isinstance(native, (MitochondrialAllele, MosaicAllele)) else native
                if type(sequence) is not SequenceAllele:
                    sequence = None
            route = "sequence" if re.fullmatch(r"allele-[0-9a-f]{24}", value.identity_id) else "native"
            geometry = (None,) * 5
            if sequence is not None:
                from genosapiens.track1_submission import _variant_kind
                geometry = (sequence.contig, sequence.position, sequence.reference,
                            sequence.alternate, _variant_kind(sequence.reference, sequence.alternate))
            yield RetainedAllele(value.retained_ordinal, value.identity_id, route, value.source_allele_ids,
                source_manifest.sha256, working_manifest.sha256, *geometry,
                _bytes({"schema": "genosapiens.wgs-retained-allele", "version": 1,
                    "classification": "controlled", "llm_policy": "no-llm", "identity_id": value.identity_id,
                    "route": route, "read_intake_sha256": self.manifest.sha256}).decode("ascii"))
