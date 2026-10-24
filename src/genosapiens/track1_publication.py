"""Acyclic immutable publication for completed and non-ready Track 1 analyses."""

from __future__ import annotations

from dataclasses import asdict, dataclass
from enum import StrEnum
import hashlib
import json
import os
import stat
import tempfile
from pathlib import Path
import re
from typing import Literal, Mapping

from genosapiens.evidence_model import _typed
from genosapiens.evidence_serialization import _loads
from genosapiens.governance import ArtifactClassification
from genosapiens.llm_guard import LLMPolicy
from genosapiens.reference_concordance import concordance_report_as_json
from genosapiens.track1_handoff import PublishedHandoff
from genosapiens.track1_artifacts import (
    AnalysisMode,
    ArtifactBinding,
    Track1ArtifactError,
    _artifact,
    _controlled_path,
    _publish_immutable,
    _publish_immutable_validated,
)
from genosapiens.wgs_reporting import (
    PreparedBoundedReports, WgsReportingError, _verify_bounded_audit_metadata,
    _AuditStreamVerifier, _open_audit_payload,
)
from genosapiens.bounded_wgs_config import BoundedWgsExecutionConfig
from genosapiens.wgs_read_rescue import MergedReadIntake, MergedIntakeIndex
from genosapiens.wgs_handoff import HandoffPublicationAuthority
from genosapiens.wgs_failures import WgsFailure
from genosapiens.wgs_prepared_case import PreparedCasePublication
from genosapiens.wgs_scientific_artifact import ScientificArtifact, scientific_artifact_complete
from genosapiens.wgs_scientific_products import (
    ScientificOrderProduct, ScientificPrimaryProduct,
    ScientificSelectionProduct, open_scientific_products, product_descriptors, scientific_report_roots,
)



OUTPUT_MANIFEST_SCHEMA = "genosapiens.track1-output-manifest"
PUBLICATION_RECEIPT_SCHEMA = "genosapiens.track1-publication-receipt"
OUTPUT_MANIFEST_PATH = "manifests/output-manifest.json"
PUBLICATION_RECEIPT_PATH = "receipts/publication.json"
HANDOFF_RECEIPT_PATH = "receipts/handoff.json"


class FinalRole(StrEnum):
    EVIDENCE = "final-evidence-document"
    PROVENANCE = "scientific-evidence"
    HYPOTHESES = "hypothesis-generation"
    RANKING = "ranking"
    SELECTION = "submission-selection"
    EPCR = "epcr-assignment"
    CSV = "track1-csv"
    PREFLIGHT = "preflight"
    AUDIT_JSON = "audit-json"
    AUDIT_MARKDOWN = "audit-markdown"
    METHODS_MARKDOWN = "methods-markdown"
    PROFILE = "analysis-profile-manifest"
    READ_INTAKE = "authorized-read-intake"


FINAL_PATHS: Mapping[FinalRole, str] = {
    FinalRole.EVIDENCE: "analysis/evidence-document.json",
    FinalRole.PROVENANCE: "analysis/scientific-evidence-manifest.json",
    FinalRole.HYPOTHESES: "analysis/hypothesis-generation.json",
    FinalRole.RANKING: "analysis/ranking.json",
    FinalRole.SELECTION: "analysis/submission-selection.json",
    FinalRole.EPCR: "analysis/epcr.json",
    FinalRole.CSV: "submission/track1.csv",
    FinalRole.PREFLIGHT: "reports/preflight.json",
    FinalRole.AUDIT_JSON: "reports/audit.json",
    FinalRole.AUDIT_MARKDOWN: "reports/audit.md",
    FinalRole.METHODS_MARKDOWN: "reports/methods.md",
    FinalRole.PROFILE: "analysis/profile/manifest.json",
    FinalRole.READ_INTAKE: "analysis/read-intake/merged-intake-manifest.json",
}
FINAL_ORDER = tuple(FinalRole)
SCIENTIFIC_ORDER = FINAL_ORDER[:6]
READY_ROLES = frozenset(FINAL_ORDER) - {FinalRole.READ_INTAKE}
SHA256 = re.compile(r"^[0-9a-f]{64}$")


class Track1PublicationError(RuntimeError):
    """Static publication failure; controlled content and paths remain undisclosed."""


def _canonical_bytes(value: object) -> bytes:
    return (
        json.dumps(
            value,
            sort_keys=True,
            separators=(",", ":"),
            ensure_ascii=True,
            allow_nan=False,
        )
        + "\n"
    ).encode("ascii")


def _sha256(content: bytes) -> str:
    return hashlib.sha256(content).hexdigest()


@dataclass(frozen=True, slots=True, repr=False)
class FinalArtifact:
    role: FinalRole
    content: bytes

    def __post_init__(self) -> None:
        if type(self.role) is not FinalRole or type(self.content) is not bytes or not self.content:
            raise Track1PublicationError("invalid final analysis artifact")

    @property
    def relative_path(self) -> str:
        return FINAL_PATHS[self.role]

    @property
    def binding(self) -> ArtifactBinding:
        return _artifact(self.role.value, self.relative_path, self.content)

    def __repr__(self) -> str:
        return "<controlled final artifact; no-llm>"


@dataclass(frozen=True, slots=True)
class FinalStageState:
    role: FinalRole
    status: Literal["available", "unavailable"]
    artifact_sha256: str | None

    def __post_init__(self) -> None:
        if type(self.role) is not FinalRole or (
            (self.status == "available") != (self.artifact_sha256 is not None)
        ) or self.artifact_sha256 is not None and SHA256.fullmatch(
            self.artifact_sha256
        ) is None:
            raise Track1PublicationError("invalid final stage state")


@dataclass(frozen=True, slots=True)
class FinalOutputManifest:
    run_id: str
    input_manifest_sha256: str
    early_snapshot_sha256: str
    mode: AnalysisMode
    outputs: tuple[ArtifactBinding, ...]
    stages: tuple[FinalStageState, ...]
    competitive_status: Literal["ready", "not-ready"]
    schema: Literal["genosapiens.track1-output-manifest"] = OUTPUT_MANIFEST_SCHEMA
    classification: ArtifactClassification = ArtifactClassification.CONTROLLED
    llm_policy: LLMPolicy = LLMPolicy.NO_LLM

    def __post_init__(self) -> None:
        for value in (self.input_manifest_sha256, self.early_snapshot_sha256):
            if SHA256.fullmatch(value) is None:
                raise Track1PublicationError("invalid final output manifest")
        available = tuple(state.role for state in self.stages if state.status == "available")
        if (
            not self.run_id.startswith("track1-run:")
            or type(self.mode) is not AnalysisMode
            or self.schema != OUTPUT_MANIFEST_SCHEMA
            or self.classification is not ArtifactClassification.CONTROLLED
            or self.llm_policy is not LLMPolicy.NO_LLM
            or tuple(state.role for state in self.stages) != FINAL_ORDER
            or tuple(item.role for item in self.outputs)
            != tuple(role.value for role in available)
            or tuple(item.relative_path for item in self.outputs)
            != tuple(FINAL_PATHS[role] for role in available)
            or any(
                state.artifact_sha256 != output.sha256
                for state, output in zip(
                    (state for state in self.stages if state.status == "available"),
                    self.outputs,
                    strict=True,
                )
            )
            or (self.competitive_status == "ready") != (READY_ROLES <= set(available))
            or self.competitive_status == "ready" and self.mode is not AnalysisMode.PRIMARY
        ):
            raise Track1PublicationError("invalid final output manifest")

    def as_bytes(self) -> bytes:
        return _canonical_bytes(asdict(self))

    @property
    def sha256(self) -> str:
        return _sha256(self.as_bytes())


@dataclass(frozen=True, slots=True)
class PublicationReceipt:
    run_id: str
    input_manifest_sha256: str
    output_manifest_sha256: str
    primary_sha256: str
    scientific_sha256: str
    competitive_status: Literal["ready", "not-ready"]
    optional_analysis_status: Literal["not_run"] = "not_run"
    schema: Literal[
        "genosapiens.track1-publication-receipt"
    ] = PUBLICATION_RECEIPT_SCHEMA
    classification: ArtifactClassification = ArtifactClassification.CONTROLLED
    llm_policy: LLMPolicy = LLMPolicy.NO_LLM

    def __post_init__(self) -> None:
        if (
            not self.run_id.startswith("track1-run:")
            or any(
                SHA256.fullmatch(value) is None
                for value in (
                    self.input_manifest_sha256,
                    self.output_manifest_sha256,
                    self.primary_sha256,
                    self.scientific_sha256,
                )
            )
            or self.competitive_status not in ("ready","not-ready")
            or self.optional_analysis_status != "not_run"
            or self.schema != PUBLICATION_RECEIPT_SCHEMA
            or self.classification is not ArtifactClassification.CONTROLLED
            or self.llm_policy is not LLMPolicy.NO_LLM
        ):
            raise Track1PublicationError("invalid publication receipt")

    def as_bytes(self) -> bytes:
        return _canonical_bytes(asdict(self))


@dataclass(frozen=True, slots=True)
class HandoffPublicationReceipt:
    run_id: str
    input_manifest_sha256: str
    output_manifest_sha256: str
    handoff_sha256: str
    handoff_relative_path: str
    schema: Literal[
        "genosapiens.track1-handoff-publication-receipt"
    ] = "genosapiens.track1-handoff-publication-receipt"
    classification: ArtifactClassification = ArtifactClassification.CONTROLLED
    llm_policy: LLMPolicy = LLMPolicy.NO_LLM

    def __post_init__(self) -> None:
        if (
            not self.run_id.startswith("track1-run:")
            or any(
                SHA256.fullmatch(value) is None
                for value in (
                    self.input_manifest_sha256,
                    self.output_manifest_sha256,
                    self.handoff_sha256,
                )
            )
            or self.handoff_relative_path
            != "controlled-work/track1/handoff/mechanism-ready-handoff.json"
            or self.schema != "genosapiens.track1-handoff-publication-receipt"
            or self.classification is not ArtifactClassification.CONTROLLED
            or self.llm_policy is not LLMPolicy.NO_LLM
        ):
            raise Track1PublicationError("invalid publication receipt")

    def as_bytes(self) -> bytes:
        return _canonical_bytes(asdict(self))


@dataclass(frozen=True, slots=True, repr=False)
class FinalPublication:
    run_root: Path
    output_manifest: FinalOutputManifest
    handoff: PublishedHandoff | None
    receipt: PublicationReceipt
    reused_outputs: bool

    def __repr__(self) -> str:
        return "<controlled final publication; no-llm>"


@dataclass(frozen=True, slots=True, repr=False)
class PrimaryPublicationAuthority:
    source: object
    prefilter: object
    working: object
    case_manifest_sha256: str
    concordance_sha256: str
    run_evidence: object
    public_reference: object
    evaluator_source_path: Path | None = None

    def __post_init__(self):
        from genosapiens.wgs_ledger import LedgerManifest,FilterManifest,WorkingSetManifest
        from genosapiens.track1_preflight import CompetitiveRunEvidence
        from genosapiens.admission_links import PinnedPublicSource
        if (type(self.source) is not LedgerManifest or type(self.prefilter) is not FilterManifest
                or type(self.working) is not WorkingSetManifest or type(self.run_evidence) is not CompetitiveRunEvidence
                or type(self.public_reference) is not PinnedPublicSource
                or SHA256.fullmatch(self.case_manifest_sha256) is None
                or SHA256.fullmatch(self.concordance_sha256) is None):
            raise Track1PublicationError("primary publication authority is invalid")


@dataclass(frozen=True, slots=True, repr=False)
class BoundedPublicationInputs:
    """Authenticated sharded leaves used by the bounded publication path."""

    scientific: ScientificArtifact
    scientific_directory: Path
    order: ScientificOrderProduct
    order_directory: Path
    primary: ScientificPrimaryProduct
    selection: ScientificSelectionProduct
    selection_directory: Path
    reports: PreparedBoundedReports
    execution_config: BoundedWgsExecutionConfig
    read_intake: MergedReadIntake | None = None
    eligibility_root: Path | None = None
    profile_root: Path | None = None

    def __post_init__(self) -> None:
        if (
            not scientific_artifact_complete(self.scientific)
            or type(self.execution_config) is not BoundedWgsExecutionConfig
            or type(self.order) is not ScientificOrderProduct
            or type(self.primary) is not ScientificPrimaryProduct
            or type(self.selection) is not ScientificSelectionProduct
            or type(self.reports) is not PreparedBoundedReports
            or any(not value.is_absolute() for value in (
                self.scientific_directory,
                self.order_directory,
                self.selection_directory,
            ))
        ):
            raise Track1PublicationError("bounded publication inputs are invalid")

    def __repr__(self) -> str:
        return "<controlled bounded publication inputs; no-llm>"


def load_final_output_manifest(content: bytes) -> FinalOutputManifest:
    """Decode canonical saved metadata, without loading scientific leaves."""
    try:
        value = _typed(_loads(content), FinalOutputManifest, decode=True)
        if value.as_bytes() != content:
            raise ValueError
        return value
    except (ValueError, TypeError, RuntimeError, AttributeError):
        raise Track1PublicationError("invalid canonical final output metadata") from None


def load_publication_receipt(content: bytes) -> PublicationReceipt:
    """Decode a canonical receipt; this creates no scientific capability."""
    try:
        value = _typed(_loads(content), PublicationReceipt, decode=True)
        if value.as_bytes() != content:
            raise ValueError
        return value
    except (ValueError, TypeError, RuntimeError, AttributeError):
        raise Track1PublicationError("invalid canonical publication receipt") from None


def check_published_inventory(root: Path, *, expected_output_manifest_sha256: str) -> dict:
    """Bounded final packaging check after normal validated publication.

    Hash small role metadata/CSV and check retained payload paths and sizes.
    Payload hashes and scientific correctness rely on the producing publication;
    this deliberately does not open a science reader or rehash large shards.
    """
    from genosapiens.wgs_shards import read_bounded_metadata
    from genosapiens.wgs_scientific_artifact import scientific_artifact_from_bytes, scientific_artifact_complete
    from genosapiens.wgs_reporting import audit_manifest_from_json, diagnostic_from_json
    from genosapiens.wgs_profile_projection import scientific_profile_component, profile_retained_files
    from genosapiens.wgs_association_projection import assessment_retained_files
    from genosapiens.wgs_analysis_eligibility import scientific_eligibility_component, eligibility_retained_files
    import csv
    import io
    try:
        root = Path(root)
        if not root.is_absolute() or root.is_symlink():
            raise ValueError
        raw = read_bounded_metadata(root / OUTPUT_MANIFEST_PATH)
        if _sha256(raw) != expected_output_manifest_sha256:
            raise ValueError
        manifest = load_final_output_manifest(raw)
        receipt = load_publication_receipt(read_bounded_metadata(root / PUBLICATION_RECEIPT_PATH))
        if (receipt.output_manifest_sha256 != expected_output_manifest_sha256
                or receipt.run_id != manifest.run_id or receipt.input_manifest_sha256 != manifest.input_manifest_sha256
                or receipt.competitive_status != manifest.competitive_status):
            raise ValueError
        roles = {}
        for output in manifest.outputs:
            content = read_bounded_metadata(root / output.relative_path)
            if (_sha256(content), len(content)) != (output.sha256, output.size_bytes):
                raise ValueError
            roles[FinalRole(output.role)] = content
        required = {FinalRole.EVIDENCE, FinalRole.PROVENANCE, FinalRole.HYPOTHESES,
            FinalRole.RANKING, FinalRole.SELECTION, FinalRole.PROFILE, FinalRole.AUDIT_JSON, FinalRole.AUDIT_MARKDOWN}
        if not required <= roles.keys():
            raise ValueError
        artifact = scientific_artifact_from_bytes(roles[FinalRole.HYPOTHESES], expected_sha256=receipt.scientific_sha256)
        profile = scientific_profile_component(artifact, root / 'analysis')
        audit = audit_manifest_from_json(roles[FinalRole.AUDIT_JSON], expected_sha256=_sha256(roles[FinalRole.AUDIT_JSON]))
        diagnostic = diagnostic_from_json(roles[FinalRole.EVIDENCE], expected_sha256=_sha256(roles[FinalRole.EVIDENCE]))
        primary_raw = read_bounded_metadata(root / 'analysis/primary-manifest.json')
        if (_sha256(primary_raw) != receipt.primary_sha256 or diagnostic.primary_manifest_sha256 != receipt.primary_sha256
                or diagnostic.scientific_manifest_sha256 != artifact.sha256
                or roles[FinalRole.PROVENANCE] != artifact.evidence.as_bytes()
                or roles[FinalRole.AUDIT_MARKDOWN] != audit.markdown_index()
                or audit.diagnostic_sha256 != diagnostic.sha256):
            raise ValueError
        checked_paths = 0
        def retained(relative, size):
            nonlocal checked_paths
            path = root / relative
            if (not relative or Path(relative).is_absolute() or '..' in Path(relative).parts
                    or any(p.is_symlink() for p in (path, *path.parents))):
                raise ValueError
            info = path.lstat()
            if not stat.S_ISREG(info.st_mode) or info.st_size != size:
                raise ValueError
            checked_paths += 1
        for shard in (*artifact.evidence.node_shards, *artifact.evidence.root_shards):
            retained('analysis/' + shard.relative_path, shard.size_bytes)
        for _, relative, _, size in profile_retained_files(root, profile):
            retained(relative, size)
        for _, relative, _, size in assessment_retained_files(root, dict(artifact.parent_bindings)):
            retained(relative, size)
        eligibility = scientific_eligibility_component(artifact, root / 'analysis')
        if eligibility is not None:
            for _, relative, _, size in eligibility_retained_files(root, eligibility):
                retained(relative, size)
        for shard in (*audit.json_shards, *audit.markdown_shards):
            retained('reports/audit/' + shard.relative_path, shard.size_bytes)
        csv_rows = None
        if FinalRole.CSV in roles:
            from genosapiens.track1_preflight import COLUMNS
            rows = list(csv.reader(io.StringIO(roles[FinalRole.CSV].decode('utf-8'))))
            if not rows or tuple(rows[0]) != COLUMNS or len(rows) > 11 or any(len(row) != len(COLUMNS) for row in rows[1:]):
                raise ValueError
            csv_rows = len(rows) - 1
            if manifest.competitive_status == 'ready' and not csv_rows:
                raise ValueError
        return dict(status='passed', verification_scope='bounded final inventory/receipt/metadata/CSV and retained path sizes',
            output_manifest_sha256=expected_output_manifest_sha256, scientific_manifest_sha256=artifact.sha256,
            profile_manifest_sha256=profile['manifest_sha256'], profile_counts=profile['counts'],
            competitive_status=manifest.competitive_status, csv_rows=csv_rows, retained_paths_checked=checked_paths,
            counts=dict(contexts=artifact.context_count, positive_contexts=artifact.positive_context_count,
                context_members=artifact.member_count, logical_pairs=artifact.pair_count,
                hypotheses=artifact.hypothesis_count, selected_groups=diagnostic.selected_group_count),
            scientific_analysis_repeated=False, deep_product_readback=False, payloads_rehashed=False)
    except (ValueError, TypeError, RuntimeError, AttributeError, KeyError, OSError, UnicodeError):
        raise Track1PublicationError('final published inventory check failed') from None


def load_handoff_receipt(content: bytes) -> HandoffPublicationReceipt:
    """Read the separate optional receipt without changing primary completion."""
    try:
        value = _typed(_loads(content), HandoffPublicationReceipt, decode=True)
        if value.as_bytes() != content:
            raise ValueError
        return value
    except (ValueError, TypeError, RuntimeError, AttributeError):
        raise Track1PublicationError("invalid canonical handoff receipt") from None


def _validate_dependency_prefix(roles: tuple[FinalRole, ...]) -> None:
    available = set(roles)
    for index, role in enumerate(SCIENTIFIC_ORDER):
        if role in available and not set(SCIENTIFIC_ORDER[:index]) <= available:
            raise Track1PublicationError("final scientific stage lacks its ancestors")
    reports = {
        FinalRole.PREFLIGHT,
        FinalRole.AUDIT_JSON,
        FinalRole.AUDIT_MARKDOWN,
        FinalRole.METHODS_MARKDOWN,
    }
    if available & reports and FinalRole.SELECTION not in available:
        raise Track1PublicationError("final report lacks submission selection")
    if FinalRole.CSV in available and not {
        FinalRole.EPCR,
        FinalRole.PREFLIGHT,
    } <= available:
        raise Track1PublicationError("competitive CSV lacks EPCR or preflight")


def _validate_artifact_envelopes(artifacts: tuple[FinalArtifact, ...]) -> None:
    for item in artifacts:
        if item.relative_path.endswith(".json"):
            try:
                document = json.loads(item.content)
                if (
                    not isinstance(document, dict)
                    or document.get("classification") != "controlled"
                    or document.get("llm_policy") != "no-llm"
                    or _canonical_bytes(document) != item.content
                ):
                    raise ValueError
            except (UnicodeError, ValueError, json.JSONDecodeError):
                raise Track1PublicationError(
                    "final JSON artifact is not canonical controlled evidence"
                ) from None
        elif item.role is FinalRole.CSV:
            try:
                text = item.content.decode("ascii", errors="strict")
            except UnicodeError:
                raise Track1PublicationError("competitive CSV is not ASCII") from None
            if not text.endswith("\n") or "," not in text.splitlines()[0]:
                raise Track1PublicationError("competitive CSV is malformed")
        else:
            try:
                text = item.content.decode("utf-8", errors="strict")
            except UnicodeError:
                raise Track1PublicationError("final report is not UTF-8") from None
            if (
                not text.endswith("\n")
                or "controlled" not in text.lower()
                or "no-llm" not in text.lower()
            ):
                raise Track1PublicationError(
                    "final report lacks controlled no-LLM labeling"
                )


def _bounded_manifest(
    parent: PreparedCasePublication,
    artifacts: tuple[FinalArtifact, ...],
) -> FinalOutputManifest:
    if type(parent) is not PreparedCasePublication:
        raise Track1PublicationError("prepared-case authority is required")
    by_role = {item.role: item for item in artifacts}
    if len(by_role) != len(artifacts):
        raise Track1PublicationError("duplicate final analysis artifact")
    _validate_artifact_envelopes(artifacts)
    ordered = tuple(by_role[role] for role in FINAL_ORDER if role in by_role)
    _validate_dependency_prefix(tuple(item.role for item in ordered))
    stages = tuple(
        FinalStageState(
            role,
            "available" if role in by_role else "unavailable",
            by_role[role].binding.sha256 if role in by_role else None,
        )
        for role in FINAL_ORDER
    )
    return FinalOutputManifest(
        run_id="track1-run:bounded-" + parent.manifest.sha256[:20],
        input_manifest_sha256=parent.manifest.sha256,
        early_snapshot_sha256=parent.early_snapshot.sha256,
        mode=AnalysisMode.PRIMARY,
        outputs=tuple(item.binding for item in ordered),
        stages=stages,
        competitive_status=(
            "ready" if READY_ROLES <= set(by_role) else "not-ready"
        ),
    )


def _validate_bounded_publication(
    value: BoundedPublicationInputs,
    artifacts: tuple[FinalArtifact, ...],
    *, evidence_parents, _scientific_reader=None,
) -> None:
    """Replay all sharded scientific and provenance leaves before publication."""

    try:
        from genosapiens.wgs_shards import read_bounded_metadata
        evidence_content = read_bounded_metadata(value.scientific_directory / "scientific-evidence-manifest.json")
        if evidence_content != value.scientific.evidence.as_bytes():
            raise ValueError
        with open_scientific_products(value.scientific, product_descriptors(value),
                directory=value.scientific_directory, expected_parent_bindings=evidence_parents,
                scratch_directory=value.scientific_directory.parent,eligibility_root=value.eligibility_root,
                profile_root=value.profile_root, _reader=_scientific_reader) as reader:
            if reader.resources.analysis_profile is None:
                raise ValueError
            if value.reports.audit.scientific_roots != scientific_report_roots(reader, product_descriptors(value)):
                raise ValueError
            domain=reader.resources.assessment_domain
            domain_sha256=reader.resources.assessment_domain_sha256
            diagnostic=value.reports.diagnostic
            if (diagnostic.assessment_domain_sha256!=domain_sha256
                    or value.reports.audit.assessment_domain_sha256!=domain_sha256
                    or diagnostic.assessment_input_count!=domain['original_inputs']
                    or diagnostic.assessment_input_dispositions!=(None if domain['input_dispositions'] is None else
                        tuple(sorted(domain['input_dispositions'].items())))):
                raise ValueError
            if value.read_intake is not None:
                with MergedIntakeIndex(value.read_intake, execution_config=value.execution_config,
                        scratch_path=value.scientific_directory / "publication-intake.sqlite") as intake_index:
                    for row in reader.iter_member_records():
                        intake_index.validate_allele(row["member"].allele)
        if (
            value.order.scientific_manifest_sha256 != value.scientific.sha256
            or value.primary.scientific_manifest_sha256 != value.scientific.sha256
            or value.primary.order_manifest_sha256 != value.order.sha256
            or value.selection.scientific_manifest_sha256 != value.scientific.sha256
            or value.selection.order_manifest_sha256 != value.order.sha256
        ):
            raise ValueError
        report = value.reports
        audit = report.audit
        diagnostic = report.diagnostic
        parents = dict(audit.parent_bindings)
        if (diagnostic.source_manifest_sha256 != evidence_parents["source_manifest"]
                or diagnostic.filter_manifest_sha256 != evidence_parents["prefilter_manifest"]
                or diagnostic.working_manifest_sha256 != evidence_parents["working_set"]
                or diagnostic.scientific_manifest_sha256 != value.scientific.sha256
                or diagnostic.order_manifest_sha256 != value.order.sha256
                or diagnostic.primary_manifest_sha256 != value.primary.sha256
                or diagnostic.selection_sha256 != value.selection.sha256
                or diagnostic.ranked_hypothesis_count != value.scientific.hypothesis_count
                or diagnostic.selected_group_count != len(value.selection.selected)
                or audit.representation_rows != value.scientific.member_count
                or audit.hypothesis_rows != value.scientific.hypothesis_count
                or audit.source_rows != diagnostic.source_record_count
                or audit.audit_rows != (audit.source_accounting_rows + audit.representation_rows
                                        + audit.detailed_hypothesis_rows + 1)
                or audit.detailed_hypothesis_rows != len(value.selection.selected)
                or any(parents.get(key) != digest for key, digest in (
                    ("source_manifest", evidence_parents["source_manifest"]),
                    ("filter_manifest", evidence_parents["prefilter_manifest"]),
                    ("working_manifest", evidence_parents["working_set"]),
                    ("scientific_manifest", value.scientific.sha256),
                    ("scientific_evidence", value.scientific.evidence.sha256),
                    ("order_manifest", value.order.sha256),
                    ("primary_manifest", value.primary.sha256),
                    ("selection", value.selection.sha256),
                    *((("read_intake", value.read_intake.manifest.sha256),)
                      if value.read_intake else ()),
                ))):
            raise ValueError
        _verify_bounded_audit_metadata(value.reports)
        expected = {
            FinalRole.EVIDENCE: value.reports.diagnostic.as_bytes(),
            FinalRole.PROVENANCE: evidence_content,
            FinalRole.HYPOTHESES: value.scientific.as_bytes(),
            FinalRole.RANKING: value.order.as_bytes(),
            FinalRole.SELECTION: value.selection.as_bytes(),
            FinalRole.AUDIT_JSON: value.reports.audit_json_index,
            FinalRole.AUDIT_MARKDOWN: value.reports.audit_markdown_index,
        }
        if value.read_intake is not None:
            expected[FinalRole.READ_INTAKE] = value.read_intake.manifest.as_bytes()
        from genosapiens.wgs_profile_projection import scientific_profile_component
        component=scientific_profile_component(value.scientific,value.scientific_directory,
            profile_root=value.profile_root)
        if component is None:
            raise ValueError
        profile_root=value.scientific_directory.parent/'analysis/profile' if value.profile_root is None else value.profile_root
        expected[FinalRole.PROFILE]=read_bounded_metadata(profile_root/'manifest.json')
        actual = {item.role: item.content for item in artifacts}
        if any(actual.get(role) != content for role, content in expected.items()):
            raise ValueError
    except WgsFailure:
        raise
    except (ValueError, TypeError, RuntimeError, OSError):
        raise Track1PublicationError(
            "bounded final publication replay failed"
        ) from None


def _publish_immutable_report_shard(
    source: Path,
    destination: Path,
    *,
    expected_sha256: str,
    expected_size: int,
    root: Path,
    audit_verifier: _AuditStreamVerifier | None = None,
    copy_progress=None,
) -> bool:
    """Copy once, or authenticate the reused destination once.

    Audit checks share these chunks. Reuse authenticates destination bytes, not
    unread source payloads. Earlier private copies may survive a late failure.
    """

    try:
        source_metadata = source.lstat()
        if (source.is_symlink() or not stat.S_ISREG(source_metadata.st_mode)
                or source_metadata.st_size != expected_size
                or any(path.is_symlink() for path in (*source.parents, *destination.parents))):
            raise ValueError
        if audit_verifier is not None and (
                source_metadata.st_uid != os.geteuid() or source_metadata.st_nlink != 1
                or stat.S_IMODE(source_metadata.st_mode) != 0o600):
            raise ValueError
        destination.parent.mkdir(mode=0o700, parents=True, exist_ok=True)
        destination.resolve().relative_to(root.resolve())
        if destination.exists():
            metadata = destination.lstat()
            if (
                destination.is_symlink()
                or not stat.S_ISREG(metadata.st_mode)
                or stat.S_IMODE(metadata.st_mode) != 0o400
                or metadata.st_size != expected_size
            ):
                raise ValueError
            digest = hashlib.sha256()
            size = 0
            reader = (_open_audit_payload(destination, expected_size=expected_size, mode=0o400)
                      if audit_verifier is not None else destination.open("rb"))
            with reader as handle:
                for content in iter(lambda: handle.read(65536), b""):
                    digest.update(content)
                    size += len(content)
                    if copy_progress is not None:
                        try: copy_progress(len(content))
                        except Exception: pass
                    if size > expected_size:
                        raise ValueError
                    if audit_verifier is not None:
                        audit_verifier.update(content)
            actual_sha256 = digest.hexdigest()
            if size != expected_size or actual_sha256 != expected_sha256:
                raise ValueError
            if audit_verifier is not None:
                audit_verifier.finish_shard(sha256=actual_sha256, size_bytes=size)
            return True
        temporary = destination.with_name(destination.name + ".partial")
        descriptor = os.open(
            temporary,
            os.O_WRONLY | os.O_CREAT | os.O_EXCL | os.O_NOFOLLOW,
            0o600,
        )
        digest = hashlib.sha256()
        size = 0
        try:
            reader = (_open_audit_payload(source, expected_size=expected_size)
                      if audit_verifier is not None else source.open("rb"))
            with os.fdopen(descriptor, "wb") as output, reader as input_:
                for content in iter(lambda: input_.read(65536), b""):
                    if size + len(content) > expected_size:
                        raise ValueError
                    output.write(content)
                    digest.update(content)
                    size += len(content)
                    if copy_progress is not None:
                        try: copy_progress(len(content))
                        except Exception: pass
                    if audit_verifier is not None:
                        audit_verifier.update(content)
                output.flush()
                os.fsync(output.fileno())
            actual_sha256 = digest.hexdigest()
            if size != expected_size or actual_sha256 != expected_sha256:
                raise ValueError
            if audit_verifier is not None:
                audit_verifier.finish_shard(sha256=actual_sha256, size_bytes=size)
            temporary.chmod(0o400)
            os.replace(temporary, destination)
        finally:
            try:
                temporary.unlink()
            except FileNotFoundError:
                pass
        return False
    except (ValueError, OSError, WgsReportingError):
        raise Track1PublicationError(
            "bounded report shard publication failed"
        ) from None


def _qualification_artifact_path(root: Path, path: Path) -> None:
    """Public qualification outputs stay beneath their already bound run root."""
    _public_artifact_path(root, path, prefix=("qualification-final",))


def _public_artifact_path(root: Path, path: Path, *, prefix: tuple[str, ...]) -> None:
    """Confine one public run without allowing another run or any symlink."""
    if not root.is_absolute() or root.resolve(strict=True) != root:
        raise ValueError
    relative = path.relative_to(root)
    if (
        ".." in relative.parts
        or not relative.parts
        or relative.parts[:len(prefix)] != prefix
        or any(part in {"controlled-work", "controlled", "data"} for part in relative.parts)
    ):
        raise ValueError
    current = root
    for part in relative.parts:
        current /= part
        if (
            current.is_symlink()
            or current.exists() and current != path and not current.is_dir()
        ):
            raise ValueError


def publish_bounded_final_analysis(
    parent: PreparedCasePublication,
    artifacts: tuple[FinalArtifact, ...],
    *,
    bounded_inputs: BoundedPublicationInputs,
    reference_root: str | Path,
    run_root: str | Path,
    base_release_bytes: bytes | None,
    handoff_authority: HandoffPublicationAuthority | None = None,
    primary_authority: PrimaryPublicationAuthority | None = None,
    public_qualification_authority=None,
    public_validation_authority_bytes: bytes | None = None,
    completed_analysis_input=None,
    expected_base_release_sha256: str | None = None,
    progress=None,
    _scientific_reader=None,
) -> FinalPublication:
    """Publish bounded outputs with original preparation and actual release parents.

    Only resolved public completed inputs permit a distinct consumer release;
    they never authorize controlled handoff or change the prepared manifest.
    """

    if (
        type(parent) is not PreparedCasePublication
        or type(artifacts) is not tuple
        or type(bounded_inputs) is not BoundedPublicationInputs
        or bounded_inputs.scientific.parent_bindings
        and ("prefilter_manifest", parent.manifest.prefilter_manifest.sha256)
        not in bounded_inputs.scientific.parent_bindings
    ):
        raise Track1PublicationError("bounded prepared publication inputs are invalid")
    from genosapiens.track1_installation import load_run_versions as load_base_release
    from genosapiens.wgs_shards import read_bounded_metadata
    evidence = bounded_inputs.scientific.evidence
    validation = None
    if parent.manifest.public_validation_authority_sha256 is not None:
        from genosapiens.wgs_validation_authority import load_public_validation_authority
        if (base_release_bytes is not None or public_qualification_authority is not None
                or handoff_authority is not None or primary_authority is not None or completed_analysis_input is not None):
            raise Track1PublicationError("public validation cannot authorize competitive publication")
        validation = load_public_validation_authority(public_validation_authority_bytes,
            expected_sha256=parent.manifest.public_validation_authority_sha256)
        policy_sha256 = dict(bounded_inputs.scientific.parent_bindings)["policy_bindings"]
        if policy_sha256 not in {row.sha256 for row in validation.candidate_bindings}:
            raise Track1PublicationError("scientific policy is outside public validation authority")
        release = None
    else:
        if public_validation_authority_bytes is not None:
            raise Track1PublicationError("public validation authority differs from prepared case")
        if completed_analysis_input is not None:
            from genosapiens.wgs_analysis_inputs import ResolvedCompletedAnalysisInput
            if (type(completed_analysis_input) is not ResolvedCompletedAnalysisInput
                    or public_qualification_authority is None or handoff_authority is not None or primary_authority is not None
                    or bounded_inputs.read_intake is not None or expected_base_release_sha256 is None):
                raise Track1PublicationError("completed input publication authority is invalid")
            release = load_base_release(base_release_bytes, expected_sha256=expected_base_release_sha256)
            completed_analysis_input.require_consumer(public_qualification_authority, release, prepared=parent)
        else:
            from genosapiens.track1_installation import InstalledVersions
            release = load_base_release(base_release_bytes,
                expected_sha256=expected_base_release_sha256 or parent.manifest.base_release_sha256)
            if (type(release) is not InstalledVersions
                    and release.sha256 != parent.manifest.base_release_sha256
                    or public_qualification_authority is not None and any(
                        row.role == "completed_analysis_input" for row in public_qualification_authority.resource_bindings)):
                raise Track1PublicationError("prepared release differs from publication")
            if (release.reference_manifest_sha256 != parent.completed_source.manifest.source.reference_manifest_sha256
                    or release.policy_bindings.common_allele_policy_sha256 != parent.manifest.prefilter_policy_sha256
                    or release.policy_bindings.bounded_wgs_execution_sha256 != parent.manifest.execution_limits_sha256):
                raise Track1PublicationError('publication preparation dependencies differ')
        policy_sha256 = release.policy_bindings.sha256
    from genosapiens.wgs_analysis_eligibility import scientific_eligibility_component,eligibility_component_parents
    eligibility=scientific_eligibility_component(bounded_inputs.scientific,bounded_inputs.scientific_directory,eligibility_root=getattr(bounded_inputs,"eligibility_root",None))
    from genosapiens.wgs_profile_projection import scientific_profile_component,profile_component_parents,profile_retained_files,verify_analysis_profile_component
    profile=scientific_profile_component(bounded_inputs.scientific,bounded_inputs.scientific_directory,
        profile_root=bounded_inputs.profile_root)
    if profile is None:
        raise Track1PublicationError('required analysis profile is missing')
    from genosapiens.wgs_association_projection import verify_assessment_projections,assessment_parents,assessment_retained_files
    assessments=verify_assessment_projections(bounded_inputs.scientific_directory.parent,
        dict(bounded_inputs.scientific.parent_bindings),full=False)
    from genosapiens.wgs_runtime import analysis_progress
    from genosapiens.wgs_work_contract import publication_work_bindings
    work_bindings = publication_work_bindings(dict(bounded_inputs.scientific.parent_bindings),
        (release if release is not None else validation).resource_manifest_sha256s)
    analysis_progress(progress, 'publication',
        'final_scientific_replay' if _scientific_reader is None else 'publication_closure')
    _validate_bounded_publication(bounded_inputs, artifacts, evidence_parents={
        "prepared_case": parent.manifest.sha256,
        "source_manifest": parent.manifest.source_manifest.sha256,
        "prefilter_manifest": parent.manifest.prefilter_manifest.sha256,
        "working_set": parent.manifest.working_sets[0].sha256,
        "execution_config": parent.manifest.execution_limits_sha256,
        "policy_bindings": policy_sha256,
        **({"read_intake": bounded_inputs.read_intake.manifest.sha256} if bounded_inputs.read_intake else {}),
        **({} if eligibility is None else eligibility_component_parents(eligibility)),
        **profile_component_parents(profile),
        **({} if assessments is None else assessment_parents(*assessments)),
        **work_bindings,
    }, _scientific_reader=_scientific_reader)
    analysis_progress(progress, 'publication', 'final_artifacts')
    output_manifest = _bounded_manifest(parent, artifacts)
    if handoff_authority is not None:
        if type(handoff_authority) is not HandoffPublicationAuthority:
            raise Track1PublicationError("handoff publication authority is invalid")
        from genosapiens.case_manifest import load_case_manifest, case_manifest_as_bytes
        derived=PrimaryPublicationAuthority(**{name:getattr(handoff_authority,name)
            for name in PrimaryPublicationAuthority.__dataclass_fields__
            if name not in {'case_manifest_sha256', 'concordance_sha256'}},
            case_manifest_sha256=hashlib.sha256(case_manifest_as_bytes(
                load_case_manifest(handoff_authority.case_manifest_path))).hexdigest(),
            concordance_sha256=hashlib.sha256(handoff_authority.concordance_content).hexdigest())
        if primary_authority is not None and primary_authority!=derived:
            raise Track1PublicationError("primary and handoff authorities differ")
        primary_authority=derived
    if (output_manifest.competitive_status == "ready") != (primary_authority is not None):
        raise Track1PublicationError("competitive primary authority is incomplete")
    if primary_authority is not None:
        _validate_primary_authority(parent=parent,inputs=bounded_inputs,artifacts={item.role:item for item in artifacts},
            authority=primary_authority,release=release,output_manifest=output_manifest,
            scientific_reader=_scientific_reader, public_qualification_authority=public_qualification_authority)
    try:
        if public_qualification_authority is None and validation is None:
            root, destination = _controlled_path(reference_root, run_root)
            publish = _publish_immutable
        else:
            from genosapiens.wgs_public_qualification import (
                PublicQualificationAuthority,
            )

            if validation is None and type(public_qualification_authority) is not PublicQualificationAuthority:
                raise ValueError
            root = Path(reference_root).resolve(strict=True)
            destination = Path(run_root)
            expected_root = Path((validation or public_qualification_authority).output_root)
            expected_destination = (expected_root / "runs" / bounded_inputs.scientific.sha256
                                    if validation is not None else expected_root / "qualification-final")
            if (
                root != expected_root.resolve(strict=True)
                or destination != expected_destination
                or destination.is_symlink()
                or not destination.resolve(strict=False).is_relative_to(root)
                or any(part in {"controlled-work", "controlled", "data"}
                       for part in destination.parts)
            ):
                raise ValueError
            def validate_public_path(root, path):
                _public_artifact_path(root, path, prefix=(
                    ("runs", bounded_inputs.scientific.sha256) if validation is not None
                    else ("qualification-final",)))

            validate_public_path(root, destination)

            def publish(path, content, root):
                return _publish_immutable_validated(
                    path, content, root, validate_public_path,
                )
        destination.mkdir(mode=0o700, parents=True, exist_ok=True)
        if (
            destination.is_symlink()
            or stat.S_IMODE(destination.lstat().st_mode) != 0o700
        ):
            raise ValueError
        from genosapiens.wgs_runtime import analysis_progress
        files_completed = bytes_completed = 0
        def observe_copy(amount):
            nonlocal bytes_completed
            bytes_completed += amount
            analysis_progress(progress, 'publication', 'output_copy',
                counts={'files_completed':files_completed,'bytes_completed':bytes_completed},
                unit='bytes',completed=bytes_completed)
        def copy_shard(*args, **kwargs):
            nonlocal files_completed
            result = _publish_immutable_report_shard(*args,copy_progress=observe_copy,**kwargs)
            files_completed += 1
            observe_copy(0)
            return result
        reused = []
        for shard in (*evidence.node_shards, *evidence.root_shards):
            reused.append(copy_shard(
                bounded_inputs.scientific_directory / shard.relative_path,
                destination / "analysis" / shard.relative_path,
                expected_sha256=shard.sha256, expected_size=shard.size_bytes, root=root))
        if eligibility is not None:
            from genosapiens.wgs_analysis_eligibility import eligibility_retained_files
            for source,relative,digest,size in eligibility_retained_files(bounded_inputs.scientific_directory.parent,eligibility,eligibility_root=getattr(bounded_inputs,"eligibility_root",None)):
                reused.append(copy_shard(source,destination/relative,
                    expected_sha256=digest,expected_size=size,root=root))
        for source,relative,digest,size in profile_retained_files(bounded_inputs.scientific_directory.parent,
                profile,profile_root=bounded_inputs.profile_root):
            if relative==FINAL_PATHS[FinalRole.PROFILE]:
                continue  # Published once, as the required final metadata role below.
            reused.append(copy_shard(source,destination/relative,
                expected_sha256=digest,expected_size=size,root=root))
        for source,relative,digest,size in assessment_retained_files(bounded_inputs.scientific_directory.parent,
                dict(bounded_inputs.scientific.parent_bindings)):
            reused.append(copy_shard(source,destination/relative,
                expected_sha256=digest,expected_size=size,root=root))
        if bounded_inputs.read_intake is not None:
            for shard in bounded_inputs.read_intake.manifest.shards:
                reused.append(copy_shard(
                    bounded_inputs.read_intake.directory / shard.relative_path,
                    destination / "analysis/read-intake" / shard.relative_path,
                    expected_sha256=shard.sha256, expected_size=shard.size_bytes, root=root))
        audit_verifier = _AuditStreamVerifier(bounded_inputs.reports)
        for shard in (
            *bounded_inputs.reports.audit.json_shards,
            *bounded_inputs.reports.audit.markdown_shards,
        ):
            audit_verifier.begin_shard(shard)
            if public_qualification_authority is not None or validation is not None:
                validate_public_path(
                    root, destination / "reports" / "audit" / shard.relative_path,
                )
            reused.append(copy_shard(
                bounded_inputs.reports.output_directory / shard.relative_path,
                destination / "reports" / "audit" / shard.relative_path,
                expected_sha256=shard.sha256,
                expected_size=shard.size_bytes,
                root=root,
                audit_verifier=audit_verifier,
            ))
        audit_verifier.finish()
        reused.append(publish(destination / "analysis/primary-manifest.json",
                              bounded_inputs.primary.as_bytes(), root))
        by_role = {item.role: item for item in artifacts}
        for role in FINAL_ORDER:
            if role in by_role:
                item = by_role[role]
                reused.append(publish(
                    destination / item.relative_path, item.content, root,
                ))
        verify_analysis_profile_component(destination,profile,full=False,
            scratch_directory=bounded_inputs.scientific_directory.parent)
        # A shared operation owner must close successfully before any final
        # publication marker. Separate/checkpoint consumers reopen normally.
        if _scientific_reader is not None:
            analysis_progress(progress, 'publication', 'publication_closure', force=True)
            _scientific_reader.close()
        receipt=PublicationReceipt(output_manifest.run_id,parent.manifest.sha256,output_manifest.sha256,
            bounded_inputs.primary.sha256,bounded_inputs.scientific.sha256,output_manifest.competitive_status)
        reused.append(publish(destination / OUTPUT_MANIFEST_PATH,output_manifest.as_bytes(),root))
        reused.append(publish(destination / PUBLICATION_RECEIPT_PATH,receipt.as_bytes(),root))
        published_handoff = None
        if handoff_authority is not None:
            published_handoff, optional_receipt, handoff_reused = _publish_current_handoff(
                parent=parent, inputs=bounded_inputs, artifacts=by_role,
                authority=handoff_authority, release=release,
                output_manifest=output_manifest, destination=destination,
                root=root, publish=publish)
            reused.extend(handoff_reused)
        return FinalPublication(
            destination,
            output_manifest,
            published_handoff,
            receipt,
            all(reused),
        )
    except Track1PublicationError:
        raise
    except (OSError, ValueError, Track1ArtifactError, WgsReportingError):
        raise Track1PublicationError(
            "bounded final publication failed"
        ) from None


def _validate_primary_authority(*,parent,inputs,artifacts,authority,release,output_manifest,
        scientific_reader=None, public_qualification_authority=None):
    from genosapiens.track1_epcr import EpcrAssignment
    from genosapiens.track1_preflight import ReleasePreflight,prepare_bounded_release_submission
    if type(authority) is not PrimaryPublicationAuthority or release is None:
        raise Track1PublicationError("competitive primary authority is invalid")
    if public_qualification_authority is not None:
        from genosapiens.wgs_public_qualification import require_public_qualification_release
        require_public_qualification_release(public_qualification_authority, release)
        if (authority.source.source.input_sha256 != public_qualification_authority.input_sha256
                or authority.source.source.input_size_bytes != public_qualification_authority.input_size_bytes
                or parent.manifest.hpo_input_sha256 != public_qualification_authority.hpo_input_sha256
                or parent.manifest.execution_limits_sha256 != public_qualification_authority.execution_config_sha256):
            raise Track1PublicationError('primary inputs differ from public qualification')
    run=authority.run_evidence
    if (authority.source.sha256!=parent.manifest.source_manifest.sha256
            or authority.prefilter.sha256!=parent.manifest.prefilter_manifest.sha256
            or authority.working.sha256!=parent.manifest.working_sets[0].sha256
            or authority.case_manifest_sha256!=parent.manifest.case_manifest_sha256
            or authority.case_manifest_sha256!=authority.source.source.case_manifest_sha256
            or authority.concordance_sha256!=parent.manifest.concordance_sha256
            or authority.concordance_sha256!=authority.source.source.concordance_sha256
            or run.run_id!=output_manifest.run_id or run.input_manifest_sha256!=parent.manifest.sha256
            or run.early_snapshot_sha256!=parent.early_snapshot.sha256
            or run.primary_snapshot_sha256!=inputs.primary.sha256
            or run.release_manifest_sha256!=release.sha256
            or run.accepted_reference_manifest_sha256!=authority.source.source.reference_manifest_sha256
            or run.concordance_report_sha256!=parent.manifest.concordance_sha256):
        raise Track1PublicationError("primary run and original input bindings differ")
    assignment = _typed(json.loads(artifacts[FinalRole.EPCR].content)["result"], EpcrAssignment, decode=True)
    if assignment.as_bytes() != artifacts[FinalRole.EPCR].content:
        raise Track1PublicationError("competitive EPCR envelope changed")
    preflight = _typed(json.loads(artifacts[FinalRole.PREFLIGHT].content), ReleasePreflight, decode=True)
    submission = prepare_bounded_release_submission(inputs.selection, assignment,
        authority.run_evidence, inputs.scientific, inputs.order,
        scientific_directory=inputs.scientific_directory,
        order_directory=inputs.order_directory, selection_directory=inputs.selection_directory,
        release_manifest_bytes=release.as_bytes(), expected_release_manifest_sha256=release.sha256,
        evaluator_source_path=authority.evaluator_source_path,
        public_reference=authority.public_reference, expected_preflight=preflight,
        read_intake_manifest=None if inputs.read_intake is None else inputs.read_intake.manifest,
        eligibility_root=inputs.eligibility_root, profile_root=inputs.profile_root,
        scientific_reader=scientific_reader, public_qualification_authority=public_qualification_authority)
    if submission.csv_bytes != artifacts[FinalRole.CSV].content:
        raise Track1PublicationError("competitive primary CSV binding changed")


def _publish_current_handoff(*, parent, inputs, artifacts, authority, release,
                             output_manifest, destination, root, publish):
    """Replay readiness and the sharded consumer before sealing publication."""
    from genosapiens.bounded_wgs_config import bounded_wgs_execution_config_as_json
    from genosapiens.case_manifest import load_case_manifest
    from genosapiens.track1_epcr import EpcrAssignment
    from genosapiens.track1_preflight import ReleasePreflight, prepare_bounded_release_submission
    from genosapiens.wgs_handoff import (
        HANDOFF_SCHEMA, HandoffArtifact, validate_handoff_bundle, write_handoff_manifest,
    )
    if type(authority) is not HandoffPublicationAuthority:
        raise Track1PublicationError("current handoff authority is required")
    assignment = _typed(json.loads(artifacts[FinalRole.EPCR].content)["result"], EpcrAssignment, decode=True)
    if assignment.as_bytes() != artifacts[FinalRole.EPCR].content:
        raise Track1PublicationError("competitive EPCR envelope changed")
    preflight = _typed(json.loads(artifacts[FinalRole.PREFLIGHT].content), ReleasePreflight, decode=True)
    submission = prepare_bounded_release_submission(inputs.selection, assignment,
        authority.run_evidence, inputs.scientific, inputs.order,
        scientific_directory=inputs.scientific_directory,
        order_directory=inputs.order_directory, selection_directory=inputs.selection_directory,
        release_manifest_bytes=release.as_bytes(), expected_release_manifest_sha256=release.sha256,
        evaluator_source_path=authority.evaluator_source_path,
        public_reference=authority.public_reference, expected_preflight=preflight,
        read_intake_manifest=None if inputs.read_intake is None else inputs.read_intake.manifest,
        eligibility_root=inputs.eligibility_root, profile_root=inputs.profile_root)
    if submission.csv_bytes != artifacts[FinalRole.CSV].content:
        raise Track1PublicationError("competitive handoff CSV binding changed")
    schema_content = authority.schema_path.read_bytes()
    contract_content = authority.contract_path.read_bytes()
    contract = json.loads(contract_content)
    canonical_path = contract["canonical_artifact_path"]
    if (contract["canonical_schema"] != HANDOFF_SCHEMA
            or hashlib.sha256(schema_content).hexdigest() != release.canonical_handoff_schema_sha256
            or hashlib.sha256(contract_content).hexdigest() != release.path_selection_policy_sha256
            or hashlib.sha256(_canonical_bytes(canonical_path)).hexdigest()
            != release.canonical_handoff_path_sha256):
        raise Track1PublicationError("handoff schema or path policy changed")
    from genosapiens.wgs_shards import read_bounded_metadata
    case_content = read_bounded_metadata(authority.case_manifest_path)
    case = load_case_manifest(authority.case_manifest_path)
    run = authority.run_evidence
    metadata = {
        "prepared": ("parents/prepared.json", parent.manifest.as_bytes()),
        "early": ("parents/early.json", parent.early_snapshot.as_bytes()),
        "source": ("parents/source.json", authority.source.as_bytes()),
        "source_accounting": ("parents/source-accounting.json",
            parent.completed_source.accounting.as_bytes()),
        "source_audit": ("parents/source-audit.json",
            parent.completed_source.audit.as_bytes()),
        "preparation_run_options": ("parents/preparation-run-options.json",
            read_bounded_metadata(parent.root.parent / "run-options.json")),
        "completed_source_input": ("parents/completed-source.json", parent.completed_source.descriptor.as_bytes()),
        "prepared_index": ("parents/prepared-index.json",
            read_bounded_metadata(parent.root.parent / "state/prepared-index.json")),
        "prefilter": ("parents/prefilter.json", authority.prefilter.as_bytes()),
        "working": ("parents/working.json", authority.working.as_bytes()),
        "execution": ("parents/execution.json", bounded_wgs_execution_config_as_json(inputs.execution_config)),
        "release": ("parents/release.json", release.as_bytes()),
        "case": ("parents/case.json", case_content),
        "concordance": ("parents/concordance.json", authority.concordance_content),
        "run_evidence": ("analysis/run-evidence.json", run.as_bytes()),
        "primary": ("analysis/primary-manifest.json", inputs.primary.as_bytes()),
        "output": (OUTPUT_MANIFEST_PATH, output_manifest.as_bytes()),
    }
    if inputs.read_intake is not None:
        metadata["read_intake"] = (FINAL_PATHS[FinalRole.READ_INTAKE], inputs.read_intake.manifest.as_bytes())
    for role, final_role in (("evidence", FinalRole.PROVENANCE), ("scientific", FinalRole.HYPOTHESES),
            ("order", FinalRole.RANKING), ("selection", FinalRole.SELECTION),
            ("preflight", FinalRole.PREFLIGHT), ("epcr", FinalRole.EPCR), ("csv", FinalRole.CSV),
            ("analysis_profile", FinalRole.PROFILE)):
        metadata[role] = (FINAL_PATHS[final_role], artifacts[final_role].content)
    bindings = tuple(HandoffArtifact(role, path, hashlib.sha256(content).hexdigest(), len(content))
                     for role, (path, content) in sorted(metadata.items()))
    run_declaration = {
        "run_id": run.run_id, "prepared_manifest_sha256": parent.manifest.sha256,
        "source_input_sha256": authority.source.source.input_sha256,
        "reference_identity_id": run.accepted_reference_identity_id,
        "reference_manifest_sha256": run.accepted_reference_manifest_sha256,
        "concordance_report_sha256": run.concordance_report_sha256,
        "configuration_sha256": case.configuration.sha256,
        "image_digest": release.runtime_image_digest,
        "release_manifest_sha256": release.sha256,
        "output_manifest_sha256": output_manifest.sha256,
    }
    evidence = inputs.scientific.evidence
    if evidence.as_bytes() != metadata["evidence"][1]:
        raise Track1PublicationError("handoff evidence declaration changed")
    with tempfile.TemporaryDirectory(prefix="handoff-stage-", dir=inputs.scientific_directory.parent) as temporary:
        stage = Path(temporary)
        for path, content in metadata.values():
            target = stage / path
            target.parent.mkdir(mode=0o700, parents=True, exist_ok=True)
            descriptor = os.open(target, os.O_WRONLY | os.O_CREAT | os.O_EXCL, 0o600)
            with os.fdopen(descriptor, "wb") as handle:
                handle.write(content)
        # Link already authenticated scientific leaves into a private rehearsal
        # bundle. Publication below streams verified bytes into immutable files.
        for source_dir, role, shards in (
            (inputs.scientific_directory, "scientific", (*evidence.node_shards, *evidence.root_shards)),
            *(((inputs.read_intake.directory, "read_intake", inputs.read_intake.manifest.shards),)
                if inputs.read_intake is not None else ()),
        ):
            target_dir = (stage / metadata[role][0]).parent
            for shard in shards:
                os.link(source_dir / shard.relative_path, target_dir / shard.relative_path)
        handoff_dir = (stage / canonical_path).parent
        from genosapiens.wgs_analysis_eligibility import scientific_eligibility_component,eligibility_retained_files
        eligibility=scientific_eligibility_component(inputs.scientific,inputs.scientific_directory,eligibility_root=getattr(inputs,"eligibility_root",None))
        if eligibility is not None:
            scope=(stage/metadata['scientific'][0]).parent.parent
            for source,relative,digest,size in eligibility_retained_files(inputs.scientific_directory.parent,eligibility,eligibility_root=getattr(inputs,"eligibility_root",None)):
                target=scope/relative
                target.parent.mkdir(mode=0o700,parents=True,exist_ok=True)
                os.link(source,target)
        from genosapiens.wgs_profile_projection import scientific_profile_component,profile_retained_files
        profile=scientific_profile_component(inputs.scientific,inputs.scientific_directory,
            profile_root=inputs.profile_root)
        if profile is None:
            raise Track1PublicationError('required handoff analysis profile is missing')
        scope=(stage/metadata['scientific'][0]).parent.parent
        for source,relative,digest,size in profile_retained_files(inputs.scientific_directory.parent,
                profile,profile_root=inputs.profile_root):
            if relative == FINAL_PATHS[FinalRole.PROFILE]:
                # The exact authenticated bytes are already a metadata role.
                continue
            target=scope/relative
            target.parent.mkdir(mode=0o700,parents=True,exist_ok=True)
            os.link(source,target)
        from genosapiens.wgs_association_projection import assessment_retained_files
        scope=(stage/metadata['scientific'][0]).parent.parent
        for source,relative,digest,size in assessment_retained_files(inputs.scientific_directory.parent,
                dict(inputs.scientific.parent_bindings)):
            target=scope/relative
            target.parent.mkdir(mode=0o700,parents=True,exist_ok=True)
            os.link(source,target)
        with open_scientific_products(inputs.scientific, product_descriptors(inputs),
                directory=inputs.scientific_directory,
                expected_parent_bindings=dict(inputs.scientific.parent_bindings),
                scratch_directory=stage,eligibility_root=inputs.eligibility_root,
                profile_root=inputs.profile_root) as reader:
            handoff = write_handoff_manifest(inputs, reader=reader, directory=handoff_dir,
                execution_config=inputs.execution_config, producer_version=authority.application_version,
                schema_sha256=hashlib.sha256(schema_content).hexdigest(),
                contract_sha256=hashlib.sha256(contract_content).hexdigest(),
                canonical_artifact_path=canonical_path, run=run_declaration, artifacts=bindings)
            validate_handoff_bundle(handoff, root=stage, directory=handoff_dir,
                scratch_path=stage / "verify.sqlite", validated_reader=reader)
        reused = []
        for path, content in metadata.values():
            reused.append(publish(destination / path, content, root))
        reused.append(publish(destination / canonical_path, handoff.as_bytes(), root))
    receipt = HandoffPublicationReceipt(
        run_id=output_manifest.run_id,
        input_manifest_sha256=output_manifest.input_manifest_sha256,
        output_manifest_sha256=output_manifest.sha256,
        handoff_sha256=handoff.sha256,
        handoff_relative_path=canonical_path,
    )
    reused.append(publish(destination / HANDOFF_RECEIPT_PATH, receipt.as_bytes(), root))
    return PublishedHandoff(destination / canonical_path, handoff.sha256, len(handoff.as_bytes())), receipt, reused
