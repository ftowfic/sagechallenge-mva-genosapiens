"""Validation-only Track 1 CSV boundary, never submission readiness or publishing.

This module renders only source-bound task 11.2 selections with matching task 11.3
assignments. Tasks 9.6, 11.6, and 13.1 must supply the missing frozen scope, run,
authorization, concordance-replay, and release contracts before a ready preflight
can exist.
"""

from __future__ import annotations

import csv
from contextlib import nullcontext
from dataclasses import dataclass
from decimal import Decimal, InvalidOperation
import hashlib
from io import StringIO
import math
from pathlib import Path
import tempfile
from typing import Literal, Mapping

from genosapiens.challenge_evaluator import load_pinned_evaluator
from genosapiens.admission_links import PinnedPublicSource
from genosapiens.evidence_model import Completeness, _Controlled, _exact, _ids, _sha, _typed
from genosapiens.evidence_serialization import _canonical
from genosapiens.track1_epcr import CONFIDENCE_LABEL, EpcrAssignment
from genosapiens.track1_installation import load_run_versions as load_base_release
from genosapiens.track1_submission import (
    PUBLIC_REFERENCE_IDENTITY,
    SubmissionAlleleTuple,
    SubmissionSelection,
    verify_sharded_selected_submission_reference,
    verify_selected_submission_reference,
)
from genosapiens.wgs_scientific_artifact import ScientificArtifact, ScientificArtifactReader, scientific_artifact_complete
from genosapiens.wgs_scientific_products import ScientificOrderProduct, ScientificSelectionProduct, scientific_products



SCHEMA = "genosapiens.track1-validation-only-preflight"
FINAL_SCHEMA = "genosapiens.track1-preflight"
COLUMNS = (
    "proband_id",
    "chrom_1",
    "pos_1",
    "ref_1",
    "alt_1",
    "chrom_2",
    "pos_2",
    "ref_2",
    "alt_2",
    "epcr",
    "finding_type",
    "notes",
)
EVALUATOR_PIN = Path(__file__).resolve().parents[2] / "config/track1/challenge-evaluator.json"
_BASE_BLOCKERS = (
    "accepted-reference-and-current-concordance-replay-unavailable",
    "base-release-freeze-unavailable",
    "frozen-assessment-scope-unavailable",
    "immutable-run-snapshot-unavailable",
    "production-epcr-policy-not-frozen",
)
_ERRORS = (
    ValueError,
    TypeError,
    RuntimeError,
    AttributeError,
    KeyError,
    IndexError,
    InvalidOperation,
    OSError,
    OverflowError,
    RecursionError,
)


class Track1PreflightError(ValueError):
    """Static diagnostics only; candidate rows and paths remain controlled."""


def _bytes(value: object) -> bytes:
    return (_canonical(value) + "\n").encode("utf-8")


@dataclass(frozen=True, slots=True, repr=False)
class Track1CsvRow(_Controlled):
    hypothesis_id: str
    first: SubmissionAlleleTuple
    second: SubmissionAlleleTuple | None
    epcr: str
    notes: str = ""

    def __post_init__(self):
        _Controlled.__post_init__(self)
        _ids((self.hypothesis_id,))
        value = _exact(self.epcr)
        if (not 0 < value <= 1 or self.second is not None and self.first == self.second
                or self.notes not in ("", "Incomplete recessive partial member; partner unresolved")
                or self.notes and self.second is not None):
            raise Track1PreflightError("invalid Track 1 CSV row")

    @property
    def identity(self) -> tuple[SubmissionAlleleTuple, ...]:
        return (self.first,) if self.second is None else (self.first, self.second)

    def fields(self) -> tuple[object, ...]:
        second = ("", "", "", "") if self.second is None else (
            self.second.contig,
            self.second.position,
            self.second.reference,
            self.second.alternate,
        )
        return (
            "PROBAND01",
            self.first.contig,
            self.first.position,
            self.first.reference,
            self.first.alternate,
            *second,
            self.epcr,
            "primary",
            self.notes,
        )


@dataclass(frozen=True, slots=True, repr=False)
class EvaluatorParserValidation(_Controlled):
    pin_sha256: str
    revision: str
    source_sha256: str
    source_git_blob_sha1: str
    source_size_bytes: int
    status: Literal["passed"] = "passed"
    private_truth_status: Literal["unavailable-not-evaluated"] = "unavailable-not-evaluated"

    def __post_init__(self):
        _Controlled.__post_init__(self)
        _sha(self.pin_sha256)
        _ids((self.revision,))
        _sha(self.source_sha256)
        if (
            len(self.source_git_blob_sha1) != 40
            or any(character not in "0123456789abcdef" for character in self.source_git_blob_sha1)
            or self.source_size_bytes <= 0
        ):
            raise Track1PreflightError("invalid evaluator parser validation")


@dataclass(frozen=True, slots=True, repr=False)
class ValidationOnlyPreflight(_Controlled):
    selection_sha256: str
    epcr_assignment_sha256: str
    rows: tuple[Track1CsvRow, ...]
    csv_sha256: str | None
    evaluator: EvaluatorParserValidation | None
    blockers: tuple[str, ...]
    warnings: tuple[str, ...]
    selected_read_evidence_present: bool
    rescue_status: Literal["run-status-unavailable", "selected-read-evidence-present-unvalidated"]
    rescue_authorization_status: Literal[
        "run-snapshot-unavailable", "task-11.6-validation-unavailable"
    ]
    confidence_label: Literal[
        "uncalibrated-relative-confidence-not-clinical-probability"
    ] = CONFIDENCE_LABEL
    schema: Literal["genosapiens.track1-validation-only-preflight"] = SCHEMA
    validation_only: Literal[True] = True
    ready: Literal[False] = False
    private_truth_status: Literal["unavailable-not-evaluated"] = "unavailable-not-evaluated"

    def __post_init__(self):
        _Controlled.__post_init__(self)
        _sha(self.selection_sha256)
        _sha(self.epcr_assignment_sha256)
        _ids(tuple(item.hypothesis_id for item in self.rows))
        _ids(self.blockers)
        _ids(self.warnings)
        if (
            not self.blockers
            or self.blockers != tuple(sorted(self.blockers))
            or self.warnings != tuple(sorted(self.warnings))
            or not 0 <= len(self.rows) <= 10
            or (self.csv_sha256 is None) != (not self.rows)
            or self.csv_sha256 is not None
            and len(self.csv_sha256) != 64
            or any(left <= right for left, right in zip(
                tuple(_exact(item.epcr) for item in self.rows),
                tuple(_exact(item.epcr) for item in self.rows)[1:],
            ))
            or self.selected_read_evidence_present
            and (
                self.rescue_status != "selected-read-evidence-present-unvalidated"
                or self.rescue_authorization_status != "task-11.6-validation-unavailable"
            )
            or not self.selected_read_evidence_present
            and (
                self.rescue_status != "run-status-unavailable"
                or self.rescue_authorization_status != "run-snapshot-unavailable"
            )
        ):
            raise Track1PreflightError("invalid validation-only preflight")
        if self.csv_sha256 is not None:
            _sha(self.csv_sha256)

    def as_bytes(self) -> bytes:
        return _bytes({
            "schema": self.schema,
            "classification": "controlled",
            "llm_policy": "no-llm",
            "validation_only": True,
            "ready": False,
            "private_truth_status": self.private_truth_status,
            "preflight": self,
        })


@dataclass(frozen=True, slots=True, repr=False)
class PreparedValidationOnlySubmission(_Controlled):
    csv_bytes: bytes | None
    preflight: ValidationOnlyPreflight

    def __post_init__(self):
        _Controlled.__post_init__(self)
        if (
            (self.csv_bytes is None) != (self.preflight.csv_sha256 is None)
            or self.csv_bytes is not None
            and hashlib.sha256(self.csv_bytes).hexdigest() != self.preflight.csv_sha256
        ):
            raise Track1PreflightError("validation-only CSV differs from preflight")


@dataclass(frozen=True, slots=True, repr=False)
class SelectedTechnicalAssessment(_Controlled):
    """One sample/source/native-allele obligation outcome with actual evidence."""

    allele_id: str
    source_allele_id: str
    sample_name: str
    source_input_sha256: str
    reason: str
    scope: str
    status: Literal["not_applicable", "satisfied", "unmet"]
    measurement_sha256: str | None

    def __post_init__(self) -> None:
        _Controlled.__post_init__(self)
        _ids((
            self.allele_id,
            self.source_allele_id,
            self.sample_name,
            self.reason,
            self.scope,
        ))
        _sha(self.source_input_sha256)
        if self.measurement_sha256 is not None:
            _sha(self.measurement_sha256)
        if (
            self.status == "satisfied" and self.measurement_sha256 is None
            or self.status == "not_applicable"
            and (
                self.measurement_sha256 is not None
                or self.reason != "no-mandatory-technical-obligation"
            )
        ):
            raise Track1PreflightError("invalid selected technical assessment")


@dataclass(frozen=True, slots=True, repr=False)
class CompetitiveRunEvidence(_Controlled):
    """Hash-bound final-run states consumed by release-aware preflight."""

    run_id: str
    input_manifest_sha256: str
    early_snapshot_sha256: str
    primary_snapshot_sha256: str
    selection_sha256: str
    epcr_assignment_sha256: str
    accepted_reference_status: Literal["accepted", "failed", "not_evaluated"]
    accepted_reference_identity_id: str
    accepted_reference_manifest_sha256: str
    concordance_status: Literal[
        "complete", "failed", "blocked_order_adjudication", "not_evaluated"
    ]
    concordance_report_sha256: str
    assessment_scope_sha256: str
    completely_enumerated_hypothesis_ids: tuple[str, ...]
    technical_assessments: tuple[SelectedTechnicalAssessment, ...]
    additional_class_status: Literal["assessed", "not_assessed"]
    rescue_recommendations: tuple[str, ...]
    rescue_execution_status: Literal["not_performed", "performed"]
    rescue_authorization_status: Literal["not_applicable", "authorized"]
    rescue_authorization_sha256: str | None
    rescue_plan_status: Literal["not_applicable", "complete", "incomplete"]
    release_manifest_sha256: str
    contradiction_warning_sha256s: tuple[str, ...]

    def __post_init__(self) -> None:
        _Controlled.__post_init__(self)
        _ids((self.run_id, self.accepted_reference_identity_id))
        _ids(self.completely_enumerated_hypothesis_ids)
        _ids(self.rescue_recommendations)
        for value in (
            self.input_manifest_sha256,
            self.early_snapshot_sha256,
            self.primary_snapshot_sha256,
            self.selection_sha256,
            self.epcr_assignment_sha256,
            self.accepted_reference_manifest_sha256,
            self.concordance_report_sha256,
            self.assessment_scope_sha256,
            self.release_manifest_sha256,
            *self.contradiction_warning_sha256s,
        ):
            _sha(value)
        if self.rescue_authorization_sha256 is not None:
            _sha(self.rescue_authorization_sha256)
        if (
            not self.run_id.startswith("track1-run")
            or self.completely_enumerated_hypothesis_ids
            != tuple(sorted(set(self.completely_enumerated_hypothesis_ids)))
            or any(type(item) is not SelectedTechnicalAssessment for item in self.technical_assessments)
            or self.technical_assessments
            != tuple(sorted(
                set(self.technical_assessments),
                key=lambda item: (
                    item.allele_id,
                    item.source_allele_id,
                    item.sample_name,
                    item.reason,
                    item.scope,
                ),
            ))
            or self.rescue_recommendations
            != tuple(sorted(set(self.rescue_recommendations)))
            or self.contradiction_warning_sha256s
            != tuple(sorted(set(self.contradiction_warning_sha256s)))
            or self.rescue_execution_status == "not_performed"
            and (
                self.rescue_authorization_status != "not_applicable"
                or self.rescue_authorization_sha256 is not None
                or self.rescue_plan_status != "not_applicable"
            )
            or self.rescue_execution_status == "performed"
            and (
                self.rescue_authorization_status != "authorized"
                or self.rescue_authorization_sha256 is None
                or self.rescue_plan_status == "not_applicable"
            )
        ):
            raise Track1PreflightError("invalid competitive run evidence")

    def as_bytes(self) -> bytes:
        return _bytes(self)


@dataclass(frozen=True, slots=True, repr=False)
class ReleasePreflight(_Controlled):
    run_evidence_sha256: str
    release_manifest_sha256: str
    selection_sha256: str
    epcr_assignment_sha256: str
    rows: tuple[Track1CsvRow, ...]
    draft_csv_sha256: str | None
    competitive_csv_sha256: str | None
    evaluator: EvaluatorParserValidation | None
    blockers: tuple[str, ...]
    warnings: tuple[str, ...]
    selected_read_evidence_present: bool
    schema: Literal["genosapiens.track1-preflight"] = FINAL_SCHEMA
    classification: Literal["controlled"] = "controlled"
    llm_policy: Literal["no-llm"] = "no-llm"
    finding_type: Literal["primary"] = "primary"
    private_truth_status: Literal[
        "unavailable-not-evaluated"
    ] = "unavailable-not-evaluated"

    def __post_init__(self) -> None:
        _Controlled.__post_init__(self)
        for value in (
            self.run_evidence_sha256,
            self.release_manifest_sha256,
            self.selection_sha256,
            self.epcr_assignment_sha256,
        ):
            _sha(value)
        for value in (self.draft_csv_sha256, self.competitive_csv_sha256):
            if value is not None:
                _sha(value)
        if (
            self.blockers != tuple(sorted(set(self.blockers)))
            or self.warnings != tuple(sorted(set(self.warnings)))
            or (self.draft_csv_sha256 is None) != (not self.rows)
            or (self.competitive_csv_sha256 is not None) != (not self.blockers)
            or self.competitive_csv_sha256 is not None
            and self.competitive_csv_sha256 != self.draft_csv_sha256
        ):
            raise Track1PreflightError("invalid release preflight")

    @property
    def ready(self) -> bool:
        return not self.blockers

    def as_bytes(self) -> bytes:
        return _bytes(self)


@dataclass(frozen=True, slots=True, repr=False)
class PreparedReleaseSubmission(_Controlled):
    csv_bytes: bytes | None
    preflight: ReleasePreflight

    def __post_init__(self) -> None:
        _Controlled.__post_init__(self)
        if (
            (self.csv_bytes is not None) != self.preflight.ready
            or self.csv_bytes is not None
            and hashlib.sha256(self.csv_bytes).hexdigest()
            != self.preflight.competitive_csv_sha256
        ):
            raise Track1PreflightError("competitive CSV differs from release preflight")


def _render_csv(rows: tuple[Track1CsvRow, ...]) -> bytes:
    stream = StringIO(newline="")
    writer = csv.writer(stream, lineterminator="\n")
    writer.writerow(COLUMNS)
    writer.writerows(item.fields() for item in rows)
    return stream.getvalue().encode("ascii")


def _validate_selected_trace_bindings(selection: SubmissionSelection) -> None:
    assessments = {item.allele_id: item for item in selection.representation_assessments}
    reference = selection.reference
    for group in selection.selected:
        for trace in group.alleles:
            assessment = assessments.get(trace.allele_id)
            if (
                assessment is None
                or assessment.native_identity_id != trace.native_identity_id
                or assessment.status != "eligible"
                or assessment.trace != trace
                or trace.submission_reference_role != reference.role
                or trace.submission_reference_identity_id != reference.identity_id
                or trace.submission_reference_manifest_sha256 != reference.manifest_sha256
                or trace.contig_aliases_sha256 != reference.contig_aliases_sha256
            ):
                raise Track1PreflightError("selected allele trace is not bound to its assessment and reference")


def _validate_evaluator_parser(
    csv_bytes: bytes,
    rows: tuple[Track1CsvRow, ...],
    source_path: str | Path,
) -> EvaluatorParserValidation:
    evaluator = load_pinned_evaluator(source_path, pin_path=EVALUATOR_PIN)
    if not callable(getattr(evaluator, "load_submission", None)):
        raise Track1PreflightError("pinned evaluator parser API is unavailable")
    with tempfile.TemporaryDirectory() as directory:
        path = Path(directory) / "validation-only-submission.csv"
        path.write_bytes(csv_bytes)
        parsed = evaluator.load_submission(str(path))
    if not isinstance(parsed, Mapping) or set(parsed) != {"PROBAND01"}:
        raise Track1PreflightError("pinned evaluator parsed an unexpected case set")
    records = tuple(parsed["PROBAND01"])
    if len(records) != len(rows):
        raise Track1PreflightError("pinned evaluator parsed an unexpected row count")
    evaluator_epcr = tuple(record.epcr for record in records)
    if (
        any(type(value) is not float or not math.isfinite(value) or not 0 < value <= 1 for value in evaluator_epcr)
        or any(left <= right for left, right in zip(evaluator_epcr, evaluator_epcr[1:]))
    ):
        raise Track1PreflightError("pinned evaluator did not preserve strict EPCR ordering")
    for rank, (row, record) in enumerate(zip(rows, records, strict=True), 1):
        identity = frozenset((
            item.contig,
            item.position,
            item.reference,
            item.alternate,
        ) for item in row.identity)
        if (
            record.variants != identity
            or record.epcr != float(row.epcr)
            or record.rank != rank
            or record.finding_type != "primary"
        ):
            raise Track1PreflightError("pinned evaluator parser changed intended row semantics")
    provenance = evaluator.__evaluator_provenance__
    return EvaluatorParserValidation(
        provenance["pin_sha256"],
        provenance["revision"],
        provenance["source_sha256"],
        provenance["source_git_blob_sha1"],
        provenance["source_size_bytes"],
    )


def prepare_validation_only_submission(
    selection: SubmissionSelection,
    assignment: EpcrAssignment,
    *,
    evaluator_source_path: str | Path | None = None,
    expected_preflight: ValidationOnlyPreflight | None = None,
) -> PreparedValidationOnlySubmission:
    """Render exact draft bytes while retaining all unavailable readiness gates."""

    try:
        selection = _typed(selection, SubmissionSelection)
        assignment = _typed(assignment, EpcrAssignment)
        _validate_selected_trace_bindings(selection)
        selection_sha256 = hashlib.sha256(selection.as_bytes()).hexdigest()
        assignment_sha256 = hashlib.sha256(assignment.as_bytes()).hexdigest()
        if assignment.selection_sha256 != selection_sha256:
            raise ValueError
        if len(selection.selected) != len(assignment.values):
            raise ValueError
        rows = []
        for group, value in zip(selection.selected, assignment.values, strict=True):
            if (
                value.hypothesis_id != group.primary_hypothesis_id
                or value.source_rank != group.source_rank
                or value.source_score != group.source_score
            ):
                raise ValueError
            rows.append(Track1CsvRow(
                group.primary_hypothesis_id,
                group.identity[0],
                None if len(group.identity) == 1 else group.identity[1],
                value.epcr,
                "Incomplete recessive partial member; partner unresolved" if group.partial_member else "",
            ))
        rows = tuple(rows)
        csv_bytes = None if not rows else _render_csv(rows)

        selected_read_evidence = any(
            source.authority == "authorized_read_evidence"
            for group in selection.selected
            for allele in group.alleles
            for source in allele.source_traces
        )
        blockers = list(_BASE_BLOCKERS)
        warnings = ["private-truth-unavailable-not-evaluated", "validation-only-not-submission-ready"]
        if selection.reference.identity_id != PUBLIC_REFERENCE_IDENTITY:
            blockers.append("production-public-reference-identity-unavailable")
        if selected_read_evidence:
            blockers.append("authorized-batch-consumer-unavailable")
        if not rows:
            blockers.append("zero-eligible-hypotheses")

        evaluator = None
        if evaluator_source_path is None:
            blockers.append("pinned-evaluator-source-unavailable")
        elif csv_bytes is None:
            warnings.append("pinned-evaluator-parser-not-applicable-no-csv")
        else:
            evaluator = _validate_evaluator_parser(csv_bytes, rows, evaluator_source_path)
        report = ValidationOnlyPreflight(
            selection_sha256,
            assignment_sha256,
            rows,
            None if csv_bytes is None else hashlib.sha256(csv_bytes).hexdigest(),
            evaluator,
            tuple(sorted(set(blockers))),
            tuple(sorted(warnings)),
            selected_read_evidence,
            "selected-read-evidence-present-unvalidated" if selected_read_evidence else "run-status-unavailable",
            "task-11.6-validation-unavailable" if selected_read_evidence else "run-snapshot-unavailable",
            assignment.confidence_label,
        )
        if expected_preflight is not None and report.as_bytes() != _typed(
            expected_preflight, ValidationOnlyPreflight,
        ).as_bytes():
            raise ValueError
        return PreparedValidationOnlySubmission(csv_bytes, report)
    except _ERRORS:
        raise Track1PreflightError("validation-only Track 1 preflight failed") from None


def _validate_bounded_selection(
    selection: ScientificSelectionProduct,
    scientific: ScientificArtifact,
    order: ScientificOrderProduct,
    *,
    scientific_directory: str | Path,
    order_directory: str | Path,
    selection_directory: str | Path,
    scientific_reader: ScientificArtifactReader | None = None,
    eligibility_root: str | Path | None = None,
    profile_root: str | Path | None = None,
) -> None:
    """Authenticate every bounded leaf and reconcile order/audit accounting."""

    if (
        type(selection) is not ScientificSelectionProduct
        or not scientific_artifact_complete(scientific)
        or type(order) is not ScientificOrderProduct
        or selection.scientific_manifest_sha256 != scientific.sha256
        or selection.order_manifest_sha256 != order.sha256
        or order.scientific_manifest_sha256 != scientific.sha256
        or selection.ranking_policy_sha256 != scientific.policy_sha256
        or order.policy_sha256 != scientific.policy_sha256
        or selection.ranked_hypothesis_count != scientific.hypothesis_count
    ):
        raise Track1PreflightError("bounded selection parent closure is invalid")
    manager = (ScientificArtifactReader(scientific, directory=scientific_directory,
            expected_parent_bindings=dict(scientific.parent_bindings),
            expected_reference_binding=selection.reference,
            scratch_directory=Path(scientific_directory).parent,
            eligibility_root=eligibility_root, profile_root=profile_root) if scientific_reader is None
            else nullcontext(scientific_reader))
    with manager as reader:
        if (type(reader) is not ScientificArtifactReader or reader.artifact != scientific
                or reader.directory != Path(scientific_directory)
                or reader.expected_reference_binding != selection.reference):
            raise Track1PreflightError("foreign bounded scientific reader")
        products = scientific_products(reader,
            expected_prefilter_manifest_sha256=dict(scientific.parent_bindings)["prefilter_manifest"])
        if products.order != order or products.selection != selection:
            raise Track1PreflightError("bounded selection does not replay from complete evidence")
    reference = selection.reference
    for group in selection.selected:
        for trace in group.alleles:
            if (
                trace.submission_reference_role != reference.role
                or trace.submission_reference_identity_id != reference.identity_id
                or trace.submission_reference_manifest_sha256
                != reference.manifest_sha256
                or trace.contig_aliases_sha256 != reference.contig_aliases_sha256
            ):
                raise Track1PreflightError(
                    "bounded selected allele trace differs from its reference"
                )


def prepare_bounded_validation_only_submission(
    selection: ScientificSelectionProduct,
    assignment: EpcrAssignment,
    scientific: ScientificArtifact,
    order: ScientificOrderProduct,
    *,
    scientific_directory: str | Path,
    order_directory: str | Path,
    selection_directory: str | Path,
    evaluator_source_path: str | Path | None = None,
    expected_preflight: ValidationOnlyPreflight | None = None,
    scientific_reader: ScientificArtifactReader | None = None,
    eligibility_root: str | Path | None = None,
    profile_root: str | Path | None = None,
) -> PreparedValidationOnlySubmission:
    """Render a draft only after complete sharded selection replay."""

    try:
        _validate_bounded_selection(
            selection,
            scientific,
            order,
            scientific_directory=scientific_directory,
            order_directory=order_directory,
            selection_directory=selection_directory,
            scientific_reader=scientific_reader,
            eligibility_root=eligibility_root,
            profile_root=profile_root,
        )
        assignment = _typed(assignment, EpcrAssignment)
        if assignment.selection_sha256 != selection.sha256:
            raise ValueError
        if len(selection.selected) != len(assignment.values):
            raise ValueError
        rows = []
        for group, value in zip(selection.selected, assignment.values, strict=True):
            if (
                value.hypothesis_id != group.primary_hypothesis_id
                or value.source_rank != group.source_rank
                or value.source_score != group.source_score
            ):
                raise ValueError
            rows.append(Track1CsvRow(
                group.primary_hypothesis_id,
                group.identity[0],
                None if len(group.identity) == 1 else group.identity[1],
                value.epcr,
                "Incomplete recessive partial member; partner unresolved" if group.partial_member else "",
            ))
        rows = tuple(rows)
        csv_bytes = None if not rows else _render_csv(rows)
        selected_read_evidence = any(
            source.authority == "authorized_read_evidence"
            for group in selection.selected
            for allele in group.alleles
            for source in allele.source_traces
        )
        blockers = list(_BASE_BLOCKERS)
        warnings = [
            "private-truth-unavailable-not-evaluated",
            "validation-only-not-submission-ready",
        ]
        if selection.reference.identity_id != PUBLIC_REFERENCE_IDENTITY:
            blockers.append("production-public-reference-identity-unavailable")
        if selected_read_evidence:
            blockers.append("authorized-batch-consumer-unavailable")
        if not rows:
            blockers.append("zero-eligible-hypotheses")
        evaluator = None
        if evaluator_source_path is None:
            blockers.append("pinned-evaluator-source-unavailable")
        elif csv_bytes is None:
            warnings.append("pinned-evaluator-parser-not-applicable-no-csv")
        else:
            evaluator = _validate_evaluator_parser(
                csv_bytes, rows, evaluator_source_path,
            )
        report = ValidationOnlyPreflight(
            selection.sha256,
            hashlib.sha256(assignment.as_bytes()).hexdigest(),
            rows,
            None if csv_bytes is None else hashlib.sha256(csv_bytes).hexdigest(),
            evaluator,
            tuple(sorted(set(blockers))),
            tuple(sorted(set(warnings))),
            selected_read_evidence,
            "selected-read-evidence-present-unvalidated"
            if selected_read_evidence else "run-status-unavailable",
            "task-11.6-validation-unavailable"
            if selected_read_evidence else "run-snapshot-unavailable",
            assignment.confidence_label,
        )
        if expected_preflight is not None and report.as_bytes() != _typed(
            expected_preflight, ValidationOnlyPreflight,
        ).as_bytes():
            raise ValueError
        return PreparedValidationOnlySubmission(csv_bytes, report)
    except _ERRORS:
        raise Track1PreflightError(
            "bounded validation-only Track 1 preflight failed"
        ) from None


def prepare_release_submission(
    selection: SubmissionSelection,
    assignment: EpcrAssignment,
    run_evidence: CompetitiveRunEvidence,
    *,
    release_manifest_bytes: bytes,
    expected_release_manifest_sha256: str,
    evaluator_source_path: str | Path | None,
    public_reference: PinnedPublicSource | None,
    expected_preflight: ReleasePreflight | None = None,
) -> PreparedReleaseSubmission:
    """Replay final bindings and emit competitive bytes only when every gate passes."""
    try:
        selection = _typed(selection, SubmissionSelection)
        assignment = _typed(assignment, EpcrAssignment)
        run_evidence = _typed(run_evidence, CompetitiveRunEvidence)
        _sha(expected_release_manifest_sha256)
        if (
            type(release_manifest_bytes) is not bytes
            or hashlib.sha256(release_manifest_bytes).hexdigest()
            != expected_release_manifest_sha256
            or run_evidence.release_manifest_sha256
            != expected_release_manifest_sha256
        ):
            raise ValueError
        release = load_base_release(
            release_manifest_bytes,
            expected_sha256=expected_release_manifest_sha256,
        )

        draft = prepare_validation_only_submission(
            selection,
            assignment,
            evaluator_source_path=evaluator_source_path,
        )
        selection_sha256 = hashlib.sha256(selection.as_bytes()).hexdigest()
        assignment_sha256 = hashlib.sha256(assignment.as_bytes()).hexdigest()
        selected_hypotheses = tuple(
            group.primary_hypothesis_id for group in selection.selected
        )
        selected_sources = {
            (
                trace.allele_id,
                source.source_allele_id,
                source.sample_name,
                source.source_input_sha256,
            )
            for group in selection.selected
            for trace in group.alleles
            for source in trace.source_traces
        }
        if (
            run_evidence.selection_sha256 != selection_sha256
            or run_evidence.epcr_assignment_sha256 != assignment_sha256
            or release.reference_identity_id
            != selection.reference.identity_id
            or any(
                trace.source_reference_identity_id != run_evidence.accepted_reference_identity_id
                or trace.source_reference_manifest_sha256 != release.reference_manifest_sha256
                for group in selection.selected for trace in group.alleles
            )
            or release.reference_manifest_sha256
            != run_evidence.accepted_reference_manifest_sha256
            or release.assessment_scope_sha256
            != run_evidence.assessment_scope_sha256
            or release.ranking_policy_sha256 != selection.ranking_policy_sha256
            or release.epcr_policy_sha256 != assignment.policy_sha256
            or assignment.selection_sha256 != selection_sha256
        ):
            raise ValueError

        blockers = []
        warnings = ["private-truth-unavailable-not-evaluated"]
        if run_evidence.accepted_reference_status != "accepted":
            blockers.append("accepted-reference-not-current")
        if run_evidence.concordance_status != "complete":
            blockers.append("complete-concordance-unavailable")
        if set(selected_hypotheses) - set(
            run_evidence.completely_enumerated_hypothesis_ids
        ):
            blockers.append("selected-required-enumeration-incomplete")
        assessment_sources = {
            (
                item.allele_id,
                item.source_allele_id,
                item.sample_name,
                item.source_input_sha256,
            )
            for item in run_evidence.technical_assessments
        }
        selected_assessments = tuple(
            item
            for item in run_evidence.technical_assessments
            if (
                item.allele_id,
                item.source_allele_id,
                item.sample_name,
                item.source_input_sha256,
            )
            in selected_sources
        )
        if (
            selected_sources - assessment_sources
            or any(item.status == "unmet" for item in selected_assessments)
        ):
            blockers.append("selected-allele-mandatory-technical-check-unmet")
        if run_evidence.additional_class_status == "not_assessed":
            warnings.append("additional-classes-not-assessed")
        if run_evidence.rescue_recommendations:
            warnings.append("rescue-recommended-not-automatically-invoked")

        selected_read_evidence = any(
            source.authority == "authorized_read_evidence"
            for group in selection.selected
            for allele in group.alleles
            for source in allele.source_traces
        )
        if selected_read_evidence:
            if (
                run_evidence.rescue_execution_status != "performed"
                or run_evidence.rescue_authorization_status != "authorized"
                or run_evidence.rescue_authorization_sha256 is None
                or run_evidence.rescue_plan_status != "complete"
                or any(item.status != "satisfied" for item in selected_assessments)
            ):
                blockers.append("selected-read-evidence-authorization-or-plan-incomplete")
        elif (
            run_evidence.rescue_execution_status != "not_performed"
            or run_evidence.rescue_authorization_status != "not_applicable"
            or run_evidence.rescue_authorization_sha256 is not None
            or run_evidence.rescue_plan_status != "not_applicable"
        ):
            blockers.append("unrescued-run-rescue-binding-invalid")

        if not selection.selected:
            blockers.append("zero-eligible-hypotheses")
        if draft.preflight.evaluator is None:
            blockers.append("pinned-evaluator-source-unavailable")
        elif draft.preflight.evaluator.pin_sha256 != release.evaluator_pin_sha256:
            raise ValueError
        if public_reference is None:
            blockers.append("public-native-reference-replay-unavailable")
        else:
            try:
                verify_selected_submission_reference(
                    selection,
                    reference=public_reference,
                    expected_public_identity_id=release.reference_identity_id,
                )
            except (ValueError, TypeError, RuntimeError, OSError):
                blockers.append("public-native-reference-replay-failed")

        draft_sha256 = (
            None
            if draft.csv_bytes is None
            else hashlib.sha256(draft.csv_bytes).hexdigest()
        )
        blockers = tuple(sorted(set(blockers)))
        preflight = ReleasePreflight(
            hashlib.sha256(run_evidence.as_bytes()).hexdigest(),
            expected_release_manifest_sha256,
            selection_sha256,
            assignment_sha256,
            draft.preflight.rows,
            draft_sha256,
            draft_sha256 if not blockers else None,
            draft.preflight.evaluator,
            blockers,
            tuple(sorted(set(warnings))),
            selected_read_evidence,
        )
        if expected_preflight is not None and preflight.as_bytes() != _typed(
            expected_preflight, ReleasePreflight
        ).as_bytes():
            raise ValueError
        return PreparedReleaseSubmission(
            draft.csv_bytes if preflight.ready else None,
            preflight,
        )
    except _ERRORS:
        raise Track1PreflightError("release Track 1 preflight failed") from None


def prepare_bounded_release_submission(
    selection: ScientificSelectionProduct,
    assignment: EpcrAssignment,
    run_evidence: CompetitiveRunEvidence,
    scientific: ScientificArtifact,
    order: ScientificOrderProduct,
    *,
    scientific_directory: str | Path,
    order_directory: str | Path,
    selection_directory: str | Path,
    release_manifest_bytes: bytes,
    expected_release_manifest_sha256: str,
    evaluator_source_path: str | Path | None,
    public_reference: PinnedPublicSource | None,
    expected_preflight: ReleasePreflight | None = None,
    read_intake_manifest=None,
    eligibility_root: str | Path | None = None,
    profile_root: str | Path | None = None,
    scientific_reader: ScientificArtifactReader | None = None,
    public_qualification_authority=None,
) -> PreparedReleaseSubmission:
    """Validate primary output under an accepted release or exact public qualification.

    Public qualification exercises the same renderer and scientific checks. It
    does not accept the release or authorize controlled input or operator use.
    """

    try:
        run_evidence = _typed(run_evidence, CompetitiveRunEvidence)
        _sha(expected_release_manifest_sha256)
        if (
            type(release_manifest_bytes) is not bytes
            or hashlib.sha256(release_manifest_bytes).hexdigest()
            != expected_release_manifest_sha256
            or run_evidence.release_manifest_sha256
            != expected_release_manifest_sha256
        ):
            raise ValueError
        release = load_base_release(
            release_manifest_bytes,
            expected_sha256=expected_release_manifest_sha256,
        )
        draft = prepare_bounded_validation_only_submission(
            selection,
            assignment,
            scientific,
            order,
            scientific_directory=scientific_directory,
            order_directory=order_directory,
            selection_directory=selection_directory,
            evaluator_source_path=evaluator_source_path,
            eligibility_root=eligibility_root,
            profile_root=profile_root,
            scientific_reader=scientific_reader,
        )
        assignment_sha256 = hashlib.sha256(assignment.as_bytes()).hexdigest()
        selected_hypotheses = tuple(
            group.primary_hypothesis_id for group in selection.selected
        )
        selected_sources = {
            (
                trace.allele_id,
                source.source_allele_id,
                source.sample_name,
                source.source_input_sha256,
            )
            for group in selection.selected
            for trace in group.alleles
            for source in trace.source_traces
        }
        if (
            run_evidence.selection_sha256 != selection.sha256
            or run_evidence.epcr_assignment_sha256 != assignment_sha256
            or release.reference_identity_id != selection.reference.identity_id
            or any(
                trace.source_reference_identity_id
                != run_evidence.accepted_reference_identity_id
                or trace.source_reference_manifest_sha256
                != release.reference_manifest_sha256
                for group in selection.selected for trace in group.alleles
            )
            or release.reference_manifest_sha256
            != run_evidence.accepted_reference_manifest_sha256
            or release.assessment_scope_sha256
            != run_evidence.assessment_scope_sha256
            or release.ranking_policy_sha256
            != selection.ranking_policy_sha256
            or release.epcr_policy_sha256 != assignment.policy_sha256
            or assignment.selection_sha256 != selection.sha256
        ):
            raise ValueError

        blockers = []
        warnings = ["private-truth-unavailable-not-evaluated"]
        if public_qualification_authority is not None:
            from genosapiens.wgs_public_qualification import (
                PublicQualificationAuthority, require_public_qualification_release,
            )
            if type(public_qualification_authority) is not PublicQualificationAuthority:
                raise ValueError
            require_public_qualification_release(public_qualification_authority, release)
            if run_evidence.rescue_execution_status != "not_performed" or read_intake_manifest is not None:
                raise ValueError
            warnings.append("public-qualification-only-no-operator-authorization")
        if run_evidence.accepted_reference_status != "accepted":
            blockers.append("accepted-reference-not-current")
        if run_evidence.concordance_status != "complete":
            blockers.append("complete-concordance-unavailable")
        if set(selected_hypotheses) - set(
            run_evidence.completely_enumerated_hypothesis_ids
        ):
            blockers.append("selected-required-enumeration-incomplete")
        assessment_sources = {
            (
                item.allele_id,
                item.source_allele_id,
                item.sample_name,
                item.source_input_sha256,
            )
            for item in run_evidence.technical_assessments
        }
        selected_assessments = tuple(
            item
            for item in run_evidence.technical_assessments
            if (
                item.allele_id,
                item.source_allele_id,
                item.sample_name,
                item.source_input_sha256,
            ) in selected_sources
        )
        if (
            selected_sources - assessment_sources
            or any(item.status == "unmet" for item in selected_assessments)
        ):
            blockers.append("selected-allele-mandatory-technical-check-unmet")
        if run_evidence.additional_class_status == "not_assessed":
            warnings.append("additional-classes-not-assessed")
        if run_evidence.rescue_recommendations:
            warnings.append("rescue-recommended-not-automatically-invoked")

        selected_read_evidence = any(
            source.authority == "authorized_read_evidence"
            for group in selection.selected
            for allele in group.alleles
            for source in allele.source_traces
        )
        if selected_read_evidence or run_evidence.rescue_execution_status == "performed":
            if (
                run_evidence.rescue_execution_status != "performed"
                or run_evidence.rescue_authorization_status != "authorized"
                or run_evidence.rescue_authorization_sha256 is None
                or run_evidence.rescue_plan_status != "complete"
                or read_intake_manifest is None
            ):
                blockers.append(
                    "selected-read-evidence-authorization-or-plan-incomplete"
                )
            if read_intake_manifest is not None:
                from genosapiens.wgs_read_rescue import validate_read_intake_run
                validate_read_intake_run(read_intake_manifest, run_evidence=run_evidence, release=release)
                read_sources = {(trace.allele_id, source.source_allele_id, source.sample_name, source.source_input_sha256)
                    for group in selection.selected for trace in group.alleles for source in trace.source_traces
                    if source.authority == "authorized_read_evidence"}
                if any(item.status != "satisfied" for item in selected_assessments
                       if (item.allele_id, item.source_allele_id, item.sample_name, item.source_input_sha256) in read_sources):
                    blockers.append("selected-read-evidence-authorization-or-plan-incomplete")
        elif (
            run_evidence.rescue_execution_status != "not_performed"
            or run_evidence.rescue_authorization_status != "not_applicable"
            or run_evidence.rescue_authorization_sha256 is not None
            or run_evidence.rescue_plan_status != "not_applicable"
            or read_intake_manifest is not None and read_intake_manifest.deltas
        ):
            blockers.append("unrescued-run-rescue-binding-invalid")
        elif read_intake_manifest is not None:
            from genosapiens.wgs_read_rescue import validate_read_intake_run
            validate_read_intake_run(read_intake_manifest, run_evidence=run_evidence, release=release)

        if not selection.selected:
            blockers.append("zero-eligible-hypotheses")
        if draft.preflight.evaluator is None:
            blockers.append("pinned-evaluator-source-unavailable")
        elif draft.preflight.evaluator.pin_sha256 != release.evaluator_pin_sha256:
            raise ValueError
        if public_reference is None:
            blockers.append("public-native-reference-replay-unavailable")
        else:
            try:
                verify_sharded_selected_submission_reference(
                    selection,
                    reference=public_reference,
                    expected_public_identity_id=release.reference_identity_id,
                )
            except (ValueError, TypeError, RuntimeError, OSError):
                blockers.append("public-native-reference-replay-failed")

        draft_sha256 = (
            None if draft.csv_bytes is None
            else hashlib.sha256(draft.csv_bytes).hexdigest()
        )
        blockers = tuple(sorted(set(blockers)))
        preflight = ReleasePreflight(
            hashlib.sha256(run_evidence.as_bytes()).hexdigest(),
            expected_release_manifest_sha256,
            selection.sha256,
            assignment_sha256,
            draft.preflight.rows,
            draft_sha256,
            draft_sha256 if not blockers else None,
            draft.preflight.evaluator,
            blockers,
            tuple(sorted(set(warnings))),
            selected_read_evidence,
        )
        if expected_preflight is not None and preflight.as_bytes() != _typed(
            expected_preflight, ReleasePreflight,
        ).as_bytes():
            raise ValueError
        return PreparedReleaseSubmission(
            draft.csv_bytes if preflight.ready else None,
            preflight,
        )
    except _ERRORS:
        raise Track1PreflightError(
            "bounded release Track 1 preflight failed"
        ) from None
