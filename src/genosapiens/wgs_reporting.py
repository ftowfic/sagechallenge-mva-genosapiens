"""Authenticated bounded diagnostics and streamed Track 1 WGS audit reports."""

from __future__ import annotations

from contextlib import ExitStack, contextmanager
from dataclasses import asdict, dataclass, field
import hashlib
import json
import os
from pathlib import Path
import re
import sqlite3
import stat
from typing import Literal

from genosapiens.bounded_wgs_config import BoundedWgsExecutionConfig
from genosapiens.evidence_model import Completeness, _exact, _typed
from genosapiens.evidence_serialization import _canonical, _loads
from genosapiens.track1_preflight import SelectedTechnicalAssessment
from genosapiens.variant_inventory import IncrementalVariantClassSummary
from genosapiens.vcf_qc import VcfQcIncrementalSummary
from genosapiens.wgs_filter import verify_working_set
from genosapiens.wgs_ledger import (
    FilterManifest,
    LedgerManifest,
    WorkingSetManifest,
    canonical,
    iter_ledger,
    read_artifact,
)
from genosapiens.wgs_preparation import (
    SourceAccountingClosure,
    WgsNormalizationSummary,
)
from genosapiens.wgs_source_routing import SourceRoutingContext
from genosapiens.wgs_evidence_store import EvidenceReference
from genosapiens.wgs_scientific_artifact import ScientificArtifact, ScientificArtifactReader
from genosapiens.wgs_scientific_products import (
    ScientificOrderProduct, ScientificPrimaryProduct,
    ScientificSelectionProduct, ScientificProducts, validate_scientific_products, scientific_report_roots,
)


DIAGNOSTIC_SCHEMA = "genosapiens.wgs-bounded-diagnostics"
AUDIT_SCHEMA = "genosapiens.wgs-bounded-audit-manifest"
VERSION = 6
AUDIT_VERSION = 6
_DIGEST = re.compile(r"^[0-9a-f]{64}$")
_ERRORS = (
    ValueError, TypeError, RuntimeError, AttributeError, KeyError, IndexError,
    OSError, sqlite3.Error, UnicodeError, OverflowError, RecursionError,
)


class WgsReportingError(RuntimeError):
    """Static bounded-report failure without controlled row content."""


def _bytes(value: object) -> bytes:
    return (_canonical(value) + "\n").encode("ascii")


def _digest(value: object) -> bool:
    return type(value) is str and _DIGEST.fullmatch(value) is not None


def _safe_relative(value: str) -> bool:
    path = Path(value)
    return (
        type(value) is str
        and value == path.as_posix()
        and not path.is_absolute()
        and path.parts
        and all(part not in ("", ".", "..") for part in path.parts)
    )


def _read_exact_private(path: Path, expected: bytes) -> bytes:
    descriptor = -1
    try:
        metadata = path.lstat()
        if (
            path.is_symlink()
            or not stat.S_ISREG(metadata.st_mode)
            or metadata.st_uid != os.geteuid()
            or stat.S_IMODE(metadata.st_mode) != 0o600
            or metadata.st_size != len(expected)
        ):
            raise OSError
        descriptor = os.open(path, os.O_RDONLY | os.O_NOFOLLOW)
        with os.fdopen(descriptor, "rb", closefd=True) as stream:
            descriptor = -1
            content = stream.read(len(expected) + 1)
        if content != expected:
            raise OSError
        return content
    except OSError:
        raise WgsReportingError("bounded WGS report replay failed") from None
    finally:
        if descriptor >= 0:
            os.close(descriptor)


@dataclass(frozen=True, slots=True, repr=False)
class ReportShard:
    role: Literal["audit_jsonl", "audit_markdown"]
    relative_path: str
    sha256: str
    size_bytes: int
    row_count: int
    first_key: str
    last_key: str

    def __post_init__(self) -> None:
        if (
            not _safe_relative(self.relative_path)
            or not _digest(self.sha256)
            or type(self.size_bytes) is not int
            or self.size_bytes < 1
            or type(self.row_count) is not int
            or self.row_count < 1
            or not self.first_key
            or not self.last_key
        ):
            raise WgsReportingError("bounded report shard is invalid")


@dataclass(frozen=True, slots=True, repr=False)
class BoundedDiagnosticReport:
    source_manifest_sha256: str
    source_accounting_sha256: str
    qc_summary_sha256: str
    normalization_summary_sha256: str
    class_summary_sha256: str
    filter_manifest_sha256: str
    working_manifest_sha256: str
    scientific_manifest_sha256: str
    order_manifest_sha256: str
    primary_manifest_sha256: str
    selection_sha256: str
    source_record_count: int
    source_physical_alternate_count: int
    retained_occurrence_count: int
    excluded_common_occurrence_count: int
    ranked_hypothesis_count: int
    selected_group_count: int
    original_class_evidence_status: Literal["complete"]
    original_technical_evidence_status: Literal["complete"]
    credible_complete_status: Literal["present", "unknown", "not_evaluated"]
    credible_complete_count: int
    selected_obligation_status: Literal["satisfied", "unmet", "not_applicable"]
    selected_technical_assessments: tuple[SelectedTechnicalAssessment, ...]
    blockers: tuple[str, ...]
    warnings: tuple[str, ...]
    assessment_domain_sha256: str
    assessment_input_count: int | None
    assessment_input_dispositions: tuple[tuple[str, int], ...] | None
    scope_policy_sha256: str = field(kw_only=True)
    scope_manifest_sha256: str = field(kw_only=True)
    scope_counts: tuple[tuple[str, int], ...] = field(kw_only=True)
    in_scope_normalized_occurrence_count: int = field(kw_only=True)
    retained_identity_count: int = field(kw_only=True)
    excluded_scope_occurrence_count: int = field(kw_only=True)
    unresolved_scope_occurrence_count: int = field(kw_only=True)
    bypassed_occurrence_count: int = field(kw_only=True)
    population_request_keys: int = field(kw_only=True)
    clinical_query_keys: int = field(kw_only=True)
    optional_analysis_status: Literal["not_run"] = "not_run"
    schema: Literal["genosapiens.wgs-bounded-diagnostics"] = DIAGNOSTIC_SCHEMA
    version: Literal[6] = VERSION
    classification: Literal["controlled"] = "controlled"
    llm_policy: Literal["no-llm"] = "no-llm"

    def __post_init__(self) -> None:
        for value in (
            self.source_manifest_sha256,
            self.source_accounting_sha256,
            self.qc_summary_sha256,
            self.normalization_summary_sha256,
            self.class_summary_sha256,
            self.filter_manifest_sha256,
            self.working_manifest_sha256,
            self.scientific_manifest_sha256,
            self.order_manifest_sha256,
            self.primary_manifest_sha256,
            self.selection_sha256,
            self.assessment_domain_sha256,
            self.scope_policy_sha256, self.scope_manifest_sha256,
        ):
            if not _digest(value):
                raise WgsReportingError("bounded diagnostic identity is invalid")
        from genosapiens.wgs_ledger import SCOPE_DISPOSITIONS
        scope = self.scope_counts
        if (type(scope) is not tuple or scope != tuple(sorted(scope))
                or len(scope) != 4 or {key for key, _ in scope} != set(SCOPE_DISPOSITIONS)
                or any(type(n) is not int or n < 0 for _, n in scope)
                or sum(n for _, n in scope) != self.source_physical_alternate_count
                or self.in_scope_normalized_occurrence_count > dict(scope)['in_scope']
                or self.population_request_keys > self.in_scope_normalized_occurrence_count
                or self.clinical_query_keys > self.in_scope_normalized_occurrence_count
                or self.retained_identity_count > self.retained_occurrence_count
                or self.retained_occurrence_count + self.excluded_common_occurrence_count
                    + self.excluded_scope_occurrence_count + self.unresolved_scope_occurrence_count
                    + self.bypassed_occurrence_count != self.source_physical_alternate_count):
            raise WgsReportingError("bounded diagnostic scope accounting is invalid")
        from genosapiens.wgs_analysis_eligibility import DISPOSITIONS
        counts = self.assessment_input_dispositions
        if ((counts is None) != (self.assessment_input_count is None)
                or self.optional_analysis_status != "not_run"
                or counts is not None and (counts != tuple(sorted(counts))
                    or len(counts) != 3 or {key for key,_ in counts} != DISPOSITIONS
                    or any(type(n) is not int or n < 0 for _,n in counts)
                    or type(self.assessment_input_count) is not int
                    or sum(n for _,n in counts) != self.assessment_input_count)):
            raise WgsReportingError("bounded diagnostic assessment domain is invalid")
        counts = (
            self.source_record_count,
            self.source_physical_alternate_count,
            self.retained_occurrence_count,
            self.excluded_common_occurrence_count,
            self.excluded_scope_occurrence_count, self.unresolved_scope_occurrence_count,
            self.bypassed_occurrence_count, self.retained_identity_count,
            self.in_scope_normalized_occurrence_count, self.population_request_keys, self.clinical_query_keys,
            self.ranked_hypothesis_count,
            self.selected_group_count,
            self.credible_complete_count,
        )
        if (
            self.schema != DIAGNOSTIC_SCHEMA or type(self.version) is not int or self.version != VERSION
            or self.classification != "controlled" or self.llm_policy != "no-llm"
            or any(type(value) is not int or value < 0 for value in counts)
            or self.selected_group_count > 10
            or self.blockers != tuple(sorted(set(self.blockers)))
            or self.warnings != tuple(sorted(set(self.warnings)))
            or self.selected_technical_assessments
            != tuple(sorted(
                set(self.selected_technical_assessments),
                key=lambda value: (
                    value.allele_id,
                    value.source_allele_id,
                    value.sample_name,
                    value.reason,
                    value.scope,
                ),
            ))
            or self.credible_complete_status == "present"
            and self.credible_complete_count == 0
            or self.credible_complete_status != "present"
            and self.credible_complete_count != 0
            or self.selected_obligation_status == "not_applicable"
            and any(
                value.status != "not_applicable"
                for value in self.selected_technical_assessments
            )
            or self.selected_obligation_status == "satisfied"
            and (
                not any(value.status == "satisfied" for value in self.selected_technical_assessments)
                or any(value.status == "unmet" for value in self.selected_technical_assessments)
            )
        ):
            raise WgsReportingError("bounded diagnostic report is invalid")

    def as_bytes(self) -> bytes:
        return _bytes({
            "classification": "controlled",
            "llm_policy": "no-llm",
            "report": self,
            "schema": self.schema,
            "version": self.version,
        })

    @property
    def sha256(self) -> str:
        return hashlib.sha256(self.as_bytes()).hexdigest()


@dataclass(frozen=True, slots=True, repr=False)
class BoundedAuditManifest:
    diagnostic_sha256: str
    parent_bindings: tuple[tuple[str, str], ...]
    json_shards: tuple[ReportShard, ...]
    markdown_shards: tuple[ReportShard, ...]
    source_rows: int
    representation_rows: int
    hypothesis_rows: int
    audit_rows: int
    maximum_rows_per_shard: int
    maximum_bytes_per_shard: int
    maximum_total_json_bytes: int
    maximum_total_markdown_bytes: int
    json_stream_sha256: str
    markdown_stream_sha256: str
    no_candidate: bool
    scientific_roots: tuple[tuple[str, EvidenceReference], ...]
    detailed_hypothesis_rows: int
    assessment_domain_sha256: str
    source_accounting_rows: int = 1
    optional_analysis_status: Literal["not_run"] = "not_run"
    schema: Literal["genosapiens.wgs-bounded-audit-manifest"] = AUDIT_SCHEMA
    version: Literal[4, 6] = AUDIT_VERSION
    classification: Literal["controlled"] = "controlled"
    llm_policy: Literal["no-llm"] = "no-llm"

    def __post_init__(self) -> None:
        if (not _digest(self.diagnostic_sha256) or not _digest(self.assessment_domain_sha256)
                or self.optional_analysis_status != "not_run"):
            raise WgsReportingError("bounded audit diagnostic binding is invalid")
        for value in (self.json_stream_sha256, self.markdown_stream_sha256):
            if not _digest(value):
                raise WgsReportingError("bounded audit stream identity is invalid")
        if (
            self.schema != AUDIT_SCHEMA or type(self.version) is not int or self.version not in (4, AUDIT_VERSION)
            or type(self.source_accounting_rows) is not int
            or self.source_accounting_rows != (self.source_rows if self.version == 4 else 1)
            or self.classification != "controlled" or self.llm_policy != "no-llm"
            or self.parent_bindings != tuple(sorted(self.parent_bindings))
            or not self.parent_bindings
            or len({name for name, _value in self.parent_bindings})
            != len(self.parent_bindings)
            or any(not name or not _digest(value)
                   for name, value in self.parent_bindings)
            or any(value.role != "audit_jsonl" for value in self.json_shards)
            or any(value.role != "audit_markdown"
                   for value in self.markdown_shards)
            or any(type(value) is not int or value < 0 for value in (
                self.source_rows,
                self.representation_rows,
                self.hypothesis_rows,
                self.audit_rows,
                self.detailed_hypothesis_rows,
            ))
            or any(type(value) is not int or value <= 0 for value in (
                self.maximum_rows_per_shard,
                self.maximum_bytes_per_shard,
                self.maximum_total_json_bytes,
                self.maximum_total_markdown_bytes,
            ))
            or self.audit_rows != sum(value.row_count for value in self.json_shards)
            or self.audit_rows != sum(value.row_count for value in self.markdown_shards)
            or self.no_candidate != (self.hypothesis_rows == 0)
            or self.detailed_hypothesis_rows > min(10, self.hypothesis_rows)
            or self.scientific_roots != tuple(sorted(self.scientific_roots))
            or tuple(name for name, _ in self.scientific_roots) != (
                "case", "primary", "primary_order",
                "reference", "resources", "selected", "selection")
            or any(type(root) is not EvidenceReference for _, root in self.scientific_roots)
            or len({root.manifest_sha256 for _, root in self.scientific_roots}) != 1
            or dict(self.parent_bindings).get("scientific_evidence")
                != self.scientific_roots[0][1].manifest_sha256
        ):
            raise WgsReportingError("bounded audit manifest is invalid")

    def as_bytes(self) -> bytes:
        manifest = asdict(self)
        if self.version == 4:
            manifest.pop('source_accounting_rows')
        return _bytes({
            "classification": "controlled",
            "llm_policy": "no-llm",
            "manifest": manifest,
            "schema": self.schema,
            "version": self.version,
        })

    @property
    def sha256(self) -> str:
        return hashlib.sha256(self.as_bytes()).hexdigest()

    def markdown_index(self) -> bytes:
        lines = [
            "# Track 1 bounded WGS audit",
            "",
            "Classification: `controlled`. LLM policy: `no-llm`.",
            "",
            f"Audit manifest SHA-256: `{self.sha256}`",
            f"Diagnostic SHA-256: `{self.diagnostic_sha256}`",
            "",
            f"Completely evaluated screened hypotheses: {self.hypothesis_rows}. Optional analyses: {self.optional_analysis_status}.",
            f"Selected hypothesis details: {self.detailed_hypothesis_rows}.",
            "",
            "Ranks apply only to evaluated contexts. The scientific evidence row states the scope and incomplete work.",
            "The normalized roots preserve evaluated scientific evidence, orders and selection dispositions:",
            "",
        ]
        lines.extend(f"- {name}: `{root.name}` — `{root.sha256}`" for name, root in self.scientific_roots)
        lines.extend([
            "",
            "These bounded shards store source accounting, member representations and selected details:",
            "",
        ])
        lines.extend(
            f"- `{value.relative_path}` — `{value.sha256}`"
            for value in (*self.json_shards, *self.markdown_shards)
        )
        return ("\n".join(lines) + "\n").encode("utf-8")


@dataclass(frozen=True, slots=True, repr=False)
class PreparedBoundedReports:
    diagnostic: BoundedDiagnosticReport
    audit: BoundedAuditManifest
    audit_json_index: bytes
    audit_markdown_index: bytes
    output_directory: Path

    def __post_init__(self) -> None:
        if (
            self.audit.diagnostic_sha256 != self.diagnostic.sha256
            or self.audit.assessment_domain_sha256 != self.diagnostic.assessment_domain_sha256
            or self.audit_json_index != self.audit.as_bytes()
            or self.audit_markdown_index != self.audit.markdown_index()
            or not self.output_directory.is_absolute()
        ):
            raise WgsReportingError("prepared bounded reports are invalid")

    def __repr__(self) -> str:
        return "<controlled bounded WGS reports; no-llm>"


def diagnostic_from_json(
    content: bytes,
    *,
    expected_sha256: str,
) -> BoundedDiagnosticReport:
    """Strict canonical bounded-diagnostic readback."""

    try:
        if (
            type(content) is not bytes
            or not _digest(expected_sha256)
            or hashlib.sha256(content).hexdigest() != expected_sha256
        ):
            raise ValueError
        data = _loads(content)
        if (
            type(data) is not dict
            or set(data) != {
                "classification", "llm_policy", "report", "schema", "version",
            }
            or data["classification"] != "controlled"
            or data["llm_policy"] != "no-llm"
            or data["schema"] != DIAGNOSTIC_SCHEMA
            or data["version"] != VERSION
        ):
            raise ValueError
        value = _typed(data["report"], BoundedDiagnosticReport, decode=True)
        if value.as_bytes() != content:
            raise ValueError
        return value
    except _ERRORS:
        raise WgsReportingError("invalid canonical bounded diagnostics") from None


def audit_manifest_from_json(
    content: bytes,
    *,
    expected_sha256: str,
) -> BoundedAuditManifest:
    """Strict canonical bounded-audit manifest readback."""

    try:
        if (
            type(content) is not bytes
            or not _digest(expected_sha256)
            or hashlib.sha256(content).hexdigest() != expected_sha256
        ):
            raise ValueError
        data = _loads(content)
        if (
            type(data) is not dict
            or set(data) != {
                "classification", "llm_policy", "manifest", "schema", "version",
            }
            or data["classification"] != "controlled"
            or data["llm_policy"] != "no-llm"
            or data["schema"] != AUDIT_SCHEMA
            or data["version"] not in (4, AUDIT_VERSION)
        ):
            raise ValueError
        manifest = dict(data['manifest'])
        if data['version'] == 4:
            if 'source_accounting_rows' in manifest:
                raise ValueError
            manifest['source_accounting_rows'] = manifest['source_rows']
        if manifest['version'] != data['version']:
            raise ValueError
        value = _typed(manifest, BoundedAuditManifest, decode=True)
        if value.as_bytes() != content:
            raise ValueError
        return value
    except _ERRORS:
        raise WgsReportingError("invalid canonical bounded audit manifest") from None


class _ReportWriter:
    def __init__(
        self,
        root: Path,
        role: Literal["audit_jsonl", "audit_markdown"],
        suffix: str,
        *,
        maximum_rows: int,
        maximum_bytes: int,
        maximum_total_bytes: int,
    ) -> None:
        self.root = root
        self.role = role
        self.suffix = suffix
        self.maximum_rows = maximum_rows
        self.maximum_bytes = maximum_bytes
        self.maximum_total_bytes = maximum_total_bytes
        self.shards: list[ReportShard] = []
        self.handle = None
        self.path = None
        self.digest = None
        self.size = self.rows = self.total = 0
        self.first = self.last = None
        self.stream = hashlib.sha256()

    def _open(self) -> None:
        ordinal = len(self.shards) + 1
        name = f"{self.role}-{ordinal:06d}.{self.suffix}.partial"
        self.path = self.root / name
        descriptor = os.open(
            self.path,
            os.O_WRONLY | os.O_CREAT | os.O_EXCL | os.O_NOFOLLOW,
            0o600,
        )
        self.handle = os.fdopen(descriptor, "wb")
        self.digest = hashlib.sha256()
        self.size = self.rows = 0
        self.first = self.last = None

    def write(self, content: bytes, *, key: str) -> None:
        if type(content) is not bytes or not content or not key:
            raise WgsReportingError("bounded report row is invalid")
        if len(content) > self.maximum_bytes:
            raise WgsReportingError("bounded report row exceeds its shard limit")
        if self.total + len(content) > self.maximum_total_bytes:
            raise WgsReportingError("bounded report exceeds its total file limit")
        if self.handle is None:
            self._open()
        if self.rows and (
            self.rows >= self.maximum_rows
            or self.size + len(content) > self.maximum_bytes
        ):
            self._close_shard()
            self._open()
        assert self.handle is not None and self.digest is not None
        self.handle.write(content)
        self.digest.update(content)
        self.stream.update(content)
        self.size += len(content)
        self.total += len(content)
        self.rows += 1
        self.first = key if self.first is None else self.first
        self.last = key

    def _close_shard(self) -> None:
        if self.handle is None:
            return
        assert self.path is not None and self.digest is not None
        self.handle.flush()
        os.fsync(self.handle.fileno())
        self.handle.close()
        final = self.path.with_suffix("")
        os.replace(self.path, final)
        self.shards.append(ReportShard(
            self.role,
            final.name,
            self.digest.hexdigest(),
            self.size,
            self.rows,
            self.first,
            self.last,
        ))
        self.handle = self.path = self.digest = None

    def close(self) -> tuple[ReportShard, ...]:
        self._close_shard()
        return tuple(self.shards)

    def abort(self) -> None:
        if self.handle is not None:
            self.handle.close()
        if self.path is not None:
            try:
                self.path.unlink()
            except FileNotFoundError:
                pass


def _report_row(kind: str, key: str, data: object) -> bytes:
    return _bytes({"data": data, "key": key, "kind": kind})


def _markdown_row(kind: str, key: str, summary: str) -> bytes:
    clean = summary.replace("\r", " ").replace("\n", " ").replace("|", "\\|")
    return f"- `{key}` | {kind} | {clean}\n".encode("utf-8")


def _source_technical_projection(row) -> dict[str, object]:
    from genosapiens.vcf_qc import VcfRecordQc, quality_acceptance_for
    payload = json.loads(row.payload_json)
    qc = payload["qc"]
    inventory = payload["inventory"]
    sample_index = payload["selected_sample_index"]
    record = _typed(qc, VcfRecordQc, decode=True)
    return {
        "record_number": row.record_number,
        "source_record_id": row.source_record_id,
        "source_allele_ids": row.source_allele_ids,
        "biological_alternate_count": row.biological_alternate_count,
        "variant_classes": qc["variant_classes"],
        "filter_value": record.filter_value,
        "filter_status": record.filter_status.value,
        "quality_acceptance": tuple(quality_acceptance_for(record, record.sample_calls[sample_index], number)
                                    for number in range(1, len(record.alternates) + 1)),
        "record_technical_flags": qc["technical_flags"],
        "alternate_technical_flags": qc["alternate_technical_flags"],
        "sample_technical_flags": qc["sample_technical_flags"][sample_index],
        "sample_alternate_technical_flags": qc[
            "sample_alternate_technical_flags"
        ][sample_index],
        "represented_classes": inventory["represented_classes"],
        "unsupported_classes": inventory["unsupported_classes"],
        "screen_calls": inventory["screen_calls"],
    }


def _verify_source_summaries(
    root: Path,
    manifest: LedgerManifest,
    closure: SourceAccountingClosure,
    qc: VcfQcIncrementalSummary,
    normalization: WgsNormalizationSummary,
    classes: IncrementalVariantClassSummary,
) -> None:
    reports = manifest.summary.reports
    if (
        closure.source_manifest_sha256 != manifest.sha256
        or closure.input_sha256 != manifest.source.input_sha256
        or len(reports) != 3
        or tuple(value.name for value in reports)
        != ("qc.json", "normalization.json", "classes.json")
        or tuple(value.sha256 for value in reports)
        != (
            closure.qc_summary_sha256,
            closure.normalization_summary_sha256,
            closure.class_summary_sha256,
        )
        or qc.input_sha256 != manifest.source.input_sha256
        or normalization.input_sha256 != manifest.source.input_sha256
        or classes.input_sha256 != manifest.source.input_sha256
        or read_artifact(root, reports[0]) != canonical(asdict(qc))
        or read_artifact(root, reports[1]) != canonical(asdict(normalization))
        or read_artifact(root, reports[2]) != canonical(asdict(classes))
        or qc.total_records != closure.total_records
        or qc.total_physical_alternate_slots
        != closure.total_physical_alternate_slots
        or classes.total_records != closure.total_records
    ):
        raise WgsReportingError("source summary closure is invalid")


def validate_profile_source_accounting(value, *, index_sha256=None, catalog_sha256=None):
    """Validate separately reported source denominators and their source bindings."""
    try:
        scope = json.loads(canonical(value))
        expected = {'schema','version','status','regulatory_index_manifest_sha256',
            'regulatory_manifest','regulatory_counts','regulatory_inventory','catalog_sha256','catalog_source_sha256',
            'catalog_records','catalog_genes','catalog_gene_disease_links','limitations'}
        if (type(scope) is not dict or set(scope) != expected
                or scope['schema'] != 'genosapiens.track1-profile-source-accounting'
                or type(scope['version']) is not int or scope['version'] != 1
                or scope['status'] != 'authenticated_source_scope'):
            raise ValueError
        for key in ('regulatory_index_manifest_sha256','catalog_sha256','catalog_source_sha256'):
            if not _digest(scope[key]):
                raise ValueError
        if (index_sha256 is not None and index_sha256 != scope['regulatory_index_manifest_sha256']
                or catalog_sha256 is not None and catalog_sha256 != scope['catalog_sha256']):
            raise ValueError
        counts, inventory = scope['regulatory_counts'], scope['regulatory_inventory']
        manifest = scope['regulatory_manifest']
        if (hashlib.sha256(canonical(manifest)).hexdigest() != scope['regulatory_index_manifest_sha256']
                or manifest.get('schema') != 'genosapiens.track1-regulatory-index'
                or manifest.get('version') != 2 or type(manifest['version']) is not int
                or manifest.get('status') != 'complete' or manifest.get('classification') != 'public'
                or manifest.get('counts') != counts or manifest.get('inventory') != inventory):
            raise ValueError
        if (type(counts) is not dict or set(counts) != {'ccre_rows','anchor_family_counts',
                're2g_rows','provenance_chunks','unique_legacy_elements','unique_source_memberships',
                'distinct_raw_targets','stable_memberships','unresolved_memberships','linked_memberships',
                'unjoined_memberships','unjoined_rows','provenance_memberships','rows_by_target_kind',
                'relations_by_kind','distinct_targets_by_kind'}
                or type(inventory) is not dict or set(inventory) != {'collection','physical_counts',
                    'stream_sha256','vocabulary_path','vocabulary_counts','target_resolution_complete'}
                or scope['limitations'] != ['source-wide counts are not candidate-level target counts',
                    'complete parsing does not establish complete target identity resolution',
                    'source links do not establish causal regulatory effects']):
            raise ValueError
        def numbers(values):
            if type(values) is not dict or any(type(v) is not int or v < 0 for v in values.values()):
                raise ValueError
        numbers({k:scope[k] for k in ('catalog_records','catalog_genes','catalog_gene_disease_links')})
        numbers({k:v for k,v in counts.items() if not isinstance(v,dict)})
        for key in ('anchor_family_counts','rows_by_target_kind','relations_by_kind','distinct_targets_by_kind'):
            numbers(counts[key])
        for key in ('rows_by_target_kind','relations_by_kind','distinct_targets_by_kind'):
            if set(counts[key]) != {'ensembl','unmapped','missing'}:
                raise ValueError
        if (not set(counts['anchor_family_counts']) <= {'EH38D','EH38F'}
                or type(inventory['stream_sha256']) is not dict
                or set(inventory['stream_sha256']) != {'ccre','re2g'}
                or any(not _digest(v) for v in inventory['stream_sha256'].values())
                or inventory['vocabulary_path'] != 'index.sqlite#vocabulary'
                or type(inventory['vocabulary_counts']) is not dict):
            raise ValueError
        for item in inventory['vocabulary_counts'].values():
            numbers(item)
            if set(item) != {'distinct_values','source_rows'}:
                raise ValueError
        if (sum(counts['anchor_family_counts'].values()) != counts['ccre_rows']
                or sum(counts['rows_by_target_kind'].values()) != counts['re2g_rows']
                or sum(counts['relations_by_kind'].values()) != counts['unique_source_memberships']
                or sum(counts['distinct_targets_by_kind'].values()) != counts['distinct_raw_targets']
                or counts['relations_by_kind']['ensembl'] != counts['stable_memberships']
                or counts['stable_memberships'] + counts['unresolved_memberships'] != counts['unique_source_memberships']
                or inventory['collection'] != 'single-validation-ingest-pass'
                or inventory['physical_counts'] != {k:counts[k] for k in
                    ('ccre_rows','re2g_rows','anchor_family_counts','rows_by_target_kind')}
                or inventory['target_resolution_complete'] != (counts['unresolved_memberships'] == 0)
                or type(inventory['target_resolution_complete']) is not bool):
            raise ValueError
        return scope
    except (ValueError, TypeError, KeyError, AttributeError):
        raise WgsReportingError('profile source accounting is incomplete or foreign') from None


def profile_source_accounting(index, catalog):
    """Summarize live authenticated source scope without reopening originals."""
    from genosapiens.regulatory_index import RegulatoryIndex
    from genosapiens.track1_analysis_profiles import ProfileCatalog
    if type(index) is not RegulatoryIndex or type(catalog) is not ProfileCatalog:
        raise WgsReportingError('profile source accounting requires authenticated owners')
    index.require_stable()
    # The independently sealed index manifest owns source-wide classification;
    # these are different denominators from a case's queried/advancing targets.
    return validate_profile_source_accounting({
        'schema': 'genosapiens.track1-profile-source-accounting', 'version': 1,
        'status': 'authenticated_source_scope',
        'regulatory_index_manifest_sha256': index.manifest_sha256,
        'regulatory_manifest': index.manifest,
        'regulatory_counts': index.manifest['counts'],
        'regulatory_inventory': index.manifest['inventory'],
        'catalog_sha256': catalog.sha256,
        'catalog_source_sha256': catalog.source_sha256,
        'catalog_records': len(catalog.records), 'catalog_genes': len(catalog.genes),
        'catalog_gene_disease_links': len(catalog.links),
        'limitations': ['source-wide counts are not candidate-level target counts',
            'complete parsing does not establish complete target identity resolution',
            'source links do not establish causal regulatory effects'],
    })


def scientific_work_counts(reader: ScientificArtifactReader) -> dict:
    """Measure serialized work using the scientific reader already in use.

    Inspect context/universe/pair headers; do not expand factorized pairs or
    open a second scientific reader just to collect telemetry.
    """
    from genosapiens.wgs_case_resources import bound_reference
    reader._require_open()
    artifact = reader.artifact
    counts = dict(contexts=artifact.context_count, positive_contexts=artifact.positive_context_count,
        context_members=artifact.member_count, logical_pairs=artifact.pair_count,
        hypotheses=artifact.hypothesis_count)
    pairs = outcomes = patterns = decisions = constructed = 0
    work_policy = reader.resources.parent_bindings.get('population_work_policy')
    work_keys = ('original_member_count','retained_member_count','excluded_member_count',
        'unresolved_member_count','original_member_occurrences','retained_member_occurrences',
        'excluded_member_occurrences')
    work_totals = dict.fromkeys(work_keys, 0)
    work_contexts = 0
    for row in reader.case.contexts():
        context = reader.reader.reconstruct(bound_reference(reader.reader, row['context'], 'scientific-case-context-v1'))
        members = reader.reader.reconstruct(bound_reference(reader.reader, context['members'], 'context-evidence-v1'))
        work = members.get('work_filter')
        if (work is None) != (work_policy is None):
            raise WgsReportingError('population work reporting policy differs')
        if work is not None:
            header = reader.reader.reconstruct(bound_reference(reader.reader, work, 'context-population-work-filter-v1'))
            if (header['policy_sha256'] != work_policy or header['resources'] != context['resources']
                    or header['context_ordinal'] != row['ordinal']
                    or header['context_id'] != context['context_id']
                    or header['retained_member_count'] != row['member_count']
                    or header['assessment_domain_sha256'] != context['assessment_domain_sha256']
                    or header['rich_excluded_hypothesis_count'] is not None
                    or header['skipped_scoring'] != 'not_evaluated'
                    or any(type(header[key]) is not int or header[key] < 0 for key in work_keys)
                    or header['original_member_count'] != header['retained_member_count'] + header['excluded_member_count']
                    or header['original_member_occurrences'] != header['retained_member_occurrences'] + header['excluded_member_occurrences']
                    or header['unresolved_member_count'] > header['retained_member_count']):
                raise WgsReportingError('population work reporting counts differ')
            for key in work_keys: work_totals[key] += header[key]
            work_contexts += 1
        universe = reader.reader.reconstruct(bound_reference(reader.reader, context['universe'], 'complete-hypothesis-universe-v1'))
        if universe['pair_evidence'] is not None:
            dependency = universe['pair_evidence']
            if type(dependency) is not dict or set(dependency) != {'name', 'sha256'}:
                raise WgsReportingError('invalid pair telemetry reference')
            reference = reader.reader.reference(dependency['name'], kind='complete-pair-evidence-v1')
            if reference.sha256 != dependency['sha256']:
                raise WgsReportingError('pair telemetry dependency differs')
            pair = reader.reader.reconstruct(reference)
            pairs += pair['pair_count']; outcomes += pair['outcome_count']; patterns += pair['pattern_count']
            decisions += pair['decision_count']; constructed += pair['constructible_count']
    if pairs != artifact.pair_count:
        raise WgsReportingError('scientific work counts differ from complete pair domain')
    counts.update(persisted_computed_pair_positions=pairs, computed_pair_outcome_classes=outcomes,
        computed_pair_patterns=patterns, computed_pair_decisions=decisions, constructed_pair_hypotheses=constructed)
    projections = reader.assessment_projections
    counts.update(discovery_targets=None if projections is None else len(projections[0].targets),
        discovery_selected_inputs=None if projections is None else projections[0].row_count)
    profile = reader.resources.analysis_profile
    if profile is not None:
        counts['retention_identities'] = profile['counts']['retention_dispositions']
        counts['retention_source_occurrences'] = profile['counts']['retention_occurrences']
    work_summary = {}
    if work_policy is not None:
        if work_contexts != artifact.context_count or work_totals['retained_member_count'] != artifact.member_count:
            raise WgsReportingError('population work reporting coverage differs')
        work_summary['population_work'] = dict(policy_sha256=work_policy, contexts=work_contexts,
            **work_totals, rich_excluded_hypothesis_count=None, skipped_scoring='not_evaluated',
            count_scope='Context member occurrences. These counts do not represent unique alleles or source records.')
    return dict(schema='genosapiens.track1-scientific-work-counts', version=2 if work_summary else 1,
        **work_summary,
        scientific_manifest_sha256=artifact.sha256,
        profile_manifest_sha256=None if profile is None else profile['manifest_sha256'],
        counts=counts, measurement_scope='Counts use existing scientific metadata without pair expansion.')


def retention_coverage_text(profile, source_scope=None):
    """Describe the evaluated universe without a negative partner inference."""
    if profile is None:
        return ''
    counts = profile['counts']
    retained = counts['retention_dispositions']['retained']
    deferred = counts['retention_dispositions']['deferred_no_retention_evidence']
    occurrences = counts['retention_occurrences']
    text = (f" The retained universe contains {retained} alleles and {occurrences['retained']} source occurrences."
        f" The selector deferred {deferred} alleles and {occurrences['deferred_no_retention_evidence']} source occurrences."
        " Deferred alleles did not enter contexts, pairs, scores, mechanisms, or drug inputs."
        " The compact deferred inventory remains in the ordered profile decision shards."
        " Selection can omit novel noncoding causes and unknown effects without qualifying evidence."
        " A deferred partner can exist in the source genome.")
    if 'finite_cohort' in profile['identity']:
        text += f" This finite input subset comes from {profile['identity']['parent']['rows']} original inputs."
    if source_scope is None:
        return text + ' Source-wide regulatory accounting remains unconfirmed.'
    return text + (f" The regulatory source contains {source_scope['regulatory_counts']['re2g_rows']} rE2G rows,"
        f" {source_scope['regulatory_counts']['distinct_raw_targets']} raw targets, and"
        f" {source_scope['regulatory_counts']['unresolved_memberships']} unresolved relations."
        f" The direct catalog contains {source_scope['catalog_genes']} genes."
        " Complete parsing does not establish target resolution or causality.")


def annotation_coverage_text(coverage):
    if coverage is None:
        return ''
    from genosapiens.wgs_annotation_coverage import REASON, POLICY, ALIAS_REASON, ALIAS_POLICY
    reasons = coverage.get('reason_counts')
    alias = reasons.get(ALIAS_REASON, 0) if type(reasons) is dict else 0
    if (coverage.get('schema') != 'genosapiens.wgs-annotation-coverage'
            or coverage.get('version') != 1 or coverage.get('policy') not in (POLICY, ALIAS_POLICY)
            or any(type(coverage.get(key)) is not int or coverage[key] < 0
                for key in ('sequence_requests', 'verified_responses', 'unassessed_requests'))
            or coverage['sequence_requests'] != coverage['verified_responses'] + coverage['unassessed_requests']
            or type(alias) is not int or not 0 <= alias <= coverage['unassessed_requests']
            or reasons != {REASON: coverage['unassessed_requests'] - alias,
                           **({ALIAS_REASON: alias} if alias else {})}
            or coverage['policy'] != (ALIAS_POLICY if alias else POLICY)
            or not _digest(coverage.get('annotation_manifest_sha256'))):
        raise WgsReportingError('annotation coverage summary is invalid')
    reason_text = ('Unassessed requests lack responses because their contigs are absent from the pinned VEP reference. '
        if not alias else
        f"{reasons[REASON]} " + ('request has a contig' if reasons[REASON] == 1 else 'requests have contigs')
        + ' absent from the pinned VEP reference. '
        + f"{alias} " + ('request has an unresolved nonprimary contig alias. ' if alias == 1
                        else 'requests have unresolved nonprimary contig aliases. ')
        +
        'A possible alias does not establish sequence equivalence or reference absence. ')
    return (f". Annotation coverage: {coverage['verified_responses']} verified responses and "
        f"{coverage['unassessed_requests']} unassessed requests. "
        + reason_text +
        'A contig is a named reference sequence. Missing annotation does not establish a negative finding.')


def build_bounded_reports(
    *,
    source_root: str | Path,
    source_manifest: LedgerManifest,
    source_accounting: SourceAccountingClosure,
    source_routing_context: SourceRoutingContext,
    qc_summary: VcfQcIncrementalSummary,
    normalization_summary: WgsNormalizationSummary,
    class_summary: IncrementalVariantClassSummary,
    filter_root: str | Path,
    filter_manifest: FilterManifest,
    working_root: str | Path,
    working_manifest: WorkingSetManifest,
    policy_sha256: str,
    scope_policy_sha256: str,
    scope_manifest_sha256: str,
    population_manifest_sha256: str,
    clinical_manifest_sha256: str,
    exception_manifest_sha256: str,
    precision_manifest_sha256: str,
    scientific: ScientificArtifact,
    scientific_directory: str | Path,
    order: ScientificOrderProduct,
    order_directory: str | Path,
    primary: ScientificPrimaryProduct,
    selection: ScientificSelectionProduct,
    selection_directory: str | Path,
    selected_technical_assessments: tuple[SelectedTechnicalAssessment, ...],
    credibility_threshold: str | None,
    execution_config: BoundedWgsExecutionConfig,
    output_directory: str | Path,
    maximum_total_json_bytes: int,
    maximum_total_markdown_bytes: int,
    scientific_reader: ScientificArtifactReader,
    read_intake=None,
    retained_filter_state_path: str | Path | None = None,
    profile_source_scope: dict | None = None,
    annotation_coverage: dict | None = None,
) -> PreparedBoundedReports:
    """Retain complete normalized science and stream source/selected audit detail."""

    json_writer = markdown_writer = None
    projection_scope = ExitStack()
    try:
        root = Path(output_directory)
        source_path = Path(source_root)
        if (
            not root.is_absolute()
            or root.exists()
            or not source_path.is_absolute()
            or type(execution_config) is not BoundedWgsExecutionConfig
            or type(maximum_total_json_bytes) is not int
            or maximum_total_json_bytes <= 0
            or type(maximum_total_markdown_bytes) is not int
            or maximum_total_markdown_bytes <= 0
        ):
            raise ValueError
        if credibility_threshold is not None:
            _exact(credibility_threshold)
        if (type(scientific_reader) is not ScientificArtifactReader
                or scientific_reader.artifact != scientific
                or scientific_reader.directory != Path(scientific_directory)):
            raise ValueError
        scientific_reader._require_open()
        coverage_text = annotation_coverage_text(annotation_coverage)
        if profile_source_scope is not None:
            component = scientific_reader.resources.analysis_profile
            if component is None:
                raise ValueError
            profile_source_scope = validate_profile_source_accounting(profile_source_scope,
                index_sha256=component['identity']['regulatory_index_sha256'],
                catalog_sha256=component['identity']['catalog_sha256'])
        projection_scope.enter_context(scientific_reader.reader.read_scope())
        validate_scientific_products(scientific_reader,
            ScientificProducts(order, primary, selection),
            expected_prefilter_manifest_sha256=filter_manifest.sha256)
        root.mkdir(mode=0o700, parents=True)
        _verify_source_summaries(
            source_path,
            source_manifest,
            source_accounting,
            qc_summary,
            normalization_summary,
            class_summary,
        )
        working_summary = verify_working_set(
            source_root=source_path,
            source_manifest=source_manifest,
            source_accounting=source_accounting,
            source_routing_context=source_routing_context,
            filter_root=filter_root,
            filter_manifest=filter_manifest,
            working_root=working_root,
            working_manifest=working_manifest,
            policy_sha256=policy_sha256,
            scope_policy_sha256=scope_policy_sha256,
            scope_manifest_sha256=scope_manifest_sha256,
            population_manifest_sha256=population_manifest_sha256,
            clinical_manifest_sha256=clinical_manifest_sha256,
            exception_manifest_sha256=exception_manifest_sha256,
            precision_manifest_sha256=precision_manifest_sha256,
            retained_state_path=retained_filter_state_path,
        )
        if (
            dict(scientific.parent_bindings).get("prefilter_manifest")
            != filter_manifest.sha256
            or order.scientific_manifest_sha256 != scientific.sha256
            or primary.scientific_manifest_sha256 != scientific.sha256
            or primary.order_manifest_sha256 != order.sha256
            or primary.prefilter_manifest_sha256 != filter_manifest.sha256
            or selection.scientific_manifest_sha256 != scientific.sha256
            or selection.order_manifest_sha256 != order.sha256
            or selection.ranking_policy_sha256 != scientific.policy_sha256
        ):
            raise ValueError
        selected_technical_assessments = _typed(
            selected_technical_assessments,
            tuple[SelectedTechnicalAssessment, ...],
        )
        selected_technical_assessments = tuple(sorted(
            selected_technical_assessments,
            key=lambda value: (
                value.allele_id,
                value.source_allele_id,
                value.sample_name,
                value.reason,
                value.scope,
            ),
        ))
        if len(set(selected_technical_assessments)) != len(
            selected_technical_assessments
        ):
            raise ValueError
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
        if any(
            source.authority == "supplied_source"
            and source.source_input_sha256 not in (
                {row.input_sha256 for row in read_intake.manifest.source_bindings} if read_intake
                else {source_manifest.source.input_sha256})
            for group in selection.selected
            for trace in group.alleles
            for source in trace.source_traces
        ):
            raise ValueError
        assessment_sources = {
            (
                value.allele_id,
                value.source_allele_id,
                value.sample_name,
                value.source_input_sha256,
            )
            for value in selected_technical_assessments
        }
        if assessment_sources - selected_sources:
            raise ValueError
        selected_obligation_status = (
            "not_applicable" if not selected_sources
            else "unmet" if selected_sources - assessment_sources or any(
                value.status == "unmet"
                for value in selected_technical_assessments
            ) else "not_applicable" if all(
                value.status == "not_applicable"
                for value in selected_technical_assessments
            ) else "satisfied"
        )

        limits = execution_config.source_shards
        json_writer = _ReportWriter(
            root,
            "audit_jsonl",
            "jsonl",
            maximum_rows=limits.max_uncompressed_rows,
            maximum_bytes=limits.max_uncompressed_bytes,
            maximum_total_bytes=maximum_total_json_bytes,
        )
        markdown_writer = _ReportWriter(
            root,
            "audit_markdown",
            "md",
            maximum_rows=limits.max_uncompressed_rows,
            maximum_bytes=limits.max_uncompressed_bytes,
            maximum_total_bytes=maximum_total_markdown_bytes,
        )

        # Completed preparation already owns record/ALT accounting and exact
        # source mappings. Keep its addressable descriptors, not two newly
        # rendered genome-wide copies of the same per-record information.
        source_rows = source_accounting.total_records
        key = 'source:accounting'
        json_writer.write(_report_row('source_accounting', key, {
            'source_manifest_sha256': source_manifest.sha256,
            'source_root': str(source_path),
            'source_manifest': asdict(source_manifest),
            'accounting': asdict(source_accounting),
            'filter_root': str(filter_root),
            'filter_manifest': asdict(filter_manifest),
            'filter_manifest_sha256': filter_manifest.sha256,
            'working_manifest_sha256': working_manifest.sha256,
            'effective_preparation_counts': asdict(working_summary),
            'scope_deferred': dict(annotation_status='not_performed',
                excluded_scope_occurrences=working_summary.excluded_scope_occurrences,
                unresolved_scope_occurrences=working_summary.unresolved_scope_occurrences),
            'request_count_scope': 'successfully committed preparation batches',
            'detail_location': 'Original source and filter ledgers preserve source allele identifiers and shard record bounds.',
            'fresh_source_replay': False,
        }), key=key)
        markdown_writer.write(_markdown_row('source accounting', key,
            f'The input contains {source_rows} original records and {source_accounting.total_physical_alternate_slots} physical ALT slots. '
            'JSON descriptors identify the source and filter ledgers. '
            'Effective scope deferrals have annotation_status=not_performed.'), key=key)

        representation_rows = 0
        for row in scientific_reader.iter_representation_records():
            representation_rows += 1
            assessment = row["assessment"]
            key = f"representation:{row['context_ordinal']}:{assessment.allele_id}"
            json_writer.write(_report_row("representation", key, row), key=key)
            markdown_writer.write(_markdown_row("representation", key,
                f"Status: {assessment.status}. Reasons: {','.join(assessment.reasons) or 'none'}"), key=key)
        if representation_rows != scientific.member_count:
            raise ValueError

        # Full logical counts refer to authenticated normalized roots. Physical
        # audit shard counts refer only to the rows actually written below.
        hypothesis_rows = scientific.hypothesis_count
        credible = (0 if credibility_threshold is None else
                    scientific_reader.credible_complete_count(credibility_threshold))
        scientific_roots = scientific_report_roots(scientific_reader,
            ScientificProducts(order, primary, selection))
        from genosapiens.wgs_science_reporting import scope_text
        scientific_scope = getattr(scientific_reader, 'scientific_scope', None)
        scope_header = None if scientific_scope is None else scientific_scope.header
        key = "scientific:complete"
        json_writer.write(_report_row("normalized_scientific_evidence", key, {
            "scientific_manifest_sha256": scientific.sha256,
            **({} if scope_header is None else {"scientific_scope": scope_header}),
            "roots": scientific_roots, "hypothesis_count": hypothesis_rows,
            "optional_analysis_status": "not_run",
            "assessment_domain": scientific_reader.resources.assessment_domain,
            **({} if annotation_coverage is None else {'annotation_coverage': annotation_coverage}),
            "products": ScientificProducts(order, primary, selection),
            **({} if scientific_reader.resources.analysis_profile is None else {
                "analysis_profile":scientific_reader.resources.analysis_profile,
                "profile_source_scope": profile_source_scope if profile_source_scope is not None else {
                    'status': 'unconfirmed', 'reason': 'source-wide accounting was not supplied'},
            }),
            **({} if scientific_reader.resources.analysis_eligibility is None else {
                "analysis_eligibility": scientific_reader.resources.analysis_eligibility,
                "analysis_domain_status": scientific_reader.resources.analysis_eligibility["analysis_domain_status"],
                "rescue_assessment_status": scientific_reader.resources.analysis_eligibility["rescue_assessment_status"],
            }),
        }), key=key)
        markdown_writer.write(_markdown_row("evaluated scientific evidence", key,
            ("" if scope_header is None else scope_text(scope_header) + " ") +
            f"The scientific output contains {hypothesis_rows} screened hypotheses. "
            "Optional analyses did not run. Complete orders and dispositions remain in normalized evidence."
            + ("" if scientific_reader.resources.analysis_eligibility is None else
                " Earlier population and RefCall decisions remain separate. Deferred rescue did not run.")
            + retention_coverage_text(scientific_reader.resources.analysis_profile, profile_source_scope)
            + coverage_text), key=key)

        selected_by_rank = {group.source_rank: group for group in selection.selected}
        selected_seen = set()
        detailed_hypothesis_rows = 0
        for row in scientific_reader.iter_selected_records():
            ordered = row["ordering"]
            group = selected_by_rank.get(ordered["rank"])
            if (group is None or group.source_rank in selected_seen
                    or group.primary_hypothesis_id != ordered["hypothesis_id"]
                    or group.source_score != row["ledger"].total
                    or group.alleles != row["allele_traces"]):
                raise ValueError
            selected_seen.add(group.source_rank)
            detailed_hypothesis_rows += 1
            key = f"hypothesis:{ordered['rank']}"
            json_writer.write(_report_row("selected_hypothesis", key, {
                "context_id": row["context_id"], "hypothesis": row["constructed"],
                "evidence": row["evidence"], "score_ledger": row["ledger"],
                "population_assessment": row["population_assessment"],
                "population_policy": row["population_policy"],
                "eligibility_reasons": row["eligibility_reasons"], "order": ordered,
                "selection": group, "partial_member_status":
                    "incomplete_partner_unresolved" if group.partial_member else "not_applicable",
            }), key=key)
            markdown_writer.write(_markdown_row("selected hypothesis", key,
                f"Submitted. Total: {ordered['total']}. "
                f"Genotype: {row['constructed'].hypothesis.completeness.value}. "
                f"partial_member={'yes' if group.partial_member else 'no'}"), key=key)
        if selected_seen != set(selected_by_rank):
            raise ValueError

        blockers = []
        if read_intake is not None:
            blockers.extend(read_intake.manifest.mandatory_stage_blockers)
            if read_intake.manifest.competitive_eligibility == "ineligible":
                blockers.append("read-intake-not-competitively-eligible")
        warnings = [
            "early-filter-can-change-candidates-and-rankings",
            "downstream-population-threshold-is-synthetic-development-separator",
            "private-truth-unavailable-not-evaluated",
        ]
        if annotation_coverage is not None and annotation_coverage['unassessed_requests']:
            warnings.append('annotation-reference-coverage-gap-unassessed-records-retained')
        if not selection.selected:
            blockers.append("zero-eligible-hypotheses")
        if selected_obligation_status == "unmet":
            blockers.append("selected-allele-mandatory-technical-check-unmet")
        diagnostic = BoundedDiagnosticReport(
            source_manifest.sha256,
            hashlib.sha256(_bytes(source_accounting)).hexdigest(),
            source_accounting.qc_summary_sha256,
            source_accounting.normalization_summary_sha256,
            source_accounting.class_summary_sha256,
            filter_manifest.sha256,
            working_manifest.sha256,
            scientific.sha256,
            order.sha256,
            primary.sha256,
            selection.sha256,
            source_accounting.total_records,
            source_accounting.total_physical_alternate_slots,
            working_summary.retained_occurrences,
            working_summary.excluded_common_occurrences,
            hypothesis_rows,
            len(selection.selected),
            "complete",
            "complete",
            "present" if credible else (
                "not_evaluated" if credibility_threshold is None else "unknown"
            ),
            credible,
            selected_obligation_status,
            selected_technical_assessments,
            tuple(sorted(blockers)),
            tuple(sorted(warnings)),
            scientific_reader.resources.assessment_domain_sha256,
            scientific_reader.resources.assessment_domain['original_inputs'],
            (None if scientific_reader.resources.assessment_domain['input_dispositions'] is None else
                tuple(sorted(scientific_reader.resources.assessment_domain['input_dispositions'].items()))),
            scope_policy_sha256=scope_policy_sha256, scope_manifest_sha256=scope_manifest_sha256,
            scope_counts=tuple(sorted((row.name, row.count) for row in filter_manifest.summary.scope_counts)),
            in_scope_normalized_occurrence_count=filter_manifest.summary.in_scope_normalized_occurrences,
            retained_identity_count=working_summary.retained_unique_identities,
            excluded_scope_occurrence_count=working_summary.excluded_scope_occurrences,
            unresolved_scope_occurrence_count=working_summary.unresolved_scope_occurrences,
            bypassed_occurrence_count=working_summary.bypassed_occurrences,
            population_request_keys=filter_manifest.summary.population_request_keys,
            clinical_query_keys=filter_manifest.summary.clinical_query_keys,
        )
        json_shards = json_writer.close()
        markdown_shards = markdown_writer.close()
        total_audit_rows = (
            1 + representation_rows + 1 + detailed_hypothesis_rows
        )
        parents = tuple(sorted({
            "class_summary": source_accounting.class_summary_sha256,
            "filter_manifest": filter_manifest.sha256,
            "order_manifest": order.sha256,
            "primary_manifest": primary.sha256,
            "scientific_manifest": scientific.sha256,
            "scientific_evidence": scientific.evidence.sha256,
            "selection": selection.sha256,
            "source_manifest": source_manifest.sha256,
            "working_manifest": working_manifest.sha256,
            **({} if annotation_coverage is None else {
                'annotation_coverage': hashlib.sha256(canonical(annotation_coverage)).hexdigest()}),
            **({"read_intake": read_intake.manifest.sha256} if read_intake else {}),
        }.items()))
        audit = BoundedAuditManifest(
            diagnostic.sha256,
            parents,
            json_shards,
            markdown_shards,
            source_rows,
            representation_rows,
            hypothesis_rows,
            total_audit_rows,
            limits.max_uncompressed_rows,
            limits.max_uncompressed_bytes,
            maximum_total_json_bytes,
            maximum_total_markdown_bytes,
            json_writer.stream.hexdigest(),
            markdown_writer.stream.hexdigest(),
            hypothesis_rows == 0,
            scientific_roots,
            detailed_hypothesis_rows,
            scientific_reader.resources.assessment_domain_sha256,
        )
        diagnostic_bytes = diagnostic.as_bytes()
        audit_bytes = audit.as_bytes()
        work = scientific_work_counts(scientific_reader)
        (root / 'runtime-work-counts.json').write_bytes(canonical(work))
        (root / "diagnostics.json").write_bytes(diagnostic_bytes)
        (root / "audit-manifest.json").write_bytes(audit_bytes)
        (root / "audit-index.md").write_bytes(audit.markdown_index())
        for path in (
            root / 'runtime-work-counts.json',
            root / "diagnostics.json",
            root / "audit-manifest.json",
            root / "audit-index.md",
        ):
            path.chmod(0o600)
        return PreparedBoundedReports(
            diagnostic,
            audit,
            audit_bytes,
            audit.markdown_index(),
            root,
        )
    except WgsReportingError:
        raise
    except _ERRORS:
        raise WgsReportingError("bounded WGS report construction failed") from None
    finally:
        if json_writer is not None:
            json_writer.abort()
        if markdown_writer is not None:
            markdown_writer.abort()
        projection_scope.close()


def _verify_bounded_audit_metadata(report: PreparedBoundedReports) -> None:
    """Authenticate small report descriptors only, never shard payloads."""
    try:
        if type(report) is not PreparedBoundedReports:
            raise ValueError
        for name, expected in (
            ("diagnostics.json", report.diagnostic.as_bytes()),
            ("audit-manifest.json", report.audit.as_bytes()),
            ("audit-index.md", report.audit.markdown_index()),
        ):
            _read_exact_private(report.output_directory / name, expected)
    except _ERRORS:
        raise WgsReportingError("bounded WGS report replay failed") from None


class _AuditStreamVerifier:
    """One copy/read owner, one current shard and ordered-stream hash.

    Private shard copies are provisional until finish closes both streams. The
    reader supplies its actual per-shard digest; this verifier never rehashes it.
    """

    def __init__(self, report: PreparedBoundedReports):
        if type(report) is not PreparedBoundedReports:
            raise WgsReportingError("bounded WGS audit stream is invalid")
        self.audit = report.audit
        self.streams = (
            (self.audit.json_shards, self.audit.json_stream_sha256, self.audit.maximum_total_json_bytes),
            (self.audit.markdown_shards, self.audit.markdown_stream_sha256, self.audit.maximum_total_markdown_bytes),
        )
        paths = [shard.relative_path for shards, _, _ in self.streams for shard in shards]
        if len(paths) != len(set(paths)):
            raise WgsReportingError("bounded WGS audit shard sequence is invalid")
        self.stream_index = self.shard_index = self.total = 0
        self.stream = hashlib.sha256()
        self.current = None
        self.size = self.rows = 0
        self.finished = False

    def _advance_streams(self):
        while self.stream_index < len(self.streams):
            shards, expected, maximum = self.streams[self.stream_index]
            if self.shard_index < len(shards):
                return
            if self.total > maximum or self.stream.hexdigest() != expected:
                raise WgsReportingError("bounded WGS audit stream closure failed")
            self.stream_index += 1
            self.shard_index = self.total = 0
            self.stream = hashlib.sha256()

    def begin_shard(self, shard: ReportShard) -> None:
        if self.finished or self.current is not None:
            raise WgsReportingError("bounded WGS audit shard sequence is invalid")
        self._advance_streams()
        if (self.stream_index == len(self.streams)
                or shard != self.streams[self.stream_index][0][self.shard_index]):
            raise WgsReportingError("bounded WGS audit shard sequence is invalid")
        self.current = shard
        self.size = self.rows = 0

    def update(self, content: bytes) -> None:
        if self.current is None or self.finished or type(content) is not bytes:
            raise WgsReportingError("bounded WGS audit stream is invalid")
        self.size += len(content)
        self.total += len(content)
        self.rows += content.count(b"\n")
        if (self.size > min(self.current.size_bytes, self.audit.maximum_bytes_per_shard)
                or self.rows > min(self.current.row_count, self.audit.maximum_rows_per_shard)
                or self.total > self.streams[self.stream_index][2]):
            raise WgsReportingError("bounded WGS audit stream exceeds capacity")
        self.stream.update(content)

    def finish_shard(self, *, sha256: str, size_bytes: int) -> None:
        shard = self.current
        if (shard is None or self.finished or type(size_bytes) is not int
                or sha256 != shard.sha256 or size_bytes != self.size
                or self.size != shard.size_bytes or self.rows != shard.row_count):
            raise WgsReportingError("bounded WGS audit shard closure failed")
        self.current = None
        self.shard_index += 1
        self._advance_streams()

    def finish(self) -> None:
        if self.finished or self.current is not None:
            raise WgsReportingError("bounded WGS audit stream is incomplete")
        self._advance_streams()
        if self.stream_index != len(self.streams):
            raise WgsReportingError("bounded WGS audit stream is incomplete")
        self.finished = True


@contextmanager
def _open_audit_payload(path: Path, *, expected_size: int, mode: int = 0o600):
    """Read a private regular file, checking named/open identity through EOF."""
    try:
        if any(parent.is_symlink() for parent in path.parents):
            raise ValueError
        descriptor = os.open(path, os.O_RDONLY | os.O_NOFOLLOW | os.O_NONBLOCK)
        with os.fdopen(descriptor, "rb") as stream:
            initial = os.fstat(stream.fileno())
            def identity(value):
                return (value.st_dev, value.st_ino, value.st_size, value.st_mtime_ns,
                        value.st_ctime_ns, value.st_mode, value.st_uid, value.st_nlink)
            if (not stat.S_ISREG(initial.st_mode) or initial.st_uid != os.geteuid()
                    or initial.st_nlink != 1 or stat.S_IMODE(initial.st_mode) != mode
                    or initial.st_size != expected_size
                    or identity(initial) != identity(path.lstat())):
                raise ValueError
            yield stream
            if identity(initial) != identity(os.fstat(stream.fileno())) or identity(initial) != identity(path.lstat()):
                raise ValueError
    except (OSError, ValueError):
        raise WgsReportingError("bounded WGS audit payload is unsafe or changed") from None


def verify_bounded_audit(report: PreparedBoundedReports) -> None:
    """Explicit standalone audit using the same checks as publication copying."""
    try:
        _verify_bounded_audit_metadata(report)
        verifier = _AuditStreamVerifier(report)
        for shards in (report.audit.json_shards, report.audit.markdown_shards):
            for shard in shards:
                verifier.begin_shard(shard)
                digest = hashlib.sha256()
                size = 0
                with _open_audit_payload(report.output_directory / shard.relative_path,
                        expected_size=shard.size_bytes) as handle:
                    for content in iter(lambda: handle.read(65536), b""):
                        digest.update(content)
                        size += len(content)
                        verifier.update(content)
                verifier.finish_shard(sha256=digest.hexdigest(), size_bytes=size)
        verifier.finish()
    except _ERRORS:
        raise WgsReportingError("bounded WGS report replay failed") from None
