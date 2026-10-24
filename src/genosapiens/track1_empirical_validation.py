"""Fail-closed loader for the public GIAB, caller, and VEP validation record."""

from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass
from fractions import Fraction
import hashlib
from pathlib import Path
import re
from typing import Any

from genosapiens.bounded_json import BoundedJSONError, load_bounded_json
from genosapiens.track1_roster import load_track1_roster


CANONICAL_REPORT_PATH = "artifacts/validation/track1-public-empirical-current-2026-09-12.json"
CANONICAL_REPORT_SHA256 = "f9afd84e5989a8e0294098c3ac95e4978bfa86d601e65e383aae3eb1c2918325"
SCHEMA = "genosapiens.track1-public-empirical-validation"
SHA256 = re.compile(r"^[0-9a-f]{64}$")


class EmpiricalValidationError(RuntimeError):
    """Raised when public validation evidence is incomplete or stale."""


@dataclass(frozen=True, slots=True)
class EmpiricalValidation:
    report_sha256: str
    giab_precision: Fraction
    giab_recall: Fraction
    caller_records: int
    roster_cases: int
    accepted: bool


@dataclass(frozen=True, slots=True)
class EmpiricalComponentEvidence:
    source_report_sha256: str
    applicable_components: tuple[str, ...]
    giab_precision: Fraction
    giab_recall: Fraction
    caller_records: int
    generating_roster_cases: int
    current_roster_annotation_accepted: bool = False


def _closed(value: object, fields: set[str], label: str) -> Mapping[str, Any]:
    if not isinstance(value, Mapping) or set(value) != fields:
        raise EmpiricalValidationError(f"{label} fields are invalid")
    return value


def _sha(value: object, label: str) -> str:
    if type(value) is not str or SHA256.fullmatch(value) is None:
        raise EmpiricalValidationError(f"{label} is not a SHA-256")
    return value


def _bound_repository_file(root: Path, path: object, digest: object) -> None:
    if type(path) is not str:
        raise EmpiricalValidationError("public resource path is invalid")
    target = root / path
    if ".." in Path(path).parts or Path(path).is_absolute() or target.is_symlink():
        raise EmpiricalValidationError("public resource path cannot be a symbolic link")
    try:
        content = target.read_bytes()
    except OSError as error:
        raise EmpiricalValidationError("public resource cannot be read") from error
    if hashlib.sha256(content).hexdigest() != _sha(digest, "resource contract SHA-256"):
        raise EmpiricalValidationError("public resource binding is stale")


def load_empirical_validation(
    path: str | Path = CANONICAL_REPORT_PATH,
    *,
    repository_root: str | Path,
    expected_sha256: str = CANONICAL_REPORT_SHA256,
    _components_only: bool = False,
) -> EmpiricalValidation:
    root = Path(repository_root).resolve(strict=True)
    source = Path(path)
    if not source.is_absolute():
        source = root / source
    try:
        document, content = load_bounded_json(
            source, reference_root=root, maximum_bytes=64 * 1024
        )
    except BoundedJSONError as error:
        raise EmpiricalValidationError("public validation report cannot be read") from error
    if hashlib.sha256(content).hexdigest() != _sha(expected_sha256, "report SHA-256"):
        raise EmpiricalValidationError("public validation report is stale")
    report = _closed(
        document,
        {
            "schema", "classification", "validation_only", "controlled_sources_used",
            "execution_date", "report_id", "runtime", "resources",
            "full_resource_validation", "giab_representation", "giab_annotation",
            "caller_preparation", "roster_annotation", "limitations",
        },
        "public validation report",
    )
    runtime = _closed(
        report["runtime"],
        {"network", "user", "references", "scratch", "controlled_mounts"},
        "validation runtime",
    )
    if (
        report["schema"] != SCHEMA
        or report["classification"] != "public"
        or report["validation_only"] is not True
        or report["controlled_sources_used"] is not False
        or runtime["network"] != "none-during-validation"
        or runtime["references"] != "read-only"
        or runtime["controlled_mounts"] != []
    ):
        raise EmpiricalValidationError("public validation boundary is invalid")

    resources = report["resources"]
    if not isinstance(resources, list) or len(resources) != 3:
        raise EmpiricalValidationError("public validation resources are incomplete")
    expected_roles = {
        "giab-hg002-nist-v4.2.1-grch38": "genomic-representation-only",
        "deepvariant-hg002-novaseq-35x-v1.0.0-grch38": "caller-fixture",
        "ensembl-vep-116-grch38": "functional-annotation",
    }
    for value in resources:
        item = _closed(
            value,
            {
                "resource_id", "contract_path", "contract_sha256",
                "installed_manifest_sha256", "release", "oracle_role",
            },
            "public validation resource",
        )
        if expected_roles.get(item["resource_id"]) != item["oracle_role"]:
            raise EmpiricalValidationError("public validation oracle roles are invalid")
        _sha(item["installed_manifest_sha256"], "installed manifest SHA-256")
        _bound_repository_file(root, item["contract_path"], item["contract_sha256"])

    full = _closed(report["full_resource_validation"], {"giab", "deepvariant"}, "full validation")
    giab_full = _closed(
        full["giab"],
        {"variant_records", "confident_intervals", "indexed_contigs", "reference_sequences", "status"},
        "GIAB full validation",
    )
    caller_full = _closed(
        full["deepvariant"],
        {"variant_records", "pass_records", "filtered_records", "indexed_contigs", "status"},
        "caller full validation",
    )
    if (
        giab_full != {
            "variant_records": 4_048_342, "confident_intervals": 481_622,
            "indexed_contigs": 22, "reference_sequences": 195, "status": "succeeded",
        }
        or caller_full["variant_records"] != caller_full["pass_records"] + caller_full["filtered_records"]
        or caller_full["status"] != "succeeded"
        or caller_full["indexed_contigs"] != 25
    ):
        raise EmpiricalValidationError("full public resource validation is incomplete")

    representation = _closed(
        report["giab_representation"],
        {"scope", "tool", "criteria", "denominators", "outcomes", "discrepancies", "artifact_sha256s"},
        "GIAB representation",
    )
    denominators = _closed(representation["denominators"], {"caller_alleles", "truth_alleles"}, "GIAB denominators")
    outcomes = _closed(
        representation["outcomes"],
        {"shared_alleles", "caller_only_alleles", "truth_only_alleles", "precision", "recall", "status"},
        "GIAB outcomes",
    )
    criteria = _closed(representation["criteria"], {"established_on", "minimum_precision", "minimum_recall"}, "GIAB criteria")
    precision = Fraction(outcomes["precision"])
    recall = Fraction(outcomes["recall"])
    if (
        outcomes["shared_alleles"] + outcomes["caller_only_alleles"] != denominators["caller_alleles"]
        or outcomes["shared_alleles"] + outcomes["truth_only_alleles"] != denominators["truth_alleles"]
        or precision != Fraction(outcomes["shared_alleles"], denominators["caller_alleles"])
        or recall != Fraction(outcomes["shared_alleles"], denominators["truth_alleles"])
        or precision < Fraction(criteria["minimum_precision"])
        or recall < Fraction(criteria["minimum_recall"])
        or outcomes["status"] != "accepted"
        or criteria["established_on"] != "development-slice-before-candidate-transition"
    ):
        raise EmpiricalValidationError("GIAB representation result is invalid")
    artifact_hashes = representation["artifact_sha256s"]
    if not isinstance(artifact_hashes, Mapping) or len(artifact_hashes) != 7:
        raise EmpiricalValidationError("GIAB artifact accounting is incomplete")
    for digest in artifact_hashes.values():
        _sha(digest, "GIAB artifact SHA-256")

    annotation = _closed(
        report["giab_annotation"],
        {"scope", "oracle_role", "input_records", "output_records", "tool", "criteria", "outcomes", "artifact_sha256s", "limitation"},
        "GIAB annotation",
    )
    annotation_outcome = _closed(annotation["outcomes"], {"expected_outputs", "semantic_replay_equal", "status"}, "annotation outcome")
    annotation_hashes = annotation["artifact_sha256s"]
    if (
        annotation["oracle_role"] != "release-specific-functional-expectations"
        or annotation["input_records"] != annotation["output_records"]
        or annotation["output_records"] != annotation_outcome["expected_outputs"]
        or annotation_outcome["semantic_replay_equal"] is not True
        or annotation_outcome["status"] != "accepted"
        or not isinstance(annotation_hashes, Mapping)
        or len(annotation_hashes) != 4
    ):
        raise EmpiricalValidationError("VEP annotation result is invalid")
    for digest in annotation_hashes.values():
        _sha(digest, "annotation artifact SHA-256")

    caller = _closed(
        report["caller_preparation"],
        {"scope", "concordance", "qc", "preparation", "missing_field_behavior", "status"},
        "caller preparation",
    )
    concordance = caller["concordance"]
    qc = caller["qc"]
    prepared = caller["preparation"]
    if (
        not isinstance(concordance, Mapping)
        or concordance.get("records_accepted") != 328
        or concordance.get("records_rejected") != 0
        or not isinstance(qc, Mapping)
        or qc.get("records") != qc.get("pass_records", 0) + qc.get("filtered_records", 0)
        or qc.get("alternate_alleles") != 330
        or not isinstance(prepared, Mapping)
        or prepared.get("normalized_alleles") + prepared.get("quarantined_alleles") != qc.get("alternate_alleles")
        or caller["status"] != "accepted-for-bounded-early-workload"
    ):
        raise EmpiricalValidationError("caller preparation result is invalid")
    _sha(prepared.get("low_balance_policy_sha256"), "low-balance policy SHA-256")

    roster_result = _closed(
        report["roster_annotation"],
        {"roster_sha256", "expectations_sha256", "causal_truth_sha256", "case_count", "development_cases", "held_out_cases", "aggregate_runs"},
        "roster annotation",
    )
    roster = None
    if not _components_only:
        try:
            roster = load_track1_roster(
                "tests/fixtures/track1/roster/roster.json", repository_root=root
            )
        except Exception:
            raise EmpiricalValidationError("current roster annotation binding is invalid") from None
    for key in ("roster_sha256", "expectations_sha256", "causal_truth_sha256"):
        _sha(roster_result[key], "generating roster annotation identity")
    runs = roster_result["aggregate_runs"]
    if (
        roster is not None and (
            roster.roster_sha256 != roster_result["roster_sha256"]
            or roster.annotation_sha256 != roster_result["expectations_sha256"]
            or roster.causal_truth_sha256 != roster_result["causal_truth_sha256"]
            or roster_result["case_count"] != len(roster.cases)
        )
        or roster_result["case_count"] != 22
        or roster_result["development_cases"] != 11
        or roster_result["held_out_cases"] != 11
        or not isinstance(runs, list)
        or len(runs) != 2
        or {run.get("split") for run in runs if isinstance(run, Mapping)} != {"development", "held_out"}
    ):
        raise EmpiricalValidationError("roster annotation binding is invalid")
    for run in runs:
        _sha(run.get("raw_output_sha256"), "roster annotation output SHA-256")
        if run.get("input_records") != run.get("output_records") + run.get("unsupported_records"):
            raise EmpiricalValidationError("roster annotation attempt accounting is incomplete")

    limitations = report["limitations"]
    if not isinstance(limitations, list) or len(limitations) < 3 or not all(
        type(item) is str and item for item in limitations
    ):
        raise EmpiricalValidationError("public validation limitations are incomplete")
    return EmpiricalValidation(
        hashlib.sha256(content).hexdigest(), precision, recall,
        int(qc["records"]), int(roster_result["case_count"]), not _components_only,
    )


def load_empirical_component_evidence(
    path: str | Path = CANONICAL_REPORT_PATH,
    *, repository_root: str | Path, expected_sha256: str = CANONICAL_REPORT_SHA256,
) -> EmpiricalComponentEvidence:
    """Authenticate applicable GIAB/caller facts without rebinding the old roster.

    The historical report keeps its generating 22-case scope. Current-roster
    annotation is assessed separately and is never accepted by this adapter.
    """
    result = load_empirical_validation(
        path, repository_root=repository_root, expected_sha256=expected_sha256,
        _components_only=True,
    )
    return EmpiricalComponentEvidence(
        result.report_sha256,
        ("giab-genomic-representation", "giab-annotation", "public-caller-preparation"),
        result.giab_precision, result.giab_recall, result.caller_records, result.roster_cases,
    )
