"""Truth-isolated Exomiser runner and PhEval-compatible output adapter.

Only public HPO and VCF bytes enter the ranking process. Causal labels are
accepted solely by the separate comparison function after raw output parsing.
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence
import csv
from dataclasses import asdict, dataclass
import hashlib
import io
import json
from pathlib import Path
import re
import subprocess
from typing import Any

from genosapiens.bounded_json import BoundedJSONError, load_bounded_json


SCHEMA = "genosapiens.exomiser-baseline"
REPORT_SCHEMA = "genosapiens.exomiser-baseline-run"
VALIDATION_SCHEMA = "genosapiens.track1-exomiser-validation"
CANONICAL_VALIDATION_PATH = "artifacts/validation/track1-exomiser-2026-09-12.json"
CANONICAL_VALIDATION_SHA256 = (
    "a216298c07d6ccc78865d299258741e542bbc6736fefe8274a86c03e00ea8bb4"
)
HPO_ID = re.compile(r"^HP:[0-9]{7}$")
SHA256 = re.compile(r"^[0-9a-f]{64}$")
IMAGE = (
    "docker.io/exomiser/exomiser-cli@"
    "sha256:a368ec91266c2719db03d5d0eb24cfc4504d4f3a60da76cff7bb550d2a2a5d57"
)
FORBIDDEN_RUNNER_KEYS = ("causal", "oracle", "recipe", "seed", "truth")
GENE_FIELDS = (
    "gene_symbol",
    "gene_identifier",
    "score",
    "grouping_id",
)
VARIANT_FIELDS = (
    "chrom",
    "start",
    "end",
    "ref",
    "alt",
    "score",
    "grouping_id",
)


class ExomiserError(RuntimeError):
    """Raised for a baseline contract or execution failure."""


def _sha256(content: bytes) -> str:
    return hashlib.sha256(content).hexdigest()


def _digest(value: object) -> str:
    return _sha256(
        (json.dumps(value, sort_keys=True, separators=(",", ":")) + "\n").encode(
            "ascii"
        )
    )


def _require_digest(value: object, label: str) -> str:
    if type(value) is not str or SHA256.fullmatch(value) is None:
        raise ExomiserError(f"{label} is not a SHA-256")
    return value


def _closed_mapping(
    value: object, fields: set[str], label: str
) -> Mapping[str, Any]:
    if not isinstance(value, Mapping) or set(value) != fields:
        raise ExomiserError(f"{label} fields are invalid")
    return value


@dataclass(frozen=True, slots=True)
class ExomiserStack:
    manifest_sha256: str
    version: str
    data_release: str
    image: str
    analysis_path: str
    analysis_sha256: str
    installed_root: str
    installed_file_manifest_sha256: str
    installed_file_count: int
    installed_bytes: int
    code_license: str
    data_use_decision: str
    pheval_version: str
    pheval_commit: str


@dataclass(frozen=True, slots=True)
class ExomiserGeneRow:
    gene_symbol: str
    gene_identifier: str
    score: float
    grouping_id: str


@dataclass(frozen=True, slots=True)
class ExomiserVariantRow:
    chrom: str
    start: int
    end: int
    ref: str
    alt: str
    score: float
    grouping_id: str


@dataclass(frozen=True, slots=True)
class ExomiserOutputs:
    genes: tuple[ExomiserGeneRow, ...]
    variants: tuple[ExomiserVariantRow, ...]


@dataclass(frozen=True, slots=True)
class BaselineComparison:
    truth_member_count: int
    supported_truth_member_count: int
    exact_variant_member_count: int
    exact_full_genotype: bool | None
    exact_variant_rank: int | None
    expected_gene_count: int
    recovered_gene_count: int
    best_gene_rank: int | None
    unsupported_truth_members: tuple[str, ...]


@dataclass(frozen=True, slots=True)
class ExomiserValidation:
    report_sha256: str
    attempts: int
    comparisons: int
    supported_comparisons: int
    unsupported_comparisons: int
    adapted_reproducible: bool


def load_exomiser_stack(
    path: str | Path, *, repository_root: str | Path
) -> ExomiserStack:
    repository = Path(repository_root).resolve(strict=True)
    manifest_path = Path(path)
    if not manifest_path.is_absolute():
        manifest_path = repository / manifest_path
    try:
        document, content = load_bounded_json(
            manifest_path, reference_root=repository
        )
    except BoundedJSONError as error:
        raise ExomiserError("Exomiser manifest cannot be read") from error
    root = _closed_mapping(
        document,
        {
            "schema",
            "classification",
            "validation_only",
            "tool",
            "runtime",
            "data",
            "analysis",
            "license_review",
            "pheval_contract",
            "supported_comparison",
        },
        "Exomiser manifest",
    )
    if (
        root["schema"] != SCHEMA
        or root["classification"] != "public"
        or root["validation_only"] is not True
    ):
        raise ExomiserError("Exomiser manifest classification is invalid")
    tool = _closed_mapping(
        root["tool"],
        {"name", "version", "release_url", "tag_commit", "archive"},
        "Exomiser tool",
    )
    runtime = _closed_mapping(
        root["runtime"],
        {"image", "image_digest", "java_version", "network", "user"},
        "Exomiser runtime",
    )
    data = _closed_mapping(
        root["data"],
        {
            "release",
            "installed_root",
            "installed_bytes",
            "installed_file_count",
            "installed_file_manifest_path",
            "installed_file_manifest_sha256",
            "archives",
        },
        "Exomiser data",
    )
    analysis = _closed_mapping(
        root["analysis"], {"path", "sha256", "adaptations"}, "Exomiser analysis"
    )
    license_review = _closed_mapping(
        root["license_review"],
        {"code_license", "code_license_url", "data_use_decision", "restrictions"},
        "Exomiser license review",
    )
    pheval = _closed_mapping(
        root["pheval_contract"],
        {
            "version",
            "commit",
            "license",
            "source_url",
            "sort_order",
            "gene_fields",
            "variant_fields",
        },
        "PhEval contract",
    )
    comparison = _closed_mapping(
        root["supported_comparison"],
        {"classes", "unsupported_classes", "metrics", "limitations"},
        "Exomiser comparison",
    )
    image_digest = _require_digest(runtime["image_digest"], "image digest")
    if (
        tool["name"] != "Exomiser"
        or tool["version"] != "15.1.0"
        or tool["tag_commit"] != "98f4e0b6f2f3f5bf790a9f955808c98968559996"
        or runtime["image"] != "docker.io/exomiser/exomiser-cli:15.1.0"
        or IMAGE != f"docker.io/exomiser/exomiser-cli@sha256:{image_digest}"
        or runtime["network"] != "none"
        or runtime["user"] != "65534:65534"
        or data["release"] != "2512"
        or license_review["code_license"] != "AGPL-3.0-or-later"
        or pheval["version"] != "0.7.14"
        or pheval["commit"] != "51d073bedbde5636c1b11c6d381f2d2fce0c74fc"
        or pheval["gene_fields"] != list(GENE_FIELDS)
        or pheval["variant_fields"] != list(VARIANT_FIELDS)
        or pheval["sort_order"] != "descending"
        or set(comparison["classes"]) != {"SNV", "indel"}
        or set(comparison["unsupported_classes"])
        != {"mixed-class", "mosaic", "regulatory", "repeat", "SV/CNV"}
    ):
        raise ExomiserError("Exomiser manifest identity or scope is invalid")
    analysis_path = repository / str(analysis["path"])
    try:
        analysis_content = analysis_path.read_bytes()
    except OSError as error:
        raise ExomiserError("Exomiser analysis configuration cannot be read") from error
    analysis_sha256 = _require_digest(analysis["sha256"], "analysis SHA-256")
    if _sha256(analysis_content) != analysis_sha256:
        raise ExomiserError("Exomiser analysis configuration is stale")
    archives = data["archives"]
    if not isinstance(archives, list) or len(archives) != 3:
        raise ExomiserError("Exomiser archives are incomplete")
    for item in archives:
        entry = _closed_mapping(
            item,
            {
                "role",
                "url",
                "size_bytes",
                "sha256",
                "publisher_md5",
                "object_generation",
            },
            "Exomiser archive",
        )
        if (
            type(entry["url"]) is not str
            or not entry["url"].startswith("https://")
            or type(entry["size_bytes"]) is not int
            or entry["size_bytes"] <= 0
        ):
            raise ExomiserError("Exomiser archive source is invalid")
        _require_digest(entry["sha256"], "archive SHA-256")
    return ExomiserStack(
        manifest_sha256=_sha256(content),
        version=str(tool["version"]),
        data_release=str(data["release"]),
        image=IMAGE,
        analysis_path=str(analysis["path"]),
        analysis_sha256=analysis_sha256,
        installed_root=str(data["installed_root"]),
        installed_file_manifest_sha256=_require_digest(
            data["installed_file_manifest_sha256"], "installed manifest SHA-256"
        ),
        installed_file_count=int(data["installed_file_count"]),
        installed_bytes=int(data["installed_bytes"]),
        code_license=str(license_review["code_license"]),
        data_use_decision=str(license_review["data_use_decision"]),
        pheval_version=str(pheval["version"]),
        pheval_commit=str(pheval["commit"]),
    )


def render_truth_free_phenopacket(
    hpo_text: bytes, *, case_name: str, created: str = "2026-09-12T00:00:00Z"
) -> bytes:
    """Convert a plain HPO list without accepting labels or evaluator metadata."""
    try:
        text = hpo_text.decode("ascii")
    except UnicodeDecodeError as error:
        raise ExomiserError("Exomiser HPO input is not ASCII") from error
    terms = text.splitlines()
    if (
        not terms
        or len(terms) > 1000
        or len(set(terms)) != len(terms)
        or any(HPO_ID.fullmatch(term) is None for term in terms)
        or re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9._-]*", case_name) is None
    ):
        raise ExomiserError("Exomiser HPO input is invalid")
    document = {
        "id": case_name,
        "subject": {"id": "SYNTHETIC", "sex": "UNKNOWN_SEX"},
        "phenotypicFeatures": tuple({"type": {"id": term}} for term in terms),
        "metaData": {
            "created": created,
            "createdBy": "genosapiens-public-validation",
            "resources": (
                {
                    "id": "hp",
                    "name": "Human Phenotype Ontology",
                    "url": "https://hpo.jax.org",
                    "version": "2512",
                    "namespacePrefix": "HP",
                    "iriPrefix": "http://purl.obolibrary.org/obo/HP_",
                },
            ),
            "phenopacketSchemaVersion": "1.0",
        },
    }
    rendered = (
        json.dumps(document, sort_keys=True, separators=(",", ":"), ensure_ascii=True)
        + "\n"
    ).encode("ascii")
    lowered = rendered.lower()
    if any(key.encode("ascii") in lowered for key in FORBIDDEN_RUNNER_KEYS):
        raise ExomiserError("Exomiser runner input crosses the truth boundary")
    return rendered


def exomiser_command(*, output_name: str) -> tuple[str, ...]:
    if re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9._-]*", output_name) is None:
        raise ExomiserError("Exomiser output name is invalid")
    return (
        "analyse",
        "--analysis",
        "/inputs/analysis.yml",
        "--vcf",
        "/inputs/case.vcf",
        "--assembly",
        "GRCh38",
        "--sample",
        "/inputs/sample.json",
        "--output-directory",
        "/results",
        "--output-filename",
        output_name,
        "--output-format",
        "TSV_GENE,TSV_VARIANT",
    )


def run_exomiser(
    stack: ExomiserStack,
    *,
    input_directory: str | Path,
    output_directory: str | Path,
    output_name: str,
    timeout: int = 900,
) -> subprocess.CompletedProcess[bytes]:
    """Run the pinned image offline with immutable inputs and reference data."""
    input_root = Path(input_directory).resolve(strict=True)
    output_root = Path(output_directory).resolve(strict=True)
    if not input_root.is_dir() or not output_root.is_dir() or timeout <= 0:
        raise ExomiserError("Exomiser execution paths or timeout are invalid")
    command = (
        "docker",
        "run",
        "--rm",
        "--pull",
        "never",
        "--network",
        "none",
        "--read-only",
        "--tmpfs",
        "/tmp:rw,nosuid,nodev,noexec,size=1g",
        "--cpus",
        "4",
        "--memory",
        "16g",
        "--pids-limit",
        "1024",
        "--user",
        "65534:65534",
        "-e",
        "EXOMISER_HG38_DATA_VERSION=2512",
        "-e",
        "EXOMISER_PHENOTYPE_DATA_VERSION=2512",
        "-e",
        "JAVA_TOOL_OPTIONS=-Xmx12g",
        "-v",
        f"{stack.installed_root}:/exomiser-data:ro",
        "-v",
        f"{input_root}:/inputs:ro",
        "-v",
        f"{output_root}:/results:rw",
        stack.image,
        *exomiser_command(output_name=output_name),
    )
    try:
        result = subprocess.run(
            command, capture_output=True, check=False, timeout=timeout
        )
    except (OSError, subprocess.TimeoutExpired) as error:
        raise ExomiserError("Exomiser execution failed or timed out") from error
    if result.returncode != 0:
        raise ExomiserError("Exomiser execution returned a non-zero status")
    return result


def _parse_score(value: str) -> float:
    try:
        score = float(value)
    except ValueError as error:
        raise ExomiserError("Exomiser output score is invalid") from error
    if score < 0 or score != score or score in (float("inf"), float("-inf")):
        raise ExomiserError("Exomiser output score is invalid")
    return score


def _tsv_rows(path: str | Path) -> tuple[dict[str, str], ...]:
    try:
        content = Path(path).read_text(encoding="utf-8")
    except (OSError, UnicodeDecodeError) as error:
        raise ExomiserError("Exomiser TSV output cannot be read") from error
    lines = content.splitlines()
    if not lines or not lines[0].startswith("#"):
        raise ExomiserError("Exomiser TSV header is invalid")
    reader = csv.DictReader([lines[0][1:], *lines[1:]], delimiter="\t")
    if reader.fieldnames is None or len(set(reader.fieldnames)) != len(reader.fieldnames):
        raise ExomiserError("Exomiser TSV fields are invalid")
    rows = tuple(dict(item) for item in reader)
    if any(None in row or any(value is None for value in row.values()) for row in rows):
        raise ExomiserError("Exomiser TSV row is malformed")
    return rows


def parse_exomiser_outputs(
    gene_tsv: str | Path, variant_tsv: str | Path
) -> ExomiserOutputs:
    genes = []
    for row in _tsv_rows(gene_tsv):
        required = {
            "GENE_SYMBOL",
            "ENTREZ_GENE_ID",
            "EXOMISER_GENE_COMBINED_SCORE",
            "ID",
        }
        if not required.issubset(row):
            raise ExomiserError("Exomiser gene fields are incomplete")
        genes.append(
            ExomiserGeneRow(
                row["GENE_SYMBOL"],
                row["ENTREZ_GENE_ID"],
                _parse_score(row["EXOMISER_GENE_COMBINED_SCORE"]),
                row["ID"],
            )
        )
    variants = []
    for row in _tsv_rows(variant_tsv):
        required = {
            "CONTIG",
            "START",
            "END",
            "REF",
            "ALT",
            "EXOMISER_GENE_COMBINED_SCORE",
            "ID",
        }
        if not required.issubset(row):
            raise ExomiserError("Exomiser variant fields are incomplete")
        try:
            start, end = int(row["START"]), int(row["END"])
        except ValueError as error:
            raise ExomiserError("Exomiser variant coordinate is invalid") from error
        if start <= 0 or end < start or not row["REF"] or not row["ALT"]:
            raise ExomiserError("Exomiser variant identity is invalid")
        variants.append(
            ExomiserVariantRow(
                row["CONTIG"].removeprefix("chr"),
                start,
                end,
                row["REF"],
                row["ALT"],
                _parse_score(row["EXOMISER_GENE_COMBINED_SCORE"]),
                row["ID"],
            )
        )
    gene_tuple = tuple(genes)
    variant_tuple = tuple(variants)
    if tuple(row.score for row in gene_tuple) != tuple(
        sorted((row.score for row in gene_tuple), reverse=True)
    ):
        raise ExomiserError("Exomiser gene results are not descending")
    return ExomiserOutputs(gene_tuple, variant_tuple)


def pheval_input_as_tsv(outputs: ExomiserOutputs, *, kind: str) -> bytes:
    if kind not in {"gene", "variant"}:
        raise ExomiserError("unknown PhEval result kind")
    fields = GENE_FIELDS if kind == "gene" else VARIANT_FIELDS
    rows: Sequence[object] = outputs.genes if kind == "gene" else outputs.variants
    stream = io.StringIO(newline="")
    writer = csv.DictWriter(stream, fieldnames=fields, delimiter="\t", lineterminator="\n")
    writer.writeheader()
    for row in rows:
        writer.writerow(asdict(row))
    return stream.getvalue().encode("utf-8")


def compare_exomiser_outputs(
    outputs: ExomiserOutputs,
    *,
    truth_members: Sequence[str],
    expected_genes: Sequence[str],
    unsupported_truth_members: Sequence[str] = (),
) -> BaselineComparison:
    """Evaluator-side comparison; never call this before the runner exits."""
    truth = tuple(truth_members)
    unsupported = tuple(unsupported_truth_members)
    if (
        not truth
        or len(set(truth)) != len(truth)
        or any(item not in truth for item in unsupported)
        or len(set(expected_genes)) != len(tuple(expected_genes))
    ):
        raise ExomiserError("Exomiser comparison truth is invalid")
    supported = tuple(item for item in truth if item not in set(unsupported))
    exact_by_group: dict[str, set[str]] = {}
    first_variant_rank: int | None = None
    for rank, row in enumerate(outputs.variants, 1):
        identity = f"{row.chrom}:{row.start}:{row.ref}:{row.alt}"
        exact_by_group.setdefault(row.grouping_id, set()).add(identity)
        if identity in supported and first_variant_rank is None:
            first_variant_rank = rank
    recovered_members = set().union(*exact_by_group.values()) if exact_by_group else set()
    exact_members = len(set(supported).intersection(recovered_members))
    exact_full = (
        None
        if not supported
        else any(set(supported).issubset(members) for members in exact_by_group.values())
    )
    expected = tuple(expected_genes)
    gene_ranks = {
        row.gene_symbol: rank for rank, row in enumerate(outputs.genes, 1)
    }
    recovered_genes = set(expected).intersection(gene_ranks)
    return BaselineComparison(
        truth_member_count=len(truth),
        supported_truth_member_count=len(supported),
        exact_variant_member_count=exact_members,
        exact_full_genotype=exact_full,
        exact_variant_rank=first_variant_rank,
        expected_gene_count=len(expected),
        recovered_gene_count=len(recovered_genes),
        best_gene_rank=(
            min(gene_ranks[gene] for gene in recovered_genes)
            if recovered_genes
            else None
        ),
        unsupported_truth_members=unsupported,
    )


def outputs_sha256(outputs: ExomiserOutputs) -> str:
    return _digest(
        {
            "genes": tuple(asdict(item) for item in outputs.genes),
            "variants": tuple(asdict(item) for item in outputs.variants),
        }
    )


def load_exomiser_validation(
    path: str | Path = CANONICAL_VALIDATION_PATH,
    *,
    repository_root: str | Path,
    expected_sha256: str = CANONICAL_VALIDATION_SHA256,
) -> ExomiserValidation:
    """Authenticate complete development/held-out comparison accounting."""
    repository = Path(repository_root).resolve(strict=True)
    source = Path(path)
    if not source.is_absolute():
        source = repository / source
    try:
        document, content = load_bounded_json(
            source, reference_root=repository, maximum_bytes=128 * 1024
        )
    except BoundedJSONError as error:
        raise ExomiserError("Exomiser validation report cannot be read") from error
    if _sha256(content) != _require_digest(expected_sha256, "validation report SHA-256"):
        raise ExomiserError("Exomiser validation report is stale")
    root = _closed_mapping(
        document,
        {
            "schema", "classification", "validation_only", "controlled_sources_used",
            "report_id", "execution_date", "stack", "provisioning", "attempts",
            "runs", "comparisons", "reproducibility", "scope", "limitations",
        },
        "Exomiser validation",
    )
    if (
        root["schema"] != VALIDATION_SCHEMA
        or root["classification"] != "public"
        or root["validation_only"] is not True
        or root["controlled_sources_used"] is not False
    ):
        raise ExomiserError("Exomiser validation boundary is invalid")
    stack_data = _closed_mapping(
        root["stack"],
        {
            "manifest_path", "manifest_sha256", "analysis_sha256", "image",
            "installed_file_manifest_sha256", "installed_file_count", "installed_bytes",
            "version", "data_release", "pheval_version", "pheval_commit",
        },
        "Exomiser validation stack",
    )
    stack = load_exomiser_stack(
        stack_data["manifest_path"], repository_root=repository
    )
    if (
        stack.manifest_sha256 != stack_data["manifest_sha256"]
        or stack.analysis_sha256 != stack_data["analysis_sha256"]
        or stack.image != stack_data["image"]
        or stack.installed_file_manifest_sha256
        != stack_data["installed_file_manifest_sha256"]
        or stack.installed_file_count != stack_data["installed_file_count"]
        or stack.installed_bytes != stack_data["installed_bytes"]
        or stack.version != stack_data["version"]
        or stack.data_release != stack_data["data_release"]
        or stack.pheval_version != stack_data["pheval_version"]
        or stack.pheval_commit != stack_data["pheval_commit"]
    ):
        raise ExomiserError("Exomiser validation stack is stale")

    provisioning = _closed_mapping(
        root["provisioning"], {"status", "archive_count", "download_bytes"},
        "Exomiser provisioning",
    )
    if provisioning != {
        "status": "succeeded", "archive_count": 3, "download_bytes": 36_861_013_169
    }:
        raise ExomiserError("Exomiser provisioning evidence is incomplete")
    attempts = root["attempts"]
    if (
        not isinstance(attempts, list)
        or len(attempts) != 3
        or [item.get("status") for item in attempts if isinstance(item, Mapping)]
        != ["failed", "succeeded-empty", "succeeded"]
        or attempts[0].get("output_published") is not False
        or any(item.get("output_published") is not True for item in attempts[1:])
    ):
        raise ExomiserError("Exomiser attempt accounting is incomplete")

    runs = root["runs"]
    if not isinstance(runs, list) or len(runs) != 2:
        raise ExomiserError("Exomiser split runs are incomplete")
    expected_outputs = {
        "development": (0, 0, "f4d4c145d2c07331a98a32cb15a60a8df7179dd1e10fe20e2e1580dd0e366a06"),
        "held_out": (2, 2, "9e0638d1544cdb67968db49fcac04a0348f6a75277c299a7b9afa7ae65a17520"),
    }
    for run in runs:
        if not isinstance(run, Mapping) or set(run) != {
            "split", "status", "input_vcf_sha256", "input_records",
            "single_alleles_loaded", "structural_variants_loaded", "genes", "variants",
            "wall_seconds", "cpu_seconds", "peak_memory_bytes", "scratch_bytes",
            "raw_gene_tsv_sha256", "raw_variant_tsv_sha256", "adapted_output_sha256",
            "pheval_gene_sha256", "pheval_variant_sha256",
        }:
            raise ExomiserError("Exomiser split run fields are invalid")
        expected = expected_outputs.get(run["split"])
        if (
            expected is None
            or (run["genes"], run["variants"], run["adapted_output_sha256"])
            != expected
            or run["single_alleles_loaded"] != 17
            or run["status"] not in {"succeeded", "succeeded-empty"}
            or run["wall_seconds"] <= 0
            or run["scratch_bytes"] <= 0
        ):
            raise ExomiserError("Exomiser split run is invalid")
        for key in (
            "input_vcf_sha256", "raw_gene_tsv_sha256", "raw_variant_tsv_sha256",
            "adapted_output_sha256", "pheval_gene_sha256", "pheval_variant_sha256",
        ):
            _require_digest(run[key], f"Exomiser {key}")

    comparisons = root["comparisons"]
    from genosapiens.track1_roster import REQUIRED_SCOPES

    if not isinstance(comparisons, list) or len(comparisons) != 22:
        raise ExomiserError("Exomiser comparison roster is incomplete")
    identities = set()
    status_counts = {"succeeded": 0, "unsupported": 0}
    for value in comparisons:
        item = _closed_mapping(
            value,
            {
                "case_id", "split", "scope", "status", "shared_aggregate_run",
                "truth_member_count", "supported_truth_member_count",
                "exact_variant_member_count", "exact_full_genotype", "exact_variant_rank",
                "expected_gene_count", "recovered_gene_count", "best_gene_rank",
                "unsupported_truth_members",
            },
            "Exomiser comparison",
        )
        identity = (item["split"], item["scope"])
        if (
            identity in identities
            or item["split"] not in {"development", "held_out"}
            or item["scope"] not in REQUIRED_SCOPES
            or item["status"] not in status_counts
            or item["shared_aggregate_run"] is not True
            or item["truth_member_count"]
            != item["supported_truth_member_count"] + len(item["unsupported_truth_members"])
            or item["exact_variant_member_count"] > item["supported_truth_member_count"]
            or item["recovered_gene_count"] > item["expected_gene_count"]
        ):
            raise ExomiserError("Exomiser comparison accounting is invalid")
        identities.add(identity)
        status_counts[item["status"]] += 1
    if identities != {
        (split, scope)
        for split in ("development", "held_out")
        for scope in REQUIRED_SCOPES
    } or status_counts != {"succeeded": 12, "unsupported": 10}:
        raise ExomiserError("Exomiser supported-scope accounting is invalid")

    reproducibility = _closed_mapping(
        root["reproducibility"],
        {
            "development_adapted_outputs_equal", "held_out_adapted_outputs_equal",
            "held_out_raw_outputs_equal", "raw_difference",
        },
        "Exomiser reproducibility",
    )
    if (
        reproducibility["development_adapted_outputs_equal"] is not True
        or reproducibility["held_out_adapted_outputs_equal"] is not True
        or reproducibility["held_out_raw_outputs_equal"] is not False
        or "P-VALUE" not in reproducibility["raw_difference"]
    ):
        raise ExomiserError("Exomiser reproducibility result is invalid")
    scope = _closed_mapping(
        root["scope"], {"supported_common_scopes", "unsupported_scopes", "comparison_gate"},
        "Exomiser scope",
    )
    if set(scope["unsupported_scopes"]) != {
        "mixed_class", "mosaic", "regulatory", "repeat", "sv_cnv"
    } or scope["comparison_gate"] != "report-only; no superiority or non-inferiority release gate":
        raise ExomiserError("Exomiser comparison scope is invalid")
    return ExomiserValidation(
        _sha256(content), len(attempts), len(comparisons),
        status_counts["succeeded"], status_counts["unsupported"], True,
    )
