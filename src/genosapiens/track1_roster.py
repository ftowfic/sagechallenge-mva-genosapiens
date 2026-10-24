"""Complete public Track 1 model-case roster and runner truth firewall."""

from __future__ import annotations

from collections import Counter
from collections.abc import Mapping
from dataclasses import dataclass
import hashlib
import json
from pathlib import Path, PurePosixPath
import re
from typing import Any

from genosapiens.bounded_json import BoundedJSONError, load_bounded_json, read_bounded_bytes
from genosapiens.track1_validation_plan import REQUIRED_SCOPES


SCHEMA = "genosapiens.track1-model-case-roster"
HPO_ID = re.compile(rb"^HP:[0-9]{7}$")
SHA256 = re.compile(r"^[0-9a-f]{64}$")
SPLITS = ("development", "held_out")
TARGET_KINDS = {"singleton", "full_pair", "locus"}
SUBMISSION_PROFILE_PATH = "config/track1/submission-acceptance-profile.json"
SUBMISSION_PROFILE_SHA256 = "5200cd15d33b813f9b21ed42acfc8bbb0df59c684ec2daf07062e49b4b86f7e4"
PROFILE_ACCEPTANCE_SHA256 = "1e9240a1d642364b5fa5eeefe64a65442039b381df0014154c7950117732b7e7"
STATUSES = ("succeeded", "failed", "timed_out", "unsupported", "truncated")
REQUIRED_MODEL_COVERAGE = (
    "clinical-unasserted-exploratory-pair",
    "completed-population-lookup-unknown-rescue",
    "independent-non-population-contradiction",
    "population-excluded-unresolved-partner",
    "unknown-count-hypothesis",
)
PUBLIC_BACKGROUND = {
    "1:602113:T:TGCCCA",
    "1:604358:G:C",
    "1:604360:T:C",
    "1:631859:CG:C",
    "1:1735008:G:GGTTT",
    "1:1735008:G:GTTTT",
}
POPULATION_GENERATOR_CASES = {
    "development-population-model-state-01": (
        "development", 101, (
            "clinical-unasserted-exploratory-pair",
            "completed-population-lookup-unknown-rescue",
            "unknown-count-hypothesis",
        ),
    ),
    "development-population-model-state-02": (
        "development", 101, (
            "completed-population-lookup-unknown-rescue",
            "population-excluded-unresolved-partner",
        ),
    ),
    "development-population-model-state-03": (
        "development", 101, ("independent-non-population-contradiction",),
    ),
    "held-out-population-model-state-01": (
        "held_out", 303, (
            "clinical-unasserted-exploratory-pair",
            "completed-population-lookup-unknown-rescue",
            "unknown-count-hypothesis",
        ),
    ),
}


class Track1RosterError(RuntimeError):
    """Raised without disclosing evaluator-only truth content."""


def _sha256(content: bytes) -> str:
    return hashlib.sha256(content).hexdigest()


def _closed(value: object, fields: set[str], label: str) -> Mapping[str, Any]:
    if not isinstance(value, Mapping) or set(value) != fields:
        raise Track1RosterError(f"{label} fields are invalid")
    return value


def _digest(value: object, label: str) -> str:
    if type(value) is not str or SHA256.fullmatch(value) is None:
        raise Track1RosterError(f"{label} is not a SHA-256")
    return value


def _safe_path(root: Path, value: object) -> Path:
    if type(value) is not str:
        raise Track1RosterError("roster artifact path is invalid")
    relative = PurePosixPath(value)
    if relative.is_absolute() or relative == PurePosixPath(".") or ".." in relative.parts:
        raise Track1RosterError("roster artifact path is unsafe")
    path = root.joinpath(*relative.parts)
    if any(component.is_symlink() for component in (path, *path.parents)):
        raise Track1RosterError("roster artifact cannot be a symbolic link")
    try:
        resolved = path.resolve(strict=True)
    except OSError as error:
        raise Track1RosterError("roster artifact cannot be resolved") from error
    if not resolved.is_relative_to(root) or not resolved.is_file():
        raise Track1RosterError("roster artifact escapes the repository")
    return resolved


def _artifact(root: Path, value: object, *, maximum: int = 2 * 1024 * 1024) -> bytes:
    entry = _closed(value, {"role", "path", "size_bytes", "sha256"}, "runner input")
    path = _safe_path(root, entry["path"])
    try:
        content = read_bounded_bytes(
            path,
            reference_root=root,
            maximum_bytes=maximum,
            expected_size_bytes=entry["size_bytes"],
        )
    except BoundedJSONError as error:
        raise Track1RosterError("roster input binding is stale") from error
    if _sha256(content) != _digest(entry["sha256"], "runner input SHA-256"):
        raise Track1RosterError("roster input binding is stale")
    return content


def _fixture_binding_bytes(repository: Path, binding: dict) -> bytes:
    if (type(binding) is not dict or not {"path", "sha256"} <= binding.keys()
            or not binding.keys() <= {"path", "sha256", "size_bytes", "audience", "role"}
            or "audience" in binding and binding["audience"] != "evaluator_only"
            or "size_bytes" in binding and (
                type(binding["size_bytes"]) is not int or binding["size_bytes"] < 0)):
        raise Track1RosterError("fixture qualification binding is invalid")
    content = read_bounded_bytes(
        _safe_path(repository, binding["path"]), reference_root=repository,
        maximum_bytes=8 * 1024 * 1024,
        expected_size_bytes=binding.get("size_bytes"),
    )
    if _sha256(content) != _digest(binding["sha256"], "fixture qualification SHA-256"):
        raise Track1RosterError("fixture qualification binding is stale")
    return content


def _retained_fixture_closures(repository, original, expanded, old_truth, new_truth,
                               old_annotations, new_annotations):
    if (len(original["cases"]) != 30 or len(expanded["cases"]) != 34
            or expanded["cases"][:30] != original["cases"]
            or {c["case_id"] for c in expanded["cases"][30:]} != set(POPULATION_GENERATOR_CASES)
            or new_truth["cases"][:len(old_truth["cases"])] != old_truth["cases"]
            or new_annotations["expectations"][:len(old_annotations["expectations"])]
               != old_annotations["expectations"]):
        raise Track1RosterError("retained fixture content differs")
    truth_by_case = {row["case_id"]: row for row in old_truth["cases"]}
    proofs = []
    for case in original["cases"]:
        for binding in case["runner_inputs"]:
            _fixture_binding_bytes(repository, binding)
        formats = []
        for relative in case["source_format_bindings"]:
            content = read_bounded_bytes(_safe_path(repository, relative), reference_root=repository,
                                         maximum_bytes=64 * 1024)
            formats.append({"path": relative, "size_bytes": len(content), "sha256": _sha256(content)})
        proofs.append({"case_id": case["case_id"], "case_sha256": _sha256(_case_bytes(case)),
            "runner_inputs": case["runner_inputs"], "source_formats": formats,
            "truth_sha256": _sha256(_case_bytes(truth_by_case[case["case_id"]])),
            "annotation_sha256": _sha256(_case_bytes(case["annotation_expectations"]))})
    return proofs


def _authenticate_population_fixture_qualification(repository, profile, singleton):
    """Verify generating provenance and retained closures for the current campaign."""
    experiment = "track1-submission-acceptance-population-models-20260917"
    report_path = "artifacts/validation/track1-population-fixture-qualification-2026-09-17.json"
    bridge_path = "artifacts/validation/track1-retained-fixture-qualification-2026-09-17-population-models.json"
    manifest_path = "config/track1/experiments/population-models-20260917-generating-source/manifest.json"
    index = {row["path"]: row for row in profile["configuration_bindings"]}
    if len(index) != len(profile["configuration_bindings"]):
        raise Track1RosterError("fixture qualification bindings are duplicated")
    report = json.loads(_fixture_binding_bytes(repository, index[report_path]))
    bridge = json.loads(_fixture_binding_bytes(repository, index[bridge_path]))
    manifest = json.loads(_fixture_binding_bytes(repository, profile["generating_source_manifest"]))
    schedule = [
        {"stage": "A", "split": "development", "case_ids": [
            "development-population-model-state-01", "development-population-model-state-02",
            "development-population-model-state-03", "development-source-backed-singleton"]},
        {"stage": "B", "split": "development", "case_ids": ["development-mixed-class", "development-repeat"]},
        {"stage": "C", "split": "held_out", "case_ids": [
            "held-out-population-model-state-01", "held-out-source-backed-singleton"]},
    ]
    names = [name for stage in schedule for name in stage["case_ids"]]
    positives = {"development-population-model-state-01", "development-source-backed-singleton",
                 "held-out-population-model-state-01", "held-out-source-backed-singleton"}
    if (profile["stage_schedule"] != schedule or profile["case_budget"] != 8
            or [row["case_id"] for row in profile["cases"]] != names
            or profile["declared_case_ids"] != sorted(names)
            or any(type(row["positive_exact_export"]) is not bool
                   or row["positive_exact_export"] is not (row["case_id"] in positives)
                   for row in profile["cases"])):
        raise Track1RosterError("population campaign matrix differs")
    if (profile["positive_export_contract"] != {
            "optional_products_required": False,
            "eligible_csv_contains_exact_genotype": True,
            "eligible_rank_at_most": 10,
            "exact_expected_genotype": True,
            "mandatory_byte_identical_replay": True,
            "required_for_each_positive_case": True,
        } or report["truth"]["sha256"] != profile["source_causal_truth_sha256"]
           or report["annotations"]["sha256"] != profile["source_annotation_sha256"]):
        raise Track1RosterError("population positive export contract differs")
    if (profile["generating_source_manifest"]["path"] != manifest_path
            or profile["generating_source_manifest"]["sha256"] != index[manifest_path]["sha256"]
            or report["schema"] != "genosapiens.track1-population-fixture-qualification"
            or report["classification"] != "public" or report["experiment_id"] != experiment
            or report["status"] != "source-and-annotation-qualified"
            or report["controlled_use_authorized"] is not False
            or report["ranking_outcomes_used"] is not False
            or report["positive_export_acceptance_established"] is not False
            or report["target_annotation_count"] != 8
            or set(report["new_case_ids"]) != set(POPULATION_GENERATOR_CASES)
            or report["source_roster"] != profile["source_roster"]
            or bridge["schema"] != "genosapiens.track1-retained-fixture-qualification"
            or bridge["classification"] != "public" or bridge["experiment_id"] != experiment
            or bridge["status"] != "exact-retention-qualified"
            or bridge["controlled_use_authorized"] is not False
            or bridge["source_roster"] != report["source_roster"]
            or profile["retained_source_roster"] != bridge["original_source_roster"]
            or singleton["source_roster"] != bridge["original_source_roster"]
            or bridge["original_source_roster"]["sha256"]
               != "1fb8e482658b1135931be20a9e6c9ba9727d55ab2ad478fdf370583ebdbf6466"
            or set(bridge["retained_case_ids"]) != {
                "development-source-backed-singleton", "held-out-source-backed-singleton",
                "development-mixed-class", "development-repeat"}):
        raise Track1RosterError("population fixture qualification is invalid")
    singleton_path = "artifacts/validation/track1-singleton-fixture-qualification-2026-09-17.json"
    if (bridge["original_singleton_qualification"]["path"] != singleton_path
            or bridge["original_singleton_qualification"]["sha256"] != index[singleton_path]["sha256"]):
        raise Track1RosterError("retained singleton qualification is invalid")
    _fixture_binding_bytes(repository, bridge["original_singleton_qualification"])
    old = json.loads(_fixture_binding_bytes(repository, bridge["original_source_roster"]))
    new = json.loads(_fixture_binding_bytes(repository, report["source_roster"]))
    old_truth = json.loads(_fixture_binding_bytes(repository, old["causal_truth"]))
    old_annotations = json.loads(_fixture_binding_bytes(repository, old["annotation_expectations"]))
    new_truth = json.loads(_fixture_binding_bytes(repository, new["causal_truth"]))
    new_annotations = json.loads(_fixture_binding_bytes(repository, new["annotation_expectations"]))
    if any(report[name]["path"] != new[field]["path"] or report[name]["sha256"] != new[field]["sha256"]
           for name, field in (("truth", "causal_truth"), ("annotations", "annotation_expectations"))):
        raise Track1RosterError("population aggregate qualification differs")
    if bridge["all_original_case_closures"] != _retained_fixture_closures(
            repository, old, new, old_truth, new_truth, old_annotations, new_annotations):
        raise Track1RosterError("retained fixture closure differs")
    freeze = json.loads(_fixture_binding_bytes(repository, report["development_numerical_freeze"]))
    selection = json.loads(_fixture_binding_bytes(repository, freeze["development_selection"]))
    qualification = json.loads(_fixture_binding_bytes(repository, freeze["development_qualification"]))
    search = json.loads(_fixture_binding_bytes(repository, selection["source_search"]))
    held = json.loads(_fixture_binding_bytes(repository, report["held_out_selection"]))
    if (freeze["schema"] != "genosapiens.track1-population-development-numerical-freeze"
            or freeze["experiment_id"] != experiment or freeze["status"] != "frozen"
            or freeze["held_out_sources_inspected"] is not False
            or qualification["status"] != "qualified"
            or qualification["selected"]["subset"] != "genomes"
            or selection["selector"] != {"subset": "genomes", "stratum": "adjusted", "ancestry": None}
            or held["development_numerical_freeze"] != report["development_numerical_freeze"]
            or held["development_settings_changed"] is not False
            or held["ranking_outcomes_used"] is not False):
        raise Track1RosterError("population numerical freeze differs")
    _authenticate_complete_enumeration_transition(repository, profile, freeze["policy"])
    _fixture_binding_bytes(repository, freeze["derivation"])
    if (manifest["schema"] != "genosapiens.track1-population-generating-source"
            or manifest["classification"] != "public"
            or manifest["scope"] != "immutable-generating-code-provenance"
            or manifest["controlled_use_authorized"] is not False
            or manifest["fixture_qualification"]["path"] != report_path
            or manifest["fixture_qualification"]["sha256"] != index[report_path]["sha256"]
            or manifest["source_search"] != selection["source_search"]):
        raise Track1RosterError("generating source provenance is invalid")
    dependencies = {}
    for binding in (*report["producer_dependencies"], *search["producer_dependencies"]):
        if binding["path"] in dependencies and dependencies[binding["path"]] != binding:
            raise Track1RosterError("generating source dependencies conflict")
        dependencies[binding["path"]] = binding
    if manifest["dependencies"] != [dependencies[path] for path in sorted(dependencies)]:
        raise Track1RosterError("generating source closure differs")
    archive = Path(manifest_path).parent.as_posix()
    for binding in manifest["dependencies"]:
        _fixture_binding_bytes(repository, {**binding, "path": archive + "/" + binding["path"]})


def _plain_hpo(content: bytes) -> None:
    lines = content.splitlines()
    if (
        not lines
        or not content.endswith(b"\n")
        or len(lines) != len(set(lines))
        or any(HPO_ID.fullmatch(line) is None for line in lines)
    ):
        raise Track1RosterError("roster HPO input is invalid")


def _vcf_identities(content: bytes) -> tuple[set[str], str]:
    identities: set[str] = set()
    sample: str | None = None
    try:
        lines = content.decode("ascii").splitlines()
    except UnicodeDecodeError as error:
        raise Track1RosterError("roster VCF is not ASCII") from error
    for line in lines:
        if line.startswith("#CHROM\t"):
            fields = line.split("\t")
            if len(fields) != 10:
                raise Track1RosterError("roster VCF sample header is invalid")
            sample = fields[9]
        elif not line.startswith("#"):
            fields = line.split("\t")
            if len(fields) != 10:
                raise Track1RosterError("roster VCF record is invalid")
            for alt in fields[4].split(","):
                identities.add(":".join((fields[0], fields[1], fields[3], alt)))
    if sample is None or not identities:
        raise Track1RosterError("roster VCF is empty")
    return identities, sample


def _repeat_member_present(member: object, contents: list[bytes]) -> bool:
    if type(member) is not str:
        return False
    encoded = member.encode("ascii")
    if any(encoded in content for content in contents):
        return True
    match = re.fullmatch(r"longtr-region:([^:]+):([0-9]+)-([0-9]+)", member)
    if match is None:
        return False
    chrom, start, end = match.groups()
    signature = f"{chrom}\t{start}\t".encode("ascii")
    end_field = f"END={end}".encode("ascii")
    return any(
        any(line.startswith(signature) and end_field in line for line in content.splitlines())
        for content in contents
    )


@dataclass(frozen=True, slots=True)
class RosterRunnerInput:
    role: str
    sha256: str
    size_bytes: int


@dataclass(frozen=True, slots=True)
class RosterRunnerView:
    inputs: tuple[RosterRunnerInput, ...]
    source_format_sha256s: tuple[str, ...]
    background_sha256: str

    def as_bytes(self) -> bytes:
        return (
            json.dumps(
                {
                    "inputs": tuple(
                        {"role": item.role, "sha256": item.sha256, "size_bytes": item.size_bytes}
                        for item in self.inputs
                    ),
                    "source_format_sha256s": self.source_format_sha256s,
                    "background_sha256": self.background_sha256,
                },
                sort_keys=True,
                separators=(",", ":"),
            )
            + "\n"
        ).encode("ascii")


@dataclass(frozen=True, slots=True)
class RosterCase:
    case_id: str
    split: str
    scope: str
    target_kind: str
    truth_members: tuple[str, ...]
    expected_gene_symbols: tuple[str, ...]
    runner_view: RosterRunnerView
    required_model_coverage: tuple[str, ...] = ()


@dataclass(frozen=True, slots=True)
class Track1Roster:
    roster_sha256: str
    annotation_sha256: str
    causal_truth_sha256: str
    cases: tuple[RosterCase, ...]
    shared_background: bool

    def runner_case(self, case_id: str, *, allow_held_out: bool = False) -> RosterRunnerView:
        matches = [case for case in self.cases if case.case_id == case_id]
        if len(matches) != 1 or matches[0].split == "held_out" and not allow_held_out:
            raise Track1RosterError("roster runner view is unavailable")
        view = matches[0].runner_view
        lowered = view.as_bytes().lower()
        if any(word in lowered for word in (b"causal", b"oracle", b"recipe", b"truth")):
            raise Track1RosterError("roster runner view crosses the truth boundary")
        return view


def load_track1_roster(path: str | Path, *, repository_root: str | Path) -> Track1Roster:
    repository = Path(repository_root).resolve(strict=True)
    roster_path = Path(path)
    if not roster_path.is_absolute():
        roster_path = repository / roster_path
    try:
        document, content = load_bounded_json(
            roster_path, reference_root=repository, maximum_bytes=512 * 1024
        )
    except BoundedJSONError as error:
        raise Track1RosterError("roster manifest cannot be read") from error
    root = _closed(
        document,
        {
            "schema",
            "classification",
            "validation_only",
            "controlled_use_authorized",
            "roster_id",
            "required_scopes",
            "required_model_coverage",
            "splits",
            "background",
            "annotation_expectations",
            "causal_truth",
            "cases",
            "metric_contract",
            "limitations",
        },
        "roster",
    )
    if (
        root["schema"] != SCHEMA
        or root["classification"] != "public"
        or root["validation_only"] is not True
        or root["controlled_use_authorized"] is not False
        or not root["required_scopes"]
        or root["required_scopes"] != sorted(set(root["required_scopes"]))
        or not set(root["required_scopes"]) <= set(REQUIRED_SCOPES)
        or root["required_model_coverage"] != sorted(set(root["required_model_coverage"]))
        or not set(root["required_model_coverage"]) <= set(REQUIRED_MODEL_COVERAGE)
        or root["splits"] != list(SPLITS)
    ):
        raise Track1RosterError("roster identity is invalid")
    background = _closed(
        root["background"],
        {
            "background_id",
            "subject_id",
            "shared_across_cases",
            "source_manifest_path",
            "disclosure",
        },
        "roster background",
    )
    background_path = _safe_path(repository, background["source_manifest_path"])
    background_sha256 = _sha256(background_path.read_bytes())
    if (
        background["background_id"] != "hg002-public-background"
        or background["subject_id"] != "HG002"
        or background["shared_across_cases"] is not True
        or "not independent patients" not in str(background["disclosure"])
    ):
        raise Track1RosterError("roster background is invalid")
    annotation_binding = _closed(
        root["annotation_expectations"], {"path", "size_bytes", "sha256"}, "annotations"
    )
    truth_binding = _closed(
        root["causal_truth"],
        {"path", "size_bytes", "sha256", "audience"},
        "causal truth",
    )
    if truth_binding["audience"] != "evaluator_only":
        raise Track1RosterError("roster truth audience is invalid")
    annotation_path = _safe_path(repository, annotation_binding["path"])
    truth_path = _safe_path(repository, truth_binding["path"])
    try:
        annotations, annotation_content = load_bounded_json(
            annotation_path,
            reference_root=repository,
            maximum_bytes=256 * 1024,
            expected_size_bytes=annotation_binding["size_bytes"],
        )
        truth, truth_content = load_bounded_json(
            truth_path,
            reference_root=repository,
            maximum_bytes=256 * 1024,
            expected_size_bytes=truth_binding["size_bytes"],
        )
    except BoundedJSONError as error:
        raise Track1RosterError("roster evaluator artifact is stale") from error
    if (
        _sha256(annotation_content)
        != _digest(annotation_binding["sha256"], "annotation SHA-256")
        or _sha256(truth_content)
        != _digest(truth_binding["sha256"], "truth SHA-256")
    ):
        raise Track1RosterError("roster evaluator artifact is stale")
    if (
        annotations.get("schema") != "genosapiens.track1-roster-annotation-expectations"
        or annotations.get("classification") != "public"
        or annotations.get("validation_only") is not True
        or truth.get("schema") != "genosapiens.track1-roster-causal-truth"
        or truth.get("audience") != "evaluator_only"
    ):
        raise Track1RosterError("roster evaluator artifact identity is invalid")
    tool = annotations.get("tool")
    if (
        not isinstance(tool, Mapping)
        or tool.get("version") != "116.2"
        or tool.get("cache_release") != "116"
        or tool.get("assembly") != "GRCh38"
        or tool.get("image_digest")
        != "f22d1f89b28c6c4caf043ddb2269fdfff8940364da589cadbbabc6212a436949"
    ):
        raise Track1RosterError("roster annotation release is invalid")
    truth_cases = truth.get("cases")
    cases = root["cases"]
    if not isinstance(cases, list) or not isinstance(truth_cases, list):
        raise Track1RosterError("roster cases are invalid")
    truth_by_id = {
        item.get("truth_id"): item for item in truth_cases if isinstance(item, Mapping)
    }
    if len(cases) < 2 * len(root["required_scopes"]) or len(truth_by_id) != len(cases):
        raise Track1RosterError("roster cardinality is invalid")
    parsed: list[RosterCase] = []
    seen_case_ids: set[str] = set()
    split_scope = Counter()
    profiles = set()
    coverage_by_split = {split: set() for split in SPLITS}
    group_splits = {}
    for case_value in cases:
        case = _closed(
            case_value,
            {
                "case_id",
                "split",
                "scope",
                "inheritance_context",
                "group_id",
                "generator",
                "background_id",
                "shared_background_disclosure",
                "runner_inputs",
                "source_format_bindings",
                "annotation_expectations",
                "truth_id",
                "target_kind",
                "metric_applicability",
                "regression_context",
                "required_model_coverage",
            },
            "roster case",
        )
        case_id, split, scope = case["case_id"], case["split"], case["scope"]
        if (
            type(case_id) is not str
            or case_id in seen_case_ids
            or split not in SPLITS
            or scope not in root["required_scopes"]
            or case["inheritance_context"] != scope
            or case["background_id"] != "hg002-public-background"
            or case["target_kind"] not in TARGET_KINDS
            or "not independent patients" not in case["shared_background_disclosure"]
        ):
            raise Track1RosterError("roster case identity is invalid")
        generator = _closed(case["generator"], {"id", "version", "seed"}, "generator")
        if (
            generator["id"] != "track1-model-roster-generator"
            or (generator["version"], generator["seed"]) not in (
                (("1", 101), ("2", 101), ("3", 101), ("4", 101)) if split == "development"
                else (("1", 202), ("2", 303), ("3", 303), ("4", 303))
            )
        ):
            raise Track1RosterError("roster generator binding is invalid")
        if generator["version"] == "3" and (
            scope != "dominant"
            or case_id != split.replace("_", "-") + "-source-backed-singleton"
        ):
            raise Track1RosterError("source-backed singleton generator scope is invalid")
        required_coverage = case["required_model_coverage"]
        if not isinstance(required_coverage, list) or required_coverage != sorted(set(required_coverage)) or not set(required_coverage) <= set(REQUIRED_MODEL_COVERAGE):
            raise Track1RosterError("roster model coverage declaration is invalid")
        if generator["version"] == "4" and (
            POPULATION_GENERATOR_CASES.get(case_id)
            != (split, generator["seed"], tuple(required_coverage))
            or scope != "compound_heterozygous"
            or case["target_kind"] != "full_pair"
        ):
            raise Track1RosterError("population generator case contract is invalid")
        coverage_by_split[split].update(required_coverage)
        group = case["group_id"]
        if type(group) is not str or not group or group in group_splits and group_splits[group] != split:
            raise Track1RosterError("related roster cases cross splits")
        group_splits[group] = split
        inputs = case["runner_inputs"]
        if not isinstance(inputs, list) or len(inputs) not in {2, 3}:
            raise Track1RosterError("roster runner inputs are incomplete")
        runner_inputs = []
        vcf_identities: set[str] = set()
        vcf_contents: list[bytes] = []
        vcf_roles_by_sha256: dict[str, str] = {}
        sample_names = set()
        roles = set()
        for input_value in inputs:
            input_entry = _closed(
                input_value, {"role", "path", "size_bytes", "sha256"}, "runner input"
            )
            role = input_entry["role"]
            if role in roles or role not in {"hpo_text", "small_variant_vcf", "repeat_native_vcf"}:
                raise Track1RosterError("roster runner input role is invalid")
            input_content = _artifact(repository, input_entry)
            if role == "hpo_text":
                _plain_hpo(input_content)
            else:
                identities, sample = _vcf_identities(input_content)
                vcf_identities.update(identities)
                vcf_contents.append(input_content)
                sample_names.add(sample)
                vcf_roles_by_sha256[input_entry["sha256"]] = str(role)
            roles.add(role)
            runner_inputs.append(
                RosterRunnerInput(
                    str(role), _digest(input_entry["sha256"], "input SHA-256"), input_entry["size_bytes"]
                )
            )
        if {"hpo_text", "small_variant_vcf"} - roles or sample_names != {"SYNTHETIC"}:
            raise Track1RosterError("roster challenge-format inputs are invalid")
        if not PUBLIC_BACKGROUND.issubset(vcf_identities):
            raise Track1RosterError("roster public background is incomplete")
        formats = case["source_format_bindings"]
        if (
            not isinstance(formats, list)
            or not formats
            or len(formats) != len(set(formats))
        ):
            raise Track1RosterError("roster source formats are incomplete")
        format_sha256s = []
        format_input_sha256s = set()
        for item in formats:
            format_path = _safe_path(repository, item)
            try:
                source_format, source_format_content = load_bounded_json(
                    format_path, reference_root=repository, maximum_bytes=64 * 1024
                )
            except BoundedJSONError as error:
                raise Track1RosterError("roster source format is invalid") from error
            source_format = _closed(
                source_format,
                {
                    "schema", "input_sha256", "reference_role",
                    "reference_identity_id", "reference_manifest_sha256",
                    "mitochondrial_contigs", "repeat_identity",
                },
                "source format",
            )
            input_sha256 = _digest(
                source_format["input_sha256"], "source format input SHA-256"
            )
            input_role = vcf_roles_by_sha256.get(input_sha256)
            expected_repeat = (
                "longtr-fd73e0d-vcf4.1/v1"
                if input_role == "repeat_native_vcf"
                else "repid-ru-repcn-pos-end"
            )
            if (
                input_role is None
                or input_sha256 in format_input_sha256s
                or source_format["schema"] != "genosapiens.supplied-vcf-format"
                or source_format["reference_role"] != "public_annotation_benchmark"
                or source_format["reference_identity_id"]
                != "grch38-no-alt-hs38d1-public-native-v1"
                or source_format["reference_manifest_sha256"]
                != "27eaf00fde35e5fa83f945c63e5dee871960374d90488cc0797b5ebea2fb1c25"
                or source_format["mitochondrial_contigs"]
                != ([] if input_role == "repeat_native_vcf" else ["M"])
                or source_format["repeat_identity"] != expected_repeat
            ):
                raise Track1RosterError("roster source format binding is invalid")
            format_input_sha256s.add(input_sha256)
            format_sha256s.append(_sha256(source_format_content))
        if format_input_sha256s != set(vcf_roles_by_sha256):
            raise Track1RosterError("roster source format coverage is incomplete")
        truth_case = truth_by_id.get(case["truth_id"])
        if not isinstance(truth_case, Mapping):
            raise Track1RosterError("roster truth binding is invalid")
        members = truth_case.get("members")
        expected_count = 2 if scope in {"compound_heterozygous", "mixed_class"} else 1
        if (
            truth_case.get("case_id") != case_id
            or truth_case.get("audience") != "evaluator_only"
            or truth_case.get("scope") != scope
            or truth_case.get("target_kind") != case["target_kind"]
            or not isinstance(members, list)
            or len(members) != expected_count
            or len(set(members)) != expected_count
            or scope != "repeat" and not set(members).issubset(vcf_identities)
        ):
            raise Track1RosterError("roster target genotype is invalid")
        if scope == "repeat" and not all(
            _repeat_member_present(member, vcf_contents) for member in members
        ):
            raise Track1RosterError("roster repeat target is absent from runner input")
        expectations = case["annotation_expectations"]
        if not isinstance(expectations, list) or len(expectations) != expected_count:
            raise Track1RosterError("roster annotation expectations are incomplete")
        expected_genes = sorted(
            {
                gene
                for expectation in expectations
                if isinstance(expectation, Mapping)
                for gene in expectation.get("gene_symbols", ())
            }
        )
        if {item.get("member") for item in expectations if isinstance(item, Mapping)} != set(members):
            raise Track1RosterError("roster annotation expectation target is invalid")
        applicability = _closed(
            case["metric_applicability"],
            {"gene_rank", "exact_variant_rank", "exact_full_pair", "repeat_locus_rank", "challenge_evaluator"},
            "metric applicability",
        )
        if (
            applicability["gene_rank"] is not bool(expected_genes)
            or applicability["exact_variant_rank"] is not (scope != "repeat")
            or applicability["exact_full_pair"] is not (expected_count == 2)
            or applicability["repeat_locus_rank"] is not (scope == "repeat")
            or applicability["challenge_evaluator"] is not (scope != "repeat")
        ):
            raise Track1RosterError("roster metric applicability is invalid")
        regression = _closed(
            case["regression_context"], {"profile", "identifiers", "mechanisms", "observations"}, "regression context"
        )
        profiles.add(regression["profile"])
        parsed.append(
            RosterCase(
                case_id,
                split,
                scope,
                str(case["target_kind"]),
                tuple(members),
                tuple(expected_genes),
                RosterRunnerView(
                    tuple(runner_inputs), tuple(format_sha256s), background_sha256
                ),
                tuple(required_coverage),
            )
        )
        seen_case_ids.add(case_id)
        split_scope[(split, scope)] += 1
    if (
        set(split_scope) != {(split, scope) for split in SPLITS for scope in root["required_scopes"]}
        or any(not set(root["required_model_coverage"]) <= value for value in coverage_by_split.values())
        or not {
            "plain_challenge_context",
            "missing_optional_identifiers",
            "conflicting_mechanisms",
            "conflicting_observations",
        } <= profiles
    ):
        raise Track1RosterError("roster split, scope, or regression coverage is incomplete")
    metrics = root["metric_contract"]
    if (
        not isinstance(metrics, Mapping)
        or metrics.get("statuses") != list(STATUSES)
        or metrics.get("non_success_accounting") != "retained-in-full-roster-denominator"
    ):
        raise Track1RosterError("roster metric contract is invalid")
    return Track1Roster(
        _sha256(content),
        _sha256(annotation_content),
        _sha256(truth_content),
        tuple(parsed),
        True,
    )


def _authenticate_complete_enumeration_transition(repository, profile, frozen_policy):
    """Authenticate the one approved execution-cap change against original bytes.

    The original population numerical freeze is immutable. This current-profile
    transition preserves its complete scientific policy, derivation and source
    selection while requiring uncapped enumeration in the bounded executor.
    """
    transition = profile["complete_enumeration_transition"]
    if (set(transition) != {"schema", "version", "decision", "historical_policy", "changed_field",
                           "previous_value", "current_value"}
            or transition["schema"] != "genosapiens.track1-complete-enumeration-policy-transition"
            or type(transition["version"]) is not int or transition["version"] != 1
            or transition["decision"] != "repair-track1-bounded-wgs-execution:D28"
            or transition["changed_field"] != "fixed_settings.maximum_pairs_per_context"
            or type(transition["previous_value"]) is not int or transition["previous_value"] != 10000
            or transition["current_value"] is not None
            or frozen_policy["path"] != profile["policy_path"]
            or frozen_policy["sha256"] != transition["historical_policy"]["sha256"]
            or frozen_policy["size_bytes"] != transition["historical_policy"]["size_bytes"]):
        raise Track1RosterError("complete enumeration transition differs")
    original = json.loads(_fixture_binding_bytes(repository, transition["historical_policy"]))
    current_bytes = read_bounded_bytes(repository / profile["policy_path"], reference_root=repository)
    current = json.loads(current_bytes)
    if (_sha256(current_bytes) != profile["policy_sha256"]
            or original["fixed_settings"]["maximum_pairs_per_context"] != 10000
            or current["fixed_settings"]["maximum_pairs_per_context"] is not None):
        raise Track1RosterError("complete enumeration policy binding differs")
    original["fixed_settings"]["maximum_pairs_per_context"] = None
    if original != current:
        raise Track1RosterError("complete enumeration changed another scientific policy field")


def profile_acceptance_digest(profile):
    """Bind the evaluator-only staged declaration independently of its roster."""
    return _sha256(_case_bytes(profile['analysis_profile_acceptance']))


def require_profile_acceptance_binding(profile, plan, *, analysis_profile):
    """Reject foreign scope/plan bindings before acquiring case resources."""
    if (analysis_profile not in {'focused', 'exhaustive'}
            or plan.analysis_profile != analysis_profile
            or plan.profile_acceptance_sha256 != PROFILE_ACCEPTANCE_SHA256
            or profile_acceptance_digest(profile) != PROFILE_ACCEPTANCE_SHA256
            or plan.submission_acceptance_profile_sha256 != SUBMISSION_PROFILE_SHA256):
        raise Track1RosterError('analysis profile acceptance binding is stale or invalid')
    if analysis_profile == 'exhaustive':
        raise Track1RosterError('profile_not_qualified')
    return profile['analysis_profile_acceptance']


def _validate_profile_acceptance(repository, profile):
    declaration = profile['analysis_profile_acceptance']
    if (profile_acceptance_digest(profile) != PROFILE_ACCEPTANCE_SHA256
            or declaration['schema'] != 'genosapiens.track1-analysis-profile-acceptance'
            or declaration['version'] != 1
            or declaration['phase_order'] != ['M1', 'M2']
            or declaration['default_profile'] != 'focused'
            or declaration['total_runtime_gate_seconds'] is not None
            or declaration['minimum_speedup'] is not None
            or declaration['optional_products_required'] is not False
            or declaration['independent_clinical_validation_claimed'] is not False):
        raise Track1RosterError('analysis profile acceptance declaration is invalid')
    cases = {row['case_id']: row for row in profile['cases']}
    if len(cases) != len(declaration['access_history']):
        raise Track1RosterError('analysis profile access history is incomplete')
    for access in declaration['access_history']:
        case = cases[access['case_id']]
        if (access['original_split'] != case['split']
                or access['runner_inputs'] != case['runner_inputs']
                or access['generator'] != case['generator']
                or access['group_id'] != case['group_id']
                or access['current_use'] != 'development_regression'
                or access['untouched_holdout_claimed'] is not False):
            raise Track1RosterError('analysis profile access history changed')
    reuse = declaration['verification_policy']['completed_input_reuse']
    if reuse != dict(authority='existing_bound_producer_completion',
            startup='bounded_metadata_bindings_and_file_stability',
            requested_records='strict_type_identity_hash_and_proof_validation',
            domain_closure='during_required_consumption', unread_payloads='producer_validation_and_retained_storage',
            full_parent_replay=False, completed_journal_replay=False, whole_resource_rehash=False,
            new_verification_checkpoint=False, skip_verification_cli=False):
        raise Track1RosterError('completed-input verification policy changed')
    pilot = declaration['pilot']
    if ([row['sample_size'] for row in pilot['strata']] != [20000, 16336]
            or pilot['seed'] != 'track1-focused-20261008-v1'
            or pilot['scientific_output_invocations'] != 1
            or pilot['deep_saved_product_readback'] is not False
            or pilot['determinism_scope'] != 'representative_small_fixtures'):
        raise Track1RosterError('analysis profile pilot binding changed')
    # These are evidence/plan locators. The pilot owner authenticates the actual
    # source bytes; static profile selection must not acquire a large sample or
    # depend on development-only OpenSpec files being in an operator package.
    for binding in (pilot['selection'], pilot['dense_selection'],
                    declaration['gate_contract']['design'], declaration['gate_contract']['scenarios']):
        _digest(binding['sha256'], 'profile evidence SHA-256')
    amendment = declaration['prospective_amendment']
    budget = amendment['finite_budget']
    if (amendment['schema'] != 'genosapiens.track1-analysis-profile-acceptance-amendment'
            or amendment['version'] != 1
            or amendment['planning_runtime_dependency'] is not False
            or amendment['original_cohorts_unchanged'] is not True
            or amendment['original_outcomes_unchanged'] is not True
            or amendment['controlled_use_authorized'] is not False
            or amendment['qualified_profiles'] != []
            or amendment['no_automatic_retries'] is not True
            or type(budget['process_seconds']) is not int
            or type(budget['host_seconds']) is not int
            or type(budget['cleanup_seconds']) is not int
            or not 0 < budget['process_seconds'] < budget['host_seconds']
            or budget['host_seconds'] - budget['process_seconds'] != budget['cleanup_seconds']
            or budget['base_seconds'] != sum(budget['allowance_seconds'].values())
            or budget['process_seconds'] != budget['base_seconds'] + budget['uncertainty_seconds'] + budget['rounded_allowance_seconds'] + budget.get('operational_reserve_seconds', 0)
            or pilot['finite_process_timeout_seconds'] != budget['process_seconds']
            or pilot['finite_host_timeout_seconds'] != budget['host_seconds']
            or budget['total_wgs_runtime_gate'] is not False):
        raise Track1RosterError('prospective profile acceptance budget or history is invalid')
    _digest(amendment['previous_profile']['sha256'], 'previous profile SHA-256')
    _digest(amendment['previous_profile_acceptance_sha256'], 'previous declaration SHA-256')
    # Planning snapshots are qualification evidence, not installed runtime inputs.
    # The pinned current declaration binds their identities without opening them.
    for binding in amendment['current_planning_snapshots']:
        _digest(binding['sha256'], 'immutable planning snapshot SHA-256')
    _validate_d41_acceptance(declaration['d41_acceptance'])


def _validate_d41_acceptance(amendment):
    """Check prospective policy/scope authority without loading test fixtures.

    The pinned declaration binds the finite matrix and immutable planning.
    Those are evaluator evidence, not dependencies of an installed runtime.
    This declaration alone cannot qualify or activate a candidate release.
    """
    from genosapiens.population_filter_candidate import WORK_POLICY_PATH, WORK_POLICY_SHA256
    from genosapiens.wgs_work_contract import PATH as WORK_CONTRACT_PATH, load_focused_work_contract
    contract=load_focused_work_contract()
    if (amendment['schema'] != 'genosapiens.track1-d41-acceptance-amendment'
            or amendment['version'] != 1
            or amendment['policy']['path'] != WORK_POLICY_PATH
            or amendment['policy']['sha256'] != WORK_POLICY_SHA256
            or amendment['execution_contract'] != dict(path=WORK_CONTRACT_PATH,
                sha256=contract.sha256,size_bytes=len(contract.content))
            or amendment['fresh_invocation_id'] != 'focused-fresh-subsample'
            or amendment['planning_runtime_dependency'] is not False
            or amendment['qualified_profiles'] != []
            or amendment['controlled_use_authorized'] is not False
            or amendment['exhaustive_full_public_required'] is not True
            or amendment['focused_qualification'] != dict(
                scope='bounded_fresh_input_end_to_end', maximum_background_source_records=512,
                positive_edge_roster='separately_frozen', execution_cap_seconds=18000,
                budget_scope='per_three_case_trial',
                termination_grace='separately_declared_no_scientific_work',
                full_public_focused_required=False, full_wgs_runtime_validated=False,
                automatic_retries=0)):
        raise Track1RosterError('D41 policy or qualification scope changed')
    for binding in (amendment['previous_profile'], amendment['finite_acceptance'],
                    *amendment['planning_snapshots']):
        _digest(binding['sha256'], 'D41 immutable acceptance evidence SHA-256')


def validate_context_restart_acceptance(inputs):
    """Check the externally frozen, exact producer/consumer compatibility pair.

    The declaration is outside the executable stack, avoiding a circular digest.
    It is authenticated by the public plan and authority before staging.
    """
    from genosapiens.wgs_profile_projection import (
        CompletedProfileDiscoveryInputs, _signature, profile_adapter_identity,
    )
    from genosapiens.wgs_shards import read_bounded_metadata
    try:
        if type(inputs) is not CompletedProfileDiscoveryInputs:
            raise ValueError
        identity = inputs.identity
        descriptor = identity['acceptance_amendment']
        path = Path(descriptor['path'])
        before = _signature(path)
        content = read_bounded_metadata(path)
        if (len(content) != descriptor['size_bytes'] or _sha256(content) != descriptor['sha256']
                or _signature(path) != before):
            raise ValueError
        compatibility = json.loads(content)['context_restart']
        if (set(compatibility) != {'inputs', 'producer_profile_acceptance_sha256',
                'consumer_profile_adapter_sha256', 'profile_schema_version',
                'association_schema_version', 'association_codec',
                'producer_discovery_policy', 'consumer_discovery_policy_sha256'}
                or compatibility['inputs'] != {k:v for k,v in identity.items() if k != 'acceptance_amendment'}
                or compatibility['profile_schema_version'] != 3
                or compatibility['association_schema_version'] not in ((4,5) if inputs.profile_only else (5,))
                or compatibility['association_codec'] != 'canonical-row-v1/zlib-1/base64'
                or compatibility['consumer_profile_adapter_sha256'] != profile_adapter_identity()):
            raise ValueError
        _digest(compatibility['producer_profile_acceptance_sha256'], 'producer acceptance SHA-256')
        _digest(compatibility['consumer_discovery_policy_sha256'], 'consumer discovery policy SHA-256')
        if type(compatibility['producer_discovery_policy']) is not dict:
            raise ValueError
        return compatibility
    except Exception:
        raise Track1RosterError('completed context input compatibility is invalid') from None


def compatible_context_discovery_policy(inputs, current_policy_bytes):
    """Bind both exact executable maps; all other discovery policy bytes agree."""
    from genosapiens.wgs_ledger import canonical
    from genosapiens.candidate_admission import _EXECUTABLE_FILES
    compatibility = validate_context_restart_acceptance(inputs)
    if inputs.profile_only:
        raise Track1RosterError('profile-only restart cannot reuse completed discovery')
    original = compatibility['producer_discovery_policy']
    current = json.loads(current_policy_bytes)
    if (_sha256(current_policy_bytes) != compatibility['consumer_discovery_policy_sha256']
            or {k:v for k,v in original.items() if k != 'executables'}
                != {k:v for k,v in current.items() if k != 'executables'}
            or set(original.get('executables', {})) != set(_EXECUTABLE_FILES)
            or set(current.get('executables', {})) != set(_EXECUTABLE_FILES)):
        raise Track1RosterError('completed discovery policy semantics differ')
    for digest in original['executables'].values():
        _digest(digest, 'producer discovery executable SHA-256')
    return _sha256(canonical(original))


def load_analysis_profile_acceptance(repository_root):
    """Authenticate runtime scope without opening evaluator fixture artifacts.

    The full profile and its declaration retain their original pinned identities.
    Callers separately authenticate the cohort and selection receipts they use.
    """
    try:
        repo = Path(repository_root).absolute()
        data, content = load_bounded_json(repo / SUBMISSION_PROFILE_PATH, reference_root=repo)
        if _sha256(content) != SUBMISSION_PROFILE_SHA256:
            raise ValueError
        _validate_profile_acceptance(repo, data)
        return data['analysis_profile_acceptance']
    except Exception:
        raise Track1RosterError("analysis profile acceptance is missing, stale or invalid") from None


def load_submission_acceptance_profile(repository_root):
    """Authenticate the predeclared current scope; never infer it from outcomes."""
    try:
        repo = Path(repository_root).absolute()
        data, content = load_bounded_json(repo / SUBMISSION_PROFILE_PATH, reference_root=repo)
        if (_sha256(content) != SUBMISSION_PROFILE_SHA256
                or data["schema"] != "genosapiens.track1-submission-acceptance-profile"
                or data["classification"] != "public"
                or data["controlled_use_authorized"] is not False
                or data["experiment_id"] != "track1-submission-acceptance-population-models-20260917"
                or data["candidate_budget"] != 1 or data["attempt_budget_per_case"] != 1
                or data["worker_timeout_seconds"] != 6480):
            raise ValueError
        policy = read_bounded_bytes(repo / data["policy_path"], reference_root=repo)
        if _sha256(policy) != data["policy_sha256"]:
            raise ValueError
        family = json.loads(policy)
        if (family["experiment_id"] != data["experiment_id"]
                or family["candidates"] != [{"policy_id": "track1-public-balanced-v1", "phenotype_weight": "1"}]
                or family["candidate_budget"] != 1):
            raise ValueError
        from genosapiens.track1_public_validation import (
            MODEL_LOCUS_DESCRIPTOR_SHA256, CLINVAR_MANIFEST_SHA256,
            GNOMAD_JOINT_MANIFEST_SHA256, GNOMAD_ALL_SITES_MANIFEST_SHA256,
        )
        if (data["model_loci_descriptor_sha256"] != MODEL_LOCUS_DESCRIPTOR_SHA256
                or data["consumer_manifest_sha256s"] != sorted((CLINVAR_MANIFEST_SHA256,
                    GNOMAD_JOINT_MANIFEST_SHA256, GNOMAD_ALL_SITES_MANIFEST_SHA256))):
            raise ValueError
        for item in data["configuration_bindings"]:
            if _sha256(read_bounded_bytes(repo / item["path"], reference_root=repo)) != item["sha256"]:
                raise ValueError
        qualification_path = "artifacts/validation/track1-singleton-fixture-qualification-2026-09-17.json"
        bindings = [item for item in data["configuration_bindings"] if item["path"] == qualification_path]
        if len(bindings) != 1:
            raise ValueError
        qualification = json.loads(read_bounded_bytes(repo / qualification_path, reference_root=repo))
        approved_singletons = {"development-source-backed-singleton", "held-out-source-backed-singleton"}
        if (qualification["schema"] != "genosapiens.track1-singleton-fixture-qualification"
                or qualification["status"] != "source-and-annotation-qualified"
                or qualification["classification"] != "public"
                or qualification["controlled_use_authorized"] is not False
                or qualification["source_roster"] != data["retained_source_roster"]
                or qualification["source_required_allele_count"] != 1
                or qualification["source_supported_allelic_model"] != "monoallelic_autosomal"
                or qualification["ranking_outcomes_used_for_selection"] is not False
                or qualification["scientific_policy_modified"] is not False
                or set(qualification["new_case_ids"]) != approved_singletons
                or {row["case_id"] for row in data["cases"]
                    if row["positive_exact_export"] and row["target_kind"] == "singleton"}
                   != approved_singletons):
            raise ValueError
        _authenticate_population_fixture_qualification(repo, data, qualification)
        _validate_profile_acceptance(repo, data)
        return data, content
    except Exception:
        raise Track1RosterError("submission acceptance profile is missing, stale or invalid") from None


def load_submission_acceptance_roster(*, repository_root):
    """Project the exact pinned selection only after reauthenticating its source."""
    try:
        repo = Path(repository_root).absolute()
        profile, content = load_submission_acceptance_profile(repo)
        binding = profile["source_roster"]
        source_path = repo / binding["path"]
        source_content = read_bounded_bytes(source_path, reference_root=repo,
            expected_size_bytes=binding["size_bytes"])
        if _sha256(source_content) != binding["sha256"]:
            raise ValueError
        source = load_track1_roster(source_path, repository_root=repo)
        if (source.annotation_sha256 != profile["source_annotation_sha256"]
                or source.causal_truth_sha256 != profile["source_causal_truth_sha256"]):
            raise ValueError
        raw = json.loads(source_content)
        original = {c["case_id"]: c for c in raw["cases"]}
        parsed = {c.case_id: c for c in source.cases}
        selected = []
        groups = {split: set() for split in SPLITS}
        genotype_hashes = {split: set() for split in SPLITS}
        seen = set()
        positives = {split: set() for split in SPLITS}
        required = set()
        for row in profile["cases"]:
            name = row["case_id"]
            case, declaration = parsed[name], original[name]
            if (name in seen or row["source_case_sha256"] != _sha256(_case_bytes(declaration))
                    or row["split"] != case.split or row["scope"] != case.scope
                    or row["group_id"] != declaration["group_id"]
                    or row["generator"] != declaration["generator"]
                    or row["runner_inputs"] != declaration["runner_inputs"]
                    or row["expected_genotype_sha256"] != _sha256(_case_bytes(case.truth_members))
                    or row["target_kind"] != case.target_kind
                    or row["required_model_coverage"] != list(case.required_model_coverage)):
                raise ValueError
            groups[case.split].add(row["group_id"])
            genotype_hashes[case.split].add(row["expected_genotype_sha256"])
            if row["positive_exact_export"]:
                if case.target_kind not in {"singleton", "full_pair"}:
                    raise ValueError
                positives[case.split].add(case.target_kind)
            if case.split == "development":
                required.update(case.required_model_coverage)
            selected.append(case); seen.add(name)
        if (seen != set(profile["declared_case_ids"])
                or list(profile["declared_case_ids"]) != sorted(seen)
                or groups["development"] & groups["held_out"]
                or genotype_hashes["development"] & genotype_hashes["held_out"]
                or any(value != {"singleton", "full_pair"} for value in positives.values())
                or not set(REQUIRED_MODEL_COVERAGE) <= required
                or profile["assessed_scopes"] != sorted({c.scope for c in selected})
                or len(selected) != profile["case_budget"]):
            raise ValueError
        return Track1Roster(_sha256(content), source.annotation_sha256, source.causal_truth_sha256,
            tuple(sorted(selected, key=lambda c: c.case_id)), source.shared_background)
    except Exception:
        raise Track1RosterError("submission acceptance roster does not replay its declared source scope") from None


def _case_bytes(value):
    return (json.dumps(value, sort_keys=True, separators=(",", ":")) + "\n").encode("ascii")
