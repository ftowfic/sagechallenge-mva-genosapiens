"""Source-bound public genotype admission and deterministic input dossier."""

from dataclasses import dataclass
from datetime import date, datetime
import hashlib
import json
from pathlib import Path
from typing import Any, Mapping

from genosapiens.json_schema import validate_json_schema
from genosapiens.track2_input import AnalysisOrigin, Track2InputError

SCHEMA_PATH = Path(__file__).resolve().parents[2] / "schemas/track2/published-patient-genotype.schema.json"


def canonical(value: Any) -> bytes:
    return (json.dumps(value, sort_keys=True, separators=(",", ":"), allow_nan=False) + "\n").encode()


def sha256(content: bytes) -> str:
    return hashlib.sha256(content).hexdigest()


def scenario_identity(data: Mapping[str, Any]) -> str:
    identity = {"genotype": data["genotype"], "disease_context": data["disease_context"],
                "source_ids": sorted(source["source_id"] for source in data["sources"])}
    return "published-genotype:" + sha256(canonical(identity))[:24]


def _source_file(root: Path, relative: str) -> Path:
    rel = Path(relative)
    path = root / rel
    if rel.is_absolute() or ".." in rel.parts or "\\" in relative or not path.is_file():
        raise Track2InputError("source document path is invalid")
    if any(parent.is_symlink() for parent in (path, *path.parents)):
        raise Track2InputError("source document cannot be a symlink")
    return path


def _indexed(records, field):
    result = {item[field]: item for item in records}
    if len(result) != len(records):
        raise Track2InputError("duplicate public record identity")
    return result


@dataclass(frozen=True)
class PublishedGenotype:
    data: Mapping[str, Any]
    input_sha256: str
    source_hashes: Mapping[str, str]
    review: Mapping[str, Any] | None
    review_sha256: str | None
    input_report_ready: bool

    def admitted_origin(self) -> AnalysisOrigin:
        if not self.input_report_ready or self.review_sha256 is None:
            raise Track2InputError("public genotype requires actual source-bound review")
        return AnalysisOrigin("published_patient_genotype", self.data["scenario_id"],
                              self.input_sha256, self.review_sha256)


def validate_genotype(data: Mapping[str, Any], root: Path, *, allow_synthetic=False,
                      retained_source_hashes=None) -> dict[str, str]:
    validate_json_schema(data, json.loads(SCHEMA_PATH.read_text()))
    if data["evidence_origin"] != "primary_literature" and not allow_synthetic:
        raise Track2InputError("synthetic evidence cannot establish a public research report")
    if data["scenario_id"] != scenario_identity(data):
        raise Track2InputError("published scenario identity mismatch")
    sources = _indexed(data["sources"], "source_id")
    evidence = _indexed(data["evidence"], "evidence_id")
    patients = _indexed(data["patients"], "label")
    models = _indexed(data["models"], "model_id")
    source_hashes = {}
    for identifier, source in sources.items():
        date.fromisoformat(source["accessed_on"])
        if not source["public_url"].startswith("https://"):
            raise Track2InputError("public source URL must use HTTPS")
        digest = (sha256(_source_file(root, source["path"]).read_bytes()) if retained_source_hashes is None
                  else retained_source_hashes.get(identifier))
        if digest != source["sha256"]:
            raise Track2InputError("source document checksum mismatch")
        source_hashes[identifier] = digest

    def locator(item):
        if item["source_id"] not in sources:
            raise Track2InputError("unresolved source locator")

    def references(ids, *, kind=None):
        if not ids or set(ids) - evidence.keys():
            raise Track2InputError("missing or unresolved genotype evidence")
        if kind and any(evidence[key]["kind"] != kind for key in ids):
            raise Track2InputError("genotype admission requires primary genetic evidence")

    kinds = {item["kind"] for item in evidence.values()}
    if not {"genetic", "functional"} <= kinds:
        raise Track2InputError("primary genetic and functional causal evidence are required")
    for item in evidence.values():
        if set(item["patient_labels"]) - patients.keys() or set(item["model_ids"]) - models.keys():
            raise Track2InputError("unresolved patient or model attribution")
        for link in item["locators"]:
            locator(link)
        if item["kind"] in {"genetic", "functional"} and not any(
                sources[link["source_id"]]["primary"] for link in item["locators"]):
            raise Track2InputError("causal evidence must resolve to primary literature")
        if item["kind"] == "genetic" and not item["patient_labels"]:
            raise Track2InputError("genetic evidence requires reported patient attribution")
        if item["kind"] == "functional" and not item["model_ids"]:
            raise Track2InputError("functional evidence requires experimental model attribution")
        for model_id in item["model_ids"]:
            label = models[model_id]["patient_label"]
            if label and label not in item["patient_labels"]:
                raise Track2InputError("assay attribution differs from its patient-derived model")
    for model in models.values():
        if model["patient_label"] is not None and model["patient_label"] not in patients:
            raise Track2InputError("unknown model patient")
        locator(model["locator"])
    genotype = data["genotype"]
    alleles = genotype["alleles"]
    _indexed(alleles, "allele_id")
    if len({allele["gene"] for allele in alleles}) != 1 or sum(a["copies"] for a in alleles) != 2:
        raise Track2InputError("complete observed diploid genotype is required")
    if genotype["zygosity"] == "homozygous":
        valid = len(alleles) == 1 and alleles[0]["copies"] == 2 and genotype["phase"] == "not_applicable_homozygous"
    else:
        valid = len(alleles) == 2 and all(a["copies"] == 1 for a in alleles) and genotype["phase"] == "confirmed_trans"
    if not valid:
        raise Track2InputError("observed allele combination or phase is incomplete")
    for allele in alleles:
        references(allele["evidence_ids"], kind="genetic")
        if allele["transcript_version"] is not None and allele["transcript_id"] is None:
            raise Track2InputError("transcript version lacks a reference identity")
    for label, patient in patients.items():
        references(patient["genotype_evidence_ids"], kind="genetic")
        if any(label not in evidence[key]["patient_labels"] for key in patient["genotype_evidence_ids"]):
            raise Track2InputError("genotype evidence does not describe this patient")
        for key in patient["assay_evidence_ids"]:
            if key not in evidence or evidence[key]["kind"] != "functional" or label not in evidence[key]["patient_labels"]:
                raise Track2InputError("assay evidence was assigned to the wrong patient")
        expected = {key for key, item in evidence.items() if item["kind"] == "functional" and label in item["patient_labels"]}
        if set(patient["assay_evidence_ids"]) != expected:
            raise Track2InputError("patient assay attribution is incomplete")
        for observation in patient["observations"]:
            locator(observation["locator"])
    if retained_source_hashes is not None and source_hashes != retained_source_hashes:
        raise Track2InputError("retained public source identities differ")
    return source_hashes


def load_published_genotype(root: str | Path, *, require_review=True, allow_synthetic=False,
                            expected_input_sha256=None, expected_review_sha256=None) -> PublishedGenotype:
    root = Path(root)
    raw = (root / "genotype.json").read_bytes()
    data = json.loads(raw)
    input_hash = sha256(raw)
    if expected_input_sha256 is not None and input_hash != expected_input_sha256:
        raise Track2InputError("frozen genotype binding mismatch")
    hashes = validate_genotype(data, root, allow_synthetic=allow_synthetic)
    review_path = root / "review.json"
    review, review_hash, ready = None, None, False
    if review_path.exists():
        raw_review = review_path.read_bytes()
        review_hash, review = sha256(raw_review), json.loads(raw_review)
        fields = {"schema", "reviewer", "reviewer_type", "independence", "reviewed_at", "input_sha256",
                  "source_hashes", "evidence_ids", "decision", "rationale", "limitations"}
        if set(review) != fields or review["schema"] != "genosapiens.published-genotype-review/v1":
            raise Track2InputError("invalid public review contract")
        if review["input_sha256"] != input_hash or review["source_hashes"] != hashes:
            raise Track2InputError("review no longer binds the source and genotype bytes")
        if set(review["evidence_ids"]) != {e["evidence_id"] for e in data["evidence"]}:
            raise Track2InputError("review does not cover the curated evidence")
        if not review["reviewer"] or not review["rationale"] or not review["limitations"]:
            raise Track2InputError("review provenance and limitations are required")
        if datetime.fromisoformat(review["reviewed_at"].replace("Z", "+00:00")).tzinfo is None:
            raise Track2InputError("review timestamp requires a time zone")
        if review["reviewer_type"] not in {"human", "ai_assisted_extraction"} or review["independence"] not in {"independent", "self_review_disclosed"}:
            raise Track2InputError("undefined review provenance")
        if review["decision"] not in {"accepted", "rejected", "pending_human_review"}:
            raise Track2InputError("undefined review decision")
        if review["reviewer_type"] != "human" and review["decision"] == "accepted":
            raise Track2InputError("AI extraction cannot impersonate manual scientific review")
        ready = review["decision"] == "accepted" and review["reviewer_type"] == "human"
    if expected_review_sha256 is not None and review_hash != expected_review_sha256:
        raise Track2InputError("frozen review binding mismatch")
    if require_review and not ready:
        raise Track2InputError("public genotype requires actual source-bound review")
    return PublishedGenotype(data, input_hash, hashes, review, review_hash, ready)


def render_input_dossier(bundle: PublishedGenotype) -> str:
    data = bundle.data
    lines = ["# Published genotype and mechanism dossier", "",
             f"Input report ready: **{str(bundle.input_report_ready).lower()}**", "",
             "Origin: published patient genotype. Challenge-case matching: not assessed. Submission ready: false.", "",
             f"Scenario: `{data['scenario_id']}`", f"Input SHA-256: `{bundle.input_sha256}`", "",
             f"Disease context: {data['disease_context']}", "", "## Observed genotype", ""]
    for allele in data["genotype"]["alleles"]:
        lines += [f"{allele['gene']}: {allele['reported_coding']} ({allele['reported_protein']}), "
                  f"{allele['copies']} copies. Zygosity: {data['genotype']['zygosity']}.", "",
                  f"Normalized display: {allele['normalized_coding']} ({allele['normalized_protein']}). "
                  + allele["normalization_rationale"], "", "Reference mapping: " + allele["reference_mapping"], ""]
    lines += ["## Published patient attribution", "", "Each label identifies a row in the source publication only. This dossier keeps each patient's observations separate.", ""]
    for patient in data["patients"]:
        lines += [f"### {patient['label']}", ""]
        for obs in patient["observations"]:
            loc = obs["locator"]
            lines += [f"- {obs['name']}: {obs['value']} ({loc['source_id']}, PDF p. {loc['page']}, {loc['section']})."]
        lines += ["", "Assays attributed in this dossier: " + (", ".join(patient["assay_evidence_ids"]) or "none documented here. This does not establish that no assay occurred") + ".", ""]
    lines += ["## Genetic, functional and mechanistic evidence", ""]
    for item in data["evidence"]:
        lines += [f"### {item['evidence_id']} ({item['kind']})", "", item["statement"], "",
                  "Patients: " + (", ".join(item["patient_labels"]) or "not patient-specific") + ". "
                  + "Models: " + (", ".join(item["model_ids"]) or "not applicable") + ".", ""]
        for loc in item["locators"]:
            lines += [f"- {loc['source_id']}, PDF p. {loc['page']}, {loc['section']}. Attribution: {loc['attribution']}."]
        lines += ["", "Limitations:", *['- ' + value for value in (item['limitations'] or ['See study-wide limitations.'])], ""]
    lines += ["## Scope and limitations", ""] + ["- " + value for value in data["limitations"]]
    lines += ["", "## Source manifest", ""]
    for source in data["sources"]:
        lines += [f"- [{source['citation']}]({source['public_url']}). Version: {source['document_version']}. SHA-256 `{source['sha256']}`. "
                  + f"Terms: {source['license_or_terms']}. Raw redistribution: {source['raw_redistribution']}."]
    lines += ["", "## Review", "", ("No scientific review supplied." if bundle.review is None else
              f"Decision: {bundle.review['decision']}. Reviewer: {bundle.review['reviewer']} ({bundle.review['reviewer_type']}). "
              f"Independence: {bundle.review['independence']}. Time: {bundle.review['reviewed_at']}.\n\nReview text (verbatim):\n\n" + bundle.review['rationale']), "",
              "This dossier supports research input review. Genetic complementation is a mechanism comparator, not medication efficacy or prescribing advice.", ""]
    if bundle.review is not None:
        lines += ["Review SHA-256: `" + bundle.review_sha256 + "`.", "",
                  "Review scope and limitations:", ""]
        lines += ["- " + item for item in bundle.review["limitations"]] + [""]
    return "\n".join(lines)
