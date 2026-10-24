"""Bounded matching-diagnosis supplement and independent arithmetic readback."""

from dataclasses import dataclass
from fractions import Fraction
from pathlib import Path

from contextlib import contextmanager

from genosapiens.track1_sensitivity import DisclosedDiagnosis
from genosapiens.wgs_failures import WgsFailure
from genosapiens.wgs_scientific_products import open_scientific_products, product_descriptors, validate_scientific_products
from genosapiens.wgs_diagnosis_derivative import write_diagnosis_derivative, open_diagnosis_derivative


CHECKS = ("fixture_source_support", "target_in_primary", "positive_prior_applied", "exact_prior_arithmetic",
    "negative_control_not_applicable", "primary_unchanged", "no_secondary_evidence_or_dispatch")


def _require(value):
    if not value:
        raise WgsFailure("bounded_products", "malformed_evidence", operation="bounded_product_validate", exception_kind="project_failure")


@dataclass(frozen=True, slots=True)
class MatchingExpectation:
    matching_diagnosis: str
    control_diagnosis: str
    gene_id: str
    gene_record_id: str
    gene_dataset: str
    target: tuple[str, int, str, str]
    hpo_aspects: tuple[tuple[str, str], ...]
    rare_disease_manifest_sha256: str
    vcf_input_sha256: str
    hpo_input_sha256: str
    hpo_manifest_sha256: str
    scientific_policy_sha256: str
    ic_sha256: str


@contextmanager
def _primary_reader(execution, scratch):
    products = execution.products
    current = getattr(execution, "scientific_reader", None)
    if current is not None:
        with current.reader.read_scope():
            _require(current.artifact == products.scientific)
            validate_scientific_products(current, product_descriptors(products),
                expected_prefilter_manifest_sha256=products.primary.prefilter_manifest_sha256)
            yield current
    else:
        with open_scientific_products(products.scientific, product_descriptors(products),
                directory=products.scientific_directory, expected_parent_bindings=dict(products.scientific.parent_bindings),
                scratch_directory=Path(scratch)) as reader:
            yield reader


def evaluate_matching_secondary(execution, *, resources, diagnosis_id, output_directory):
    """One authorized post-primary sensitivity over unchanged complete roots."""
    products = execution.products
    before = tuple(getattr(products, role).as_bytes() for role in (
        "scientific", "order", "primary", "selection"))
    output = Path(output_directory)
    output.mkdir(mode=0o700, parents=True, exist_ok=True)
    with _primary_reader(execution, output.parent) as reader:
        matching = write_diagnosis_derivative(reader, products.primary, DisclosedDiagnosis(diagnosis_id), output_directory=output)
        # The requested negative control has a separate optional owner as well.
        write_diagnosis_derivative(reader, products.primary, reader.resources.diagnosis, output_directory=output / "control")
    _require(before == tuple(getattr(products, role).as_bytes() for role in (
        "scientific", "order", "primary", "selection")))
    return matching


def check_matching_science(execution, matching, *, matching_directory, control, control_directory, resources, expected, scratch_path):
    """Authenticate complete prior/order replay and independently expected sources."""
    products, primary = execution.products, execution.products.primary
    _require(primary.prior_policy.strength == "1"
        and primary.prior_policy.derivation_sha256 == expected.scientific_policy_sha256
        and primary.ic_sha256 == expected.ic_sha256
        and matching.status == "applied" and matching.reason is None
        and matching.diagnosis == DisclosedDiagnosis(expected.matching_diagnosis)
        and control.diagnosis == DisclosedDiagnosis(expected.control_diagnosis)
        and control.status == "not_applicable" and control.reason == "no_prior_match")
    scratch = Path(scratch_path).parent
    scratch.mkdir(mode=0o700, parents=True, exist_ok=True)
    with _primary_reader(execution, scratch) as reader:
        genes = reader.resources.genes.prepared.evidence
        gene_supported = False
        for owner in genes:
            if owner.public_manifest_sha256 != expected.rare_disease_manifest_sha256:
                continue
            record = owner.store.records.get(expected.gene_record_id)
            if record is None:
                continue
            fields = record.source_fields
            _require(fields.get("source_disease_id") == expected.matching_diagnosis
                and fields.get("ensembl_gene_id") == expected.gene_id and fields.get("source_dataset") == expected.gene_dataset)
            gene_supported |= any(assertion.record_id == record.record_id and assertion.gene_id == expected.gene_id and any(
                row.required_allele_count == 1 and row.allelic_requirement == "monoallelic_autosomal"
                for row in assertion.requirements) for assertion in owner.store.assertions)
        terms = set()
        profile = reader.resources.phenotypes.profile
        _require(profile is not None and profile.public_manifest_sha256 == expected.rare_disease_manifest_sha256)
        for row in profile.store.annotations:
            fields = row.source_fields
            if fields.get("source_disease_id") != expected.matching_diagnosis:
                continue
            term = fields.get("source_hpo_id")
            if term in dict(expected.hpo_aspects):
                _require(fields.get("aspect") == dict(expected.hpo_aspects)[term]
                    and row.polarity == "positive" and row.resolution.canonical_id == term
                    and all(fields.get(key) for key in ("source_release", "source_dataset", "assertion_key", "lineage_key", "dependency_key")))
                terms.add(term)
        _require(gene_supported and terms == set(dict(expected.hpo_aspects)))
        observations = reader.resources.phenotypes.observations
        _require(len(observations) == len(expected.hpo_aspects)
            and {r.scoring.source_hpo_id for r in observations} == set(dict(expected.hpo_aspects)))
        for row in observations:
            _require(row.scoring.scope == "proband" and row.scoring.status == "present"
                and row.scoring.provenance.input_sha256 == expected.hpo_input_sha256
                and row.ontology_manifest_sha256 == expected.hpo_manifest_sha256)
        db = reader.reader.connection
        db.execute("CREATE TABLE matching_target(allele TEXT PRIMARY KEY)")
        try:
            for member in reader.iter_member_records():
                allele = member["member"].allele
                for identity in allele.identities:
                    native = identity.identity
                    if (identity.stage == "native" and (getattr(native, "contig", None), getattr(native, "position", None),
                            getattr(native, "reference", None), getattr(native, "alternate", None)) == expected.target):
                        _require(any((call.source.source_fields or {}).get("source_input_sha256") == expected.vcf_input_sha256
                            and any(v not in (".", "0") for v in call.sample.genotype.replace("|", "/").split("/")) for call in allele.calls))
                        db.execute("INSERT OR IGNORE INTO matching_target VALUES(?)", (allele.allele_id,))
            _require(db.execute("SELECT count(*) FROM matching_target").fetchone()[0] > 0)
            target_proven = False
            with open_diagnosis_derivative(matching, reader, primary, directory=matching_directory, scratch_directory=scratch) as derivative:
                for ordered in derivative.ordered():
                    prior = derivative.decision(ordered["locator"][0])
                    primary_total = str(Fraction(ordered["total"]) - Fraction(prior["contribution"]))
                    rank = reader.case.rank_of(ordered["locator"], total=primary_total, ties=ordered["ties"])
                    original = reader.case.ordered_at(rank)
                    _require(original["hypothesis_id"] == ordered["hypothesis_id"]
                        and Fraction(ordered["total"]) == Fraction(original["total"]) + Fraction(prior["contribution"]))
                    if any(db.execute("SELECT 1 FROM matching_target WHERE allele=?", (key,)).fetchone()
                            for key in reader.case.identities.subject_ids(ordered["locator"])):
                        target_proven |= any(proof["match"] == "exact_native"
                            and proof["native_disease_id"] == expected.matching_diagnosis
                            and proof["target_gene_id"] == expected.gene_id for proof in prior["proofs"])
            _require(target_proven)
            # Replay the explicitly published control against unchanged primary roots.
            _require(reader.resources.diagnosis == control.diagnosis)
            with open_diagnosis_derivative(control, reader, primary, directory=control_directory, scratch_directory=scratch) as derivative:
                for original, altered in zip(reader.ordered(), derivative.ordered(), strict=True):
                    _require({k:v for k,v in original.items() if k != "mode"} == altered)
        finally:
            db.execute("DROP TABLE matching_target")
    hashes = tuple(sorted({matching.sha256, control.sha256, *(value for _, value in execution.manifest.output_bindings)}))
    return {name: {"status": "passed", "supporting_artifact_sha256s": list(hashes)} for name in CHECKS}
