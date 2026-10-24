"""Bounded public-only Track 1 source-to-ranking execution.

This module materializes only the hashes present in a runner-safe case view.  It
does not accept causal truth, case recipes, or a scope label.  Development and
held-out coordination remains in :mod:`track1_validation_runner`; this module
only prepares and executes one already-authorized public view.
"""

from __future__ import annotations

from collections.abc import Mapping
from dataclasses import asdict, dataclass
from fractions import Fraction
import hashlib
import json
from pathlib import Path
import re
import resource
import socket
import time
from typing import get_args

from genosapiens import candidate_admission as admission
from genosapiens import evidence_fusion as fusion
from genosapiens.challenge_evaluator import load_pinned_evaluator
from genosapiens.admission_inventory import (
    SuppliedCallInputs,
    SuppliedSourceSetInputs,
    build_admission_inventory,
    source_set_manifest_bytes,
)
from genosapiens.admission_links import PinnedPublicSource
from genosapiens.disease_gene_evidence import (
    PUBLIC_TABLES,
    AggregationPolicy,
    DiseaseGeneEvidenceStore,
    build_disease_gene_evidence,
    read_pinned_disease_gene_evidence,
)
from genosapiens.disease_profiles import (
    DiseaseProfileStore,
    build_disease_profiles,
    read_pinned_disease_profiles,
)
from genosapiens.evidence_fusion import (
    ComponentPolicy,
    ExactAffineClamp,
    FusionPolicy,
    PenaltyPolicy,
    SourceFactMapping,
    TieRule,
    policy_as_json,
)
from genosapiens.evidence_model import (
    ComponentName,
    DiseaseGeneEvidence,
    DiseaseProfile,
    SourcedRecord,
    scoring_view,
    _validation_session,
)
from genosapiens.evidence_serialization import from_json
from genosapiens.genotype_hypotheses import (
    HypothesisPopulationInputs,
    HypothesisPolicy,
    HypothesisReplay,
    hypothesis_policy_bytes,
)
from genosapiens import population_annotation as population
from genosapiens import population_policy as frequency
from genosapiens.grch38_reference import IndexedFasta
from genosapiens.hpo_ic import HpoICStore, load_completed_hpo_ic
from genosapiens.model_requirements import ModelRequirementPolicy
from genosapiens.model_locus_resources import (
    ModelLocusResource, authenticated_model_loci, load_model_locus_resource,
)
from genosapiens.partner_qualification import (
    ClinicalQualificationPolicy,
    ClinicalStrengthPolicy,
)
from genosapiens.phenotype_ingestion import ingest_phenotypes
from genosapiens.reference_concordance import (
    ReferenceIdentityContract,
    ReferenceRoleBinding,
    evaluate_public_validation_vcf_concordance,
    load_construction_report,
    load_reference_identity,
    open_installed_reference,
)
from genosapiens.submission_rendering import load_submission_contig_renderer
from genosapiens.repeat_catalog_admission import (
    RepeatReferenceInput,
    load_pinned_repeat_catalog,
)
from genosapiens.track1_diagnostics import DiagnosticPolicy, SourceDiagnosticEvidence
from genosapiens.track1_epcr import EpcrPolicy, policy_from_json as epcr_policy_from_json
from genosapiens.track1_policy_bindings import PolicyBindings
from genosapiens.track1_sensitivity import (
    DOMINANT_TIE,
    DisclosedDiagnosis,
    PhenotypeTransform,
    PriorPolicy,
    _AuthenticatedSensitivityResources,
    _authenticate_loaded_ic_store,
)
from genosapiens.track1_validation_plan import (
    InputKind,
    RankerCaseView,
    RankerResource,
    ranker_case_view_as_json,
)
from genosapiens.track1_validation_runner import (
    ExperimentProtocol,
    PolicyCandidate,
    ResourceMeasurements,
)
from genosapiens.variant_inventory import build_variant_class_inventory
from genosapiens.variant_normalization import (
    normalize_public_validation_vcf_small_variants,
)
from genosapiens.vep_annotation import (
    OpenTargetsConstraintReference,
    VepReferenceCompatibility,
    VepRequest,
    run_vep_annotation,
)
from genosapiens.vcf_qc import summarize_vcf_genotype_qc


RUNTIME_IMAGE = (
    "genosapiens-track1-pipeline@sha256:"
    "640ba2fa699cc8ed677ed19f57aa39d256f05ff8377cf745f8f12131b95b0b92"
)
REFERENCE_MANIFEST_SHA256 = (
    "27eaf00fde35e5fa83f945c63e5dee871960374d90488cc0797b5ebea2fb1c25"
)
VEP_BUNDLE_SHA256 = (
    "3071600d1ab2db89a9d0988465d137c0221378ebfdedb021b5392d4a8deb075a"
)
OPEN_TARGETS_MANIFEST_SHA256 = (
    "c6edfc9b8726fd05f5f5cadc19d09ae49d6866d4742c8e9b52b038796637e55e"
)
RARE_DISEASE_BUNDLE_SHA256 = (
    "12c2bf613424137e3c8e7aec7f69ace0cf83a84a5a3bd0bcca81fc63970823d0"
)
HPO_STORE_SHA256 = (
    "f376198b4faa577910e18f1f75bdf482cd7f9810b232e29c6c44f0d1edd54efd"
)
HPO_MANIFEST_SHA256 = (
    "81bb8fe18512f05ba06a9286e6a4c0514a49cf4883d066ef31779e2cf4954232"
)
REPEAT_MANIFEST_SHA256 = (
    "7058790cfa8b76bcdc6684f18cfd3152b70a51026412641955262c35e83a95c6"
)
CLINVAR_MANIFEST_SHA256 = "78527bd4e5b7dfb3debb88e66c08ce1e3b66ab23692a258544416dffdea4e9e5"
CLINVAR_MANIFEST_PATH = "manifests/clinvar/clinvar-grch38-2026-08-22.19fc8e9cd788.json"
GNOMAD_JOINT_MANIFEST_SHA256 = "a453a49f573e4a5999b38bcd02744ab389f0ef1a4e878a9adf7fa656e59985af"
GNOMAD_ALL_SITES_MANIFEST_SHA256 = "8f1df94310ec81af4118ab6760694997a932ceaeb10469b7b91e4c73d6400305"
COMMON_ALLELE_EXCEPTION_MANIFEST_SHA256 = "99a3da5f52cae81e097634147c79d2154d20c10646277482cbe3e6c271bcce82"
COMMON_ALLELE_PUBLIC_INDEX_MANIFEST_SHA256 = "e65a6a5b6c236c794991130ab2a9cb4ffac5fc326e847e964b78a269eb3efb1d"
GNOMAD_PRECISION_REPORT_SHA256 = "4fb171a7fb5a23ed8ed24dc81b35880f84fa79948b25393048cbcdc991410303"
COMMON_ALLELE_EXCEPTION_MANIFEST_PATH = "manifests/track1-common-allele-exceptions/79596dff01ac7291604c7715fdd3fa6d7ea89257ee8172e7c2b926402184586e.json"
COMMON_ALLELE_PUBLIC_INDEX_MANIFEST_PATH = "manifests/common-allele-public-indexes/c8e7740ccc3b024c8502b035529ad2fc185a5d9dd18a32f9b7cc421d7d8fac4f-cache-bounded-runtime.json"
GNOMAD_PRECISION_REPORT_PATH = "manifests/gnomad/4.1/faf95-precision/4fb171a7fb5a23ed8ed24dc81b35880f84fa79948b25393048cbcdc991410303.json"
EVALUATOR_SOURCE_SHA256 = (
    "6d18b581e65a45e1ccc120071d588e740c2e42e983ff50704c60a40232b19180"
)
SHA256 = re.compile(r"^[0-9a-f]{64}$")
POLICY_PATH = Path("config/track1/public-ranking-development-policy-family.json")
REFERENCE_IDENTITY_PATH = Path(
    "config/reference-identities/challenge-grch38-masked-v2-hs38d1.json"
)
VEP_CONTRACT_PATH = Path("config/public-resources/ensembl-vep-116-grch38.json")
OPEN_TARGETS_CONTRACT_PATH = Path(
    "config/public-resources/open-targets-platform-26.06.json"
)
RARE_DISEASE_CONTRACT_PATH = Path(
    "config/public-resources/rare-disease-evidence-2026-08-30.json"
)
ONTOLOGY_CONTRACT_PATH = Path("config/public-resources/core-ontologies.json")
REPEAT_CONTRACT_PATH = Path("config/public-resources/variant-class-resources.json")
MODEL_LOCUS_CONTRACT_PATH = Path("config/public-resources/grch38-model-loci.json")
MODEL_LOCUS_DESCRIPTOR_SHA256 = "5b9c9760bdb5fa901e0690cf8e6ec5e5a4f6c9fa7382086652b52ca29b02d61e"
MODEL_LOCUS_MANIFEST_PATH = "manifests/model-loci/grch38-authoritative-model-loci.22b69fc15bbb.json"
MODEL_LOCUS_MANIFEST_SHA256 = "b8442cac9dcf9068ef263deab3f453a5fb8425eb8b9fd758236efd5ab5ab1eb5"


class PublicTrack1ValidationError(RuntimeError):
    """Static failure text; source rows and runner inputs are never disclosed."""

    def __init__(self, message: str, *, preparation_stage: str | None = None,
                 execution_stage: str | None = None):
        from genosapiens.track1_validation_runner import PREPARATION_STAGES, EXECUTION_STAGES
        if preparation_stage is not None and (
            type(preparation_stage) is not str or preparation_stage not in PREPARATION_STAGES
        ):
            raise ValueError("invalid public preparation stage") from None
        if execution_stage is not None and (
            type(execution_stage) is not str or execution_stage not in EXECUTION_STAGES
        ):
            raise ValueError("invalid public execution stage") from None
        if preparation_stage is not None and execution_stage is not None:
            raise ValueError("ambiguous public failure stage") from None
        super().__init__(message)
        self.preparation_stage = preparation_stage
        self.execution_stage = execution_stage


def _sha256(content: bytes) -> str:
    return hashlib.sha256(content).hexdigest()


def _canonical(value: object) -> bytes:
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


def _bound_bytes(path: Path, expected_sha256: str) -> bytes:
    try:
        if path.is_symlink() or not path.is_file():
            raise ValueError
        content = path.read_bytes()
        if _sha256(content) != expected_sha256:
            raise ValueError
        return content
    except (OSError, ValueError):
        raise PublicTrack1ValidationError("public resource binding is stale") from None


@dataclass(frozen=True, slots=True)
class PublicPolicyRuntime:
    candidate: PolicyCandidate
    fusion_policy: FusionPolicy
    fusion_policy_bytes: bytes
    prior_policy: PriorPolicy
    diagnostic_policy: DiagnosticPolicy


@dataclass(frozen=True, slots=True)
class PublicPolicyFamily:
    policy_bytes: bytes
    policy_sha256: str
    protocol: ExperimentProtocol
    hypothesis_policy: HypothesisPolicy
    hypothesis_scientific_sha256: str
    runtimes: tuple[PublicPolicyRuntime, ...]
    disclosed_diagnosis: DisclosedDiagnosis
    epcr_policy: EpcrPolicy
    epcr_policy_bytes: bytes
    epcr_policy_sha256: str
    epcr_derivation_bytes: bytes


def _policy_source(policy_sha256: str, record_id: str) -> SourcedRecord:
    return SourcedRecord(
        "Genosapiens public Track 1 development policy family",
        "1",
        record_id,
        policy_sha256,
        policy_sha256,
        {
            "schema": "genosapiens.track1-public-ranking-development-policy-family",
            "selection": "prespecified-before-development-execution",
            "classification": "public",
        },
    )


def _fusion_policy(
    policy_id: str,
    phenotype_weight: str,
    *,
    derivation_sha256: str,
) -> FusionPolicy:
    aggregations = {
        "disease_gene": "shared",
        "molecular": "member_mean",
        "inheritance": "shared",
        "population": "shared",
        "clinical": "member_mean",
        "technical_and_conflict_penalties": "shared_plus_member_mean",
    }
    components = tuple(
        ComponentPolicy(
            name,
            ExactAffineClamp(
                "1",
                "0",
                "-10" if name == "population" else "0",
                "10",
            ),
            "0",
            "0",
            "1",
            aggregations[name],
        )
        for name in get_args(ComponentName)[1:]
    )
    values = {
        "supported": "1",
        "contradicted": "0",
        "match": "1",
        "no_match": "0",
        "positive": "1",
        "negative": "0",
        "conflict": "1",
        "read_validation_required": "1",
        "truncated": "1",
        "failed": "1",
    }
    source_facts = tuple(
        SourceFactMapping(name, values.get(name.rsplit(".", 1)[-1]))
        for name in fusion._FACTS
    )
    return FusionPolicy(
        policy_id,
        PhenotypeTransform("1", "0", "0", "10", "0", phenotype_weight),
        components,
        source_facts,
        AggregationPolicy(),
        PenaltyPolicy("1", "1", "5", "2", "1", "2", "3"),
        (
            TieRule("molecular", "descending"),
            TieRule("technical_and_conflict_penalties", "ascending"),
            TieRule("completeness", "descending"),
            TieRule("phase", "descending"),
            TieRule("native_disease_id", "ascending"),
            TieRule("hypothesis_id", "ascending"),
        ),
        derivation_sha256,
        "public",
    )


def load_public_model_loci(repository_root, reference_root) -> ModelLocusResource:
    """Qualify the selected small source before any dependent annotation."""
    try:
        return load_model_locus_resource(
            Path(repository_root) / MODEL_LOCUS_CONTRACT_PATH,
            repository_root=repository_root, reference_root=reference_root,
            expected_descriptor_sha256=MODEL_LOCUS_DESCRIPTOR_SHA256,
        )
    except Exception:
        raise PublicTrack1ValidationError("public model-locus qualification failed",
            preparation_stage="manifest-qualification") from None


def _authenticated_policy_loci(root, model_loci):
    regions = authenticated_model_loci(model_loci)
    _bound_bytes(Path(root) / MODEL_LOCUS_CONTRACT_PATH, MODEL_LOCUS_DESCRIPTOR_SHA256)
    contract = load_reference_identity(Path(root) / REFERENCE_IDENTITY_PATH)
    roles = {role.role: role.identity_id for role in contract.reference_roles}
    if (not regions or model_loci.descriptor_sha256 != MODEL_LOCUS_DESCRIPTOR_SHA256
            or model_loci.manifest_sha256 != MODEL_LOCUS_MANIFEST_SHA256
            or {r.reference_role.value for r in regions} != set(roles)
            or any(r.reference_identity_id != roles[r.reference_role.value]
                or r.reference_manifest_sha256 != REFERENCE_MANIFEST_SHA256 for r in regions)):
        raise ValueError
    return regions


def load_public_policy_family(repository_root: str | Path, *, model_loci: ModelLocusResource) -> PublicPolicyFamily:
    """Bind the single declared policy only after authenticating source regions."""
    try:
        root = Path(repository_root).resolve(strict=True)
        from genosapiens.installed_public_index import load_installed_index_pin
        index_pin = load_installed_index_pin(root)
        index_manifest_sha = (index_pin['manifest_sha256'] if index_pin
                              else COMMON_ALLELE_PUBLIC_INDEX_MANIFEST_SHA256)
        loci = _authenticated_policy_loci(root, model_loci)
        policy_bytes = (root / POLICY_PATH).read_bytes()
        policy_sha256 = _sha256(policy_bytes)
        data = json.loads(policy_bytes)
        fixed = data["fixed_settings"]
        _population_template(policy_bytes)
        derivation_bytes = _bound_bytes(
            root / "config/track1/public-population-development-derivation.json",
            fixed["population"]["derivation_sha256"],
        )
        derivation = json.loads(derivation_bytes)
        if (
            derivation.get("schema") != "genosapiens.track1-public-population-development-derivation"
            or derivation.get("classification") != "public"
            or derivation.get("validation_only") is not True
            or derivation.get("controlled_use_authorized") is not False
            or derivation.get("source_manifest_sha256s") != sorted((
                GNOMAD_JOINT_MANIFEST_SHA256, OPEN_TARGETS_MANIFEST_SHA256,
                RARE_DISEASE_BUNDLE_SHA256,
            ))
            or derivation.get("selected_frequency_threshold") != fixed["population"]["frequency_rule"]["threshold"]
            or derivation.get("minimum_observed_development_an") != fixed["population"]["frequency_rule"]["minimum_an"]
            or (derivation.get("selected_subset"), derivation.get("selected_stratum"),
                derivation.get("selected_ancestry")) != ("genomes", "adjusted", None)
            or tuple(fixed["population"]["frequency_rule"][key] for key in ("subset", "stratum", "ancestry"))
               != ("genomes", "adjusted", None)
            or len(derivation["selected_low_identity_sha256s"]) != 2
            or len(set((*derivation["selected_low_identity_sha256s"], derivation["selected_high_identity_sha256"]))) != 3
            or any(SHA256.fullmatch(value) is None for value in (*derivation["selected_low_identity_sha256s"], derivation["selected_high_identity_sha256"]))
            or not Fraction(derivation["development_selected_low_frequency_upper_bound"]) < Fraction(derivation["selected_frequency_threshold"]) < Fraction(derivation["development_selected_high_frequency"])
        ):
            raise ValueError
        if (
            data["schema"]
            != "genosapiens.track1-public-ranking-development-policy-family"
            or data["classification"] != "public"
            or data["validation_only"] is not True
            or data["controlled_use_authorized"] is not False
            or data["candidate_budget"] != 1
            or data["attempt_budget_per_case"] != 1
            or data["stopping_condition"]
            != "complete-declared-family-or-budget-exhausted"
            or data["candidates"] != [{"policy_id": "track1-public-balanced-v1", "phenotype_weight": "1"}]
            or fixed["model_loci_descriptor_sha256"] != MODEL_LOCUS_DESCRIPTOR_SHA256
            or fixed["maximum_pairs_per_context"] is not None
            or fixed["dominant_term_tie"] != DOMINANT_TIE
        ):
            raise ValueError
        source = _policy_source(policy_sha256, "model-requirement-policy")
        requirements = ModelRequirementPolicy(source)
        positive = ClinicalStrengthPolicy(
            tuple(fixed["clinical_positive_significances"]),
            tuple(fixed["clinical_review_statuses"]),
            tuple(fixed["clinical_confidences"]),
        )
        negative = ClinicalStrengthPolicy(
            tuple(fixed["clinical_negative_significances"]),
            tuple(fixed["clinical_review_statuses"]),
            tuple(fixed["clinical_confidences"]),
        )
        clinical = ClinicalQualificationPolicy(
            "track1-public-clinical-strength-v1",
            _policy_source(policy_sha256, "clinical-qualification-policy"),
            positive,
            negative,
        )
        hypothesis = HypothesisPolicy(
            requirements,
            clinical,
            fixed["incompleteness_penalty"],
            maximum_pairs_per_context=fixed["maximum_pairs_per_context"],
            loci=loci,
        )
        hypothesis_sha256 = _sha256(
            admission._bytes({"hypothesis_policy": hypothesis})
        )
        prior = PriorPolicy(
            "track1-public-diagnosis-prior-v1",
            fixed["diagnosis_prior_strength"],
            policy_sha256,
            "public",
        )
        prior_sha256 = _sha256(_canonical(asdict(prior)))
        diagnostic = DiagnosticPolicy(
            "track1-public-observational-diagnostics-v1",
            None,
            None,
            None,
            (),
        )
        dominant_sha256 = _sha256(_canonical({
            "dominant_term_tie": DOMINANT_TIE,
            "phenotype_transform": fixed["phenotype_transform"],
        }))
        epcr_bytes = (root / fixed["epcr_policy_path"]).read_bytes()
        epcr_sha256 = _sha256(epcr_bytes)
        epcr = epcr_policy_from_json(epcr_bytes, expected_sha256=epcr_sha256)
        epcr_derivation = (root / fixed["epcr_suitability_path"]).read_bytes()
        if _sha256(epcr_derivation) != epcr.derivation_sha256:
            raise ValueError
        runtimes = []
        for value in data["candidates"]:
            selected = _fusion_policy(
                value["policy_id"],
                value["phenotype_weight"],
                derivation_sha256=policy_sha256,
            )
            selected_bytes = policy_as_json(selected)
            selected_sha256 = _sha256(selected_bytes)
            bindings = PolicyBindings(
                admission.POLICY_SHA256,
                hypothesis_sha256,
                selected_sha256,
                prior_sha256,
                dominant_sha256,
                diagnostic.sha256,
                epcr_sha256,
                _sha256((root / "config/track1/common-allele-prefilter-v1.json").read_bytes()),
                COMMON_ALLELE_EXCEPTION_MANIFEST_SHA256,
                index_manifest_sha,
                GNOMAD_PRECISION_REPORT_SHA256,
                _sha256((root / "config/track1/bounded-wgs-execution-v5.json").read_bytes()),
            )
            runtimes.append(PublicPolicyRuntime(
                PolicyCandidate(
                    value["policy_id"],
                    bindings,
                    tuple(sorted((
                        "clinical-qualification-source-vocabulary",
                        "complete-supplied-call-accounting",
                        "diagnosis-blind-primary",
                        "epcr-numerical-policy",
                        "population-eligibility-separate-from-qualification",
                    ))),
                    (),
                ),
                selected,
                selected_bytes,
                prior,
                diagnostic,
            ))
        if tuple(item.candidate.policy_id for item in runtimes) != tuple(sorted(
            item.candidate.policy_id for item in runtimes
        )):
            raise ValueError
        protocol = ExperimentProtocol(
            data["experiment_id"],
            data["policy_family_id"],
            tuple(data["objective_metrics"]),
            tuple(data["deterministic_selection_ties"]),
            data["candidate_budget"],
            data["attempt_budget_per_case"],
            data["stopping_condition"],
            "last-attempt-is-terminal-all-attempts-retained",
        )
        return PublicPolicyFamily(
            policy_bytes,
            policy_sha256,
            protocol,
            hypothesis,
            hypothesis_sha256,
            tuple(runtimes),
            DisclosedDiagnosis(fixed["diagnosis_sensitivity_control"]),
            epcr,
            epcr_bytes,
            epcr_sha256,
            epcr_derivation,
        )
    except Exception:
        raise PublicTrack1ValidationError("public policy family is invalid") from None


@dataclass(frozen=True, slots=True, repr=False)
class PublicResourceContext:
    repository_root: Path
    reference_root: Path
    scratch_root: Path
    reference_contract: ReferenceIdentityContract
    reference_manifest_path: Path
    reference: IndexedFasta
    source_to_reference_contigs: Mapping[str, str]
    construction_report: object
    hpo_store: HpoICStore
    disease_profiles: DiseaseProfileStore
    disease_genes: DiseaseGeneEvidenceStore
    repeat_catalog: DiseaseGeneEvidence
    repeat_reference: RepeatReferenceInput
    vep_constraint: OpenTargetsConstraintReference
    vep_compatibility: VepReferenceCompatibility
    submission_reference: PinnedPublicSource
    evaluator_source_path: Path
    common_resource_bytes: tuple[bytes, ...]
    common_runner_resources: tuple[RankerResource, ...]
    sensitivity_resources: _AuthenticatedSensitivityResources | None = None
    clinical_open_targets: PinnedPublicSource | None = None
    clinvar: PinnedPublicSource | None = None
    population_joint: population.GnomadJointReference | None = None
    population_all_sites: population.GnomadAllSitesReference | None = None
    population_reference: population.PopulationReferenceCompatibility | None = None
    model_loci: ModelLocusResource | None = None
    executable_stack_content: bytes | None = None
    public_output_root: Path | None = None
    profile_catalog: object | None = None
    regulatory_index: object | None = None
    regulatory_sources: tuple | None = None

    def __repr__(self) -> str:
        return "<public Track 1 validation resources>"


@dataclass(frozen=True, slots=True, repr=False)
class PublicPreparationResourceContext:
    """Authenticated resources whose live objects preparation actually uses."""

    repository_root: Path
    reference_root: Path
    scratch_root: Path
    reference_contract: ReferenceIdentityContract
    reference_manifest_path: Path
    reference: IndexedFasta
    source_to_reference_contigs: Mapping[str, str]
    hpo_store: HpoICStore
    clinical_open_targets: PinnedPublicSource
    clinvar: PinnedPublicSource
    population_joint: population.GnomadJointReference
    population_all_sites: population.GnomadAllSitesReference
    population_reference: population.PopulationReferenceCompatibility
    authenticated_descriptors: tuple[RankerResource, ...]

    def __repr__(self) -> str:
        return "<public Track 1 preparation resources>"


@dataclass(frozen=True, slots=True)
class _StableFileIdentity:
    device: int
    inode: int
    size: int
    modified_ns: int
    changed_ns: int


class _StageBoundBytes:
    """Reuse exact bytes only while their stage-local file identity is stable."""

    def __init__(self) -> None:
        self._values: dict[Path, tuple[_StableFileIdentity, str, bytes]] = {}

    @staticmethod
    def _identity(path: Path) -> _StableFileIdentity:
        if path.is_symlink() or not path.is_file():
            raise PublicTrack1ValidationError("public resource binding is stale") from None
        value = path.stat()
        return _StableFileIdentity(
            value.st_dev, value.st_ino, value.st_size,
            value.st_mtime_ns, value.st_ctime_ns,
        )

    def read(self, path: Path, expected_sha256: str | None = None) -> bytes:
        path = Path(path)
        before = self._identity(path)
        saved = self._values.get(path)
        if saved is not None:
            identity, digest, content = saved
            if before != identity or (
                expected_sha256 is not None and expected_sha256 != digest
            ):
                raise PublicTrack1ValidationError("public resource binding is stale") from None
            return content
        try:
            content = path.read_bytes()
        except OSError:
            raise PublicTrack1ValidationError("public resource binding is stale") from None
        after = self._identity(path)
        digest = _sha256(content)
        if before != after or (
            expected_sha256 is not None and digest != expected_sha256
        ):
            raise PublicTrack1ValidationError("public resource binding is stale") from None
        self._values[path] = (after, digest, content)
        return content


def _authenticated_public_reference_aliases(
    *,
    identity_path: Path,
    reference_manifest_path: Path,
    reference_root: Path,
    reference_contract: ReferenceIdentityContract,
    reference: IndexedFasta,
) -> Mapping[str, str]:
    """Load the exact native/canonical map from authenticated reference data."""

    role = reference_contract.reference_role("public_annotation_benchmark")
    renderer = load_submission_contig_renderer(
        identity_path=identity_path,
        reference_manifest_path=reference_manifest_path,
        reference_root=reference_root,
        reference_binding=ReferenceRoleBinding(
            "public_annotation_benchmark", role.identity_id,
        ),
    )
    aliases = renderer.reference_aliases
    if (
        not aliases
        or set(reference.index) - set(aliases.values())
        or any(target not in reference.index for target in aliases.values())
    ):
        raise PublicTrack1ValidationError(
            "public reference alias binding is invalid"
        )
    return aliases


def _clinical_source_pin(
    contract_path: Path,
    references: Path,
    manifest_path: Path,
    manifest_sha256: str,
) -> PinnedPublicSource:
    """Clinical contracts resolve manifest and output paths from one namespace."""
    from genosapiens.vep_annotation import _local
    _local(contract_path)
    return PinnedPublicSource(
        contract_path, _sha256(contract_path.read_bytes()), references,
        manifest_path, manifest_sha256,
    )


def _public_consumer_descriptors(repository, references):
    from genosapiens.gnomad_publication import load_gnomad_publication_contract
    from genosapiens.gnomad_all_sites import load_gnomad_all_sites_contract
    from genosapiens.vep_annotation import _local
    open_targets_contract = repository / OPEN_TARGETS_CONTRACT_PATH
    open_targets_manifest = references / "manifests/open-targets/open-targets-platform-26.06.0cae4540cfea.json"
    clinvar_contract = repository / "config/public-resources/clinvar-grch38-2026-08-22.json"
    population_contracts = tuple(repository / "config/public-resources" / name for name in (
        "gnomad-joint-projection.json", "gnomad-joint-publication.json",
        "gnomad-chr21-chr22-validation.json", "gnomad-all-sites-v4.1.json",
    ))
    for path in population_contracts:
        _local(path)
    publication = load_gnomad_publication_contract(population_contracts[1])
    all_sites = load_gnomad_all_sites_contract(population_contracts[3])
    return (
        _clinical_source_pin(open_targets_contract, references, open_targets_manifest, OPEN_TARGETS_MANIFEST_SHA256),
        _clinical_source_pin(clinvar_contract, references, references / CLINVAR_MANIFEST_PATH, CLINVAR_MANIFEST_SHA256),
        population.GnomadJointReference(
            *population_contracts[:3], references / publication.manifest_path,
            GNOMAD_JOINT_MANIFEST_SHA256, references / publication.output_directory,
        ),
        population.GnomadAllSitesReference(population_contracts[3], references, GNOMAD_ALL_SITES_MANIFEST_SHA256),
    )


def _preflight_public_consumers(references, clinical_open_targets, clinvar, joint, all_sites):
    """Reject unusable current metadata before authenticating whole resources.

    Consumers reuse pinned publication reconciliation with metadata and live
    query checks. An explicit strict verifier still performs whole-file replay.
    """
    from contextlib import ExitStack
    from genosapiens import admission_links as links, gnomad_all_sites as sites
    from genosapiens import gnomad_publication as publication, gnomad_validation as validation
    from genosapiens.gnomad_projection import load_gnomad_projection_contract
    from genosapiens.manifest import PublicResourceManifest
    from genosapiens.vep_annotation import _local
    stage = "resource-topology"
    try:
        references = _local(references)
        for pin in (clinical_open_targets, clinvar):
            if _local(pin.root) != references or not _local(pin.manifest_path).is_relative_to(references):
                raise ValueError
        if (_local(all_sites.reference_root) != references
                or not _local(joint.manifest_path).is_relative_to(references)
                or not _local(joint.release_directory).is_relative_to(references)):
            raise ValueError
        stage = "manifest-qualification"
        links.preflight_clinical_source(clinical_open_targets, "open_targets")
        links.preflight_clinical_source(clinvar, "clinvar")
        with ExitStack() as stack:
            files = links._Files(stack)
            projection = load_gnomad_projection_contract(_local(joint.projection_contract_path))
            contract = publication.load_gnomad_publication_contract(_local(joint.publication_contract_path))
            policy = validation.load_gnomad_validation_contract(_local(joint.validation_contract_path))
            population.validate_joint_contract_bindings(
                projection, contract, policy, _sha256(joint.projection_contract_path.read_bytes()),
            )
            manifest = PublicResourceManifest.from_mapping(files.json(joint.manifest_path, joint.manifest_sha256))
            publication._verify_manifest_contract(
                manifest, contract, projection, publication._transformation_parameters_sha256(
                    joint.publication_contract_path, joint.projection_contract_path, joint.validation_contract_path, policy,
                ),
            )
            if (_local(joint.release_directory) != references / contract.output_directory
                    or sum(o.size_bytes for o in manifest.outputs) > contract.aggregate_output_cap_bytes):
                raise ValueError
            for output in manifest.outputs:
                path = _local(references / output.path)
                if not path.is_relative_to(references) or not path.is_file() or path.stat().st_size != output.size_bytes:
                    raise ValueError
            all_contract = sites.load_gnomad_all_sites_contract(_local(all_sites.contract_path))
            all_sha = _sha256(all_sites.contract_path.read_bytes())
            bundle = files.json(references / all_contract.bundle_manifest_path, all_sites.bundle_manifest_sha256)
            entries = sites.validate_bundle_metadata(bundle, all_contract, all_sha)
            output_bytes = 0
            for source, (entry, _transform) in zip(all_contract.sources, entries, strict=True):
                current = PublicResourceManifest.from_mapping(files.json(
                    references / entry["source_manifest_path"], entry["source_manifest_sha256"],
                ))
                sites._verify_source_manifest_contract(current, all_contract, source, all_sha)
                for output in current.outputs:
                    path = _local(references / output.path)
                    if not path.is_relative_to(references) or not path.is_file() or path.stat().st_size != output.size_bytes:
                        raise ValueError
                    output_bytes += output.size_bytes
            if output_bytes > all_contract.maximum_derived_bytes:
                raise ValueError
    except Exception:
        raise PublicTrack1ValidationError("public consumer resource preflight failed", preparation_stage=stage) from None


def _load_public_resources(
    repository_root: str | Path,
    reference_root: str | Path,
    scratch_root: str | Path,
    *,
    materialize_analysis: bool,
    analysis_hpo_loader=None,
) -> PublicResourceContext | PublicPreparationResourceContext:
    """Authenticate one stage's resources through shared verified boundaries."""
    stage = "environment"
    try:
        if any(name != "lo" for _index, name in socket.if_nameindex()):
            raise ValueError
        repository = Path(repository_root).resolve(strict=True)
        references = Path(reference_root).resolve(strict=True)
        scratch = Path(scratch_root).resolve(strict=True)
        if not scratch.is_relative_to(references / "scratch/public-build"):
            raise ValueError
        stable_bytes = _StageBoundBytes()
        stage = "resource-topology"
        model_loci = None
        if materialize_analysis:
            model_loci = load_public_model_loci(repository, references)
            _authenticated_policy_loci(repository, model_loci)
        clinical_open_targets, clinvar, population_joint, population_all_sites = _public_consumer_descriptors(repository, references)
        _preflight_public_consumers(references, clinical_open_targets, clinvar, population_joint, population_all_sites)
        stage = "reference"
        identity_path = repository / REFERENCE_IDENTITY_PATH
        identity_bytes = stable_bytes.read(identity_path)
        identity_sha256 = _sha256(identity_bytes)
        contract = load_reference_identity(identity_path)
        reference_manifest = references / (
            "manifests/grch38/"
            "challenge-grch38-masked-v2-plus-hs38d1.a92c50f61c1b.json"
        )
        reference_manifest_bytes = stable_bytes.read(
            reference_manifest,
            REFERENCE_MANIFEST_SHA256,
        )
        with open_installed_reference(
            contract,
            identity_path=identity_path,
            reference_manifest_path=reference_manifest,
            reference_root=references,
            reuse_published_validation=True,
        ) as (manifest_sha256, _fetch, reference_contigs):
            if (
                manifest_sha256 != REFERENCE_MANIFEST_SHA256
                or len(reference_contigs) != contract.sequence_count
            ):
                raise ValueError
        public_role = contract.reference_role("public_annotation_benchmark")
        public_fasta = references / contract.expected_output(
            public_role.fasta_output_role
        ).path
        public_fai = references / contract.expected_output(
            public_role.fai_output_role
        ).path
        reference = IndexedFasta(public_fasta, public_fai)
        source_to_reference_contigs = _authenticated_public_reference_aliases(
            identity_path=identity_path,
            reference_manifest_path=reference_manifest,
            reference_root=references,
            reference_contract=contract,
            reference=reference,
        )
        construction = load_construction_report(
            references / contract.expected_output("construction_report").path
        )
        if (
            not construction.accepted
            or construction.public_annotation_benchmark_identity_id
            != public_role.identity_id
            or construction.total_records != len(reference.index)
        ):
            raise ValueError

        stage = "phenotype-and-disease"
        hpo_path = references / (
            "derived/phenotype/hpo-ic/v1/"
            "060cf7c49671f9a98f9a8578c35f0fb8eeae62166d79cb40c2bcc81855abcc89/"
            "store.json"
        )
        hpo_loader = (load_completed_hpo_ic if not materialize_analysis or analysis_hpo_loader is None
            else analysis_hpo_loader)
        hpo = hpo_loader(reference_root=references,
            manifest_path=references / 'manifests/phenotype/hpo-ic-060cf7c49671f9a98f9a8578c35f0fb8eeae62166d79cb40c2bcc81855abcc89.json',
            expected_manifest_sha256=HPO_MANIFEST_SHA256,expected_store_sha256=HPO_STORE_SHA256)
        rare_contract = repository / RARE_DISEASE_CONTRACT_PATH
        ontology_contract = repository / ONTOLOGY_CONTRACT_PATH
        sensitivity_resources = None
        profiles = None
        genes = None
        if materialize_analysis:
            sensitivity_resources = getattr(hpo_loader, 'bind_sensitivity',
                _authenticate_loaded_ic_store)(
                hpo,
                expected_store_sha256=HPO_STORE_SHA256,
            )
            profiles = read_pinned_disease_profiles(
                contract_path=rare_contract,
                ontology_contract_path=ontology_contract,
                reference_root=references,
                expected_bundle_sha256=RARE_DISEASE_BUNDLE_SHA256,
            )
            genes = read_pinned_disease_gene_evidence(
                contract_path=rare_contract,
                reference_root=references,
                expected_bundle_sha256=RARE_DISEASE_BUNDLE_SHA256,
            )

        stage = "repeat"
        repeat_contract = repository / REPEAT_CONTRACT_PATH
        repeat_manifest = references / (
            "manifests/variant-class-resources/"
            "broad-str-analysis-disease-catalog-grch38-d41c85c.json"
        )
        repeat_catalog_path = references / (
            "releases/repeats/broad-str-analysis/"
            "d41c85c97da3cd8c8ac9d7d47883570c95069874/"
            "variant_catalog_without_offtargets.GRCh38.json"
        )
        repeat_contract_bytes = stable_bytes.read(repeat_contract)
        repeat_manifest_bytes = stable_bytes.read(
            repeat_manifest,
            REPEAT_MANIFEST_SHA256,
        )
        aliases_path = references / contract.expected_output("contig_aliases").path
        repeat_catalog = None
        repeat_reference = None
        if materialize_analysis:
            repeat_catalog = load_pinned_repeat_catalog(
                stable_bytes.read(repeat_catalog_path),
                contract_bytes=repeat_contract_bytes,
                manifest_bytes=repeat_manifest_bytes,
                expected_contract_sha256=_sha256(repeat_contract_bytes),
                expected_manifest_sha256=REPEAT_MANIFEST_SHA256,
                evidence_id="track1-public-repeat-catalog",
            )
            repeat_reference = RepeatReferenceInput(
                identity_bytes,
                reference_manifest_bytes,
                stable_bytes.read(aliases_path),
                identity_sha256,
                REFERENCE_MANIFEST_SHA256,
            )

        stage = "annotation"
        vep_contract = repository / VEP_CONTRACT_PATH
        vep_bundle = references / (
            "manifests/ensembl/vep/"
            "ensembl-vep-116-grch38.f1c47e44cd07.bundle.json"
        )
        open_targets_contract = repository / OPEN_TARGETS_CONTRACT_PATH
        open_targets_manifest = references / (
            "manifests/open-targets/"
            "open-targets-platform-26.06.0cae4540cfea.json"
        )
        constraint = None
        compatibility = None
        submission_reference = None
        if materialize_analysis:
            constraint = OpenTargetsConstraintReference(
                open_targets_contract,
                references / "releases/open-targets/26.06",
                open_targets_manifest,
                OPEN_TARGETS_MANIFEST_SHA256,
            )
            compatibility = VepReferenceCompatibility(
                identity_path,
                identity_sha256,
                references,
                reference_manifest,
                REFERENCE_MANIFEST_SHA256,
            )
            submission_reference = PinnedPublicSource(
                identity_path,
                identity_sha256,
                references,
                reference_manifest,
                REFERENCE_MANIFEST_SHA256,
            )
        stage = "clinical-and-population-descriptors"
        clinvar_contract, clinvar_manifest = clinvar.contract_path, clinvar.manifest_path
        joint_manifest = population_joint.manifest_path
        all_sites_contract = json.loads(stable_bytes.read(population_all_sites.contract_path))
        all_sites_manifest = references / all_sites_contract["bundle_manifest_path"]
        population_contracts = (
            population_joint.projection_contract_path, population_joint.publication_contract_path,
            population_joint.validation_contract_path, population_all_sites.contract_path,
        )
        population_reference = population.PopulationReferenceCompatibility(
            identity_path, identity_sha256, references, reference_manifest, REFERENCE_MANIFEST_SHA256,
        )

        stage = "evaluator"
        evaluator_source = references / (
            "releases/challenge-evaluator/"
            "SageBio-rare-disease-real-kid-mva-hackathon-2026/"
            "1c761cc23d90aebe6a011fd5b0b99517df42408c/evaluation.py"
        )
        stable_bytes.read(evaluator_source, EVALUATOR_SOURCE_SHA256)
        if materialize_analysis:
            load_pinned_evaluator(
                evaluator_source,
                pin_path=repository / "config/track1/challenge-evaluator.json",
            )

        stage = "descriptors"
        from genosapiens.installed_public_index import load_installed_index_pin
        index_pin = load_installed_index_pin(repository) or dict(
            manifest_path=COMMON_ALLELE_PUBLIC_INDEX_MANIFEST_PATH,
            manifest_sha256=COMMON_ALLELE_PUBLIC_INDEX_MANIFEST_SHA256)
        descriptor_paths = (
            ("runner-model-loci-contract", repository / MODEL_LOCUS_CONTRACT_PATH, MODEL_LOCUS_DESCRIPTOR_SHA256),
            ("runner-model-loci-manifest", references / MODEL_LOCUS_MANIFEST_PATH, MODEL_LOCUS_MANIFEST_SHA256),
            ("runner-reference-contract", identity_path, identity_sha256),
            ("runner-reference-manifest", reference_manifest, REFERENCE_MANIFEST_SHA256),
            ("runner-vep-contract", vep_contract, _sha256(stable_bytes.read(vep_contract))),
            ("runner-vep-manifest", vep_bundle, VEP_BUNDLE_SHA256),
            (
                "runner-open-targets-contract",
                open_targets_contract,
                _sha256(stable_bytes.read(open_targets_contract)),
            ),
            (
                "runner-open-targets-manifest",
                open_targets_manifest,
                OPEN_TARGETS_MANIFEST_SHA256,
            ),
            ("runner-rare-disease-contract", rare_contract, _sha256(stable_bytes.read(rare_contract))),
            (
                "runner-rare-disease-manifest",
                references / (
                    "manifests/rare-disease-evidence/"
                    "rare-disease-evidence-2026-08-30-v6.cded6e513755.bundle.json"
                ),
                RARE_DISEASE_BUNDLE_SHA256,
            ),
            (
                "runner-hpo-manifest",
                references / (
                    "manifests/phenotype/"
                    "hpo-ic-060cf7c49671f9a98f9a8578c35f0fb8eeae62166d79cb40c2bcc81855abcc89.json"
                ),
                HPO_MANIFEST_SHA256,
            ),
            ("runner-repeat-contract", repeat_contract, _sha256(repeat_contract_bytes)),
            ("runner-repeat-manifest", repeat_manifest, REPEAT_MANIFEST_SHA256),
            ("runner-clinvar-contract", clinvar_contract, _sha256(stable_bytes.read(clinvar_contract))),
            ("runner-clinvar-manifest", clinvar_manifest, CLINVAR_MANIFEST_SHA256),
            ("runner-gnomad-joint-manifest", joint_manifest, GNOMAD_JOINT_MANIFEST_SHA256),
            ("runner-gnomad-all-sites-manifest", all_sites_manifest, GNOMAD_ALL_SITES_MANIFEST_SHA256),
            ("runner-common-allele-exceptions", references / COMMON_ALLELE_EXCEPTION_MANIFEST_PATH, COMMON_ALLELE_EXCEPTION_MANIFEST_SHA256),
            ("runner-common-allele-public-index", references / index_pin['manifest_path'], index_pin['manifest_sha256']),
            ("runner-gnomad-precision", references / GNOMAD_PRECISION_REPORT_PATH, GNOMAD_PRECISION_REPORT_SHA256),
            *(("runner-population-contract", path, _sha256(stable_bytes.read(path))) for path in population_contracts),
            ("runner-evaluator-source", evaluator_source, EVALUATOR_SOURCE_SHA256),
        )
        contents = tuple(stable_bytes.read(path, digest) for _kind, path, digest in descriptor_paths)
        runner_resources = tuple(
            RankerResource(kind, _sha256(content))
            for (kind, _path, _digest), content in zip(descriptor_paths, contents, strict=True)
        )
        if not materialize_analysis:
            return PublicPreparationResourceContext(
                repository,
                references,
                scratch,
                contract,
                reference_manifest,
                reference,
                source_to_reference_contigs,
                hpo,
                clinical_open_targets,
                clinvar,
                population_joint,
                population_all_sites,
                population_reference,
                runner_resources,
            )
        from genosapiens.wgs_profile_projection import load_installed_profile_resources
        stage='profile-resources'
        profile_catalog,regulatory_index,regulatory_sources=load_installed_profile_resources(
            repository,references,gene_store=genes,reference=submission_reference)
        return PublicResourceContext(
            repository,
            references,
            scratch,
            contract,
            reference_manifest,
            reference,
            source_to_reference_contigs,
            construction,
            hpo,
            profiles,
            genes,
            repeat_catalog,
            repeat_reference,
            constraint,
            compatibility,
            submission_reference,
            evaluator_source,
            contents,
            runner_resources,
            sensitivity_resources,
            clinical_open_targets, clinvar, population_joint, population_all_sites, population_reference,
            model_loci,
            profile_catalog=profile_catalog,regulatory_index=regulatory_index,
            regulatory_sources=regulatory_sources,
        )
    except PublicTrack1ValidationError:
        raise
    except Exception:
        from genosapiens.track1_validation_runner import PREPARATION_STAGES
        raise PublicTrack1ValidationError(
            f"public resource context is invalid at {stage} stage",
            preparation_stage=stage if stage in PREPARATION_STAGES else "preparation",
        ) from None


def load_public_preparation_resources(
    repository_root: str | Path,
    reference_root: str | Path,
    scratch_root: str | Path,
) -> PublicPreparationResourceContext:
    """Load the authenticated live objects required by public preparation."""
    result = _load_public_resources(
        repository_root, reference_root, scratch_root,
        materialize_analysis=False,
    )
    if type(result) is not PublicPreparationResourceContext:
        raise AssertionError("invalid preparation resource loader result")
    return result


def load_public_resource_context(
    repository_root: str | Path,
    reference_root: str | Path,
    scratch_root: str | Path,
    *, analysis_hpo_loader=None,
) -> PublicResourceContext:
    """Load the complete authenticated resources required by public analysis."""
    result = _load_public_resources(
        repository_root, reference_root, scratch_root,
        materialize_analysis=True,
        analysis_hpo_loader=analysis_hpo_loader,
    )
    if type(result) is not PublicResourceContext:
        raise AssertionError("invalid analysis resource loader result")
    return result


@dataclass(frozen=True, slots=True, repr=False)
class MaterializedRunnerView:
    input_paths: tuple[Path, ...]
    input_sha256s: tuple[str, ...]
    input_sizes: tuple[int, ...]
    resource_bytes: tuple[bytes, ...]

    def __repr__(self) -> str:
        return "<materialized public runner view>"


def _content_index(paths: tuple[Path, ...]) -> Mapping[str, tuple[Path, int]]:
    import os
    import stat
    index: dict[str, tuple[Path, int]] = {}
    for path in paths:
        if any(parent.is_symlink() for parent in (path, *path.parents)) or not path.is_file():
            continue
        with path.open("rb") as stream:
            before = os.fstat(stream.fileno())
            if not stat.S_ISREG(before.st_mode):
                raise PublicTrack1ValidationError("runner artifact is not a regular file")
            digest = hashlib.sha256()
            for block in iter(lambda: stream.read(65536), b""):
                digest.update(block)
            after = os.fstat(stream.fileno())
        identity = lambda value: (value.st_dev, value.st_ino, value.st_size, value.st_mtime_ns, value.st_ctime_ns)
        if identity(before) != identity(after):
            raise PublicTrack1ValidationError("runner artifact changed during authentication")
        previous = index.get(digest.hexdigest())
        if previous is not None and previous[1] != before.st_size:
            raise PublicTrack1ValidationError("runner artifact digest collision")
        index[digest.hexdigest()] = (path, before.st_size)
    return index


def materialize_runner_view(
    view: RankerCaseView,
    resources: PublicResourceContext,
) -> MaterializedRunnerView:
    """Resolve only hashes in a truth-free view from closed public fixture roots."""
    try:
        if type(view) is not RankerCaseView or type(resources) is not PublicResourceContext:
            raise ValueError
        roots = (
            resources.repository_root / "tests/fixtures/track1/roster",
            resources.repository_root / "tests/fixtures/track1/materialized",
        )
        input_paths = tuple(
            path
            for root in roots
            for path in root.rglob("*")
            if path.name.endswith(".vcf") or path.name.endswith(".hpo.txt")
        )
        descriptor_paths = (
            *tuple(sorted(
                path
                for path in (
                    resources.repository_root
                    / "tests/fixtures/track1/roster/source-formats"
                ).glob("*.json")
                if path.is_file() and not path.is_symlink()
            )),
            *tuple(sorted(
                path
                for path in (
                    resources.repository_root
                    / "tests/fixtures/track1/materialized"
                ).glob("*.source-format.json")
                if path.is_file() and not path.is_symlink()
            )),
            resources.repository_root
            / "tests/fixtures/track1/materialized/hg002-public-background.json",
        )
        inputs = _content_index(input_paths)
        descriptors = _content_index(descriptor_paths)
        resolved_paths = []
        resolved_sizes = []
        for item in view.inputs:
            path, size = inputs[item.sha256]
            if (
                item.kind is InputKind.HPO_TEXT
                and not path.name.endswith(".hpo.txt")
                or item.kind is InputKind.VCF
                and not path.name.endswith(".vcf")
            ):
                raise ValueError
            resolved_paths.append(path)
            resolved_sizes.append(size)
        common = resources.common_resource_bytes
        common_hashes = tuple(_sha256(value) for value in common)
        requested = tuple(item.descriptor_sha256 for item in view.resources)
        if requested[:len(common)] != common_hashes:
            raise ValueError
        from genosapiens.bounded_json import read_bounded_bytes
        remaining = tuple(read_bounded_bytes(descriptors[value][0], maximum_bytes=2 * 1024 * 1024)
                          for value in requested[len(common):])
        if tuple(_sha256(value) for value in remaining) != requested[len(common):]:
            raise ValueError
        return MaterializedRunnerView(
            tuple(resolved_paths),
            tuple(item.sha256 for item in view.inputs),
            tuple(resolved_sizes),
            (*common, *remaining),
        )
    except Exception:
        raise PublicTrack1ValidationError("runner view cannot be materialized") from None


def _prepare_call_input(
    path: Path,
    source_format_bytes: bytes,
    resources: PublicResourceContext,
) -> tuple[SuppliedCallInputs, SourceDiagnosticEvidence]:
    """Finite upstream fixture probe for the shared external VCF adapter.

    Case execution uses prepare_public_case and source ledgers. This helper is
    restricted to the roster's existing 2 MiB plain-VCF fixture contract; its
    callers qualify invented source fixtures before any ranking experiment.
    """
    from genosapiens.bounded_json import read_bounded_bytes
    path = Path(path).absolute()
    if (path.suffix != ".vcf" or not any(path.parts[index:index + 3] == ("tests", "fixtures", "track1")
            for index in range(len(path.parts) - 2))):
        raise PublicTrack1ValidationError("fixture VCF probe is outside its source contract")
    content = read_bounded_bytes(path, maximum_bytes=2 * 1024 * 1024, reference_root=path.parent)
    if type(source_format_bytes) is not bytes or len(source_format_bytes) > 64 * 1024:
        raise PublicTrack1ValidationError("fixture format exceeds its source contract")
    report = evaluate_public_validation_vcf_concordance(
        path,
        resources.reference.fetch,
        {name: item.length for name, item in resources.reference.index.items()},
        resources.source_to_reference_contigs,
        contract=resources.reference_contract,
        reference_manifest_sha256=REFERENCE_MANIFEST_SHA256,
    )
    qc = summarize_vcf_genotype_qc(path, concordance_report=report)
    normalized = normalize_public_validation_vcf_small_variants(
        path,
        concordance_report=report,
        fetch_reference=resources.reference.fetch,
        source_to_reference_contigs=resources.source_to_reference_contigs,
    )
    inventory = build_variant_class_inventory(qc, allele_balance_screen=True)
    source_format_sha256 = _sha256(source_format_bytes)
    supplied = SuppliedCallInputs(
        content,
        qc.input_sha256,
        qc,
        normalized,
        inventory,
        source_format_bytes=source_format_bytes,
        source_format_sha256=source_format_sha256,
    )
    diagnostic = SourceDiagnosticEvidence(
        supplied.input_sha256,
        qc.sample_names[0],
        report,
        qc,
        inventory,
    )
    return supplied, diagnostic


def _annotation_artifacts(
    requests: tuple[VepRequest, ...],
    *,
    resources: PublicResourceContext,
    publication_id: str,
) -> Mapping[str, bytes]:
    """Bounded external-tool artifacts for upstream fixture qualification only.

    No prepared-case, execution, handoff or saved-worker reader calls this
    adapter. Its flat VEP response is an internal shared annotation-kernel
    input and cannot authorize a current scientific execution.
    """
    from genosapiens.bounded_json import read_bounded_bytes
    from genosapiens.bounded_wgs_config import load_bounded_wgs_execution_config
    config_content = read_bounded_bytes(resources.repository_root / "config/track1/bounded-wgs-execution-v5.json")
    config = load_bounded_wgs_execution_config(config_content, _sha256(config_content))
    if (type(requests) is not tuple or len(requests) > config.annotation.max_partition_requests
            or type(publication_id) is not str or not re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9._-]*", publication_id)):
        raise PublicTrack1ValidationError("fixture annotation exceeds its source contract")
    destination = resources.scratch_root / publication_id
    if not destination.exists():
        run_vep_annotation(
            requests,
            contract_path=resources.repository_root / VEP_CONTRACT_PATH,
            reference_root=resources.reference_root,
            bundle_manifest_sha256=VEP_BUNDLE_SHA256,
            runtime_image=RUNTIME_IMAGE,
            upstream_image_digest=(
                "sha256:f22d1f89b28c6c4caf043ddb2269fdfff8940364da589cadbbabc6212a436949"
            ),
            output_directory=destination,
            constraint_reference=resources.vep_constraint,
            reference_compatibility=resources.vep_compatibility,
        )
    if destination.is_symlink() or not destination.is_dir():
        raise PublicTrack1ValidationError("annotation publication is invalid")
    artifacts = {}
    remaining = config.annotation.max_shard_bytes
    for path in sorted(destination.iterdir()):
        if path.is_symlink() or not path.is_file():
            raise PublicTrack1ValidationError("fixture annotation publication is invalid")
        content = read_bounded_bytes(path, maximum_bytes=remaining, reference_root=destination)
        remaining -= len(content)
        artifacts[path.name] = content
    vep_data = json.loads(read_bounded_bytes(resources.repository_root / VEP_CONTRACT_PATH))
    artifacts["reference.fai"] = read_bounded_bytes(
        resources.reference_root / (vep_data["fasta"]["output_path"] + ".fai"),
        maximum_bytes=remaining, reference_root=resources.reference_root)
    remaining -= len(artifacts["reference.fai"])
    binding = json.loads(artifacts["binding.json"])
    if "native_annotation" in binding:
        artifacts["cache.inventory.json"] = read_bounded_bytes(
            resources.reference_root
            / vep_data["cache"]["output_root"]
            / "cache.inventory.json", maximum_bytes=remaining, reference_root=resources.reference_root)
    return artifacts


def _subset_gene_store(
    full: DiseaseGeneEvidenceStore,
    gene_ids: frozenset[str],
    condition_ids: frozenset[str] = frozenset(),
) -> tuple[DiseaseGeneEvidenceStore | None, frozenset[str]]:
    if not gene_ids and not condition_ids:
        return None, frozenset()
    # VEP emits versioned Ensembl gene accessions for this pinned cache while
    # the public disease-gene publication records the stable accession.  The
    # version suffix is annotation metadata, not part of the Ensembl gene
    # identity used by the disease-gene join.
    stable_gene_ids = frozenset(value.split(".", 1)[0] for value in gene_ids)
    rows = {
        table: [
            record.source_fields
            for record in full.records.values()
            if record.table == table
        ]
        for table in PUBLIC_TABLES
    }
    # Resource closure includes the exact alias and subtype proof records.
    # This broad connected subset is only a resource selection; the admission
    # kernel still applies directed reachability and source model requirements.
    from collections import defaultdict,deque
    neighbours=defaultdict(set)
    for relation in full.relations:
        if relation.subject_id and relation.object_id:
            neighbours[relation.subject_id].add(relation.object_id)
            neighbours[relation.object_id].add(relation.subject_id)
    def closure(initial):
        reached=set(initial);pending=deque(reached)
        while pending:
            for value in neighbours[pending.popleft()]-reached:
                reached.add(value);pending.append(value)
        return reached
    related=closure(condition_ids)
    selected_gene_rows = [
        row for row in rows["gene_assertion"]
        if any(
            type(value) is str and value.split(".", 1)[0] in stable_gene_ids
            for value in (row.get("ensembl_gene_id"), row.get("gene_id"))
        ) or row.get('source_disease_id') in related or row.get('mapped_disease_id') in related
    ]
    diseases = closure(related | {
        value
        for row in selected_gene_rows
        for value in (row.get("source_disease_id"), row.get("mapped_disease_id"))
        if type(value) is str and value
    })
    relations = rows["disease_relation"]
    selected_relations = [
        row for row in relations
        if row.get("root_disease_id") in diseases
        and row.get("target_disease_id") in diseases
    ]
    store = build_disease_gene_evidence(
        selected_gene_rows,
        source_binding=full.source_binding,
        inheritance_rows=[
            row for row in rows["disease_inheritance"]
            if row.get("disease_id") in diseases
        ],
        disease_rows=[
            row for row in rows["disease"]
            if row.get("disease_id") in diseases
        ],
        alias_rows=[
            row for row in rows["disease_alias"]
            if row.get("disease_id") in diseases
        ],
        relation_rows=selected_relations,
    )
    return store, frozenset(diseases)


def _subset_profile_store(
    full: DiseaseProfileStore,
    diseases: frozenset[str],
    hpo: HpoICStore,
) -> DiseaseProfileStore | None:
    if not diseases:
        return None
    rows = [
        item.source_fields
        for item in full.annotations
        if item.source_fields.get("source_disease_id") in diseases
        or item.source_fields.get("mapped_disease_id") in diseases
    ]
    return build_disease_profiles(
        hpo.graph,
        rows,
        source_binding=full.source_binding,
    )


@dataclass(slots=True, repr=False)
class PreparedPublicCase:
    case_id: str
    view: RankerCaseView
    materialized: MaterializedRunnerView
    public_case: object
    clinical_indexes: object
    preparation_wall_seconds: str
    preparation_cpu_seconds: str
    screening_arguments: dict
    annotation_source: object | None = None
    invocation_count: int = 0

    def __repr__(self) -> str:
        return "<prepared public Track 1 case>"


def _population_template(policy_bytes: bytes):
    fixed = json.loads(policy_bytes)["fixed_settings"]
    template = fixed["population"]
    if set(template) != {"frequency_rule", "transform", "rescues", "derivation_sha256"} or SHA256.fullmatch(template["derivation_sha256"]) is None:
        raise ValueError
    rule = dict(template["frequency_rule"])
    for name in ("allowed_filters", "region_flags", "unreliable_signals"):
        rule[name] = tuple(rule[name])
    rescues = dict(template["rescues"])
    for name in ("pathogenic_significances", "pathogenic_review_statuses", "technical_signals"):
        rescues[name] = tuple(rescues[name])
    return frequency.FrequencyRule(**rule), frequency.PenaltyTransform(**template["transform"]), frequency.RescueRules(**rescues)


def preflight_public_population_policy(repository_root: str | Path) -> None:
    """Reject impossible source requirements before launching case workers.

    Historical policies remain loadable for inspection. This execution gate is
    additional to manifest authentication and actual per-record quality checks.
    """
    from genosapiens.gnomad_projection import load_gnomad_projection_contract
    root = Path(repository_root)
    rule, _, _ = _population_template((root / POLICY_PATH).read_bytes())
    projection = load_gnomad_projection_contract(
        root / "config/public-resources/gnomad-joint-projection.json")
    if "filters_unknown" in rule.unreliable_signals and not projection.dataset(rule.subset).has_filters:
        raise PublicTrack1ValidationError(
            "population policy requires source filters absent from the selected dataset",
            preparation_stage="population-policy-binding",
        )


def authenticate_population_development_derivation(experiment, policies, repository_root, *, recovery_resolver):
    """Check source-selection arithmetic against actual completed attempts.

    The recorded source observations are independently reconstructed during
    checkpoint recovery. This check cannot turn an operational pilot, caller
    declaration or historical result into development evidence.
    """
    fixed=json.loads(policies.policy_bytes)["fixed_settings"]
    content=_bound_bytes(Path(repository_root)/"config/track1/public-population-development-derivation.json",fixed["population"]["derivation_sha256"])
    derivation=json.loads(content)
    low_keys=tuple(derivation["selected_low_identity_sha256s"])
    high_key=derivation["selected_high_identity_sha256"]
    for candidate in experiment.candidates:
        observed={key:set() for key in (*low_keys,high_key)}
        for attempt in experiment.attempts:
            if attempt.split!="development" or attempt.policy_id!=candidate.policy_id:continue
            for key,value,an in recovery_resolver.frequency_observations(attempt, (*low_keys, high_key)):
                if key in observed:observed[key].add((Fraction(value),an))
        if any(len(values)!=1 for values in observed.values()):
            raise PublicTrack1ValidationError("population development observations are incomplete or inconsistent")
        rows={key:next(iter(values)) for key,values in observed.items()}
        if (
            max(rows[key][0] for key in low_keys)!=Fraction(derivation["development_selected_low_frequency_upper_bound"])
            or rows[high_key][0]!=Fraction(derivation["development_selected_high_frequency"])
            or min(value[1] for value in rows.values())!=derivation["minimum_observed_development_an"]
        ):
            raise PublicTrack1ValidationError("population development derivation does not match source replay")



def prepare_public_case(case_id, view, *, resources, policies, output_directory=None) -> PreparedPublicCase:
    """Stream all supplied sources into the single current prepared contract."""
    from genosapiens.bounded_wgs_config import load_bounded_wgs_execution_config
    from genosapiens.case_manifest import (CaseInput, CaseFileClass, ConfigurationBinding,
        ReferenceBundleBinding, _input_id, create_case_manifest, case_manifest_as_bytes)
    from genosapiens.common_allele_filter import load_common_allele_filter_policy
    from genosapiens.common_allele_public_indexes import verify_common_allele_public_indexes
    from genosapiens.gnomad_faf95_adapter import (load_gnomad_faf95_precision_contract,
        verify_gnomad_faf95_precision_qualification)
    from genosapiens.population_clinical_session import PopulationClinicalResourceSession
    from genosapiens.track1_public_validation_coordinator import executable_stack_bytes, publish_validation_file
    from genosapiens.vcf_qc import iter_vcf_qc_records
    from genosapiens.wgs_admission import NativeSourceContract, ReferenceRole
    from genosapiens.wgs_ledger import canonical, publish_bytes, load_model, LedgerManifest
    from genosapiens.wgs_preparation import (
        PreparationRunOptions, prepare_bounded_case, source_worker_reference,
    )
    from genosapiens.wgs_public_case import publish_public_case
    from genosapiens.wgs_public_filter import PublicEarlyFilterDecisionProvider, PublicFilterWorkerConfig
    from genosapiens.wgs_read_rescue import load_prepared_intake_source, merge_prepared_intake, MergedReadIntake
    from genosapiens.wgs_shards import read_bounded_metadata
    from genosapiens.wgs_validation_authority import PublicValidationAuthority, PublicValidationInput

    wall_start, cpu_start = time.monotonic_ns(), time.process_time_ns()
    children_start = resource.getrusage(resource.RUSAGE_CHILDREN)
    stage = "manifest-qualification"
    try:
        loci = _authenticated_policy_loci(resources.repository_root, resources.model_loci)
        if policies.hypothesis_policy.loci != loci:
            raise ValueError
        stage = "materialization"
        materialized = materialize_runner_view(view, resources)
        stage = "input-descriptors"
        formats = {}
        for item, content in zip(view.resources, materialized.resource_bytes, strict=True):
            if item.kind == "runner-source-format":
                descriptor = json.loads(content)
                if (set(descriptor) != {"schema", "input_sha256", "reference_role", "reference_identity_id",
                        "reference_manifest_sha256", "mitochondrial_contigs", "repeat_identity"}
                        or descriptor["schema"] != "genosapiens.supplied-vcf-format"
                        or descriptor["input_sha256"] in formats):
                    raise ValueError
                formats[descriptor["input_sha256"]] = descriptor
        items = tuple(PublicValidationInput(item.kind.value, str(path), item.sha256, size)
            for item, path, size in zip(view.inputs, materialized.input_paths, materialized.input_sizes, strict=True))
        vcfs = tuple(sorted((item for item in items if item.kind == "vcf"), key=lambda item: item.sha256))
        hpo, = (item for item in items if item.kind == "hpo_text")
        if not vcfs or set(formats) != {item.sha256 for item in vcfs}:
            raise ValueError
        reference = resources.reference_contract.reference_role("public_annotation_benchmark")
        for row in formats.values():
            if (row["reference_role"] != "public_annotation_benchmark"
                    or row["reference_identity_id"] != reference.identity_id
                    or row["reference_manifest_sha256"] != REFERENCE_MANIFEST_SHA256):
                raise ValueError
        stack = resources.executable_stack_content
        if stack is None:
            stack = executable_stack_bytes(resources.repository_root)
        # A thin worker must receive and verify a source/config projection of
        # the controller's full stack; it must never fabricate a smaller stack.
        from genosapiens.track1_public_worker import verify_scientific_stack_projection
        verify_scientific_stack_projection(resources.repository_root, stack)
        stack_sha = _sha256(stack)
        bindings = policies.runtimes[0].candidate.bindings
        config_content = _bound_bytes(resources.repository_root / "config/track1/bounded-wgs-execution-v5.json",
                                     bindings.bounded_wgs_execution_sha256)
        config = load_bounded_wgs_execution_config(config_content, bindings.bounded_wgs_execution_sha256)
        output_root = resources.public_output_root or resources.scratch_root / "scratch/public-build/bounded-public"
        authority = PublicValidationAuthority(case_id, _sha256(ranker_case_view_as_json(view)), view, items,
            tuple(item.descriptor_sha256 for item in view.resources), policies.policy_sha256,
            tuple(sorted((runtime.candidate.bindings for runtime in policies.runtimes), key=lambda row: row.sha256)),
            reference.identity_id, REFERENCE_MANIFEST_SHA256,
            tuple(sorted({*(item.descriptor_sha256 for item in view.resources),
                GNOMAD_JOINT_MANIFEST_SHA256, CLINVAR_MANIFEST_SHA256, HPO_MANIFEST_SHA256,
                bindings.common_allele_exception_manifest_sha256, bindings.common_allele_public_index_manifest_sha256,
                bindings.gnomad_projection_precision_sha256})),
            stack_sha, RUNTIME_IMAGE.split("@", 1)[1], str(output_root),
            PreparationRunOptions(False).sha256)
        output = Path(output_directory) if output_directory is not None else Path(output_root) / ("case-" + authority.sha256[:24])
        authority.require_execution(original_input_path=Path(vcfs[0].path), hpo_input_sha256=hpo.sha256,
            output_path=output, executable_stack_sha256=stack_sha, runtime_image_digest=authority.runtime_image_digest)
        if output.exists():
            raise ValueError
        output.mkdir(mode=0o700, parents=True)
        publish_bytes(output, "public-authority.json", authority.as_bytes())
        private_hpo = output / "phenotypes.hpo.txt"
        publish_validation_file(private_hpo, Path(hpo.path), expected_sha256=hpo.sha256, expected_size_bytes=hpo.size_bytes)
        # HPO ingestion requires a private input even for a public fixture.
        private_hpo.chmod(0o600)
        case_inputs = tuple(CaseInput(_input_id(role, kind, item.sha256), "public-" + role,
            role, kind, item.sha256, item.size_bytes) for role, kind, item in (
                *(("vcf-" + str(index), CaseFileClass.VCF, row) for index, row in enumerate(vcfs, 1)),
                ("phenotype", CaseFileClass.HPO_LIST, hpo)))
        case = create_case_manifest(sample_id="public-validation-proband", inputs=case_inputs,
            configuration=ConfigurationBinding("genosapiens.public-policy-family", policies.policy_sha256),
            reference_bundle=ReferenceBundleBinding("public-annotation-benchmark", REFERENCE_MANIFEST_SHA256))
        case_content = case_manifest_as_bytes(case)
        publish_bytes(output, "case-manifest.json", case_content)
        case_path = output / "case-manifest.json"
        filter_policy = load_common_allele_filter_policy(_bound_bytes(
            resources.repository_root / "config/track1/common-allele-prefilter-v1.json",
            bindings.common_allele_policy_sha256), bindings.common_allele_policy_sha256)
        precision_content = read_bounded_metadata(resources.repository_root / "config/public-resources/gnomad-joint-faf95-precision-v1.json")
        precision_contract = load_gnomad_faf95_precision_contract(precision_content, _sha256(precision_content))
        precision = verify_gnomad_faf95_precision_qualification(_bound_bytes(
            resources.reference_root / GNOMAD_PRECISION_REPORT_PATH, bindings.gnomad_projection_precision_sha256),
            bindings.gnomad_projection_precision_sha256, precision_contract)
        from genosapiens.installed_public_index import load_installed_index_pin
        installed_pin = load_installed_index_pin(resources.repository_root)
        indexes = verify_common_allele_public_indexes(reference_root=resources.reference_root,
            manifest_path=installed_pin['manifest_path'] if installed_pin else COMMON_ALLELE_PUBLIC_INDEX_MANIFEST_PATH,
            expected_manifest_sha256=bindings.common_allele_public_index_manifest_sha256,
            verify_sources=False, reuse_qualified_reconciliation=True)
        limits = population.PopulationQueryLimits(config.public_queries.batch_rows, config.public_queries.batch_rows,
            65536, 1_000_000, config.public_queries.max_response_bytes, config.public_queries.timeout_seconds)
        renderer = load_submission_contig_renderer(identity_path=resources.population_reference.identity_path,
            reference_manifest_path=resources.population_reference.manifest_path,
            reference_root=resources.population_reference.reference_root,
            reference_binding=ReferenceRoleBinding("public_annotation_benchmark", reference.identity_id))
        def fetch(contig, position, length):
            return resources.reference.fetch(resources.source_to_reference_contigs[contig], position, length)
        from genosapiens.track1_early_scope import EarlyScopeSession, validate_scope_pin
        from genosapiens.bounded_json import load_bounded_json
        scope_pin, _ = load_bounded_json(resources.repository_root / 'track1-early-scope-pin.json',
            reference_root=resources.repository_root, maximum_bytes=65536)
        scope_pin = validate_scope_pin(scope_pin)
        sources = []
        stage = "call-preparation"
        for ordinal, item in enumerate(vcfs, 1):
            path = Path(item.path)
            concordance = evaluate_public_validation_vcf_concordance(path, resources.reference.fetch,
                {name: entry.length for name, entry in resources.reference.index.items()}, resources.source_to_reference_contigs,
                contract=resources.reference_contract, reference_manifest_sha256=REFERENCE_MANIFEST_SHA256)
            stream = iter_vcf_qc_records(path, concordance_report=concordance,
                max_record_bytes=config.input.max_record_bytes, max_header_bytes=config.input.max_header_bytes)
            try:
                sample = stream.header.sample_names[0]
            finally:
                stream.close()
            descriptor = formats[item.sha256]
            def native(parent):
                return NativeSourceContract(parent, ReferenceRole.PUBLIC_ANNOTATION_BENCHMARK,
                    reference.identity_id, REFERENCE_MANIFEST_SHA256, tuple(descriptor["mitochondrial_contigs"]),
                    descriptor["repeat_identity"])
            from genosapiens.wgs_source_routing import SourceRoutingContext
            def routing(context):
                return SourceRoutingContext(native(context.source_manifest.sha256), renderer)
            def decisions(context):
                route = routing(context)
                return PublicEarlyFilterDecisionProvider(filter_policy, precision, route,
                    PublicFilterWorkerConfig(resources.population_joint, resources.population_all_sites,
                        resources.population_reference, RUNTIME_IMAGE, limits, config.public_queries,
                        indexes, output / ("source-" + str(ordinal)) / "prefilter-staging"), scope_session)(context)
            worker_role = resources.reference_contract.reference_role("public_annotation_benchmark")
            worker_reference = source_worker_reference(
                resources.reference_root / resources.reference_contract.expected_output(
                    worker_role.fasta_output_role
                ).path,
                resources.reference_root / resources.reference_contract.expected_output(
                    worker_role.fai_output_role
                ).path,
                aliases=tuple(sorted(resources.source_to_reference_contigs.items())),
            )
            with EarlyScopeSession(resources.reference_root, scope_pin, renderer=renderer) as scope_session:
                result = prepare_bounded_case(output_root=output / ("source-" + str(ordinal)), original_input_path=path,
                    hpo_input_path=private_hpo, hpo_input_sha256=hpo.sha256, hpo_input_size_bytes=hpo.size_bytes,
                    hpo_ontology=resources.hpo_store.graph.ontology, hpo_reference_root=output_root,
                    concordance_report=concordance, fetch_reference=fetch,
                    source_worker_reference=worker_reference, selected_sample_name=sample,
                    selected_sample_sha256=_sha256(canonical({"schema": "genosapiens.public-source-sample",
                        "input_sha256": item.sha256, "sample_name": sample})), case_manifest_sha256=_sha256(case_content),
                    ontology_manifest_sha256=HPO_MANIFEST_SHA256, executable_stack_sha256=stack_sha,
                    runtime_image_digest=authority.runtime_image_digest, execution_config=config,
                    run_options=PreparationRunOptions(False),
                    policy_sha256=bindings.common_allele_policy_sha256, population_manifest_sha256=GNOMAD_JOINT_MANIFEST_SHA256,
                    clinical_manifest_sha256=CLINVAR_MANIFEST_SHA256, exception_manifest_sha256=bindings.common_allele_exception_manifest_sha256,
                    precision_manifest_sha256=bindings.gnomad_projection_precision_sha256,
                    scope_policy_sha256=scope_session.policy_sha256,
                    scope_manifest_sha256=scope_session.manifest.sha256,
                    decision_provider=decisions, source_routing_context_provider=routing,
                    base_release_bytes=None, expected_base_release_sha256=None, public_validation_authority_bytes=authority.as_bytes(),
                    expected_public_validation_authority_sha256=authority.sha256)
            if result.receipt.status != "passed":
                raise ValueError
            ledger = load_model(read_bounded_metadata(result.root / "source/manifest.json"), LedgerManifest)
            sources.append(load_prepared_intake_source(preparation_root=result.root, original_input_path=path,
                expected_prepared_sha256=result.prepared_case.manifest.sha256,
                expected_early_snapshot_sha256=result.prepared_case.early_snapshot.sha256,
                source_routing_context=SourceRoutingContext(native(ledger.sha256), renderer),
                execution_config=config, base_release_bytes=None,
                public_validation_authority_bytes=authority.as_bytes()))
        sources = tuple(sources)
        intake_root = output / "intake"
        intake = MergedReadIntake(merge_prepared_intake((), sources=sources, case_manifest_path=case_path,
            reference_root=resources.reference_root, output_directory=intake_root), intake_root)
        saved = publish_public_case(directory=output, authority=authority, sources=sources, intake=intake,
            case_manifest_path=case_path, execution_config=config, resources=resources)
        children = resource.getrusage(resource.RUSAGE_CHILDREN)
        from genosapiens.wgs_public_analysis import analysis_screening_arguments
        screening=analysis_screening_arguments(resources,bindings,indexes,
            common_policy=filter_policy,precision=precision,
            dependencies=(('merged_intake_manifest',str((intake_root/'merged-intake-manifest.json').absolute()),
                intake.manifest.sha256,len(intake.manifest.as_bytes())),))
        return PreparedPublicCase(case_id, view, materialized, saved, indexes,
            str(Fraction(time.monotonic_ns() - wall_start, 1_000_000_000)),
            str(Fraction(time.process_time_ns() - cpu_start, 1_000_000_000)
                + Fraction(str(children.ru_utime + children.ru_stime - children_start.ru_utime - children_start.ru_stime))),
            screening)
    except Exception:
        raise PublicTrack1ValidationError("public case preparation failed", preparation_stage=stage) from None


def execute_prepared_public_case(prepared, runtime, *, resources, policies, output_directory=None):
    """Run one bounded scientific invocation, reusing exact saved annotation."""
    from genosapiens.wgs_public_analysis import PublicBoundedAnalysisBuilder
    from genosapiens.wgs_public_execution import execute_public_case
    wall_start, cpu_start = time.monotonic_ns(), time.process_time_ns()
    children_start = resource.getrusage(resource.RUSAGE_CHILDREN)
    stage = "candidate-authentication"
    try:
        if (runtime.candidate.bindings not in prepared.public_case.authority.candidate_bindings
                or policies.policy_sha256 != prepared.public_case.authority.policy_family_sha256):
            raise ValueError
        stage = "bounded-scientific-execution"
        prepared.invocation_count += 1
        authority = prepared.public_case.authority
        builder = PublicBoundedAnalysisBuilder(resources, policies, runtime, prepared.public_case.sources[0].native_contract,
            authority.executable_stack_sha256, authority.runtime_image_digest, prepared.public_case.intake.manifest.sample_id,
            prepared.clinical_indexes, annotation_source=prepared.annotation_source,
            **prepared.screening_arguments)
        with _validation_session():
            result = execute_public_case(prepared.public_case, builder=builder,
                output_root=(Path(output_directory) if output_directory is not None else
                    prepared.public_case.case_manifest_path.parent / ("execution-" + str(prepared.invocation_count))))
        prepared.annotation_source = result.products.annotation_source
    except Exception:
        raise PublicTrack1ValidationError("public scientific execution failed", execution_stage=stage) from None
    children = resource.getrusage(resource.RUSAGE_CHILDREN)
    measurements = ResourceMeasurements(
        str(Fraction(prepared.preparation_wall_seconds) + Fraction(time.monotonic_ns() - wall_start, 1_000_000_000)),
        str(Fraction(prepared.preparation_cpu_seconds) + Fraction(time.process_time_ns() - cpu_start, 1_000_000_000)
            + Fraction(str(children.ru_utime + children.ru_stime - children_start.ru_utime - children_start.ru_stime))), None, None)
    return result, measurements
