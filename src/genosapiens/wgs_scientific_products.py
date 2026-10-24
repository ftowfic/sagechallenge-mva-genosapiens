"""Current normalized scientific product metadata, derived from one case owner.

Primary order and selection live in the complete primary artifact.
These small descriptors bind their roots and presentation metadata; consumers
validate them against one open ScientificArtifactReader before using them.
There are no alternate rich-row codecs or saved private-index authority here.
"""
from dataclasses import dataclass, fields
from contextlib import contextmanager
import hashlib
from pathlib import Path
from typing import ClassVar, Literal, Mapping

from genosapiens.evidence_fusion import PhenotypeTransform, ranking_ties
from genosapiens.evidence_model import _sha, _typed, _exact
from genosapiens.evidence_serialization import _canonical, _loads
from genosapiens.track1_sensitivity import DisclosedDiagnosis, PriorPolicy, RankTie
from genosapiens.track1_submission import SubmissionGroup, SubmissionReferenceBinding
from genosapiens.wgs_case_resources import bound_reference
from genosapiens.wgs_evidence_store import EvidenceReference, EvidenceStoreError
from genosapiens.wgs_phenotype_evidence import _hash

MAX_PRODUCT_BYTES = 8 * 1024**2


@dataclass(frozen=True,slots=True,repr=False)
class _Product:
    scientific_manifest_sha256: str
    root: EvidenceReference
    schema: ClassVar[str]
    root_kind: ClassVar[str]
    version: ClassVar[int] = 1

    def __post_init__(self):
        if type(self.root) is not EvidenceReference or self.root.kind != self.root_kind:
            raise EvidenceStoreError("foreign normalized scientific product root")
        for field in fields(self):
            value = getattr(self,field.name)
            if field.name.endswith("_sha256"):
                _sha(value)
            elif field.name.endswith("_count") and (type(value) is not int or not 0 <= value < 2**63):
                raise EvidenceStoreError("invalid normalized scientific product count")
        self.as_bytes()

    @property
    def evidence_manifest_sha256(self):
        return self.root.manifest_sha256

    def as_bytes(self):
        content = (_canonical({"schema":self.schema,"version":self.version,
            "classification":"controlled","llm_policy":"no-llm","product":self})+"\n").encode("ascii")
        if len(content) > MAX_PRODUCT_BYTES:
            raise EvidenceStoreError("normalized scientific product metadata exceeds capacity",capacity=True)
        return content

    @property
    def sha256(self):
        return hashlib.sha256(self.as_bytes()).hexdigest()


@dataclass(frozen=True,slots=True,repr=False)
class ScientificOrderProduct(_Product):
    policy_sha256: str
    tie_rules: tuple[RankTie,...]
    hypothesis_count: int
    schema: ClassVar[str] = "genosapiens.wgs-normalized-order-product"
    root_kind: ClassVar[str] = "complete-case-order-v1"


@dataclass(frozen=True,slots=True,repr=False)
class ScientificPrimaryProduct(_Product):
    order_manifest_sha256: str
    prefilter_manifest_sha256: str
    candidate_universe_sha256: str
    ranking_policy_id: str
    ranking_policy_sha256: str
    derivation_sha256: str
    scoring_view_sha256: str
    disease_profile_sha256: str
    disease_gene_evidence_sha256: str
    ic_sha256: str
    profile_evidence_id: str | None
    transform: PhenotypeTransform
    tie_rules: tuple[RankTie,...]
    prior_policy: PriorPolicy
    hypothesis_count: int
    schema: ClassVar[str] = "genosapiens.wgs-normalized-primary-product"
    root_kind: ClassVar[str] = "scientific-primary-evidence-v1"


@dataclass(frozen=True,slots=True,repr=False)
class ScientificSensitivityProduct(_Product):
    mode: Literal["dominant_term_sensitivity","diagnosis_informed_sensitivity"]
    status: Literal["applied","not_applicable","not_supplied"]
    reason: str | None
    primary_manifest_sha256: str
    prefilter_manifest_sha256: str
    candidate_universe_sha256: str
    ranking_policy_sha256: str
    hypothesis_count: int
    order_root: EvidenceReference
    source_root: EvidenceReference
    dominant_term: Mapping[str,object] | None
    removed_observation_ids: tuple[str,...]
    diagnosis: DisclosedDiagnosis | None
    prior_policy: PriorPolicy
    schema: ClassVar[str] = "genosapiens.wgs-normalized-sensitivity-product"
    root_kind: ClassVar[str] = "scientific-case-views-v1"

    def __post_init__(self):
        _Product.__post_init__(self)
        kind = "case-ablation-evidence-v1" if self.mode == "dominant_term_sensitivity" else "case-diagnosis-evidence-v1"
        if (self.mode not in ("dominant_term_sensitivity","diagnosis_informed_sensitivity")
                or type(self.order_root) is not EvidenceReference or self.order_root.kind != "complete-case-order-v1"
                or type(self.source_root) is not EvidenceReference or self.source_root.kind != kind
                or any(r.manifest_sha256 != self.evidence_manifest_sha256 for r in (self.order_root,self.source_root))
                or self.status not in ("applied","not_applicable","not_supplied")
                or (self.status == "applied") != (self.reason is None)
                or self.mode == "dominant_term_sensitivity" and self.diagnosis is not None
                or self.mode == "diagnosis_informed_sensitivity" and (self.dominant_term is not None or self.removed_observation_ids)
                or self.status == "not_supplied" and (self.mode != "diagnosis_informed_sensitivity"
                    or self.diagnosis is not None or self.reason != "diagnosis_not_supplied")):
            raise EvidenceStoreError("invalid normalized sensitivity product")


@dataclass(frozen=True,slots=True,repr=False)
class ScientificSelectionProduct(_Product):
    order_manifest_sha256: str
    ranking_policy_sha256: str
    reference: SubmissionReferenceBinding
    selected: tuple[SubmissionGroup,...]
    selected_root: EvidenceReference
    ranked_hypothesis_count: int
    submitted_group_count: int
    eligible_tail_group_count: int
    excluded_hypothesis_count: int
    alternate_hypothesis_count: int
    schema: ClassVar[str] = "genosapiens.wgs-normalized-selection-product"
    root_kind: ClassVar[str] = "complete-submission-selection-v1"
    maximum_groups: ClassVar[int] = 10
    epcr_assigned: ClassVar[bool] = False
    csv_emitted: ClassVar[bool] = False
    controlled_use_ready: ClassVar[bool] = False

    def __post_init__(self):
        _Product.__post_init__(self)
        if (len(self.selected) != self.submitted_group_count or len(self.selected) > self.maximum_groups
                or any(g.alternative_hypothesis_ids for g in self.selected)
                or tuple(g.source_rank for g in self.selected) != tuple(sorted({g.source_rank for g in self.selected}))
                or self.ranked_hypothesis_count != sum((self.submitted_group_count,self.eligible_tail_group_count,
                    self.excluded_hypothesis_count,self.alternate_hypothesis_count))
                or type(self.selected_root) is not EvidenceReference
                or self.selected_root.kind != "scientific-selected-groups-v1"
                or self.selected_root.manifest_sha256 != self.evidence_manifest_sha256):
            raise EvidenceStoreError("invalid normalized selection product")


PRODUCT_TYPES = (ScientificOrderProduct,ScientificPrimaryProduct,ScientificSensitivityProduct,ScientificSelectionProduct)


def scientific_product_from_bytes(content, *, expected_sha256, product_type):
    try:
        _sha(expected_sha256)
        if (product_type not in PRODUCT_TYPES or type(content) is not bytes or len(content) > MAX_PRODUCT_BYTES
                or hashlib.sha256(content).hexdigest() != expected_sha256):
            raise EvidenceStoreError("normalized scientific product identity changed")
        data = _loads(content)
        if (type(data) is not dict or set(data) != {"schema","version","classification","llm_policy","product"}
                or data["schema"] != product_type.schema or type(data["version"]) is not int or data["version"] != 1
                or data["classification"] != "controlled" or data["llm_policy"] != "no-llm"):
            raise EvidenceStoreError("obsolete or malformed scientific product")
        value = _typed(data["product"],product_type,decode=True)
        if value.as_bytes() != content:
            raise EvidenceStoreError("noncanonical normalized scientific product")
        return value
    except (ValueError,TypeError,KeyError,AttributeError):
        raise EvidenceStoreError("invalid normalized scientific product") from None


@dataclass(frozen=True,slots=True,repr=False)
class ScientificProducts:
    order: ScientificOrderProduct
    primary: ScientificPrimaryProduct
    selection: ScientificSelectionProduct


@dataclass(frozen=True,slots=True,repr=False)
class ScientificSensitivityProducts:
    dominant: ScientificSensitivityProduct
    diagnosis: ScientificSensitivityProduct


def scientific_products(reader, *, expected_prefilter_manifest_sha256):
    """Derive every small descriptor from the operation's verified complete case."""
    from genosapiens.wgs_scientific_artifact import ScientificArtifactReader, ScopedScientificArtifact
    if type(reader) is not ScientificArtifactReader:
        raise EvidenceStoreError("verified complete scientific artifact owner required")
    reader._require_open()
    if type(reader.artifact) is ScopedScientificArtifact and reader.artifact.completion_status == 'incomplete':
        raise EvidenceStoreError('incomplete science cannot supply complete submission products')
    _sha(expected_prefilter_manifest_sha256)
    if reader.parents.get("prefilter_manifest") != expected_prefilter_manifest_sha256:
        raise EvidenceStoreError("normalized scientific prefilter parent changed")
    with reader.reader.read_scope():
        case,resources,artifact = reader.case,reader.resources,reader.artifact
        def root(value):
            return bound_reference(reader.reader,value,value["kind"])
        order = ScientificOrderProduct(artifact.sha256,root(case.orders["primary"].reference),
            resources.fusion_policy_sha256,ranking_ties(resources.fusion_policy),case.hypothesis_count)
        primary = ScientificPrimaryProduct(artifact.sha256,root(case.primary_reference),order.sha256,
            expected_prefilter_manifest_sha256,case.identities.sha256,resources.fusion_policy.policy_id,
            resources.fusion_policy_sha256,resources.fusion_policy.derivation_sha256,
            resources.header["view_sha256"],resources.header["profile_sha256"],
            resources.gene_logical_sha256,resources.header['ic_sha256'],resources.profile_evidence_id,
            resources.fusion_policy.phenotype,order.tie_rules,resources.prior_policy,case.hypothesis_count)
        selection = ScientificSelectionProduct(artifact.sha256,root(case.primary_header["selection"]),order.sha256,
            resources.fusion_policy_sha256,resources.reference_binding,case.selected_groups,root(case.header["selected"]),
            case.hypothesis_count,len(case.selection.selected),case.selection.tail_count,
            case.selection.excluded_count,case.selection.alternative_count)
        return ScientificProducts(order,primary,selection)


def validate_scientific_products(reader, products, *, expected_prefilter_manifest_sha256):
    if type(products) is not ScientificProducts or products != scientific_products(
            reader,expected_prefilter_manifest_sha256=expected_prefilter_manifest_sha256):
        raise EvidenceStoreError("normalized scientific products do not replay")


@contextmanager
def open_scientific_products(artifact, products, *, directory, expected_parent_bindings, scratch_directory,
                             eligibility_root=None, profile_root=None, _reader=None):
    """One replay owner shared by every consumer in a public operation."""
    from genosapiens.wgs_scientific_artifact import ScientificArtifactReader
    if type(products) is not ScientificProducts:
        raise EvidenceStoreError("normalized scientific products required")
    if _reader is not None:
        directory = Path(directory)
        if (type(_reader) is not ScientificArtifactReader
                or _reader.artifact != artifact or _reader.directory != directory
                or _reader.parents != expected_parent_bindings
                or _reader.expected_reference_binding != products.selection.reference
                or _reader.eligibility_root != (directory.parent / 'eligibility' if eligibility_root is None else Path(eligibility_root))
                or _reader.profile_root != (directory.parent / 'analysis/profile' if profile_root is None else Path(profile_root))):
            raise EvidenceStoreError('foreign live scientific publication reader')
        _reader._require_open()
        if (_reader.reader._read_scope_depth or _reader.reader._active_derivative is not None
                or _reader.reader._source_writer is not None):
            raise EvidenceStoreError('scientific publication requires a closed input scope')
        with _reader.reader.read_scope():
            validate_scientific_products(_reader, products,
                expected_prefilter_manifest_sha256=products.primary.prefilter_manifest_sha256)
            yield _reader
        return
    with ScientificArtifactReader(artifact, directory=directory,
            expected_parent_bindings=expected_parent_bindings,
            expected_reference_binding=products.selection.reference,
            scratch_directory=scratch_directory,eligibility_root=eligibility_root,
            profile_root=profile_root) as reader, reader.reader.read_scope():
        validate_scientific_products(reader, products,
            expected_prefilter_manifest_sha256=products.primary.prefilter_manifest_sha256)
        yield reader


def product_descriptors(value):
    """Extract the three current primary metadata products."""
    return ScientificProducts(value.order, value.primary, value.selection)


def scientific_report_roots(reader, products):
    """Root closure retained by the current audit presentation."""
    reader._require_open()
    artifact, evidence = reader.artifact, reader.reader
    return tuple(sorted({
        "case": evidence.reference(artifact.case_name, kind="complete-scientific-case-v1"),
        "resources": evidence.reference(artifact.resources_name, kind="scientific-case-resources-v1"),
        "reference": evidence.reference(artifact.reference_name, kind="submission-reference-metadata-v1"),
        "primary": products.primary.root, "primary_order": products.order.root,
        "selection": products.selection.root, "selected": products.selection.selected_root,
    }.items()))
