"""Versioned conservative common-allele prefilter contracts.

The pure decision engine is deliberately separate from resource I/O.  Loaders
authenticate canonical policy bytes; adapters must provide typed, authenticated
facts and may never pass phenotype or outcome-dependent values into this API.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass
from fractions import Fraction
import hashlib
import json
import os
from pathlib import Path
import re
import sqlite3
from typing import Any, Literal

from genosapiens.governance import ArtifactClassification
from genosapiens.llm_guard import LLMPolicy


POLICY_SCHEMA = "genosapiens.common-allele-prefilter-policy"
POLICY_VERSION = 1
POLICY_ID = "track1-common-allele-prefilter-v1"
_SHA256 = re.compile(r"^[0-9a-f]{64}$")
_ID = re.compile(r"^[A-Za-z0-9][A-Za-z0-9._/-]{0,127}$")
_ERRORS = (ValueError, TypeError, KeyError, UnicodeError, RecursionError, OverflowError)


class CommonAlleleFilterError(ValueError):
    """A static policy/decision error that never includes controlled values."""


def _fail(message: str = "invalid common-allele filter policy") -> None:
    raise CommonAlleleFilterError(message)


def _pairs(items: list[tuple[str, Any]]) -> dict[str, Any]:
    result: dict[str, Any] = {}
    for key, value in items:
        if key in result:
            _fail()
        result[key] = value
    return result


def _constant(_value: str) -> None:
    _fail()


def _canonical(value: object) -> bytes:
    try:
        return (json.dumps(value, sort_keys=True, separators=(",", ":"),
                           ensure_ascii=True, allow_nan=False) + "\n").encode("ascii")
    except _ERRORS:
        _fail("cannot encode common-allele filter policy")


def _object(value: object, fields: set[str]) -> dict[str, Any]:
    if type(value) is not dict or set(value) != fields:
        _fail()
    return value


def _text(value: object) -> str:
    if type(value) is not str or _ID.fullmatch(value) is None:
        _fail()
    return value


def _texts(value: object, *, nonempty: bool = True) -> tuple[str, ...]:
    if type(value) is not list or (nonempty and not value):
        _fail()
    result = tuple(_text(item) for item in value)
    if len(result) != len(set(result)):
        _fail()
    return result


def _literal(value: object, expected: object) -> None:
    if type(value) is not type(expected) or value != expected:
        _fail()


@dataclass(frozen=True, slots=True)
class RationalPolicyValue:
    numerator: int
    denominator: int

    def __post_init__(self) -> None:
        if (type(self.numerator) is not int or type(self.denominator) is not int
                or self.numerator <= 0 or self.denominator <= 0
                or Fraction(self.numerator, self.denominator).numerator != self.numerator
                or Fraction(self.numerator, self.denominator).denominator != self.denominator):
            _fail()

    @property
    def fraction(self) -> Fraction:
        return Fraction(self.numerator, self.denominator)


@dataclass(frozen=True, slots=True)
class FrequencyPolicy:
    dataset: Literal["gnomad-joint-v4.1"]
    subset: Literal["joint"]
    metric: Literal["FAF95"]
    maximum_frequency_field: Literal["JT_FAFMAX95"]
    maximum_ancestry_field: Literal["JT_FAFMAX95_ANC"]
    ancestry_an_field: Literal["JT_AN_ANC"]
    ancestry_ac_field: Literal["JT_AC_ANC"]
    ancestry_order_binding: Literal["projection.major_genetic_ancestries"]
    threshold: RationalPolicyValue
    comparison: Literal["strictly_greater"]
    minimum_matching_ancestry_an: int
    allowed_maximum_ancestries: tuple[str, ...]
    prohibited_substitutes: tuple[str, ...]

    def __post_init__(self) -> None:
        if (self.dataset != "gnomad-joint-v4.1" or self.subset != "joint"
                or self.metric != "FAF95" or self.maximum_frequency_field != "JT_FAFMAX95"
                or self.maximum_ancestry_field != "JT_FAFMAX95_ANC"
                or self.ancestry_an_field != "JT_AN_ANC" or self.ancestry_ac_field != "JT_AC_ANC"
                or self.ancestry_order_binding != "projection.major_genetic_ancestries"
                or self.threshold.fraction != Fraction(1, 20)
                or self.comparison != "strictly_greater"
                or type(self.minimum_matching_ancestry_an) is not int
                or self.minimum_matching_ancestry_an != 2000
                or self.allowed_maximum_ancestries != ("afr", "amr", "eas", "nfe", "sas")
                or not {"JT_GRPMAX_AF", "JT_GRPMAX_AN", "JT_AN", "genome_only_af",
                        "callset_info_af", "patient_allele_balance"}.issubset(self.prohibited_substitutes)):
            _fail()


@dataclass(frozen=True, slots=True)
class PrecisionPolicy:
    source_storage: Literal["float32_bcf_with_text_export"]
    comparison_interval: Literal["lower_bound_must_exceed_threshold"]
    require_evidence: Literal[True]

    def __post_init__(self) -> None:
        if (self.source_storage != "float32_bcf_with_text_export"
                or self.comparison_interval != "lower_bound_must_exceed_threshold"
                or self.require_evidence is not True):
            _fail()


@dataclass(frozen=True, slots=True)
class PopulationQualityPolicy:
    projection_filter_is_source_pass: Literal[False]
    exome_filters_defined_field: str
    exome_filters_field: str
    genome_filters_defined_field: str
    genome_filters_field: str
    exome_interval_qc_field: str
    exome_not_called_field: str
    genome_not_called_field: str
    joint_discordance_field: str
    joint_discordance_p_field: str
    discordance_minimum_p: RationalPolicyValue
    unknown_warning_behavior: Literal["retain_unknown"]

    def __post_init__(self) -> None:
        expected = (
            "EX_FILTERS_DEFINED", "EX_FILTERS", "GN_FILTERS_DEFINED", "GN_FILTERS",
            "REGION_FAIL_INTERVAL_QC", "REGION_NOT_CALLED_IN_EXOMES",
            "REGION_NOT_CALLED_IN_GENOMES", "DISC", "DISC_UNION_P",
        )
        actual = (
            self.exome_filters_defined_field, self.exome_filters_field,
            self.genome_filters_defined_field, self.genome_filters_field,
            self.exome_interval_qc_field, self.exome_not_called_field,
            self.genome_not_called_field, self.joint_discordance_field,
            self.joint_discordance_p_field,
        )
        if (self.projection_filter_is_source_pass is not False or actual != expected
                or self.discordance_minimum_p.fraction != Fraction(1, 10000)
                or self.unknown_warning_behavior != "retain_unknown"):
            _fail()


@dataclass(frozen=True, slots=True)
class ClinicalProtectionPolicy:
    protected_significance: tuple[str, ...]
    unprotected_significance: tuple[str, ...]
    conflict_markers: tuple[str, ...]
    vus_alone_protects: Literal[False]
    unrecognized_behavior: Literal["retain_unknown"]

    def __post_init__(self) -> None:
        if (self.vus_alone_protects is not False or self.unrecognized_behavior != "retain_unknown"
                or "pathogenic" not in self.protected_significance
                or "likely_pathogenic" not in self.protected_significance
                or "uncertain_significance" not in self.unprotected_significance
                or not self.conflict_markers
                or set(self.protected_significance) & set(self.unprotected_significance)):
            _fail()


@dataclass(frozen=True, slots=True)
class ProtectionPolicy:
    curated_exact_allele: Literal[True]
    curated_gene_interval: Literal[True]
    gene_interval_includes_introns: Literal[True]
    explicit_phase: Literal[True]
    unresolved_technical_obligation: Literal[True]
    representation_classes: tuple[str, ...]
    native_classes: tuple[str, ...]

    def __post_init__(self) -> None:
        required_representations = {"mnv", "complex_delins", "ambiguous_composite",
                                    "overlapping_composite"}
        required_native = {"x_linked", "y_linked", "mitochondrial", "structural",
                           "repeat", "mosaic"}
        if (self.curated_exact_allele is not True or self.curated_gene_interval is not True
                or self.gene_interval_includes_introns is not True or self.explicit_phase is not True
                or self.unresolved_technical_obligation is not True
                or not required_representations.issubset(self.representation_classes)
                or not required_native.issubset(self.native_classes)):
            _fail()


@dataclass(frozen=True, slots=True)
class CommonAlleleFilterPolicy:
    frequency: FrequencyPolicy
    precision: PrecisionPolicy
    population_quality: PopulationQualityPolicy
    clinical_protection: ClinicalProtectionPolicy
    protections: ProtectionPolicy
    required_resource_roles: tuple[str, ...]
    prohibited_inputs: tuple[str, ...]
    unknown_behavior: Literal["retain_unknown"] = "retain_unknown"
    policy_id: Literal["track1-common-allele-prefilter-v1"] = POLICY_ID
    schema: Literal["genosapiens.common-allele-prefilter-policy"] = POLICY_SCHEMA
    version: Literal[1] = POLICY_VERSION

    def __post_init__(self) -> None:
        required_inputs = {"hpo", "phenotype", "diagnosis", "score", "rank", "selected_rows",
                           "validation_truth", "previous_campaign_outcome"}
        required_resources = {"gnomad_joint_projection", "gnomad_all_sites", "clinvar",
                              "open_targets_variant", "open_targets_eva",
                              "common_allele_exceptions", "source_reference",
                              "precision_qualification"}
        if (self.schema != POLICY_SCHEMA or self.version != POLICY_VERSION
                or self.policy_id != POLICY_ID or self.unknown_behavior != "retain_unknown"
                or not required_inputs.issubset(self.prohibited_inputs)
                or not required_resources.issubset(self.required_resource_roles)):
            _fail()

    @property
    def threshold(self) -> Fraction:
        return self.frequency.threshold.fraction


def _rational(value: object) -> RationalPolicyValue:
    data = _object(value, {"numerator", "denominator"})
    return RationalPolicyValue(data["numerator"], data["denominator"])


def load_common_allele_filter_policy(
    content: bytes, expected_sha256: str
) -> CommonAlleleFilterPolicy:
    """Load exact canonical bytes bound to an independently supplied digest."""

    try:
        if (type(content) is not bytes or type(expected_sha256) is not str
                or _SHA256.fullmatch(expected_sha256) is None
                or hashlib.sha256(content).hexdigest() != expected_sha256):
            _fail()
        data = json.loads(content.decode("ascii"), object_pairs_hook=_pairs,
                          parse_constant=_constant)
        if _canonical(data) != content:
            _fail()
        top = _object(data, {
            "schema", "version", "policy_id", "frequency", "precision",
            "population_quality", "clinical_protection", "protections",
            "required_resource_roles", "outcome_independent_inputs", "unknown_behavior",
        })
        frequency = _object(top["frequency"], {
            "dataset", "subset", "metric", "maximum_frequency_field",
            "maximum_ancestry_field", "ancestry_an_field", "ancestry_ac_field",
            "ancestry_order_binding", "threshold", "comparison",
            "minimum_matching_ancestry_an", "allowed_maximum_ancestries",
            "prohibited_substitutes",
        })
        precision = _object(top["precision"], {
            "source_storage", "comparison_interval", "require_evidence",
        })
        quality = _object(top["population_quality"], {
            "projection_filter_is_source_pass", "exome_filters_defined_field",
            "exome_filters_field", "genome_filters_defined_field", "genome_filters_field",
            "exome_interval_qc_field", "exome_not_called_field", "genome_not_called_field",
            "joint_discordance_field", "joint_discordance_p_field",
            "discordance_minimum_p", "unknown_warning_behavior",
        })
        clinical = _object(top["clinical_protection"], {
            "protected_significance", "unprotected_significance", "conflict_markers",
            "vus_alone_protects", "unrecognized_behavior",
        })
        protections = _object(top["protections"], {
            "curated_exact_allele", "curated_gene_interval", "gene_interval_includes_introns",
            "explicit_phase", "unresolved_technical_obligation", "representation_classes",
            "native_classes",
        })
        inputs = _object(top["outcome_independent_inputs"], {"prohibited"})
        return CommonAlleleFilterPolicy(
            frequency=FrequencyPolicy(
                **{key: frequency[key] for key in frequency if key not in {
                    "threshold", "allowed_maximum_ancestries", "prohibited_substitutes"}},
                threshold=_rational(frequency["threshold"]),
                allowed_maximum_ancestries=_texts(frequency["allowed_maximum_ancestries"]),
                prohibited_substitutes=_texts(frequency["prohibited_substitutes"]),
            ),
            precision=PrecisionPolicy(**precision),
            population_quality=PopulationQualityPolicy(
                **{key: quality[key] for key in quality if key != "discordance_minimum_p"},
                discordance_minimum_p=_rational(quality["discordance_minimum_p"]),
            ),
            clinical_protection=ClinicalProtectionPolicy(
                protected_significance=_texts(clinical["protected_significance"]),
                unprotected_significance=_texts(clinical["unprotected_significance"]),
                conflict_markers=_texts(clinical["conflict_markers"]),
                vus_alone_protects=clinical["vus_alone_protects"],
                unrecognized_behavior=clinical["unrecognized_behavior"],
            ),
            protections=ProtectionPolicy(
                **{key: protections[key] for key in protections if key not in {
                    "representation_classes", "native_classes"}},
                representation_classes=_texts(protections["representation_classes"]),
                native_classes=_texts(protections["native_classes"]),
            ),
            required_resource_roles=_texts(top["required_resource_roles"]),
            prohibited_inputs=_texts(inputs["prohibited"]),
            unknown_behavior=top["unknown_behavior"],
            policy_id=top["policy_id"],
            schema=top["schema"],
            version=top["version"],
        )
    except CommonAlleleFilterError:
        raise
    except _ERRORS:
        _fail()


def common_allele_filter_policy_as_json(policy: CommonAlleleFilterPolicy) -> bytes:
    """Return the exact v1 canonical policy representation."""

    if type(policy) is not CommonAlleleFilterPolicy:
        _fail()
    value = asdict(policy)
    value["frequency"]["allowed_maximum_ancestries"] = list(
        value["frequency"]["allowed_maximum_ancestries"]
    )
    value["frequency"]["prohibited_substitutes"] = list(
        value["frequency"]["prohibited_substitutes"]
    )
    value["clinical_protection"]["protected_significance"] = list(
        value["clinical_protection"]["protected_significance"]
    )
    value["clinical_protection"]["unprotected_significance"] = list(
        value["clinical_protection"]["unprotected_significance"]
    )
    value["clinical_protection"]["conflict_markers"] = list(
        value["clinical_protection"]["conflict_markers"]
    )
    value["protections"]["representation_classes"] = list(
        value["protections"]["representation_classes"]
    )
    value["protections"]["native_classes"] = list(value["protections"]["native_classes"])
    value["required_resource_roles"] = list(value.pop("required_resource_roles"))
    prohibited = list(value.pop("prohibited_inputs"))
    value["outcome_independent_inputs"] = {"prohibited": prohibited}
    return _canonical(value)


SOURCE_RECORD_ID = re.compile(r"^source-record-[0-9a-f]{24}$")
SOURCE_ALLELE_ID = re.compile(r"^source-allele-[0-9a-f]{24}$")
NORMALIZED_ALLELE_ID = re.compile(r"^allele-[0-9a-f]{24}$")
DECISION_SCHEMA = "genosapiens.common-allele-prefilter-decision"
CLINICAL_ROLES = (
    "clinvar", "open_targets_variant", "open_targets_eva",
    "common_allele_exceptions",
)
PREPARATION_DISPOSITIONS = ("normalized", "native_class", "quarantined", "non_variant")
SOURCE_CLASSES = (
    "snv", "mnv", "insertion", "deletion", "delins", "symbolic",
    "breakend", "spanning-deletion", "non-variant", "other",
)
NORMALIZED_KINDS = ("snv", "mnv", "insertion", "deletion", "delins")
CONTIG_CLASSES = ("autosomal", "x_linked", "y_linked", "mitochondrial", "other")
PHASE_EVIDENCE = ("phased_gt", "phase_set", "pgt_pid", "phase_conflicting",
                  "phase_incomplete", "phase_unresolved")
REPRESENTATION_FLAGS = (
    "mnv", "complex_delins", "ambiguous_composite", "overlapping_composite",
)
NATIVE_CLASSES = ("x_linked", "y_linked", "mitochondrial", "structural", "repeat", "mosaic")
DECISION_DISPOSITIONS = ("retain", "retain_unknown", "exclude_common")


def _decision_fail() -> None:
    raise CommonAlleleFilterError("invalid common-allele filter decision input")


def _decision_text(value: object) -> bool:
    return (
        isinstance(value, str)
        and bool(value)
        and value.strip() == value
        and not any(character in value for character in "\t\r\n")
    )


def _decision_digest(value: object) -> bool:
    return isinstance(value, str) and _SHA256.fullmatch(value) is not None


def _decision_count(value: object) -> bool:
    return isinstance(value, int) and not isinstance(value, bool) and value >= 0


def _ordered_values(values: object, allowed: set[str] | None = None) -> bool:
    return (
        isinstance(values, tuple)
        and values == tuple(sorted(set(values)))
        and all(_decision_text(value) for value in values)
        and (allowed is None or set(values) <= allowed)
    )


@dataclass(frozen=True, slots=True)
class RationalEvidence:
    numerator: int
    denominator: int

    def __post_init__(self) -> None:
        if (
            not _decision_count(self.numerator)
            or not isinstance(self.denominator, int)
            or isinstance(self.denominator, bool)
            or self.denominator <= 0
            or Fraction(self.numerator, self.denominator).numerator != self.numerator
            or Fraction(self.numerator, self.denominator).denominator != self.denominator
        ):
            _decision_fail()

    @property
    def fraction(self) -> Fraction:
        return Fraction(self.numerator, self.denominator)


@dataclass(frozen=True, slots=True)
class PrefilterOccurrence:
    source_record_id: str
    source_allele_id: str
    record_number: int
    allele_number: int
    normalized_allele_id: str | None
    preparation_disposition: Literal[
        "normalized", "native_class", "quarantined", "non_variant"
    ]
    contig_class: Literal[
        "autosomal", "x_linked", "y_linked", "mitochondrial", "other"
    ]
    source_variant_class: str
    normalized_kind: str | None
    reference_equivalence: Literal["exact", "uncertain", "unavailable"]
    phase_evidence: tuple[str, ...] = ()
    representation_flags: tuple[str, ...] = ()
    native_classes: tuple[str, ...] = ()
    technical_obligations: tuple[str, ...] = ()

    def __post_init__(self) -> None:
        normalized = self.preparation_disposition == "normalized"
        if (
            SOURCE_RECORD_ID.fullmatch(self.source_record_id) is None
            or SOURCE_ALLELE_ID.fullmatch(self.source_allele_id) is None
            or not _decision_count(self.record_number)
            or self.record_number <= 0
            or not _decision_count(self.allele_number)
            or self.allele_number <= 0
            or self.preparation_disposition not in PREPARATION_DISPOSITIONS
            or self.contig_class not in CONTIG_CLASSES
            or self.source_variant_class not in SOURCE_CLASSES
            or self.normalized_kind not in (*NORMALIZED_KINDS, None)
            or normalized != (
                isinstance(self.normalized_allele_id, str)
                and NORMALIZED_ALLELE_ID.fullmatch(self.normalized_allele_id) is not None
            )
            or normalized != (self.normalized_kind is not None)
            or self.reference_equivalence not in {"exact", "uncertain", "unavailable"}
            or not _ordered_values(self.phase_evidence, set(PHASE_EVIDENCE))
            or not _ordered_values(
                self.representation_flags, set(REPRESENTATION_FLAGS)
            )
            or not _ordered_values(self.native_classes, set(NATIVE_CLASSES))
            or not _ordered_values(self.technical_obligations)
        ):
            _decision_fail()


@dataclass(frozen=True, slots=True)
class PopulationNonContributionProof:
    """Supplemental authenticated all-sites witness; never variant count data.

    Producers authenticate the complete resource/query and reference geometry.
    This pure value preserves those bindings for strict eligibility replay.
    """
    subset: Literal['exomes']
    assembly: Literal['GRCh38']
    reference_identity_id: str
    reference_manifest_sha256: str
    normalized_allele_id: str
    contig: str
    position: int
    reference: str
    alternate: str
    checked_start: int
    checked_end: int
    coordinate_convention: Literal['one_based_inclusive']
    all_sites_resource_sha256: str
    all_sites_contract_sha256: str
    all_sites_index_sha256: str
    query_identity_sha256: str
    query_status: Literal['matched', 'no_match']
    query_output_sha256: str
    proof_basis: Literal['explicit_zero_callable', 'absent_from_complete_all_sites']
    raw_variant_ac: int | None
    raw_variant_an: int | None
    source_fields_sha256: str
    schema: Literal['genosapiens.population-noncontribution-proof'] = 'genosapiens.population-noncontribution-proof'
    version: Literal[1] = 1

    def __post_init__(self):
        digests = (self.reference_manifest_sha256, self.all_sites_resource_sha256,
                   self.all_sites_contract_sha256, self.all_sites_index_sha256,
                   self.query_identity_sha256, self.query_output_sha256, self.source_fields_sha256)
        if (self.schema != 'genosapiens.population-noncontribution-proof'
                or type(self.version) is not int or self.version != 1
                or self.subset != 'exomes' or self.assembly != 'GRCh38'
                or not _decision_text(self.reference_identity_id)
                or not _decision_text(self.normalized_allele_id)
                or NORMALIZED_ALLELE_ID.fullmatch(self.normalized_allele_id) is None
                or type(self.contig) is not str or re.fullmatch(r'(?:chr)?(?:[1-9]|1[0-9]|2[0-2])', self.contig) is None
                or any(not _decision_count(v) or v < 1 for v in
                       (self.position, self.checked_start, self.checked_end))
                or self.coordinate_convention != 'one_based_inclusive'
                or type(self.reference) is not str or re.fullmatch('[ACGT]+', self.reference) is None
                or type(self.alternate) is not str or re.fullmatch('[ACGT]+', self.alternate) is None
                or self.reference == self.alternate
                or not (len(self.reference) == len(self.alternate) == 1
                        or len(self.reference) == 1 < len(self.alternate) and self.alternate.startswith(self.reference)
                        or len(self.alternate) == 1 < len(self.reference) and self.reference.startswith(self.alternate))
                or self.checked_start != self.position
                or self.checked_end != self.position+len(self.reference)-1
                or any(not _decision_digest(v) for v in digests)
                or self.proof_basis not in ('explicit_zero_callable', 'absent_from_complete_all_sites')
                or self.query_status != ('matched' if self.proof_basis == 'explicit_zero_callable' else 'no_match')
                or any(v is not None and not _decision_count(v) for v in (self.raw_variant_ac, self.raw_variant_an))):
            _decision_fail()

    def as_bytes(self):
        return _canonical(asdict(self))

    @property
    def sha256(self):
        return hashlib.sha256(self.as_bytes()).hexdigest()


def load_population_noncontribution_proof(content: bytes, expected_sha256: str):
    try:
        if type(content) is not bytes or len(content) > 8192 or not _decision_digest(expected_sha256):
            _decision_fail()
        if hashlib.sha256(content).hexdigest() != expected_sha256:
            _decision_fail()
        data = json.loads(content, object_pairs_hook=_pairs, parse_constant=_constant)
        _object(data, set(PopulationNonContributionProof.__dataclass_fields__))
        result = PopulationNonContributionProof(**data)
        if result.as_bytes() != content:
            _decision_fail()
        return result
    except _ERRORS:
        _decision_fail()


@dataclass(frozen=True, slots=True)
class PopulationSubsetFacts:
    subset: Literal["exomes", "genomes"]
    allele_count: int | None
    allele_number: int | None
    count_status: Literal["valid", "missing", "inconsistent"]
    contribution: Literal["observed", "proven_no_callable", "unknown"]
    source_filters_defined: bool | None
    source_filters: tuple[str | None, ...] | None
    variant_not_called: bool | None
    interval_qc_pass: bool | None
    noncontribution_proof: PopulationNonContributionProof | None = None

    def __post_init__(self) -> None:
        if (
            self.subset not in {"exomes", "genomes"}
            or self.allele_count is not None and not _decision_count(self.allele_count)
            or self.allele_number is not None and not _decision_count(self.allele_number)
            or self.count_status not in {"valid", "missing", "inconsistent"}
            or self.contribution not in {"observed", "proven_no_callable", "unknown"}
            or self.source_filters_defined not in {True, False, None}
            or self.source_filters is not None and (
                not isinstance(self.source_filters, tuple)
                or any(value is not None and not _decision_text(value)
                       for value in self.source_filters)
            )
            or self.source_filters_defined is True and self.source_filters is None
            or self.source_filters_defined is not True and self.source_filters is not None
            or self.variant_not_called not in {True, False, None}
            or self.interval_qc_pass not in {True, False, None}
            or self.subset == "genomes" and self.interval_qc_pass is not None
            or self.count_status == "valid" and (
                self.allele_count is None or self.allele_number is None
            )
            or self.count_status == "missing" and (
                self.allele_count is not None and self.allele_number is not None
            )
            or self.contribution == "observed" and (
                self.allele_number is None or self.allele_number <= 0
            )
            or self.noncontribution_proof is not None and (
                type(self.noncontribution_proof) is not PopulationNonContributionProof
                or self.subset != 'exomes'
                or self.noncontribution_proof.raw_variant_ac != self.allele_count
                or self.noncontribution_proof.raw_variant_an != self.allele_number
            )
            or self.contribution == "proven_no_callable" and not (
                self.noncontribution_proof is not None
                and self.count_status == 'missing'
                and self.allele_count in (None, 0) and self.allele_number in (None, 0)
            ) and (
                self.count_status != "valid" or self.allele_number != 0
                or self.allele_count != 0
            )
        ):
            _decision_fail()


@dataclass(frozen=True, slots=True)
class PopulationFilterFacts:
    query_status: Literal[
        "exact_match", "no_match", "unassessable", "unavailable", "conflicting"
    ]
    maximum_raw_token: str | None
    maximum_ancestry: str | None
    precision_lower_bound: RationalEvidence | None
    precision_upper_bound: RationalEvidence | None
    precision_status: Literal["qualified", "missing", "unassessable"]
    matching_ancestry_ac: int | None
    matching_ancestry_an: int | None
    joint_ac: int | None
    joint_an: int | None
    exomes: PopulationSubsetFacts | None
    genomes: PopulationSubsetFacts | None
    discordant: bool | None
    discordance_p: RationalEvidence | None
    comparison_status: Literal["applicable", "proved_not_applicable", "unknown"]
    technical_warnings: tuple[str, ...]
    population_manifest_sha256: str
    precision_manifest_sha256: str
    dataset: Literal["gnomad-joint-v4.1"] = "gnomad-joint-v4.1"
    subset: Literal["joint"] = "joint"
    metric: Literal["FAF95"] = "FAF95"

    def __post_init__(self) -> None:
        has_interval = (
            self.precision_lower_bound is not None
            and self.precision_upper_bound is not None
        )
        if (
            self.query_status not in {
                "exact_match", "no_match", "unassessable", "unavailable", "conflicting"
            }
            or self.dataset != "gnomad-joint-v4.1"
            or self.subset != "joint"
            or self.metric != "FAF95"
            or self.maximum_raw_token is not None
            and not _decision_text(self.maximum_raw_token)
            or self.maximum_ancestry is not None
            and not _decision_text(self.maximum_ancestry)
            or self.precision_status not in {"qualified", "missing", "unassessable"}
            or (self.precision_lower_bound is None)
            != (self.precision_upper_bound is None)
            or (self.precision_status == "qualified") != has_interval
            or has_interval and self.precision_lower_bound.fraction
            > self.precision_upper_bound.fraction
            or self.matching_ancestry_ac is not None
            and not _decision_count(self.matching_ancestry_ac)
            or self.matching_ancestry_an is not None
            and not _decision_count(self.matching_ancestry_an)
            or self.joint_ac is not None and not _decision_count(self.joint_ac)
            or self.joint_an is not None and not _decision_count(self.joint_an)
            or self.exomes is not None and not isinstance(
                self.exomes, PopulationSubsetFacts
            )
            or self.genomes is not None and not isinstance(
                self.genomes, PopulationSubsetFacts
            )
            or self.exomes is not None and self.exomes.subset != "exomes"
            or self.genomes is not None and self.genomes.subset != "genomes"
            or self.discordant not in {True, False, None}
            or self.comparison_status not in {
                "applicable", "proved_not_applicable", "unknown"
            }
            or not _ordered_values(self.technical_warnings)
            or not _decision_digest(self.population_manifest_sha256)
            or not _decision_digest(self.precision_manifest_sha256)
        ):
            _decision_fail()


@dataclass(frozen=True, slots=True)
class ClinicalProtectionFacts:
    query_status: Literal["matches", "no_match", "unassessable", "unavailable"]
    significance_tokens: tuple[str, ...]
    conflict_tokens: tuple[str, ...]
    curated_exact_match: bool
    curated_interval_match: bool
    curated_overlap_status: Literal["complete", "uncertain", "unavailable"]
    searched_roles: tuple[str, ...]
    clinical_manifest_sha256: str
    exception_manifest_sha256: str
    assertion_scopes: tuple[str, ...] = ()

    def __post_init__(self) -> None:
        if (
            self.query_status not in {"matches", "no_match", "unassessable", "unavailable"}
            or not _ordered_values(self.significance_tokens)
            or not _ordered_values(self.conflict_tokens)
            or type(self.curated_exact_match) is not bool
            or type(self.curated_interval_match) is not bool
            or self.curated_overlap_status not in {"complete", "uncertain", "unavailable"}
            or self.searched_roles != CLINICAL_ROLES
            or not _decision_digest(self.clinical_manifest_sha256)
            or not _decision_digest(self.exception_manifest_sha256)
            or not _ordered_values(
                self.assertion_scopes, {"aggregate", "individual"}
            )
            or self.query_status == "no_match" and (
                self.significance_tokens or self.conflict_tokens
                or self.curated_exact_match or self.curated_interval_match
            )
        ):
            _decision_fail()


@dataclass(frozen=True, slots=True, repr=False)
class PrefilterProof:
    population_query_status: str
    dataset: str
    subset: str
    metric: str
    maximum_raw_token: str | None
    maximum_ancestry: str | None
    precision_lower_bound: RationalEvidence | None
    precision_upper_bound: RationalEvidence | None
    matching_ancestry_ac: int | None
    matching_ancestry_an: int | None
    protection_query_status: str
    protected_significance_tokens: tuple[str, ...]
    conflict_tokens: tuple[str, ...]
    assertion_scopes: tuple[str, ...]
    curated_exact_match: bool
    curated_interval_match: bool

    def __post_init__(self) -> None:
        if (
            self.population_query_status not in {
                "exact_match", "no_match", "unassessable", "unavailable", "conflicting"
            }
            or self.dataset != "gnomad-joint-v4.1"
            or self.subset != "joint"
            or self.metric != "FAF95"
            or self.maximum_raw_token is not None
            and not _decision_text(self.maximum_raw_token)
            or self.maximum_ancestry is not None
            and not _decision_text(self.maximum_ancestry)
            or (self.precision_lower_bound is None)
            != (self.precision_upper_bound is None)
            or self.matching_ancestry_ac is not None
            and not _decision_count(self.matching_ancestry_ac)
            or self.matching_ancestry_an is not None
            and not _decision_count(self.matching_ancestry_an)
            or self.protection_query_status
            not in {"matches", "no_match", "unassessable", "unavailable"}
            or not _ordered_values(self.protected_significance_tokens)
            or not _ordered_values(self.conflict_tokens)
            or not _ordered_values(
                self.assertion_scopes, {"aggregate", "individual"}
            )
            or type(self.curated_exact_match) is not bool
            or type(self.curated_interval_match) is not bool
        ):
            _decision_fail()


@dataclass(frozen=True, slots=True, repr=False)
class PrefilterDecision:
    source_allele_id: str
    record_number: int
    allele_number: int
    normalized_allele_id: str | None
    disposition: Literal["retain", "retain_unknown", "exclude_common"]
    primary_reason: str
    reason_codes: tuple[str, ...]
    proof: PrefilterProof
    policy_sha256: str
    population_manifest_sha256: str
    clinical_manifest_sha256: str
    exception_manifest_sha256: str
    precision_manifest_sha256: str
    classification: ArtifactClassification = ArtifactClassification.CONTROLLED
    llm_policy: LLMPolicy = LLMPolicy.NO_LLM
    schema: str = DECISION_SCHEMA

    def __post_init__(self) -> None:
        if (
            self.schema != DECISION_SCHEMA
            or SOURCE_ALLELE_ID.fullmatch(self.source_allele_id) is None
            or not _decision_count(self.record_number)
            or self.record_number <= 0
            or not _decision_count(self.allele_number)
            or self.allele_number <= 0
            or self.normalized_allele_id is not None
            and NORMALIZED_ALLELE_ID.fullmatch(self.normalized_allele_id) is None
            or self.disposition not in DECISION_DISPOSITIONS
            or not _decision_text(self.primary_reason)
            or not self.reason_codes
            or self.reason_codes[0] != self.primary_reason
            or len(set(self.reason_codes)) != len(self.reason_codes)
            or any(not _decision_text(reason) for reason in self.reason_codes)
            or not isinstance(self.proof, PrefilterProof)
            or any(not _decision_digest(value) for value in (
                self.policy_sha256, self.population_manifest_sha256,
                self.clinical_manifest_sha256, self.exception_manifest_sha256,
                self.precision_manifest_sha256,
            ))
            or self.classification is not ArtifactClassification.CONTROLLED
            or self.llm_policy is not LLMPolicy.NO_LLM
        ):
            _decision_fail()


def _decision(
    occurrence: PrefilterOccurrence,
    population: PopulationFilterFacts,
    protection: ClinicalProtectionFacts,
    policy: CommonAlleleFilterPolicy,
    disposition: str,
    reasons: list[str],
) -> PrefilterDecision:
    if not reasons:
        _decision_fail()
    return PrefilterDecision(
        occurrence.source_allele_id,
        occurrence.record_number,
        occurrence.allele_number,
        occurrence.normalized_allele_id,
        disposition,
        reasons[0],
        tuple(dict.fromkeys(reasons)),
        PrefilterProof(
            population.query_status,
            population.dataset,
            population.subset,
            population.metric,
            population.maximum_raw_token,
            population.maximum_ancestry,
            population.precision_lower_bound,
            population.precision_upper_bound,
            population.matching_ancestry_ac,
            population.matching_ancestry_an,
            protection.query_status,
            protection.significance_tokens,
            protection.conflict_tokens,
            protection.assertion_scopes,
            protection.curated_exact_match,
            protection.curated_interval_match,
        ),
        hashlib.sha256(common_allele_filter_policy_as_json(policy)).hexdigest(),
        population.population_manifest_sha256,
        protection.clinical_manifest_sha256,
        protection.exception_manifest_sha256,
        population.precision_manifest_sha256,
    )


def _protection_reasons(
    facts: ClinicalProtectionFacts, policy: CommonAlleleFilterPolicy,
) -> tuple[str, list[str]]:
    reasons: list[str] = []
    if facts.curated_exact_match:
        reasons.append("curated_exact_allele_protection")
    if facts.curated_interval_match:
        reasons.append("curated_interval_protection")
    recognized_conflicts = set(facts.conflict_tokens).intersection(
        policy.clinical_protection.conflict_markers
    )
    if recognized_conflicts:
        reasons.append("explicit_clinical_conflict_protection")
    protected = set(facts.significance_tokens).intersection(
        policy.clinical_protection.protected_significance
    )
    if protected:
        reasons.append("pathogenic_or_likely_pathogenic_protection")
    if reasons:
        return "retain", reasons
    if facts.query_status in {"unassessable", "unavailable"}:
        return "retain_unknown", ["clinical_or_curated_lookup_incomplete"]
    if facts.curated_overlap_status in {"uncertain", "unavailable"}:
        return "retain_unknown", ["curated_interval_overlap_unresolved"]
    if facts.conflict_tokens and not set(facts.conflict_tokens) <= set(
        policy.clinical_protection.conflict_markers
    ):
        return "retain_unknown", ["clinical_conflict_marker_uninterpretable"]
    recognized = set(policy.clinical_protection.protected_significance).union(
        policy.clinical_protection.unprotected_significance
    )
    if facts.query_status == "matches" and (
        not facts.significance_tokens
        or not set(facts.significance_tokens) <= recognized
    ):
        return "retain_unknown", ["clinical_assertion_uninterpretable"]
    return "continue", ["no_protection_found_in_pinned_sources"]


def _subset_quality_reasons(subset: PopulationSubsetFacts) -> list[str]:
    reasons: list[str] = []
    if subset.noncontribution_proof is not None and (
            subset.allele_count not in (None, 0) or subset.allele_number not in (None, 0)):
        reasons.append('exomes_noncontribution_conflicts_with_variant_counts')
    if (subset.subset == 'exomes' and subset.contribution == 'proven_no_callable'
            and subset.noncontribution_proof is not None
            and subset.allele_count in (None, 0) and subset.allele_number in (None, 0)):
        return reasons
    if subset.count_status != "valid":
        reasons.append(f"{subset.subset}_counts_{subset.count_status}")
        return reasons
    if subset.allele_count > subset.allele_number:
        reasons.append(f"{subset.subset}_counts_inconsistent")
    if subset.contribution == "unknown":
        reasons.append(f"{subset.subset}_contribution_unknown")
    if subset.allele_count > 0:
        if subset.source_filters_defined is not True:
            reasons.append(f"{subset.subset}_source_filters_undefined")
        elif any(value not in {None, "PASS"} for value in subset.source_filters):
            reasons.append(f"{subset.subset}_source_filtered")
        if subset.variant_not_called is not False:
            reasons.append(f"{subset.subset}_positive_alternate_not_called_or_unknown")
    if subset.subset == "exomes" and subset.contribution == "observed":
        if subset.interval_qc_pass is not True:
            reasons.append("exome_interval_qc_failed_or_unknown")
    if subset.allele_count > 0 and subset.variant_not_called is True:
        reasons.append(f"{subset.subset}_not_called_contradicts_positive_count")
    return reasons


def _population_reasons(
    facts: PopulationFilterFacts, policy: CommonAlleleFilterPolicy,
) -> tuple[str, list[str]]:
    if facts.query_status != "exact_match":
        return "retain_unknown", [
            {
                "no_match": "gnomad_exact_match_absent",
                "unassessable": "gnomad_evidence_unassessable",
                "unavailable": "gnomad_evidence_unavailable",
                "conflicting": "gnomad_exact_matches_conflict",
            }[facts.query_status]
        ]
    reasons: list[str] = []
    if facts.maximum_ancestry not in policy.frequency.allowed_maximum_ancestries:
        reasons.append("maximum_ancestry_missing_or_unsupported")
    if facts.maximum_raw_token is None:
        reasons.append("maximum_faf95_missing")
    if facts.precision_status != "qualified":
        reasons.append("faf95_precision_evidence_missing_or_unassessable")
    if facts.matching_ancestry_ac is None or facts.matching_ancestry_an is None:
        reasons.append("matching_ancestry_counts_missing")
    elif (
        facts.matching_ancestry_ac <= 0
        or facts.matching_ancestry_ac > facts.matching_ancestry_an
    ):
        reasons.append("matching_ancestry_counts_inconsistent")
    elif facts.matching_ancestry_an < policy.frequency.minimum_matching_ancestry_an:
        reasons.append("matching_ancestry_an_below_2000")
    if facts.joint_ac is not None or facts.joint_an is not None:
        if (
            facts.joint_ac is None or facts.joint_an is None
            or facts.joint_ac > facts.joint_an
        ):
            reasons.append("joint_counts_missing_or_inconsistent")
    for subset in (facts.exomes, facts.genomes):
        if subset is None:
            reasons.append("subset_quality_evidence_missing")
        else:
            reasons.extend(_subset_quality_reasons(subset))
    if facts.exomes is not None and facts.genomes is not None:
        exomes, genomes = facts.exomes, facts.genomes
        missing_exome_exception = (exomes.count_status == 'missing'
            and exomes.contribution == 'proven_no_callable' and exomes.noncontribution_proof is not None)
        if missing_exome_exception and (
                genomes.count_status != 'valid' or genomes.contribution != 'observed'
                or genomes.allele_count is None or genomes.allele_count <= 0
                or genomes.allele_number is None or genomes.allele_number <= 0
                or facts.joint_ac != genomes.allele_count or facts.joint_an != genomes.allele_number):
            reasons.append('noncontribution_joint_genome_counts_inconsistent')
        if (
            facts.joint_ac is not None and facts.joint_an is not None
            and exomes.count_status == genomes.count_status == "valid"
            and (
                exomes.allele_count + genomes.allele_count != facts.joint_ac
                or exomes.allele_number + genomes.allele_number != facts.joint_an
            )
        ):
            reasons.append("joint_subset_counts_inconsistent")
        comparison_applicable = (
            exomes.contribution == "observed"
            and genomes.contribution == "observed"
        )
        if comparison_applicable:
            if facts.comparison_status != "applicable":
                reasons.append("discordance_comparison_missing")
            elif facts.discordant is not False:
                reasons.append("joint_technology_discordant_or_unknown")
            elif (
                facts.discordance_p is None
                or facts.discordance_p.fraction
                < policy.population_quality.discordance_minimum_p.fraction
            ):
                reasons.append("discordance_p_missing_or_below_minimum")
        elif facts.comparison_status != "proved_not_applicable":
            reasons.append("discordance_not_applicable_proof_missing")
    if facts.technical_warnings:
        reasons.append("population_warning_unrecognized")
    if reasons:
        return "retain_unknown", reasons
    lower = facts.precision_lower_bound.fraction
    upper = facts.precision_upper_bound.fraction
    if lower > policy.threshold:
        return "exclude_common", ["joint_faf95_lower_bound_above_1_over_20"]
    if upper < policy.threshold:
        return "retain", ["joint_faf95_upper_bound_below_1_over_20"]
    return "retain_unknown", ["joint_faf95_precision_touches_or_crosses_1_over_20"]


def prefilter_bypass_reasons(occurrence: PrefilterOccurrence) -> tuple[str, ...]:
    """Shared pre-evidence routing; skipped lookups cannot alter these reasons."""
    if not isinstance(occurrence, PrefilterOccurrence):
        _decision_fail()
    if occurrence.preparation_disposition != "normalized":
        return (f"preparation_{occurrence.preparation_disposition}",)
    bypass: list[str] = []
    if occurrence.contig_class != "autosomal":
        bypass.append("non_autosomal_bypass")
    if occurrence.source_variant_class not in {"snv", "insertion", "deletion"}:
        bypass.append("source_representation_class_bypass")
    if occurrence.normalized_kind not in {"snv", "insertion", "deletion"}:
        bypass.append("normalized_representation_class_bypass")
    bypass.extend(f"explicit_phase_{value}" for value in occurrence.phase_evidence)
    bypass.extend(f"representation_{value}" for value in occurrence.representation_flags)
    bypass.extend(f"native_class_{value}" for value in occurrence.native_classes)
    bypass.extend(
        f"unresolved_technical_obligation_{value}"
        for value in occurrence.technical_obligations
    )
    return tuple(bypass)


def evaluate_common_allele_filter(
    occurrence: PrefilterOccurrence,
    population_facts: PopulationFilterFacts,
    protection_facts: ClinicalProtectionFacts,
    *,
    policy: CommonAlleleFilterPolicy,
) -> PrefilterDecision:
    """Apply v1 precedence to authenticated, typed facts only."""
    if (
        not isinstance(occurrence, PrefilterOccurrence)
        or not isinstance(population_facts, PopulationFilterFacts)
        or not isinstance(protection_facts, ClinicalProtectionFacts)
        or not isinstance(policy, CommonAlleleFilterPolicy)
    ):
        _decision_fail()
    for subset in (population_facts.exomes, population_facts.genomes):
        if subset is not None and subset.noncontribution_proof is not None:
            if subset.noncontribution_proof.normalized_allele_id != occurrence.normalized_allele_id:
                _decision_fail()
    bypass = prefilter_bypass_reasons(occurrence)
    if bypass:
        return _decision(
            occurrence, population_facts, protection_facts, policy, "retain", bypass,
        )
    if occurrence.reference_equivalence != "exact":
        return _decision(
            occurrence, population_facts, protection_facts, policy,
            "retain_unknown", ["source_reference_equivalence_uncertain"],
        )

    protection_status, protection_reasons = _protection_reasons(
        protection_facts, policy,
    )
    if protection_status != "continue":
        return _decision(
            occurrence, population_facts, protection_facts, policy,
            protection_status, protection_reasons,
        )
    population_status, population_reasons = _population_reasons(
        population_facts, policy,
    )
    return _decision(
        occurrence, population_facts, protection_facts, policy,
        population_status, [*population_reasons, *protection_reasons],
    )


@dataclass(frozen=True, slots=True, repr=False)
class NormalizedAlleleClosure:
    normalized_allele_id: str
    disposition: Literal["retain", "exclude_common"]
    total_occurrences: int
    retained_occurrences: int
    unknown_occurrences: int
    excluded_occurrences: int
    first_record_number: int
    first_allele_number: int

    def __post_init__(self) -> None:
        if (
            NORMALIZED_ALLELE_ID.fullmatch(self.normalized_allele_id) is None
            or self.disposition not in {"retain", "exclude_common"}
            or not _decision_count(self.total_occurrences)
            or self.total_occurrences <= 0
            or any(not _decision_count(value) for value in (
                self.retained_occurrences, self.unknown_occurrences,
                self.excluded_occurrences,
            ))
            or self.retained_occurrences + self.unknown_occurrences
            + self.excluded_occurrences != self.total_occurrences
            or (self.disposition == "exclude_common")
            != (self.excluded_occurrences == self.total_occurrences)
            or not _decision_count(self.first_record_number)
            or self.first_record_number <= 0
            or not _decision_count(self.first_allele_number)
            or self.first_allele_number <= 0
        ):
            _decision_fail()


def _sqlite_rows(cursor, batch_size: int = 1024):
    while True:
        rows = cursor.fetchmany(batch_size)
        if not rows:
            return
        yield from rows


class NormalizedRetentionAccumulator:
    """Disk-backed all-occurrence closure across arbitrary source shards."""

    __slots__ = (
        "path", "expected_occurrences", "_connection", "_count", "_previous",
        "_sealed", "_closed",
    )

    def __init__(self, scratch_path: str | Path, *, expected_occurrences: int) -> None:
        if not _decision_count(expected_occurrences) or expected_occurrences <= 0:
            _decision_fail()
        self.path = Path(scratch_path)
        self.expected_occurrences = expected_occurrences
        descriptor = -1
        try:
            descriptor = os.open(
                self.path,
                os.O_WRONLY | os.O_CREAT | os.O_EXCL | os.O_NOFOLLOW,
                0o600,
            )
            os.fchmod(descriptor, 0o600)
            os.close(descriptor)
            descriptor = -1
            self._connection = sqlite3.connect(self.path)
            self._connection.enable_load_extension(False)
            self._connection.execute("PRAGMA foreign_keys=ON")
            self._connection.execute("PRAGMA cache_size=-8192")
            self._connection.execute("PRAGMA temp_store=FILE")
            self._connection.execute("PRAGMA mmap_size=0")
            self._connection.execute("PRAGMA trusted_schema=OFF")
            self._connection.executescript("""
                CREATE TABLE normalized_identities(
                    normalized_allele_id TEXT PRIMARY KEY,
                    first_record_number INTEGER NOT NULL,
                    first_allele_number INTEGER NOT NULL
                );
                CREATE TABLE occurrences(
                    source_allele_id TEXT PRIMARY KEY,
                    record_number INTEGER NOT NULL,
                    allele_number INTEGER NOT NULL,
                    normalized_allele_id TEXT NOT NULL
                        REFERENCES normalized_identities(normalized_allele_id),
                    disposition TEXT NOT NULL
                        CHECK(disposition IN ('retain','retain_unknown','exclude_common')),
                    UNIQUE(record_number,allele_number)
                );
                CREATE INDEX occurrence_normalized_order ON occurrences(
                    normalized_allele_id,record_number,allele_number
                );
            """)
            self._connection.commit()
        except (OSError, sqlite3.Error):
            if descriptor >= 0:
                os.close(descriptor)
            raise CommonAlleleFilterError(
                "cannot create common-allele retention state"
            ) from None
        self._count = 0
        self._previous = (0, 0)
        self._sealed = False
        self._closed = False

    def add(self, decision: PrefilterDecision) -> None:
        if (
            self._closed
            or self._sealed
            or not isinstance(decision, PrefilterDecision)
            or decision.normalized_allele_id is None
            or (decision.record_number, decision.allele_number) <= self._previous
        ):
            raise CommonAlleleFilterError(
                "common-allele retention source order is invalid"
            )
        try:
            self._connection.execute(
                "INSERT OR IGNORE INTO normalized_identities VALUES (?,?,?)",
                (
                    decision.normalized_allele_id,
                    decision.record_number,
                    decision.allele_number,
                ),
            )
            self._connection.execute(
                "INSERT INTO occurrences VALUES (?,?,?,?,?)",
                (
                    decision.source_allele_id,
                    decision.record_number,
                    decision.allele_number,
                    decision.normalized_allele_id,
                    decision.disposition,
                ),
            )
        except sqlite3.Error:
            raise CommonAlleleFilterError(
                "cannot update common-allele retention state"
            ) from None
        self._previous = (decision.record_number, decision.allele_number)
        self._count += 1
        if self._count % 1024 == 0:
            self._connection.commit()

    def seal(self) -> None:
        if self._closed or self._sealed or self._count != self.expected_occurrences:
            raise CommonAlleleFilterError(
                "common-allele retention accounting is invalid"
            )
        self._connection.commit()
        self._sealed = True

    def iter_closure(self):
        if self._closed or not self._sealed:
            raise CommonAlleleFilterError("common-allele retention state is unavailable")
        for row in _sqlite_rows(self._connection.execute("""
            SELECT i.normalized_allele_id,
                   CASE WHEN sum(o.disposition='exclude_common')=count(*)
                        THEN 'exclude_common' ELSE 'retain' END,
                   count(*),
                   sum(o.disposition='retain'),
                   sum(o.disposition='retain_unknown'),
                   sum(o.disposition='exclude_common'),
                   i.first_record_number,i.first_allele_number
            FROM normalized_identities i
            JOIN occurrences o USING(normalized_allele_id)
            GROUP BY i.normalized_allele_id
            ORDER BY i.first_record_number,i.first_allele_number,i.normalized_allele_id
        """)):
            yield NormalizedAlleleClosure(*row)

    def iter_occurrences(self, normalized_allele_id: str):
        if (
            self._closed
            or not self._sealed
            or not isinstance(normalized_allele_id, str)
            or NORMALIZED_ALLELE_ID.fullmatch(normalized_allele_id) is None
        ):
            _decision_fail()
        yield from _sqlite_rows(self._connection.execute("""
            SELECT source_allele_id,record_number,allele_number,disposition
            FROM occurrences WHERE normalized_allele_id=?
            ORDER BY record_number,allele_number
        """, (normalized_allele_id,)))

    def close(self) -> None:
        if not self._closed:
            self._connection.close()
            self._closed = True

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        incomplete = exc_type is None and not self._sealed
        self.close()
        if incomplete:
            raise CommonAlleleFilterError(
                "common-allele retention state was not sealed"
            )
        return False
