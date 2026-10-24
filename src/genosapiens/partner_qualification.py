"""Decision 7 member qualification kernel, NOT source authentication.

The consumer must replay admission, original clinical/annotation sources, typed
requirements and supplied genotypes against independently held inputs and pins.
Hashes below bind inputs, not their truth. No population decision, enumeration,
completeness claim, technical discharge, file access or callback occurs here.
"""

from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass, replace
import re
from typing import Literal, get_args

from genosapiens import candidate_admission as admission
from genosapiens.admission_links import parse_open_targets_variant_row
from genosapiens.clinical_annotation import _CONFLICT_MARKERS, _polarities, _projection, _ref
from genosapiens.disease_gene_evidence import GeneRequirements, _G2P_ALLELES, _requirements_disagree
from genosapiens.evidence_model import (
    AlleleEvidence, DiseaseGeneEvidence, EvidenceReference, HypothesisContext,
    HypothesisPhase, SourcedRecord, _Controlled, _ids, _typed,
)
from genosapiens.genotype_properties import (
    POLICY_SHA256 as GENOTYPE_POLICY_SHA256, AlleleGenotypeProperties,
    PairPhaseProperties, SuppliedGenotypeProperties, selected_source_calls,
)
from genosapiens.genotype_model import (
    MONOALLELIC_HOMOZYGOTE_REASON, ModelDecisionGenotype, requires_monoallelic_homozygote_evidence,
)
from genosapiens.model_requirements import (
    RequirementAssessment, RequirementField, _hash, _identifier_match,
)


SCHEMA = "genosapiens.partner-qualification"
QualificationState = Literal["qualifying", "nonqualifying", "unresolved"]
_POSITIVE = {"pathogenic", "likely pathogenic", "pathogenic/likely pathogenic"}
_NEGATIVE = {"benign", "likely benign", "benign/likely benign"}
_MODELS = {
    "dominant", "homozygous_recessive", "compound_heterozygous", "x_linked",
    "incomplete", "biallelic", "mitochondrial", "mosaic", "structural", "repeat",
    "source-requirement-context-only",
}


class PartnerQualificationError(ValueError):
    """Static diagnostics; malformed bindings are not biological unknowns."""


@dataclass(frozen=True, slots=True, repr=False)
class ClinicalStrengthPolicy(_Controlled):
    """Exact publisher spellings, not stars or a numerical confidence threshold.

    Empty vocabularies accept nothing. Null never satisfies a vocabulary.
    Review and confidence are separately explicit even when the adapter projects
    both from the same source field. Significance is existential within a row;
    all opposing signs and affirmative conflict markers are still retained.
    """

    accepted_significances: tuple[str, ...]
    accepted_review_statuses: tuple[str, ...]
    accepted_confidences: tuple[str, ...]

    def __post_init__(self):
        _Controlled.__post_init__(self)
        for values in (self.accepted_significances, self.accepted_review_statuses, self.accepted_confidences):
            _ids(values)


@dataclass(frozen=True, slots=True, repr=False)
class ClinicalQualificationPolicy(_Controlled):
    """Caller-selected, independently pinned settings; no production defaults.

    The current condition-linked reader is EVA 26.06. ClinVar VCF aggregates
    remain raw evidence, never a condition/significance Cartesian product.
    """

    policy_id: str
    source: SourcedRecord
    positive: ClinicalStrengthPolicy
    negative: ClinicalStrengthPolicy
    adapter: Literal["open-targets-26.06-eva"] = "open-targets-26.06-eva"

    def __post_init__(self):
        _Controlled.__post_init__(self)
        _ids((self.policy_id,))
        if self.source.public_manifest_sha256 is None or self.source.source_fields is None:
            raise PartnerQualificationError("source-bound clinical policy required")
        for policy, allowed in ((self.positive, _POSITIVE), (self.negative, _NEGATIVE)):
            if any(value.lower().replace("_", " ") not in allowed for value in policy.accepted_significances):
                raise PartnerQualificationError("invalid clinical significance vocabulary")

    @property
    def policy_sha256(self) -> str:
        return _hash({"schema": SCHEMA, "policy": self})


@dataclass(frozen=True, slots=True, repr=False)
class ClinicalQualificationBasis(_Controlled):
    source_ref: EvidenceReference
    source: SourcedRecord
    condition_id: str | None
    target_id: str | None
    applicability: Literal["exact", "equivalent", "broader", "unresolved", "other_target", "aggregate"]
    positive: bool
    negative: bool
    self_conflict: bool
    relation_refs: tuple[EvidenceReference, ...]
    conflict_refs: tuple[EvidenceReference, ...]
    superseded_by: tuple[EvidenceReference, ...]
    dependency_group: str
    reasons: tuple[str, ...]


@dataclass(frozen=True, slots=True, repr=False)
class PartnerQualificationAssessment(_Controlled):
    context_id: str
    allele_id: str
    native_identity_id: str
    state: QualificationState
    clinical: tuple[ClinicalQualificationBasis, ...]
    clinical_support_refs: tuple[EvidenceReference, ...]
    clinical_contradiction_refs: tuple[EvidenceReference, ...]
    clinical_conflict_refs: tuple[EvidenceReference, ...]
    precedence_refs: tuple[EvidenceReference, ...]
    source_refs: tuple[SourcedRecord, ...]
    evidence_refs: tuple[EvidenceReference, ...]
    requirement_assessments: tuple[RequirementAssessment, ...]
    genotype: AlleleGenotypeProperties
    phase: PairPhaseProperties | None
    reasons: tuple[str, ...]
    binding: Mapping[str, object]
    authentication: Literal["consumer_must_replay_sources"] = "consumer_must_replay_sources"

    def as_bytes(self) -> bytes:
        """Controlled/no-LLM audit only; not a replay or no-partner certificate."""
        return admission._bytes({"schema": SCHEMA, "classification": "controlled",
                                 "llm_policy": "no-llm", "assessment": self})


def _clinical(context, allele, policy, evidence, *, prepared=None):
    sources = admission._unique(c.source for c in allele.clinical or ())
    if not sources:
        return ()
    reached = (admission._reach(context.native_disease_id or context.disease_id, admission._relation_edges(evidence))
               if prepared is None else prepared.clinical_reach(context))
    records = {}
    variants = {}
    for source in sources:
        ref = _ref(source)
        if ref in records and records[ref] != source:
            raise ValueError
        records[ref] = source
        if source.source_dataset == "Open Targets variant" and source.source_release == "26.06":
            row = source.source_fields["record"]
            parse_open_targets_variant_row(row)
            variants[source.public_manifest_sha256, row["variantId"]] = source
    bases, tokens = [], []
    for source in sources:
        ref = _ref(source)
        releases = {"ClinVar VCF": "2026-08-22", "Open Targets variant": "26.06", "Open Targets evidence_eva": "26.06"}
        if releases.get(source.source_dataset) != source.source_release:
            bases.append(ClinicalQualificationBasis(ref, source, None, None, "unresolved", False, False,
                False, (), (), (), "", ("unrecognized-clinical-source",)))
            tokens.append(set())
            continue
        original = _projection(source)
        row = source.source_fields["record"]
        condition = next(iter(original.condition_ids or ()), None)
        distance, relation_refs = reached.get(condition, (None, ()))
        applicability = "unresolved"
        reasons = []
        if not original.explicit_condition_link:
            applicability = "aggregate" if source.source_dataset == "ClinVar VCF" else "unresolved"
            reasons.append("no-source-linked-condition-assertion; not-benign")
        elif (context.locus_id is not None or _identifier_match(original.target_id, context.gene_id,
                stable=context.gene_id is not None and "." not in context.gene_id) is not True):
            applicability = "other_target" if original.target_id is not None else "unresolved"
            reasons.append("source-target-not-established-for-native-context")
        elif distance is not None:
            applicability = ("exact" if condition == (context.native_disease_id or context.disease_id) else
                             "equivalent" if distance == 0 else "broader")
        else:
            reasons.append("native-condition-not-linked; mapped-view-not-equivalence")
        if original.assertion_id is not None:
            # Exact allele / alias / REF binding is authenticated by main replay,
            # not by decoding variantId or by trusting enriched ClinicalEvidence.
            if (source.public_manifest_sha256, row["variantId"]) not in variants:
                raise ValueError
            if any(row.get(field) for field in ("allelicRequirements", "alleleOrigins")):
                # Preserve publisher restrictions until a source-specific
                # interpretation exists; they are not wildcard applicability.
                applicability = "unresolved"
                reasons.append("source-assertion-model-restriction-interpretation-unavailable")
        signs = _polarities(original.significance)
        conflict = len(signs) == 2
        fields = row["info"] if source.source_dataset == "ClinVar VCF" else row
        conflict |= source.source_dataset == "ClinVar VCF" and fields.get("CLNSIGCONF") not in (None, ".", "")
        for field, vocabulary in _CONFLICT_MARKERS.get(source.source_dataset, {}).items():
            values = fields.get(field)
            if type(values) is str:
                values = values.split("|") if field == "CLNSIG" else (values,)
            conflict |= any(v is not None and v.lower() in vocabulary for v in values or ())
        strengths = tuple(original.explicit_condition_link is True
            and original.review_status in strength.accepted_review_statuses
            and original.confidence in strength.accepted_confidences
            and any(v in strength.accepted_significances for v in original.significance or ())
            for strength in (policy.positive, policy.negative))
        if not any(strengths):
            reasons.append("clinical-strength-policy-not-met; not-negative-by-default")
        if conflict:
            reasons.append("affirmative-source-self-conflict")
        bases.append(ClinicalQualificationBasis(ref, source, condition, original.target_id, applicability,
            *strengths, bool(conflict), relation_refs, (ref,) if conflict else (), (), "", tuple(reasons)))
        # Same-family / shared-publication dependence follows the existing source
        # rule. It never deletes opposing assertions or creates independent votes.
        family = "clinvar" if source.source_dataset in ("ClinVar VCF", "Open Targets evidence_eva") else source.source_dataset
        group = {"family:" + family}
        for publication in row.get("literature") or ():
            if publication and (match := re.fullmatch(r"(?:PMID:)?([1-9][0-9]*)", publication, re.I)):
                group.add("PMID:" + match[1])
        tokens.append(group)
    specific = tuple(b.source_ref for b in bases if b.applicability in ("exact", "equivalent")
                     and (b.positive or b.negative or b.self_conflict))
    for index, basis in enumerate(bases):
        if basis.applicability == "broader" and specific:
            bases[index] = replace(basis, superseded_by=specific,
                reasons=(*basis.reasons, "umbrella-evidence-superseded"))
    pending = set(range(len(bases)))
    while pending:
        members, queue = set(), [min(pending)]
        while queue:
            index = queue.pop()
            if index in members:
                continue
            members.add(index)
            queue.extend(j for j in pending - members if tokens[index] & tokens[j])
        group = "clinical-dependency-" + _hash({
            "tokens": sorted({t for i in members for t in tokens[i]}),
            "source_refs": sorted(
                (bases[i].source_ref for i in members),
                key=lambda ref: (ref.evidence_id, ref.record_id, ref.public_manifest_sha256),
            ),
        })
        for index in members:
            basis = bases[index]
            conflicts = list(basis.conflict_refs)
            if basis.applicability in ("exact", "equivalent"):
                conflicts.extend(peer.source_ref for peer in bases if peer.applicability in ("exact", "equivalent")
                    and (basis.positive and peer.negative or basis.negative and peer.positive))
            bases[index] = replace(basis, dependency_group=group, conflict_refs=admission._unique(conflicts))
        pending -= members
    return tuple(bases)


@dataclass(frozen=True, slots=True, repr=False)
class QualificationBasisDecision:
    state: str
    has_references: bool


@dataclass(frozen=True, slots=True, repr=False)
class QualificationRequirementDecision:
    requirement: GeneRequirements
    field: str
    state: str
    has_rule: bool
    bases: tuple[QualificationBasisDecision, ...]


@dataclass(frozen=True, slots=True, repr=False)
class QualificationDecisionInputs:
    model: str
    context_requirements: tuple[GeneRequirements, ...]
    genotype: ModelDecisionGenotype
    read_validation_required: bool
    requirements: tuple[QualificationRequirementDecision, ...]
    phase: HypothesisPhase | None
    support: bool
    contradiction: bool
    conflicts: bool


def qualification_decision_inputs(context, genotype, requirements, phase, *, support, contradiction, conflicts):
    """Exact projection of every state/reason dependency; audit evidence stays separate."""
    return QualificationDecisionInputs(context.model, context.requirements,
        ModelDecisionGenotype(genotype.alternate_present, genotype.minimum_alt_copies,
            genotype.maximum_alt_copies, genotype.reported_ploidies, genotype.conflicting_source_genotypes,
            genotype.screened_low_balance_candidate, bool(genotype.phase_conflict_bases)),
        genotype.read_validation_required,
        tuple(QualificationRequirementDecision(r.requirement, r.field, r.state, r.rule is not None,
            tuple(QualificationBasisDecision(b.state, bool(b.source_refs or b.evidence_refs)) for b in r.bases))
            for r in requirements), None if phase is None else phase.phase,
        bool(support), bool(contradiction), bool(conflicts))


def evaluate_qualification_decision(inputs):
    """The sole state/reason kernel; source identity construction is separate.

    Phase contributes only its categorical relation here. Full pair/source bases
    remain mandatory in the returned assessment. This separation allows bounded
    callers to evaluate equivalent decision inputs once without conflating their
    distinct evidence, bindings or members.
    """
    genotype, requirements, phase = inputs.genotype, inputs.requirements, inputs.phase
    support, contradiction, conflicts = inputs.support, inputs.contradiction, inputs.conflicts
    reasons, negative, unavailable, material_conflict = [], bool(contradiction), False, bool(conflicts)
    if support:
        reasons.append("applicable-clinical-positive-policy-met")
    if contradiction:
        reasons.append("applicable-clinical-negative-policy-met")
    if not support and not contradiction:
        reasons.append("no-qualifying-clinical-basis; absence-not-benign")
    if not inputs.context_requirements or inputs.model not in _MODELS:
        unavailable = True
        reasons.append("model-interpretation-unavailable")
    if genotype.conflicting_source_genotypes:
        material_conflict = True
        reasons.append("material-source-genotype-conflict")
    elif genotype.screened_low_balance_candidate:
        unavailable = True
        reasons.append("screen-candidate-not-confirmed-carriage; source-GT-not-rewritten")
    elif genotype.alternate_present is False:
        negative = True
        reasons.append("evaluated-source-GT-excludes-member-carriage")
    elif genotype.alternate_present is not True or genotype.minimum_alt_copies is None or genotype.minimum_alt_copies < 1:
        unavailable = True
        reasons.append("member-carriage-unresolved")
    counts = {r.required_allele_count for r in inputs.context_requirements}
    if not counts or None in counts:
        reasons.append("required-model-count-unknown; no-inheritance-count-inference")
    if inputs.model == "homozygous_recessive":
        if len(counts) != 1 or None in counts:
            unavailable = True
            reasons.append("homozygous-copy-requirement-unavailable")
        else:
            count, = counts
            if (genotype.maximum_alt_copies is not None and genotype.maximum_alt_copies < count
                    and not genotype.screened_low_balance_candidate):
                negative = True
                reasons.append("evaluated-copies-contradict-explicit-homozygous-model")
            elif (genotype.minimum_alt_copies != count or genotype.maximum_alt_copies != count
                  or genotype.reported_ploidies != (count,)):
                unavailable = True
                reasons.append("homozygous-copy-and-ploidy-match-unresolved")
    trans_required = inputs.model == "compound_heterozygous"
    if trans_required:
        if genotype.phase_conflict_bases or phase is HypothesisPhase.CONFLICTING:
            material_conflict = True
            reasons.append("material-required-phase-conflict")
        elif phase is HypothesisPhase.IN_CIS:
            negative = True
            reasons.append("established-cis-contradicts-compound-heterozygous-model")
        elif phase is None or phase is HypothesisPhase.UNKNOWN:
            reasons.append("ordinary-unknown-phase-retained; not-trans-or-veto")
    for result in requirements:
        value = getattr(result.requirement, result.field)
        if result.field in ("inheritance", "required_allele_count") or value is None:
            continue
        if result.field == "allelic_requirement":
            kind = result.requirement.allelic_requirement_kind
            if requires_monoallelic_homozygote_evidence(result.requirement, genotype):
                unavailable = True
                reasons.append(MONOALLELIC_HOMOZYGOTE_REASON)
                continue
            if kind == "literal_allelic_model" and value in (
                "monoallelic_X_hemizygous", "monoallelic_Y_hemizygous", "monoallelic_X_heterozygous",
            ):
                count = result.requirement.required_allele_count
                ploidy = 2 if value == "monoallelic_X_heterozygous" else 1
                if count is None or genotype.screened_low_balance_candidate:
                    unavailable = True
                    reasons.append("source-zygosity-requirement-unresolved")
                elif (genotype.maximum_alt_copies is not None and genotype.maximum_alt_copies < count
                      or genotype.minimum_alt_copies is not None and genotype.minimum_alt_copies > count
                      or genotype.reported_ploidies and ploidy not in genotype.reported_ploidies):
                    negative = True
                    reasons.append("evaluated-copies-or-ploidy-contradict-source-zygosity")
                elif (genotype.minimum_alt_copies != count or genotype.maximum_alt_copies != count
                      or genotype.reported_ploidies != (ploidy,)):
                    unavailable = True
                    reasons.append("source-zygosity-requirement-unresolved")
                continue
            if kind == "mode_of_inheritance" or kind == "literal_allelic_model" and value in _G2P_ALLELES:
                reasons.append("source-allelic-label-retained; whole-model-check-not-member-copy-count")
                continue
        if result.state == "contradicted" and result.has_rule and result.bases and all(
            b.state == "contradicted" and b.has_references for b in result.bases
        ):
            negative = True
            reasons.append("applicable-allele-property-contradiction:" + result.field)
        elif result.state != "supported":
            unavailable = True
            reasons.append("required-allele-property-unresolved:" + result.field)
            if {b.state for b in result.bases} >= {"supported", "contradicted"}:
                material_conflict = True
                reasons.append("material-required-property-conflict:" + result.field)
    if inputs.read_validation_required:
        reasons.append("mandatory-read-validation-retained; not-discharged-or-negative-by-itself")
    if material_conflict:
        state = "unresolved"
        reasons.append("material-qualification-conflict-after-precedence")
    elif negative:
        state = "nonqualifying"
    elif support and not unavailable:
        state = "qualifying"
    else:
        state = "unresolved"
    return state, tuple(sorted(set(reasons)))


def _qualification_decision(context, genotype, requirements, phase, *, support, contradiction, conflicts):
    return evaluate_qualification_decision(qualification_decision_inputs(
        context, genotype, requirements, phase, support=support, contradiction=contradiction, conflicts=conflicts))


def assess_partner_qualification(
    context: HypothesisContext, allele: AlleleEvidence, *, clinical_policy: ClinicalQualificationPolicy,
    requirement_assessments: tuple[RequirementAssessment, ...], genotype_properties: SuppliedGenotypeProperties,
    phase: PairPhaseProperties | None = None, disease_gene_evidence: tuple[DiseaseGeneEvidence, ...] = (),
    prepared=None,
) -> PartnerQualificationAssessment:
    """Assess a member, optionally in a particular evaluated pair, without I/O.

    Only context.model == 'compound_heterozygous' requires trans here; there is
    no caller trans_required flag. Unknown phase is not a veto. A heterozygous
    member of a general biallelic/incomplete context need not itself carry two
    copies. Explicit homozygous_recessive models instead use the source-required
    count and evaluated copy bounds. Inheritance labels supply no copy counts.
    """
    try:
        from genosapiens.genotype_evaluation import prepare_genotype_evaluation
        if type(clinical_policy) is not ClinicalQualificationPolicy:
            raise ValueError
        if prepared is None:
            prepared = prepare_genotype_evaluation(context=context, alleles=(allele,),
                genotypes=genotype_properties, gene_evidence=disease_gene_evidence,
                policies={'clinical': clinical_policy}, pairs=())
        prepared.require(context=context, alleles=(allele,), genotypes=genotype_properties,
            gene_evidence=disease_gene_evidence, policies={'clinical': clinical_policy})
        cached_requirements = prepared.requirement_hashes.get(id(requirement_assessments))
        if cached_requirements is not None and cached_requirements[0] is requirement_assessments:
            requirements = requirement_assessments
            requirement_hash = cached_requirements[1]
        else:
            requirements = _typed(requirement_assessments, tuple[RequirementAssessment, ...])
            requirement_hash = _hash(requirements)
        properties, evidence = prepared.genotypes, prepared.gene_evidence
        phase = _typed(phase, PairPhaseProperties | None)
        prepared.validate_sources(context)
        inputs = prepared.context_allele_hash(context, allele)
        expected = {(i, field) for i in range(len(context.requirements)) for field in get_args(RequirementField)}
        if (len(requirements) != len(expected) or {(r.requirement_index, r.field) for r in requirements} != expected
                or len({r.policy_sha256 for r in requirements}) > 1):
            raise ValueError
        for result in requirements:
            req = context.requirements[result.requirement_index]
            refs = tuple(r for r in context.evidence_refs if r.record_id in req.record_ids)
            if ((result.context_id, result.allele_id, result.native_identity_id, result.input_sha256)
                    != (context.context_id, allele.allele_id, allele.native_identity_id, inputs)
                    or result.requirement != req or set(result.requirement_refs) != set(refs)):
                raise ValueError
        genotype = prepared.properties[allele.allele_id]
        native = next(i for i in allele.identities if i.identity_id == allele.native_identity_id)
        selected_calls = prepared.calls_for(allele)
        source_ids = {c.source.source_fields["source_allele_id"] for c in selected_calls}
        if (genotype.native_identity_id != allele.native_identity_id
                or native.reference_manifest_sha256 != properties.reference_manifest_sha256
                or not selected_calls or source_ids != set(genotype.source_allele_ids)):
            raise ValueError
        if phase is not None and (prepared.phases.get(tuple(sorted(phase.allele_ids))) != phase
                                  or allele.allele_id not in phase.allele_ids):
            raise ValueError
        clinical_key = (prepared.context_key(context), allele.allele_id)
        if clinical_key not in prepared.clinical:
            prepared.clinical[clinical_key] = _clinical(context, allele, clinical_policy, evidence, prepared=prepared)
        clinical = prepared.clinical[clinical_key]
        applicable = tuple(b for b in clinical if b.applicability in ("exact", "equivalent") and not b.superseded_by)
        support = admission._unique(b.source_ref for b in applicable if b.positive)
        contradiction = admission._unique(b.source_ref for b in applicable if b.negative)
        conflicts = admission._unique(r for b in applicable for r in b.conflict_refs)
        state, reasons = _qualification_decision(context, genotype, requirements, phase,
            support=support, contradiction=contradiction, conflicts=conflicts)
        binding = {"schema": SCHEMA, "context_allele_sha256": inputs,
            "clinical_policy_sha256": clinical_policy.policy_sha256,
            "requirement_assessments_sha256": requirement_hash,
            **prepared.genotype_binding_fields(), "phase_sha256": _hash(phase),
            "disease_gene_evidence_sha256": prepared.gene_evidence_sha256, "input_sha256": properties.input_sha256,
            "inventory_sha256": properties.inventory_sha256,
            "reference_manifest_sha256": properties.reference_manifest_sha256,
            "scope": "member-qualification-only; mandatory-upstream-source-replay; no-population-policy"}
        return PartnerQualificationAssessment(context.context_id, allele.allele_id, allele.native_identity_id, state,
            clinical, support, contradiction, conflicts, admission._unique(r for b in clinical for r in b.superseded_by),
            admission._unique((clinical_policy.source, native.source, *(c.source for c in selected_calls),
                *(b.source for b in clinical), *(s for r in requirements for b in r.bases for s in b.source_refs))),
            admission._unique((*context.evidence_refs, *(b.source_ref for b in clinical),
                *(r for b in clinical for r in b.relation_refs),
                *(ref for r in requirements for b in r.bases for ref in b.evidence_refs))),
            requirements, genotype, phase, reasons, binding)
    except (ValueError, TypeError, KeyError, IndexError, AttributeError, RuntimeError, StopIteration, OverflowError, RecursionError):
        raise PartnerQualificationError("invalid partner qualification inputs or source bindings") from None
