"""Replay supplied evidence before constructing context-specific research models.

This development API is not a release freeze, reference acceptance, score ledger,
or submission gate. Input pins come from the independently selected candidate
stack, never from a persisted result offered for consumption.
"""

from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass, fields, replace
import hashlib
from collections import Counter
from itertools import combinations, islice
from pathlib import Path
from typing import Literal

from genosapiens import candidate_admission as admission
from genosapiens.admission_inventory import SuppliedInputs
from genosapiens.allele_identity import MitochondrialAllele, MosaicAllele, RepeatAllele, StructuralAllele
from genosapiens.clinical_annotation import ClinicalAnnotationResult, annotate_clinical
from genosapiens.evidence_model import (
    AlleleEvidence, AlleleReference, Completeness, DiseaseGeneEvidence, DiseaseProfile,
    EnumerationStatus, GenotypeHypothesis, HypothesisContext, HypothesisPhase, IdentityStage,
    ScoringView, SourcedRecord, _Controlled, _exact, _ids, _typed,
)
from genosapiens.genotype_properties import SuppliedGenotypeProperties, _allele_event, assess_supplied_genotypes
from genosapiens import longtr_admission as longtr
from genosapiens.genotype_model import GenotypeModelAssessment, ModelLocusRegion, assess_genotype_model
from genosapiens.evidence_serialization import _plain
from genosapiens.hpo_ic import HpoICStore
from genosapiens.model_requirements import ModelRequirementPolicy, assess_model_requirements
from genosapiens.partner_qualification import (
    ClinicalQualificationPolicy, PartnerQualificationAssessment, assess_partner_qualification,
)
from genosapiens.partner_relevance import ContextRelevanceAssessment, assess_context_relevance
from genosapiens import population_annotation as population, population_policy as frequency


SCHEMA = "genosapiens.genotype-hypotheses"
EVALUABLE_DOMAIN_RULES = {
    "schema": "genosapiens.source-evaluable-partner-domain",
    "inventory": "complete-original-source-accounting-with-explicit-assessment-dispositions",
    "members": "source-supported-assessed-context-target-associations-with-designated-sample-properties",
    "pairs": "all-distinct-source-events; retain-native-and-cross-file-members",
    "pruning": "fixed-source-screening-before-contexts;no-within-assessed-context-score-or-qualification-pruning",
    "unresolved": "retain-all-source-slots-in-broader-coverage; never-exclusion-proof",
    "selected": "independent-identity-carriage-and-whole-model-checks-required",
    "recessive_singleton_shape": "source-recessive-context-and-exact-diploid-two-copy-GT-select-homozygous-model;unknown-required-count-and-completeness-remain-unknown",
    "incomplete_recessive_singleton": "eligible-supported-member-with-incomplete-disclosure-and-no-partner-proof",
}
_EXECUTABLES = (
    "genotype_hypotheses.py", "genotype_properties.py", "genotype_evaluation.py", "genotype_model.py", "model_requirements.py",
    "partner_qualification.py", "partner_relevance.py", "clinical_annotation.py",
    "population_annotation.py", "population_policy.py",
    "evidence_model.py", "evidence_serialization.py", "longtr_admission.py",
)


class GenotypeHypothesesError(ValueError):
    """Static diagnostics; invalid input is not a smaller successful universe."""


@dataclass(frozen=True, slots=True, repr=False)
class HypothesisPolicy(_Controlled):
    requirements: ModelRequirementPolicy
    clinical: ClinicalQualificationPolicy
    incompleteness_penalty: str
    loci: tuple[ModelLocusRegion, ...] = ()
    maximum_pairs_per_context: int | None = None

    def __post_init__(self):
        _Controlled.__post_init__(self)
        if _exact(self.incompleteness_penalty) <= 0 or (
            self.maximum_pairs_per_context is not None and self.maximum_pairs_per_context < 0
        ):
            raise GenotypeHypothesesError("invalid hypothesis policy")


@dataclass(frozen=True, slots=True, repr=False)
class HypothesisPopulationInputs:
    """Independent resource and policy pins; never a persisted exclusion label.

    Sources are queried once for the complete allele inventory. Thresholds and
    exceptions are evaluated only after context-specific hypotheses exist.
    """

    joint: population.GnomadJointReference
    all_sites: population.GnomadAllSitesReference
    reference: population.PopulationReferenceCompatibility
    runtime_image: str
    limits: population.PopulationQueryLimits
    policy_bytes: bytes
    policy_sha256: str
    public_curation_bytes: Mapping[str, bytes]


@dataclass(frozen=True, slots=True, repr=False)
class RepresentationEligibilityAssessment(_Controlled):
    """Negative CSV-representation metadata, never a partner-accounting input."""

    allele_id: str
    native_identity_id: str
    status: Literal["audit-only"] = "audit-only"
    reason: Literal["identity-not-faithfully-representable"] = "identity-not-faithfully-representable"
    scope: Literal["competitive-csv-representation"] = "competitive-csv-representation"

    def __post_init__(self):
        _Controlled.__post_init__(self)
        _ids((self.allele_id, self.native_identity_id))


@dataclass(frozen=True, slots=True, repr=False)
class PartnerSearchAssessment(_Controlled):
    context_id: str
    relevance: ContextRelevanceAssessment
    member_assessments: tuple[PartnerQualificationAssessment, ...]
    pair_assessments: tuple[tuple[PartnerQualificationAssessment, PartnerQualificationAssessment], ...]
    pair_models: tuple[GenotypeModelAssessment, ...]
    total_distinct_pairs: int
    evaluated_pairs: int
    traversal: Literal["complete", "truncated"]
    unresolved_source_allele_ids: tuple[str, ...]
    unresolved_relevance_source_allele_ids: tuple[str, ...]
    enumeration_status: EnumerationStatus
    limits: tuple[str, ...]
    source_inventory_status: Literal["complete", "assessed-domain-complete"]
    evaluable_allele_ids: tuple[str, ...]
    evaluable_domain_binding: SourcedRecord


@dataclass(frozen=True, slots=True, repr=False)
class ConstructedHypothesis(_Controlled):
    hypothesis: GenotypeHypothesis
    admission_context_id: str
    members: tuple[PartnerQualificationAssessment, ...]
    model: GenotypeModelAssessment
    conditional_singleton_search_satisfied: bool
    incompleteness_penalty: str | None
    read_validation_required: bool
    limits: tuple[str, ...]


@dataclass(frozen=True, slots=True, repr=False)
class HypothesisGeneration(_Controlled):
    admission: admission.CandidateAdmission
    clinical: ClinicalAnnotationResult | None
    genotypes: SuppliedGenotypeProperties
    searches: tuple[PartnerSearchAssessment, ...]
    hypotheses: tuple[ConstructedHypothesis, ...]
    binding: SourcedRecord
    population_annotation: population.PopulationAnnotationResult | None = None
    population_assessments: tuple[frequency.PopulationPolicyResult, ...] = ()
    representation_eligibility: tuple[RepresentationEligibilityAssessment, ...] = ()

    def __post_init__(self):
        _Controlled.__post_init__(self)
        alleles = {allele.allele_id: allele for allele in self.alleles}
        assessments = self.representation_eligibility
        if (assessments != tuple(sorted(assessments, key=lambda item: item.allele_id))
                or len({item.allele_id for item in assessments}) != len(assessments)
                or any(item.allele_id not in alleles
                       or item.native_identity_id != alleles[item.allele_id].native_identity_id
                       for item in assessments)):
            raise GenotypeHypothesesError("invalid representation eligibility binding")

    @property
    def alleles(self) -> tuple[AlleleEvidence, ...]:
        """Complete evidence from the latest performed stage, including empty output."""
        if self.clinical is not None:
            return self.clinical.alleles
        if self.population_annotation is not None:
            return self.population_annotation.alleles
        return self.admission.inventory.alleles

    def as_bytes(self) -> bytes:
        return admission._bytes({"schema": SCHEMA, "classification": "controlled", "llm_policy": "no-llm", "result": self})


@dataclass(frozen=True, slots=True, repr=False)
class HypothesisReplay:
    """Complete independently pinned inputs for exact hypothesis replay."""

    inputs: SuppliedInputs
    sources: admission.AdmissionSources
    view: ScoringView
    profiles: DiseaseProfile | None
    genes: tuple[DiseaseGeneEvidence, ...] | None
    ic_store: HpoICStore | None
    admission_policy_sha256: str
    sample_name: str
    policy: HypothesisPolicy
    expected_policy_sha256: str
    expected_result: HypothesisGeneration | None = None
    population_inputs: HypothesisPopulationInputs | None = None
    representation_eligibility: tuple[RepresentationEligibilityAssessment, ...] = ()


def hypothesis_policy_bytes(policy: HypothesisPolicy, *, admission_policy_sha256: str) -> bytes:
    """Pin executable rules and all independent public policy/source projections.

    As with the existing public gene/IC stores, independently held source-bound
    policy projections are inputs, not imported from a purported result. Their
    production provenance and empirical selection belong to the candidate stack.
    This pin alone never authorizes controlled use.
    """
    try:
        policy = _typed(policy, HypothesisPolicy)
        if len(admission_policy_sha256) != 64 or any(c not in "0123456789abcdef" for c in admission_policy_sha256):
            raise ValueError
        return admission._bytes({"schema": SCHEMA, "policy": policy,
            "evaluable_domain_rules": EVALUABLE_DOMAIN_RULES,
            "admission_policy_sha256": admission_policy_sha256,
            "executables": {name: hashlib.sha256(Path(__file__).with_name(name).read_bytes()).hexdigest()
                            for name in _EXECUTABLES}})
    except (ValueError, TypeError, OSError, AttributeError):
        raise GenotypeHypothesesError("invalid hypothesis policy binding") from None


def _context(seed, model):
    return replace(seed, model=model, context_id="genotype-context-" + admission._hash([seed.context_id, model]))


_PAIR_SEAL = object()


@dataclass(frozen=True, slots=True, repr=False)
class _PairUniverse:
    total: int
    pairs: tuple[tuple[str, str], ...]
    owners: tuple[AlleleEvidence, ...]
    limit: int | None
    seal: object

    def __iter__(self):
        yield self.total
        yield self.pairs

    def require(self, keys, alleles, limit):
        if (self.seal is not _PAIR_SEAL or self.limit != limit
                or len(self.owners) != len(keys)
                or any(owner is not alleles[key] for owner, key in zip(self.owners, keys, strict=True))):
            raise GenotypeHypothesesError("complete pair owner changed")


def _pairs(keys, alleles, limit):
    owners = tuple(alleles[key] for key in keys)
    keys = [k for k in keys if not longtr.is_parent(alleles[k])]
    events = {key: _allele_event(alleles[key]) for key in keys}
    counts = Counter(events.values())
    total = len(keys) * (len(keys) - 1) // 2 - sum(n * (n - 1) // 2 for n in counts.values())
    pairs = (p for p in combinations(keys, 2) if events[p[0]] != events[p[1]])
    return _PairUniverse(total, tuple(pairs if limit is None else islice(pairs, limit)), owners, limit, _PAIR_SEAL)


def _model(seed, allele, prop):
    native = admission._native(allele).identity
    if isinstance(native, RepeatAllele) and not longtr.is_projection(allele):
        return "repeat"
    if isinstance(native, MitochondrialAllele) or seed.inheritance == "HP:0001427":
        return "mitochondrial"
    if isinstance(native, MosaicAllele):
        return "mosaic"
    counts = {r.required_allele_count for r in seed.requirements}
    # Model shape and required count are independent. A source recessive context
    # plus two observed diploid copies supports a homozygous model even when
    # its required count is unknown; the model checker preserves that unknown.
    if counts == {2} or seed.inheritance == "HP:0000007":
        return "homozygous_recessive" if (prop.minimum_alt_copies, prop.maximum_alt_copies, prop.reported_ploidies) == (2, 2, (2,)) else "incomplete"
    if seed.inheritance in ("HP:0001417", "HP:0001419", "HP:0001423",
                            "g2p:monoallelic_X_hemizygous", "g2p:monoallelic_X_heterozygous"):
        return "x_linked"
    if seed.inheritance in ("HP:0001450", "g2p:monoallelic_Y_hemizygous"):
        return "y_linked"
    if seed.inheritance in ("HP:0000006", "g2p:monoallelic_PAR"):
        return "dominant"
    return "structural" if isinstance(native, StructuralAllele) else "source-requirement-context-only"


def _requires_pairs(context):
    counts = {r.required_allele_count for r in context.requirements}
    return counts == {2} or counts != {1} and context.inheritance in ("HP:0000007", "HP:0001419", "g2p:biallelic_PAR")


def consuming_model_contexts(seed, allele, properties):
    """Exact singleton and potential pair contexts before rich construction.

    Presence/qualification can later rule out construction; retaining those
    alternatives here is conservative. LongTR parents never consume a pair.
    """
    single = _context(seed, _model(seed, allele, properties))
    return ((single, _context(seed, 'compound_heterozygous'))
            if _requires_pairs(seed) and not longtr.is_parent(allele) else (single,))


def _record(kind, value):
    digest = admission._hash(value)
    return SourcedRecord(SCHEMA, "current", kind + "-" + digest, digest, None, {"kind": kind, "value": _plain(value)})


def qualify_whole_model(state, reasons, model_state):
    """Apply whole-model precedence without changing the member's audit facts."""
    if model_state == "unresolved" and state == "qualifying":
        state = "unresolved"
    elif model_state == "contradicted" and not any("material-" in reason for reason in reasons):
        state = "nonqualifying"
    return state, (*reasons, "whole-model-" + model_state)


def construction_carriage_supported(properties, *, locus_call=False):
    return locus_call or not any(value.alternate_present is not True
        and not value.screened_low_balance_candidate for value in properties)


def hypothesis_id_for_members(context_id, member_ids, sample_name):
    """The original identity expression, shared by scalar and rich consumers."""
    return "hypothesis-" + admission._hash([context_id,member_ids,sample_name])


def construct_evaluated_hypothesis(*, context, admission_context_id, references,
        properties, assessments, phase, phase_evidence, model, enumeration,
        enumeration_evidence, sample_name, observation_ids, incompleteness_penalty,
        outside_context_partner_uncertainty, conditional=False, locus_call=False):
    """Bind one evaluated singleton/pair through the shared construction kernel.

    Upstream owners establish complete enumeration, qualification and source
    closure. This function retains their facts, IDs, uncertainty and penalties;
    it does not decide whether a nonqualifying assessment enters the universe.
    """
    if not construction_carriage_supported(properties,locus_call=locus_call):
        return None
    complete = model.completeness
    read_required = any(value.read_validation_required for value in properties)
    result_limits = ["research-hypothesis-not-solved-genotype","parental-origin-unknown"]
    if phase is HypothesisPhase.UNKNOWN:
        result_limits.append("trans-unconfirmed")
    if any(value.state == "unresolved" for value in assessments):
        result_limits.append("exploratory-qualification-unresolved")
    if complete is Completeness.INCOMPLETE:
        result_limits.append("outside-supplied-scope-partners-unknown")
    if locus_call:
        result_limits.append("repeat-locus-evidence; individual-allele-correspondence-unresolved")
    ids = tuple(ref.allele_id for ref in references)
    hypothesis = GenotypeHypothesis(
        hypothesis_id_for_members(context.context_id,ids,sample_name),context,() if locus_call else references,
        complete,phase,enumeration,context.evidence_refs,tuple(sorted(observation_ids)),
        completeness_evidence=(_record("genotype-properties",properties),),
        phase_evidence=(_record("phase",phase_evidence),) if len(ids) == 2 else (),
        enumeration_evidence=enumeration_evidence,locus_call_refs=references if locus_call else ())
    conditional = (conditional and not outside_context_partner_uncertainty
                   and complete is Completeness.INCOMPLETE and model.state == "supported")
    return ConstructedHypothesis(hypothesis,admission_context_id,assessments,model,conditional,
        incompleteness_penalty if complete is Completeness.INCOMPLETE else None,read_required,tuple(result_limits))


def evaluate_complete_context(
    *,
    context,
    relevance: ContextRelevanceAssessment,
    keys: tuple[str, ...],
    alleles: Mapping[str, AlleleEvidence],
    genotypes: SuppliedGenotypeProperties,
    policy: HypothesisPolicy,
    disease_gene_evidence: tuple[DiseaseGeneEvidence, ...],
    sample_name: str,
    observation_ids: tuple[str, ...],
    input_sha256: str,
    inventory_sha256: str,
    policy_sha256: str,
    pair_universe: _PairUniverse | None,
    source_inventory_status: Literal["complete", "assessed-domain-complete"] = "complete",
    outside_context_partner_uncertainty: bool = False,
) -> tuple[PartnerSearchAssessment, tuple[ConstructedHypothesis, ...]]:
    """Run the shared hypothesis kernel for one already closed context.

    Callers own source/admission authentication and context closure.  The kernel
    never treats a partition as a pairing universe; ``keys`` must name the whole
    complete context.  A filtered or otherwise wider source universe is carried
    as explicit partner uncertainty and cannot satisfy the conditional-singleton
    search.
    """

    if (
        type(relevance) is not ContextRelevanceAssessment
        or relevance.context != context
        or type(keys) is not tuple
        or keys != tuple(sorted(set(keys)))
        or any(key not in alleles for key in keys)
        or type(genotypes) is not SuppliedGenotypeProperties
        or source_inventory_status not in {"complete", "assessed-domain-complete"}
        or type(outside_context_partner_uncertainty) is not bool
    ):
        raise GenotypeHypothesesError("complete context inputs are invalid")
    from genosapiens.genotype_evaluation import prepare_genotype_evaluation
    if _requires_pairs(context):
        if type(pair_universe) is not _PairUniverse:
            raise GenotypeHypothesesError("complete pair universe is missing")
        pair_universe.require(keys, alleles, policy.maximum_pairs_per_context)
        total, selected = pair_universe
    else:
        if pair_universe is not None:
            raise GenotypeHypothesesError("unexpected pair universe")
        total, selected = 0, ()
    prepared = prepare_genotype_evaluation(context=context, alleles=alleles,
        genotypes=genotypes, gene_evidence=disease_gene_evidence,
        policies={'requirements': policy.requirements, 'clinical': policy.clinical, 'loci': policy.loci},
        pairs=selected)
    props, phases = prepared.properties, prepared.phases
    if set(keys) - set(props) or any(pair not in phases for pair in selected):
        raise GenotypeHypothesesError("complete context genotype coverage is invalid")
    model_contexts = {}
    def context_for(model):
        if model not in model_contexts:
            model_contexts[model] = _context(context, model)
        return model_contexts[model]

    def qualify(ctx, key, phase=None, model=None):
        allele = alleles[key]
        requirements = assess_model_requirements(
            ctx, allele, policy=policy.requirements, prepared=prepared,
        )
        result = assess_partner_qualification(
            ctx, allele, clinical_policy=policy.clinical,
            requirement_assessments=requirements,
            genotype_properties=genotypes, phase=phase,
            disease_gene_evidence=disease_gene_evidence, prepared=prepared,
        )
        if model is not None:
            state, reasons = qualify_whole_model(result.state, result.reasons, model.state)
            result = replace(
                result, state=state,
                reasons=reasons,
                binding={
                    **result.binding,
                    "whole_model_sha256": admission._hash(model),
                },
            )
        return result

    pair_context = context_for("compound_heterozygous")
    pair_models = tuple(
        assess_genotype_model(
            pair_context, tuple(alleles[value] for value in pair),
            genotypes, loci=policy.loci, prepared=prepared,
        )
        for pair in selected
    )
    pair_assessments = tuple(
        tuple(
            qualify(pair_context, member, phases[pair], model)
            for member in pair
        )
        for pair, model in zip(selected, pair_models, strict=True)
    )
    # No pair-wide rescan per member. The only partner facts needed by the
    # conditional singleton certificate are count and all-nonqualifying.
    partner_counts, qualifying_partners = Counter(), set()
    if not outside_context_partner_uncertainty:
        for left, right in pair_assessments:
            partner_counts[left.allele_id] += 1
            partner_counts[right.allele_id] += 1
            if right.state != "nonqualifying": qualifying_partners.add(left.allele_id)
            if left.state != "nonqualifying": qualifying_partners.add(right.allele_id)
    events = {key: _allele_event(alleles[key]) for key in keys}
    event_counts = Counter(events[key] for key in keys if not longtr.is_parent(alleles[key]))
    distinct_members = sum(event_counts.values())
    members = tuple(qualify(context, member) for member in keys)
    unresolved_source = tuple(
        entry.source_allele_id for entry in relevance.entries
        if entry.state != "excluded_by_evidence" and (
            entry.unresolved_entry_ids
            or entry.state == "relevant" and any(
                allele_id in props
                and props[allele_id].alternate_present is None
                and not (
                    longtr.is_parent(alleles[allele_id])
                    and longtr.metadata(alleles[allele_id]).get(
                        "source_alt_scope_complete"
                    ) is True
                )
                for allele_id in entry.derived_allele_ids
            )
        )
    )
    unresolved_relevance = tuple(
        entry.source_allele_id for entry in relevance.entries
        if entry.state == "unresolved_relevance"
    )
    truncated = len(selected) != total
    enumeration = (
        EnumerationStatus.TRUNCATED if truncated else EnumerationStatus.COMPLETE
    )
    domain_value = {
        "rules": EVALUABLE_DOMAIN_RULES,
        "policy_sha256": policy_sha256,
        "context_id": context.context_id,
        "input_sha256": input_sha256,
        "inventory_sha256": inventory_sha256,
        "source_allele_ids": tuple(
            entry.source_allele_id for entry in relevance.entries
        ),
        "evaluable_allele_ids": keys,
        "association_proofs_sha256": admission._hash(relevance),
    }
    if source_inventory_status != "complete" or outside_context_partner_uncertainty:
        domain_value.update({
            "source_inventory_status": source_inventory_status,
            "outside_context_partner_uncertainty":
                outside_context_partner_uncertainty,
        })
    domain = _record("source-evaluable-partner-domain", domain_value)
    limits = [
        "outside-supplied-scope-partners-unknown",
        "qualification-not-population-eligibility",
        "enumeration-is-source-evaluable-domain-not-exhaustive-discovery",
    ]
    if outside_context_partner_uncertainty:
        limits.append("excluded-or-unresolved-source-partners-not-absence-evidence")
    search = PartnerSearchAssessment(
        context.context_id, relevance, members, pair_assessments, pair_models,
        total, len(selected), "truncated" if truncated else "complete",
        unresolved_source, unresolved_relevance, enumeration, tuple(limits),
        source_inventory_status, keys, domain,
    )
    enumeration_evidence = (_record("partner-search", {
        "context_id": context.context_id,
        "search_sha256": admission._hash(search),
    }),)
    results = []

    def construct(ctx, ids, assessments, phase, model, conditional=False, locus_call=False):
        result = construct_evaluated_hypothesis(context=ctx,admission_context_id=context.context_id,
            references=tuple(AlleleReference(key,alleles[key].native_identity_id) for key in ids),
            properties=tuple(props[key] for key in ids),assessments=assessments,phase=phase,
            phase_evidence=phases[ids] if len(ids) == 2 else None,model=model,enumeration=enumeration,
            enumeration_evidence=enumeration_evidence,sample_name=sample_name,observation_ids=observation_ids,
            incompleteness_penalty=policy.incompleteness_penalty,
            outside_context_partner_uncertainty=outside_context_partner_uncertainty,
            conditional=conditional,locus_call=locus_call)
        if result is not None:
            results.append(result)

    for ids, assessments, model in zip(
        selected, pair_assessments, pair_models, strict=True,
    ):
        if not any(value.state == "nonqualifying" for value in assessments):
            construct(
                pair_context, ids, assessments, phases[ids].phase, model,
            )
    for member in members:
        allele_id = member.allele_id
        ctx = context_for(_model(context, alleles[allele_id], props[allele_id]))
        model = assess_genotype_model(
            ctx, (alleles[allele_id],), genotypes, loci=policy.loci, prepared=prepared,
        )
        own = qualify(ctx, allele_id, model=model)
        no_partner = (
            enumeration is EnumerationStatus.COMPLETE
            and not unresolved_relevance
            and not unresolved_source
            and not outside_context_partner_uncertainty
            and own.state == "qualifying"
            and allele_id not in qualifying_partners
            and partner_counts[allele_id] == distinct_members - event_counts[events[allele_id]]
        )
        if own.state != "nonqualifying":
            native = admission._native(alleles[allele_id])
            interpretation = (
                native.source.source_fields or {}
            ).get("native_interpretation") or {}
            locus_call = (
                isinstance(native.identity, RepeatAllele)
                and interpretation.get("repeat_representation") == "locus_call"
            )
            construct(
                ctx, (allele_id,), (own,), HypothesisPhase.NOT_APPLICABLE,
                model,
                conditional=(
                    ctx.model == "incomplete"
                    and no_partner
                    and not props[allele_id].read_validation_required
                ),
                locus_call=locus_call,
            )
    return search, tuple(results)


def replay_genotype_hypotheses(replay: HypothesisReplay) -> HypothesisGeneration:
    """Use the one exact public replay path for diagnostics, fusion, and consumers."""
    if type(replay) is not HypothesisReplay:
        raise GenotypeHypothesesError("typed hypothesis replay required")
    return generate_genotype_hypotheses(**{
        item.name: getattr(replay, item.name) for item in fields(HypothesisReplay)
    })


def generate_genotype_hypotheses(
    *, inputs: SuppliedInputs, sources: admission.AdmissionSources, view: ScoringView,
    profiles: DiseaseProfile | None, genes: tuple[DiseaseGeneEvidence, ...] | None,
    ic_store: HpoICStore | None, admission_policy_sha256: str, sample_name: str,
    policy: HypothesisPolicy, expected_policy_sha256: str,
    expected_result: HypothesisGeneration | None = None,
    population_inputs: HypothesisPopulationInputs | None = None,
    representation_eligibility: tuple[RepresentationEligibilityAssessment, ...] = (),
) -> HypothesisGeneration:
    """Replay full source accounting and enumerate without population pruning.

    Persisted results can only be consumed by full replay with expected_result.
    Qualification uncertainty is not construction failure or resource truncation.
    No missing native record, unmatched rescue, or audit-only identity is dropped
    to obtain a conditional singleton. Mandatory technical checks are not cleared.
    """
    try:
        pin = hypothesis_policy_bytes(policy, admission_policy_sha256=admission_policy_sha256)
        if hashlib.sha256(pin).hexdigest() != expected_policy_sha256:
            raise ValueError
        candidates = admission.generate_candidate_admission(inputs=inputs, sources=sources, view=view,
            profiles=profiles, genes=genes, ic_store=ic_store, expected_policy_sha256=admission_policy_sha256)
        alleles = candidates.inventory.alleles
        population_result = None
        if population_inputs is not None:
            if type(population_inputs) is not HypothesisPopulationInputs or any(a.population for a in alleles):
                raise ValueError
            selected_population_policy = frequency.prepare_population_policy(
                population_inputs.policy_bytes, expected_sha256=population_inputs.policy_sha256)
            frequency.validate_population_curation(selected_population_policy.policy, population_inputs.public_curation_bytes)
            population_requests = []
            for allele in alleles:
                normalized = tuple(i.identity_id for i in allele.identities if i.stage is IdentityStage.NORMALIZED)
                if len(normalized) > 1:
                    raise ValueError
                population_requests.append(population.PopulationRequest(allele, normalized[0] if normalized else None))
            population_result = population.run_population_annotation(
                tuple(population_requests),
                joint=population_inputs.joint, all_sites=population_inputs.all_sites,
                reference=population_inputs.reference, runtime_image=population_inputs.runtime_image,
                limits=population_inputs.limits)
            alleles = population_result.alleles
        clinical = None
        if sources.clinvar is not None or sources.open_targets is not None:
            if sources.reference is None:
                raise ValueError
            clinical = annotate_clinical(alleles, clinvar=sources.clinvar,
                open_targets=sources.open_targets, reference=sources.reference)
            alleles = clinical.alleles
        elif any(a.clinical for a in alleles):
            raise ValueError
        by_id = {a.allele_id: a for a in alleles}
        representation_eligibility = tuple(sorted(
            _typed(representation_eligibility, tuple[RepresentationEligibilityAssessment, ...]),
            key=lambda item: item.allele_id,
        ))
        if (len({item.allele_id for item in representation_eligibility}) != len(representation_eligibility)
                or any(item.allele_id not in by_id
                       or item.native_identity_id != by_id[item.allele_id].native_identity_id
                       for item in representation_eligibility)):
            raise ValueError
        genotypes = assess_supplied_genotypes(inputs, sample_name=sample_name)
        props = {p.allele_id: p for p in genotypes.allele_properties}
        relevance, keys_by_context, pair_keys, pair_universes = {}, {}, set(), {}
        for seed in candidates.contexts:
            key = seed.context.context_id
            relevance[key] = assess_context_relevance(candidates, seed.context)
            # Include actual source target associations, not just variant-first
            # seeds. Different sample projections remain in source accounting.
            keys = sorted({r.subject_id for e in relevance[key].entries for r in e.routes
                           if r.matches_context_target and r.subject_id in props})
            keys_by_context[key] = keys
            if _requires_pairs(seed.context):
                pair_universes[key] = _pairs(keys, by_id, policy.maximum_pairs_per_context)
                pair_keys.update(pair_universes[key].pairs)
            else:
                pair_universes[key] = None
        genotypes = assess_supplied_genotypes(inputs, sample_name=sample_name, pairs=tuple(sorted(pair_keys)))
        searches, results = [], []
        for seed in candidates.contexts:
            key = seed.context.context_id
            search, hypotheses = evaluate_complete_context(
                context=seed.context,
                relevance=relevance[key],
                keys=tuple(keys_by_context[key]),
                alleles=by_id,
                genotypes=genotypes,
                policy=policy,
                disease_gene_evidence=candidates.disease_gene_evidence,
                sample_name=sample_name,
                observation_ids=tuple(
                    observation.observation_id
                    for observation in candidates.view.observations
                ),
                input_sha256=inputs.input_sha256,
                inventory_sha256=candidates.binding["inventory_sha256"],
                policy_sha256=expected_policy_sha256,
                pair_universe=pair_universes[key],
            )
            searches.append(search)
            results.extend(hypotheses)
        results = tuple(sorted(results, key=lambda r: r.hypothesis.hypothesis_id))
        population_assessments = () if population_inputs is None else tuple(
            frequency.evaluate_population_policy(
                item.hypothesis, tuple(by_id[key] for key in item.hypothesis.evidence_subject_ids),
                prepared_policy=selected_population_policy)
            for item in results)
        binding = _record("generation-binding", {"schema": SCHEMA, "policy_sha256": expected_policy_sha256,
            "admission_sha256": hashlib.sha256(candidates.as_bytes()).hexdigest(),
            "genotype_sha256": admission._hash(genotypes), "clinical_sha256": admission._hash(clinical),
            "population_annotation_sha256": admission._hash(population_result),
            "population_policy_sha256": None if population_inputs is None else population_inputs.policy_sha256,
            "population_assessments_sha256": admission._hash(population_assessments),
            "representation_eligibility_sha256": admission._hash(representation_eligibility),
            "sample_name": sample_name, "input_sha256": inputs.input_sha256,
            "source_set": _plain(candidates.inventory.source_set),
            "scope": "supplied-call-development; no-competitive-or-controlled-use-authorization"})
        result = HypothesisGeneration(candidates, clinical, genotypes, tuple(searches),
            results, binding, population_result, population_assessments, representation_eligibility)
        if expected_result is not None and result.as_bytes() != _typed(expected_result, HypothesisGeneration).as_bytes():
            raise ValueError
        return result
    except Exception:
        raise GenotypeHypothesesError("genotype hypothesis replay failed") from None
