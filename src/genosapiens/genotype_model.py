"""Source-bound whole-model checks, not source authentication or qualification.

The consumer replays the context, supplied genotype properties and independently
pinned public locus projection before calling. This kernel performs no I/O and
does not enumerate partners, infer sex, select clinical evidence or waive checks.
"""

from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass
import hashlib
from typing import Literal

from genosapiens.admission_inventory import supplied_structural_envelope
from genosapiens.allele_identity import (
    GenomicInterval, MitochondrialAllele, MosaicAllele, ReferenceRole, RepeatAllele,
    SequenceAllele, StructuralAllele, _reference,
)
from genosapiens.disease_gene_evidence import GeneRequirements, _G2P_ALLELES
from genosapiens.evidence_model import (
    AlleleEvidence, AlleleReference, Completeness, EvidenceReference, HypothesisContext,
    HypothesisPhase, IdentityMapping, SourcedRecord, _Controlled, _ids, _sha, _typed,
)
from genosapiens.evidence_serialization import _canonical
from genosapiens.genotype_properties import (
    POLICY_SHA256 as GENOTYPE_POLICY_SHA256, AlleleGenotypeProperties, PairPhaseProperties,
    SourceGenotypeProperties, SuppliedGenotypeProperties, _allele_event, selected_source_calls, source_sample_name,
)
from genosapiens import longtr_admission as longtr
from genosapiens.vep_native_annotation import native_interval_precision


SCHEMA = "genosapiens.genotype-model"
LOCUS_SCHEMA = "genosapiens.model-locus-region"
State = Literal["supported", "contradicted", "unresolved"]
LocusKind = Literal["autosome", "x_nonpar", "y_nonpar", "par", "mitochondrial"]
MONOALLELIC_HOMOZYGOTE_REASON = "monoallelic-autosomal-homozygote-applicability-unresolved"

# These are the source adapter's actual tokens, not spelling/contig heuristics.
_LITERAL_RULES = {
    "monoallelic_autosomal": ("dominant", "autosome", (2,)),
    "biallelic_autosomal": ("recessive", "autosome", (2,)),
    "monoallelic_PAR": ("dominant", "par", (2,)),
    "biallelic_PAR": ("recessive", "par", (2,)),
    "monoallelic_X_hemizygous": ("x_linked", "x_nonpar", (1,)),
    "monoallelic_X_heterozygous": ("x_linked", "x_nonpar", (2,)),
    "monoallelic_Y_hemizygous": ("y_linked", "y_nonpar", (1,)),
    "mitochondrial": ("mitochondrial", "mitochondrial", ()),
}
_INHERITANCE_RULES = {
    "HP:0000006": _LITERAL_RULES["monoallelic_autosomal"],
    "HP:0000007": _LITERAL_RULES["biallelic_autosomal"],
    "HP:0001417": ("x_linked", "x_nonpar", (1, 2)),
    "HP:0001419": ("x_recessive", "x_nonpar", (1, 2)),
    "HP:0001423": ("x_linked", "x_nonpar", (1, 2)),
    "HP:0001450": _LITERAL_RULES["monoallelic_Y_hemizygous"],
    "HP:0001427": _LITERAL_RULES["mitochondrial"],
    **{_G2P_ALLELES[token][0]: rule for token, rule in _LITERAL_RULES.items()
       if _G2P_ALLELES[token][0].startswith("g2p:")},
}
_MODELS = {
    "dominant": ("dominant",),
    "recessive": ("homozygous_recessive", "compound_heterozygous", "biallelic", "incomplete"),
    "x_linked": ("x_linked",),
    "x_recessive": ("x_linked", "homozygous_recessive", "compound_heterozygous", "incomplete"),
    "y_linked": ("y_linked",),
    "mitochondrial": ("mitochondrial",),
}
POLICY_BYTES = (_canonical({
    "schema": SCHEMA, "literal_rules": _LITERAL_RULES, "inheritance_rules": _INHERITANCE_RULES,
    "models": _MODELS, "source_count_vocabulary": _G2P_ALLELES,
    "counts": "source-required-only; exact-supported-copies-and-model-ploidy; no-duplicate-events",
    "monoallelic_autosomal_homozygote": "retain-source-copies-and-candidate; applicability-unresolved-without-explicit-source-model; neither-minimum-count-nor-exact-count-contradiction",
    "loci": "public-source-projection; exact-reference; complete-closed-span; no-aliases-or-sex",
    "native": "all-source-bounded-geometry-or-supplied-endpoint-CI-envelope; locus-envelope-not-affected-bases; unavailable-locus-does-not-erase-independent-autosomal-copies",
    "phase": "recessive-pairs-require-trans; unknown-permitted; conflict-unresolved; copies-independent",
    "repeat": "explicit-source-locus-context-can-support-locus-model-with-unknown-genotype",
    "longtr_projections": longtr.POLICY_SHA256,
    "limits": "unprojected-repeat-copies-unknown; screened-candidates-unresolved; mitochondrial-GT-not-organism-copies",
}) + "\n").encode("ascii")
POLICY_SHA256 = hashlib.sha256(POLICY_BYTES).hexdigest()


class GenotypeModelError(ValueError):
    """Static diagnostics; malformed source bindings are not biological unknowns."""


def _hash(value):
    return hashlib.sha256((_canonical(value) + "\n").encode("ascii")).hexdigest()


def requires_monoallelic_homozygote_evidence(requirement: GeneRequirements, genotype: AlleleGenotypeProperties) -> bool:
    """Identify the approved source-specific uncertainty, not a contradiction.

    The current source vocabulary's monoallelic_autosomal label supplies no
    homozygous-model applicability assertion. Ordinary clinical pathogenicity,
    an unrelated biallelic context, or an extra caller flag cannot supply it.
    Copies come from replayed source ALT indices, never literal GT spelling.
    """
    return (
        requirement.allelic_requirement_kind == "literal_allelic_model"
        and requirement.allelic_requirement == "monoallelic_autosomal"
        and not genotype.conflicting_source_genotypes
        and not genotype.screened_low_balance_candidate
        and genotype.alternate_present is True
        and (genotype.minimum_alt_copies, genotype.maximum_alt_copies, genotype.reported_ploidies) == (2, 2, (2,))
    )


@dataclass(frozen=True, slots=True, repr=False)
class ModelLocusRegion(_Controlled):
    """Public-derived 1-based closed geometry, independently pinned by the consumer.

    source_fields['model_locus_region'] must contain this exact projection, with
    schema=LOCUS_SCHEMA and coordinate_system='1-based-closed'. Original public
    fields/derivation evidence may also be retained. Matching hashes/constructors
    do not authenticate a projection or authorize a caller to invent its kind.
    """

    reference_role: ReferenceRole
    reference_identity_id: str
    reference_manifest_sha256: str
    contig: str
    start: int
    end: int
    kind: LocusKind
    source: SourcedRecord

    def __post_init__(self):
        _Controlled.__post_init__(self)
        _reference(self.reference_role, self.reference_identity_id)
        _sha(self.reference_manifest_sha256)
        GenomicInterval(self.contig, self.start, self.end)
        expected = {"schema": LOCUS_SCHEMA, "coordinate_system": "1-based-closed",
            "reference_role": self.reference_role.value, "reference_identity_id": self.reference_identity_id,
            "reference_manifest_sha256": self.reference_manifest_sha256, "contig": self.contig,
            "start": self.start, "end": self.end, "kind": self.kind}
        if (self.source.public_manifest_sha256 is None or self.source.source_fields is None
                or _canonical(self.source.source_fields.get("model_locus_region")) != _canonical(expected)):
            raise GenotypeModelError("public source-bound locus projection required")


@dataclass(frozen=True, slots=True, repr=False)
class ModelLocusFacts(_Controlled):
    kind: LocusKind | None
    geometry: GenomicInterval | None
    regions: tuple[ModelLocusRegion, ...]
    sources: tuple[SourcedRecord, ...]
    reasons: tuple[str, ...]


@dataclass(frozen=True, slots=True, repr=False)
class ModelMemberFacts(_Controlled):
    allele: AlleleReference
    native: IdentityMapping
    genotype: AlleleGenotypeProperties
    source_genotypes: tuple[SourceGenotypeProperties, ...]
    locus: ModelLocusFacts


@dataclass(frozen=True, slots=True, repr=False)
class ModelRequirementFacts(_Controlled):
    requirement: GeneRequirements
    source_refs: tuple[EvidenceReference, ...]
    family: str | None
    required_locus_kind: LocusKind | None
    state: State
    completeness: Completeness
    reasons: tuple[str, ...]


@dataclass(frozen=True, slots=True, repr=False)
class GenotypeModelAssessment(_Controlled):
    context_id: str
    model: str
    state: State
    completeness: Completeness
    reasons: tuple[str, ...]
    members: tuple[ModelMemberFacts, ...]
    requirement_facts: tuple[ModelRequirementFacts, ...]
    phase: HypothesisPhase
    phase_evidence: PairPhaseProperties | None
    read_validation_required: bool
    screened_low_balance_candidate: bool
    binding: SourcedRecord
    authentication: Literal["consumer_must_replay_sources"] = "consumer_must_replay_sources"


def _locus(allele, native, loci):
    identity = native.identity
    if isinstance(identity, (MitochondrialAllele, MosaicAllele)):
        identity = identity.sequence
    sources, geometry, geometry_limits = (native.source,), None, ()
    if isinstance(identity, SequenceAllele):
        geometry = GenomicInterval(identity.contig, identity.position,
                                   identity.position + len(identity.reference) - 1)
    elif isinstance(identity, StructuralAllele):
        precision = native_interval_precision(allele)
        sources += (precision,)
        contributors = precision.source_fields["contributors"]
        if identity.interval is not None and contributors and all(
            c["recipe"] is not None and c["reason"] in (None, "source-endpoints-imprecise") for c in contributors
        ):
            recipes = [c["recipe"] for c in contributors]
            smin, smax = max(r["start"]["min"] for r in recipes), min(r["start"]["max"] for r in recipes)
            emin, emax = max(r["end"]["min"] for r in recipes), min(r["end"]["max"] for r in recipes)
            if smin <= identity.interval.start <= smax and emin <= identity.interval.end <= emax:
                geometry = GenomicInterval(identity.interval.contig, smin, emax)
        if geometry is None:
            geometry = supplied_structural_envelope(allele)
            if geometry is not None:
                sources += tuple(c.source for c in allele.calls or ())
                geometry_limits = ("supplied-endpoint-CI-locus-envelope; not-affected-bases-or-molecular-effect",)
    elif isinstance(identity, RepeatAllele):
        geometry = identity.interval
    if geometry is None:
        return ModelLocusFacts(None, None, (), sources, ("source-native-geometry-unavailable",))
    regions = tuple(r for r in loci if (
        r.reference_role, r.reference_identity_id, r.reference_manifest_sha256, r.contig
    ) == (identity.reference_role, identity.reference_identity_id, native.reference_manifest_sha256, geometry.contig)
        and r.start <= geometry.end and r.end >= geometry.start)
    cursor = geometry.start
    for region in sorted(regions, key=lambda r: (r.start, r.end)):
        if region.start > cursor:
            break
        cursor = max(cursor, region.end + 1)
    kinds = {r.kind for r in regions}
    known = cursor > geometry.end and len(kinds) == 1
    reason = ("complete-source-bound-locus-span" if known else "locus-kind-conflict-or-boundary-crossing"
              if len(kinds) > 1 else "complete-reference-bound-locus-coverage-unavailable")
    return ModelLocusFacts(next(iter(kinds)) if known else None, geometry, regions, sources, (reason, *geometry_limits))


@dataclass(frozen=True, slots=True, repr=False)
class ModelDecisionGenotype:
    """All genotype inputs read by the model state/completeness decision."""
    alternate_present: bool | None
    minimum_alt_copies: int | None
    maximum_alt_copies: int | None
    reported_ploidies: tuple[int, ...]
    conflicting_source_genotypes: bool
    screened_low_balance_candidate: bool
    phase_conflict_bases: bool


@dataclass(frozen=True, slots=True, repr=False)
class ModelDecisionMember:
    genotype: ModelDecisionGenotype
    locus_kind: LocusKind | None
    locus_reasons: tuple[str, ...]
    repeat_locus_id: str | None
    repeat_locus_call: bool | None
    unprojected_repeat: bool
    nonorganism_copy_count: bool


@dataclass(frozen=True, slots=True, repr=False)
class ModelRequirementDecision:
    family: str | None
    required_locus_kind: LocusKind | None
    state: State
    completeness: Completeness
    reasons: tuple[str, ...]


def model_decision_member(member, *, projected):
    """Project decision dependencies, retaining the full member separately.

    Source identities, geometry and phase bases stay in ModelMemberFacts. Their
    interpreted locus kind/reasons, native class and genotype facts are the only
    dependencies of the state kernel below; none is replaced by an inferred fact.
    """
    prop, native = member.genotype, member.native
    repeat = isinstance(native.identity, RepeatAllele)
    # Only the repeat-locus branch reads this metadata. An unrelated source
    # field on another native class must not become a new decision dependency.
    interpretation = ((native.source.source_fields or {}).get("native_interpretation") or {}
                      ) if repeat else {}
    return ModelDecisionMember(
        ModelDecisionGenotype(prop.alternate_present, prop.minimum_alt_copies,
            prop.maximum_alt_copies, prop.reported_ploidies, prop.conflicting_source_genotypes,
            prop.screened_low_balance_candidate, bool(prop.phase_conflict_bases)),
        member.locus.kind, member.locus.reasons,
        native.identity.locus_id if repeat else None,
        (repeat and interpretation.get("repeat_representation") == "locus_call")
            if isinstance(interpretation, Mapping) else None,
        repeat and member.allele.allele_id not in projected,
        isinstance(native.identity, (MitochondrialAllele, MosaicAllele)),
    )


def evaluate_model_requirement(*, model, locus_id, has_locus_binding, requirement, members, phase):
    """The single model requirement decision, over explicit immutable inputs.

    Bounded evaluation may group identical tuples of these inputs. Every member
    of such a group has exactly this already evaluated result; source records
    and result bindings are attached independently, never borrowed from a peer.
    """
    token = requirement.allelic_requirement
    rule = _INHERITANCE_RULES.get(requirement.inheritance)
    if requirement.allelic_requirement_kind == "literal_allelic_model":
        rule = _LITERAL_RULES.get(token)
        expected = _G2P_ALLELES.get(token)
        if expected is None or requirement.inheritance != expected[0] or (
            requirement.required_allele_count is not None and requirement.required_allele_count != expected[1]
        ):
            rule = None
    elif token is not None and requirement.allelic_requirement_kind != "mode_of_inheritance":
        rule = None
    if not requirement.record_ids or rule is None or requirement.required_allele_count not in (None, 1, 2):
        return ModelRequirementDecision(None, None, "unresolved", Completeness.UNKNOWN,
                                     ("source-model-interpretation-unavailable",))
    family, locus_kind, ploidies = rule
    if model == "repeat" and len(members) == 1:
        member, = members
        if (member.repeat_locus_id is not None and locus_id == member.repeat_locus_id
                and has_locus_binding and member.repeat_locus_call is None):
            # Match the expanded kernel's refusal exactly at its original read
            # boundary; irrelevant malformed metadata is not a new condition.
            raise AttributeError("invalid native repeat interpretation")
        if (member.repeat_locus_id is not None and locus_id == member.repeat_locus_id
                and has_locus_binding and member.repeat_locus_call):
            state = ("unresolved" if member.locus_kind is None else
                     "supported" if member.locus_kind == locus_kind else "contradicted")
            return ModelRequirementDecision(family, locus_kind, state, Completeness.UNKNOWN,
                (*member.locus_reasons, "source-native-repeat-locus-model", "individual-allele-genotype-unresolved",
                 "locus-support-not-carriage-count-phase-or-clinical-qualification"))
    if model not in _MODELS[family]:
        # Native class labels and source-requirement-context-only are not models.
        state = "contradicted" if model in {m for models in _MODELS.values() for m in models} else "unresolved"
        return ModelRequirementDecision(family, locus_kind, state, Completeness.UNKNOWN,
                                     ("requested-model-not-supported-by-source-requirement",))

    reasons, negative, unavailable, conflict = [], False, False, False
    pair = len(members) == 2
    cardinality_ok = pair if model == "compound_heterozygous" else (
        not pair or model == "biallelic")
    if not cardinality_ok:
        negative = True
        reasons.append("distinct-member-cardinality-contradicts-model")
    expected_ploidies = (2,) if pair or model in ("homozygous_recessive", "incomplete") else ploidies
    count_usable = cardinality_ok and family != "mitochondrial"
    for member in members:
        p = member.genotype
        reasons.extend(member.locus_reasons)
        if member.locus_kind is None:
            unavailable = True
            # Autosomal copy accounting can be independently known from GT, but
            # X/Y/PAR applicability defines which genotype model is being counted.
            if locus_kind != "autosome":
                count_usable = False
        elif member.locus_kind != locus_kind:
            negative, count_usable = True, False
            reasons.append("bound-locus-contradicts-source-model")
        if p.conflicting_source_genotypes:
            conflict, count_usable = True, False
            reasons.append("material-source-genotype-conflict")
        elif p.screened_low_balance_candidate:
            unavailable, count_usable = True, False
            reasons.append("screen-candidate-not-confirmed-carriage; source-GT-unchanged")
        elif member.unprojected_repeat:
            unavailable, count_usable = True, False
            reasons.append("source-ALT-not-bound-to-native-repeat-sizes")
        elif p.alternate_present is False:
            negative = True
            reasons.append("source-GT-excludes-member-carriage")
        elif p.alternate_present is not True:
            unavailable = True
            reasons.append("member-carriage-unresolved")
        if member.nonorganism_copy_count:
            count_usable = False
        if expected_ploidies and not p.conflicting_source_genotypes:
            if p.reported_ploidies and not set(p.reported_ploidies) <= set(expected_ploidies):
                negative, count_usable = True, False
                reasons.append("source-ploidy-contradicts-model")
            elif len(p.reported_ploidies) != 1:
                unavailable, count_usable = True, False
                reasons.append("model-ploidy-unresolved")
        exact_copies = 2 if model == "homozygous_recessive" else 1 if pair or (
            requirement.inheritance in ("g2p:monoallelic_X_hemizygous", "g2p:monoallelic_X_heterozygous",
                                        "g2p:monoallelic_Y_hemizygous")
        ) else None
        if not p.conflicting_source_genotypes and not p.screened_low_balance_candidate:
            if exact_copies is not None:
                if (p.maximum_alt_copies is not None and p.maximum_alt_copies < exact_copies
                        or p.minimum_alt_copies is not None and p.minimum_alt_copies > exact_copies):
                    negative = True
                    reasons.append("source-copies-contradict-model-shape")
                elif (p.minimum_alt_copies, p.maximum_alt_copies) != (exact_copies, exact_copies):
                    unavailable = True
                    reasons.append("model-copy-shape-unresolved")
            if model == "incomplete" and p.minimum_alt_copies is not None and p.minimum_alt_copies >= 2:
                negative = True
                reasons.append("supported-two-copy-singleton-is-not-incomplete-shape")
        if p.minimum_alt_copies is None or p.minimum_alt_copies != p.maximum_alt_copies:
            count_usable = False
            reasons.append("full-source-copy-count-unresolved")
        if pair and (p.minimum_alt_copies, p.maximum_alt_copies) != (1, 1):
            count_usable = False
        if requires_monoallelic_homozygote_evidence(requirement, p):
            unavailable, count_usable = True, False
            reasons.append(MONOALLELIC_HOMOZYGOTE_REASON)
    completeness = Completeness.UNKNOWN
    required = requirement.required_allele_count
    if required is None:
        reasons.append("required-model-count-unknown; no-inheritance-count-inference")
    elif count_usable:
        copies = sum(m.genotype.minimum_alt_copies for m in members)
        if copies == required:
            completeness = Completeness.COMPLETE
        elif copies < required:
            completeness = Completeness.INCOMPLETE
        reasons.append("independent-source-copy-accounting")
    if family == "mitochondrial":
        reasons.append("mitochondrial-GT-not-organism-copy-count")
    if pair and cardinality_ok:
        if phase is HypothesisPhase.CONFLICTING or any(m.genotype.phase_conflict_bases for m in members):
            conflict = True
            reasons.append("material-required-phase-conflict")
        elif phase is HypothesisPhase.IN_CIS:
            negative = True
            reasons.append("established-cis-contradicts-pair-model")
        elif phase is HypothesisPhase.UNKNOWN:
            reasons.append("ordinary-unknown-phase-retained; not-trans-or-veto")
    state = "unresolved" if conflict else "contradicted" if negative else "unresolved" if unavailable else "supported"
    return ModelRequirementDecision(family, locus_kind, state, completeness,
                                 tuple(sorted(set(reasons) | {"source-model-shape-evaluated"})))


def _requirement(context, requirement, members, phase, projected):
    refs = tuple(r for r in context.evidence_refs if r.record_id in requirement.record_ids)
    decision = evaluate_model_requirement(model=context.model, locus_id=context.locus_id,
        has_locus_binding=bool(context.locus_bindings), requirement=requirement,
        members=tuple(model_decision_member(member, projected=projected) for member in members), phase=phase)
    return ModelRequirementFacts(requirement, refs, decision.family, decision.required_locus_kind,
        decision.state, decision.completeness, decision.reasons)


def summarize_model_requirements(facts, requirements):
    """Combine evaluated source alternatives without inferring missing counts."""
    states = {r.state for r in facts}
    counts = {r.required_allele_count for r in requirements}
    complete = {r.completeness for r in facts}
    state = next(iter(states)) if len(states) == 1 else "unresolved"
    completeness = next(iter(complete)) if len(complete) == 1 and len(counts) == 1 else Completeness.UNKNOWN
    reasons = {reason for r in facts for reason in r.reasons}
    if not facts:
        reasons.add("source-model-requirements-unavailable")
    if len(states) > 1 or len(counts - {None}) > 1 or len({(f.family, f.required_locus_kind) for f in facts}) > 1:
        state = "unresolved"
        reasons.add("source-model-alternatives-not-combined")
    return state, completeness, tuple(sorted(reasons))


def assess_genotype_model(
    context: HypothesisContext, alleles: tuple[AlleleEvidence, ...], genotypes: SuppliedGenotypeProperties,
    *, loci: tuple[ModelLocusRegion, ...], prepared=None,
) -> GenotypeModelAssessment:
    """Assess one singleton or distinct unordered pair before qualification.

    State is model construction support, NOT clinical qualification. Completeness
    is independent copy accounting and may remain COMPLETE with contradicted cis
    or unresolved native locus geometry. Always gate positive/complete/no-partner
    claims on state as well. Unknown phase/count do not fabricate missing facts.
    Source authentication and executable/public-policy pin verification are upstream.
    """
    try:
        from genosapiens.genotype_evaluation import prepare_genotype_evaluation
        if prepared is None:
            _typed(loci, tuple[ModelLocusRegion, ...])
            prepared = prepare_genotype_evaluation(context=context, alleles=alleles,
                genotypes=genotypes, gene_evidence=(), policies={'loci': loci}, pairs=())
        prepared.require(context=context, alleles=alleles, genotypes=genotypes, policies={'loci': loci})
        alleles = tuple(sorted(alleles, key=lambda a: a.allele_id))
        if len(alleles) not in (1, 2):
            raise ValueError
        _ids(tuple(a.allele_id for a in alleles))
        loci = prepared.loci
        props, source_props = prepared.properties, prepared.source_properties
        members = []
        for allele in alleles:
            if allele.allele_id in prepared.model_members:
                members.append(prepared.model_members[allele.allele_id])
                continue
            native = next(i for i in allele.identities if i.identity_id == allele.native_identity_id)
            prop = props[allele.allele_id]
            _ids(prop.source_allele_ids)
            calls = prepared.calls_for(allele)
            if (prop.native_identity_id != native.identity_id
                    or native.reference_manifest_sha256 != genotypes.reference_manifest_sha256 or not calls
                    or {c.source.source_fields["source_allele_id"] for c in calls} != set(prop.source_allele_ids)):
                raise ValueError
            sources = tuple(source_props[key] for key in prop.source_allele_ids)
            # A replay-established native derivative may coexist with unresolved
            # measurements. Membership, not source-wide success, binds the member.
            if any(s.sample.sample_name != source_sample_name(genotypes, s.input_sha256 or genotypes.input_sha256)
                   or allele.allele_id not in s.derived_allele_ids
                   or s.source_status == "non_variant" for s in sources) or any(
                c.sample != source_props[c.source.source_fields["source_allele_id"]].sample for c in calls
            ):
                raise ValueError
            event_sources = tuple(source_props[key] for key in (prop.event_source_allele_ids or prop.source_allele_ids))
            if not set(prop.source_allele_ids) <= {s.source_allele_id for s in event_sources}:
                raise ValueError
            members.append(ModelMemberFacts(AlleleReference(allele.allele_id, native.identity_id), native,
                                            prop, event_sources, _locus(allele, native, loci)))
            prepared.model_members[allele.allele_id] = members[-1]
        members = tuple(members)
        if len({_allele_event(a) for a in alleles}) != len(members) or len(members) == 2 and any(longtr.is_parent(a) for a in alleles):
            raise ValueError
        ids = tuple(a.allele_id for a in alleles)
        phase_evidence = prepared.phases.get(ids)
        if phase_evidence is not None and phase_evidence.phase is HypothesisPhase.NOT_APPLICABLE:
            raise ValueError
        phase = (phase_evidence.phase if phase_evidence else HypothesisPhase.UNKNOWN
                 ) if len(members) == 2 else HypothesisPhase.NOT_APPLICABLE
        projected = {a.allele_id for a in alleles if longtr.is_projection(a)}
        facts = tuple(_requirement(context, r, members, phase, projected) for r in context.requirements)
        state, completeness, reasons = summarize_model_requirements(facts, context.requirements)
        binding_fields = {"schema": SCHEMA, "policy_sha256": POLICY_SHA256,
            "context_sha256": prepared.context_key(context), "alleles_sha256": prepared.alleles_hash(alleles),
            **prepared.genotype_binding_fields(), "locus_policy_sha256": prepared.policy_hashes["loci"],
            "input_sha256": genotypes.input_sha256, "inventory_sha256": genotypes.inventory_sha256,
            "reference_manifest_sha256": genotypes.reference_manifest_sha256, "sample_name": genotypes.sample_name}
        digest = _hash(binding_fields)
        binding = SourcedRecord(SCHEMA, "current", "model-binding-" + digest, digest, None, binding_fields)
        return GenotypeModelAssessment(context.context_id, context.model, state, completeness, tuple(sorted(reasons)),
            members, facts, phase, phase_evidence, any(m.genotype.read_validation_required for m in members),
            any(m.genotype.screened_low_balance_candidate for m in members), binding)
    except (ValueError, TypeError, KeyError, AttributeError, RuntimeError, OverflowError, RecursionError):
        raise GenotypeModelError("invalid genotype model inputs or source bindings") from None
