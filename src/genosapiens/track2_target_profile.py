"""Target views over existing mechanism, source and candidate decisions."""

from dataclasses import asdict, dataclass

from genosapiens.track2_target_hypothesis import TargetHypothesis


@dataclass(frozen=True)
class KnownDrugHit:
    drug_id: str
    name: str
    reported_action: str
    action_direction: str
    maximum_clinical_phase: str
    candidate_id: str | None
    evidence_ids: tuple[str, ...]

    def __post_init__(self):
        if (not self.drug_id or not self.name or not self.reported_action or not self.evidence_ids
                or self.action_direction not in {"increase", "decrease", "mixed", "unknown"}):
            raise ValueError("known-drug identity/action provenance is incomplete")


@dataclass(frozen=True)
class KnownDrugLookup:
    source: str
    query_id: str
    state: str
    maximum_results: int
    hits: tuple[KnownDrugHit, ...]
    evidence_ids: tuple[str, ...]
    limitation: str

    def __post_init__(self):
        if (self.state not in {"success", "no_result", "unavailable", "invalid", "truncated", "screened_exclusion"}
                or not self.source or not self.query_id or not self.evidence_ids or not self.limitation
                or type(self.maximum_results) is not int or self.maximum_results < 1
                or len(self.hits) > self.maximum_results
                or (self.state in {"no_result", "unavailable", "invalid", "screened_exclusion"} and self.hits)
                or (self.state == "success" and not self.hits)):
            raise ValueError("known-drug lookup status or bounds are inconsistent")


@dataclass(frozen=True)
class ExpressionObservation:
    tissue_or_cell_type: str
    assay: str
    population: str
    result: str
    evidence_ids: tuple[str, ...]

    def __post_init__(self):
        if not all((self.tissue_or_cell_type, self.assay, self.population, self.result, self.evidence_ids)):
            raise ValueError("baseline expression requires tissue, assay, population and source")


@dataclass(frozen=True)
class TargetProfileInput:
    target_id: str
    residual_functional_product: str
    residual_product_context: str
    residual_product_evidence_ids: tuple[str, ...]
    separately_supported_restoration: bool
    restoration_evidence_ids: tuple[str, ...]
    known_drug_lookups: tuple[KnownDrugLookup, ...]
    expression_state: str
    expression_evidence_ids: tuple[str, ...]
    expression_observations: tuple[ExpressionObservation, ...]
    human_genetic_context: str
    human_genetic_evidence_ids: tuple[str, ...]

    def __post_init__(self):
        if (not self.target_id or self.residual_functional_product not in {"present", "absent", "unknown"}
                or not self.residual_product_context or not self.residual_product_evidence_ids
                or not self.known_drug_lookups
                or self.expression_state not in {"success", "no_result", "unavailable", "invalid", "truncated"}
                or not self.expression_evidence_ids or not self.human_genetic_context
                or (self.separately_supported_restoration and not self.restoration_evidence_ids)
                or (self.expression_state == "success" and not self.expression_observations)
                or (self.expression_state in {"no_result", "unavailable", "invalid"} and self.expression_observations)):
            raise ValueError("target profile lacks source context or honest lookup coverage")


def action_comparison(target: TargetHypothesis, hit: KnownDrugHit, profile: TargetProfileInput) -> str:
    desired = target.required_modulation.removeprefix("partial_")
    if hit.action_direction in {"mixed", "unknown"}:
        return "direction_unknown"
    if hit.action_direction != desired:
        return "opposed"
    if (target.causal_role == "direct" and desired == "increase"
            and profile.residual_functional_product != "present"
            and not profile.separately_supported_restoration):
        return "functional_product_not_established"
    return "direction_aligned_requires_independent_action_evidence"


def build_target_profile(target, profile, *, genotype_effect, candidate_decisions):
    if profile.target_id != target.target_id or not target.ready_for_drug_retrieval:
        raise ValueError("target profile has no qualified target rationale")
    hits = []
    for lookup in profile.known_drug_lookups:
        for hit in lookup.hits:
            if hit.candidate_id is not None and hit.candidate_id not in candidate_decisions:
                raise ValueError("known drug references an unevaluated candidate")
            decision = candidate_decisions.get(hit.candidate_id)
            comparison = action_comparison(target, hit, profile)
            if decision and decision["eligible"] and comparison != "direction_aligned_requires_independent_action_evidence":
                raise ValueError("eligible candidate contradicts source-bound action or functional product")
            hits.append({**asdict(hit), "source": lookup.source, "lookup_state": lookup.state,
                         "action_comparison": comparison,
                         "scientific_tier": decision["scientific_tier"] if decision else "not_assessed",
                         "research_eligible": decision["eligible"] if decision else False,
                         "eligibility_reason": decision["reason_codes"] if decision else ["no_verified_medication_candidate"],
                         "approval_from_known_drug_metadata": False})
    return {"target_id": target.target_id, "genotype_effect": genotype_effect,
            "required_modulation": target.required_modulation, "signed_rationale": asdict(target),
            "context": asdict(profile), "retrieved_drugs": hits,
            "eligible_candidate_ids": sorted({h["candidate_id"] for h in hits if h["research_eligible"]}),
            "expression_interpretation": "Reference expression is context, not patient product availability, rescue or safety.",
            "human_genetic_interpretation": "This context is optional. No genetic association or lifelong effect establishes a medication benefit.",
            "coverage_interpretation": "The source results have a bounded scope. No-hit, failed, truncated and excluded records remain distinct."}


def render_target_profiles(profiles):
    lines = ["## Intervention target profiles", ""]
    for profile in profiles:
        context = profile["context"]
        lines += ["### " + profile["target_id"], "",
                  f"Genotype effect: {profile['genotype_effect']}. Required modulation: {profile['required_modulation']}.", "",
                  f"Residual functional product: {context['residual_functional_product']}. {context['residual_product_context'].rstrip('.')}.", "",
                  "Searches for known medications:", ""]
        for record in context['known_drug_lookups']:
            lines += [f"- {record['source']}: {record['state']}."]
        lines += ['']
        for hit in profile["retrieved_drugs"]:
            lines += [f"- {hit['name']} ({hit['drug_id']}): {hit['reported_action'].rstrip('.')}.",
                      f"  Direction: {hit['action_comparison']}. Tier: {hit['scientific_tier']}. "
                      f"Research eligible: {str(hit['research_eligible']).lower()}.",
                      f"  Evidence: {', '.join(hit['evidence_ids'])}."]
        if not profile["retrieved_drugs"]:
            lines += ["The recorded searches found no usable medication result. This does not establish an exhaustive absence of therapies."]
        lines += ["", f"Baseline expression: {context['expression_state']}. " + profile["expression_interpretation"], ""]
        for row in context["expression_observations"]:
            lines += [f"- {row['tissue_or_cell_type']}: {row['result']}.",
                      f"  Assay: {row['assay']}. Population: {row['population']}.",
                      f"  Evidence: {', '.join(row['evidence_ids'])}."]
        lines += ["", "Human-genetic context: " + context["human_genetic_context"].rstrip('.') + ". " + profile["human_genetic_interpretation"], "",
                  "Hazards and conflicts:", ""]
        for statement in (*profile['signed_rationale']['harm_conditions'], *profile['signed_rationale']['conflicts']):
            lines += ['- ' + statement.rstrip('.') + '.']
        lines += ['']
    return "\n".join(lines)
