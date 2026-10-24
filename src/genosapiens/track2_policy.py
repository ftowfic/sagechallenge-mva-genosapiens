"""Strict source-matrix and G0-G7 policy validation."""

from __future__ import annotations

from dataclasses import dataclass
from datetime import date
import json
from pathlib import Path
from typing import Any, Mapping


SOURCE_SCHEMA = "genosapiens.track2-source-domain-matrix/v1"
DECISION_SCHEMA = "genosapiens.track2-decision-policy/v1"
EXPANSION_SCHEMA = "genosapiens.track2-public-expansion-policy/v1"
GATES = {f"G{number}" for number in range(8)}
TIERS = {"supported", "conditional", "exploratory", "excluded", "contraindicated"}
GATE_OUTCOMES = {"pass", "bounded_uncertainty", "not_evaluable", "fail", "contraindicating"}
ROUTES = {
    "direct_target_modulation",
    "product_restoration_or_stabilization",
    "transcript_rescue",
    "upstream_regulation",
    "downstream_compensation",
    "substrate_or_cofactor_modulation",
    "pathway_bypass",
    "qualified_target_drug_linkage",
    "intervention_literature",
    "lincs_cell_specific_expression",
}


class Track2PolicyError(ValueError):
    """Raised when a required policy name or transition is undefined."""


def _load(path: str | Path) -> dict[str, Any]:
    try:
        data = json.loads(Path(path).read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as error:
        raise Track2PolicyError(f"cannot load policy {path}: {error}") from error
    if not isinstance(data, dict):
        raise Track2PolicyError("policy root must be an object")
    return data


def _string_list(value: Any, field: str, *, allow_empty: bool = False) -> list[str]:
    if (
        not isinstance(value, list)
        or (not allow_empty and not value)
        or any(not isinstance(item, str) or not item for item in value)
        or len(value) != len(set(value))
    ):
        raise Track2PolicyError(f"{field} must be a unique string list")
    return value


@dataclass(frozen=True, slots=True)
class Track2Policies:
    source_matrix: Mapping[str, Any]
    decision_policy: Mapping[str, Any]


def validate_public_expansion_policy(
    data: Mapping[str, Any], source_names: set[str]
) -> None:
    required = {
        "date_range",
        "entity_types",
        "global_bounds",
        "mechanism_branches",
        "policy_version",
        "relations",
        "schema",
        "screening_profiles",
        "default_submission_routes",
    }
    if set(data) != required or data.get("schema") != EXPANSION_SCHEMA:
        raise Track2PolicyError("public expansion policy is incomplete or incompatible")
    defaults = _string_list(data.get("default_submission_routes"), "default_submission_routes")
    if set(defaults) - ROUTES:
        raise Track2PolicyError("default submission routes are undefined")
    entity_types = set(_string_list(data.get("entity_types"), "entity_types"))
    screening_profiles = set(
        _string_list(data.get("screening_profiles"), "screening_profiles")
    )
    _string_list(data.get("mechanism_branches"), "mechanism_branches")
    bounds = data.get("global_bounds")
    if not isinstance(bounds, Mapping) or set(bounds) != {
        "maximum_graph_depth",
        "maximum_results_per_seed_source_relation",
        "stopping",
    }:
        raise Track2PolicyError("public expansion global bounds are incomplete")
    for field in ("maximum_graph_depth", "maximum_results_per_seed_source_relation"):
        value = bounds[field]
        if not isinstance(value, int) or isinstance(value, bool) or value <= 0:
            raise Track2PolicyError(f"public expansion {field} must be positive")
    if bounds["stopping"] != "configured_bound_or_source_exhausted":
        raise Track2PolicyError("public expansion stopping criterion is undefined")
    date_range = data.get("date_range")
    if not isinstance(date_range, Mapping) or set(date_range) != {
        "end",
        "start",
        "undated_record_handling",
    }:
        raise Track2PolicyError("public expansion date range is incomplete")
    try:
        start = date.fromisoformat(date_range["start"])
        end = date.fromisoformat(date_range["end"])
    except (TypeError, ValueError) as error:
        raise Track2PolicyError("public expansion date range must use ISO dates") from error
    if start > end or date_range["undated_record_handling"] != "retain_and_flag":
        raise Track2PolicyError("public expansion date range is invalid")
    relations = data.get("relations")
    if not isinstance(relations, Mapping) or not relations:
        raise Track2PolicyError("public expansion relations are required")
    for name, relation in relations.items():
        if not isinstance(relation, Mapping) or set(relation) != {
            "from",
            "maximum_depth",
            "screening",
            "sources",
            "to",
        }:
            raise Track2PolicyError(f"public expansion relation {name} is incomplete")
        from_types = set(_string_list(relation["from"], f"{name}.from"))
        to_types = set(_string_list(relation["to"], f"{name}.to"))
        relation_sources = set(_string_list(relation["sources"], f"{name}.sources"))
        relation_screening = set(
            _string_list(relation["screening"], f"{name}.screening")
        )
        if (from_types | to_types) - entity_types:
            unknown = sorted((from_types | to_types) - entity_types)[0]
            raise Track2PolicyError(f"public expansion relation uses unknown entity: {unknown}")
        if relation_sources - source_names:
            unknown = sorted(relation_sources - source_names)[0]
            raise Track2PolicyError(f"public expansion relation uses unknown source: {unknown}")
        if relation_screening - screening_profiles:
            unknown = sorted(relation_screening - screening_profiles)[0]
            raise Track2PolicyError(f"public expansion relation uses unknown screening: {unknown}")
        depth = relation["maximum_depth"]
        if (
            not isinstance(depth, int)
            or isinstance(depth, bool)
            or depth <= 0
            or depth > bounds["maximum_graph_depth"]
        ):
            raise Track2PolicyError(f"public expansion relation {name} has invalid depth")


def load_public_expansion_policy(
    path: str | Path, source_matrix_path: str | Path
) -> dict[str, Any]:
    source_matrix = _load(source_matrix_path)
    validate_source_matrix(source_matrix)
    policy = _load(path)
    validate_public_expansion_policy(policy, set(source_matrix["sources"]))
    return policy


def validate_public_expansion_request(
    request: Mapping[str, Any], policy: Mapping[str, Any]
) -> None:
    required = {
        "date_end",
        "date_start",
        "from_entity_type",
        "graph_depth",
        "max_results",
        "relation_type",
        "screening",
        "seed_id",
        "source",
        "stopping",
        "to_entity_type",
    }
    if set(request) != required:
        raise Track2PolicyError("public expansion request is incomplete or unbounded")
    relation = policy["relations"].get(request.get("relation_type"))
    if not isinstance(relation, Mapping):
        raise Track2PolicyError("public expansion edge is not allowed")
    if (
        request.get("from_entity_type") not in relation["from"]
        or request.get("to_entity_type") not in relation["to"]
        or request.get("source") not in relation["sources"]
        or request.get("screening") not in relation["screening"]
    ):
        raise Track2PolicyError("public expansion edge is not allowed")
    depth = request.get("graph_depth")
    result_count = request.get("max_results")
    if (
        not isinstance(depth, int)
        or isinstance(depth, bool)
        or not 1 <= depth <= relation["maximum_depth"]
        or not isinstance(result_count, int)
        or isinstance(result_count, bool)
        or not 1
        <= result_count
        <= policy["global_bounds"]["maximum_results_per_seed_source_relation"]
    ):
        raise Track2PolicyError("public expansion request is unbounded")
    if request.get("date_start") != policy["date_range"]["start"] or request.get(
        "date_end"
    ) != policy["date_range"]["end"]:
        raise Track2PolicyError("public expansion request is outside the fixed date range")
    if request.get("stopping") != policy["global_bounds"]["stopping"]:
        raise Track2PolicyError("public expansion request is unbounded")
    if not isinstance(request.get("seed_id"), str) or not request["seed_id"]:
        raise Track2PolicyError("public expansion seed identity is required")


def validate_source_matrix(data: Mapping[str, Any]) -> set[str]:
    required = {
        "domains",
        "jurisdiction_exception",
        "policy_version",
        "primary_jurisdiction",
        "schema",
        "sources",
    }
    if set(data) != required or data.get("schema") != SOURCE_SCHEMA:
        raise Track2PolicyError("source-domain matrix is incomplete or incompatible")
    if data.get("primary_jurisdiction") != "US":
        raise Track2PolicyError("the Track 2 primary jurisdiction must be US")
    exception = data.get("jurisdiction_exception")
    if not isinstance(exception, Mapping) or set(exception) != {
        "allowed",
        "requires_current_marketing_status",
        "requires_manual_review",
        "requires_named_jurisdiction",
        "requires_official_regulator",
    } or any(value is not True for value in exception.values()):
        raise Track2PolicyError("named-jurisdiction exception requirements are incomplete")
    sources = data.get("sources")
    domains = data.get("domains")
    if not isinstance(sources, Mapping) or not sources:
        raise Track2PolicyError("source definitions are required")
    if not isinstance(domains, Mapping) or not domains:
        raise Track2PolicyError("source domains are required")
    source_names = set(sources)
    for name, source in sources.items():
        if not isinstance(source, Mapping) or set(source) != {"authority", "kind"}:
            raise Track2PolicyError(f"source {name} is incomplete")
        if any(not isinstance(value, str) or not value for value in source.values()):
            raise Track2PolicyError(f"source {name} values must be non-empty strings")
    for name, domain in domains.items():
        if not isinstance(domain, Mapping):
            raise Track2PolicyError(f"source domain {name} must be an object")
        required_fields = {
            "alternate",
            "freshness_days",
            "manual_import",
            "missing",
            "primary",
            "redistribution",
        }
        optional_fields = {"optional", "supporting_only"}
        if not required_fields.issubset(domain) or set(domain) - required_fields - optional_fields:
            raise Track2PolicyError(f"source domain {name} fields are incomplete or unknown")
        referenced: list[str] = []
        for field in ("primary", "alternate", "optional", "supporting_only"):
            if field in domain:
                referenced.extend(
                    _string_list(domain[field], f"{name}.{field}", allow_empty=True)
                )
        undefined = sorted(set(referenced) - source_names)
        if undefined:
            raise Track2PolicyError(
                f"source domain {name} references undefined required source: {undefined[0]}"
            )
        freshness = domain["freshness_days"]
        if freshness is not None and (
            not isinstance(freshness, int) or isinstance(freshness, bool) or freshness <= 0
        ):
            raise Track2PolicyError(f"{name}.freshness_days must be null or positive")
    return set(domains)


def validate_decision_policy(data: Mapping[str, Any], source_domains: set[str]) -> None:
    required = {
        "blocking_precedence",
        "candidate_generation_constraints",
        "evidence_maturity",
        "g5_exposure_policy",
        "gates",
        "input_definitions",
        "missing_data_default",
        "policy_version",
        "reason_codes",
        "review",
        "research_submission",
        "route_bounds",
        "route_completeness",
        "schema",
        "source_matrix_schema",
        "tier_transitions",
    }
    if set(data) != required or data.get("schema") != DECISION_SCHEMA:
        raise Track2PolicyError("decision policy is incomplete or incompatible")
    if data.get("source_matrix_schema") != SOURCE_SCHEMA:
        raise Track2PolicyError("decision policy references the wrong source matrix")
    if data.get("missing_data_default") != "fail_closed":
        raise Track2PolicyError("missing decision data must fail closed")
    inputs = data.get("input_definitions")
    gates = data.get("gates")
    reason_codes = set(_string_list(data.get("reason_codes"), "reason_codes"))
    if not isinstance(inputs, Mapping) or not inputs:
        raise Track2PolicyError("gate input definitions are required")
    if any(not isinstance(value, str) or not value for value in inputs.values()):
        raise Track2PolicyError("gate input definitions must be named types")
    if not isinstance(gates, Mapping) or set(gates) != GATES:
        raise Track2PolicyError("decision policy must define exactly G0-G7")
    for gate_id, gate in gates.items():
        if not isinstance(gate, Mapping) or set(gate) != {
            "name",
            "reason_codes",
            "required_inputs",
            "required_sources",
        }:
            raise Track2PolicyError(f"gate {gate_id} fields are incomplete or unknown")
        gate_inputs = set(_string_list(gate["required_inputs"], f"{gate_id}.required_inputs"))
        undefined_inputs = sorted(gate_inputs - set(inputs))
        if undefined_inputs:
            raise Track2PolicyError(
                f"gate {gate_id} references undefined gate input: {undefined_inputs[0]}"
            )
        gate_sources = set(_string_list(gate["required_sources"], f"{gate_id}.required_sources"))
        undefined_sources = sorted(gate_sources - source_domains)
        if undefined_sources:
            raise Track2PolicyError(
                f"gate {gate_id} references undefined required source: {undefined_sources[0]}"
            )
        undefined_reasons = sorted(
            set(_string_list(gate["reason_codes"], f"{gate_id}.reason_codes"))
            - reason_codes
        )
        if undefined_reasons:
            raise Track2PolicyError(
                f"gate {gate_id} references undefined reason code: {undefined_reasons[0]}"
            )
    transitions = data.get("tier_transitions")
    if not isinstance(transitions, Mapping) or set(transitions) != TIERS:
        raise Track2PolicyError("decision policy has an undefined state transition")
    serialized = json.dumps(transitions, sort_keys=True)
    for token in re_tokens(serialized):
        if token.startswith("G") and token[1:].isdigit() and token not in GATES:
            raise Track2PolicyError(f"state transition references undefined gate: {token}")
        if token in {"pass", "bounded_uncertainty", "not_evaluable", "fail", "contraindicating"}:
            if token not in GATE_OUTCOMES:
                raise Track2PolicyError(f"state transition references undefined outcome: {token}")
    conditional = transitions["conditional"]
    if set(conditional.get("forbidden_reasons", [])) - reason_codes:
        raise Track2PolicyError("state transition references undefined reason code")
    precedence = _string_list(data.get("blocking_precedence"), "blocking_precedence")
    for item in precedence:
        try:
            gate_id, outcome = item.split(":", 1)
        except ValueError as error:
            raise Track2PolicyError("blocking precedence entries must be gate:outcome") from error
        if gate_id not in GATES or outcome not in GATE_OUTCOMES:
            raise Track2PolicyError("blocking precedence has an undefined state transition")
    routes = data.get("route_bounds")
    if not isinstance(routes, Mapping) or set(routes) != ROUTES:
        raise Track2PolicyError("route policy must define every bounded route")
    for route, bounds in routes.items():
        if not isinstance(bounds, Mapping) or set(bounds) != {
            "lookback_years", "max_depth", "max_results", "screening", "stopping"
        }:
            raise Track2PolicyError(f"route {route} bounds are incomplete")
        for field in ("lookback_years", "max_depth", "max_results"):
            value = bounds[field]
            if not isinstance(value, int) or isinstance(value, bool) or value <= 0:
                raise Track2PolicyError(f"route {route}.{field} must be positive")
    completeness = data.get("route_completeness")
    if not isinstance(completeness, Mapping) or set(completeness) != {
        "affected_candidate_tier_ceiling",
        "coverage_limiting_states",
        "limited_coverage_claim",
        "limited_no_candidate_claim",
        "milestone_blocking_states",
        "scope",
    }:
        raise Track2PolicyError("route completeness policy is incomplete")
    if (
        set(_string_list(completeness["milestone_blocking_states"], "milestone_blocking_states"))
        != {"deferred", "unattempted", "unimplemented"}
        or set(_string_list(completeness["coverage_limiting_states"], "coverage_limiting_states"))
        != {"truncated", "unavailable"}
        or completeness["affected_candidate_tier_ceiling"] != "exploratory"
        or completeness["limited_coverage_claim"] != "incomplete"
        or completeness["limited_no_candidate_claim"] != "not_exhaustive"
        or completeness["scope"] != "predeclared_selected_routes"
    ):
        raise Track2PolicyError("route completeness policy has undefined states")
    research = data.get("research_submission")
    if research != {
        "required_pass": ["G0", "G1", "G3", "G4", "G6", "G7"],
        "allowed_mechanism_outcomes": ["pass", "bounded_uncertainty"],
        "allowed_exploratory_reasons": ["target_site_exposure_unknown", "rescue_magnitude_unknown"],
        "required_review_domains": sorted([
            "approval", "mechanism", "drug_action", "current_warnings", "negative_evidence",
            "harmful_alternatives", "target_liabilities", "claims", "validation",
        ]),
        "allowed_gap_fields": ["target_site_exposure", "rescue_magnitude", "pediatric_applicability"],
        "tier_override_allowed": False,
        "lincs_alternative": "lincs-research-eligibility/v1:repeated_consistent",
    }:
        raise Track2PolicyError("research-submission policy has undefined requirements or transitions")
    constraints = data.get("candidate_generation_constraints")
    if not isinstance(constraints, Mapping) or set(constraints) != {
        "cannot_be_sole_primary_support",
        "external_computational_tier_ceiling",
        "prohibited_internal_methods",
        "required_external_metadata",
        "score_substitution_allowed",
    }:
        raise Track2PolicyError("candidate-generation constraints are incomplete")
    if (
        set(_string_list(constraints["prohibited_internal_methods"], "prohibited_internal_methods"))
        != {
            "broad_signature_reversal",
            "generalized_network_proximity",
            "language_model_candidate_generation",
        }
        or set(_string_list(constraints["required_external_metadata"], "required_external_metadata"))
        != {
            "input_provenance",
            "model_name",
            "model_parameters",
            "model_version",
            "output_provenance",
        }
        or constraints["external_computational_tier_ceiling"] != "exploratory"
        or constraints["cannot_be_sole_primary_support"] is not True
        or constraints["score_substitution_allowed"] is not False
    ):
        raise Track2PolicyError("candidate-generation constraints are undefined")
    exposure_policy = data.get("g5_exposure_policy")
    if not isinstance(exposure_policy, Mapping) or set(exposure_policy) != {
        "bounded_secondary_uncertainties",
        "material_unknowns",
        "minimum_potency_exposure_ratio",
    }:
        raise Track2PolicyError("G5 exposure policy is incomplete")
    minimum_ratio = exposure_policy["minimum_potency_exposure_ratio"]
    if (
        not isinstance(minimum_ratio, (int, float))
        or isinstance(minimum_ratio, bool)
        or minimum_ratio <= 0
        or set(
            _string_list(
                exposure_policy["bounded_secondary_uncertainties"],
                "bounded_secondary_uncertainties",
            )
        )
        != {
            "interindividual_variability",
            "protein_binding_estimate",
            "tissue_partition_estimate",
        }
        or set(_string_list(exposure_policy["material_unknowns"], "material_unknowns"))
        != {"rescue_magnitude", "target_site_exposure"}
    ):
        raise Track2PolicyError("G5 exposure policy has undefined states")


def re_tokens(value: str) -> tuple[str, ...]:
    return tuple(
        token.strip('"{},:[] ') for token in value.replace("\n", " ").split()
    )


def load_track2_policies(
    source_matrix_path: str | Path,
    decision_policy_path: str | Path,
) -> Track2Policies:
    source = _load(source_matrix_path)
    domains = validate_source_matrix(source)
    decision = _load(decision_policy_path)
    validate_decision_policy(decision, domains)
    return Track2Policies(source, decision)
