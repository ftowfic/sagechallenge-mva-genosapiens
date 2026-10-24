"""Single offline scientific composer with source-bound public G0 admission.

Reviews describe scientific judgments. Gates, coverage and readiness are replayed
from those judgments; accepting the input never approves a medication.
"""
from dataclasses import asdict, fields, replace
from datetime import datetime
import json
from pathlib import Path
from typing import Any, Mapping
from genosapiens.track2_scientific_storage import json_identity,file_identity
from genosapiens.track2_scientific_reader import read_scientific

from genosapiens.track2_codec import decode
from genosapiens.track2_candidate_union import build_candidate_union
from genosapiens.track2_endpoint import select_endpoint
from genosapiens.track2_exposure import ActivityExposureComparison
from genosapiens.track2_freeze import (SynthesisRequest, load_frozen_public_envelope,
    run_frozen_synthesis, validate_public_research_scope, freeze_public_envelope)
from genosapiens.track2_g5 import ExposureUncertainty, evaluate_g5
from genosapiens.track2_g6 import MVAHazardAssessment, evaluate_mva_g6
from genosapiens.track2_gates import evaluate_g0_g4
from genosapiens.track2_manual_import import ManualEvidenceImport
from genosapiens.track2_mechanism import assess_mechanism_readiness
from genosapiens.track2_milestone import assess_research_eligibility, assess_route_milestone
from genosapiens.track2_no_candidate import build_objective_candidate_outcome
from genosapiens.track2_policy import load_public_expansion_policy, load_track2_policies
from genosapiens.track2_population_applicability import assess_population_applicability
from genosapiens.track2_public_query import build_public_seed_registry
from genosapiens.track2_published_genotype import (canonical, sha256, _source_file,
    load_published_genotype, render_input_dossier)
from genosapiens.track2_input import (load_analysis_input, input_evidence,
    assessed_report_origins as report_origins, assessed_mechanism_inputs as source_mechanism_inputs)
from genosapiens.track2_report_claims import ReportClaim, build_report_claim_bundle
from genosapiens.track2_report_tables import (ReportAssessment, ReportTableInputs,
    build_report_tables, report_tables_as_json, report_tables_as_markdown)
from genosapiens.track2_routes import RouteAssessment, build_route_portfolio
from genosapiens.track2_assessment import assess_plan
from genosapiens.track2_target_hypothesis import build_target_hypothesis
from genosapiens.track2_target_profile import (TargetProfileInput, action_comparison,
    build_target_profile, render_target_profiles)
from genosapiens.track2_tiers import decide_candidate_tier
from genosapiens.track2_validation_plan import ValidationPlan, review_structural_g7

ROOT = Path(__file__).resolve().parents[2]
SCHEMA = "genosapiens.track2-synthesis-plan/v4"
NATIVE_SCHEMA = 'genosapiens.track2-synthesis-plan/v5'
PLAN_FIELDS = {"schema", "classification", "input_origin", "tables", "mechanism_inputs",
    "candidate_inputs", "routes", "requests", "objective_candidates", "candidate_route_dependencies",
    "profile_inputs", "claims", "methods", "evidence_imports", "scope_sha256", "lincs", "mechanism_assessment", "frozen_assessment", "assessment_scope"}


def assessment_bindings(plan):
    return {**({'mechanism_assessment': json_identity(plan['mechanism_assessment'])['sha256']}
            if plan['mechanism_assessment'] is not None else {}),
            **({'frozen_assessment':plan['frozen_assessment']['sha256']}
               if plan['frozen_assessment'] is not None else {}),
            **({'discovery_resources':plan['exploratory_discovery']['resources']['sha256'],
                'discovery_policy':plan['exploratory_discovery']['policy']['sha256']}
               if plan.get('schema')==NATIVE_SCHEMA else {})}


class Track2SynthesisError(ValueError):
    def __init__(self, message="scientific input and replayed decision disagree"):
        super().__init__(message)


def code_identity():
    """Bind study production without unrelated operator or finalization code."""
    from genosapiens.track2_construction_dependencies import python_dependencies
    from genosapiens.track2_report_finalization import VALIDATORS
    dependencies = python_dependencies(ROOT, [('genosapiens.track2_synthesis', 'synthesize')],
        ignored=VALIDATORS | {('genosapiens.track2_synthesis', 'code_identity')})
    for name in ('read_track2_lincs.R', 'render_track2_lincs.R'):
        dependencies['script:'+name] = file_identity(ROOT/'scripts'/name)['sha256']
    return sha256(canonical(dependencies))


def plan_acceptance_identity():
    """Bind the technical plan decision independently of report presentation."""
    from genosapiens.track2_construction_dependencies import python_dependencies
    return sha256(canonical(python_dependencies(ROOT,
        [('genosapiens.track2_synthesis', 'plan_acceptance')])))


def plan_acceptance(plan, input_sha256, *, timestamp, image):
    from genosapiens.track2_adjudication import automatic_decision
    scope=plan['assessment_scope'] or {sid:{'required':True} for sid in plan['mechanism_inputs']}
    selected={sid for sid,row in scope.items() if row['required']}
    mechanisms={sid:assess_mechanism_readiness(row) for sid,row in plan['mechanism_inputs'].items()}
    if plan.get('schema')==NATIVE_SCHEMA:
        native=plan['exploratory_discovery']
        objectives={o['scenario_id'] for o in native['objectives']}
        checks={k:{'state':'accepted' if value else 'unresolved','reasons':[reason]} for k,value,reason in (
            ('integrity',True,'validated_native_plan_and_source_bindings'),
            ('mechanisms',bool(selected) and all(mechanisms[s].permits_planning for s in selected),'accepted_required_mechanisms'),
            ('hypotheses',bool(selected) and selected==objectives and all(o['kind']=='restoration_or_compensation_hypothesis'
                for o in native['objectives']),'native_restoration_hypotheses'),
            ('prospective_scope',set(scope)==set(mechanisms),'all_scenarios_retained_before_results'))}
        return automatic_decision(domain='discovery_plan',subject=plan,
            evidence={'input_sha256':input_sha256,'assessment':plan['mechanism_assessment'],
                'scope':plan['assessment_scope'],'discovery':native},checks=checks,
            execution={'code_sha256':plan_acceptance_identity(),'image':image},timestamp=timestamp)
    bound={r['scenario_id'] for r in (plan['frozen_assessment'] or {}).get('scenario_bindings',[])}
    checks={k:{'state':'accepted' if value else 'unresolved','reasons':[reason]} for k,value,reason in (
        ('integrity',True,'validated_plan_and_source_bindings'),
        ('mechanisms',bool(selected) and all(mechanisms[s].permits_planning for s in selected),'accepted_required_mechanisms'),
        ('applicability',bool(selected) and selected==bound,'applicable_required_objectives'),
        ('prospective_scope',set(scope)==set(mechanisms),'all_scenarios_retained_before_results'))}
    return automatic_decision(domain='plan',subject=plan,
        evidence={'input_sha256':input_sha256,'assessment':plan['mechanism_assessment'],
                  'scope':plan['assessment_scope'],'frozen_assessment':plan['frozen_assessment']},
        checks=checks,execution={'code_sha256':plan_acceptance_identity(),'image':image},timestamp=timestamp)


def load_plan_acceptance(path, plan, bundle):
    target=Path(path).with_name('plan-acceptance.json')
    if not target.exists(): return None
    raw=target.read_bytes();record=json.loads(raw)
    if (bundle.data['input_origin']!='diagnosis_guided_vcf' or raw!=canonical(record)
        or record!=plan_acceptance(plan,bundle.input_sha256,timestamp=record['timestamp'],
                                  image=bundle.data['builder']['runtime']['image'])):
        raise Track2SynthesisError('Automatic plan acceptance changed or has stale sources.')
    return record


def load_plan(bundle_root, plan_path, *, require_review=True, allow_synthetic=False, discovery_reference_root='/discovery-resources', retained_science=None):
    root, path = Path(bundle_root), Path(plan_path)
    plan = read_scientific(path)
    plan_digest = file_identity(path)['sha256']
    native=plan.get('schema')==NATIVE_SCHEMA
    if not native:
        # The existing public study contract keeps its current typed composer.
        # Native discovery uses disk-backed scientific records throughout.
        plan=json.loads(path.read_bytes())
    bundle = load_analysis_input(root,input_origin=plan.get('input_origin'),allow_synthetic=allow_synthetic,
        require_review=require_review or plan.get('input_origin')!='diagnosis_guided_vcf')
    if (set(plan) != PLAN_FIELDS | ({'exploratory_discovery'} if native else set()) or plan["schema"] not in {SCHEMA,NATIVE_SCHEMA}
            or plan["classification"] != bundle.data['classification']
            or plan["input_origin"] != bundle.data["input_origin"]):
        raise Track2SynthesisError("unsupported synthesis plan or origin")
    if (plan['input_origin']=='diagnosis_guided_vcf') != (plan['mechanism_assessment'] is not None):
        raise Track2SynthesisError('Diagnosis planning requires its downstream mechanism assessment')
    source_mechanism_inputs(bundle,plan['mechanism_inputs'],draft=not require_review,
                            assessment=plan['mechanism_assessment'])
    if plan['assessment_scope'] is not None:
        from genosapiens.track2_mechanism_review import validate_scope
        validate_scope(plan['mechanism_assessment'],plan['assessment_scope'])
    if native:
        from genosapiens.track2_exploratory_planning import load_plan_resources,validate_native_plan
        resources=load_plan_resources(path,plan['exploratory_discovery'],discovery_reference_root)
        try:
            validate_native_plan(plan,resources)
        finally:
            resources['connection'].close()
    if plan['frozen_assessment'] is not None:
        from genosapiens.track2_scientific_planning import frozen_bindings
        frozen_bindings(path.parent,plan,bundle)
        if plan['lincs'] is not None:
            raise Track2SynthesisError('Frozen assessment and fresh LINCS execution are separate input modes')
    query_ids = set()
    sources = {**bundle.source_hashes, **assessment_bindings(plan)}
    if plan['lincs'] is not None:
        from genosapiens.track2_lincs_integration import bound_inputs
        from genosapiens.track2_report_science_reuse import retained_study_science
        retained=retained_study_science(retained_science,ROOT,plan,bundle.input_sha256)
        if retained is None:
            _, _, _, lincs_bindings = bound_inputs(root, plan['lincs'])
        else:
            lincs_bindings=retained[2]
        sources.update(lincs_bindings)
    evidence = {e["evidence_id"] for e in input_evidence(bundle)}
    if plan['mechanism_assessment'] is not None:
        additions = {e['evidence_id'] for e in plan['mechanism_assessment']['evidence']}
        if evidence & additions:
            raise Track2SynthesisError('Assessment evidence duplicates construction evidence')
        evidence.update(additions)
    for source in plan["evidence_imports"]:
        if set(source) != {"record", "path", "evidence_ids", "license_or_terms", "redistribution", "query_ids"}:
            raise Track2SynthesisError("source import contract is incomplete")
        record = decode(ManualEvidenceImport, source["record"])
        source_root=path.parent if plan['input_origin']=='diagnosis_guided_vcf' else root
        digest = sha256(_source_file(source_root, source["path"]).read_bytes())
        if digest != record.document_checksum or not source["license_or_terms"] or not source["redistribution"]:
            raise Track2SynthesisError("source import checksum or terms missing")
        if record.import_id in sources or not source["evidence_ids"] or evidence & set(source["evidence_ids"]):
            raise Track2SynthesisError("duplicate source or evidence identity")
        query_ids.update(source["query_ids"])
        sources[record.import_id] = digest
        evidence.update(source["evidence_ids"])
    attempted_queries = {q for row in plan["routes"] if row["attempt"] and row["attempt"]["executed"]
                         for q in row["attempt"]["query_ids"]}
    if attempted_queries - query_ids:
        raise Track2SynthesisError("route attempts must resolve to retained source imports")
    if native:
        from typing import get_type_hints,get_args
        from genosapiens.track2_codec import json_value
        models=get_type_hints(ReportTableInputs)
        if set(plan['tables'])!=set(models):raise Track2SynthesisError('Native input table fields changed.')
        catalog=set()
        for field,model in models.items():
            for row in plan['tables'][field]:
                record=decode(get_args(model)[0],json_value(row))
                if field=='evidence_catalog':catalog.add(record.evidence_id)
    else:
        tables = decode(ReportTableInputs, plan["tables"])
        catalog = {e.evidence_id for e in tables.evidence_catalog}
    if catalog != evidence:
        raise Track2SynthesisError("report catalog must resolve to the actual input or source imports")
    if plan['lincs'] is not None and 'lincs-evidence:'+plan['lincs']['sources_sha256'] not in evidence:
        raise Track2SynthesisError('LINCS evidence must be bound to the existing reviewed source catalog')
    _resolve_evidence(plan, evidence)
    review,review_hash = _load_plan_review(path,plan_digest,sources,bundle,
        require_review=require_review,allow_synthetic=allow_synthetic)
    acceptance=load_plan_acceptance(path,plan,bundle)
    bindings = {"input": bundle.input_sha256, "input_review": bundle.review_sha256,
                "plan": plan_digest, "plan_review": review_hash, "code": native_plan_identity() if native else code_identity(),
                "input_schema": sha256((ROOT / ('schemas/track2/published-patient-genotype.schema.json'
                    if plan['input_origin']=='published_patient_genotype' else 'schemas/track2/diagnosis-guided-package.schema.json')).read_bytes()),
                **sources}
    if plan['input_origin']=='diagnosis_guided_vcf':
        bindings['automatic_input']=bundle.automatic_acceptance_sha256
    if acceptance is not None: bindings['automatic_plan']=sha256(canonical(acceptance))
    return bundle, plan, review, bindings


def native_plan_identity():
    from genosapiens.track2_construction_dependencies import python_dependencies
    return sha256(canonical(python_dependencies(ROOT,[
        ('genosapiens.track2_mechanism_assessment','validate_assessment'),
        ('genosapiens.track2_exploratory_planning','validate_native_plan')])))


def _load_plan_review(path,plan_digest,sources,bundle,*,require_review,allow_synthetic):
    review_path=path.with_name('plan-review.json')
    if (not require_review and bundle.data['input_origin']=='diagnosis_guided_vcf'
            and not review_path.exists()):
        return None,None
    review_bytes = review_path.read_bytes()
    review = json.loads(review_bytes)
    expected_fields = {"schema", "reviewer", "reviewer_type", "independence", "decision",
                       "rationale", "reviewed_at", "plan_sha256", "source_hashes"}
    if (set(review) != expected_fields or review["schema"] != "genosapiens.track2-plan-review/v1"
            or review["plan_sha256"] != plan_digest or review["source_hashes"] != sources
            or not review["reviewer"] or not review["rationale"]
            or review["reviewer_type"] not in {"human", "ai_assisted_extraction", "synthetic_fixture"}
            or review["independence"] not in {"independent", "self_review_disclosed"}
            or review["decision"] not in {"accepted", "pending_human_review", "rejected"}
            or datetime.fromisoformat(review["reviewed_at"].replace("Z", "+00:00")).tzinfo is None):
        raise Track2SynthesisError("source-bound plan review is incomplete or stale")
    synthetic_review=(review['reviewer_type']=='synthetic_fixture' and allow_synthetic
        and bundle.data.get('evidence_origin')=='independently_synthetic')
    if review['reviewer_type']=='synthetic_fixture' and not synthetic_review:
        raise Track2SynthesisError('synthetic review cannot admit actual scientific inputs')
    accepted = (review["reviewer_type"] == "human" or synthetic_review) and review["decision"] == "accepted"
    if (require_review and not accepted) or (review["decision"] == "accepted" and not accepted):
        raise Track2SynthesisError("actual scientific plan review is required")
    return review,sha256(review_bytes)


def _resolve_evidence(value, evidence):
    if isinstance(value, dict):
        for key, child in value.items():
            if key.endswith("evidence_id") and child is not None and child not in evidence:
                raise Track2SynthesisError("unresolved evidence reference: " + str(child))
            if key.endswith("evidence_ids") and set(child) - evidence:
                raise Track2SynthesisError("unresolved evidence references")
            _resolve_evidence(child, evidence)
    elif isinstance(value, list):
        for child in value:
            _resolve_evidence(child, evidence)


def _scenario_readiness(plan,bundle,policy,*,draft=False):
    raw_tables = decode(ReportTableInputs, plan["tables"])
    origins={o.scenario_id:o for o in source_mechanism_inputs(bundle,plan['mechanism_inputs'],draft=draft,
        assessment=plan['mechanism_assessment'])}
    records = (*raw_tables.mechanism_records, *raw_tables.endpoint_records,
               *raw_tables.target_records, *raw_tables.candidate_records)
    if (any(row.context.origin != origins.get(row.context.origin.scenario_id) or row.context.policy_version != policy["policy_version"]
                   for row in records)
            or {row.readiness.scenario_id for row in raw_tables.mechanism_records} != set(origins)
            or set(plan["candidate_inputs"]) != {row.candidate.active_ingredient_id for row in raw_tables.candidate_records}):
        raise Track2SynthesisError("report origin must match admitted input and review")
    readiness = {scenario: assess_mechanism_readiness(facts) for scenario,facts in plan['mechanism_inputs'].items()}
    if any(row.readiness != readiness[row.context.origin.scenario_id] for row in raw_tables.mechanism_records):
        raise Track2SynthesisError("mechanism readiness was not replayed")
    return raw_tables,readiness


def _compose(plan, bundle, frozen, policy, *, draft=False, scenario_id=None):
    if plan['schema']==NATIVE_SCHEMA:
        from genosapiens.track2_discovery_report import native_composition
        return native_composition(plan,frozen)
    raw_tables,readiness=_scenario_readiness(plan,bundle,policy,draft=draft)
    diagnosis=plan['input_origin']=='diagnosis_guided_vcf'
    if diagnosis and scenario_id is None and len(readiness)>1:
        from genosapiens.track2_report_tables import Track2ReportTables, TABLE_NAMES
        from genosapiens.track2_assessment import assess_scope
        branches=[_compose(plan,bundle,frozen,policy,draft=draft,scenario_id=sid) for sid in sorted(readiness)]
        tables=Track2ReportTables(raw_tables.evidence_catalog,(policy['policy_version'],),
            **{name:tuple(row for branch in branches for row in getattr(branch[0],name)) for name in TABLE_NAMES})
        outcomes=tuple(row for branch in branches for row in branch[1])
        milestones={key:value for branch in branches for key,value in branch[2].items()}
        profiles=tuple({**row,'scenario_id':sid} for sid,branch in zip(sorted(readiness),branches) for row in branch[3])
        claims=build_report_claim_bundle(tuple(decode(ReportClaim,c) for c in plan['claims']),report_tables=tables)
        assessment=assess_scope([s for branch in branches for s in branch[5].scenarios],
            [o for branch in branches for o in branch[5].objectives])
        return tables,outcomes,milestones,profiles,claims,assessment,[r for branch in branches for r in branch[6]]
    if scenario_id is not None:
        import copy
        plan=copy.deepcopy(plan)
        raw_tables=replace(raw_tables,**{name:tuple(row for row in getattr(raw_tables,name)
            if row.context.origin.scenario_id==scenario_id) for name in
            ('mechanism_records','endpoint_records','target_records','candidate_records')})
        readiness={scenario_id:readiness[scenario_id]}
        if plan['assessment_scope'] is not None:
            plan['assessment_scope']={scenario_id:plan['assessment_scope'][scenario_id]}
        objectives={r.objective.objective_id for r in raw_tables.target_records}
        candidates={r.candidate.active_ingredient_id for r in raw_tables.candidate_records}
        targets={r.target.target_id for r in raw_tables.target_records}
        plan['objective_candidates']={k:v for k,v in plan['objective_candidates'].items() if k in objectives}
        plan['candidate_inputs']={k:v for k,v in plan['candidate_inputs'].items() if k in candidates}
        plan['candidate_route_dependencies']={k:[d for d in v if d[0] in objectives]
            for k,v in plan['candidate_route_dependencies'].items() if k in candidates}
        plan['profile_inputs']=[p for p in plan['profile_inputs'] if p['target_id'] in targets]
        plan['routes']=[r for r in plan['routes'] if r['objective_id'] in objectives]
        plan['requests']=[r for r in plan['requests'] if r['objective_id'] in objectives]
    if diagnosis and not bundle.input_report_ready and (
            raw_tables.endpoint_records or raw_tables.target_records or raw_tables.candidate_records
            or plan['requests'] or plan['routes'] or plan['objective_candidates'] or plan['lincs'] is not None):
        raise Track2SynthesisError('unreviewed diagnosis reports require the input-only plan')
    targets = {row.target.target_id: row.target for row in raw_tables.target_records}
    profile_inputs = {p.target_id: p for p in (decode(TargetProfileInput, p) for p in plan["profile_inputs"])}
    if set(targets) != set(profile_inputs) or len(profile_inputs) != len(plan["profile_inputs"]):
        raise Track2SynthesisError("every qualified target requires one source-bound profile")
    for record in raw_tables.mechanism_records:
        if (record.readiness.exact_disease_validity == "exact_supported"
                and record.context.subtype_applicability != "exact_subtype"):
            raise Track2SynthesisError()
    endpoints = {item.endpoint.endpoint_id: item.endpoint for item in raw_tables.endpoint_records}
    for record in raw_tables.endpoint_records:
        endpoint = record.endpoint
        replayed = select_endpoint(**{item.name: getattr(endpoint, item.name)
            for item in fields(endpoint) if item.name not in {"relevance", "modifiability"}})
        if replayed != endpoint:
            raise Track2SynthesisError()
    for record in raw_tables.target_records:
        target, objective = record.target, record.objective
        if not readiness[record.context.origin.scenario_id].permits_planning:
            raise Track2SynthesisError('direction-dependent targets require a supported mechanism')
        endpoint = endpoints[objective.endpoint_id]
        if (objective.direction != endpoint.expected_direction
                or objective.tissue != endpoint.tissue or objective.stage != endpoint.stage
                or endpoint.relevance != "supported"
                or not all((objective.desired_change, objective.compartment, objective.assumptions,
                    objective.harm_conditions, objective.mechanism_review_id, objective.endpoint_review_id,
                    objective.evidence_ids))):
            raise Track2SynthesisError()
        replayed = build_target_hypothesis(objective=objective, **{
            item.name: getattr(target, item.name) for item in fields(target)
            if item.name not in {"objective_id", "ready_for_drug_retrieval"}})
        if replayed != target:
            raise Track2SynthesisError()
    for record in raw_tables.candidate_records:
        if build_candidate_union(record.candidate.rationales).candidates != (record.candidate,):
            raise Track2SynthesisError()
    assessments = [decode(RouteAssessment, item) for item in plan["routes"]]
    objectives = plan["objective_candidates"]
    candidate_ids = {row.candidate.active_ingredient_id for row in raw_tables.candidate_records}
    if ((not objectives and not diagnosis) or set().union(*(set(ids) for ids in objectives.values())) != candidate_ids
            or set(plan["candidate_route_dependencies"]) != candidate_ids):
        raise Track2SynthesisError()
    if diagnosis and set(objectives)!={r.objective.objective_id for r in raw_tables.target_records}:
        raise Track2SynthesisError('diagnosis-guided objectives require qualified source-bound targets')
    requests = {}
    affected = set()
    coverage = []
    coverage_by_dependency = {}
    for item in plan["requests"]:
        if set(item) != {"objective_id", "request"}:
            raise Track2SynthesisError()
        if item['objective_id'] not in objectives:
            raise Track2SynthesisError('query has no qualified intervention objective')
        request = decode(SynthesisRequest, item["request"])
        key = (item["objective_id"], request.route)
        if key in requests:
            raise Track2SynthesisError()
        result = run_frozen_synthesis(request, frozen, lambda _: True)
        requests[key] = result
        coverage_by_dependency[key] = result.coverage
        coverage.extend(result.coverage)
        if not result.authorized:
            affected.add(key)
    expected = {(item.objective_id, item.route) for item in assessments
                if item.route in frozen["public_scope"]["selected_routes"] and item.attempt is not None}
    if set(requests) != expected:
        raise Track2SynthesisError()
    assessments = [replace(item, attempt=replace(item.attempt, state="unavailable"))
                   if (item.objective_id, item.route) in affected else item for item in assessments]
    portfolios = {objective: build_route_portfolio(
        (objective,), (item for item in assessments if item.objective_id == objective), policy,
        public_scope=frozen["public_scope"],
    ) for objective in objectives}
    dependencies = {key: tuple(tuple(item) for item in items)
                    for key, items in plan["candidate_route_dependencies"].items()}
    limited = {(item.objective_id, item.route) for item in assessments if item.attempt is not None
               and item.attempt.state in policy["route_completeness"]["coverage_limiting_states"]}
    milestones = {objective: assess_route_milestone(
        portfolio, candidate_route_dependencies={key: dependencies[key] for key in objectives[objective]},
        decision_policy=policy,
    ) for objective, portfolio in portfolios.items()}
    decisions, gates_by_id, reviews, validation, candidates = {}, {}, {}, {}, []
    for row in raw_tables.candidate_records:
        candidate_id = row.candidate.active_ingredient_id
        origin = row.context.origin
        facts = plan["candidate_inputs"][candidate_id]
        if set(facts) != {"g0_g4", "g5", "g6", "g7", "evidence_gaps", "resolution_experiments"}:
            raise Track2SynthesisError()
        g0g4 = dict(facts["g0_g4"])
        automatic=next((s for s in (plan['mechanism_assessment'] or {}).get('scenarios',[])
            if s['scenario_id']==origin.scenario_id and 'applicability' in s),None)
        if automatic is not None:
            g0g4={**g0g4,'input_values':{**g0g4['input_values'],
                'mechanism_readiness':readiness[origin.scenario_id].readiness_state,
                'translational_maturity':'computational_inference','lesion_qualification':'unknown',
                'independent_rescue_lineage':'missing'}}
        if g0g4["input_values"]["input_integrity"] != "valid":
            raise Track2SynthesisError()
        if g0g4["input_values"]["mechanism_readiness"] != readiness[origin.scenario_id].readiness_state:
            raise Track2SynthesisError()
        if (g0g4["input_values"]["exact_subtype_applicability"] != row.context.subtype_applicability
                or g0g4["input_values"]["translational_maturity"] != row.context.translational_maturity):
            raise Track2SynthesisError()
        for target_id in row.target_hypothesis_ids:
            target = targets[target_id]
            profile = profile_inputs[target_id]
            if (target.causal_role == "direct" and target.required_modulation in {"increase", "partial_increase"}
                    and profile.residual_functional_product != "present"
                    and not profile.separately_supported_restoration
                    and g0g4["input_values"]["drug_action_direction"] == "aligned"):
                raise Track2SynthesisError("absent or unknown functional product cannot establish drug action")
            matches = [hit for lookup in profile.known_drug_lookups for hit in lookup.hits
                       if hit.candidate_id == candidate_id]
            if any(action_comparison(target, hit, profile) !=
                   "direction_aligned_requires_independent_action_evidence" for hit in matches):
                if g0g4["input_values"]["drug_action_direction"] == "aligned":
                    raise Track2SynthesisError("candidate action contradicts the source-bound target profile")
        population = assess_population_applicability(
            row.population_applicability.child_applicability.dimensions,
            row.population_applicability.general_mva_applicability.dimensions,
        )
        if population != row.population_applicability:
            raise Track2SynthesisError()
        g0g4["evidence_ids_by_gate"] = {key: tuple(ids) for key, ids in g0g4["evidence_ids_by_gate"].items()}
        first = evaluate_g0_g4(**g0g4, decision_policy=policy).results
        g5 = dict(facts["g5"])
        g5["comparison"] = decode(ActivityExposureComparison | None, g5["comparison"])
        if (g5["comparison"] is not None and g5["comparison"].active_ingredient_id != candidate_id
                or g5["time_above_threshold_supported"] is not None
                and type(g5["time_above_threshold_supported"]) is not bool):
            raise Track2SynthesisError()
        g5["uncertainties"] = tuple(decode(ExposureUncertainty, item) for item in g5["uncertainties"])
        g5["evidence_ids"] = tuple(g5["evidence_ids"])
        exposure = evaluate_g5(**g5, decision_policy=policy)
        g6 = dict(facts["g6"])
        if type(g6["claimed_cancer_prevention_from_proximal_biomarker"]) is not bool:
            raise Track2SynthesisError()
        if (g6["child_applicability_state"] != population.child_applicability.state
                or g6["general_mva_applicability_state"] != population.general_mva_applicability.state):
            raise Track2SynthesisError()
        g6["assessments"] = tuple(decode(MVAHazardAssessment, item) for item in g6["assessments"])
        for field in ("cancer_surveillance_requirements", "candidate_benefit_evidence_ids"):
            g6[field] = tuple(g6[field])
        safety = evaluate_mva_g6(**g6, decision_policy=policy)
        g7 = dict(facts["g7"])
        g7["plan"] = decode(ValidationPlan, g7["plan"])
        if g7["plan"].active_ingredient_id != candidate_id:
            raise Track2SynthesisError()
        structural = review_structural_g7(**g7, decision_policy=policy)
        gates = (*first, exposure, safety.gate_result, structural.gate_result)
        if any(key in affected for key in dependencies[candidate_id]):
            gates = tuple(replace(gate, coverage_state="incomplete") for gate in gates)
            structural = replace(structural, gate_result=gates[7])
            safety = replace(safety, gate_result=gates[6])
        if all(key in limited for key in dependencies[candidate_id]) and gates[4].outcome == "pass":
            # Recall confined to unavailable/truncated evidence cannot establish
            # complete independent action/rescue support or a primary tier.
            gates = (*gates[:4], replace(gates[4], outcome="not_evaluable",
                      reason_codes=("insufficient_directness",), coverage_state="incomplete"), *gates[5:])
        decision = decide_candidate_tier(gates, policy)
        review = row.research_submission_review
        eligibility = assess_research_eligibility(candidate_id, gates, review, policy)
        decisions[candidate_id], gates_by_id[candidate_id], reviews[candidate_id] = decision, gates, review
        validation[candidate_id] = structural
        context = replace(row.context,
            child_applicability=population.child_applicability.state,
            general_mva_applicability=population.general_mva_applicability.state,
            exposure=ReportAssessment(exposure.outcome, "; ".join(exposure.reason_codes), exposure.evidence_ids),
            effect_magnitude=(g5["comparison"].effect_magnitude_state if g5["comparison"] else "magnitude_unknown"),
            safety=ReportAssessment(safety.subtype_safety_conclusion, "; ".join(safety.gate_result.reason_codes), safety.gate_result.evidence_ids),
            decision_reasons=decision.decision_reasons or ("pass",),
        )
        # Discovery membership resolves to an actual target/objective and a
        # declared route dependency, rather than a separate hand-written list.
        if any((rationale.objective_id, rationale.route) not in dependencies[candidate_id]
               for rationale in row.candidate.rationales):
            raise Track2SynthesisError()
        affected_routes = {route for objective, route in dependencies[candidate_id]
                           if (objective, route) in affected}
        candidate_coverage = tuple(item for item in row.route_coverage if item.route not in affected_routes)
        candidate_coverage += tuple(item for key in dependencies[candidate_id] if key in affected
                                    for item in coverage_by_dependency[key])
        candidates.append(replace(row, context=context, route_coverage=candidate_coverage,
            exposure_comparison=g5["comparison"], safety_review=safety,
            population_applicability=population, gate_results=gates,
            tier_decision=decision, research_eligibility=eligibility, exclusion=None,
            evidence_dimensions={**row.evidence_dimensions,
                "exact_subtype": context.subtype_applicability,
                "translational_maturity": context.translational_maturity,
                "mechanism_readiness": readiness[origin.scenario_id].readiness_state,
                "pediatric_practicality": population.child_applicability.state,
                "exposure_feasibility": gates[5].outcome, "safety": gates[6].outcome},
        ))
    objective_scenarios = {row.objective.objective_id: row.context.origin.scenario_id
                          for row in raw_tables.target_records}
    if set(objective_scenarios) != set(objectives):
        raise Track2SynthesisError('Assessment objectives require source-bound target scenarios')
    assessment = assess_plan(readiness, objective_scenarios, portfolios,plan['assessment_scope'])
    objective_assessments = {row['objective_id']: row for row in assessment.objectives}
    outcomes = tuple(build_objective_candidate_outcome(
        objective_id=objective, decisions={key: decisions[key] for key in ids},
        gate_results={key: gates_by_id[key] for key in ids},
        evidence_gaps={key: tuple(plan["candidate_inputs"][key]["evidence_gaps"]) for key in ids},
        resolution_experiments={key: tuple(plan["candidate_inputs"][key]["resolution_experiments"]) for key in ids},
        other_submission_requirements_met=False, research_reviews={key: reviews[key] for key in ids},
        decision_policy=policy, assessment=objective_assessments[objective],
    ) for objective, ids in sorted(objectives.items()))
    exclusions = {item.active_ingredient_id: item for outcome in outcomes for item in outcome.exclusions}
    tables = build_report_tables(replace(raw_tables,
        candidate_records=tuple(replace(row, exclusion=exclusions.get(row.candidate.active_ingredient_id))
                                for row in sorted(candidates, key=lambda item: item.candidate.active_ingredient_id))),
        decision_policies={policy["policy_version"]: policy},allow_unqualified=diagnosis,
    )
    candidate_decisions = {d.candidate_id: asdict(d) for outcome in outcomes for d in outcome.research_decisions}
    target_origins={row.target.target_id:row.context.origin.scenario_id for row in raw_tables.target_records}
    profiles = tuple(build_target_profile(targets[key], profile_inputs[key],
        genotype_effect=readiness[target_origins[key]].direction, candidate_decisions=candidate_decisions)
        for key in sorted(targets))
    claims = build_report_claim_bundle(tuple(decode(ReportClaim, c) for c in plan["claims"]), report_tables=tables)
    if plan['classification']=='public' and any(c.publication_classification != "public" for c in claims.claims):
        raise Track2SynthesisError("public report contains a nonpublic claim")
    return tables, outcomes, milestones, profiles, claims, assessment, assessments


def policies():
    policy = load_track2_policies(ROOT / "config/track2/source-domain-matrix-v1.json",
                                 ROOT / "config/track2/decision-policy-g0-g7-v1.json").decision_policy
    expansion = load_public_expansion_policy(ROOT / "config/track2/public-expansion-policy-v1.json",
        ROOT / "config/track2/source-domain-matrix-v1.json")
    return policy, expansion


def validate_research_report(output):
    """Readiness consumers must verify the immutable report and required figures."""
    from genosapiens.track2_report_metadata import checked_report_metadata
    return checked_report_metadata(output).readiness


def validate_study_report(output, *, presentation=True):
    """Replay the current study report's retained scientific evidence."""
    root = Path(output)
    manifest = json.loads((root/'manifest.json').read_bytes())
    if manifest.get('schema') != 'genosapiens.track2-research-manifest/v1':
        raise Track2SynthesisError('The study report contract is obsolete or unsupported.')
    if 'report.md' not in manifest['outputs']:
        raise Track2SynthesisError('Report lacks its historical presentation.')
    detail_name = 'scientific-details.md' if 'presentation.json' in manifest['outputs'] else 'report.md'
    validate_study_science(root, manifest, detail_name=detail_name)
    if presentation and 'presentation.json' in manifest['outputs']:
        from genosapiens.track2_report_checks import check_paper
        check_paper(root)
    return manifest['readiness']


def validate_study_science(root, manifest, *, detail_name='scientific-details.md', verify_outputs=True):
    """Replay retained study decisions without requiring a report presentation."""
    root = Path(root)
    from genosapiens.track2_report_captions import scientific_text
    def scientific_prose(name):
        return scientific_text(_source_file(root,name).read_text())
    from genosapiens.governance import _file_identity
    exclusion_files={'filter-exclusions/records.jsonl','filter-exclusions/manifest.json'}
    if exclusion_files & set(manifest['outputs']):
        from genosapiens.track2_exclusions import check_exclusions
        if not exclusion_files <= set(manifest['outputs']):
            raise Track2SynthesisError('Report exclusion artifact coverage is incomplete.')
        descriptor_raw=_source_file(root,'filter-exclusions/manifest.json').read_bytes()
        descriptor=json.loads(descriptor_raw)
        if (sha256(descriptor_raw)!=manifest['outputs']['filter-exclusions/manifest.json']
                or descriptor['sha256']!=manifest['outputs']['filter-exclusions/records.jsonl']):
            raise Track2SynthesisError('Report exclusion artifact binding differs.')
        check_exclusions(root/'filter-exclusions',descriptor)
    for name, digest in (manifest['outputs'].items() if verify_outputs else ()):
        if name in exclusion_files:
            continue
        if _file_identity(_source_file(root,name))[0] != digest:
            raise Track2SynthesisError('report output is missing or changed')
    if (root/'lincs').exists():
        from genosapiens.track2_lincs_report import validate_figures
        if 'lincs/results.json' not in manifest['outputs']:
            raise Track2SynthesisError('LINCS outputs missing from report manifest')
        validate_figures(root,json.loads((root/'lincs/results.json').read_bytes()))
    if json.loads((root/'readiness.json').read_bytes()) != manifest['readiness']:
        raise Track2SynthesisError('report readiness disagrees with manifest')
    required = {'assessment.json', 'assessment-inputs.json', 'plan-source.json', 'outcomes.json', 'coverage.json', detail_name}
    if not required <= set(manifest['outputs']):
        raise Track2SynthesisError('Report lacks retained assessment completion evidence')
    retained = json.loads(_source_file(root, 'assessment-inputs.json').read_bytes())
    plan = retained['plan']
    if plan['classification'] != manifest['classification']:
        raise Track2SynthesisError('The study report classification differs from its retained plan.')
    if plan['input_origin']=='diagnosis_guided_vcf':
        from genosapiens.track2_exclusions import validate_descriptor
        if not (exclusion_files|{'input-package-source.json','selection-evidence.json','source-evidence.json'}) <= set(manifest['outputs']):
            raise Track2SynthesisError('Report lacks complete exclusion evidence.')
        package_raw=_source_file(root,'input-package-source.json').read_bytes()
        package=json.loads(package_raw)
        selection_raw=_source_file(root,'selection-evidence.json').read_bytes()
        selection=json.loads(selection_raw)
        if (sha256(package_raw)!=manifest['bindings']['input']
                or sha256(selection_raw)!=package['artifacts']['selection']['sha256']
                or selection['coverage']['filter_exclusions']!=descriptor):
            raise Track2SynthesisError('Report selection differs from the bound input package.')
        validate_descriptor(descriptor,selection=selection)
        from genosapiens.track2_deferred import check_inventory, DIRECTORY as DEFERRED_DIRECTORY
        deferred_files={DEFERRED_DIRECTORY+'/records.jsonl',DEFERRED_DIRECTORY+'/manifest.json'}
        if not deferred_files <= set(manifest['outputs']):
            raise Track2SynthesisError('Report lacks complete deferred selection evidence.')
        check_inventory(root/DEFERRED_DIRECTORY,selection)
        from genosapiens.track2_selection_report import render_selection
        for name in (detail_name,'input-dossier.md'):
            if render_selection(selection,deferred_path='deferred-alleles/records.jsonl') not in scientific_prose(name):
                raise Track2SynthesisError('Report selection limits differ from the bound input package.')
        source_raw=_source_file(root,'source-evidence.json').read_bytes()
        if sha256(source_raw)!=package['artifacts']['evidence']['sha256']:
            raise Track2SynthesisError('Report source evidence differs from its package.')
        from genosapiens.track2_vcf_evidence import render as render_vcf_evidence, validate_calls
        vcf_parts={'selection':selection,'evidence':json.loads(source_raw)}
        validate_calls(vcf_parts)
        if '\n'.join(render_vcf_evidence(vcf_parts)) not in scientific_prose('input-dossier.md'):
            raise Track2SynthesisError('Dossier VCF tables differ from source evidence.')
    if (sha256(_source_file(root, 'plan-source.json').read_bytes()) != manifest['bindings']['plan'] or
        json.loads(_source_file(root, 'plan-source.json').read_bytes()) != plan):
        raise Track2SynthesisError('Retained assessment plan differs from its bound source')
    if plan['frozen_assessment'] is not None:
        from genosapiens.track2_scientific_planning import validate_packet
        config=plan['frozen_assessment']
        if ('frozen-assessment.json' not in manifest['outputs'] or
            manifest['outputs']['frozen-assessment.json'] != config['sha256'] or
            manifest['bindings'].get('frozen_assessment') != config['sha256']):
            raise Track2SynthesisError('Report lacks its bound frozen assessment evidence')
        packet=json.loads(_source_file(root,'frozen-assessment.json').read_bytes())
        _,source_coverage=validate_packet(packet,config['public_evidence_sha256'])
        from genosapiens.track2_scientific_planning import validate_frozen_routes
        validate_frozen_routes(plan,source_coverage)
        if json.loads(_source_file(root,'coverage.json').read_bytes())['routes'] != plan['routes']:
            raise Track2SynthesisError('Report routes differ from the bound frozen assessment')
    if retained['plan_review'] is not None:
        if ('plan-review-source.json' not in manifest['outputs'] or
            sha256(_source_file(root, 'plan-review-source.json').read_bytes()) != manifest['bindings']['plan_review'] or
            json.loads(_source_file(root, 'plan-review-source.json').read_bytes()) != retained['plan_review']):
            raise Track2SynthesisError('Retained assessment review differs from its bound source')
    mechanism_states = {sid: assess_mechanism_readiness(facts) for sid, facts in plan['mechanism_inputs'].items()}
    coverage = json.loads(_source_file(root, 'coverage.json').read_bytes())
    if coverage['scope']['content_sha256'] != plan['scope_sha256']:
        raise Track2SynthesisError('Assessment scope differs from its reviewed plan')
    policy, _ = policies()
    if sha256(canonical(policy)) != manifest['bindings']['policy']:
        raise Track2SynthesisError('Assessment policy differs from the current decision contract')
    portfolios = {oid: build_route_portfolio([oid],
        [decode(RouteAssessment, r) for r in coverage['routes'] if r['objective_id'] == oid], policy,
        public_scope=coverage['scope']) for oid in plan['objective_candidates']}
    objective_scenarios = {r['objective']['objective_id']:r['context']['origin']['scenario_id']
                          for r in plan['tables']['target_records']}
    if set(objective_scenarios) != set(portfolios):
        raise Track2SynthesisError('Assessment objectives lack bound scenario evidence')
    if plan['assessment_scope'] is not None:
        from genosapiens.track2_mechanism_review import validate_scope
        validate_scope(plan['mechanism_assessment'],plan['assessment_scope'])
    assessment = assess_plan(mechanism_states, objective_scenarios, portfolios,plan['assessment_scope'])
    actual = json.loads(_source_file(root, 'assessment.json').read_bytes())
    if canonical(assessment.as_dict()) != canonical(actual):
        raise Track2SynthesisError('Retained assessment status disagrees with scenario or route evidence')
    review = retained['plan_review']
    reviewed = bool(review and review['decision'] == 'accepted' and
                    review['reviewer_type'] in {'human', 'synthetic_fixture'})
    if 'automatic_plan' in manifest['bindings']:
        raw=_source_file(root,'plan-acceptance.json').read_bytes();automatic=json.loads(raw)
        if (sha256(raw)!=manifest['bindings']['automatic_plan'] or automatic!=plan_acceptance(
                plan,manifest['bindings']['input'],timestamp=automatic['timestamp'],image=manifest['runtime']['image'])):
            raise Track2SynthesisError('Retained automatic plan decision is stale or changed.')
        reviewed=reviewed or bool(automatic['outcome']=='accepted' and (review is None or review['decision']!='rejected'))
        if automatic_summary(plan,automatic) not in scientific_prose(detail_name):
            raise Track2SynthesisError('Automatic decision wording differs from the retained evidence.')
    if 'automatic_input' in manifest['bindings']:
        raw=_source_file(root,'input-acceptance.json').read_bytes()
        accepted_input=json.loads(raw)
        if (sha256(raw)!=manifest['bindings']['automatic_input'] or accepted_input['authority']!='machine'
                or accepted_input['subject_sha256']!=manifest['bindings']['input']):
            raise Track2SynthesisError('Retained automatic input decision changed.')
    readiness = manifest['readiness']
    if (readiness['assessment_status'] != assessment.status or
        readiness['research_report_ready'] != (assessment.complete and readiness['input_report_ready'] and reviewed)):
        raise Track2SynthesisError('Research readiness disagrees with retained assessment evidence')
    outcomes = json.loads(_source_file(root, 'outcomes.json').read_bytes())
    by_objective = {o['objective_id']:o for o in assessment.objectives}
    if len(outcomes) != len(by_objective) or {o['objective_id'] for o in outcomes} != set(by_objective):
        raise Track2SynthesisError('Objective outcomes omit required assessment evidence')
    for outcome in outcomes:
        evidence = by_objective[outcome['objective_id']]
        expected_status = ('research_hypothesis_identified' if outcome['research_candidate_ids'] else
                           'no_plausible_candidate_identified' if evidence['status'] == 'complete' else
                           'assessment_' + evidence['status'])
        if (outcome['assessment_status'] != evidence['status'] or
            outcome['assessment_complete'] != (evidence['status'] == 'complete') or
            outcome['assessment_reasons'] != evidence['reasons'] or outcome['track2_objective_status'] != expected_status):
            raise Track2SynthesisError('Objective outcome disagrees with retained assessment evidence')
    if readiness['candidate_objective_met'] != any(o['research_candidate_ids'] for o in outcomes):
        raise Track2SynthesisError('Candidate objective disagrees with retained decisions')
    expected_summary = assessment_summary(assessment)
    if expected_summary not in scientific_prose(detail_name):
        raise Track2SynthesisError('Report wording disagrees with retained assessment evidence')
    if 'analysis_origins' in manifest:
        from genosapiens.track2_structural import validate_report_section
        name = 'structural-evidence.json'
        if name not in manifest['outputs'] or manifest['outputs'][name] != manifest['bindings']['structural']:
            raise Track2SynthesisError('report structural evidence binding differs')
        structural = json.loads(_source_file(root, name).read_bytes())
        for report in ('input-dossier.md', detail_name):
            validate_report_section(scientific_prose(report), structural)
    return manifest['readiness']


def assessment_summary(assessment):
    text = {'not_performed': 'Medication assessment was not performed. No completed negative result is available.',
            'incomplete': 'Medication assessment is incomplete. Required work remains unresolved.',
            'complete': 'Medication assessment is complete within the declared scope.'}[assessment.status]
    reasons = '\n'.join('- ' + r.code + ': ' + ' / '.join(v for v in
                        (r.scenario_id, r.objective_id, r.route) if v) for r in assessment.reasons)
    excluded='\n'.join('- '+s['scenario_id']+': '+s['rationale'] for s in assessment.scenarios if not s['required'])
    return ('## Medication assessment\n\n' + text
        + ('\n\nReasons (exact identifiers):\n\n' + reasons if reasons else '')
        + ('\n\nProspectively excluded scenarios (exact identifiers):\n\n'+excluded if excluded else ''))


def automatic_summary(plan, decision):
    if decision is None: return ''
    lines=['## Automatic acceptance','',
        'Authority: machine. Policy: '+decision['policy']+'.',
        'Plan outcome: '+decision['outcome']+'.',
        'The automatic decision does not record human review or submission authorization.',
        'Prediction acceptance permits exploratory research. It does not establish curated function or experimental rescue.','',
        '| Scenario | Mechanism acceptance | Exploratory applicability | Required assessment |',
        '|---|---|---|---|']
    scope=plan['assessment_scope'] or {}
    for row in plan['mechanism_assessment']['scenarios']:
        acceptance=row['mechanism_inputs'].get('prediction_acceptance')
        state=acceptance['decision']['outcome'] if acceptance else 'actual_review_path'
        lines.append('| '+' | '.join((row['scenario_id'],state,
            row.get('applicability',{}).get('outcome','actual_review_path'),
            str(scope.get(row['scenario_id'],{}).get('required',True)).lower()))+' |')
    functional_rows=[]
    for row in plan['mechanism_assessment']['scenarios']:
        acceptance=row['mechanism_inputs'].get('prediction_acceptance')
        proof=(acceptance or {}).get('evidence',{}).get('scientific_proof')
        if proof is None: continue
        checks=list(json.loads(proof).get('functional',{}).values())
        if not checks: continue
        native=sorted({check['native']['state'] for check in checks})
        imported=sorted({check['reviewed_import']['state'] for check in checks})
        limits={gap for check in checks for source in check['native'].get('checks',[])
            for gap in source.get('gaps',[])}
        limits.update(source['coverage'] for check in checks for source in check['native'].get('sources',[])
            if source.get('coverage')=='gene_not_in_snapshot')
        if any(check.get('optional_import_gap') for check in checks): limits.add('optional_study_mapping_unavailable')
        functional_rows.append('| '+' | '.join((row['scenario_id'],', '.join(native),', '.join(imported),
            ', '.join(sorted(limits)) or 'No recorded lookup gap'))+' |')
    if functional_rows:
        lines+=['','| Scenario | Functional catalog | Public study import | Source limits |',
            '|---|---|---|---|',*functional_rows]
    lines+=['','Functional evidence applies only within the checked frozen sources.',
        'A complete lookup without contrary records found no contrary evidence within those sources.',
        'Unavailable or invalid required sources leave acceptance unresolved.',
        'Optional study imports can lack mappings. The assessment preserves these gaps separately.',
        'Protein annotations apply only to the matched protein alteration and source experimental context.',
        'They do not establish nucleotide-specific splicing or whole-genotype effects.',
        'Exploratory reuse transfers no variant assay, patient observation, or source eligibility to the new genotype.',
        'Each exploratory objective retains its source tissue, stage, assumptions, hazards, and alternatives.',
        'Unknown residual product cannot support direct activation.',
        'The report identifies the assessed subset. Excluded alternatives do not receive medication assessment.']
    return '\n'.join(lines)


def freeze(bundle_root, plan_path, scope_path, runtime, reference_root, frozen_on, *, draft=False):
    bundle, plan, review, bindings = load_plan(bundle_root, plan_path, require_review=not draft)
    if plan['input_origin']!='published_patient_genotype':
        raise Track2SynthesisError('diagnosis-guided synthesis requires an independently frozen public envelope')
    policy, expansion = policies()
    scope = json.loads(Path(scope_path).read_bytes())
    validate_public_research_scope(scope, decision_policy=policy, expansion_policy=expansion)
    if plan["scope_sha256"] != scope["content_sha256"]:
        raise Track2SynthesisError("plan review does not bind declared public scope")
    if ('lincs_cell_specific_expression' in scope['selected_routes']) != (plan['lincs'] is not None):
        raise Track2SynthesisError('selected LINCS route requires explicit source and query bindings')
    bindings["policy"] = sha256(canonical(policy))
    registry = build_public_seed_registry(bundle)
    publication = freeze_public_envelope(registry, source_manifest_sha256s=bindings, runtime=runtime,
        coverage_records=(), expansion_policy=expansion, frozen_on=frozen_on,
        reference_root=reference_root, public_scope=scope)
    return {"path": str(publication.path), "sha256": publication.sha256,
            "envelope_id": publication.envelope_id, "draft": draft}


def freeze_input_snapshot(bundle_root,plan_path,envelope_path,output,*,runtime,draft=False,allow_synthetic=False,
                          controlled_store=None,discovery_reference_root='/discovery-resources'):
    """Private input pin separate from the independently frozen public envelope."""
    bundle,plan,review,bindings=load_plan(bundle_root,plan_path,require_review=not draft,allow_synthetic=allow_synthetic,
        discovery_reference_root=discovery_reference_root)
    if plan['input_origin']!='diagnosis_guided_vcf':
        raise Track2SynthesisError('private snapshot is for the diagnosis-guided package origin')
    policy,expansion=policies()
    bindings['policy']=sha256(canonical(policy))
    if plan['schema']==NATIVE_SCHEMA:
        from genosapiens.track2_discovery_sources import native_snapshot
        descriptor=json.loads(_source_file(Path(plan_path).parent,plan['exploratory_discovery']['resources']['path']).read_bytes())
        snapshot=native_snapshot(plan,bindings,runtime,descriptor)
    else:
        frozen=load_frozen_public_envelope(envelope_path)
        validate_public_research_scope(frozen['public_scope'],decision_policy=policy,expansion_policy=expansion)
        if plan['scope_sha256']!=frozen['public_scope']['content_sha256']:
            raise Track2SynthesisError('plan and independently frozen public scope differ')
        if (set(runtime)!={'image','code_revision'} or not runtime['image'].startswith('genosapiens/')
            or '@sha256:' not in runtime['image'] or not runtime['code_revision']):
            raise Track2SynthesisError('immutable synthesis runtime binding is required')
        snapshot={'schema':'genosapiens.track2-private-input-snapshot/v1','classification':plan['classification'],
            'bindings':bindings,'envelope_sha256':sha256(Path(envelope_path).read_bytes()),'runtime':runtime}
    raw=canonical(snapshot)
    from genosapiens.track2_storage import output_storage
    with output_storage(output, classification=plan['classification'], controlled_store=controlled_store,
                        kind='track2-input-snapshot') as target:
        target.parent.mkdir(mode=0o700,parents=True,exist_ok=True)
        with target.open('xb') as handle:handle.write(raw)
        target.chmod(0o600)
    return {'path':str(target),'sha256':sha256(raw)}


def synthesize(bundle_root, plan_path, envelope_path, output, *, draft=False, allow_synthetic=False,
               input_snapshot_path=None,expected_snapshot_sha256=None,controlled_store=None,
               public_evidence_path=None,discovery_path=None,discovery_reference_root='/discovery-resources',
               discovery_figures='/discovery-figures', retained_science=None,_evidence_only=False):
    bundle, plan, review, bindings = load_plan(bundle_root, plan_path, require_review=not draft,
        allow_synthetic=allow_synthetic,discovery_reference_root=discovery_reference_root,retained_science=retained_science)
    if plan['schema']==NATIVE_SCHEMA:
        from genosapiens.track2_discovery_report import synthesize_native
        return synthesize_native(bundle,plan,review,bindings,plan_path,output,
            input_snapshot_path=input_snapshot_path,expected_snapshot_sha256=expected_snapshot_sha256,
            controlled_store=controlled_store,discovery_path=discovery_path,
            reference_root=discovery_reference_root,figures_root=discovery_figures,_evidence_only=_evidence_only)
    frozen = load_frozen_public_envelope(envelope_path)
    policy, expansion = policies()
    bindings["policy"] = sha256(canonical(policy))
    diagnosis=plan['input_origin']=='diagnosis_guided_vcf'
    matched_public = None
    if diagnosis:
        from genosapiens.track2_public_matching import (retained_matching, load_public_evidence,
            verify_matching, render_matching, copy_matched_report)
        matched_public = retained_matching(bundle.components)
        if matched_public is not None:
            if public_evidence_path is None:
                raise Track2SynthesisError('Public matching requires its installed source evidence.')
            public = load_public_evidence(public_evidence_path,allow_synthetic=allow_synthetic)
            verify_matching(bundle.components, public)
    runtime=frozen['runtime']
    if diagnosis:
        if input_snapshot_path is None or expected_snapshot_sha256 is None:
            raise Track2SynthesisError('diagnosis synthesis requires its private frozen input snapshot')
        raw=Path(input_snapshot_path).read_bytes();snapshot=json.loads(raw)
        if (sha256(raw)!=expected_snapshot_sha256 or snapshot.get('schema')!='genosapiens.track2-private-input-snapshot/v1'
            or snapshot['classification']!=plan['classification'] or snapshot['bindings']!=bindings
            or snapshot['envelope_sha256']!=sha256(Path(envelope_path).read_bytes())):
            raise Track2SynthesisError('private frozen input, source, review, code, schema or envelope changed')
        runtime=snapshot['runtime']
    elif bindings != frozen["manifest_checksums"]:
        raise Track2SynthesisError("frozen input, source, review, code, schema or policy changed")
    validate_public_research_scope(frozen["public_scope"], decision_policy=policy, expansion_policy=expansion)
    if plan["scope_sha256"] != frozen["public_scope"]["content_sha256"]:
        raise Track2SynthesisError("frozen public scope differs from reviewed plan")
    lincs_result = None
    retained_lincs_figures = None
    if plan['frozen_assessment'] is not None:
        from genosapiens.track2_scientific_planning import frozen_bindings
        packet,_=frozen_bindings(Path(plan_path).parent,plan,bundle)
        if 'report/lincs/results.json' in packet['files']:
            retained_result=json.loads(packet['files']['report/lincs/results.json'])
            from genosapiens.track2_lincs_integration import category_eligibility
            category_eligibility(retained_result)
            lincs_result={k:retained_result[k] for k in ('schema','queries','cells','cross_cell',
                'identity_exclusions','source_sha256','terms','redistribution')}
            lincs_result['case_objective_bindings']=plan['frozen_assessment']['scenario_bindings']
    if plan['lincs'] is not None:
        _scenario_readiness(plan,bundle,policy,draft=draft)
        if diagnosis and not bundle.input_report_ready:
            raise Track2SynthesisError('unreviewed diagnosis reports require the input-only plan')
        from genosapiens.track2_lincs import ROUTE
        from genosapiens.track2_lincs_integration import execute, route_attempt
        from genosapiens.track2_target_hypothesis import TargetHypothesis
        targets = [decode(TargetHypothesis, r['target']) for r in plan['tables']['target_records']]
        requests = [r for r in plan['requests'] if r['request']['route'] == ROUTE]
        routes = [r for r in plan['routes'] if r['route'] == ROUTE]
        if len(targets) != 1 or len(requests) != 1 or len(routes) != 1:
            raise Track2SynthesisError('bounded LINCS target/route must be uniquely selected')
        scenario=plan['tables']['target_records'][0]['context']['origin']['scenario_id']
        if not assess_mechanism_readiness(plan['mechanism_inputs'][scenario]).permits_planning:
            raise Track2SynthesisError('direction-dependent LINCS queries require a supported mechanism')
        if not run_frozen_synthesis(decode(SynthesisRequest,requests[0]['request']),frozen,lambda _: True).authorized:
            raise Track2SynthesisError('LINCS execution outside frozen public scope')
        from genosapiens.track2_report_science_reuse import retained_study_science
        retained=retained_study_science(retained_science,ROOT,plan,bundle.input_sha256)
        if retained is None:
            lincs_result = execute(bundle_root,plan['lincs'],target=targets[0],scenario_id=scenario)
        else:
            lincs_result,retained_lincs_figures=retained[:2]
        plan = {**plan,'routes':[route_attempt(lincs_result,r) if r['route']==ROUTE else r for r in plan['routes']]}
    tables, outcomes, milestones, profiles, claims, assessment, route_records = _compose(plan, bundle, frozen, policy,draft=draft)
    if lincs_result is not None:
        from genosapiens.track2_lincs_integration import discovery_union
        raw = decode(ReportTableInputs, plan['tables'])
        from genosapiens.track2_no_candidate import include_lincs_research_membership
        from genosapiens.track2_report_summary import include_lincs_review
        if 'case_objective_bindings' in lincs_result:
            from genosapiens.track2_report_tables import TABLE_NAMES
            scenario_decisions=[]
            for binding in lincs_result['case_objective_bindings']:
                sid=binding['scenario_id']
                subset=replace(tables,**{name:tuple(r for r in getattr(tables,name) if r.origin.scenario_id==sid)
                                       for name in TABLE_NAMES})
                joined={**lincs_result,**discovery_union(lincs_result,
                    [r.candidate for r in raw.candidate_records if r.scenario_id==sid],subset.candidates)}
                for decision in joined['candidate_decisions']:
                    decision['objective_ids']=[binding['objectives'][oid] for oid in decision['objective_ids']]
                    scenario_decisions.append({'scenario_id':sid,'decision':decision})
                outcomes=include_lincs_research_membership(outcomes,joined)
                subset=include_lincs_review(subset,joined)
                tables=replace(tables,**{name:tuple(r for r in getattr(tables,name) if r.origin.scenario_id!=sid)
                    +getattr(subset,name) for name in TABLE_NAMES})
            # Figures summarize one frozen experiment. Per-scenario scientific
            # decisions remain separate and never replace another scenario's gates.
            lincs_result.update(discovery_union(lincs_result,[],[]))
            for decision in lincs_result['candidate_decisions']:
                decision['objective_ids']=sorted({b['objectives'][oid] for b in lincs_result['case_objective_bindings']
                                                 for oid in decision['objective_ids']})
            lincs_result['scenario_candidate_decisions']=scenario_decisions
        else:
            lincs_result.update(discovery_union(lincs_result,[r.candidate for r in raw.candidate_records],tables.candidates))
            outcomes=include_lincs_research_membership(outcomes,lincs_result)
            tables=include_lincs_review(tables,lincs_result)
    accepted = review is not None and review["decision"] == "accepted" and (review["reviewer_type"] == "human" or
        allow_synthetic and review['reviewer_type']=='synthetic_fixture' and
        bundle.data.get('evidence_origin')=='independently_synthetic')
    automatic=load_plan_acceptance(plan_path,json.loads(Path(plan_path).read_bytes()),bundle) if diagnosis else None
    accepted=accepted or bool(automatic and automatic['outcome']=='accepted'
                             and (review is None or review['decision']!='rejected'))
    readiness = {"input_report_ready": bundle.input_report_ready,
        "research_report_ready": bundle.input_report_ready and accepted and assessment.complete,
        "assessment_status": assessment.status,
        "candidate_objective_met": any(o.research_candidate_ids for o in outcomes),
        "case_matching": "not_assessed", "submission_ready": False}
    methods = plan["methods"]
    if (set(methods) != {"description", "llm_disclosure", "effort", "limitations", "validation_followup"}
            or any(not value for value in methods.values())):
        raise Track2SynthesisError("methods, disclosure, actual effort and limitations are required")
    from genosapiens.track2_report_summary import target_and_candidates, route_coverage, review_summary
    report = ["# Diagnosis-guided repurposing research report" if diagnosis else
              "# TRIP13 repurposing research report" if len(profiles) == 1 and profiles[0]['target_id'] == 'TRIP13'
              else "# Published-genotype repurposing research report", "",
        ("This report evaluates the frozen input within the diagnosis scope. It retains alternative and unresolved mechanisms. "
         "Case linkage and submission readiness require separate assessment." if diagnosis else
         "This report searches for medication leads for a published genetic condition. "
         "It does not identify an established treatment or match a challenge patient."), "",
        "Research report ready: **" + str(readiness["research_report_ready"]).lower() + "**. "
        "Candidate objective met: **" + str(readiness["candidate_objective_met"]).lower() + "**. "
        "Submission ready: **false**.", "",
        assessment_summary(assessment), "",
        automatic_summary(plan,automatic), "",
        target_and_candidates(profiles, tables, lincs_result, assessment_status=assessment.status), "",
        route_coverage(plan, profiles, lincs_result), ""]
    if diagnosis:
        if not bundle.input_report_ready or not accepted:
            report[2:2] = ['**Draft report. Input or plan acceptance is incomplete.**',
                'Report generation does not record human acceptance. Review the report before use.','']
        from genosapiens.track2_phenotype import render_phenotype_context
        if bundle.data['evidence_origin']=='independently_synthetic':
            report += ['**This report uses independently synthetic evidence. It is not a patient analysis or a human scientific review.**','']
        report += [render_phenotype_context(bundle.components['phenotype'])]
        from genosapiens.track2_structural import render_structural
        report += [render_structural(bundle.components['structural'])]
        from genosapiens.track2_selection_report import render_selection
        report += [render_selection(bundle.components['selection'],deferred_path='deferred-alleles/records.jsonl')]
        report += [render_matching(bundle.components, report_link=True)]
        if not profiles:
            report += ['## Mechanism and candidate gaps','',
                'No qualified intervention objective is available from the retained mechanisms. '
                'The pipeline did not start a medication query that requires a known direction. The candidate objective remains unmet.','']
    if lincs_result is not None:
        report += [review_summary(tables, lincs_result), ""]
    report += ["## Next steps", "", methods["validation_followup"], "",
               "## Review", "",
               (f"Decision: {review['decision'].replace('_', ' ')}. Reviewer: {review['reviewer']}.\n\nReview text (verbatim):\n\n{review['rationale']}"
                if review is not None else 'This report contains no human plan review.'), "",
               "The plan review and final report review are separate. Final report reviews identify the exact files by their hashes."]
    supplements = ["", "## Supplementary material", "",
        "Detailed methods, source records and complete evidence tables follow. "
        "[Input dossier](input-dossier.md) · [Structured evidence tables](tables.json) · [Standalone tables](tables.md).", "",
        "### Methods and disclosure", "", methods["description"], "", methods["llm_disclosure"], "",
        "Effort: " + methods["effort"], "", "Limitations: " + methods["limitations"], "",
        "### Scientific interpretation and evidence links", ""]
    for claim in claims.claims:
        supplements += [f"- [{claim.claim_kind}] {claim.text if claim.text is not None else '`' + json.dumps(claim.value, sort_keys=True) + '`'}",
                        "  Evidence: " + ", ".join(claim.evidence_ids) + "."]
    supplements += ["", render_target_profiles(profiles).replace("### ", "#### ").replace("## Intervention", "### Intervention"),
                    "", "### Assessed scientific outcomes", ""]
    for outcome in outcomes:
        supplements += [f"- {outcome.objective_id}: {outcome.track2_objective_status}. Scientific outcome: {outcome.scientific_outcome}. "
                        + "Research candidates: " + (", ".join(outcome.research_candidate_ids) or "none") + "."]
    supplements += ["", "<details>", "<summary>Full evidence tables: candidates, exclusions and G0–G7</summary>", "",
                    report_tables_as_markdown(tables, heading_level=3), "", "</details>", "",
                    "### Source catalog", "", "<details>", "<summary>Show all source citations and provenance</summary>", ""]
    for e in tables.evidence_catalog:
        supplements += [f"- {e.evidence_id}: {e.citation}", f"  Source details: {e.provenance}"]
    supplements += ["", "</details>", ""]
    from genosapiens.track2_report_language import report_terms, report_acknowledgement
    supplements += [report_terms(lincs=lincs_result is not None, diagnosis=diagnosis), "",
                    report_acknowledgement(), ""]
    from genosapiens.track2_diagnosis_package import render_diagnosis_dossier
    contents = {
        "input-dossier.md": (render_diagnosis_dossier(bundle,exclusion_path='filter-exclusions/records.jsonl',
            deferred_path='deferred-alleles/records.jsonl')
            if diagnosis else render_input_dossier(bundle)).encode(),
        "report.md": "\n".join(report).encode(),
        "tables.json": report_tables_as_json(tables).encode(),
        "tables.md": report_tables_as_markdown(tables).encode(),
        "target-profiles.json": canonical(profiles),
        "claims.json": canonical(claims.as_dict()),
        "outcomes.json": canonical([asdict(o) for o in outcomes]),
        "assessment.json": canonical(assessment.as_dict()),
        "assessment-inputs.json": canonical({'plan': json.loads(Path(plan_path).read_bytes()), 'plan_review':review}),
        "plan-source.json": Path(plan_path).read_bytes(),
        "coverage.json": canonical({"scope": frozen["public_scope"], "routes":
                                   [asdict(r) for r in route_records],
                                   "milestones": {k: asdict(v) for k, v in milestones.items()}}),
        "methods.json": canonical(methods), "readiness.json": canonical(readiness),
    }
    if review is not None:
        contents['plan-review-source.json'] = Path(plan_path).with_name('plan-review.json').read_bytes()
    if automatic is not None:
        contents['plan-acceptance.json']=canonical(automatic)
    if diagnosis:
        contents['input-acceptance.json']=canonical(bundle.automatic_acceptance)
        contents['input-package-source.json']=canonical(bundle.data)
        contents['selection-evidence.json']=canonical(bundle.components['selection'])
        contents['source-evidence.json']=canonical(bundle.components['evidence'])
    if plan['frozen_assessment'] is not None:
        contents['frozen-assessment.json']=_source_file(Path(plan_path).parent,plan['frozen_assessment']['path']).read_bytes()
    if diagnosis:
        contents['phenotype.json']=canonical(bundle.components['phenotype'])
        contents['structural-evidence.json']=canonical(bundle.components['structural'])
        if matched_public is not None:
            contents['public-evidence-matching.json']=canonical(matched_public)
        supported=[key for key,facts in plan['mechanism_inputs'].items()
                   if assess_mechanism_readiness(facts).permits_planning]
        contents['mechanism-outcomes.json']=canonical({'mechanism_count':len(plan['mechanism_inputs']),
            'supported_mechanism_ids':sorted(supported),'qualified_objective_count':len(plan['objective_candidates']),
            'curated_ready_ids':sorted(k for k,v in plan['mechanism_inputs'].items()
                                      if assess_mechanism_readiness(v).readiness_state=='ready'),
            'prediction_accepted_ids':sorted(k for k,v in plan['mechanism_inputs'].items()
                                      if assess_mechanism_readiness(v).readiness_state=='prediction_accepted'),
            'status':'no_qualified_mechanism' if not supported else 'no_qualified_intervention_objective'
                if not plan['objective_candidates'] else 'intervention_objectives_assessed',
            'candidate_objective_met':readiness['candidate_objective_met'],
            'scope_limitations':bundle.components['scope']['limitations']})
    manifest = {"schema": "genosapiens.track2-research-manifest/v1", "classification": plan['classification'],
                **({'analysis_origins':[asdict(o) for o in report_origins(bundle,draft=draft,assessment=plan['mechanism_assessment'])],
                    'input_snapshot_sha256':expected_snapshot_sha256} if diagnosis else
                   {'analysis_origin':asdict(bundle.admitted_origin())}), "runtime": runtime,
                "bindings": bindings, "frozen_envelope_sha256": sha256(Path(envelope_path).read_bytes()),
                "outputs": {name: sha256(content) for name, content in sorted(contents.items())},
                "readiness": readiness}
    contents["manifest.json"] = canonical(manifest)
    from genosapiens.track2_storage import output_storage
    with output_storage(output, classification=plan['classification'], controlled_store=controlled_store,
                        kind='track2-research-report') as target:
        target.mkdir(parents=True, exist_ok=False, mode=0o700)
        for name, content in contents.items():
            if name not in {'readiness.json','manifest.json'} and not (_evidence_only and name=='report.md'):
                (target / name).write_bytes(content)
        if diagnosis:
            from genosapiens.track2_exclusions import copy_exclusions, DIRECTORY, ROWS, MANIFEST
            descriptor=bundle.components['selection']['coverage']['filter_exclusions']
            copy_exclusions(bundle.root/DIRECTORY,target/DIRECTORY,descriptor)
            manifest['outputs'].update({DIRECTORY+'/'+ROWS:descriptor['sha256'],
                DIRECTORY+'/'+MANIFEST:sha256(canonical(descriptor))})
            from genosapiens.track2_deferred import copy_inventory, DIRECTORY as DEFERRED_DIRECTORY
            descriptor=bundle.components['selection']['coverage']['deferred_inventory']
            copy_inventory(bundle.root/DEFERRED_DIRECTORY,target/DEFERRED_DIRECTORY,bundle.components['selection'])
            manifest['outputs'].update({DEFERRED_DIRECTORY+'/'+ROWS:descriptor['sha256'],
                DEFERRED_DIRECTORY+'/'+MANIFEST:sha256(canonical(descriptor))})
        if matched_public is not None and matched_public['matched_public_result'] is not None:
            manifest['outputs'].update(copy_matched_report(public_evidence_path,target))
        if lincs_result is not None:
            from genosapiens.track2_lincs_report import render_lincs, validate_figures
            section, assets = render_lincs(lincs_result, target,
                retained_figures_root=Path(public_evidence_path)/'report' if plan['frozen_assessment'] is not None else retained_lincs_figures)
            validate_figures(target, lincs_result)
            report += [section]
            manifest['outputs'].update({name:sha256((target/name).read_bytes()) for name in assets})
        (target/'scientific-details.md').write_text('\n'.join(report) + '\n'.join(supplements))
        from genosapiens.track2_report_adapters import study_inputs
        from genosapiens.track2_report_paper import compose_paper
        from genosapiens.track2_report_gene_metadata import load_gene_metadata
        from genosapiens.track2_report_context import retained_gene_ids
        objectives,nominations,coverage,cells=study_inputs(bundle,plan,profiles,route_records,lincs_result,tables)
        if _evidence_only:
            from genosapiens.track2_report_paper import export_scientific_supplements
            from genosapiens.track2_report_evidence import publish_study_evidence
            model=export_scientific_supplements(target,bundle=bundle,plan=plan,objectives=objectives,
                nominations=nominations,coverage=coverage,selected_cells=cells,
                gene_metadata=load_gene_metadata(ROOT,retained_gene_ids(bundle,plan,objectives)))
            (target/'readiness.json').write_bytes(contents['readiness.json'])
            return publish_study_evidence(target,model=model,scientific_metadata=manifest)
        contents['report.md'] = compose_paper(target,bundle=bundle,plan=plan,objectives=objectives,
            nominations=nominations,coverage=coverage,selected_cells=cells,readiness=readiness,
            title=report[0],gene_metadata=load_gene_metadata(ROOT,retained_gene_ids(bundle, plan, objectives)))
        (target/'report.md').write_bytes(contents['report.md'])
        from genosapiens.track2_scientific_storage import file_identity as report_file_identity
        presentation_files=['scientific-details.md','presentation.json'] + [str(p.relative_to(target)) for p in (target/'supplements').iterdir()]
        manifest['outputs'].update({name:report_file_identity(target/name)['sha256'] for name in presentation_files})
        for name in manifest['outputs']:
            if (name.endswith('.md') and (target/name).parent in {target,target/'lincs'}) or name=='lincs/manifest.json':
                manifest['outputs'][name]=report_file_identity(target/name)['sha256']
        manifest['outputs']['report.md'] = sha256(contents['report.md'])
        contents['manifest.json'] = canonical(manifest)
        (target/'readiness.json').write_bytes(contents['readiness.json'])
        (target/'manifest.json').write_bytes(contents['manifest.json'])
    return {"output": str(target), "manifest_sha256": sha256(contents["manifest.json"]), **readiness}


def export_study_evidence(*args, **kwargs):
    """Produce retained study evidence without writing the historical paper."""
    if '_evidence_only' in kwargs:
        raise Track2SynthesisError('The evidence export mode is fixed.')
    return synthesize(*args,**kwargs,_evidence_only=True)
