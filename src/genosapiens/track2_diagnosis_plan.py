"""Source-bound input report plan, before reviewed intervention objectives exist."""
from dataclasses import asdict
import copy

from genosapiens.track2_mechanism import assess_mechanism_readiness
from genosapiens.track2_input import assessed_report_origins as report_origins
from genosapiens.track2_report_claims import make_report_claim
from genosapiens.track2_report_tables import (
    ReportEvidenceReference, ReportAssessment, ResolvedReportContext, MechanismTableRecord, ReportTableInputs,
)
from genosapiens.track2_synthesis import SCHEMA, policies


def input_only_plan(bundle, *, public_scope_sha256, methods, draft=False, assessment=None, streamed=False):
    """Construct a reviewable plan; no review, intervention or query is invented."""
    from genosapiens.track2_mechanism_assessment import build_assessment
    assessment = build_assessment(bundle) if assessment is None else assessment
    # report_origins replays the assessment before it admits these scenarios.
    origins={o.scenario_id:o for o in report_origins(bundle,draft=draft,assessment=assessment)}
    policy,_=policies()
    catalog=tuple(ReportEvidenceReference(e['evidence_id'],e['source_identity'],
        'Retained diagnosis-guided evidence',e['lineage_key']) for e in
        (*bundle.components['evidence'], *assessment['evidence']))
    mechanisms=[]
    for row in assessment['scenarios']:
        facts=row['mechanism_inputs'];readiness=assess_mechanism_readiness(facts)
        evidence=tuple(row['evidence_ids'])
        context=ResolvedReportContext(origins[row['scenario_id']],
            'exact_subtype' if readiness.exact_disease_validity=='exact_supported' else 'unresolved',
            'computational_inference','magnitude_unknown','not_assessed','not_assessed',evidence,policy['policy_version'],
            ReportAssessment('not_assessed','No candidate exposure assessed',evidence),
            ReportAssessment('not_assessed','No candidate safety assessed',evidence),('intervention_objectives_not_assessed',))
        mechanisms.append(MechanismTableRecord(context,readiness))
    tables=ReportTableInputs(catalog,tuple(mechanisms),(),(),())
    claims=[asdict(make_report_claim(source_identity=bundle.input_sha256,role='scope_summary',
        claim_kind='computed_result',evidence_ids=tuple(e.evidence_id for e in catalog),
        publication_classification=bundle.data['classification'],
        value={'source_alleles':len(bundle.components['selection']['selected_alleles']),
            'mechanism_alternatives':len(mechanisms),'supported_mechanisms':sum(m.readiness.readiness_state=='ready' for m in mechanisms),
            'prediction_accepted_mechanisms':sum(m.readiness.readiness_state=='prediction_accepted' for m in mechanisms),
            'qualified_intervention_objectives':0,'candidate_objective_met':False}))]
    return {'schema':SCHEMA,'classification':bundle.data['classification'],'input_origin':'diagnosis_guided_vcf',
        'mechanism_assessment':assessment,'assessment_scope':None,
        'tables':asdict(tables),'mechanism_inputs':{row['scenario_id']:
            row['mechanism_inputs'] if streamed else copy.deepcopy(row['mechanism_inputs']) for row in assessment['scenarios']},
        'candidate_inputs':{},'routes':[],'requests':[],'objective_candidates':{},'candidate_route_dependencies':{},
        'profile_inputs':[],'claims':claims,'methods':methods,'evidence_imports':[],
        'scope_sha256':public_scope_sha256,'lincs':None,'frozen_assessment':None}
