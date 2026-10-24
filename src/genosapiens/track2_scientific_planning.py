"""Prepare case plans from reviewed mechanisms and authenticated frozen studies.

The source study supplies scientific inputs. A separate local review establishes
case mechanism and objective applicability. Matching alone grants neither.
"""
import copy
from dataclasses import asdict, fields, replace
import json
from pathlib import Path

from genosapiens.track2_codec import decode, json_value
from genosapiens.track2_diagnosis_plan import input_only_plan
from genosapiens.track2_endpoint import EndpointSelection, select_endpoint
from genosapiens.track2_mechanism import assess_mechanism_readiness, build_mechanism_model
from genosapiens.track2_operator_inputs import require, small_json, write_new
from genosapiens.track2_published_genotype import canonical, sha256, _source_file
from genosapiens.track2_report_tables import (ReportTableInputs, ReportEvidenceReference,
    EndpointTableRecord, TargetHypothesisTableRecord)
from genosapiens.track2_target_hypothesis import derive_intervention_objective, build_target_hypothesis

PACKET_SCHEMA = 'genosapiens.track2-frozen-assessment-evidence/v1'


def source_tables(plan):
    """Project the added optional readiness field without changing source bytes."""
    tables=copy.deepcopy(plan['tables'])
    for row in tables['mechanism_records']:
        row['readiness'].setdefault('prediction_acceptance',None)
    return decode(ReportTableInputs,tables)


def exploratory_applicability(bundle, scenario, packet, public_sha256):
    """Evaluate reuse as a research hypothesis in the declared source context."""
    source,_=validate_packet(packet,public_sha256)
    mapping=json.loads(packet['files']['mapping.json']) if 'mapping.json' in packet['files'] else None
    readiness=assess_mechanism_readiness(scenario['mechanism_inputs'])
    tables=source_tables(source)
    endpoints={r.endpoint.endpoint_id:r.endpoint for r in tables.endpoint_records}
    directions={r['direction'] for r in source['mechanism_inputs'].values()}
    same_model=(mapping is not None and (mapping['gene_id'],mapping['disease_id'])==
                (scenario['gene_id'],scenario['disease_id']))
    rows=[]
    for record in tables.target_records:
        target,objective=record.target,record.objective
        endpoint=endpoints[objective.endpoint_id]
        values={f.name:getattr(target,f.name) for f in fields(target)
                if f.name not in {'objective_id','ready_for_drug_retrieval'}}
        build_target_hypothesis(objective=objective,**values)
        residual=[a['residual_product'] for a in (scenario['model'] or {}).get('allele_models',[])]
        conditions={
            'signed_endpoint':endpoint.relevance=='supported' and endpoint.expected_direction in {'increase','decrease'},
            'research_endpoint':endpoint.endpoint_class not in {'fixed_developmental_injury','cancer_risk_outcome'},
            'tissue_stage':target.tissue==endpoint.tissue and target.stage==endpoint.stage
                and all(v.lower() not in {'unknown','unresolved','not_assessed'} for v in (target.tissue,target.stage)),
            'residual_product':not (target.causal_role=='direct' and target.required_modulation in {'increase','partial_increase'})
                or bool(residual) and all(v=='present' for v in residual),
            'assumptions_hazards_alternatives':bool(objective.assumptions and objective.harm_conditions
                and target.harm_conditions and target.alternatives),
            'no_target_conflicts':not target.conflicts}
        rows.append({'objective_id':objective.objective_id,'target_id':target.target_id,
            'conditions':conditions,'tissue':target.tissue,'stage':target.stage,
            'assumptions':list(objective.assumptions),'harm_conditions':list(target.harm_conditions),
            'alternatives':list(target.alternatives),'conflicts':list(target.conflicts),
            'endpoint_id':endpoint.endpoint_id,'source_reference_ids':sorted(set(target.evidence_ids+objective.evidence_ids))})
    reasons=[]
    if not readiness.permits_planning: reasons.append('mechanism_not_accepted')
    if not same_model: reasons.append('gene_or_specific_disease_not_supported')
    if readiness.direction not in directions: reasons.append('source_mechanism_direction_differs')
    if not rows: reasons.append('no_source_objectives')
    reasons.extend(r['objective_id']+':'+k for r in rows for k,v in r['conditions'].items() if not v)
    from genosapiens.track2_public_matching import retained_matching
    matching=retained_matching(bundle.components)
    exact=bool(matching and any(d['scenario_id']==scenario['construction_scenario_id'] and d['status']=='matched'
                              for d in matching['decisions']))
    return {'schema':'genosapiens.track2-exploratory-applicability/v1','authority':'machine',
        'construction_sha256':bundle.input_sha256,'scenario_id':scenario['scenario_id'],
        'public_evidence_sha256':public_sha256,'exact_genotype_match':exact,
        'outcome':'accepted' if not reasons else 'unresolved','reasons':reasons or ['compatible_frozen_research_context'],
        'records':rows,'case_matched':False,'variant_specific_evidence_transferred':False,
        'limitations':['Research applies to the stated source tissue and stage.',
            'The new genotype has no transferred assay result or observed treatment response.',
            'Unknown residual product cannot support direct activation.']}


def frozen_packet(public_root):
    root=Path(public_root)
    names=['evidence.json','plan.json','plan-review.json','report/manifest.json','report/coverage.json']
    for name in ('input/genotype.json','input/review.json','mapping.json','mapping-review.json'):
        if (root/name).is_file(): names.append(name)
    if (root/'report/lincs/results.json').is_file():
        names.append('report/lincs/results.json')
    return {'schema':PACKET_SCHEMA,'files':{name:_source_file(root,name).read_text() for name in names}}


def validate_packet(packet, expected):
    require(set(packet)=={'schema','files'} and packet['schema']==PACKET_SCHEMA,
            'Frozen assessment evidence has an unsupported contract.')
    files=packet['files']; export=json.loads(files['evidence.json'])
    require(sha256(files['evidence.json'].encode())==expected,'Frozen assessment belongs to different public evidence.')
    for name,raw in files.items():
        if name!='evidence.json':
            require(export['files'].get(name)==sha256(raw.encode()),'Frozen assessment source bytes changed.')
    plan=json.loads(files['plan.json']); review=json.loads(files['plan-review.json'])
    manifest=json.loads(files['report/manifest.json'])
    require(review['decision']=='accepted' and review['reviewer_type']=='human'
        and review['plan_sha256']==sha256(files['plan.json'].encode())==manifest['bindings']['plan']
        and manifest['bindings']['plan_review']==sha256(files['plan-review.json'].encode()),
        'Frozen assessment lacks its actual source review.')
    for name in ('coverage.json','lincs/results.json'):
        if 'report/'+name in files:
            require(manifest['outputs'].get(name)==sha256(files['report/'+name].encode()),
                    'Frozen assessment report evidence differs from its source manifest.')
    if 'report/lincs/results.json' in files:
        result=json.loads(files['report/lincs/results.json'])
        config=plan['lincs']
        require(config is not None and config['sources_sha256']==result['source_sha256']
            and manifest['bindings']['lincs:sources']==result['source_sha256']
            and set(config['queries'])==set(result['queries'])
            and all(config['queries'][cell]['sha256']==sha256(canonical(query))==
                    manifest['bindings']['lincs:query:'+cell] for cell,query in result['queries'].items()),
            'Frozen LINCS evidence is outside its reviewed source plan.')
    return plan,json.loads(files['report/coverage.json'])


def validate_frozen_routes(plan, coverage):
    """Frozen completion is valid only for the exact retained route evidence."""
    expected=[]
    objectives=set()
    for binding in plan['frozen_assessment']['scenario_bindings']:
        objectives.update(binding['objectives'].values())
        expected.extend({**row,'objective_id':binding['objectives'][row['objective_id']]}
                        for row in coverage['routes'])
    key=lambda r:(r['objective_id'],r['route'])
    require(objectives==set(plan['objective_candidates'])
        and canonical(sorted(expected,key=key))==canonical(sorted(plan['routes'],key=key)),
        'Frozen routes differ from retained source evidence. New searches require a new public assessment.')


def frozen_bindings(root, plan, bundle):
    """Verify frozen reuse against immutable matching and new case review evidence."""
    config=plan['frozen_assessment']
    if config is None:
        return None,{}
    require(set(config)=={'path','sha256','public_evidence_sha256','scenario_bindings'},
            'Frozen assessment binding is incomplete.')
    raw=_source_file(Path(root),config['path']).read_bytes()
    require(sha256(raw)==config['sha256'],'Frozen assessment evidence changed.')
    packet=json.loads(raw)
    source,coverage=validate_packet(packet,config['public_evidence_sha256'])
    from genosapiens.track2_public_matching import retained_matching
    matching=retained_matching(bundle.components)
    require(matching is not None and matching['public_evidence_sha256']==config['public_evidence_sha256'],
            'Frozen assessment has no source binding.')
    matched={r['scenario_id'] for r in matching['decisions'] if r['status']=='matched'}
    scenarios={r['scenario_id']:r for r in plan['mechanism_assessment']['scenarios']}
    source_directions={r['direction'] for r in source['mechanism_inputs'].values()}
    require(config['scenario_bindings'],'Frozen assessment has no reviewed case objectives.')
    for row in config['scenario_bindings']:
        scenario=scenarios[row['scenario_id']]
        if 'applicability' in scenario:
            expected=exploratory_applicability(bundle,scenario,packet,config['public_evidence_sha256'])
            require(scenario['applicability']==expected and expected['outcome']=='accepted'
                and row['review_evidence_id']==scenario['applicability_evidence_id']
                and row['applicability_path']=='exploratory', 'Exploratory applicability differs from its sources.')
            require(set(row['objectives'])==set(source['objective_candidates'])
                and set(row['objectives'].values())<=set(plan['objective_candidates']),
                'Exploratory objective mapping differs from its source.')
            continue
        evidence={e['evidence_id']:e for e in plan['mechanism_assessment']['evidence']}
        require(scenario['construction_scenario_id'] in matched
            and scenario['effect_subject_id'] is not None
            and scenario['mechanism_inputs']['direction'] in source_directions
            and assess_mechanism_readiness(scenario['mechanism_inputs']).readiness_state=='ready'
            and row['review_evidence_id'] in evidence,
            'Frozen assessment requires separate supported case mechanism and objective review.')
        endorsement=json.loads(evidence[row['review_evidence_id']]['assertion'])
        require(endorsement['rows'] and all(r['objective_applicability']=='supported'
            and r['decision']=='accepted' and r['construction_sha256']==bundle.input_sha256
            and r['public_evidence_sha256']==config['public_evidence_sha256']
            and r['construction_scenario_id']==scenario['construction_scenario_id'] for r in endorsement['rows']),
            'Frozen objectives lack their actual source-bound case review.')
        require(set(row['objectives'])==set(source['objective_candidates'])
            and set(row['objectives'].values()) <= set(plan['objective_candidates']),
            'Frozen assessment objective mapping differs from the reviewed plan.')
    validate_frozen_routes(plan,coverage)
    return packet,{'frozen_assessment':config['sha256']}


def prepare_scientific_plan(bundle, *, assessment, public, public_root, frozen, endorsements,
                            methods, output, draft=False, assessment_scope=None, discovery_resources=None):
    """Generate plans with existing endpoint, objective, target and route owners."""
    if discovery_resources is not None:
        from genosapiens.track2_exploratory_planning import prepare_native_plan
        return prepare_native_plan(bundle,assessment=assessment,resources=discovery_resources,
            methods=methods,output=output,draft=draft)
    plan=input_only_plan(bundle,assessment=assessment,public_scope_sha256=frozen['public_scope']['content_sha256'],
                         methods=methods,draft=draft)
    plan=json_value(plan)
    if assessment_scope is not None:
        from genosapiens.track2_mechanism_review import validate_scope
        validate_scope(assessment,assessment_scope)
        plan['assessment_scope']=assessment_scope
    plan['frozen_assessment']=None
    gaps=[]
    ready={s['scenario_id']:s for s in assessment['scenarios']
           if assess_mechanism_readiness(s['mechanism_inputs']).permits_planning}
    if not ready or public['summary'] is None:
        return plan,[{'scenario_id':s['scenario_id'],'reason':'mechanism_unresolved' if s['scenario_id'] not in ready
                     else 'no_frozen_objective_evidence'} for s in assessment['scenarios']]
    from genosapiens.track2_public_matching import matching_decisions
    matches=matching_decisions(bundle.components,public)
    matched={r['scenario_id'] for r in matches['decisions'] if r['status']=='matched'}
    packet=frozen_packet(public_root)
    source,coverage=validate_packet(packet,public['sha256'])
    require(source['scope_sha256']==plan['scope_sha256'],'Frozen therapeutic scope differs from the installed scope.')
    source_records=source_tables(source)
    source_direction={r['direction'] for r in source['mechanism_inputs'].values()}
    tables=decode(ReportTableInputs,plan['tables'])
    contexts={r.context.origin.scenario_id:r.context for r in tables.mechanism_records}
    endpoints=[];targets=[];candidates=[];profiles={};bindings=[]
    for scenario in assessment['scenarios']:
        sid=scenario['scenario_id']; review=endorsements.get(sid)
        exploratory='applicability' in scenario
        if exploratory and scenario['applicability']['outcome']=='accepted':
            review={'decision':'accepted','objective_applicability':'supported',
                    'review_evidence_id':scenario['applicability_evidence_id']}
        if assessment_scope is not None and not assessment_scope[sid]['required']:
            gaps.append({'scenario_id':sid,'reason':'prospectively_excluded','rationale':assessment_scope[sid]['rationale']})
            continue
        if (sid not in ready or not exploratory and scenario['construction_scenario_id'] not in matched
            or scenario['mechanism_inputs']['direction'] not in source_direction or review is None):
            gaps.append({'scenario_id':sid,'reason':'mechanism_unresolved' if sid not in ready else
                'no_exact_frozen_objective_evidence' if scenario['construction_scenario_id'] not in matched else
                'objective_applicability_review_missing'})
            continue
        require(review['decision']=='accepted' and review['objective_applicability']=='supported',
                'Scientific planning requires actual objective applicability review.')
        context=contexts[sid]
        model=build_mechanism_model(scenario['model'])
        suffix=sha256(sid.encode())[:16]
        endpoint_map={r.endpoint.endpoint_id:r.endpoint.endpoint_id+':'+suffix for r in source_records.endpoint_records}
        objective_map={oid:oid+':'+suffix for oid in source['objective_candidates']}
        generated={}
        for row in source_records.endpoint_records:
            values={f.name:getattr(row.endpoint,f.name) for f in fields(EndpointSelection)
                    if f.name not in {'relevance','modifiability'}}
            endpoint=select_endpoint(**{**values,'endpoint_id':endpoint_map[row.endpoint.endpoint_id]})
            generated[row.endpoint.endpoint_id]=endpoint
            endpoints.append(EndpointTableRecord(context if exploratory else replace(row.context,origin=context.origin),sid,endpoint))
        for row in source_records.target_records:
            original=row.objective
            objective=derive_intervention_objective(objective_id=objective_map[original.objective_id],
                mechanism=model,endpoint=generated[original.endpoint_id],
                desired_change=original.desired_change,compartment=original.compartment,
                assumptions=original.assumptions,harm_conditions=original.harm_conditions,
                mechanism_review_id=(scenario['applicability_evidence_id'] if exploratory else scenario['effect_subject_id']),
                endpoint_review_id=review['review_evidence_id'],
                evidence_ids=(*original.evidence_ids,review['review_evidence_id']))
            values={f.name:getattr(row.target,f.name) for f in fields(row.target)
                    if f.name not in {'objective_id','ready_for_drug_retrieval'}}
            values['signed_rescue_path']=tuple(replace(e,to_id=endpoint_map.get(e.to_id,e.to_id),
                from_id=endpoint_map.get(e.from_id,e.from_id)) for e in row.target.signed_rescue_path)
            target=build_target_hypothesis(objective=objective,**values)
            targets.append(TargetHypothesisTableRecord(context if exploratory else replace(row.context,origin=context.origin),objective,target))
        for row in source_records.candidate_records:
            from genosapiens.track2_candidate_union import build_candidate_union
            candidate=build_candidate_union(tuple(replace(r,objective_id=objective_map[r.objective_id])
                for r in row.candidate.rationales)).candidates[0]
            # Source acceptance never becomes a case candidate review by matching.
            candidate_context=(replace(row.context,origin=context.origin,translational_maturity='computational_inference',
                effect_magnitude='magnitude_unknown',decision_reasons=('exploratory_source_context_only',))
                if exploratory else replace(row.context,origin=context.origin))
            candidates.append(replace(row,context=candidate_context,scenario_id=sid,
                selected_endpoint_id=endpoint_map[row.selected_endpoint_id],candidate=candidate,
                research_submission_review=None,research_eligibility=None))
        for row in source['profile_inputs']:
            if exploratory:
                row={**row,'residual_functional_product':'unknown',
                    'residual_product_context':'The accepted prediction does not establish residual functional product in this genotype.',
                    'residual_product_evidence_ids':sorted(s['applicability_evidence_id'] for s in assessment['scenarios']
                        if s.get('applicability',{}).get('outcome')=='accepted'),
                    'separately_supported_restoration':False,'restoration_evidence_ids':[]}
            require(profiles.setdefault(row['target_id'],row)==row,'Frozen target profiles disagree.')
        for key,value in source['candidate_inputs'].items():
            require(plan['candidate_inputs'].setdefault(key,value)==value,'Scenario candidate evidence disagrees.')
        for old,new in objective_map.items():
            plan['objective_candidates'][new]=source['objective_candidates'][old]
        for key,dependencies in source['candidate_route_dependencies'].items():
            plan['candidate_route_dependencies'].setdefault(key,[]).extend(
                [[objective_map[oid],route] for oid,route in dependencies])
        plan['routes'].extend({**row,'objective_id':objective_map[row['objective_id']]} for row in coverage['routes'])
        plan['requests'].extend({**row,'objective_id':objective_map[row['objective_id']]} for row in source['requests'])
        bindings.append({'scenario_id':sid,'objectives':objective_map,'review_evidence_id':review['review_evidence_id'],
            **({'applicability_path':'exploratory'} if exploratory else {})})
    if not bindings:
        return plan,gaps
    raw=canonical(packet);path='frozen-assessment.json'
    write_new(Path(output)/path,raw)
    plan['frozen_assessment']={'path':path,'sha256':sha256(raw),'public_evidence_sha256':public['sha256'],
                               'scenario_bindings':bindings}
    # A reviewed frozen assessment is a retained dataset record. It preserves the
    # original review, raw plan bytes, route evidence, and report output bindings.
    from genosapiens.track2_manual_import import make_manual_evidence_import, DocumentLocator, ExtractedField
    source_review=json.loads(packet['files']['plan-review.json'])
    record=make_manual_evidence_import(document_type='dataset_record',document_identity=public['sha256'],
        document_version='frozen-reviewed-assessment',document_checksum=sha256(raw),
        document_url=small_json(Path(public_root)/'input/genotype.json')['sources'][0]['public_url'],
        accessed_on=source_review['reviewed_at'][:10],locator=DocumentLocator(None,'Original reviewed plan and report evidence'),
        extracted_fields=(ExtractedField('assessment','frozen source records','See retained original bytes',None,
            'Case applicability requires separate review','Source limitations remain applicable'),),
        reviewer_role='scientific_reviewer',review_record_id=source_review['plan_sha256'],
        review_time=source_review['reviewed_at'],structured_interface_status='reviewed_manual_extraction',
        overall_transcription_uncertainty='Retained original source bytes',
        overall_interpretation_uncertainty='Case applicability requires separate review')
    plan['evidence_imports']=[{'record':asdict(record),'path':path,
        'evidence_ids':[e.evidence_id for e in source_records.evidence_catalog],
        'license_or_terms':'Original source terms remain in the frozen public export',
        'redistribution':'Private local evidence reuse only',
        'query_ids':sorted({q for r in coverage['routes'] if r['attempt'] and r['attempt']['executed']
                            for q in r['attempt']['query_ids']})}]
    tables=replace(tables,evidence_catalog=(*tables.evidence_catalog,*source_records.evidence_catalog),
        endpoint_records=tuple(endpoints),target_records=tuple(targets),candidate_records=tuple(candidates))
    plan.update(tables=asdict(tables),profile_inputs=list(profiles.values()))
    from genosapiens.track2_report_claims import make_report_claim
    plan['claims']=[asdict(make_report_claim(source_identity=bundle.input_sha256,role='scope_summary',
        claim_kind='computed_result',evidence_ids=tuple(e.evidence_id for e in tables.evidence_catalog),
        publication_classification=bundle.data['classification'],value={
            'mechanism_scenarios':len(assessment['scenarios']),
            'qualified_intervention_objectives':len(plan['objective_candidates']),
            'frozen_source_sha256':public['sha256'], 'planning_gaps':gaps}))]
    return json_value(plan),gaps
