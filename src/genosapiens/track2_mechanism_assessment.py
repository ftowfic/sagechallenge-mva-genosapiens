"""Immutable downstream effects over current-format construction evidence."""
import copy
from dataclasses import asdict, dataclass, replace
import json

from genosapiens.evidence_model import _typed, AlleleEvidence
from genosapiens.evidence_serialization import _plain
from genosapiens.genotype_hypotheses import ConstructedHypothesis
from genosapiens.track2_adjudication import AdjudicationRecord, automatic_decision
from genosapiens.track2_codec import decode, json_value
from genosapiens.track2_diagnosis_builder import MODEL_SCHEMA, scenario_model
from genosapiens.track2_mechanism import INPUT_FIELDS, build_mechanism_model
from genosapiens.track2_published_genotype import canonical, sha256
from genosapiens.track2_scoped_interpretation import require
from genosapiens.track2_scoped_mechanisms import ReviewedEffect, assemble_scoped_mechanisms, allele_prediction, derived_evidence

SCHEMA = 'genosapiens.track2-mechanism-assessment/v2'
NATIVE_SCHEMA = 'genosapiens.track2-mechanism-assessment/v3'


def automatic_effect(bundle, scenario, inputs, *, fact_cache=None, public_packet_sha256=None):
    """Replay automatic research acceptance from retained scientific producers."""
    from genosapiens.model_requirements import ModelRequirementPolicy, TranscriptScope
    from genosapiens.track2_functional_sources import functional_check
    from genosapiens.track2_public_matching import retained_matching, reviewed_functional_check
    from genosapiens.track2_vcf_selection import replay_sequence_quality
    if not any(e.startswith('scoped-model:') for e in scenario['evidence_ids']):
        from genosapiens.track2_adjudication import AUTOMATIC_CHECKS
        facts=copy.deepcopy(scenario['mechanism_inputs'])
        checks={k:{'state':'unresolved','reasons':['source_scoped_proof_unavailable']}
                for k in AUTOMATIC_CHECKS['mechanism']}
        details={'checks':checks,'construction_sha256':bundle.input_sha256}
        facts['prediction_acceptance']={'decision':automatic_decision(domain='mechanism',subject=facts,
            evidence=details,checks=checks,execution=inputs['execution'],timestamp=inputs['timestamp']),
            'evidence':details}
        return facts,None,derived_evidence('automatic-mechanism-evidence',details,scenario['genotype_id'])
    constructed,mapping=construction_proof(bundle,scenario)
    context=constructed.hypothesis.context
    selected=[e for e in bundle.components['evidence'] if e['evidence_id'] in scenario['evidence_ids']]
    def records(prefix):
        result=[]
        for row in selected:
            if row['evidence_id'].startswith(prefix+':'):
                require(sha256(row['assertion'].encode())==row['source_identity'],'Automatic source proof changed.')
                result.append(json.loads(row['assertion']))
        return result
    alleles=[_typed(r,AlleleEvidence,decode=True) for r in records('scoped-allele')]
    require({a.allele_id for a in alleles}==set(mapping),'Automatic allele evidence domain differs.')
    policies=records('scoped-functional-sources')
    require(len(policies)<=1,'Automatic functional policy is ambiguous.')
    policy=_typed(policies[0],ModelRequirementPolicy,decode=True) if policies else None
    public=retained_matching(bundle.components)
    packet=inputs['public_packet']
    require(packet is None or public is not None,'Frozen functional import lacks its construction binding.')
    predictions={};qualities={};functional={};per_allele={}
    policy_key = [row['source_identity'] for row in selected
                  if row['evidence_id'].startswith('scoped-functional-sources:')]
    allele_keys = {json.loads(row['assertion'])['allele_id']: row['source_identity']
                  for row in selected if row['evidence_id'].startswith('scoped-allele:')}
    context_key = sha256(canonical(_plain(context)))
    for allele in alleles:
        key=mapping[allele.allele_id]
        def calculate():
            prediction=allele_prediction(allele,context,policy)
            quality=[replay_sequence_quality(call) for call in allele.calls or ()]
            tx=tuple(t['transcript_id'] for t in prediction['transcripts'])
            native=(functional_check(context,allele,TranscriptScope(tx),policy) if policy is not None else
                {'state':'unavailable','sources':[],'matches':[],'contrary':False})
            imported=reviewed_functional_check(packet,public['public_evidence_sha256'] if public else None,context,allele,tx)
            return {'prediction':prediction,'quality':quality,'native':native,'imported':imported}
        if fact_cache is None:
            result = calculate()
        else:
            require(public_packet_sha256 is not None, 'Allele fact reuse requires its public source binding.')
            cache_key = sha256(canonical({'allele':allele_keys[allele.allele_id], 'context':context_key,
                'policy':policy_key, 'public_packet':public_packet_sha256,
                'public_evidence':public['public_evidence_sha256'] if public else None}))
            result = fact_cache.resolve(cache_key, calculate)
        prediction=result['prediction'];qualities[key]=result['quality']
        native=result['native'];imported=result['imported']
        available=([native] if native['state']!='unavailable' else [])+(
            [imported] if imported['state']!='unavailable' else [])
        functional[key]={'native':native,'reviewed_import':imported,
            'complete':bool(available) and all(r['state']=='complete' for r in available),
            'contrary':any(r['contrary'] for r in available),
            'optional_import_gap':packet is not None and imported['state']=='unavailable'}
        quality=qualities[key]
        state=('conflicting' if any(q and q['state']=='excluded' for q in quality) or functional[key]['contrary'] else
            'accepted' if quality and all(q and q['state']=='accepted' for q in quality)
            and prediction['state']=='predicted' and functional[key]['complete'] else 'unresolved')
        per_allele[key]={'outcome':state,'direction':prediction['direction'],
            'strength':'policy_prediction','residual_product':'unknown','transcript_fate':'unknown'}
        predictions[key]=prediction
    reqs=context.requirements
    conflicts=any('track2_specific_requirement_conflict' in r.diagnostics for r in reqs)
    incompatible=any(r.mechanism_direction not in (None,'loss_of_function') for r in reqs)
    specific=bool(reqs) and bool(context.native_disease_id) and bool(context.evidence_refs) and all(
        r.required_allele_count is not None and 'disease_level_inheritance_not_gene_specific' not in r.diagnostics for r in reqs)
    def check(ok,reason,conflict=False):
        return {'state':'conflicting' if conflict else 'accepted' if ok else 'unresolved','reasons':[reason]}
    checks={
        'genotype':check(constructed.model.state=='supported' and constructed.hypothesis.completeness.value=='complete'
            and constructed.hypothesis.phase.value in {'in_trans','not_applicable'},
            'complete_compatible_model_required',constructed.model.state=='contradicted'),
        'quality':check(all(q and all(r and r['state']=='accepted' for r in q) for q in qualities.values()),
            'selected_sample_quality_required',any(r and r['state']=='excluded' for q in qualities.values() for r in q)),
        'inheritance':check(specific and not conflicts,'specific_source_model_required',conflicts),
        'transcripts':check(all(r['state']=='predicted' for r in predictions.values()),'relevant_transcript_predictions_required'),
        'direction':check(all(r['direction']=='loss_of_function' for r in predictions.values()) and not incompatible,
            'compatible_joint_loss_of_function_assumption_required',incompatible),
        'functional_sources':check(all(r['complete'] and not r['contrary'] for r in functional.values()),
            'completed_applicable_functional_lookup_required',any(r['contrary'] for r in functional.values()))}
    accepted=all(c['state']=='accepted' for c in checks.values())
    facts=copy.deepcopy(scenario['mechanism_inputs'])
    model=json.loads(canonical(asdict(scenario_model(bundle.components,scenario))))
    proof={'construction_sha256':bundle.input_sha256,'source_hashes':bundle.source_hashes,
        'genotype':_plain(constructed),'predictions':predictions,'quality':qualities,'functional':functional,
        'allele_acceptance':per_allele,'checks':checks}
    # Native source references have their own namespaces and replay owner.
    # Preserve their original bytes without treating them as report catalog IDs.
    details={'checks':checks,'scientific_proof':canonical(proof).decode('ascii')}
    evidence=derived_evidence('automatic-mechanism-evidence',details,scenario['genotype_id'])
    for item in model['allele_models']:
        key=item['allele_id']; p=predictions[key]
        if per_allele[key]['outcome']=='accepted':
            item.update(transcript_id='&'.join(t['transcript_id'] for t in p['transcripts']),
                molecular_function='loss_of_function',evidence_support_types=['policy_prediction'],
                evidence_ids=[evidence['evidence_id']])
    facts.update(technical_support='supported' if checks['quality']['state']=='accepted' else
        'unsupported' if checks['quality']['state']=='conflicting' else 'uncertain',
        exact_disease_validity='conflicting' if conflicts else 'exact_supported' if specific else 'unknown',
        transcript_relevance='relevant' if checks['transcripts']['state']=='accepted' else 'unknown')
    if accepted:
        facts.update(direction='loss_of_function',joint_evidence='prediction_accepted',inheritance_consistency='consistent')
        model['joint_model'].update(joint_effect='loss_of_function',evidence_support_types=['policy_prediction'],
            evidence_ids=[evidence['evidence_id']])
    else:
        facts['conflicts']=sorted(set(facts['conflicts'])|{k for k,v in checks.items() if v['state']=='conflicting'})
    for item in facts['allele_evidence']:
        p=predictions[item['allele_id']]
        item.update(predicted_consequence='&'.join(sorted({c for t in p['transcripts'] for c in t['consequences']})) or None,
            curated_consequence=None,curation_status='predicted' if p['transcripts'] else 'unknown',
            evidence_ids=[evidence['evidence_id']])
    decision=automatic_decision(domain='mechanism',subject=facts,evidence=details,checks=checks,
        execution=inputs['execution'],timestamp=inputs['timestamp'])
    facts['prediction_acceptance']={'decision':decision,'evidence':details}
    return facts,model,evidence


@dataclass(frozen=True)
class MechanismProjection:
    alleles: tuple


@dataclass(frozen=True)
class GenotypeProjection:
    hypotheses: tuple
    binding: dict


def assess_effect(bundle, scenario, constructed, effect):
    """Project retained allele/genotype records into the existing mechanism owner.

    This projection has its own binding. It never impersonates the original
    interpretation or changes its source hash. No annotation runs again.
    """
    selected=set(scenario['evidence_ids'])
    alleles=[]
    for row in bundle.components['evidence']:
        if row['evidence_id'] in selected and row['evidence_id'].startswith('scoped-allele:'):
            require(sha256(row['assertion'].encode()) == row['source_identity'], 'Observed allele proof changed.')
            alleles.append(_typed(json.loads(row['assertion']), AlleleEvidence, decode=True))
    require({a.allele_id for a in alleles} == set(constructed.hypothesis.allele_ids),
            'Mechanism projection lacks the observed allele domain.')
    if effect.assessments is not None and effect.assessments['technical_support']=='supported':
        from genosapiens.track2_vcf_selection import replay_sequence_quality
        qualities=[replay_sequence_quality(call) for a in alleles for call in a.calls or ()]
        require(qualities and all(q is not None and q['state']=='accepted' for q in qualities),
                'A generic mechanism review cannot waive failed or missing sample quality.')
    projection=MechanismProjection(tuple(alleles))
    genotypes=GenotypeProjection((constructed,), {
        'construction_sha256':bundle.input_sha256,
        'interpretation_sha256':sha256(canonical(_plain(projection)))})
    results=assemble_scoped_mechanisms(genotypes,projection,evidence_origin=bundle.data['evidence_origin'],
                                     reviewed_effects=(effect,))
    require(len(results)==1, 'Reviewed effect has an ambiguous mechanism projection.')
    return results[0]


def construction_proof(bundle, scenario):
    records = [r for r in bundle.components['evidence']
               if r['evidence_id'] in scenario['evidence_ids'] and r['evidence_id'].startswith('scoped-model:')]
    require(len(records) == 1, 'Mechanism assessment requires its construction proof.')
    record = records[0]
    require(sha256(record['assertion'].encode()) == record['source_identity'], 'Construction proof changed.')
    proof = json.loads(record['assertion'])
    require(proof['schema'] == MODEL_SCHEMA and proof['mechanism']['genotype_id'] == scenario['genotype_id'],
            'Mechanism assessment has a foreign genotype.')
    constructed = _typed(proof['genotype'], ConstructedHypothesis, decode=True)
    require(constructed.hypothesis.hypothesis_id == scenario['genotype_id'], 'Constructed genotype identity differs.')
    mapping = proof['normalized_to_package_allele_ids']
    require(set(mapping) == set(constructed.hypothesis.allele_ids)
            and set(mapping.values()) == set(scenario['allele_ids']), 'Constructed allele mapping differs.')
    return constructed, mapping


def effect_from_record(record):
    require(set(record) == {'genotype_id', 'model', 'evidence', 'adjudication', 'assessments'},
            'Reviewed effect fields differ from the current contract.')
    return ReviewedEffect(record['genotype_id'], record['model'], tuple(record['evidence']),
                          decode(AdjudicationRecord, record['adjudication']), record['assessments'])


def _assessment_header(bundle, effects, automatic_inputs):
    require(bundle.data['input_origin'] == 'diagnosis_guided_vcf', 'Mechanism assessment requires constructed input.')
    if automatic_inputs is not None:
        require(set(automatic_inputs) in ({'execution','timestamp','public_packet'},
            {'execution','timestamp','public_packet','discovery'}),'Automatic assessment inputs are incomplete.')
        if 'discovery' in automatic_inputs:
            from genosapiens.track2_exploratory_planning import installed_policy
            binding=automatic_inputs['discovery']
            require(set(binding)=={'resources_sha256','policy_sha256'}
                and binding['policy_sha256']==sha256(canonical(installed_policy(binding['policy_sha256'])))
                and isinstance(binding['resources_sha256'],str) and len(binding['resources_sha256'])==64,
                'Native mechanism assessment requires its discovery policy and resource bindings.')
    native=automatic_inputs is not None and 'discovery' in automatic_inputs
    result = {'schema': NATIVE_SCHEMA if native else SCHEMA, 'classification': bundle.data['classification'],
        'construction_sha256': bundle.input_sha256, 'construction_review_sha256': bundle.review_sha256,
        'construction_components': bundle.source_hashes,
        'effects': sorted(effects, key=lambda r: (r['construction_scenario_id'], sha256(canonical(r['effect']))))}
    if automatic_inputs is not None:
        result['automatic_inputs'] = automatic_inputs
    return json_value(result)


def _assessment_rows(bundle, effects, automatic_inputs):
    """Replay one scenario and its evidence at a time with the scientific owners."""
    source = {r['scenario_id']: r for r in bundle.components['scenarios']}
    # The package loader already authenticates the complete evidence catalog.
    # Each scientific owner below needs this scenario's cited records and the
    # shared public matching record. Index once instead of scanning the catalog
    # for every allele, construction proof, mechanism model, and alternative.
    indexed = {row['evidence_id']: (position, row)
               for position, row in enumerate(bundle.components['evidence'])}
    shared = {key for key in indexed if key.startswith('public-match:')}
    by_source = {sid: [] for sid in source}
    for item in effects:
        require(set(item) == {'construction_scenario_id', 'effect'} and item['construction_scenario_id'] in source,
                'Reviewed effect is outside the construction scope.')
        by_source[item['construction_scenario_id']].append(effect_from_record(item['effect']))
    if automatic_inputs is not None:
        from genosapiens.track2_allele_facts import AlleleFactCache
        from genosapiens.track2_scientific_planning import exploratory_applicability
        from genosapiens.track2_public_matching import retained_matching
        matching = retained_matching(bundle.components)
        fact_cache = AlleleFactCache()
        packet_sha256 = sha256(canonical(automatic_inputs['public_packet']))
    for sid, scenario in source.items():
        selected = sorted((indexed[key] for key in set(scenario['evidence_ids']) | shared if key in indexed),
                          key=lambda item: item[0])
        scoped = replace(bundle, components={**bundle.components, 'evidence': [row for _, row in selected]})
        for effect in by_source[sid] or (None,):
            evidence = []
            application_record = None
            if effect is None:
                assessment_id = sid
                if automatic_inputs is not None:
                    facts,model,automatic_evidence=automatic_effect(scoped,scenario,automatic_inputs,
                        fact_cache=fact_cache,public_packet_sha256=packet_sha256)
                    evidence.append(automatic_evidence)
                else:
                    facts = copy.deepcopy(scenario['mechanism_inputs'])
                    model = (asdict(scenario_model(scoped.components, scenario))
                        if any(e.startswith('scoped-model:') for e in scenario['evidence_ids']) else None)
            else:
                require(bundle.input_report_ready, 'Reviewed construction is required before mechanism assessment.')
                constructed, mapping = construction_proof(scoped, scenario)
                result = assess_effect(scoped, scenario, constructed, effect)
                checked = result.model
                require(not effect.adjudication.blocking_fact, 'Blocking mechanism evidence cannot establish readiness.')
                readiness = result.readiness
                facts = asdict(readiness)
                facts = {k: v for k, v in facts.items() if k in INPUT_FIELDS}
                assessment_id = 'diagnosis-genotype:' + sha256(canonical({
                    'construction_scenario_id': sid, 'effect_sha256': sha256(effect.reviewed_bytes())}))[:24]
                model = asdict(checked)
                for allele in model['allele_models']:
                    allele['allele_id'] = mapping[allele['allele_id']]
                model['joint_model']['allele_ids'] = [mapping[a] for a in model['joint_model']['allele_ids']]
                for allele in facts['allele_evidence']:
                    allele['allele_id'] = mapping[allele['allele_id']]
                evidence.extend(effect.evidence)
            facts['scenario_id'] = assessment_id
            if model is not None:
                model['scenario_id'] = assessment_id
            assessed = {'scenario_id': assessment_id, 'construction_scenario_id': sid,
                'gene_id': scenario['gene_id'], 'genotype_id': scenario['genotype_id'],
                'disease_id': scenario['disease_id'], 'mechanism_inputs': facts, 'model': model,
                'effect_subject_id': effect.subject_id if effect else None,
                'evidence_ids': sorted(set(scenario['evidence_ids']) | (
                    {e['evidence_id'] for e in effect.evidence} if effect else
                    {automatic_evidence['evidence_id']} if automatic_inputs is not None else set()))}
            if automatic_inputs is not None:
                application=(exploratory_applicability(scoped,assessed,automatic_inputs['public_packet'],
                    matching['public_evidence_sha256']) if automatic_inputs['public_packet'] is not None else
                    {'outcome':'unresolved','reasons':['frozen_objective_evidence_unavailable'],
                     'case_matched':False,'variant_specific_evidence_transferred':False})
                row=derived_evidence('exploratory-applicability',application,assessed['scenario_id'])
                application_record = row
                assessed.update(applicability=application,applicability_evidence_id=row['evidence_id'],
                    evidence_ids=sorted([*assessed['evidence_ids'],row['evidence_id']]))
                if 'discovery' in automatic_inputs:
                    from genosapiens.track2_exploratory_planning import discovery_eligibility, installed_policy
                    from genosapiens.track2_mechanism import assess_mechanism_readiness
                    eligibility=discovery_eligibility(assessed,assess_mechanism_readiness(facts),
                        installed_policy(automatic_inputs['discovery']['policy_sha256']))
                    native_record=derived_evidence('discovery-eligibility',eligibility,assessed['scenario_id'])
                    evidence.append(native_record)
                    assessed.update(discovery_eligibility=eligibility,discovery_eligibility_evidence_id=native_record['evidence_id'],
                        evidence_ids=sorted([*assessed['evidence_ids'],native_record['evidence_id']]))
            yield json_value(assessed), json_value(evidence), effect is not None, application_record


def _merge_evidence(evidence, records, reviewed, *, digest=False):
    """Preserve the existing reviewed-conflict and derived-record semantics."""
    for row in records:
        key = row['evidence_id']
        value = sha256(canonical(row)) if digest else row
        if reviewed:
            require(evidence.setdefault(key, value) == value,
                    'Mechanism evidence identity has conflicting content.')
        else:
            evidence[key] = value


def build_assessment(bundle, effects=(), *, automatic_inputs=None, streamed=False):
    """Bind effects to unchanged construction with one copy of each source text."""
    result = _assessment_header(bundle, effects, automatic_inputs)
    scenarios, evidence, applications = {}, {}, {}
    rows=_assessment_rows(bundle,effects,automatic_inputs)
    if streamed:
        from genosapiens.track2_scientific_records import stored_json
        from genosapiens.track2_scientific_storage import Records
        rows=stored_json(Records(rows))
    for scenario, records, reviewed, application in rows:
        sid = scenario['scenario_id']
        require(sid not in scenarios, 'Duplicate mechanism assessment scenario.')
        scenarios[sid] = scenario
        _merge_evidence(evidence, records, reviewed)
        if application is not None:
            applications[application['evidence_id']] = application
    # The original contract inserts applicability records after all mechanisms.
    evidence.update(applications)
    return {**result, 'scenarios': [scenarios[k] for k in sorted(scenarios)],
            'evidence': [evidence[k] for k in sorted(evidence)]}


def validate_assessment(bundle, assessment):
    """Replay every decision without retaining a second complete assessment."""
    require(isinstance(assessment, dict) and assessment.get('schema') in {SCHEMA,NATIVE_SCHEMA},
            'A current downstream mechanism assessment is required.')
    reason = 'Mechanism assessment disagrees with construction, evidence, or actual review.'
    effects, inputs = assessment['effects'], assessment.get('automatic_inputs')
    header = _assessment_header(bundle, effects, inputs)
    require(set(assessment) == set(header) | {'scenarios', 'evidence'}
            and all(assessment[k] == v for k, v in header.items()), reason)
    def index(rows, key):
        result = {row[key]: row for row in rows}
        require([row[key] for row in rows] == sorted(result), reason)
        return result
    scenarios = index(assessment['scenarios'], 'scenario_id')
    evidence = index(assessment['evidence'], 'evidence_id')
    seen_scenarios, expected_evidence, applications = set(), {}, {}
    for scenario, records, reviewed, application in _assessment_rows(bundle, effects, inputs):
        sid = scenario['scenario_id']
        require(sid not in seen_scenarios and scenarios.get(sid) == scenario, reason)
        seen_scenarios.add(sid)
        _merge_evidence(expected_evidence, records, reviewed, digest=True)
        if application is not None:
            applications[application['evidence_id']] = sha256(canonical(application))
    expected_evidence.update(applications)
    require(seen_scenarios == set(scenarios) and set(expected_evidence) == set(evidence), reason)
    require(all(sha256(canonical(evidence[key])) == digest for key, digest in expected_evidence.items()), reason)
    return assessment
