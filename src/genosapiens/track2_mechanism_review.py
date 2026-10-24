"""Generated local worksheets for case-specific mechanism and objective review."""
import copy
import csv
from dataclasses import asdict, replace
from datetime import datetime, timezone
import io
import json
import os
from pathlib import Path

from genosapiens.bounded_json import read_bounded_bytes
from genosapiens.track2_adjudication import AdjudicationRecord
from genosapiens.track2_mechanism_assessment import construction_proof, build_assessment
from genosapiens.track2_operator_inputs import require, private_path
from genosapiens.track2_published_genotype import canonical, sha256
from genosapiens.track2_scoped_mechanisms import ReviewedEffect

FIELDS=('construction_sha256','public_evidence_sha256','construction_scenario_id','allele_id',
    'source_direction','alternative_id','allele_effect','joint_effect','technical_support','exact_disease_validity',
    'transcript_relevance','technical_evidence_ids','disease_evidence_ids','transcript_evidence_ids',
    'objective_applicability','scope','scope_rationale','evidence_ids','decision','reviewer',
    'independence','rationale','limitations')


def worksheet(bundle, public):
    if public['summary'] is None:
        return None
    from genosapiens.track2_public_matching import matching_decisions
    matching=matching_decisions(bundle.components,public)
    matched={r['scenario_id'] for r in matching['decisions'] if r['status']=='matched'}
    if not matched or public['summary'] is None:
        return None
    directions={r['direction'] for r in public['summary']['mechanism_inputs'].values()}
    require(len(directions)==1,'Installed source mechanism alternatives require separate reviewed planning evidence.')
    stream=io.StringIO(); writer=csv.DictWriter(stream,fieldnames=FIELDS,delimiter='\t',lineterminator='\n')
    writer.writeheader()
    for scenario in bundle.components['scenarios']:
        constructed,_=construction_proof(bundle,scenario)
        for allele in constructed.hypothesis.allele_ids:
            writer.writerow(dict(construction_sha256=bundle.input_sha256,public_evidence_sha256=public['sha256'],
                construction_scenario_id=scenario['scenario_id'],allele_id=allele,
                source_direction=next(iter(directions)) if scenario['scenario_id'] in matched else 'unknown',
                alternative_id='primary',
                allele_effect='unknown',joint_effect='unknown',technical_support='uncertain',
                exact_disease_validity='unknown',transcript_relevance='unknown',objective_applicability='unassessed',
                scope='required',scope_rationale='All retained scenarios are in the requested assessment scope.',
                decision='pending_human_review'))
    return stream.getvalue().encode()


def read_mechanism_review(path,bundle,public,*,governance,host_uid):
    from genosapiens.participant_governance import authorize_controlled_access
    raw=read_bounded_bytes(private_path(path),maximum_bytes=2*1024*1024)
    reader=csv.DictReader(io.StringIO(raw.decode('utf-8')),delimiter='\t')
    require(tuple(reader.fieldnames or ())==FIELDS,'Use the generated mechanism worksheet columns.')
    rows=list(reader)
    template=worksheet(bundle,public)
    require(template is not None,'No exact frozen evidence supports this mechanism worksheet.')
    expected=list(csv.DictReader(io.StringIO(template.decode()),delimiter='\t'))
    identity=FIELDS[:5]
    expected_identities={tuple(r[k] for k in identity) for r in expected}
    require(rows and {tuple(r[k] for k in identity) for r in rows}==expected_identities
            and all(r.get('alternative_id') for r in rows),
            'Mechanism worksheet source, scenario, or allele bindings changed.')
    by_scenario={s['scenario_id']:s for s in bundle.components['scenarios']}
    known={e['evidence_id'] for e in bundle.components['evidence']} | {
        e['evidence_id'] for e in public['summary']['source_catalog']}
    effects=[]; endorsements={}; submitted=datetime.now(timezone.utc).isoformat()
    for sid,alternative in sorted({(r['construction_scenario_id'],r['alternative_id']) for r in rows}):
        selected=[r for r in rows if (r['construction_scenario_id'],r['alternative_id'])==(sid,alternative)]
        require(len(selected)==len({r['allele_id'] for r in selected})
            and {tuple(r[k] for k in identity) for r in selected}==
                {i for i in expected_identities if i[2]==sid},
            'Each mechanism alternative must cover the complete observed genotype.')
        require(all(set(r)==set(FIELDS) and r['decision'] in {'accepted','rejected','pending_human_review'}
                    for r in selected),'Mechanism worksheet decision is invalid.')
        if any(r['decision']!='accepted' for r in selected):
            continue
        require(bundle.input_report_ready,'Review construction before mechanism evidence.')
        for row in selected:
            require(all(isinstance(row[k],str) and row[k].strip() for k in FIELDS),
                    'Complete every field before accepting a mechanism.')
            authorize_controlled_access(governance,row['reviewer'],host_uid)
            require(row['independence']=='self_review_disclosed',
                    'This single-review worksheet requires self_review_disclosed.')
            refs=row['evidence_ids'].split('|')
            require(refs and set(refs)<=known,'Mechanism review cites evidence outside the retained sources.')
            technical=row['technical_evidence_ids'].split('|')
            require(set(technical)<=set(by_scenario[sid]['evidence_ids'])
                and any(e.startswith('scoped-allele:') for e in technical),
                'Technical review must cite this construction scenario and its observed allele evidence.')
            require(all(set(row[k].split('|'))<=known for k in ('disease_evidence_ids','transcript_evidence_ids')),
                'Disease and transcript review must cite retained evidence.')
        shared=('joint_effect','technical_support','exact_disease_validity','transcript_relevance',
                'objective_applicability','scope','scope_rationale','reviewer','independence','rationale','limitations')
        require(all(len({r[k] for r in selected})==1 for k in shared),
                'Allele rows disagree on their shared genotype assessment.')
        row=selected[0]; scenario=by_scenario[sid]
        proof=next(json.loads(e['assertion']) for e in bundle.components['evidence']
            if e['evidence_id'] in scenario['evidence_ids'] and e['evidence_id'].startswith('scoped-model:'))
        model=copy.deepcopy(proof['mechanism']['model']);model['scenario_id']=scenario['genotype_id']
        assertion=canonical({'rows':selected,'submitted_at':submitted,'worksheet_sha256':sha256(raw)}).decode()
        eid='mechanism-review:'+sha256(assertion.encode())[:24]
        evidence={'evidence_id':eid,'source_identity':sha256(assertion.encode()),
            'lineage_key':public['sha256'],'assertion':assertion,'evidence_role':'curated'}
        values={r['allele_id']:r for r in selected}
        for allele in model['allele_models']:
            allele.update(molecular_function=values[allele['allele_id']]['allele_effect'],
                evidence_support_types=['curated_functional_evidence'],evidence_ids=[eid])
        model['joint_model'].update(joint_effect=row['joint_effect'],
            evidence_support_types=['curated_functional_evidence'],evidence_ids=[eid])
        effect=ReviewedEffect(scenario['genotype_id'],model,(evidence,),None,
            {k:row[k] for k in ('technical_support','exact_disease_validity','transcript_relevance')})
        review=AdjudicationRecord('adjudication:'+sha256(assertion.encode())[:24],effect.subject_id,
            'mechanism','initial','scientific_reviewer',1,row['independence'],'unreviewed','accepted',
            row['rationale'],submitted,(eid,),None,False,sha256(effect.reviewed_bytes()))
        effects.append({'construction_scenario_id':sid,'effect':asdict(replace(effect,adjudication=review))})
        endorsements[effect.subject_id]={'decision':row['decision'],
            'objective_applicability':row['objective_applicability'],'review_evidence_id':eid}
    assessment=build_assessment(bundle,json.loads(canonical(effects)))
    resolved={s['scenario_id']:endorsements[s['effect_subject_id']] for s in assessment['scenarios']
              if s['effect_subject_id'] in endorsements
              and endorsements[s['effect_subject_id']]['objective_applicability']=='supported'}
    scope=scope_from_assessment(assessment)
    return {'assessment':assessment,'endorsements':resolved,'assessment_scope':scope,'worksheet_sha256':sha256(raw),
            'submitted_at':submitted,'all_decisions_accepted':all(r['decision']=='accepted' for r in rows)}


def scope_from_assessment(assessment):
    evidence={e['evidence_id']:e for e in assessment['evidence']}
    result={}
    for scenario in assessment['scenarios']:
        if assessment['schema']=='genosapiens.track2-mechanism-assessment/v3':
            eligibility=scenario['discovery_eligibility']
            result[scenario['scenario_id']]={'required':eligibility['outcome']=='eligible',
                'rationale':'Native prospective scope: '+', '.join(eligibility['reasons']),
                'objective_ids':eligibility['objective_ids']}
            continue
        scope={'required':True,'rationale':'All retained mechanism scenarios are in the assessment scope.'}
        if 'applicability' in scenario:
            application=scenario['applicability']
            scope={'required':application['outcome']=='accepted',
                'rationale':'Automatic prospective scope: '+', '.join(application['reasons'])}
        records=[evidence[e] for e in scenario['evidence_ids'] if e in evidence and e.startswith('mechanism-review:')]
        if records:
            require(len(records)==1,'Scenario scope has ambiguous mechanism review evidence.')
            rows=json.loads(records[0]['assertion'])['rows'];row=rows[0]
            require(row['scope'] in {'required','excluded'} and row['scope_rationale'].strip(),
                    'Declare required or excluded scope with its prospective rationale.')
            scope={'required':row['scope']=='required','rationale':row['scope_rationale']}
        result[scenario['scenario_id']]=scope
    return result


def validate_scope(assessment,scope):
    require(assessment is not None and scope==scope_from_assessment(assessment),
            'Assessment scope differs from its source-bound prospective review.')
