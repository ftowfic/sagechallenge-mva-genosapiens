"""Separate actual genotypes from explicitly supported mechanism alternatives."""
from dataclasses import dataclass
import copy
import re

from genosapiens.track2_adjudication import AdjudicationRecord
from genosapiens.track2_mechanism import (
    NON_VEP_SUPPORT, assess_mechanism_readiness, build_mechanism_model,
)
from genosapiens.track2_published_genotype import canonical, sha256
from genosapiens.track2_scoped_interpretation import require


JOINT_FIELDS = ('cellular_lesion_class', 'discriminating_baseline_lesion', 'cell_process',
    'tissue', 'developmental_stage', 'disease_stage', 'natural_history', 'temporal_phase',
    'response_shape', 'feedback')

PREDICTION_POLICY = 'track2-vep-direction-assumption/1.0.0'


def allele_prediction(allele, context, requirement_policy=None):
    """Select relevant transcripts before applying the declared damage heuristic."""
    from genosapiens.model_requirements import _transcripts, TranscriptScope
    from genosapiens.evidence_serialization import _plain
    available,gaps=_transcripts(context,allele,TranscriptScope())
    scopes=[]
    for rule in getattr(requirement_policy,'rules',()):
        scope=getattr(rule,'transcripts',None)
        # Projection changes copy derivation, not the original source identity.
        if (scope is not None and scope.transcript_ids
                and set(rule.selector.source_refs)<=set(context.evidence_refs) and any(
                rule.selector.requirement.record_ids==r.record_ids
                and rule.selector.requirement.mechanism_direction==r.mechanism_direction
                and rule.selector.requirement.consequence_requirement==r.consequence_requirement
                for r in context.requirements)):
            scopes.append(scope)
    if scopes:
        selected=[];missing=[]
        for scope in scopes:
            rows,errors=_transcripts(context,allele,scope)
            selected.extend(rows);missing.extend(errors)
        reason='applicable_disease_transcripts'
    else:
        selected=[r for r in available if any(s in (r[0].transcript_status or ())
            for s in ('MANE_Select','MANE_Plus_Clinical'))]
        reason='mane'
        if not selected:
            selected=[r for r in available if r[0].canonical is True]
            reason='canonical_fallback'
        missing=gaps
    unique={row[2]:row for row in selected}
    if len(unique)!=len(selected) and any(
        _plain(row[0])!=_plain(unique[row[2]][0]) for row in selected):
        missing=[*missing,'conflicting_transcript_annotations']
    selected=[unique[k] for k in sorted(unique)]
    predictions=[]
    for tx,raw,identifier,terms in selected:
        terms=tuple(terms or ())
        source_raw=(tx.source.source_fields or {}).get('csq',raw)
        sift=source_raw.get('SIFT');polyphen=source_raw.get('PolyPhen')
        def label(value):
            if value in (None,'','.'): return None
            match=re.fullmatch(r'([a-z_]+)(?:\(([0-9]+(?:\.[0-9]+)?(?:e-[0-9]+)?)\))?',value)
            if not match or match[2] is not None and not 0<=float(match[2])<=1:
                return 'malformed'
            return match[1]
        s,p=label(sift),label(polyphen)
        rule=None
        if set(terms)&{'stop_gained','frameshift_variant'}: rule='truncating_consequence_assumption'
        elif 'missense_variant' in terms:
            if s=='deleterious': rule='sift_deleterious_assumption'
            elif s is None and p=='probably_damaging': rule='absent_sift_polyphen_fallback'
        predictions.append({'transcript_id':identifier,'source':_plain(tx.source),'consequences':list(terms),
            'sift':sift,'polyphen':polyphen,'rule':rule,'direction':'loss_of_function' if rule else 'unknown',
            'predictor_disagreement':s=='deleterious' and p not in (None,'probably_damaging')})
    accepted=bool(predictions) and not missing and all(p['direction']=='loss_of_function' for p in predictions)
    return {'policy':PREDICTION_POLICY,'allele_id':allele.allele_id,'gene_id':context.gene_id,
        'selection_reason':reason,'transcripts':predictions,
        'alternatives':[_plain(tx) for tx,_,identifier,_ in available if identifier not in unique],
        'state':'predicted' if accepted else 'unresolved','direction':'loss_of_function' if accepted else 'unknown',
        'reasons':(['research_loss_of_function_assumption'] if accepted else
            ['missing_or_disagreeing_relevant_transcript_predictions']),
        'selection_gaps':[_plain(gap) for gap in missing],
        'transcript_fate':'unknown','residual_product':'unknown'}


@dataclass(frozen=True)
class ReviewedEffect:
    genotype_id: str
    model: dict
    evidence: tuple[dict, ...]
    adjudication: AdjudicationRecord
    assessments: dict | None = None

    def reviewed_bytes(self):
        return canonical({'genotype_id':self.genotype_id,'model':self.model,'evidence':list(self.evidence),
                          'assessments':self.assessments})

    @property
    def subject_id(self):
        return 'scoped-effect:'+sha256(self.reviewed_bytes())[:24]


@dataclass(frozen=True)
class ScopedMechanism:
    mechanism_id: str
    genotype_id: str
    model: object
    readiness: object
    evidence: tuple[dict, ...]
    adjudication: AdjudicationRecord | None
    grouping_key: dict
    limitations: tuple[str, ...]


def _evidence(allele):
    # One record per biological allele; physical duplicate source observations
    # remain inside its assertion and never count as independent functional data.
    from genosapiens.evidence_serialization import _plain
    assertion=canonical(_plain(allele)).decode('ascii')
    return {'evidence_id':'scoped-allele:'+sha256(assertion.encode())[:24],
        'source_identity':sha256(assertion.encode()),'lineage_key':allele.allele_id,
        'assertion':assertion,'evidence_role':'source_observation'}


def derived_evidence(prefix, value, lineage, role='prediction'):
    assertion=canonical(value).decode('ascii')
    return {'evidence_id':prefix+':'+sha256(assertion.encode())[:24],
        'source_identity':sha256(assertion.encode()),'lineage_key':lineage,
        'assertion':assertion,'evidence_role':role}


def unknown_model(hypothesis, alleles, evidence):
    """Record predictions without asserting transcript fate or drug direction."""
    models=[]
    gene=hypothesis.context.gene_id.split('.')[0]
    for key in hypothesis.allele_ids:
        transcripts=[t for t in alleles[key].transcript_consequences or ()
            if (t.source.source_fields or {}).get('csq',{}).get('Gene','').split('.')[0]==gene]
        consequences=sorted({c for t in transcripts for c in t.consequences or ()})
        models.append({'allele_id':key,'transcript_id':'unresolved',
            'vep_consequence':'&'.join(consequences) or 'unknown','vep_severity':'unknown',
            'transcript_fate':'unknown','residual_product':'unknown','molecular_function':'unknown',
            'evidence_support_types':['vep_prediction'] if consequences else ['source_observation'],
            'evidence_ids':[evidence[key]['evidence_id']]})
    joint={'allele_ids':list(hypothesis.allele_ids),'joint_effect':'unknown',
        **dict.fromkeys(JOINT_FIELDS,'unknown'),'evidence_support_types':['source_observation'],
        'evidence_ids':[evidence[key]['evidence_id'] for key in hypothesis.allele_ids]}
    return {'scenario_id':hypothesis.hypothesis_id,'allele_models':models,'joint_model':joint,
        'explicit_unknowns':['joint_model.'+name for name in ('joint_effect',*JOINT_FIELDS)]
            +[f'allele_models[{i}].'+name for i in range(len(models))
              for name in ('transcript_fate','residual_product','molecular_function')]}


def _reviewed(effect, hypothesis, evidence_origin):
    require(type(effect) is ReviewedEffect and effect.genotype_id==hypothesis.hypothesis_id
            and effect.model['scenario_id']==effect.genotype_id,'mechanism review has a foreign genotype')
    review=effect.adjudication
    require(type(review) is AdjudicationRecord and review.subject_id==effect.subject_id
            and review.domain=='mechanism' and review.decision=='accepted'
            and review.reviewed_content_sha256==sha256(effect.reviewed_bytes()),
            'mechanism effect lacks source-bound adjudication')
    model=build_mechanism_model(effect.model)
    require(effect.assessments is None or set(effect.assessments)=={
        'technical_support','exact_disease_validity','transcript_relevance'},
        'reviewed mechanism assessment fields differ from the current contract')
    require(set(model.joint_model.allele_ids)==set(hypothesis.allele_ids),
            'mechanism model does not match the actual genotype')
    evidence={e['evidence_id']:e for e in effect.evidence}
    require(len(evidence)==len(effect.evidence) and evidence and set(review.evidence_ids)==set(evidence))
    for row in evidence.values():
        require(set(row)=={'evidence_id','source_identity','lineage_key','assertion','evidence_role'}
            and all(isinstance(v,str) and v for v in row.values())
            and row['evidence_role'] in {'curated','synthetic'},'reviewed mechanism evidence is incomplete')
        require(row['evidence_role']!='synthetic' or evidence_origin=='independently_synthetic',
                'synthetic mechanism evidence cannot support actual case science')
    for record in (*model.allele_models,model.joint_model):
        require(record.evidence_ids and set(record.evidence_ids)<=set(evidence),
                'mechanism effect evidence lineage is incomplete')
    return model


def assemble_scoped_mechanisms(genotypes, interpreted, *, evidence_origin, reviewed_effects=()):
    """Preserve genotype/requirement/effect alternatives and existing readiness gates.

    Review can supply molecular effects; it cannot rewrite the observed copy,
    phase or model results. VCF QC alone does not establish technical validation.
    """
    from genosapiens.evidence_serialization import _plain
    require(genotypes.binding['interpretation_sha256']==sha256(canonical(_plain(interpreted))))
    by_id={a.allele_id:a for a in interpreted.alleles}
    observations={key:_evidence(allele) for key,allele in by_id.items()}
    generated={h.hypothesis.hypothesis_id:h for h in genotypes.hypotheses}
    alternatives={key:[] for key in generated}
    for effect in reviewed_effects:
        require(type(effect) is ReviewedEffect and effect.genotype_id in generated,
                'reviewed mechanism is outside the constructed genotype domain')
        _reviewed(effect,generated[effect.genotype_id].hypothesis,evidence_origin)
        alternatives[effect.genotype_id].append(effect)
    results={}
    for key,constructed in generated.items():
        hypothesis=constructed.hypothesis
        # Without supplied reviewed effects retain an explicit unknown branch.
        for effect in alternatives[key] or (None,):
            raw=unknown_model(hypothesis,by_id,observations) if effect is None else copy.deepcopy(effect.model)
            evidence=tuple(observations[a] for a in hypothesis.allele_ids)
            policy=getattr(genotypes,'requirement_policy',None)
            if policy is not None:
                from genosapiens.track2_functional_sources import scoped_policy
                policy=scoped_policy(policy,hypothesis.context.gene_id)
            evidence+=tuple(derived_evidence('scoped-prediction',
                allele_prediction(by_id[a],hypothesis.context,policy),a) for a in hypothesis.allele_ids)
            if policy is not None:
                evidence+=(derived_evidence('scoped-functional-sources',_plain(policy),
                    policy.policy_sha256,'source_observation'),)
            if effect is not None:evidence+=effect.evidence
            # Identity covers source requirements, genotype and distinct effects;
            # HPO, scores, ranks and processing partition order are absent.
            grouping={'disease_id':hypothesis.context.native_disease_id,'gene_id':hypothesis.context.gene_id,
                'requirements':_plain(hypothesis.context.requirements),'genotype_id':key,
                'effect_sha256':sha256(canonical(raw))}
            mechanism_id='scoped-mechanism:'+sha256(canonical(grouping))[:24]
            raw['scenario_id']=mechanism_id
            model=build_mechanism_model(raw)
            allele_facts=[]
            for item in model.allele_models:
                curated=bool(set(item.evidence_support_types)&NON_VEP_SUPPORT) and item.molecular_function!='unknown'
                predicted=item.vep_consequence!='unknown'
                allele_facts.append({'allele_id':item.allele_id,
                    'predicted_consequence':item.vep_consequence if predicted else None,
                    'curated_consequence':item.molecular_function if curated else None,
                    'curation_status':'curated' if curated else 'predicted' if predicted else 'unknown',
                    'evidence_ids':list(item.evidence_ids)})
            phase={'in_trans':'confirmed_trans','in_cis':'confirmed_cis',
                'not_applicable':'not_applicable'}.get(hypothesis.phase.value,'unknown')
            readiness=assess_mechanism_readiness({'scenario_id':mechanism_id,
                'genotype_completeness':hypothesis.completeness.value,'phase':phase,
                'technical_support':'uncertain','exact_disease_validity':'unknown',
                'inheritance_consistency':{'supported':'consistent','contradicted':'inconsistent'}.get(constructed.model.state,'unknown'),
                'transcript_relevance':'unknown','allele_evidence':allele_facts,
                'joint_evidence':'supported' if effect is not None and model.joint_model.joint_effect!='unknown' else 'unknown',
                'conflicts':['source phase conflict'] if hypothesis.phase.value=='conflicting' else [],
                'direction':model.joint_model.joint_effect,
                **({} if effect is None or effect.assessments is None else effect.assessments)})
            result=ScopedMechanism(mechanism_id,key,model,readiness,evidence,
                None if effect is None else effect.adjudication,grouping,
                (*constructed.limits,'Outside-scope and unsupported-class partners remain unknown',
                 *(('Technical validation, disease validity and transcript relevance are unassessed',)
                   if effect is None or effect.assessments is None else ())))
            # Identical effects are one branch, while different source attestations
            # cannot silently overwrite each other as the same reviewed result.
            if mechanism_id in results:
                require(results[mechanism_id]==result,'duplicate mechanism identity has conflicting provenance')
            results[mechanism_id]=result
    return tuple(results[key] for key in sorted(results))
