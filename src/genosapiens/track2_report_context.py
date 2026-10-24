"""Source-bound opening context for both report origins, without new inference."""
import json
import re

from genosapiens.track2_operator_inputs import require
from genosapiens.track2_scientific_storage import json_identity


def anchor(identity):
    return 'context-' + json_identity(identity)['sha256'][:20]


def display(value):
    if value is None or value == '' or value == []:
        return 'Unavailable'
    if isinstance(value, (tuple, list)):
        return ', '.join(display(v) for v in value)
    return str(value).replace('&', '&amp;').replace('<', '&lt;').replace('>', '&gt;').replace('|', '&#124;').replace('\n', ' ')


def _records(evidence, prefix):
    return [(row, json.loads(row['assertion'])) for row in evidence
            if row['evidence_id'].startswith(prefix + ':')]


def predicted_effect(prediction, transcript):
    terms = set(transcript.get('consequences') or [])
    source = (transcript.get('source') or {}).get('source_fields') or {}
    csq = source.get('csq', source)
    position = csq.get('Protein_position')
    if not position:
        match = re.search(r'(\d+)(?:Ter|\*)', csq.get('HGVSp') or '')
        position = match[1] if match else None
    if 'stop_gained' in terms:
        text = 'The prediction introduces an early stop' + (f' at amino-acid position {display(position)}.' if position else '.')
    elif 'frameshift_variant' in terms:
        text = 'The prediction changes the protein reading frame.'
    elif 'missense_variant' in terms:
        text = 'The prediction changes an amino acid.'
    else:
        text = 'The functional effect remains unresolved.'
    return text + ' This prediction does not establish transcript decay, protein absence, or measured loss of function.'


def _detected_observations(parts):
    evidence = parts['evidence']
    variants = {}
    for variant in parts['variants']:
        for eid in variant['evidence_ids']:
            variants.setdefault(eid, []).append(variant)
    predictions = {}
    for record, value in _records(evidence, 'scoped-prediction'):
        predictions.setdefault(value['allele_id'], []).append(dict(evidence_id=record['evidence_id'], **value))
    observed = []
    for record, allele in _records(evidence, 'scoped-allele'):
        related = variants.get(record['evidence_id'], [])
        selected = predictions.get(allele['allele_id'], [])
        transcripts = []
        for prediction in selected:
            for tx in prediction['transcripts']:
                raw = (tx.get('source') or {}).get('source_fields') or {}
                csq = raw.get('csq', raw)
                transcripts.append(dict(transcript_id=tx['transcript_id'], coding=csq.get('HGVSc'),
                    protein=csq.get('HGVSp'), exon=csq.get('EXON'), intron=csq.get('INTRON'),
                    amino_acid_position=csq.get('Protein_position'), selection_reason=prediction['selection_reason'],
                    evidence_id=prediction['evidence_id'], source=tx.get('source'),
                    interpretation_state=prediction['state'], direction=tx['direction'],
                    plain_effect=predicted_effect(prediction, tx)))
        observed.append(dict(observation_id=allele['allele_id'], label='Detected variant',
            source_record_ids=[v['source_record_id'] for v in related],
            source_allele_ids=[v['source_allele_id'] for v in related],
            package_allele_ids=[v['allele_id'] for v in related],
            gene_ids=sorted({g for v in related for g in v['gene_ids']}),
            evidence_id=record['evidence_id'], identities=allele['identities'],
            mapping_state='retained' if allele['identities'] else 'unavailable',
            calls=allele.get('calls'), transcripts=transcripts,
            all_transcript_consequences=allele.get('transcript_consequences'), predictions=selected,
            quality=[(c['source'].get('source_fields') or {}).get('automatic_quality') for c in allele.get('calls') or []],
            source=record))
    return observed


def retained_scenarios(bundle, plan):
    """Select exactly the hypotheses displayed by the context projection."""
    assessment = plan.get('mechanism_assessment')
    if assessment:
        return assessment['scenarios']
    if bundle.data['input_origin'] == 'diagnosis_guided_vcf':
        return bundle.components['scenarios']
    require(bundle.data['input_origin'] == 'published_patient_genotype',
            'The report context origin is unsupported.')
    return [dict(scenario_id=sid, mechanism_inputs=facts)
            for sid, facts in plan['mechanism_inputs'].items()]


def retained_gene_ids(bundle, plan, objectives):
    """Use exact scenario and objective identifiers, without alias inference."""
    return sorted({row['gene_id'] for row in [*retained_scenarios(bundle, plan), *objectives]
                   if row.get('gene_id')})


def mechanism_blockers(scenario, source_evidence):
    """Project retained checks without interpreting missing evidence as a search."""
    facts = scenario['mechanism_inputs']
    acceptance = facts.get('prediction_acceptance') or {}
    details = acceptance.get('evidence') or {}
    checks = details.get('checks') or {}
    proof = json.loads(details['scientific_proof']) if details.get('scientific_proof') else {}
    evidence_ids = [key for key in scenario.get('evidence_ids', []) if key in source_evidence]
    rows = []
    def add(category, state, reasons, path):
        if state != 'accepted':
            rows.append(dict(category=category, state=state, reasons=sorted(set(reasons)),
                evidence_ids=evidence_ids, evidence_path=path))
    for category, key, fallback in (
            ('disease_scope', 'inheritance', 'exact_disease_validity'),
            ('genotype_completeness', 'genotype', 'genotype_completeness')):
        check = checks.get(key)
        if check:
            add(category, check['state'], check.get('reasons', []), 'prediction_acceptance.evidence.checks.'+key)
        elif facts.get(fallback) not in {'exact_supported', 'complete'}:
            add(category, 'unresolved', ['retained_'+fallback+'_'+str(facts.get(fallback) or 'unavailable')], fallback)
    functional = proof.get('functional') or {}
    for allele, value in functional.items():
        for origin in ('native', 'reviewed_import'):
            result = value.get(origin) or {}
            for gap in [gap for check in [result, *result.get('checks', [])]
                        for gap in check.get('gaps') or []]:
                category = ('functional_interpretation' if gap in {
                    'functional_assertion_requires_interpretation', 'effect_scope_or_negation_unresolved',
                    'experimental_source_citation_unavailable'
                } else 'protein_identity')
                add(category, 'unresolved', [gap], 'scientific_proof.functional.'+allele+'.'+origin)
            if result.get('contrary'):
                add('functional_interpretation', 'conflicting', ['contrary_functional_evidence'],
                    'scientific_proof.functional.'+allele+'.'+origin)
    if not functional:
        add('protein_identity', 'unavailable', ['retained_protein_identity_check_unavailable'], 'prediction_acceptance.evidence')
        add('functional_interpretation', 'unavailable', ['retained_functional_interpretation_unavailable'], 'prediction_acceptance.evidence')
    elif checks.get('functional_sources', {}).get('state') != 'accepted' and not any(
            row['category'] == 'functional_interpretation' for row in rows):
        add('functional_interpretation', 'unresolved', ['functional_check_incomplete'],
            'prediction_acceptance.evidence.checks.functional_sources')
    return rows


def project_context(bundle, plan, objectives, gene_metadata=None):
    """Project authenticated inputs. A missing optional mapping remains missing."""
    origin = bundle.data['input_origin']
    metadata = gene_metadata or {}
    scenarios = retained_scenarios(bundle, plan)
    if origin == 'diagnosis_guided_vcf':
        parts = bundle.components
        observations = _detected_observations(parts)
        constructions = {s['scenario_id']: s for s in parts['scenarios']}
        scope = parts['scope']
        source_evidence = {e['evidence_id']: e for e in parts['evidence']}
        patients = []
    else:
        require(origin == 'published_patient_genotype', 'The report context origin is unsupported.')
        observations = []
        genotype = bundle.data['genotype']
        for allele in genotype['alleles']:
            observations.append(dict(observation_id=allele['allele_id'], label='Reported variant',
                gene_symbol=allele['gene'], gene_ids=[], reported_coding=allele['reported_coding'],
                reported_protein=allele['reported_protein'], normalized_coding=allele['normalized_coding'],
                normalized_protein=allele['normalized_protein'], transcript_id=allele['transcript_id'],
                transcript_version=allele['transcript_version'], reference_mapping=allele['reference_mapping'],
                normalization_rationale=allele['normalization_rationale'], copies=allele['copies'],
                zygosity=genotype['zygosity'], phase=genotype['phase'], evidence_ids=allele['evidence_ids'],
                identities=[], mapping_state='not_supplied_by_published_contract',
                quality=None, source=allele))
        constructions = {}
        scope = None
        source_evidence = {e['evidence_id']: e for e in bundle.data['evidence']}
        patients = bundle.data['patients']
    source_evidence.update({e['evidence_id']: e for e in
                            (plan.get('mechanism_assessment') or {}).get('evidence', [])})
    hypotheses = []
    for scenario in scenarios:
        sid = scenario['scenario_id']
        local_objectives = [o for o in objectives if o['scenario_id'] == sid]
        construction = constructions.get(scenario.get('construction_scenario_id', sid))
        gene_ids = sorted({o['gene_id'] for o in local_objectives if o.get('gene_id')} |
                          ({scenario['gene_id']} if scenario.get('gene_id') else set()))
        genes = [dict(gene_id=gid, reference_metadata=metadata.get(gid),
                      metadata_state=metadata[gid].get('mapping_state', 'available') if gid in metadata else 'unavailable',
                      symbol=(metadata.get(gid) or {}).get('symbol') or
                        next((o.get('gene_symbol') for o in local_objectives if o.get('gene_id') == gid), None))
                 for gid in gene_ids]
        proof_records = [source_evidence[key] for key in (construction or {}).get('evidence_ids', [])
                         if key.startswith('scoped-model:')]
        proofs = [json.loads(r['assertion']) for r in proof_records]
        allele_ids = set((construction or {}).get('allele_ids', []))
        linked = [o['observation_id'] for o in observations if origin == 'published_patient_genotype'
                  or allele_ids.intersection(o['package_allele_ids'])]
        hypotheses.append(dict(scenario_id=sid, anchor=anchor(sid), label='Hypothesis '+str(len(hypotheses)+1),
            construction_scenario_id=(construction or {}).get('scenario_id'), genes=genes,
            observation_ids=linked, mechanism_inputs=scenario['mechanism_inputs'],
            model=scenario.get('model'), construction_proofs=proofs,
            construction_evidence_ids=[r['evidence_id'] for r in proof_records],
            objectives=local_objectives, limitations=scenario.get('limitations', [])))
        hypotheses[-1]['mechanism_blockers'] = mechanism_blockers(scenario, source_evidence)
    return dict(schema='genosapiens.track2-report-context/v1', input_origin=origin,
        input_sha256=bundle.input_sha256, observations=observations, hypotheses=hypotheses,
        scope=scope, patient_attribution=patients,
        selection_coverage=bundle.components['selection']['coverage'] if origin == 'diagnosis_guided_vcf' else None,
        source_evidence=list(source_evidence.values()))


def _context_facts(context):
    """Render the opening summary from projections, keeping full details linked."""
    yield '## Gene and variant context\n\n'
    yield '[Complete allele, transcript, quality, and hypothesis records](supplements/context.json).\n\n'
    observations = {o['observation_id']: o for o in context['observations']}
    displayed_observations = set()
    for hypothesis in context['hypotheses']:
        yield '<a id="'+hypothesis['anchor']+'"></a>\n\n### '+hypothesis['label']+'\n\n'
        for gene in hypothesis['genes']:
            metadata = gene['reference_metadata'] or {}
            yield 'Gene: '+display(gene['symbol'])+'. Stable identifier: `'+display(gene['gene_id'])+'`.\n\n'
            yield 'Full gene name: '+display(metadata.get('name'))+'. HGNC identifier: '+display(metadata.get('hgnc_id'))+'.\n\n'
            locus = metadata.get('locus')
            yield 'Gene metadata mapping: '+display(gene['metadata_state'])+'.\n\n'
            if locus:
                yield 'Gene locus: '+display(locus['chromosome'])+':'+display(locus['start'])+'–'+display(locus['end'])+'.\n\n'
                yield 'Gene assembly: '+display(locus['assembly'])+'. Coordinate convention: '+display(locus['coordinate_system'])+'.\n\n'
            else:
                yield 'Gene locus: unavailable.\n\n'
            yield 'The gene locus does not identify the allele location.\n\n'
        for oid in hypothesis['observation_ids']:
            row = observations[oid]
            if oid in displayed_observations:
                yield '[Shared '+row['label'].lower()+'](#'+anchor(oid)+'). This hypothesis retains its own interpretation.\n\n'
                continue
            displayed_observations.add(oid)
            yield '<a id="'+anchor(oid)+'"></a>\n\n'
            yield row['label']+': `'+display(oid)+'`.\n\n'
            if row['label'] == 'Reported variant':
                yield display(row['gene_symbol'])+': '+display(row['reported_coding'])+' ('+display(row['reported_protein'])+').\n\n'
                yield 'Transcript: '+display(row['transcript_id'])+'. Transcript version: '+display(row['transcript_version'])+'.\n\n'
                yield 'Reported copies: '+display(row['copies'])+'. Phase: '+display(row['phase'])+'. Detection quality: unavailable.\n\n'
                yield 'Reference mapping (retained source text): '+display(row['reference_mapping'])+'\n\n'
            else:
                for identity in row['identities']:
                    location = identity['identity']
                    yield 'Allele location ('+display(identity['stage'])+'): '+display(location.get('contig'))+':'+display(location.get('position'))+'.\n\n'
                    yield 'Reference: '+display(location.get('reference'))+'. Alternate: '+display(location.get('alternate'))+'. Coordinates use one-based closed positions.\n\n'
                    yield 'Reference assembly identity: `'+display(location.get('reference_identity_id'))+'`. Reference role: `'+display(location.get('reference_role'))+'`.\n\n'
                for tx in row['transcripts']:
                    yield 'Transcript: `'+display(tx['transcript_id'])+'`. Selection rule: `'+display(tx['selection_reason'])+'`.\n\n'
                    yield 'Coding change: '+display(tx['coding'])+'. Protein change: '+display(tx['protein'])+'.\n\n'
                    yield 'Exon: '+display(tx['exon'])+'. Intron: '+display(tx['intron'])+'.\n\n'
                    yield tx['plain_effect']+'\n\n'
                if not row['transcripts']:
                    yield 'Transcript effect: unavailable. This report does not select a replacement transcript.\n\n'
                for quality in row['quality']:
                    quality = quality or {}
                    yield 'Quality state: '+display(quality.get('state'))+'. Genotype: '+display(quality.get('genotype'))+'.\n\n'
                    yield 'Genotype quality: '+display(quality.get('genotype_quality'))+'. Depth: '+display(quality.get('depth'))+'. Allele balance: '+display(quality.get('allele_balance'))+'.\n\n'
        facts = hypothesis['mechanism_inputs']
        yield 'Mechanism direction: '+display(facts.get('direction'))+'. Evidence status: '+display(facts.get('joint_evidence'))+'.\n\n'
        yield 'Phase: '+display(facts.get('phase'))+'. Genotype completeness: '+display(facts.get('genotype_completeness'))+'.\n\n'
        for blocker in hypothesis.get('mechanism_blockers', []):
            yield 'Mechanism blocker: '+display(blocker['category'].replace('_', ' '))+'. State: '+display(blocker['state'])+'.\n\n'
            yield 'Retained reason: `'+display(blocker['reasons'])+'`. [Retained checks and source evidence](supplements/context.json).\n\n'
        if hypothesis.get('mechanism_blockers'):
            yield 'Unavailable evidence does not establish a completed negative search.\n\n'
        for proof in hypothesis['construction_proofs']:
            genotype = proof['genotype']
            requirements = genotype['hypothesis']['context']['requirements']
            yield 'Required affected copies: '+display([r.get('required_allele_count') for r in requirements])+'.\n\n'
            for member in genotype['model']['members']:
                prop = member['genotype']
                yield 'Observed ALT copies: '+display(prop.get('minimum_alt_copies'))+' to '+display(prop.get('maximum_alt_copies'))+'.\n\n'
        for objective in hypothesis['objectives']:
            yield 'Intervention objective: `'+display(objective['objective_id'])+'`.\n\n'
            yield 'Desired modulation: '+display(objective.get('desired_modulation'))+'. This objective is a research hypothesis.\n\n'
        if not hypothesis['objectives']:
            yield 'No retained intervention objective exists for this hypothesis.\n\n'
    yield 'Alternative hypotheses remain separate. Their shared observations do not establish multiple patients or independent evidence.\n\n'


def context_markdown(context):
    """Keep facts in tight labeled groups without changing their content."""
    group = None
    for chunk in _context_facts(context):
        if chunk.startswith(('## ', '<a ', '[Complete', 'Alternative hypotheses')):
            if group is not None:
                yield '\n'
            group = None
            yield chunk
            continue
        if chunk.startswith(('Gene:', 'Full gene name:', 'Gene metadata', 'Gene locus:', 'Gene assembly:', 'The gene locus')):
            current = 'Gene and location'
        elif chunk.startswith(('Mechanism direction:', 'Phase:', 'Required affected', 'Observed ALT')):
            current = 'Genotype and interpretation'
        elif chunk.startswith(('Intervention objective:', 'Desired modulation:', 'No retained intervention')):
            current = 'Objective'
        else:
            current = 'Variant observations'
        if current != group:
            yield '\n**'+current+'**\n\n'
            group = current
        yield '- '+chunk.rstrip()+'\n'
    yield '\n'
