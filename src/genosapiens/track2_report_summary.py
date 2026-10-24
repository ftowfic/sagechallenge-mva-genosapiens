"""Readable report views over retained discovery and scientific assessments."""
from dataclasses import replace

from genosapiens.track2_lincs_priority import CAVEAT, GROUPS, md, priority_table
from genosapiens.track2_lincs_integration import research_priority
from genosapiens.track2_report_tables import ReportAssessment


GATE_LABELS = {
    'G0': 'Input identity and integrity', 'G1': 'Medication identity and current approval',
    'G2': 'Disease mechanism and model match', 'G3': 'Drug action in the required direction',
    'G4': 'Independent drug and rescue evidence', 'G5': 'Achievable exposure',
    'G6': 'Safety and hazards', 'G7': 'Validation plan',
}
ROUTE_LABELS = {
    'direct_target_modulation': 'Drugs acting directly on the target',
    'product_restoration_or_stabilization': 'Restoring or stabilizing the missing protein',
    'transcript_rescue': 'Repairing the gene message (RNA)',
    'upstream_regulation': 'Changing signals that control the target',
    'downstream_compensation': 'Compensating for effects downstream',
    'substrate_or_cofactor_modulation': 'Changing helper molecules or substrates',
    'pathway_bypass': 'Using an alternative biological pathway',
    'qualified_target_drug_linkage': 'Target–drug databases (Pharos)',
    'intervention_literature': 'Published drug studies',
    'lincs_cell_specific_expression': 'Cell response patterns (LINCS)',
}


def include_lincs_review(tables, result):
    """Expose discovery in the shared tables without fabricating gate decisions.

    Discovery gate rows are report-only assessments, explicitly not GateResult
    inputs to the scientific tier engine. Existing assessed rows are preserved.
    """
    reference = tables.target_hypotheses[0]
    diagnosis=reference.origin.input_origin=='diagnosis_guided_vcf'
    evidence = ('lincs-evidence:' + result['source_sha256'],)
    if not set(evidence) <= {r.evidence_id for r in tables.evidence_catalog}:
        raise ValueError('LINCS review tables lack source evidence')
    candidates, exclusions, gates = list(tables.candidates), list(tables.exclusions), list(tables.gates)
    decisions = {r['medication_entity_key']: r for r in result.get('candidate_decisions', [])}
    priorities = research_priority(result)
    unknown = ReportAssessment('not_assessed', 'Drug-specific review has not been completed.', evidence)
    for priority in priorities:
        entity = priority['medication_entity_key']
        decision = decisions.get(entity)
        candidate_id = decision['candidate_id'] if decision else 'lincs-screen:' + entity
        details = {'record_kind': 'discovery_review', 'name': priority['name'],
                   'active_ingredient_id': candidate_id, 'medication_entity_key': entity,
                   'route': 'lincs_cell_specific_expression', 'review_priority': priority,
                   'scientific_assessment': decision['scientific_assessment'] if decision else 'not_assessed',
                   'scientific_tier': decision['scientific_tier'] if decision else None,
                   'research_eligible': decision['research_eligible'] if decision else False,
                   'eligibility_bases': decision['eligibility_bases'] if decision else [],
                   'category_eligibility': decision['category_eligibility'] if decision else None,
                   'review_state': decision['review_state'] if decision else 'not_locally_nominated',
                   'case_matching': 'not_assessed', 'caveat': CAVEAT}
        base = replace(reference, row_id=candidate_id, subtype_applicability='unresolved',
                       translational_maturity='computational_inference', effect_magnitude='magnitude_unknown',
                       child_applicability='unresolved', general_mva_applicability='unresolved',
                       evidence_ids=evidence, exposure=unknown, safety=unknown,
                       decision_reasons=tuple(decision['reason_codes']) if decision else ('not_locally_nominated',), details=details)
        if decision:
            existing = next((i for i, row in enumerate(candidates) if row.row_id == candidate_id), None)
            if existing is None:
                candidates.append(base)
                for gate_id in GATE_LABELS:
                    gates.append(replace(base, row_id=f'discovery-gate:{candidate_id}:{gate_id}', details={
                        'record_kind': 'discovery_review', 'name': priority['name'],
                        'active_ingredient_id': candidate_id, 'medication_entity_key': entity,
                        'route': 'lincs_cell_specific_expression', 'gate_id': gate_id,
                        'gate_outcome': 'pass' if gate_id == 'G0' else 'not_assessed',
                        'gate_review_state': ('validated_diagnosis_input' if diagnosis else 'validated_public_input') if gate_id == 'G0' else 'pending_drug_specific_review',
                        'gate_reason': (('The diagnosis-guided input passed validation. This does not establish medication benefit.' if diagnosis else
                                        'The published-genotype input passed validation. This does not imply a match to a challenge case.') if gate_id == 'G0'
                                        else 'The report retains the FDA snapshot link. Review of the current product and warnings remains incomplete.' if gate_id == 'G1'
                                        else 'LINCS correlation cannot establish this requirement.'),
                        'coverage_state': 'complete' if gate_id == 'G0' else 'incomplete',
                    }))
            else:
                assessed = candidates[existing]
                candidates[existing] = replace(assessed, details={**assessed.details, 'lincs_review_priority': priority,
                    'combined_research_eligible':decision['research_eligible'],
                    'eligibility_bases':decision['eligibility_bases'],
                    'category_eligibility':decision['category_eligibility']})
        # Local exclusions and missing research assessments are distinct. A
        # discovery queue entry is not a contraindication or a scientific tier.
        if not details['research_eligible']:
            exclusions.append(replace(base, row_id='lincs-status:' + entity, details={**details,
                'exclusion_scope': 'research_shortlist_pending_review' if decision else 'local_expression_nomination',
                'exclusion_summary': ('The expression checks found a match. Independent medication review remains missing.' if decision and decision['scientific_assessment'] == 'not_assessed'
                                      else '; '.join(decision['reason_codes']) if decision
                                      else 'No assessed cell line passed both directional checks. There is no LINCS nomination.'),
            }))
    return replace(tables, candidates=tuple(candidates), exclusions=tuple(exclusions), gates=tuple(gates))


def route_coverage(plan, profiles, lincs):
    lines = ['## Route coverage', '',
             'Green means complete or evidence found. Amber means partial. Gray means not assessed. '
             'Evidence does not establish a suitable treatment.', '',
             '| Evidence type | Executed? | Any result? | What we found |', '|---|---|---|---|']
    for route in plan['routes']:
        attempt = route['attempt'] or {}
        state = attempt.get('state', 'not_applicable')
        executed = attempt.get('executed', False)
        execution = '🟡 Partial' if executed and state in {'truncated', 'unavailable', 'invalid'} else '🟢 Yes' if executed else '⚪ No'
        if not executed:
            obtained, summary = '⚪ Not assessed', 'Not included in this search.'
        elif route['route'] == 'lincs_cell_specific_expression' and lincs is not None:
            count = len(lincs.get('candidate_decisions', []))
            obtained = '🟢 Yes' if count else '⚪ No matches'
            eligible = sum(r['research_eligible'] for r in lincs.get('candidate_decisions', []))
            summary = f'{count} medications matched in at least one cell line. {eligible} are eligible for research. Scientific review of each medication remains separate.'
        elif route['route'] == 'intervention_literature':
            hits = {h['drug_id']: h for p in profiles for h in p['retrieved_drugs']}
            opposed = sum(h['action_comparison'] == 'opposed' for h in hits.values())
            obtained = '🟢 Yes' if hits else '🟡 Uncertain' if state in {'truncated', 'unavailable', 'invalid'} else '⚪ No hits'
            summary = f'The search found {len(hits)} named compounds. {opposed} act against the required direction. '
            summary += 'The literature search was partial.' if state == 'truncated' else 'See the source records for scope.'
        elif state == 'no_result':
            obtained, summary = '⚪ No hits', 'This source returned no medication linked to the target.'
        elif state in {'unavailable', 'invalid'}:
            obtained, summary = '🟡 Uncertain', 'The search did not return usable evidence.'
        else:
            obtained, summary = ('🟡 Partial' if state == 'truncated' else '🟢 Yes'), route['biological_rationale']
        lines.append(f"| {ROUTE_LABELS.get(route['route'], md(route['route']))} | {execution} | {obtained} | {md(summary)} |")
    return '\n'.join(lines)


def target_and_candidates(profiles, tables, lincs, *, assessment_status):
    lines = ['## Desired target, direction and identified candidates', '']
    for profile in profiles:
        direction = profile['required_modulation'].replace('_', ' ')
        lines += [f"**Target: {profile['target_id']}. Desired direction: {direction} functional activity.**", '',
                  'The goal is to restore the affected biological function. An expression pattern resembling more '
                  'target activity does not prove that a drug restores the protein or its function.', '',
                  'Functional protein: ' + profile['context']['residual_product_context'], '']
    assessed = [r for r in tables.candidates if r.details.get('record_kind') != 'discovery_review']
    eligible_ids = {r.row_id for r in assessed if r.details['research_eligibility'].eligible}
    if lincs is not None:
        eligible_ids.update(r['candidate_id'] for r in lincs.get('candidate_decisions', []) if r['research_eligible'])
    eligible = len(eligible_ids)
    if assessment_status == 'not_performed':
        lines += ['**Medication assessment was not performed.**',
                  'No completed negative result is available.', '']
    elif assessment_status == 'incomplete':
        lines += ['**Medication assessment is incomplete.**',
                  f'Observed research-eligible candidates so far: {eligible}. Required work remains incomplete.', '']
    if lincs is not None:
        priorities = [r for r in research_priority(lincs) if r['passing_cells']]
        repeated = sum(r['passing_cells'] >= 2 for r in priorities)
        lines += [f'**{len(priorities)} LINCS discovery candidates**, including {repeated} that pass in multiple cell lines. '
                  f'**{eligible} research-eligible candidates** across both eligibility routes.', '', CAVEAT, '',
                  'The table shows the first 20 discovery candidates in review order. '
                  '[Complete ranked list, including non-matches](lincs/review-priority.md) · '
                  '[Download the ranking](lincs/review-priority.csv).', '',
                  'Colors show agreement with the desired direction. **Review score = passing cell lines × sum of finite local joint scores**, '
                  'including negative scores. Agreement groups determine the order before score size.\n\n✓ means a local pass. × means a local failure. '
                  '“Direction +” counts cells that match the desired pattern alone. A pass requires both checks. '
                  'Missing tests stay visible in the coverage count. [Full ordering explanation](#how-to-read-the-medication-order).', '',
                  priority_table(priorities[:20], sorted(lincs['cells'])), '']
    elif assessment_status == 'complete':
        lines += [f'**{eligible} research-eligible candidates** in the assessed evidence.', '']
    lines += ['### Candidates and compounds found through other evidence', '',
              '| Name | Evidence type | Direction / result | Research status |', '|---|---|---|---|']
    seen = set()
    lincs_decisions = {r['candidate_id']: r for r in (lincs or {}).get('candidate_decisions', [])}
    for profile in profiles:
        for hit in profile['retrieved_drugs']:
            if hit['drug_id'] in seen: continue
            seen.add(hit['drug_id'])
            direction = {'opposed': 'Opposes the desired direction', 'direction_unknown': 'Direction unclear',
                         'functional_product_not_established': 'Direction matches. Functional protein remains unconfirmed',
                         'direction_aligned_requires_independent_action_evidence': 'Direction matches. Independent evidence remains necessary'}[hit['action_comparison']]
            combined = lincs_decisions.get(hit.get('candidate_id'))
            status = ('Eligible by LINCS agreement. Scientific assessment remains separate'
                      if combined and 'lincs_repeated_consistent' in combined['eligibility_bases']
                      else 'Eligible for research' if hit['research_eligible'] else 'Not eligible on current evidence')
            lines.append(f"| {md(hit['name'])} | Published study | {direction}: {md(hit['reported_action'])} | {status} |")
    for row in assessed:
        if any(h.get('candidate_id') == row.row_id for p in profiles for h in p['retrieved_drugs']): continue
        status = ('Eligible by LINCS agreement. Scientific assessment remains separate'
                  if 'lincs_repeated_consistent' in row.details.get('eligibility_bases', [])
                  else 'Eligible for research' if row.details['research_eligibility'].eligible else 'Not eligible on current evidence')
        lines.append(f"| {md(row.row_id)} | Assessed drug evidence | {md(row.details['priority_tier'].tier)} | {status} |")
    if not seen and not assessed:
        lines.append('| None | Other assessed routes | No additional candidate found | — |'
            if assessment_status == 'complete' else
            '| Not established | Medication assessment | Assessment remains unfinished | Not a negative result |')
    return '\n'.join(lines)


def exploratory_summary(plan,result,assessment):
    candidates=result['concordance']
    return '\n\n'.join(['# Exploratory repurposing research',
        'Assessment status: **'+assessment.status.replace('_',' ')+ '**.',
        str(len(candidates))+' medication and objective combinations meet a discovery nomination rule.',
        'These nominations are research hypotheses. They do not establish treatment benefit, patient rescue, or clinical safety.',
        'Both routes use the same ingredient, gene, scenario, and objective to establish overlap.',
        'Overlap does not establish independent confirmation. Single-route and single-cell nominations remain visible.'])


def review_summary(tables, lincs):
    lines = ['## Candidates, exclusions and G0–G7 checks', '', CAVEAT, '',
             'A discovery candidate is a lead to investigate. “Not assessed” means evidence is missing. '
             'It does not mean a medication failed a safety test.', '']
    if lincs is not None:
        rows = research_priority(lincs)
        lines += ['| LINCS result | Medications | What happens next |', '|---|---:|---|']
        for group, (color, label) in GROUPS.items():
            matching = [r for r in rows if r['agreement_group'] == group]
            if matching:
                lines.append(f'| {color} {label} | {len(matching)} | '
                             + ('Eligible by expression agreement. Scientific review remains necessary' if group == 'repeated_consistent'
                                else 'Drug-specific evidence review' if matching[0]['passing_cells'] else 'Outside the LINCS nomination queue') + ' |')
    lines += ['', '| Check | LINCS discovery candidates without drug-specific review |', '|---|---|']
    for gate, label in GATE_LABELS.items():
        summary = ('🟢 Public input passed validation. Matching to a challenge case remains unassessed.' if gate == 'G0'
                   else '🟡 The report links the FDA snapshot. Full review of current eligibility remains pending.' if gate == 'G1'
                   else '⚪ Not assessed for these drugs.')
        lines.append(f'| {gate}: {label} | {summary} |')
    lines += ['', 'Existing drug-specific gate decisions, including exclusions or hazards, remain unchanged and visible. '
              'Category-only membership does not override a scientific finding. '
              'The supplementary Candidates, Exclusions and G0–G7 tables include every LINCS nomination, '
              'its source, local results and missing assessments. [Structured tables](tables.json).']
    return '\n'.join(lines)
