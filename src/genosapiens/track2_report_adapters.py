"""Read retained producer contracts into the common presentation contract."""
from dataclasses import asdict, is_dataclass

from genosapiens.track2_report_ranking import ROUTES


def plain(value):
    return asdict(value) if is_dataclass(value) else value


def study_inputs(bundle, plan, profiles, routes, lincs, evaluated_tables=None):
    """Keep target identity distinct from disease gene and scenario identity."""
    scenarios = {s['scenario_id']: s for s in (plan.get('mechanism_assessment') or {}).get('scenarios', [])}
    objectives = []
    for record in plan['tables']['target_records']:
        sid = record['context']['origin']['scenario_id']
        target = record['target']
        if bundle.data['input_origin'] == 'published_patient_genotype':
            symbols = sorted({a['gene'] for a in bundle.data['genotype']['alleles']})
            gene = symbols[0] if len(symbols) == 1 else None
        else:
            gene = scenarios.get(sid, {}).get('gene_id')
            if gene is None:
                construction = next((s for s in bundle.components['scenarios'] if s['scenario_id'] == sid), {})
                genes = construction.get('gene_ids', [])
                gene = genes[0] if len(genes) == 1 else None
        # An unresolved disease mapping must not become the intervention target.
        gene = gene or 'unresolved:'+sid
        objectives.append(dict(objective_id=record['objective']['objective_id'], scenario_id=sid,
            gene_id=gene, gene_symbol=gene if not gene.startswith(('ENSG', 'unresolved:')) else None,
            gene_identity_state='retained_stable_id' if gene.startswith('ENSG') else 'source_symbol_or_unresolved',
            target_ids=[target['target_id']], desired_modulation=target['required_modulation'],
            source_target=target, source_objective=record['objective']))
    by_id = {o['objective_id']:o for o in objectives}
    coverage = []
    for item in routes:
        item = plain(item)
        state = (item.get('attempt') or {}).get('state', 'not_applicable' if item.get('applicability') == 'not_applicable' else 'unavailable')
        coverage.append(dict(objective_id=item['objective_id'], route=item['route'],
            state={'success':'completed_with_candidates', 'attempted_complete':'completed_with_candidates', 'no_result':'completed_no_candidates',
                   'screened_exclusion':'completed_no_candidates'}.get(state, state), source_record=item))
    rows = {}
    def row(cid, oid, name, source):
        obj = by_id[oid]
        return rows.setdefault((cid, oid), dict(canonical_active_ingredient_id=cid, ingredient_key=cid,
            name=name, **{k:obj[k] for k in ('gene_id','gene_symbol','scenario_id','objective_id')},
            mechanism_nominated=False, lincs_nominated=False, cells=[], mechanisms=[],
            gates=[], eligibility=[], source_reference=source))
    for i, record in enumerate(plan['tables']['candidate_records']):
        candidate = record['candidate']
        cid = candidate['active_ingredient_id']
        for rationale in candidate['rationales']:
            oid = rationale['objective_id']
            r = row(cid, oid, cid, 'plan-source.json#/tables/candidate_records/'+str(i))
            r['mechanism_nominated'] |= rationale['route'] == ROUTES['mechanism']
            r['lincs_nominated'] |= rationale['route'] == ROUTES['lincs']
            r['gates'] = record.get('gate_results', [])
            r['eligibility'] = {k:record[k] for k in ('tier_decision','research_eligibility','conflicts') if k in record}
    for i, profile in enumerate(profiles):
        oid = profile['signed_rationale']['objective_id']
        if oid not in by_id:
            continue
        for j, hit in enumerate(profile['retrieved_drugs']):
            cid = hit['candidate_id']
            if (cid, oid) not in rows:
                continue
            r = rows[cid, oid]
            r['name'] = hit['name']
            r['mechanisms'].append(dict(mechanism_id=hit['drug_id'], targets=[profile['target_id']],
                action=hit['reported_action'], direction_class={'increase':'reported_positive','decrease':'reported_negative'}.get(
                    hit['action_direction'], 'unclassified'), references=hit['evidence_ids'], group_warning=None,
                source_reference='target-profiles.json#/'+str(i)+'/retrieved_drugs/'+str(j), source_record=hit))
    if evaluated_tables is not None:
        for index, candidate in enumerate(evaluated_tables.candidates):
            details = candidate.as_dict()['details']
            for (cid, oid), value in rows.items():
                if cid != candidate.row_id or by_id[oid]['scenario_id'] != candidate.origin.scenario_id:
                    continue
                value['gates'] = details['g0_g7']
                value['eligibility'] = {key:details[key] for key in ('priority_tier','research_eligibility','conflicts')}
                value['applicability_conflicts'] = details['conflicts']
                value['source_reference'] = 'tables.json#/tables/candidates/'+str(index)
    selected_cells = {}
    if lincs:
        decisions = lincs.get('scenario_candidate_decisions') or [dict(decision=d) for d in lincs['candidate_decisions']]
        for entry in decisions:
            decision = entry['decision']
            for oid in decision['objective_ids']:
                obj = by_id[oid]
                r = row(decision['candidate_id'], oid, decision.get('name') or decision['candidate_id'], 'lincs/results.json')
                r['lincs_nominated'] = True
                r['eligibility'].append(decision) if isinstance(r['eligibility'], list) else r['eligibility'].update(lincs=decision)
                bindings = [b for b in lincs.get('case_objective_bindings', []) if b['scenario_id'] == obj['scenario_id']]
                source_oids = {old for b in bindings for old, new in b['objectives'].items() if new == oid} or {oid}
                for cid, query in lincs['queries'].items():
                    if query['objective_id'] not in source_oids:
                        continue
                    method = query.get('method', 'paired_oe_kd')
                    selected_cells.setdefault(oid, []).append(dict(cell_id=cid, method=method,
                        source_reference='lincs/results.json#/queries/'+cid))
                    for index, medication in enumerate(lincs['cells'][cid]['medications']):
                        if medication['medication_entity_key'] == decision['medication_entity_key']:
                            r['name'] = medication['name']
                            r['cells'].append(dict(medication, cell_id=cid, method=method,
                                query_gene_count=sum(bool(g['selected']) for g in query['genes']),
                                query_sha256=query['query_sha256'],
                                source_reference='lincs/results.json#/cells/'+cid+'/medications/'+str(index)))
    projected = []
    for index, value in enumerate(rows.values()):
        projected.append(dict(value, source_record_reference=value['source_reference'],
            source_reference='supplements/nomination-inputs.json#/'+str(index)))
    return objectives, projected, coverage, selected_cells
