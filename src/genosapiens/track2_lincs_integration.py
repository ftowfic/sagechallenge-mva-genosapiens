"""LINCS discovery joins the existing route and review boundaries without granting gates."""
from dataclasses import asdict, replace
import json
from pathlib import Path

from genosapiens.track2_candidate_union import CandidateRouteObservation, build_candidate_union
from genosapiens.track2_drug_identity import DrugIdentityAssertion, resolve_drug_identity
from genosapiens.track2_lincs import ROUTE, load_source_bundle, load_cell, read_csv, resolve_allowlist
from genosapiens.track2_lincs_scoring import refine_query, score_cell, cross_cell_summary
from genosapiens.track2_published_genotype import _source_file, canonical, sha256
from genosapiens.track2_lincs_priority import ELIGIBILITY_POLICY, review_priority


def category_eligibility(result):
    """Replay membership from local results, never from a supplied category flag.

    Call only after source/query admission and local scoring. The shared
    composer verifies those inputs before this research-membership boundary.
    """
    from genosapiens.track2_lincs_scoring import cross_cell_summary, summarize_conditions
    if result['cross_cell'] != cross_cell_summary(result['cells']):
        raise ValueError('LINCS cross-cell evidence differs from local results')
    queries = {cell: q['query_sha256'] for cell, q in sorted(result['queries'].items())}
    if set(queries) != set(result['cells']):
        raise ValueError('LINCS eligibility query coverage is incomplete')
    for cell, local in result['cells'].items():
        if (local['query_sha256'] != queries[cell]
                or local['medications'] != summarize_conditions(local['conditions'])):
            raise ValueError('LINCS eligibility differs from originating condition evidence')
    return {row['medication_entity_key']: {
        'eligible': row['agreement_group'] == ELIGIBILITY_POLICY['agreement_group'],
        'basis': ELIGIBILITY_POLICY['basis'], 'policy': ELIGIBILITY_POLICY,
        'agreement_group': row['agreement_group'], 'source_sha256': result['source_sha256'],
        'query_sha256s': queries, 'passing_cells': row['passing_cells'],
        'failed_cells': row['failed_cells'], 'conflicting_cells': row['conflicting_cells'],
        'missing_cells': row['missing_cells'], 'tested_cell_count': row['tested_cell_count'],
        'selected_cell_count': row['selected_cell_count'],
    } for row in review_priority(result)}


def research_priority(result):
    """One export view over replayed category and independent evidence decisions."""
    category = category_eligibility(result)
    decisions = {row['medication_entity_key']:row for row in result.get('candidate_decisions', [])}
    rows = []
    for priority in review_priority(result):
        entity = priority['medication_entity_key']
        decision = decisions.get(entity)
        scientific = decision.get('evidence_review_eligibility') if decision else None
        bases = ([ELIGIBILITY_POLICY['basis']] if category[entity]['eligible'] else [])
        if scientific and scientific['eligible']:
            bases.append('reviewed_scientific_evidence')
        if decision and (decision['category_eligibility'] != category[entity]
                         or decision['eligibility_bases'] != bases
                         or decision['research_eligible'] != bool(bases)):
            raise ValueError('LINCS research eligibility does not reproduce local evidence')
        rows.append({**priority,'category_eligibility':category[entity],
                     'research_eligible':bool(bases),'eligibility_bases':bases})
    return rows


def bound_inputs(root, config):
    if set(config) != {'sources_path','sources_sha256','queries'}:
        raise ValueError('incomplete LINCS plan binding')
    path = _source_file(Path(root), config['sources_path'])
    if sha256(path.read_bytes()) != config['sources_sha256']:
        raise ValueError('LINCS source manifest changed')
    manifest, paths = load_source_bundle(root, path)
    if set(config['queries']) != {c['cell_id'] for c in manifest['cells']}:
        raise ValueError('LINCS selected cells and query bindings disagree')
    bindings = {'lincs:sources':config['sources_sha256']}
    bindings.update({'lincs:'+key:value['sha256'] for key,value in manifest['inputs'].items()})
    queries = {}
    for cell, binding in config['queries'].items():
        if set(binding) != {'path','sha256'}:
            raise ValueError('incomplete frozen query binding')
        raw = _source_file(Path(root), binding['path']).read_bytes()
        if sha256(raw) != binding['sha256']:
            raise ValueError('frozen LINCS query changed')
        queries[cell] = json.loads(raw)
        bindings['lincs:query:'+cell] = binding['sha256']
    return manifest, paths, queries, bindings


def execute(root, config, *, target, scenario_id):
    manifest, paths, queries, _ = bound_inputs(root, config)
    allowed, excluded = resolve_allowlist(read_csv(paths['allowlist']),read_csv(paths['products']),read_csv(paths['identity_audit']))
    results = {}
    for cell_config in manifest['cells']:
        cell = load_cell(cell_config, paths, allowed)
        replay = refine_query(cell, cell_config, scenario_id=scenario_id, target=target,
                              source_sha256=config['sources_sha256'])
        if replay != queries[cell['cell_id']]:
            raise ValueError('frozen local query disagrees with actual source/target/policy')
        results[cell['cell_id']] = score_cell(replay, cell)
    return {'schema':'genosapiens.track2-lincs-result/v1','queries':queries,'cells':results,
            'cross_cell':cross_cell_summary(results), 'identity_exclusions':excluded,
            'source_sha256':config['sources_sha256'], 'terms':manifest['terms'],
            'redistribution':manifest['redistribution']}


def route_attempt(result, row):
    """Replace declared unexecuted state with actual local execution accounting."""
    if row['route'] != ROUTE or row['attempt']['state'] != 'unattempted' or row['attempt']['executed']:
        raise ValueError('LINCS route must start as a declared unattempted route')
    queries = result['queries']
    if {q['objective_id'] for q in queries.values()} != {row['objective_id']}:
        raise ValueError('LINCS route objective mismatch')
    incomplete = any(q['state'] != 'ready' for q in queries.values())
    nominated = any(r['state']=='pass' for c in result['cells'].values() for r in c['medications'])
    return {**row,'attempt':{**row['attempt'], 'executed':True,
        'state':'unavailable' if incomplete else 'attempted_complete' if nominated else 'no_result',
        'query_ids':['lincs-query:'+q['query_sha256'] for _,q in sorted(queries.items())],
        'screening_exclusions':[f"{c}: {len(v['exclusions'])} QC/identity exclusions; "
            f"{sum(r['state']=='fail' for r in v['medications'])} failed and "
            f"{sum(r['state']=='unassessed' for r in v['medications'])} unassessed medication summaries"
            for c,v in sorted(result['cells'].items())]}}


def discovery_union(result, existing_candidates, assessed_rows):
    categorical = category_eligibility(result)
    observations = [r for c in existing_candidates for r in c.rationales]
    entity_ids = {}
    for cell_id, cell in sorted(result['cells'].items()):
        query = result['queries'][cell_id]
        nominated = {r['medication_entity_key'] for r in cell['medications'] if r['state']=='pass'}
        for condition in cell['conditions']:
            entity = condition['medication_entity_key']
            if entity not in nominated: continue
            source = condition['identity']
            assertion = DrugIdentityAssertion('drugs-at-fda',entity,source['fda_ingredient_names'],
                'active_ingredient',exact_substance_key=entity,active_moiety_keys=(entity,))
            identity = resolve_drug_identity('drugs-at-fda:'+entity,(assertion,))
            identity = replace(identity,
                source_identifiers=tuple(sorted({*identity.source_identifiers,'lincs:'+condition['pert_id']})),
                product_identifiers=tuple(sorted({p['application']+':'+p['product_number'] for p in source['products']})))
            entity_ids[entity] = identity.active_ingredient_ids[0]
            # One discovery rationale per local query and entity; same release
            # is one lineage, not one independent experiment per dose or cell.
            observations.append(CandidateRouteObservation(identity,query['objective_id'],ROUTE,
                'lincs-rationale:'+query['query_sha256']+':'+entity,
                'Exploratory same-cell expression match; independent drug action and rescue unassessed.',
                (query['target_id'],'lincs-query:'+query['query_sha256'],entity),
                ('lincs-evidence:'+result['source_sha256'],),('GSE92742',)))
    union = build_candidate_union(observations)
    assessed = {row.row_id:row for row in assessed_rows}
    decisions = []
    for entity,candidate_id in sorted(entity_ids.items()):
        row = assessed.get(candidate_id)
        decision = row.details['research_eligibility'] if row else None
        category = categorical[entity]
        bases = ([ELIGIBILITY_POLICY['basis']] if category['eligible'] else [])
        if decision and decision.eligible:
            bases.append('reviewed_scientific_evidence')
        decisions.append({'medication_entity_key':entity,'candidate_id':candidate_id,
            'objective_ids':sorted({result['queries'][c]['objective_id'] for c in result['queries']}),
            'local_nominated_cells':[r['cell_id'] for c in result['cells'].values() for r in c['medications']
                                     if r['medication_entity_key']==entity and r['state']=='pass'],
            'scientific_tier':decision.scientific_tier if decision else None,
            'scientific_assessment':'assessed' if row else 'not_assessed',
            'research_eligible':bool(bases), 'eligibility_bases':bases, 'category_eligibility':category,
            'evidence_review_eligibility':asdict(decision) if decision else None,
            'reason_codes':(['lincs_repeated_agreement_no_observed_conflict'] if category['eligible'] else [])
                + (list(decision.reason_codes) if decision else [
                'independent_mechanism_and_drug_action_missing','current_warning_review_missing',
                'exposure_and_hazard_review_missing','validation_plan_missing']),
            'case_matching':'not_assessed', 'submission_ready':False,
            'gate_results':[asdict(g) for g in row.details['g0_g7']] if row else [],
            'review_state':'existing_source_bound_assessment' if row else 'pending_drug_specific_review'})
    return {'candidate_union':asdict(union),'candidate_decisions':decisions,
            'research_eligibility_policy':ELIGIBILITY_POLICY}
