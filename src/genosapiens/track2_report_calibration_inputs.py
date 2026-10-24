"""Calibration scopes from admitted complete evaluated evidence."""
import json
from pathlib import Path

from genosapiens.track2_operator_inputs import require
from genosapiens.track2_scientific_reader import read_scientific
from genosapiens.track2_scientific_storage import json_identity
from genosapiens.track2_derived_strength import METHODS,aggregate_strength
from genosapiens.track2_report_derived import condition_index,report_conditions
from genosapiens.track2_calibration_bank import scientific_scope
from genosapiens.track2_frozen_calibration import identity


def retained_scopes(root, *, eligibility, versions, classification, source_manifest_sha256):
    """Return complete preparation scopes and explicit rows for every evaluation.

    The caller admits the source package first. This extractor also checks query
    content, duplicate observations, objective scope, and evaluated membership.
    """
    root=Path(root)
    if eligibility['state']!='available':return dict(inputs=[],evaluations=[])
    require(versions['calibration']['state']=='available', 'Calibration lacks original source versions.')
    from genosapiens.track2_tissue_records import validate_tissue_record
    queries={}
    for record in read_scientific(root/'discovery/tissue-queries.json'):
        if record['result'] is None:continue
        query=record['result']['query'];validate_tissue_record(query,'query')
        key=query['query_sha256']
        require(key not in queries or json_identity(queries[key])==json_identity(query),
                'A retained query identity has conflicting content.')
        queries[key]=query
    objectives={row['objective_id']:row for row in eligibility['objectives']}
    expected={(oid,row['medication_entity_key']) for oid,obj in objectives.items() for row in obj['evaluated']}
    seen=set();evaluations=[];groups={}
    with condition_index(report_conditions(root)) as index:
        for score in read_scientific(root/'discovery/tissue-scores.json'):
            oid=score['objective_id'];ingredient=score['medication_entity_key'];key=(oid,ingredient)
            require(key in expected and key not in seen, 'Calibration evaluated scope differs from frozen eligibility.')
            seen.add(key);objective=objectives[oid]
            require(all(score[k]==objective[k] for k in ('gene_id','scenario_id')),
                    'Calibration cannot combine different objective genes.')
            for method in METHODS:
                coverage=[row for row in score['coverage'] if row['method']==method]
                if not coverage:continue
                tissues=sorted({row['tissue_id'] for row in coverage})
                selected={row['query_sha256'] for row in coverage if row.get('query_sha256')}
                rows=[]
                for (raw,) in index.execute('SELECT record FROM conditions WHERE ingredient=? AND gene=? AND method=? ORDER BY query,signature',
                        (ingredient,objective['gene_id'],method)):
                    row=json.loads(raw)
                    if row['query_sha256'] in selected and row['tissue_id'] in tissues:
                        validate_tissue_record(row,'condition');rows.append(row)
                incomplete=any(row.get('incomplete_jobs') for row in coverage)
                summary=aggregate_strength(rows,method=method,direction=objective['direction'],
                    declared_tissues=tissues,incomplete=incomplete)
                for field in ('n_conditions','n_valid','n_undefined'):
                    if all(field in row for row in coverage) and (field=='n_conditions' or objective['direction']=='increase'):
                        require(summary[field]==sum(row[field] for row in coverage),
                                'Calibration conditions differ from retained coverage.')
                evaluation=dict(objective_id=oid,scenario_id=objective['scenario_id'],gene_id=objective['gene_id'],
                    medication_entity_key=ingredient,method=method,direction=objective['direction'],
                    status=summary['status'],nominated=score['nominated'],provisional=score['provisional'],
                    declared_tissues=tissues,missing_tissues=summary['missing_tissues'],scope_sha256=None)
                evaluations.append(evaluation)
                if summary['status']!='complete':continue
                needed={row['query_sha256'] for row in rows}
                require(needed<=queries.keys(), 'A calibration condition lacks its retained query.')
                for row in rows:
                    query=queries[row['query_sha256']]
                    require((query['gene_id'],query['tissue_id'],query['method'])==
                        (row['gene_id'],row['tissue_id'],row['method']),
                        'Calibration conditions differ from their scientific query scope.')
                group=groups.setdefault((oid,method),dict(queries={},conditions=[],evaluations=[]))
                group['queries'].update({key:queries[key] for key in needed})
                group['conditions'].extend(rows);group['evaluations'].append(evaluation)
    require(seen==expected, 'Calibration is missing evaluated ingredients.')
    deduplicated={}
    for (oid,method),group in groups.items():
        objective=objectives[oid]
        scope=scientific_scope(list(group['queries'].values()),group['conditions'],direction=objective['direction'],
            eligible_ingredients=eligibility['ingredients'],
            fda_snapshot_sha256=versions['calibration']['fda_snapshot_sha256'],
            source_release=versions['calibration']['source_release'])
        digest=identity(scope)
        for evaluation in group['evaluations']:evaluation['scope_sha256']=digest
        provenance=dict(objective_id=oid,scenario_id=objective['scenario_id'],gene_id=objective['gene_id'],method=method)
        if digest in deduplicated:
            deduplicated[digest]['provenance']['objectives'].append(provenance)
            continue
        deduplicated[digest]=dict(schema='genosapiens.track2-calibration-inputs/v1',
            queries=list(group['queries'].values()),conditions=group['conditions'],direction=objective['direction'],
            eligible_ingredients=eligibility['ingredients'],fda_snapshot_sha256=scope['fda_snapshot_sha256'],
            source_release=scope['source_release'],classification=classification,
            provenance=dict(classification=classification,source_manifest_sha256=source_manifest_sha256,
                objectives=[provenance]))
    return dict(inputs=[deduplicated[key] for key in sorted(deduplicated)],evaluations=evaluations)


def extract_calibration_inputs(root, manifest):
    from genosapiens.track2_report_evidence import check_evidence,scientific_root
    from genosapiens.track2_report_sources import check_source_records
    root=Path(root);check_evidence(root,manifest)
    science=scientific_root(root,manifest)
    versions,eligibility=check_source_records(root,manifest,scientific_root=science)
    result=retained_scopes(science,eligibility=eligibility,versions=versions,
        classification=manifest['classification'],source_manifest_sha256=json_identity(manifest)['sha256'])
    result['applicability']=calibration_applicability(result['evaluations'],eligibility=eligibility,
        classification=manifest['classification'],scientific_outcome=manifest['readiness']['assessment_status'])
    return result


def calibration_applicability(evaluations, *, eligibility, classification, scientific_outcome):
    """Keep scientific limits separate from operational preparation failures."""
    reasons=dict(complete=None,provisional='incomplete_required_conditions',
        unsupported_direction='unsupported_direction',unavailable='no_evaluable_conditions')
    scopes=[];covered=set()
    for evaluation in evaluations:
        require(evaluation['status'] in reasons, 'Calibration received an unknown scientific status.')
        applicable=evaluation['status']=='complete'
        require((evaluation['scope_sha256'] is not None)==applicable,
                'Calibration applicability differs from extracted preparation scope.')
        scopes.append(dict(evaluation,applicability='applicable' if applicable else 'not_applicable',
            reason=reasons[evaluation['status']]))
        covered.add(evaluation['objective_id'])
    for objective in eligibility.get('objectives',[]):
        if objective['objective_id'] in covered:continue
        scopes.append(dict(objective_id=objective['objective_id'],scenario_id=objective['scenario_id'],
            gene_id=objective['gene_id'],direction=objective['direction'],method=None,
            medication_entity_key=None,scope_sha256=None,applicability='not_applicable',reason='empty_condition_scope'))
    for hypothesis in eligibility.get('hypotheses',[]):
        if hypothesis['objective_ids']:continue
        scopes.append(dict(objective_id=None,scenario_id=hypothesis['scenario_id'],gene_id=None,
            direction=None,method=None,medication_entity_key=None,scope_sha256=None,
            applicability='not_applicable',reason='no_retained_intervention_objective'))
    if eligibility['state']!='available':
        scopes.append(dict(objective_id=None,scenario_id=None,gene_id=None,direction=None,method=None,
            medication_entity_key=None,scope_sha256=None,applicability='not_applicable',
            reason=eligibility.get('reason','no_retained_tissue_evaluation')))
    elif not scopes:
        scopes.append(dict(objective_id=None,scenario_id=None,gene_id=None,direction=None,method=None,
            medication_entity_key=None,scope_sha256=None,applicability='not_applicable',reason='empty_evaluated_scope'))
    return dict(schema='genosapiens.track2-calibration-applicability/v1',classification=classification,
        scientific_outcome=scientific_outcome,scopes=scopes,
        applicable_count=sum(row['applicability']=='applicable' for row in scopes))
