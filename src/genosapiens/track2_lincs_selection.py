"""Prospective gene and cell selection from public metadata, before drug scoring."""
import json
import math

from genosapiens.track2_published_genotype import canonical, sha256


def resolve_perturbation_gene(record, gene_index):
    if record['pert_type'] not in {'trt_oe', 'trt_sh.cgs'}:
        return dict(gene_id=None, reason='not_a_selected_genetic_type')
    symbol = record['pert_iname'].strip().upper()
    approved = gene_index['approved'].get(symbol, set())
    matches = approved or gene_index['aliases'].get(symbol, set())
    return dict(gene_id=next(iter(matches)) if len(matches)==1 else None,
        reason=('approved_symbol' if approved else 'source_alias') if len(matches)==1
        else 'ambiguous_gene_alias' if matches else 'unmapped_gene_name')


def select_cell_method(selected_signatures, policy):
    directions = set(selected_signatures) & {'oe', 'kd'}
    key = 'oe+kd' if directions == {'oe', 'kd'} else next(iter(directions), '')
    return policy['lincs']['method_selection'].get(key)


def _rows(connection, sql, parameters):
    cursor = connection.execute(sql, parameters)
    fields = [field[0] for field in cursor.description]
    return [dict(zip(fields, row)) for row in cursor]


def _qc(row, thresholds):
    return all(isinstance(row.get(k), (float,int)) and math.isfinite(row[k]) and row[k]>=limit
        for k, limit in thresholds.items())


def select_genetic_signatures(connection, gene_id, policy):
    rows = _rows(connection, 'SELECT * FROM lincs_signatures WHERE gene_id=? AND pert_type IN (?,?)',
        (gene_id,'trt_oe','trt_sh.cgs'))
    results = []
    for cell in sorted({row['cell_id'] for row in rows}):
        selected, exclusions = {}, []
        for direction, kind in policy['lincs']['genetic_types'].items():
            local = [row for row in rows if row['cell_id']==cell and row['pert_type']==kind]
            passing = sorted((row for row in local if _qc(row, policy['lincs']['genetic_qc'])),
                key=lambda r:(-r['tas'],-r['distil_cc_q75'],-r['distil_nsample'],r['sig_id']))
            if passing:
                selected[direction] = {k:v for k,v in passing[0].items() if k!='source_json'}
            exclusions.extend(dict(sig_id=row['sig_id'],direction=direction,
                reason='lower_genetic_priority' if _qc(row,policy['lincs']['genetic_qc']) else 'failed_or_missing_genetic_QC')
                for row in local if not passing or row['sig_id']!=passing[0]['sig_id'])
            if not local:
                exclusions.append(dict(sig_id=None,direction=direction,reason='missing_genetic_direction'))
        annotation = connection.execute('SELECT annotations_json FROM lincs_cells WHERE cell_id=?',(cell,)).fetchone()
        results.append(dict(gene_id=gene_id,cell_id=cell,selected=selected,method=select_cell_method(selected,policy),
            exclusions=sorted(exclusions,key=lambda r:(r['direction'],r['sig_id'] or '')),
            cell_metadata=json.loads(annotation[0]) if annotation else None))
    return results


def select_drug_conditions(connection, cell_id, medication_catalog, policy):
    # The catalog's provider identity is independent of Open Targets admission.
    rows = _rows(connection, '''SELECT s.*, c.ingredient_key, m.data AS medication
        FROM lincs_signatures s LEFT JOIN compound_identity c ON c.provider='lincs' AND c.compound_id=s.pert_id
        LEFT JOIN medications m ON m.ingredient_key=c.ingredient_key
        WHERE s.cell_id=? AND s.pert_type='trt_cp' ORDER BY s.sig_id''', (cell_id,))
    conditions, exclusions = [], []
    for row in rows:
        reason = ('unresolved_medication_identity' if row['ingredient_key'] is None else
            'unsupported_drug_time' if row['pert_time_unit']!='h' or row['pert_time'] not in policy['lincs']['drug_times_hours'] else
            'failed_or_missing_drug_QC' if not _qc(row,policy['lincs']['drug_qc']) else None)
        if reason:
            exclusions.append(dict(sig_id=row['sig_id'],reason=reason))
        else:
            conditions.append({k:v for k,v in row.items() if k not in {'source_json','medication'}} | {
                'medication_entity_key':row['ingredient_key'],'medication':json.loads(row['medication'])})
    return dict(conditions=conditions, exclusions=exclusions)


def discovery_job_key(selection, resource_bindings, policy_sha256):
    from genosapiens.track2_discovery import scientific_identity
    return 'job-'+sha256(canonical(dict(gene_id=selection['gene_id'],direction='increase',
        cell_id=selection['cell_id'],method=selection['method'],selected=selection['selected'],
        resources=resource_bindings,policy_sha256=policy_sha256,scientific_code_sha256=scientific_identity())))


def select_tissue_genetic_signatures(connection, gene_id, policy, membership):
    """Group the existing prospective per-cell choices without additional weight."""
    from genosapiens.track2_operator_inputs import require
    from genosapiens.track2_tissue_records import validate_tissue_record
    cells = {}
    for row in membership:
        validate_tissue_record(row, 'membership')
        require(row['cell_id'] not in cells, 'Tissue membership has duplicate cells.')
        cells[row['cell_id']] = row
    selections = select_genetic_signatures(connection, gene_id, policy)
    require(all(row['cell_id'] in cells for row in selections), 'A genetic cell lacks a membership record.')
    results = []
    for tissue in sorted({row['tissue_id'] for row in membership if row['tissue_id'] is not None}):
        members = sorted((row for row in membership if row['tissue_id']==tissue), key=lambda r:r['cell_id'])
        local = [row for row in selections if cells[row['cell_id']]['tissue_id']==tissue]
        contributors = {direction:[row['selected'][direction] for row in local if direction in row['selected']]
            for direction in ('oe', 'kd')}
        excluded = []
        for row in local:
            for exclusion in row['exclusions']:
                source = _rows(connection, 'SELECT * FROM lincs_signatures WHERE sig_id=?',
                    (exclusion['sig_id'],)) if exclusion['sig_id'] else []
                excluded.append(dict(exclusion, cell_id=row['cell_id'],
                    source={k:v for k,v in source[0].items() if k!='source_json'} if source else None))
        results.append(dict(gene_id=gene_id, tissue_id=tissue, membership=members,
            contributors=contributors, excluded_signatures=excluded,
            method=select_cell_method({d:rows for d,rows in contributors.items() if rows}, policy)))
    unresolved = [dict(membership=row, genetic_selections=[s for s in selections if s['cell_id']==row['cell_id']])
        for row in sorted(membership, key=lambda r:r['cell_id']) if row['tissue_id'] is None]
    return dict(tissues=results, unresolved_cells=unresolved)
