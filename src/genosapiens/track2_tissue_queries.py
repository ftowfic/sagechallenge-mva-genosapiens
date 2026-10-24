"""Prospective equal-cell medians and immutable genetic query refinement."""
from statistics import median

from genosapiens.track2_lincs import number
from genosapiens.track2_lincs_scoring import spearman
from genosapiens.track2_lincs_selection import select_cell_method
from genosapiens.track2_operator_inputs import require
from genosapiens.track2_tissue_records import identity, schema_name, validate_tissue_record


def pool_tissue_values(selection, genetic, features, policy):
    ids = [row['pr_gene_id'] for row in features]
    require(len(ids)==len(set(ids)) and all(row['measured_landmark']==1 for row in features),
        'Tissue queries require unique measured landmark features.')
    members = {r['cell_id'] for r in selection['membership'] if r['tissue_id']==selection['tissue_id']}
    contributors = {}
    for direction in ('oe', 'kd'):
        selected = selection['contributors'][direction]
        require(len({r['cell_id'] for r in selected})==len(selected),
            'A genetic cell cannot contribute more than once per direction.')
        contributors[direction] = []
        for source in selected:
            require(source['cell_id'] in members and source['gene_id']==selection['gene_id']
                and source['pert_type']==policy['lincs']['genetic_types'][direction],
                'A genetic contributor differs from the tissue query scope.')
            require(source['sig_id'] in genetic and isinstance(genetic[source['sig_id']], dict),
                'A selected genetic matrix column is missing or lacks feature identities.')
            values = {gene:number(genetic[source['sig_id']].get(gene)) for gene in ids}
            contributors[direction].append(dict(signature=dict(source), values=values))
    genes = []
    for feature in features:
        gene = feature['pr_gene_id']
        values = {d:[dict(cell_id=r['signature']['cell_id'], sig_id=r['signature']['sig_id'], value=r['values'][gene])
            for r in sources] for d,sources in contributors.items()}
        missing = [dict(direction=d,cell_id=r['cell_id'],sig_id=r['sig_id'],reason='missing_or_nonfinite_value')
            for d,rows in values.items() for r in rows if r['value'] is None]
        pooled = {d:median(r['value'] for r in rows) if rows and all(r['value'] is not None for r in rows) else None
            for d,rows in values.items()}
        genes.append(dict(gene_id=gene, gene_symbol=feature['symbol'], source_values=values,
            missing_contributors=missing, **pooled))
    return contributors, genes


def refine_tissue_query(selection, genetic, features, *, source_sha256, policy, scientific_code_sha256):
    contributors, genes = pool_tissue_values(selection, genetic, features, policy)
    method = select_cell_method({d:rows for d,rows in contributors.items() if rows}, policy)
    require(method==selection['method'], 'The prospective tissue method changed.')
    paired = method=='paired_oe_kd'
    minimum = policy['lincs']['minimum_absolute_modz']
    for row in genes:
        values = [row[d] for d, sources in contributors.items() if sources]
        reason = ('missing_genetic_directions' if method is None else
            'missing_contributor_value' if row['missing_contributors'] else
            'not_opposing' if paired and row['oe']*row['kd']>=0 else
            'below_magnitude' if any(abs(v)<minimum for v in values) else None)
        desired = -row['kd'] if method=='kd_only_reversal' and row['kd'] is not None else row['oe']
        row.update(selected=reason is None, reason=reason, desired=desired)
    selected = [row for row in genes if row['selected']]
    desired = [row['desired'] for row in selected]
    reason = 'missing_genetic_directions' if method is None else None
    if method and (len(selected)<policy['lincs']['minimum_features'] or spearman(desired,desired) is None or
            paired and spearman([r['kd'] for r in selected],[r['kd'] for r in selected]) is None):
        reason = 'inadequate_paired_query' if paired else 'inadequate_single_direction_query'
    query = dict(schema=schema_name('query'), gene_id=selection['gene_id'], tissue_id=selection['tissue_id'],
        method=method, required_modulation='increase', membership_sha256=identity(selection['membership']),
        source_sha256=source_sha256, policy_sha256=identity(policy), scientific_code_sha256=scientific_code_sha256,
        contributors=contributors, excluded_signatures=selection['excluded_signatures'], genes=genes,
        gene_count=len(selected), state='unavailable' if reason else 'ready', reason=reason)
    query['query_sha256'] = identity(query)
    return validate_tissue_record(query, 'query')
