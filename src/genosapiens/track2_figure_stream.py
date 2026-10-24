"""One native figure page at a time, with disk-backed condition ordering."""
import itertools
import json
from pathlib import Path
import sqlite3
import tempfile

from genosapiens.track2_lincs_scoring import POLICY
from genosapiens.track2_operator_inputs import require


def condition_rows(path):
    with Path(path).open() as stream:
        for line in stream:
            yield json.loads(line)


def figure_pages(query, result, conditions_path):
    """Bound memory by one query, five medications, and one 24-by-eight page.

    SQLite uses a 2 MiB cache and disk-backed sorting. Temporary files stay
    private and the worker's ephemeral storage owns their lifecycle.
    """
    if query['schema']=='genosapiens.track2-tissue-query/v1':
        from genosapiens.track2_tissue_figures import tissue_condition_pages
        yield from tissue_condition_pages(query,result,conditions_path)
        return
    selected=[g for g in query['genes'] if g['selected']]
    displayed=list(itertools.islice((r for r in result['medications'] if r['state']=='pass'),POLICY['top_n_display']))
    order={row['medication_entity_key']:i for i,row in enumerate(displayed)}
    limit=max([1.]+[abs(g[k]) for g in selected for k in ('oe','kd') if g[k] is not None])
    with tempfile.TemporaryDirectory(prefix='track2-figure-') as temporary:
        connection=sqlite3.connect(Path(temporary)/'conditions.sqlite')
        try:
            connection.execute('PRAGMA cache_size=-2048')
            (Path(temporary)/'conditions.sqlite').chmod(0o600)
            connection.execute('PRAGMA temp_store=FILE')
            connection.execute('CREATE TABLE conditions (position INTEGER, time REAL, identity TEXT, data TEXT)')
            count=0
            for row in condition_rows(conditions_path):
                if row['medication_entity_key'] not in order:continue
                require(len(row['values'])==len(selected),'Figure gene and condition values differ.')
                limit=max(limit,max((abs(v) for v in row['values'] if v is not None),default=1.))
                connection.execute('INSERT INTO conditions VALUES (?,?,?,?)',
                    (order[row['medication_entity_key']],float(row['pert_time']),row['sig_id'],json.dumps(row)))
                count+=1
            connection.commit()
            common=dict(schema='genosapiens.track2-native-figure-page/v1',cell_id=query['cell_id'],method=query['method'],gene_id=query['gene_id'],
                query_sha256=query['query_sha256'],gene_ids=[g['gene_id'] for g in selected],
                top_n_display=POLICY['top_n_display'],palette=['#2166AC','#F7F7F7','#D97706'],
                missing='gray labeled NA',expression_limit=limit,clipping=False)
            caption=('The figure shows '+query['method']+' in '+query['cell_id']+'. '
                'Scores apply only to this gene, cell, and method. Expression matching does not establish rescue. '
                + ('Both genetic directions contribute to the query. Selected opposition does not provide independent validation.'
                    if query['method']=='paired_oe_kd' else
                    'The opposite genetic direction is missing. Opposite correlation and joint score remain unavailable.'))
            # The full source query remains separately authenticated. Rendering
            # needs its feature records and state, not duplicated source metadata.
            projection={k:query[k] for k in ('gene_id','cell_id','method','genes','state','reason','absent_direction')}
            yield dict(common,id=query['cell_id']+'-selection',kind='selection',query=projection,
                caption=caption,alt=query['cell_id']+' selection and retained feature counts.')
            for gene_start in range(0,len(selected),24):
                cursor=connection.execute('SELECT data FROM conditions ORDER BY position,time,identity')
                for condition_start in range(0,max(1,count),8):
                    rows=[json.loads(row[0]) for row in cursor.fetchmany(8)]
                    subset=[dict(row,values=row['values'][gene_start:gene_start+24],
                        gene_ids=row['gene_ids'][gene_start:gene_start+24]) for row in rows]
                    yield dict(common,id=f"{query['cell_id']}-heatmap-{gene_start//24+1}-{condition_start//8+1}",
                        kind='heatmap',genes=selected[gene_start:gene_start+24],gene_offset=0,
                        conditions=subset,condition_ids=[row['sig_id'] for row in rows],caption=caption,
                        alt=query['cell_id']+' gene and condition responses with desired directions.')
            yield dict(common,id=query['cell_id']+'-ranking',kind='ranking',medications=displayed,
                n_tested=len(result['medications']),n_pass=sum(r['state']=='pass' for r in result['medications']),
                n_fail=sum(r['state']=='fail' for r in result['medications']),
                n_unassessed=sum(r['state']=='unassessed' for r in result['medications']),
                caption=caption,alt=query['cell_id']+' local expression ranks and signed scores.')
        finally:connection.close()


def check_figure_pages(folder, query, result, conditions_path, *, context=None):
    from genosapiens.track2_checked_pages import CheckedPageStreams
    if context is None:
        owned=CheckedPageStreams()
        try:return owned.check(folder,query,result,conditions_path)
        finally:owned.close()
    require(isinstance(context,CheckedPageStreams),'Page checks require an invocation-owned context.')
    return context.check(folder,query,result,conditions_path)
