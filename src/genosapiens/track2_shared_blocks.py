"""Revision-owned bounded matrix blocks and exact query projections."""
from pathlib import Path
import math

from genosapiens.track2_operator_inputs import require
from genosapiens.track2_tissue_records import identity
from genosapiens.track2_scientific_storage import file_identity, publish_json
from genosapiens.track2_scientific_completion import (
    SCHEMA, artifact_path, reference, check_reference, check_completion, publish_completion)

REPRESENTATION = 'genosapiens.track2-matrix-block/json-numbers-v1'


def plan_shared_blocks(jobs, queries, resources):
    """Plan only after every genetic query attempt reaches a terminal state."""
    require(all(q['state'] in {'completed','failed'} for q in queries.values()),
        'Drug extraction requires terminal genetic query attempts.')
    limit = resources['policy']['lincs']['column_block_size']
    require(type(limit) is int and 1 <= limit <= 256, 'Invalid shared block column bound.')
    groups = {}
    for job in sorted(jobs,key=lambda j:j['job_id']):
        outcome = queries[job['query_job_id']]
        if outcome['state'] != 'completed' or outcome['query']['state'] != 'ready':
            continue
        query = outcome['query']
        ids = tuple(r['gene_id'] for r in query['genes'] if r['selected'])
        signatures = tuple(job['drug_signature_ids'])
        require(len(signatures)==len(set(signatures)), 'Shared extraction has duplicate signatures.')
        key = (job['cell_id'], signatures)
        group = groups.setdefault(key, dict(rows=set(),consumers=[]))
        group['rows'].update(ids)
        group['consumers'].append(job['job_id'])
    reader = file_identity(Path(__file__).resolve().parents[2]/'scripts/read_track2_lincs.R')
    blocks = []
    for (cell,signatures),group in sorted(groups.items()):
        for start in range(0,len(signatures),limit):
            binding = dict(matrix=resources['descriptor']['matrix'],reader=reader,
                row_ids=sorted(group['rows']),signature_ids=list(signatures[start:start+limit]),representation=REPRESENTATION)
            blocks.append(dict(block_id=identity(binding),binding=binding,consumers=group['consumers']))
    return blocks


def project_block(block, gene_ids):
    positions = {gene:i for i,gene in enumerate(block['row_ids'])}
    require(len(positions)==len(block['row_ids']) and len(gene_ids)==len(set(gene_ids)) and
        all(gene in positions for gene in gene_ids), 'Shared block feature inventory differs.')
    require(len(block['signature_ids'])==len(set(block['signature_ids'])) and set(block['values'])==set(block['signature_ids']),
        'Shared block signature inventory differs.')
    result = {}
    for sig in block['signature_ids']:
        values = block['values'][sig]
        require(len(values)==len(positions), 'Shared block dimensions differ.')
        require(all(v is None or type(v) in (int,float) and math.isfinite(v) for v in values),
            'Shared block has invalid numeric values.')
        result[sig] = [values[positions[gene]] for gene in gene_ids]
    return result


def block_reference(root, block):
    key=block['block_id'];binding=block['binding']
    require(key==identity(binding), 'Shared block dependency identity changed.')
    return dict(block_id=key,binding=binding,completion=reference(root,
        'extractions/'+key+'/completion.json',role='shared-block-completion',schema=SCHEMA,dependencies=key))


def check_shared_reference(root, record, *, checked=None):
    key=record['block_id'];binding=record['binding']
    require(key==identity(binding) and record['completion']['path']=='extractions/'+key+'/completion.json',
        'A shared block reference has incompatible dependencies or location.')
    # Resolve paths even when this validation session already checked the bytes.
    completion=artifact_path(root,record['completion']['path'])
    cache_key=identity(record)
    if checked is None or cache_key not in checked:
        check_reference(root,record['completion'],role='shared-block-completion',schema=SCHEMA,dependencies=key)
        check_completion(completion.parent,stage='tissue-block',binding=binding)
        if checked is not None:checked.add(cache_key)
    return completion.parent/'block.json'


def publish_shared_block(root, block, resources):
    from genosapiens.track2_discovery import _extract
    from genosapiens.track2_discovery_outcome import record_failure
    from genosapiens.track2_lincs import number
    from genosapiens.track2_run_state import files
    root=Path(root);key=block['block_id'];binding=block['binding']
    require(key==identity(binding) and binding['representation']==REPRESENTATION and
        binding['matrix']==resources['descriptor']['matrix'] and
        binding['reader']==file_identity(Path(__file__).resolve().parents[2]/'scripts/read_track2_lincs.R'),
        'Shared extraction dependencies changed.')
    require(0<len(binding['signature_ids'])<=resources['policy']['lincs']['column_block_size']<=256,
        'Shared extraction exceeds its column bound.')
    final=artifact_path(root,'extractions/'+key)
    counters=resources.setdefault('execution_counters',{})
    if final.exists():
        record=block_reference(root,block)
        check_shared_reference(root,record)
        counters['shared_block_reuses']=counters.get('shared_block_reuses',0)+1
        return record
    pending=artifact_path(root,'extractions/'+key+'.pending')
    if pending.exists():
        closure=files(pending)
        retained=pending.with_name(key+'-incomplete-'+identity(closure))
        require(not retained.exists(), 'An incomplete shared block retention path already exists.')
        pending.rename(retained)
    pending.mkdir(parents=True,mode=0o700)
    working=dict(resources,extraction_cache={})
    try:
        values={}
        for part in _extract(binding['row_ids'],binding['signature_ids'],working,pending/'raw'):
            require(not set(values)&set(part) and len(values)+len(part)<=len(binding['signature_ids']),
                'Shared extraction exceeded its requested inventory.')
            values.update(part)
        require(list(values)==binding['signature_ids'],
            'Shared extraction returned a different signature inventory.')
        data=dict(row_ids=binding['row_ids'],signature_ids=binding['signature_ids'],
            values={sig:[number(v) for v in vector] for sig,vector in values.items()})
        project_block(data,binding['row_ids'])
        publish_json(pending/'block.json',data)
        publish_completion(pending,stage='tissue-block',binding=binding,result={})
        pending.rename(final)
        counters['shared_block_publications']=counters.get('shared_block_publications',0)+1
        return block_reference(root,block)
    except BaseException as error:
        record_failure(pending,stage='tissue-block',role='shared-extraction',error=error)
        raise
    finally:
        if working.get('reader_session') is not None:
            resources['reader_session']=working['reader_session']


def read_shared_block(root, record, gene_ids, *, checked=None):
    from genosapiens.track2_discovery import _read
    path=check_shared_reference(root,record,checked=checked)
    data=_read(path)
    require(data['row_ids']==record['binding']['row_ids'] and data['signature_ids']==record['binding']['signature_ids'],
        'Shared block identifiers differ from their binding.')
    return project_block(data,gene_ids)


def iter_shared_blocks(root, records, gene_ids, *, checked=None):
    """Decode only the reference currently requested by the scoring consumer."""
    for record in records:
        yield read_shared_block(root,record,gene_ids,checked=checked)
