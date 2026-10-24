"""Dependent query and drug-cell jobs with bounded matrix extraction."""
from pathlib import Path
import shutil

from genosapiens.track2_discovery_resources import load_discovery_resources, resource_bindings
from genosapiens.track2_operator_inputs import require
from genosapiens.track2_scientific_completion import check_completion, publish_completion
from genosapiens.track2_scientific_storage import publish_json, publish_chunks, json_chunks
from genosapiens.track2_tissue_records import identity, validate_tissue_record
from genosapiens.track2_tissue_admission import TissueAdmission


def _admission(resources):
    saved = resources.get('tissue_admission')
    if saved is not None:
        saved.check(resources['descriptor'], resources['policy'])
        return saved, None
    connection = load_discovery_resources(resources['descriptor'], resources['reference_root'])
    try:
        context = TissueAdmission(resources['descriptor'], resources['policy'], connection,
            counters=resources.setdefault('execution_counters', {}))
        return context, connection
    except BaseException:
        connection.close()
        raise


def admit_tissue_drug_job(job, query, resources, admission):
    admission.check(resources['descriptor'],resources['policy'])
    require(query['source_sha256']==resources['descriptor']['matrix']['content_sha256'],
        'The tissue query matrix source changed.')
    require(query['policy_sha256']==identity(resources['policy']), 'The tissue query policy changed.')
    conditions=admission.drugs(job['cell_id'])
    require(conditions==dict(conditions=job['conditions'],exclusions=job['drug_exclusions']),
        'Prospective drug conditions changed.')
    require(job['drug_signature_ids']==[r['sig_id'] for r in conditions['conditions']],
        'Prospective drug signature order changed.')


def tissue_drug_binding(job,query,resources):
    from genosapiens.track2_discovery import scientific_identity
    return dict(selection_sha256=identity(job),query_sha256=query['query_sha256'],
        resources=resource_bindings(resources['descriptor']),code=scientific_identity(),runtime=resources['runtime'])


def _begin(output, stage, binding, *, force_restart=False):
    from genosapiens.track2_discovery import _read
    from genosapiens.track2_run_state import files, atomic_record
    if not force_restart and (output/'completion.json').exists():
        try:
            check_completion(output, stage=stage, binding=binding)
            return _read(output/'result.json')
        except (OSError, ValueError, KeyError, TypeError):
            pass
    if output.exists():
        # Preserve failed and incompatible evidence before replacing this attempt.
        closure=files(output)
        archive=output.parent/(output.name+'-retained-'+identity(closure))
        if not archive.exists(): shutil.copytree(output,archive)
        else: require(files(archive)==closure, 'Retained scientific job evidence changed.')
        shutil.rmtree(output)
    output.mkdir(parents=True, mode=0o700)
    atomic_record(output/'progress.json',dict(status='running',binding=binding))
    return None


def execute_tissue_query(selection, resources, output):
    from genosapiens.track2_discovery import _extract, scientific_identity, _read
    from genosapiens.track2_tissue_queries import refine_tissue_query
    from genosapiens.track2_discovery_outcome import record_failure
    output=Path(output); connection=None
    binding=dict(selection_sha256=identity(selection), resources=resource_bindings(resources['descriptor']),
        code=scientific_identity(), runtime=resources['runtime'])
    try:
        if saved:=_begin(output,'tissue-query',binding): return saved
        publish_json(output/'selection.json',selection)
        admission,connection=_admission(resources)
        features=admission.features()
        ids=[r['pr_gene_id'] for r in features]
        signatures=sorted(r['sig_id'] for rows in selection['contributors'].values() for r in rows)
        genetic={}
        if signatures:
            for block in _extract(ids,signatures,resources,output/'genetic'):
                genetic.update({sig:dict(zip(ids,values)) for sig,values in block.items()})
        query=refine_tissue_query(selection,genetic,features,
            source_sha256=resources['descriptor']['matrix']['content_sha256'], policy=resources['policy'],
            scientific_code_sha256=scientific_identity())
        publish_json(output/'query.json',query)
        result=dict(schema='genosapiens.track2-tissue-query-job-result/v1',query_job_id=selection['query_job_id'],
            gene_id=query['gene_id'],tissue_id=query['tissue_id'],state='completed',query=query)
        publish_json(output/'result.json',result)
        publish_completion(output,stage='tissue-query',binding=binding,result={})
        return _read(output/'result.json')
    except BaseException as error:
        record_failure(output,stage='tissue-query',role='query-artifacts',error=error)
        raise
    finally:
        if connection is not None: connection.close()


def execute_tissue_drug_job(job, query, resources, output):
    from genosapiens.track2_discovery import _extract, scientific_identity, _read
    from genosapiens.track2_tissue_scoring import prepare_tissue_query, score_prepared_tissue_condition, summarize_tissue_conditions
    from genosapiens.track2_figure_stream import condition_rows
    from genosapiens.track2_discovery_outcome import record_failure
    output=Path(output); connection=None
    validate_tissue_record(query,'query')
    require(query['gene_id']==job['gene_id'] and query['tissue_id']==job['tissue_id'] and query['method']==job['method'],
        'A drug-cell job refers to a different tissue query.')
    binding=tissue_drug_binding(job,query,resources)
    try:
        admission,connection=_admission(resources)
        admit_tissue_drug_job(job,query,resources,admission)
        shared=resources.get('shared_blocks_by_job',{}).get(job['job_id'])
        if shared is not None:
            check_job_blocks(resources['revision_root'],shared,query,job,
                checked=resources.setdefault('checked_shared_blocks',set()))
        if saved:=_begin(output,'tissue-drug',binding): return saved
        publish_json(output/'selection.json',job); publish_json(output/'query.json',query)
        metadata={r['sig_id']:r for r in job['conditions']}
        counters=resources.setdefault('execution_counters', {})
        context=prepare_tissue_query(query, counters=counters)
        ids=list(context.gene_ids)
        def drug_blocks():
            if shared is None:
                yield from _extract(ids,job['drug_signature_ids'],resources,output/'drugs')
            else:
                from genosapiens.track2_shared_blocks import iter_shared_blocks
                yield from iter_shared_blocks(resources['revision_root'],shared,ids,
                    checked=resources['checked_shared_blocks'])
        def scored():
            if query['state']!='ready':
                for row in job['conditions']:
                    yield from json_chunks(score_prepared_tissue_condition(context,row,{},counters=counters))
            else:
                for block in drug_blocks():
                    for sig, values in block.items():
                        require(len(values)==len(ids), 'The extracted condition features are misaligned.')
                        yield from json_chunks(score_prepared_tissue_condition(context,metadata[sig],dict(zip(ids,values)),counters=counters))
                    del block
        publish_chunks(output/'conditions.jsonl',scored())
        medications=summarize_tissue_conditions(condition_rows(output/'conditions.jsonl'))
        result=dict(schema='genosapiens.track2-tissue-drug-result/v1',job_id=job['job_id'],query_job_id=job['query_job_id'],
            gene_id=job['gene_id'],tissue_id=job['tissue_id'],cell_id=job['cell_id'],method=job['method'],
            query_sha256=query['query_sha256'],gene_count=query['gene_count'],
            state='completed' if query['state']=='ready' else 'unavailable',reason=query['reason'],
            medications=medications,exclusions=job['drug_exclusions'],conditions_path='conditions.jsonl',
            query_path='query.json',source_sha256=query['source_sha256'])
        if shared is not None:result['shared_blocks']=shared
        publish_json(output/'result.json',result)
        publish_completion(output,stage='tissue-drug',binding=binding,result={})
        return _read(output/'result.json')
    except BaseException as error:
        record_failure(output,stage='tissue-drug',role='job-artifacts',error=error)
        raise
    finally:
        if connection is not None: connection.close()


def check_job_blocks(root, records, query, job, *, checked=None):
    from genosapiens.track2_shared_blocks import check_shared_reference, REPRESENTATION
    signatures=[]
    ids={r['gene_id'] for r in query['genes'] if r['selected']}
    for record in records:
        binding=record['binding']
        require(binding['matrix']['content_sha256']==query['source_sha256'] and
            binding['representation']==REPRESENTATION and ids<=set(binding['row_ids']),
            'A shared block differs from its consumer query.')
        check_shared_reference(root,record,checked=checked)
        signatures.extend(binding['signature_ids'])
    require(signatures==(job['drug_signature_ids'] if query['state']=='ready' else []),
        'The shared blocks differ from the planned signature inventory.')


def validate_tissue_job(path, *, counters=None, revision_root=None, shared_checked=None):
    from genosapiens.track2_discovery import _read
    from genosapiens.track2_figure_stream import condition_rows
    from genosapiens.track2_tissue_scoring import summarize_tissue_conditions
    path=Path(path);receipt=_read(path/'completion.json')
    check_completion(path,stage='tissue-drug',binding=receipt['binding'])
    result=_read(path/'result.json');query=_read(path/'query.json');selected=_read(path/'selection.json')
    validate_tissue_record(query,'query')
    require(receipt['binding']['selection_sha256']==identity(selected)
        and receipt['binding']['query_sha256']==query['query_sha256']==result['query_sha256'],
        'The tissue job dependencies changed.')
    if 'shared_blocks' in result:
        check_job_blocks(revision_root or path.parent.parent,result['shared_blocks'],query,selected,checked=shared_checked)
    metadata={r['sig_id']:r for r in selected['conditions']}
    require(len(metadata)==len(selected['conditions']) and
        sorted(metadata)==sorted(selected['drug_signature_ids']), 'The planned condition inventory changed.')
    fields=('medication_entity_key','sig_id','pert_id','pert_iname','cell_id',
        'pert_dose','pert_dose_unit','pert_time','pert_time_unit')
    seen=set()
    ids=[r['gene_id'] for r in query['genes'] if r['selected']]
    def checked_conditions():
        if counters is not None:
            counters['condition_semantic_passes']=counters.get('condition_semantic_passes',0)+1
        for row in condition_rows(path/'conditions.jsonl'):
            sig=row['sig_id']
            require(sig in metadata and sig not in seen, 'A tissue job omitted, duplicated, or added admitted conditions.')
            seen.add(sig)
            require(all(row[k]==query[k] for k in ('query_sha256','gene_id','tissue_id','method'))
                and row['cell_id']==selected['cell_id'] and row['gene_ids']==ids,
                'A condition has a different query, cell, or feature binding.')
            require(all(row[k]==metadata[sig][k] for k in fields) and
                row['quality']=={k:metadata[sig][k] for k in ('tas','distil_cc_q75','distil_nsample')},
                'A retained condition differs from prospective admission.')
            yield row
        require(seen==set(metadata), 'A tissue job omitted or added admitted conditions.')
    require(result['medications']==summarize_tissue_conditions(checked_conditions()),
        'Tissue job summaries differ from retained conditions.')
    return result
