"""Independent calibration journals with preserved attempts and entry reuse."""
import json
import os
import time
import uuid
from pathlib import Path

from genosapiens.track2_operator_inputs import require
from genosapiens.track2_calibration_bank import prepare_bank,scientific_scope,entry_scope,_path,_write,admit_entry
from genosapiens.track2_frozen_calibration import identity,canonical


def prepare_report_calibration(operation, inputs, *, classification, source_manifest_sha256,
                               batch_size=8, on_entry=None, retained_calibration=()):
    operation=Path(operation).absolute()
    _path(operation,'request.json');operation.mkdir(mode=0o700,parents=True,exist_ok=True)
    require(operation.stat().st_uid==os.geteuid() and not operation.stat().st_mode & 0o077,
            'Calibration state requires a private owned directory.')
    scopes=[];provenance={};scope_by_entry={}
    for source in inputs:
        require(source['classification']==classification and
            source['provenance']['source_manifest_sha256']==source_manifest_sha256,
            'Calibration inputs have another source or classification.')
        scope=scientific_scope(source['queries'],source['conditions'],direction=source['direction'],
            fda_snapshot_sha256=source['fda_snapshot_sha256'],eligible_ingredients=source['eligible_ingredients'],
            source_release=source['source_release'])
        scopes.append(identity(scope))
        for pair in sorted({(r['medication_entity_key'],r['method']) for r in scope['conditions'] if r['state']=='valid'}):
            entry=entry_scope(scope,*pair);key=identity(entry);scope_by_entry[key]=entry
            record=provenance.setdefault(key,dict(classification=classification,
                source_manifest_sha256=source_manifest_sha256,objectives=[]))
            for objective in source['provenance'].get('objectives',[]):
                if objective not in record['objectives']:record['objectives'].append(objective)
    request=dict(schema='genosapiens.track2-report-calibration-request/v1',classification=classification,
        source_manifest_sha256=source_manifest_sha256,inputs_sha256=[identity(source) for source in inputs],scopes=scopes,
        retained_calibration=[record['expected_manifest_sha256'] for record in retained_calibration])
    request_path=operation/'request.json'
    if request_path.exists():
        require(request_path.read_bytes()==canonical(request), 'Calibration dependencies changed. Preserve the original operation.')
    else:_write(request_path,request)
    expected_root=_path(operation,'expected');expected_root.mkdir(mode=0o700,exist_ok=True)
    attempts=_path(operation,'attempts');attempts.mkdir(mode=0o700,exist_ok=True)
    attempt=attempts/uuid.uuid4().hex;attempt.mkdir(mode=0o700)
    _write(attempt/'started.json',dict(status='incomplete',request_sha256=identity(request),classification=classification))
    expected={}
    for path in sorted(expected_root.glob('*.json')):
        binding=json.loads(_path(expected_root,path.name).read_bytes())
        require(path.name==binding['entry_key']+'.json', 'A calibration expectation has another entry identity.')
        expected[binding['entry_key']]=binding
    bank=_path(operation,'bank')
    initial={key for key in expected if _path(bank/'entries',key).exists()}
    completed_keys=set(initial);imported=set();counts=dict(query_permutations=0,aggregate_iterations=0)
    entries={};start=time.monotonic()
    def ready(binding):
        key=binding['entry_key'];path=expected_root/(key+'.json')
        if path.exists():
            require(path.read_bytes()==canonical(binding), 'A completed calibration identity changed.')
        else:_write(path,binding)
        expected[key]=binding
    def completed(binding):
        completed_keys.add(binding['entry_key'])
        if on_entry is not None:on_entry(binding)
    try:
        from genosapiens.track2_report_calibration import check_calibration
        from genosapiens.track2_calibration_bank import publish_entry
        for retained in retained_calibration:
            component=check_calibration(retained['root'],expected_manifest_sha256=retained['expected_manifest_sha256'])
            require(classification!='public' or component['classification']=='public',
                    'Controlled calibration cannot enter a public operation.')
            for original in component['entries']:
                key=original['entry_key']
                if key not in scope_by_entry or key in initial or key in imported:continue
                background=admit_entry(retained['root'],expected_entry_key=key,
                    expected_distribution_sha256=original['distribution_sha256'])
                content={k:original[k] for k in ('entry_key','scope_sha256','distribution_sha256')}
                binding=dict(content,classification=classification,
                    provenance_sha256=identity(dict(content,classification=classification,source=provenance[key])))
                ready(binding)
                published=publish_entry(bank,scope_by_entry[key],list(background.draws),
                    provenance=provenance[key],classification=classification)
                require(published==binding, 'Retained calibration changed during independent publication.')
                imported.add(key);completed(binding)
        for source in inputs:
            resource=prepare_bank(bank,source,admitted_entries=list(expected.values()),counters=counts,
                batch_size=batch_size,on_entry=completed,on_ready=ready,provenance_by_entry=provenance)
            for binding in resource['entries']:
                prior=entries.setdefault(binding['entry_key'],binding)
                require(prior==binding, 'Calibration scopes have conflicting provenance.')
        resource=dict(schema='genosapiens.track2-calibration-resource/v1',scope_sha256=identity(scopes),
            classification=classification,entries=[entries[key] for key in sorted(entries)])
        for binding in resource['entries']:
            admit_entry(bank,expected_entry_key=binding['entry_key'],expected_distribution_sha256=binding['distribution_sha256'])
        path=operation/'resource.json'
        if path.exists():require(path.read_bytes()==canonical(resource), 'The retained calibration resource changed.')
        else:_write(path,resource)
    except BaseException as error:
        _write(attempt/'outcome.json',dict(status='incomplete',error_type=type(error).__name__,
            elapsed_seconds=time.monotonic()-start,operation_counts=counts,
            entries_new=len(completed_keys-initial-imported),entries_reused=len((initial|imported) & set(entries)),
            entries_imported=len(imported)))
        raise
    result=dict(status='complete',resource_contract_sha256=identity(resource),
        elapsed_seconds=time.monotonic()-start,operation_counts=counts,
        entries_new=len(completed_keys-initial-imported),entries_reused=len((initial|imported) & set(entries)),
        entries_imported=len(imported),entry_count=len(entries))
    _write(attempt/'outcome.json',result)
    return dict(bank=bank,resource_contract=path,expected_resource_sha256=identity(resource),operations=result)
