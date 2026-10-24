"""Immutable calibration entries with separately bound source provenance."""
from bisect import bisect_left, bisect_right
from dataclasses import dataclass
import hashlib
import json
import os
from pathlib import Path
import shutil
import uuid

from genosapiens.track2_derived_strength import BACKGROUND_SCHEMA, coefficient, policy, validate_contract
from genosapiens.track2_frozen_calibration import canonical, identity, query_content

CONDITION_FIELDS = ('gene_id','tissue_id','method','medication_entity_key','sig_id','pert_id',
    'cell_id','pert_dose','pert_dose_unit','pert_time','pert_time_unit','quality','gene_ids','values',
    'r_oe','r_kd','state','reason')


def scientific_scope(queries, conditions, *, direction, fda_snapshot_sha256, eligible_ingredients, source_release):
    """Bind scientific dependencies without patient, run, or objective labels."""
    _digest(fda_snapshot_sha256)
    if direction!='increase':raise ValueError('Unsupported calibration direction.')
    query_map={};contents={}
    for query in queries:
        content=query_content(query);key=identity(content)
        previous=query_map.setdefault(query['query_sha256'],key)
        if previous!=key:raise ValueError('A query identity has conflicting content.')
        contents[key]=content
    eligible=sorted(set(eligible_ingredients));rows={}
    for row in conditions:
        if row['medication_entity_key'] not in eligible:continue
        value={field:row[field] for field in CONDITION_FIELDS}
        value['query_key']=query_map[row['query_sha256']]
        key=(value['query_key'],value['sig_id'])
        if key in rows and rows[key]!=value:raise ValueError('Repeated scientific observations disagree.')
        rows[key]=value
    root=Path(__file__).resolve().parent
    code={name:hashlib.sha256((root/name).read_bytes()).hexdigest() for name in
        ('track2_derived_strength.py','track2_frozen_calibration.py','track2_lincs_scoring.py')}
    return dict(schema='genosapiens.track2-calibration-scope/v1',direction=direction,
        fda_snapshot_sha256=fda_snapshot_sha256,eligible_ingredients=eligible,source_release=source_release,
        queries={key:contents[key] for key in sorted(contents)},conditions=[rows[key] for key in sorted(rows)],
        policy=policy(),scoring_code=code)


def entry_scope(scope, ingredient, method):
    if ingredient not in scope['eligible_ingredients']:raise ValueError('Ingredient is outside the frozen FDA scope.')
    conditions=[r for r in scope['conditions'] if r['medication_entity_key']==ingredient and r['method']==method]
    if not any(r['state']=='valid' for r in conditions):raise ValueError('An entry requires evaluable eligible conditions.')
    keys={r['query_key'] for r in conditions}
    return dict(scope,queries={k:scope['queries'][k] for k in sorted(keys)},conditions=conditions,
        medication_entity_key=ingredient,method=method)


def _digest(value):
    if not isinstance(value,str) or len(value)!=64 or any(c not in '0123456789abcdef' for c in value):
        raise ValueError('An expected SHA-256 identity is required.')
    return value


def _path(root, name):
    root=Path(root).absolute();path=root/name
    if root.resolve()!=root or path.resolve()!=path or path.is_symlink():
        raise ValueError('Calibration links are unsupported.')
    return path


def _write(path, value):
    with path.open('xb') as stream:
        stream.write(canonical(value));stream.flush();os.fsync(stream.fileno())
    path.chmod(0o400)


def publish_entry(bank, scope, draws, *, provenance, classification):
    """Keep failed private attempts and publish only complete immutable entries."""
    if classification not in {'public','controlled'}:raise ValueError('Unknown calibration classification.')
    if provenance.get('classification')!=classification:raise ValueError('Source classification differs.')
    _digest(provenance.get('source_manifest_sha256'))
    bank=Path(bank).absolute();key=identity(scope)
    _path(bank,'entries');bank.mkdir(mode=0o700,parents=True,exist_ok=True)
    if bank.stat().st_uid!=os.geteuid() or bank.stat().st_mode & 0o077:
        raise ValueError('Calibration publication requires a private owned bank.')
    for name in ('attempts','entries','provenance'):
        _path(bank,name).mkdir(mode=0o700,exist_ok=True)
    attempt=bank/'attempts'/uuid.uuid4().hex;attempt.mkdir(mode=0o700)
    _write(attempt/'started.json',dict(status='incomplete',entry_key=key))
    record=dict(schema=BACKGROUND_SCHEMA,method=scope['method'],direction=scope['direction'],
        entry_key=key,draw_count=10000,draws=draws,
        numerical_contract=policy()['calibration']['numerical_contract'])
    validate_contract(record,'frozen-background')
    _write(attempt/'distribution.json',record)
    _write(attempt/'scope.json',scope)
    binding=dict(entry_key=key,distribution_sha256=identity(record),scope_sha256=identity(scope))
    # A scientific distribution can be shared. Every consuming provenance remains separate.
    provenance_record=dict(binding,classification=classification,source=provenance)
    _write(attempt/'provenance.json',provenance_record)
    destination=_path(bank,'entries')/key
    payload=attempt/'entry';payload.mkdir(mode=0o700)
    for name in ('distribution.json','scope.json'):
        shutil.copyfile(attempt/name,payload/name)
        (payload/name).chmod(0o400)
    _write(payload/'binding.json',binding)
    provenance_path=_path(bank,'provenance')/(identity(provenance_record)+'.json')
    try:_write(provenance_path,provenance_record)
    except FileExistsError:
        if provenance_path.read_bytes()!=canonical(provenance_record):raise ValueError('Existing provenance changed.')
    # Atomic rename refuses an existing nonempty directory. Never replace an entry.
    try:
        if destination.exists() or destination.is_symlink():
            raise FileExistsError('The calibration entry destination already exists.')
        os.rename(payload,destination)
    except OSError:
        if not destination.is_dir():raise
        if {p.name for p in destination.iterdir()}!={'distribution.json','scope.json','binding.json'}:
            raise ValueError('An existing calibration entry has an incomplete inventory.')
        for name in ('distribution.json','scope.json','binding.json'):
            if _path(destination,name).read_bytes()!=_path(payload,name).read_bytes():
                raise ValueError('An existing immutable calibration entry differs.')
    _write(attempt/'complete.json',dict(status='complete',**binding,provenance_sha256=identity(provenance_record)))
    return dict(binding,provenance_sha256=identity(provenance_record),classification=classification)


@dataclass(frozen=True)
class AdmittedBackground:
    entry_key: str
    distribution_sha256: str
    method: str
    draws: tuple


def admit_entry(bank, *, expected_entry_key, expected_distribution_sha256):
    """Expect identities from an independently admitted resource contract."""
    _digest(expected_entry_key);_digest(expected_distribution_sha256)
    root=_path(bank,'entries')/expected_entry_key
    if {p.name for p in root.iterdir()}!={'distribution.json','scope.json','binding.json'}:
        raise ValueError('Calibration entry inventory differs.')
    raw=_path(root,'distribution.json').read_bytes()
    if hashlib.sha256(raw).hexdigest()!=expected_distribution_sha256:
        raise ValueError('Calibration distribution differs from its admitted identity.')
    record=json.loads(raw);validate_contract(record,'frozen-background')
    scope=json.loads(_path(root,'scope.json').read_bytes())
    binding=json.loads(_path(root,'binding.json').read_bytes())
    if identity(scope)!=expected_entry_key or record['entry_key']!=expected_entry_key:
        raise ValueError('Calibration scientific scope differs.')
    if binding!=dict(entry_key=expected_entry_key,distribution_sha256=expected_distribution_sha256,scope_sha256=expected_entry_key):
        raise ValueError('Calibration entry binding differs.')
    if scope['method']!=record['method'] or scope['direction']!=record['direction']:
        raise ValueError('Calibration method or direction differs.')
    return AdmittedBackground(expected_entry_key,expected_distribution_sha256,record['method'],tuple(record['draws']))


def percentile(background, observed):
    if not isinstance(background,AdmittedBackground):raise ValueError('Percentile lookup requires admitted evidence.')
    coefficient(observed)
    if observed is None:raise ValueError('Percentile lookup requires an observed strength.')
    draws=background.draws
    lower=bisect_left(draws,observed);equal=bisect_right(draws,observed)-lower
    numerator=2*lower+equal;denominator=2*len(draws)
    return dict(method=background.method,strength=observed,lower_count=lower,equal_count=equal,
        numerator=numerator,denominator=denominator,percentile=100*numerator/denominator,
        degenerate_background=draws[0]==draws[-1],entry_key=background.entry_key,
        distribution_sha256=background.distribution_sha256)


def retain_entries(bank, destination, bindings, *, classification):
    """Copy each consumed entry once and check it without the original bank."""
    if classification not in {'public','controlled'}:raise ValueError('Unknown report classification.')
    destination=Path(destination).absolute()
    _path(destination,'entries')
    destination.mkdir(mode=0o700,parents=True,exist_ok=False)
    (destination/'entries').mkdir(mode=0o700);(destination/'provenance').mkdir(mode=0o700)
    copied=set();retained=[]
    for binding in bindings:
        entry=admit_entry(bank,expected_entry_key=binding['entry_key'],
            expected_distribution_sha256=binding['distribution_sha256'])
        provenance_hash=_digest(binding['provenance_sha256'])
        provenance_raw=_path(bank,'provenance/'+provenance_hash+'.json').read_bytes()
        if hashlib.sha256(provenance_raw).hexdigest()!=provenance_hash:
            raise ValueError('Calibration provenance differs from its admitted identity.')
        provenance=json.loads(provenance_raw)
        if (provenance['entry_key']!=entry.entry_key or provenance['distribution_sha256']!=entry.distribution_sha256
                or provenance['scope_sha256']!=entry.entry_key
                or provenance['classification']!=binding['classification']
                or provenance['source']['classification']!=provenance['classification']):
            raise ValueError('Calibration provenance and distribution bindings differ.')
        if classification=='public' and provenance['classification']!='public':
            raise ValueError('Controlled calibration cannot enter a public report.')
        if entry.entry_key not in copied:
            target=destination/'entries'/entry.entry_key;target.mkdir(mode=0o700)
            for name in ('distribution.json','scope.json','binding.json'):
                shutil.copyfile(_path(Path(bank)/'entries'/entry.entry_key,name),target/name)
            copied.add(entry.entry_key)
        target=destination/'provenance'/(provenance_hash+'.json')
        if not target.exists():_write(target,provenance)
        retained.append(dict(binding))
        admit_entry(destination,expected_entry_key=entry.entry_key,expected_distribution_sha256=entry.distribution_sha256)
    _write(destination/'consumed.json',dict(classification=classification,entries=retained))
    return retained


def _prepare_bank(bank, inputs, *, admitted_entries=(), counters=None, batch_size=8, on_entry=None,
                  on_ready=None, provenance_by_entry=None):
    """Explicit preparation with independently admitted identities for reuse."""
    from genosapiens.track2_frozen_calibration import prepare_distributions
    from genosapiens.track2_tissue_records import validate_tissue_record
    if inputs.get('schema')!='genosapiens.track2-calibration-inputs/v1':
        raise ValueError('Unsupported calibration preparation inputs.')
    if type(batch_size) is not int or not 1<=batch_size<=64:
        raise ValueError('Calibration batches require one through 64 entries.')
    for query in inputs['queries']:validate_tissue_record(query,'query')
    for row in inputs['conditions']:validate_tissue_record(row,'condition')
    counters={} if counters is None else counters
    counters.setdefault('query_permutations',0)
    counters.setdefault('aggregate_iterations',0)
    scope=scientific_scope(inputs['queries'],inputs['conditions'],direction=inputs['direction'],
        fda_snapshot_sha256=inputs['fda_snapshot_sha256'],eligible_ingredients=inputs['eligible_ingredients'],
        source_release=inputs['source_release'])
    pairs=sorted({(r['medication_entity_key'],r['method']) for r in scope['conditions'] if r['state']=='valid'})
    expected={}
    for binding in admitted_entries:
        if binding['entry_key'] in expected and expected[binding['entry_key']]!=binding:
            raise ValueError('Admitted calibration bindings conflict.')
        expected[binding['entry_key']]=binding
    entries={pair:entry_scope(scope,*pair) for pair in pairs}
    reused={};missing=[]
    for pair,entry in entries.items():
        key=identity(entry)
        present=_path(Path(bank)/'entries',key).exists()
        if key not in expected:
            if present:raise ValueError('An existing entry lacks an independently admitted identity.')
            missing.append(pair);continue
        if not present:
            missing.append(pair);continue
        saved=expected[key]
        reused[pair]=admit_entry(bank,expected_entry_key=key,
            expected_distribution_sha256=saved['distribution_sha256'])
    bindings={}
    def publish(pair,values):
        scope=entries[pair];key=identity(scope)
        provenance=inputs['provenance'] if provenance_by_entry is None else provenance_by_entry[key]
        record=dict(schema=BACKGROUND_SCHEMA,method=scope['method'],direction=scope['direction'],
            entry_key=key,draw_count=10000,draws=values,numerical_contract=policy()['calibration']['numerical_contract'])
        validate_contract(record,'frozen-background')
        content=dict(entry_key=key,distribution_sha256=identity(record),scope_sha256=key)
        prepared=dict(content,provenance_sha256=identity(dict(content,classification=inputs['classification'],source=provenance)),
            classification=inputs['classification'])
        if key in expected and expected[key]['distribution_sha256']!=prepared['distribution_sha256']:
            raise ValueError('Prepared calibration differs from its independently expected identity.')
        if on_ready is not None:on_ready(prepared)
        binding=publish_entry(bank,scope,values,provenance=provenance,classification=inputs['classification'])
        if binding!=prepared:raise ValueError('Calibration publication changed its prepared identity.')
        bindings[pair]=binding
        if on_entry is not None:on_entry(binding)
    for pair,background in reused.items():publish(pair,list(background.draws))
    for start in range(0,len(missing),batch_size):
        batch=missing[start:start+batch_size]
        # Each query starts the same deterministic stream in every batch.
        # A complete ingredient/method never splits across batches.
        selected=[row for row in inputs['conditions']
            if row['state']=='valid' and (row['medication_entity_key'],row['method']) in batch]
        query_ids={row['query_sha256'] for row in selected}
        draws=prepare_distributions([q for q in inputs['queries'] if q['query_sha256'] in query_ids],selected,
            direction=inputs['direction'],eligible_ingredients=inputs['eligible_ingredients'],counters=counters)
        for pair in batch:publish(pair,draws[pair])
    return dict(schema='genosapiens.track2-calibration-resource/v1',scope_sha256=identity(scope),
        classification=inputs['classification'],entries=[bindings[pair] for pair in entries])


def prepare_bank(bank, inputs, *, admitted_entries=(), counters=None, batch_size=8, on_entry=None,
                 on_ready=None, provenance_by_entry=None):
    import time
    bank=Path(bank).absolute();_path(bank,'preparation-attempts')
    bank.mkdir(mode=0o700,parents=True,exist_ok=True)
    if bank.stat().st_uid!=os.geteuid() or bank.stat().st_mode & 0o077:
        raise ValueError('Calibration preparation requires a private owned bank.')
    parent=_path(bank,'preparation-attempts');parent.mkdir(mode=0o700,exist_ok=True)
    attempt=parent/uuid.uuid4().hex;attempt.mkdir(mode=0o700)
    _write(attempt/'started.json',dict(status='incomplete',inputs_sha256=identity(inputs),
        classification=inputs.get('classification')))
    counters={} if counters is None else counters
    start=time.monotonic()
    try:
        result=_prepare_bank(bank,inputs,admitted_entries=admitted_entries,counters=counters,
            batch_size=batch_size,on_entry=on_entry,on_ready=on_ready,provenance_by_entry=provenance_by_entry)
    except BaseException as error:
        _write(attempt/'outcome.json',dict(status='incomplete',error_type=type(error).__name__,
            elapsed_seconds=time.monotonic()-start,operation_counts=counters))
        raise
    _write(attempt/'outcome.json',dict(status='complete',resource_contract_sha256=identity(result),
        elapsed_seconds=time.monotonic()-start,operation_counts=counters))
    return result


def main():
    import argparse
    import time
    parser=argparse.ArgumentParser(description='Explicit preparation of frozen calibration from retained input records.')
    parser.add_argument('inputs',type=Path);parser.add_argument('bank',type=Path)
    parser.add_argument('--expected-input-sha256',required=True)
    parser.add_argument('--resource-contract',type=Path,required=True)
    parser.add_argument('--admitted-contract',type=Path)
    parser.add_argument('--expected-contract-sha256')
    args=parser.parse_args()
    if args.resource_contract.exists():parser.error('Choose a new resource contract destination.')
    raw=_path(args.inputs.parent,args.inputs.name).read_bytes()
    if hashlib.sha256(raw).hexdigest()!=_digest(args.expected_input_sha256):
        parser.error('Preparation inputs differ from their expected identity.')
    admitted=[]
    if args.admitted_contract is not None:
        prior=_path(args.admitted_contract.parent,args.admitted_contract.name).read_bytes()
        if hashlib.sha256(prior).hexdigest()!=_digest(args.expected_contract_sha256):
            parser.error('The resource contract differs from its independently expected identity.')
        contract=json.loads(prior)
        if contract['schema']!='genosapiens.track2-calibration-resource/v1':parser.error('Unsupported calibration resource contract.')
        admitted=contract['entries']
    elif args.expected_contract_sha256 is not None:parser.error('An expected contract digest requires a contract path.')
    start=time.monotonic();counters={}
    resource=prepare_bank(args.bank,json.loads(raw),admitted_entries=admitted,counters=counters)
    _write(args.resource_contract,resource)
    print(json.dumps(dict(status='complete',resource_contract_sha256=identity(resource),
        elapsed_seconds=time.monotonic()-start,operation_counts=counters)))


if __name__=='__main__':main()
