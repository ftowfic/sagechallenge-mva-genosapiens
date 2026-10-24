"""Public discovery acquisition and immutable, offline catalog admission."""
import hashlib
import json
import os
from pathlib import Path
import sqlite3
import urllib.request

from genosapiens.track2_operator_inputs import require
from genosapiens.track2_published_genotype import canonical, sha256
from genosapiens.track2_run_state import atomic_record


def public_path(root, relative):
    root = Path(root).absolute()
    path = root / relative
    require(not Path(relative).is_absolute() and '..' not in Path(relative).parts
        and path.resolve() == path and path.is_relative_to(root),
        'The discovery resource path is invalid.')
    return path


def file_digest(path):
    with Path(path).open('rb') as stream:
        return hashlib.file_digest(stream, 'sha256').hexdigest()


def file_state(path):
    info = Path(path).stat()
    return dict(size=info.st_size, mtime_ns=info.st_mtime_ns, ctime_ns=info.st_ctime_ns,
        inode=info.st_ino, device=info.st_dev)


def validate_partition_inventory(contract, integrity, croissant):
    """Require the complete declared datasets from retained release metadata."""
    datasets = set(contract['required_fields'])
    expected = {p['path']: p['sha1'] for p in contract['partitions']}
    require(len(expected) == len(contract['partitions']), 'Duplicate source partitions.')
    observed = {}
    for line in integrity.splitlines():
        digest, name = line.split(maxsplit=1)
        name = name.removeprefix('./output/')
        if name.split('/')[0] in datasets and name.endswith('.parquet'):
            require(name not in observed, 'Duplicate publisher partitions.')
            observed[name] = digest
    require(expected == observed and all(any(p.startswith(d+'/') for p in expected)
        for d in datasets), 'The declared Open Targets partition inventory is incomplete or changed.')
    schemas = {r['name']: {f['name'] for f in r['field']} for r in croissant['recordSet']}
    for dataset, fields in contract['required_fields'].items():
        require(set(fields) <= schemas.get(dataset, set()),
            'The retained Open Targets schema differs: '+dataset)
    return expected


def acquire_sources(contract, reference_root):
    """Acquire only the declared drug partitions, with resumable producer receipts."""
    root = Path(reference_root).absolute()
    prefix = contract['output_prefix']
    metadata = {}
    for entry in contract['metadata']:
        path = public_path(root, prefix+'/'+entry['name'])
        require(path.is_file(), 'Restore retained Open Targets release metadata: '+str(path))
        raw = path.read_bytes()
        require(len(raw) == entry['size_bytes'] and sha256(raw) == entry['sha256'],
            'The retained Open Targets metadata changed: '+entry['name'])
        metadata[entry['name']] = raw
    inventory = validate_partition_inventory(contract, metadata['release_data_integrity'].decode(),
        json.loads(metadata['croissant.json']))
    contract_digest = sha256(canonical(contract))
    receipt = public_path(root, 'manifests/track2/discovery/acquisition-'+contract_digest+'.json')
    receipt.parent.mkdir(parents=True, exist_ok=True)
    import fcntl
    with (receipt.parent/('.acquisition-'+contract_digest+'.lock')).open('a') as lock:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        previous = json.loads(receipt.read_bytes()) if receipt.exists() else {}
        require(not previous or previous.get('source_contract_sha256') == contract_digest,
            'The acquisition receipt has a different source binding.')
        old = {r['partition']: r for r in previous.get('files', [])}
        if previous.get('status') == 'complete' and any(
                not public_path(root,row['path']).is_file()
                or file_state(public_path(root,row['path'])) != row['file_state'] for row in old.values()):
            from genosapiens.portable_resources import preserve_record
            preserve_record(receipt)
        state = dict(schema='genosapiens.track2-discovery-acquisition/v1', classification='public',
            status='incomplete', release=contract['release'], source_contract_sha256=contract_digest,
            schema_sha256=sha256(canonical(contract['required_fields'])),
            metadata=contract['metadata'], files=[])
        atomic_record(receipt, state)
        try:
            for name, expected in sorted(inventory.items()):
                path = public_path(root, prefix+'/'+name)
                row = old.get(name)
                if row and path.is_file() and row['sha1'] == expected and row['file_state'] == file_state(path):
                    state['files'].append(row)
                    continue
                path.parent.mkdir(parents=True, exist_ok=True)
                pending = public_path(root, prefix+'/'+name+'.partial')
                h1, h256, size = hashlib.sha1(), hashlib.sha256(), 0
                # A complete existing file without a receipt needs one admission read.
                if path.exists():
                    source = path.open('rb')
                    destination = None
                else:
                    source = urllib.request.urlopen(contract['base_url']+name, timeout=120)
                    destination = pending.open('wb')
                try:
                    with source:
                        while block := source.read(1024*1024):
                            h1.update(block)
                            h256.update(block)
                            size += len(block)
                            if destination is not None:
                                destination.write(block)
                    if destination is not None:
                        destination.flush()
                        os.fsync(destination.fileno())
                finally:
                    if destination is not None:
                        destination.close()
                require(h1.hexdigest() == expected, 'An Open Targets partition failed its publisher checksum: '+name)
                if destination is not None:
                    pending.chmod(0o444)
                    pending.replace(path)
                row = dict(partition=name, path=path.relative_to(root).as_posix(),
                    url=contract['base_url']+name, sha1=expected, sha256=h256.hexdigest(),
                    size=size, file_state=file_state(path))
                state['files'].append(row)
                atomic_record(receipt, state)
            state['status'] = 'complete'
            atomic_record(receipt, state)
        except BaseException:
            state['status'] = 'incomplete'
            atomic_record(receipt, state)
            raise
    return dict(path=receipt.relative_to(root).as_posix(), sha256=file_digest(receipt), **state)


def _json(path):
    with Path(path).open('rb') as stream:
        return json.load(stream)


def retained_matrix(contract, root):
    """Admit the producer's matrix identity without rereading matrix content."""
    source = contract['lincs']
    receipt = public_path(root, source['root']+'/'+source['matrix_receipt'])
    transfer = public_path(root, source['root']+'/'+source['transfer_receipt'])
    require(receipt.is_file() and transfer.is_file(),
        'Restore the retained LINCS matrix-integrity.json and matrix-transfer-state.json receipts.')
    data, progress = _json(receipt), _json(transfer)
    require(data.get('publisher_match') is True and data['decompressed_file'] == source['matrix_file']
        and data['decompressed_bytes'] == source['matrix_size'] and data['decompressed_sha256'] == source['matrix_sha256']
        and progress['status'].startswith('download_complete') and progress['total_bytes'] == data['compressed_bytes'],
        'The retained LINCS matrix receipts conflict with the source contract.')
    matrix = public_path(root, source['root']+'/'+source['matrix_file'])
    require(matrix.is_file() and matrix.stat().st_size == source['matrix_size'],
        'Restore the installed LINCS matrix from its retained producer output. Automatic matrix download is disabled.')
    return dict(path=matrix.relative_to(root).as_posix(), size=source['matrix_size'],
        producer_receipt_path=receipt.relative_to(root).as_posix(), producer_receipt_sha256=file_digest(receipt),
        content_sha256=source['matrix_sha256'])


def discovery_input_contract(contract):
    source = contract['lincs']
    target = contract['open_targets']
    return dict(matrix=dict(path=source['root']+'/'+source['matrix_file'],
        sha256=source['matrix_sha256'], size_bytes=source['matrix_size']),
        target=dict(path=target['target_manifest_path'], sha256=target['target_manifest_sha256']))


def admit_discovery_inputs(contract, reference_root):
    """Explicitly authenticate installed or copied inputs before catalog production."""
    from genosapiens.portable_resources import authenticate, completion, reusable
    root = Path(reference_root).absolute()
    bindings = discovery_input_contract(contract)
    if saved := reusable(root, 'discovery-inputs', bindings): return saved
    try:
        return require_discovery_admission(contract, root)
    except ValueError:
        pass
    retained_matrix(contract, root)
    progress = public_path(root, 'manifests/installation/discovery-inputs-progress.json')
    atomic_record(progress, dict(status='incomplete'))
    matrix = bindings['matrix']
    path = public_path(root, matrix['path'])
    # Matrix preparation already authenticates these bytes. Reuse its sealed
    # local completion when this is the standard production layout.
    receipt = _json(public_path(root, contract['lincs']['root']+'/'+contract['lincs']['matrix_receipt']))
    keys = ('compressed_file', 'compressed_bytes', 'sha512', 'decompressed_file',
            'decompressed_bytes', 'decompressed_sha256')
    prepared = None
    if contract['lincs']['root'] == 'lincs/GSE92742' and all(k in receipt for k in keys):
        prepared = reusable(root, 'lincs-matrix', {k:receipt[k] for k in keys})
    if prepared is None:
        authenticate(path, matrix['sha256'], matrix['size_bytes'])
    paths = [path]
    target = public_path(root, bindings['target']['path'])
    require(target.stat().st_size <= 8*1024**2, 'The target manifest exceeds its metadata limit.')
    authenticate(target, bindings['target']['sha256'], target.stat().st_size)
    paths.append(target)
    for row in _json(target)['outputs']:
        path = public_path(root, row['path'])
        authenticate(path, row['sha256'], row['size_bytes'])
        paths.append(path)
    record = completion(root, 'discovery-inputs', bindings, paths)
    atomic_record(progress, dict(status='complete'))
    return record


def require_discovery_admission(contract, root):
    """Use local producer evidence. Never authenticate large inputs at startup."""
    from genosapiens.portable_resources import reusable
    bindings = discovery_input_contract(contract)
    if record := reusable(root, 'discovery-inputs', bindings):
        return dict(local_completion_sha256=sha256(canonical(record)))
    # Existing catalogs retain their original trust contract. They can support
    # a software-only rebuild only while their authenticated local files match.
    for path in sorted(public_path(root, 'derived/track2/discovery/v1').glob('*/manifest.json')):
        try:
            descriptor = _json(path)
            matrix = descriptor['matrix']
            if (matrix['path'] != bindings['matrix']['path']
                    or matrix['content_sha256'] != bindings['matrix']['sha256']
                    or matrix['size'] != bindings['matrix']['size_bytes']): continue
            if not any(row['path'] == bindings['target']['path']
                    and row['sha256'] == bindings['target']['sha256']
                    for row in descriptor['source_manifests']): continue
            load_discovery_resources(descriptor, root).close()
            prior = _json(path.parent/'completion.json')['bindings']
            states = prior.get('input_file_states', {})
            if any(file_state(public_path(root, name)) != state for name,state in states.items()): continue
            target = public_path(root, bindings['target']['path'])
            current_states = {row['path']:file_state(public_path(root, row['path']))
                for row in _json(target)['outputs']}
            # Bind inputs rather than the catalog identity. Otherwise each
            # software rebuild would change its own next preparation request.
            return dict(prior_inputs_sha256=sha256(canonical(dict(contract=bindings,
                matrix_file_state=prior['matrix_file_state'], target_states=current_states))))
        except (OSError, ValueError, KeyError, TypeError):
            continue
    raise ValueError('Public source identities changed or lack local admission. Run '
        'scripts/prepare_public_resources.py discovery --admit-copied --reference-root PATH '
        'to authenticate the exact inputs before catalog construction.')


def _bound_inputs(contract, root, drug_contract):
    """Check small manifests and retained source identities during construction."""
    target = public_path(root, contract['open_targets']['target_manifest_path'])
    require(target.is_file() and file_digest(target) == contract['open_targets']['target_manifest_sha256'],
        'Restore the declared Open Targets target manifest.')
    sources = [dict(source='open_targets_targets', release=contract['open_targets']['release'],
        path=target.relative_to(root).as_posix(), sha256=file_digest(target))]
    targets = []
    for row in _json(target)['outputs']:
        if '/target/' in row['path'] and row['path'].endswith('.parquet'):
            path = public_path(root, row['path'])
            require(path.is_file() and path.stat().st_size == row['size_bytes'], 'Restore the declared Open Targets target partition: '+row['path'])
            targets.append(path)
    require(targets, 'The target manifest has no target partitions.')
    acquisition = public_path(root, 'manifests/track2/discovery/acquisition-'+sha256(canonical(drug_contract))+'.json')
    require(acquisition.is_file(), 'Prepare the declared Open Targets drug partitions before catalog construction.')
    acquired = _json(acquisition)
    require(acquired.get('status') == 'complete' and acquired['source_contract_sha256'] == sha256(canonical(drug_contract))
        and {r['partition']: r['sha1'] for r in acquired['files']} == {r['path']:r['sha1'] for r in drug_contract['partitions']},
        'Open Targets drug acquisition is incomplete.')
    partitions = {'drug_molecule': [], 'drug_mechanism_of_action': []}
    for row in acquired['files']:
        path = public_path(root, row['path'])
        require(path.is_file() and file_state(path) == row['file_state'], 'An acquired Open Targets partition changed.')
        partitions[row['partition'].split('/')[0]].append(path)
    sources.append(dict(source='open_targets_drugs', release=drug_contract['release'],
        path=acquisition.relative_to(root).as_posix(), sha256=file_digest(acquisition)))
    paths = {}
    for name, manifest_key, identity_key in [('lincs','metadata_manifest','name'), ('fda','manifest','file')]:
        settings = contract[name]
        manifest = public_path(root, settings['root']+'/'+settings[manifest_key])
        require(manifest.is_file(), 'Restore the retained public source manifest: '+str(manifest))
        entries = {row[identity_key]:row for row in _json(manifest)}
        paths[name] = {}
        for kind, filename in settings['metadata_files' if name=='lincs' else 'files'].items():
            require(filename in entries, 'A required public asset has no retained acquisition identity: '+filename)
            row = entries[filename]
            path = public_path(root, settings['root']+'/'+filename)
            require(path.is_file() and path.stat().st_size == row['size' if name=='lincs' else 'bytes'],
                'Restore the declared public source asset: '+str(path))
            digest = file_digest(path)
            require(digest == row['sha256'] and (name!='lincs' or digest == settings['metadata_sha256'][filename]),
                'A public source asset differs from its retained receipt: '+filename)
            paths[name][kind] = path
        sources.append(dict(source=name, release=settings.get('release',settings.get('snapshot')),
            path=manifest.relative_to(root).as_posix(), sha256=file_digest(manifest)))
    return dict(targets=targets, partitions=partitions, paths=paths, sources=sources)


def build_catalog(contract, reference_root, output, *, policy=None, builder=None, drug_contract=None):
    from genosapiens.track2_discovery_catalog import ingest_catalog, parquet_rows
    from genosapiens.track2_medication_catalog import build_fda_catalog
    repository = Path(__file__).resolve().parents[2]
    policy = policy or _json(repository/'config/track2/exploratory-policy-v1.json')
    builder = builder or _json(repository/'config/track2/discovery-builder-runtime.json')
    drug_contract = drug_contract or _json(repository/contract['open_targets']['contract_path'])
    root, staging = Path(reference_root).absolute(), Path(output).absolute()
    admission = require_discovery_admission(contract, root)
    require(not staging.exists() or not any(staging.iterdir()), 'Catalog construction requires empty staging output.')
    staging.mkdir(parents=True, exist_ok=True)
    atomic_record(staging/'progress.json', dict(status='incomplete'))
    inputs = _bound_inputs(contract, root, drug_contract)
    matrix = retained_matrix(contract, root)
    code_files = ['track2_discovery_resources.py','track2_discovery_catalog.py','track2_medication_catalog.py','track2_lincs_selection.py']
    code_digest = sha256(canonical({name:file_digest(Path(__file__).parent/name) for name in code_files}))
    binding = dict(source_contract_sha256=sha256(canonical(contract)), source_manifests=inputs['sources'],
        matrix=matrix, matrix_file_state=file_state(public_path(root,matrix['path'])),
        policy_sha256=sha256(canonical(policy)), identity_policy_sha256=sha256(canonical(policy['identity'])),
        catalog_schema=1, builder=dict(code_sha256=code_digest, image=builder['image']))
    binding['input_admission'] = admission
    binding['input_file_states'] = {p.relative_to(root).as_posix():file_state(p) for p in inputs['targets']}
    build_id = sha256(canonical(binding))
    fda = build_fda_catalog(inputs['paths']['fda'], policy['identity'])
    sources = [dict(source_id='open_targets' if r['source']=='open_targets_drugs' else r['source'],
        release=r['release'], manifest_sha256=r['sha256'], schema=r['source'],
        lineage_id='open_targets_chembl:'+r['release'] if r['source'].startswith('open_targets') else r['source']+':'+r['release'])
        for r in inputs['sources']]
    # The first source fixes the release used for mechanism record identity.
    coverage = ingest_catalog(staging/'catalog.sqlite',
        genes=parquet_rows(inputs['targets'], ['id','approvedSymbol','symbolSynonyms','obsoleteSymbols']),
        molecules=parquet_rows(inputs['partitions']['drug_molecule'], drug_contract['required_fields']['drug_molecule']),
        mechanisms=parquet_rows(inputs['partitions']['drug_mechanism_of_action'], drug_contract['required_fields']['drug_mechanism_of_action']),
        lincs=inputs['paths']['lincs'], fda=fda, policy=policy, sources=sources, audit_path=staging/'identity-audit.jsonl')
    atomic_record(staging/'coverage.json', dict(schema='genosapiens.track2-discovery-coverage/v1', exclusions=coverage,
        limitations=['The FDA snapshot covers resolved active single ingredients.',
            'Missing structural identity excludes some biologics.', 'Cell metadata does not establish patient tissue relevance.']))
    artifacts = {}
    for name in ['catalog.sqlite','identity-audit.jsonl','coverage.json']:
        path = staging/name
        path.chmod(0o444)
        artifacts[name] = dict(sha256=file_digest(path), size=path.stat().st_size, file_state=file_state(path))
    atomic_record(staging/'completion.json', dict(schema='genosapiens.track2-discovery-catalog-completion/v1',
        status='complete', build_id=build_id, bindings=binding, artifacts=artifacts))
    prefix = 'derived/track2/discovery/v1/'+build_id+'/'
    manifest = dict(schema='genosapiens.track2-discovery-resources/v1', classification='public', build_id=build_id,
        catalog=dict(path=prefix+'catalog.sqlite', sha256=artifacts['catalog.sqlite']['sha256'],
            size=artifacts['catalog.sqlite']['size'], schema_version=1), source_manifests=inputs['sources'], matrix=matrix,
        policy_sha256=binding['policy_sha256'], identity_policy_sha256=binding['identity_policy_sha256'], builder=binding['builder'],
        coverage_path=prefix+'coverage.json', coverage_sha256=artifacts['coverage.json']['sha256'],
        completion_sha256=file_digest(staging/'completion.json'))
    atomic_record(staging/'manifest.json', manifest)
    atomic_record(staging/'progress.json', dict(status='complete', build_id=build_id))
    return manifest


def publish_catalog(staging, destination):
    staging, destination = Path(staging), Path(destination)
    manifest, completion = _json(staging/'manifest.json'), _json(staging/'completion.json')
    require(completion['status']=='complete' and completion['build_id']==manifest['build_id']
        and file_digest(staging/'completion.json')==manifest['completion_sha256'], 'Catalog publication is incomplete.')
    for name, record in completion['artifacts'].items():
        require(file_state(public_path(staging,name))==record['file_state'], 'Catalog output changed before publication.')
    require(not destination.exists() and destination.name==manifest['build_id'], 'Catalog publication requires a new immutable build destination.')
    destination.parent.mkdir(parents=True, exist_ok=True)
    staging.rename(destination)
    for name in ['manifest.json','completion.json','progress.json']:
        (destination/name).chmod(0o444)
    return destination/'manifest.json'


def load_discovery_resources(descriptor, reference_root, *, matrix_required=True):
    """Check completed descriptors and retained file identities, then open read-only."""
    data = _json(descriptor) if isinstance(descriptor, (str,Path)) else descriptor
    if data.get('schema') == 'genosapiens.track2-discovery-resources/v2':
        from genosapiens.track2_tissue_resources import load_tissue_resources
        return load_tissue_resources(data, reference_root, matrix_required=matrix_required)
    from genosapiens.json_schema import validate_json_schema
    validate_json_schema(data, _json(Path(__file__).resolve().parents[2]/'schemas/track2/discovery-resources.schema.json'))
    expected = {'schema','classification','build_id','catalog','source_manifests','matrix','policy_sha256',
        'identity_policy_sha256','builder','coverage_path','coverage_sha256','completion_sha256'}
    require(set(data)==expected and data['schema']=='genosapiens.track2-discovery-resources/v1'
        and data['classification']=='public' and data['catalog']['schema_version']==1, 'Invalid discovery descriptor.')
    root = Path(reference_root).absolute()
    database = public_path(root, data['catalog']['path'])
    completion_path = database.parent/'completion.json'
    require(file_digest(completion_path)==data['completion_sha256'], 'The discovery completion identity changed.')
    completion = _json(completion_path)
    require(completion['status']=='complete' and completion['build_id']==data['build_id']
        and sha256(canonical(completion['bindings']))==data['build_id'], 'The discovery catalog is incomplete.')
    bindings = completion['bindings']
    require(all(data[k]==bindings[k] for k in ['source_manifests','matrix','policy_sha256','identity_policy_sha256','builder']),
        'The discovery descriptor bindings changed.')
    require(data['catalog']['sha256']==completion['artifacts']['catalog.sqlite']['sha256']
        and data['catalog']['size']==completion['artifacts']['catalog.sqlite']['size']
        and data['coverage_sha256']==completion['artifacts']['coverage.json']['sha256'], 'Discovery artifact bindings changed.')
    for name, record in completion['artifacts'].items():
        path = public_path(root, (database.parent/name).relative_to(root).as_posix())
        require(path.is_file() and file_state(path)==record['file_state'], 'A completed discovery artifact changed: '+name)
    matrix = public_path(root, data['matrix']['path'])
    receipt = public_path(root, data['matrix']['producer_receipt_path'])
    require(not matrix_required or file_digest(receipt)==data['matrix']['producer_receipt_sha256'] and matrix.is_file()
        and matrix.stat().st_size==data['matrix']['size'] and file_state(matrix)==bindings['matrix_file_state'],
        'The selected matrix differs from its retained identity.')
    connection = sqlite3.connect(database.as_uri()+'?mode=ro&immutable=1', uri=True)
    connection.row_factory = sqlite3.Row
    try:
        require(connection.execute('PRAGMA user_version').fetchone()[0]==1, 'The discovery catalog schema changed.')
        connection.execute('SELECT source_id FROM sources LIMIT 1').fetchone()
    except BaseException:
        connection.close()
        raise
    return connection


def resource_bindings(descriptor):
    data = _json(descriptor) if isinstance(descriptor, (str,Path)) else descriptor
    if data.get('schema') == 'genosapiens.track2-discovery-resources/v2':
        return resource_bindings(data['producer']) | dict(descriptor_sha256=sha256(canonical(data)),
            query_policy_sha256=sha256(canonical(data['query_policy'])),
            compatibility_sha256=sha256(canonical(data['compatibility'])),
            tissue_membership_sha256=data['tissue_companion']['sha256'])
    return dict(build_id=data['build_id'], descriptor_sha256=sha256(canonical(data)),
        catalog_sha256=data['catalog']['sha256'], matrix_sha256=data['matrix']['content_sha256'],
        policy_sha256=data['policy_sha256'], completion_sha256=data['completion_sha256'])


TISSUE_QUERY_FIELDS = frozenset(('schema', 'version', 'lincs.schema', 'lincs.aggregation',
    'lincs.undefined_condition', 'lincs.tissue_membership', 'lincs.genetic_pooling',
    'lincs.missing_feature', 'lincs.query_scope', 'lincs.drug_scope', 'lincs.score',
    'lincs.paired_research_eligibility'))


def _policy_fields(policy, prefix=''):
    result = {}
    for key, value in policy.items():
        name = prefix+key
        if isinstance(value, dict):
            result.update(_policy_fields(value, name+'.'))
        else:
            result[name] = value
    return result


def query_policy_compatibility(descriptor, reference_root, producer_policy, query_policy, *, matrix_required=True):
    """Admit one declared query transition without changing producer authority."""
    require(descriptor['policy_sha256']==sha256(canonical(producer_policy)),
        'The supplied producer policy does not match the retained catalog.')
    require(descriptor['identity_policy_sha256']==sha256(canonical(producer_policy['identity']))
        ==sha256(canonical(query_policy['identity'])), 'The medication identity policy changed.')
    root = Path(__file__).resolve().parents[2]
    legacy = _json(root/'config/track2/exploratory-policy-v1.json')
    tissue = _json(root/'config/track2/exploratory-policy-v2.json')
    require(producer_policy in (legacy, tissue) and query_policy==tissue,
        'Unsupported producer-to-query policy transition.')
    before, after = _policy_fields(producer_policy), _policy_fields(query_policy)
    changed = sorted(k for k in set(before)|set(after) if before.get(k)!=after.get(k))
    require(set(changed)<=TISSUE_QUERY_FIELDS, 'An unknown or catalog-producing policy field changed.')
    connection = load_discovery_resources(descriptor, reference_root, matrix_required=matrix_required)
    try:
        required = {
            'lincs_cells': {'cell_id','annotations_json'},
            'lincs_signatures': {'sig_id','pert_id','gene_id','cell_id','pert_type','pert_time',
                'pert_time_unit','pert_dose','pert_dose_unit','distil_nsample','distil_cc_q75','tas'},
            'lincs_features': {'pr_gene_id','symbol','measured_landmark'},
            'compound_identity': {'provider','compound_id','ingredient_key'},
            'medications': {'ingredient_key','data'},
        }
        for table, names in required.items():
            require(names <= {row[1] for row in connection.execute('PRAGMA table_info('+table+')')},
                'The retained catalog lacks required tissue query fields: '+table)
    finally:
        connection.close()
    return dict(schema='genosapiens.track2-query-compatibility/v1',
        transition='native-cell-v1-to-tissue-v2' if producer_policy==legacy else 'tissue-v2-identity',
        producer_descriptor_sha256=sha256(canonical(descriptor)),
        producer_completion_sha256=descriptor['completion_sha256'], catalog_sha256=descriptor['catalog']['sha256'],
        source_manifests=descriptor['source_manifests'], matrix=descriptor['matrix'],
        producer_policy_sha256=sha256(canonical(producer_policy)), query_policy_sha256=sha256(canonical(query_policy)),
        identity_policy_sha256=descriptor['identity_policy_sha256'], changed_query_fields=changed,
        unchanged_catalog_fields_sha256=sha256(canonical({k:v for k,v in before.items() if k not in TISSUE_QUERY_FIELDS})))


def gene_coverage(connection, gene_id, policy):
    from genosapiens.track2_lincs_selection import select_genetic_signatures, select_drug_conditions
    selections = select_genetic_signatures(connection,gene_id,policy)
    for selection in selections:
        drugs = select_drug_conditions(connection,selection['cell_id'],None,policy)
        selection['drug_signature_ids'] = [r['sig_id'] for r in drugs['conditions']]
        selection['drug_exclusions'] = drugs['exclusions']
    return dict(gene_id=gene_id,selections=selections,
        reason_codes=[] if any(s['method'] and s['drug_signature_ids'] for s in selections) else
            ['no_admitted_drug_conditions'] if any(s['method'] for s in selections) else ['no_usable_genetic_signatures'])
