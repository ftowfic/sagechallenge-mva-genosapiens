"""Self-contained calibration artifacts for the challenge package contract."""
import hashlib
import json
from pathlib import Path

from genosapiens.track2_calibration_bank import (
    _digest, _path, _write, admit_entry, percentile, retain_entries,
)
from genosapiens.track2_derived_strength import aggregate_strength
from genosapiens.track2_frozen_calibration import identity

SCHEMA = 'genosapiens.track2-report-calibration/v1'


def _read(path, expected):
    path = Path(path)
    raw = _path(path.parent, path.name).read_bytes()
    if hashlib.sha256(raw).hexdigest() != _digest(expected):
        raise ValueError('A report calibration artifact differs from its expected identity.')
    return json.loads(raw)


def _resource(path, expected):
    contract = _read(path, expected)
    if contract['schema'] != 'genosapiens.track2-calibration-resource/v1':
        raise ValueError('The calibration resource contract is unsupported.')
    entries = {}
    for binding in contract['entries']:
        if binding['entry_key'] in entries:
            raise ValueError('The calibration resource repeats an entry.')
        entries[binding['entry_key']] = binding
    return contract, entries


def aggregate_key(condition_scope_sha256, row):
    return (condition_scope_sha256,row['method'],row['direction'],tuple(row['declared_tissues']),row['status'])


def _check_observation(bank, binding, row, *, aggregates=None):
    background = admit_entry(bank, expected_entry_key=binding['entry_key'],
        expected_distribution_sha256=binding['distribution_sha256'])
    scope = json.loads(_path(Path(bank)/'entries'/binding['entry_key'], 'scope.json').read_bytes())
    if (row['medication_entity_key'] != scope['medication_entity_key'] or row['method'] != scope['method']
            or row['direction'] != scope['direction']):
        raise ValueError('The observed score differs from its frozen scientific scope.')
    if row.get('condition_scope_sha256') is not None and row['condition_scope_sha256'] != identity(scope['conditions']):
        raise ValueError('The observed conditions and queries differ from the frozen scientific scope.')
    aggregates={} if aggregates is None else aggregates
    key=aggregate_key(identity(scope['conditions']),row)
    if key not in aggregates:
        conditions = (dict(condition, query_sha256=condition['query_key']) for condition in scope['conditions'])
        aggregates[key] = aggregate_strength(conditions, method=row['method'], direction=row['direction'],
            declared_tissues=row['declared_tissues'], incomplete=row['status']=='provisional')
    expected=aggregates[key]
    # Undefined references use scientific query keys in the resource scope. Source query
    # labels remain in the original derived export and do not change score arithmetic.
    for field in ('status','strength','n_conditions','n_valid','n_undefined','n_cells','n_times',
                  'contributing_tissues','missing_tissues','tissue_strengths'):
        if row[field] != expected[field]:
            raise ValueError('The observed aggregate differs from retained calibration conditions.')
    if row['status'] != 'complete':
        raise ValueError('A calibrated report score requires complete scientific execution.')
    calculated = percentile(background, row['strength'])
    if any(row.get(key) != value for key,value in calculated.items()):
        raise ValueError('A reported percentile differs from its retained distribution.')
    return calculated


def publish_calibration(destination, *, bank, resource_contract, expected_resource_sha256,
                        observations, classification, context=None):
    """Retain consumed distributions and independently check percentile counts.

    A current producer may seed the owned context with checked summaries.
    Saved rankings provide no seed authority. Unseeded keys calculate their
    aggregates normally, including standalone readers and resumed composition.
    """
    if context is None:
        from genosapiens.track2_checked_evidence import CheckedEvidenceContext
        with CheckedEvidenceContext() as owned:
            return publish_calibration(destination,bank=bank,resource_contract=resource_contract,
                expected_resource_sha256=expected_resource_sha256,observations=observations,
                classification=classification,context=owned)
    contract, entries = _resource(resource_contract, expected_resource_sha256)
    if classification not in {'public','controlled'} or (
            classification=='public' and contract['classification']!='public'):
        raise ValueError('The report cannot reduce calibration classification.')
    from genosapiens.track2_scientific_records import stored_json
    from genosapiens.track2_scientific_storage import Records
    observations = stored_json(Records(observations))
    consumed = {}; seen = set()
    aggregates={} if context is None else context._calibration_aggregates
    if context is not None:context._check_active()
    for row in observations:
        key = row['medication_entity_key'],row['objective_id'],row['method']
        if key in seen:raise ValueError('A report calibration observation is duplicated.')
        seen.add(key)
        binding = entries.get(row['entry_key'])
        if binding is None:raise ValueError('An observed score has no independently admitted calibration entry.')
        _check_observation(bank, binding, row,aggregates=aggregates)
        consumed[binding['entry_key']] = binding
    destination = Path(destination)
    bindings = [consumed[key] for key in sorted(consumed)]
    retain_entries(bank, destination, bindings, classification=classification)
    # Retain the independently admitted contract, including identities for unused entries.
    # Only consumed distribution bytes are required for a self-contained report.
    with (destination/'resource-contract.json').open('xb') as stream:
        stream.write(_path(Path(resource_contract).parent, Path(resource_contract).name).read_bytes())
    (destination/'resource-contract.json').chmod(0o400)
    _write_records(destination/'observations.json', observations)
    manifest = dict(schema=SCHEMA, classification=classification, resource_contract_sha256=expected_resource_sha256,
        observations_sha256=_file_digest(destination/'observations.json'), consumed_sha256=hashlib.sha256(
            (destination/'consumed.json').read_bytes()).hexdigest(), entries=bindings)
    _write(destination/'manifest.json', manifest)
    if context is None:
        _check_calibration(destination,expected_manifest_sha256=identity(manifest),aggregates=aggregates)
    else:context.checked_calibration(destination,expected_manifest_sha256=identity(manifest))
    return dict(path='calibration/manifest.json', sha256=identity(manifest))


def check_calibration(root, *, expected_manifest_sha256):
    """Verify frozen lookup using only artifacts inside this report component."""
    from genosapiens.track2_disk_memo import DiskMemo
    with DiskMemo() as aggregates:
        return _check_calibration(root,expected_manifest_sha256=expected_manifest_sha256,aggregates=aggregates)


def _check_calibration(root, *, expected_manifest_sha256, aggregates):
    root = Path(root)
    manifest = _read(root/'manifest.json', expected_manifest_sha256)
    if set(manifest) != {'schema','classification','resource_contract_sha256','observations_sha256',
                        'consumed_sha256','entries'} or manifest['schema'] != SCHEMA:
        raise ValueError('The report calibration contract is unsupported.')
    contract, admitted = _resource(root/'resource-contract.json', manifest['resource_contract_sha256'])
    classification = manifest['classification']
    if classification not in {'public','controlled'} or (classification=='public' and contract['classification']!='public'):
        raise ValueError('The report calibration classification differs.')
    consumed = _read(root/'consumed.json', manifest['consumed_sha256'])
    if consumed != dict(classification=classification, entries=manifest['entries']):
        raise ValueError('The consumed calibration inventory differs.')
    bindings = {}
    for binding in manifest['entries']:
        key = binding['entry_key']
        if key in bindings or admitted.get(key) != binding:
            raise ValueError('A consumed entry differs from the admitted resource.')
        provenance = _read(root/'provenance'/(binding['provenance_sha256']+'.json'), binding['provenance_sha256'])
        if (provenance['entry_key'] != key or provenance['scope_sha256'] != key
                or provenance['distribution_sha256'] != binding['distribution_sha256']
                or provenance['classification'] != binding['classification']
                or provenance['source']['classification'] != binding['classification']
                or classification=='public' and binding['classification']!='public'):
            raise ValueError('A retained provenance binding differs.')
        bindings[key] = binding
    from genosapiens.track2_scientific_reader import read_scientific
    if _file_digest(root/'observations.json')!=manifest['observations_sha256']:
        raise ValueError('A report calibration artifact differs from its expected identity.')
    observations = read_scientific(root/'observations.json')
    used = set(); seen = set()
    for row in observations:
        key = row['entry_key']
        score_key = row['medication_entity_key'], row['objective_id'], row['method']
        if key not in bindings or score_key in seen:
            raise ValueError('A retained calibration observation has an invalid binding.')
        seen.add(score_key); used.add(key)
        _check_observation(root, bindings[key], row,aggregates=aggregates)
    if used != set(bindings):raise ValueError('The calibration package retains an unconsumed distribution.')
    expected = {'manifest.json','resource-contract.json','observations.json','consumed.json'}
    for binding in bindings.values():
        expected.update('entries/'+binding['entry_key']+'/'+name for name in ('binding.json','distribution.json','scope.json'))
        expected.add('provenance/'+binding['provenance_sha256']+'.json')
    actual = set()
    for path in root.rglob('*'):
        if path.is_symlink():raise ValueError('A calibration artifact contains a symbolic link.')
        if path.is_file():actual.add(path.relative_to(root).as_posix())
        elif not path.is_dir():raise ValueError('A calibration artifact has an unsupported file type.')
    if actual != expected:raise ValueError('The calibration package inventory differs.')
    return manifest


def _file_digest(path):
    from genosapiens.track2_scientific_storage import file_identity
    return file_identity(path)['sha256']


def _write_records(path,value):
    """Stream the established compact JSON contract without its final newline."""
    import os
    from genosapiens.track2_scientific_storage import json_chunks
    with Path(path).open('xb') as stream:
        pending=None
        for chunk in json_chunks(value):
            if pending is not None:stream.write(pending)
            pending=chunk
        if pending!=b'\n':raise ValueError('Unexpected canonical record terminator.')
        stream.flush();os.fsync(stream.fileno())
    Path(path).chmod(0o400)
