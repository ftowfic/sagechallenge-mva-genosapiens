"""Compact source-bound selection inventory, separate from scientific evidence."""
import hashlib
import json
import os
from pathlib import Path

from genosapiens.track2_exclusions import _identity, DIGEST, MAX_LINE_BYTES, BUFFER_BYTES
from genosapiens.track2_input import Track2InputError
from genosapiens.track2_published_genotype import canonical

SCHEMA = 'genosapiens.track2-deferred-inventory/v1'
DIRECTORY = 'deferred-alleles'
ROWS = 'records.jsonl'
MANIFEST = 'manifest.json'
FIELDS = {'identity_id', 'identity', 'source_occurrences', 'annotation_summary',
          'reason', 'policy_version', 'policy_sha256', 'source_versions'}


def require(value):
    if not value:
        raise Track2InputError('Deferred inventory is incomplete or differs from its selection.')


def binding(selection, retention):
    return {**{k: selection[k] for k in ('source_sha256', 'reference_sha256', 'scope_sha256',
        'sample', 'classification', 'evidence_origin')},
        **{k: retention[k] for k in ('policy_sha256', 'policy_version', 'catalog_sha256')}}


def check_row(row, bound):
    require(type(row) is dict and set(row) == FIELDS)
    require(isinstance(row['identity_id'], str) and row['identity_id']
        and type(row['identity']) is dict and row['identity']
        and isinstance(row['source_occurrences'], list) and row['source_occurrences']
        and all(isinstance(v, str) and v for v in row['source_occurrences'])
        and len(set(row['source_occurrences'])) == len(row['source_occurrences'])
        and isinstance(row['annotation_summary'], list)
        and all(isinstance(v, str) for v in row['annotation_summary'])
        and row['reason'] == 'deferred_no_retention_evidence'
        and row['policy_sha256'] == bound['policy_sha256']
        and row['policy_version'] == bound['policy_version']
        and isinstance(row['source_versions'], list))
    from genosapiens.track2_codec import decode
    from genosapiens.variant_retention import LookupCoverage
    for value in row['source_versions']:
        decode(LookupCoverage, value)


def describe(rows, selection, retention):
    bound = binding(selection, retention)
    digest = hashlib.sha256()
    size = count = occurrences = 0
    seen, source_ids = set(), set()
    retained = {v['source_allele_id'] for v in selection['selected_alleles']}
    for row in rows:
        check_row(row, bound)
        require(row['identity_id'] not in seen
            and not source_ids.intersection(row['source_occurrences'])
            and not retained.intersection(row['source_occurrences']))
        seen.add(row['identity_id'])
        source_ids.update(row['source_occurrences'])
        raw = canonical(row)
        require(len(raw) <= MAX_LINE_BYTES)
        digest.update(raw)
        size += len(raw)
        count += 1
        occurrences += len(row['source_occurrences'])
    require(count == retention['counts']['deferred']
        and occurrences == retention['occurrences']['deferred'])
    return {'schema': SCHEMA, 'state': 'complete', 'binding': bound, 'path': ROWS,
        'sha256': digest.hexdigest(), 'bytes': size,
        'normalized_sequence_alleles': count, 'source_ALT_occurrences': occurrences}


def validate_descriptor(value, selection):
    retention = selection['coverage']['retention']
    require(type(value) is dict and set(value) == {'schema','state','binding','path','sha256','bytes',
        'normalized_sequence_alleles','source_ALT_occurrences'})
    require(value['schema'] == SCHEMA and value['state'] == 'complete' and value['path'] == ROWS
        and value['binding'] == binding(selection, retention)
        and isinstance(value['sha256'], str) and DIGEST.fullmatch(value['sha256'])
        and all(type(value[k]) is int and value[k] >= 0 for k in
            ('bytes','normalized_sequence_alleles','source_ALT_occurrences'))
        and value['normalized_sequence_alleles'] == retention['counts']['deferred']
        and value['source_ALT_occurrences'] == retention['occurrences']['deferred'])


def write_inventory(root, rows, selection):
    """Publish the completion manifest only after writing and checking all rows."""
    root = Path(root)
    expected = selection['coverage']['deferred_inventory']
    validate_descriptor(expected, selection)
    require(not root.exists() and not any(p.is_symlink() for p in (root, *root.parents)))
    root.mkdir(mode=0o700)
    with (root / ROWS).open('xb', buffering=BUFFER_BYTES) as stream:
        os.fchmod(stream.fileno(), 0o600)
        def emitted():
            for row in rows:
                raw = canonical(row)
                require(len(raw) <= MAX_LINE_BYTES)
                stream.write(raw)
                yield row
        actual = describe(emitted(), selection, selection['coverage']['retention'])
        require(actual == expected)
        stream.flush()
        os.fsync(stream.fileno())
    partial = root / (MANIFEST + '.partial')
    with partial.open('xb') as stream:
        os.fchmod(stream.fileno(), 0o600)
        stream.write(canonical(actual))
        stream.flush()
        os.fsync(stream.fileno())
    partial.rename(root / MANIFEST)


def check_inventory(root, selection):
    root = Path(root)
    expected = selection['coverage']['deferred_inventory']
    validate_descriptor(expected, selection)
    require(root.is_dir() and {p.name for p in root.iterdir()} == {ROWS, MANIFEST})
    before = {p: _identity(root / p) for p in (ROWS, MANIFEST)}
    with (root / MANIFEST).open('rb') as stream:
        require(stream.read(MAX_LINE_BYTES + 1) == canonical(expected))
    with (root / ROWS).open('rb', buffering=BUFFER_BYTES) as stream:
        def records():
            while raw := stream.readline(MAX_LINE_BYTES + 1):
                require(len(raw) <= MAX_LINE_BYTES and raw.endswith(b'\n'))
                row = json.loads(raw)
                require(raw == canonical(row))
                yield row
        require(describe(records(), selection, selection['coverage']['retention']) == expected)
    require(before == {p: _identity(root / p) for p in (ROWS, MANIFEST)})


def copy_inventory(source, target, selection):
    import shutil
    source, target = Path(source), Path(target)
    require(source.is_dir() and {p.name for p in source.iterdir()} == {ROWS, MANIFEST}
        and not any(p.is_symlink() for p in (target, *target.parents)))
    before = {p: _identity(source / p) for p in (ROWS, MANIFEST)}
    shutil.copytree(source, target)
    check_inventory(target, selection)
    require(before == {p: _identity(source / p) for p in (ROWS, MANIFEST)})
