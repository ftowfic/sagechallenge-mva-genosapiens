"""Bounded, authenticated exact-PASS exclusion artifacts.

One artifact preserves every excluded source occurrence. The completion footer
and descriptor are published only after the source stream completes.
"""
import hashlib
import json
import os
from pathlib import Path
import re

from genosapiens.track2_input import Track2InputError
from genosapiens.track2_published_genotype import canonical, sha256


SCHEMA = 'genosapiens.track2-filter-exclusions/v1'
DIRECTORY = 'filter-exclusions'
ROWS = 'records.jsonl'
MANIFEST = 'manifest.json'
MAX_LINE_BYTES = 4 * 1024 * 1024
BUFFER_BYTES = 64 * 1024
DIGEST = re.compile(r'[0-9a-f]{64}')
ROW_FIELDS = {'source_record_id', 'record_number', 'raw_line_sha256',
              'filter_value', 'alternate_count', 'reason'}
BINDING_FIELDS = ('source_sha256', 'sample', 'reference_sha256', 'scope_sha256',
                  'classification', 'evidence_origin', 'selection_policy', 'code_sha256','vcf_policy_sha256')


def require(condition, message='FILTER exclusion evidence is incomplete or invalid.'):
    if not condition:
        raise Track2InputError(message)


def exclusion_binding(bound):
    return {key: bound[key] for key in BINDING_FIELDS}


def _identity(path):
    path = Path(path)
    require(not any(p.is_symlink() for p in (path, *path.parents)))
    s = path.stat()
    return s.st_dev, s.st_ino, s.st_size, s.st_mtime_ns, s.st_ctime_ns


def _row(row, previous, total):
    require(type(row) is dict and set(row) == ROW_FIELDS)
    require(type(row['record_number']) is int and previous < row['record_number'] <= total
            and type(row['alternate_count']) is int and row['alternate_count'] > 0
            and isinstance(row['source_record_id'], str) and bool(row['source_record_id'])
            and isinstance(row['raw_line_sha256'], str) and DIGEST.fullmatch(row['raw_line_sha256'])
            and isinstance(row['filter_value'], str) and row['filter_value'] not in ('', 'PASS')
            and not any(c in row['filter_value'] for c in '\t\r\n')
            and row['reason'] == 'filter_not_exact_pass')


class ExclusionWriter:
    """Write one row at a time with a fixed byte buffer and incremental digest."""

    def __init__(self, root, *, binding, source_records):
        self.root = Path(root)
        self.binding = exclusion_binding(binding)
        self.source_records = source_records
        require(type(source_records) is int and source_records >= 0)
        require(not self.root.exists() and not self.root.is_symlink())
        require(not any(p.is_symlink() for p in self.root.parents))
        self.root.mkdir(mode=0o700)
        self.stream = (self.root / ROWS).open('xb', buffering=BUFFER_BYTES)
        (self.root / ROWS).chmod(0o600)
        self.digest = hashlib.sha256()
        self.bytes = self.count = self.alternates = self.last = 0
        self.complete = False
        try:
            self._write({'schema': SCHEMA, 'binding': self.binding, 'source_record_count': source_records})
        except BaseException:
            self.stream.close()
            raise

    def _write(self, value):
        raw = canonical(value)
        require(len(raw) <= MAX_LINE_BYTES, 'FILTER exclusion row exceeds its byte limit.')
        self.stream.write(raw)
        self.digest.update(raw)
        self.bytes += len(raw)

    def append(self, row):
        require(not self.complete)
        _row(row, self.last, self.source_records)
        self._write(row)
        self.count += 1
        self.alternates += row['alternate_count']
        self.last = row['record_number']

    def finish(self, *, source_records, source_alternates):
        require(not self.complete and source_records == self.source_records
                and type(source_alternates) is int and source_alternates >= self.alternates)
        footer = {'state': 'complete', 'excluded_record_count': self.count,
                  'excluded_alternate_count': self.alternates,
                  'source_record_count': source_records, 'source_alternate_count': source_alternates}
        self._write(footer)
        self.stream.flush()
        os.fsync(self.stream.fileno())
        self.stream.close()
        descriptor = {'schema': SCHEMA, 'binding': self.binding, 'path': ROWS,
                      'sha256': self.digest.hexdigest(), 'bytes': self.bytes, **footer}
        partial = self.root / (MANIFEST + '.partial')
        with partial.open('xb') as output:
            output.write(canonical(descriptor))
            output.flush()
            os.fsync(output.fileno())
        partial.chmod(0o600)
        partial.rename(self.root / MANIFEST)
        self.complete = True
        return descriptor

    def __enter__(self):
        return self

    def __exit__(self, kind, value, traceback):
        self.stream.close()
        if kind is None:
            require(self.complete, 'FILTER exclusion writing did not complete.')


def validate_descriptor(value, *, selection=None):
    require(type(value) is dict and set(value) == {'schema', 'binding', 'path', 'sha256', 'bytes',
        'state', 'excluded_record_count', 'excluded_alternate_count', 'source_record_count', 'source_alternate_count'})
    require(value['schema'] == SCHEMA and value['path'] == ROWS and value['state'] == 'complete'
            and isinstance(value['sha256'], str) and DIGEST.fullmatch(value['sha256'])
            and type(value['binding']) is dict and set(value['binding']) == set(BINDING_FIELDS)
            and all(type(value[k]) is int and value[k] >= 0 for k in ('bytes', 'excluded_record_count',
                'excluded_alternate_count', 'source_record_count', 'source_alternate_count'))
            and value['bytes'] > 0 and value['excluded_record_count'] <= value['source_record_count']
            and value['excluded_record_count'] <= value['excluded_alternate_count'] <= value['source_alternate_count'])
    binding = value['binding']
    from genosapiens.vcf_contracts import require_current_policy
    require_current_policy(binding['vcf_policy_sha256'])
    require(all(isinstance(binding[k], str) and DIGEST.fullmatch(binding[k]) for k in
                ('source_sha256', 'reference_sha256', 'scope_sha256', 'code_sha256'))
            and isinstance(binding['sample'], str) and bool(binding['sample'])
            and isinstance(binding['selection_policy'], str) and bool(binding['selection_policy'])
            and (binding['classification'], binding['evidence_origin']) in {
                ('controlled', 'controlled_case'), ('public', 'public_source'), ('public', 'independently_synthetic')})
    if selection is not None:
        coverage = selection['coverage']
        require(all(binding[k] == selection[k] for k in BINDING_FIELDS[:6])
                and binding['vcf_policy_sha256']==selection['vcf_policy_sha256']
                and binding['selection_policy'] == coverage['filter_policy']
                and value['source_record_count'] == coverage['source_record_count']
                and value['source_alternate_count'] == sum(coverage['source_alternate_class_counts'].values())
                and value['excluded_record_count'] == coverage['filter_excluded_record_count']
                and value['excluded_alternate_count'] == coverage['filter_excluded_alternate_count'])


def iter_exclusions(root, descriptor):
    """Authenticate complete coverage while yielding bounded individual rows.

    Callers must exhaust this iterator before accepting its evidence.
    """
    validate_descriptor(descriptor)
    root = Path(root)
    require(root.is_dir() and {p.name for p in root.iterdir()} == {ROWS, MANIFEST})
    manifest = root / MANIFEST
    manifest_identity = _identity(manifest)
    expected = canonical(descriptor)
    require(len(expected) <= MAX_LINE_BYTES)
    with manifest.open('rb') as source:
        raw = source.read(len(expected) + 1)
    require(raw == expected)
    path = root / ROWS
    before = _identity(path)
    digest = hashlib.sha256()
    size = count = alternates = previous = 0
    with path.open('rb', buffering=BUFFER_BYTES) as source:
        def read():
            nonlocal size
            raw = source.readline(MAX_LINE_BYTES + 1)
            require(raw and len(raw) <= MAX_LINE_BYTES and raw.endswith(b'\n'))
            try:
                value = json.loads(raw)
            except (ValueError, UnicodeError):
                raise Track2InputError('FILTER exclusion serialization is invalid.') from None
            require(canonical(value) == raw)
            digest.update(raw)
            size += len(raw)
            return value
        require(read() == {'schema': SCHEMA, 'binding': descriptor['binding'],
                           'source_record_count': descriptor['source_record_count']})
        while True:
            row = read()
            if isinstance(row, dict) and row.get('state') == 'complete':
                break
            _row(row, previous, descriptor['source_record_count'])
            previous = row['record_number']
            count += 1
            alternates += row['alternate_count']
            require(count <= descriptor['excluded_record_count'] and alternates <= descriptor['excluded_alternate_count'])
            yield row
        require(row == {k: descriptor[k] for k in ('state', 'excluded_record_count',
                'excluded_alternate_count', 'source_record_count', 'source_alternate_count')}
                and source.read(1) == b'' and size == descriptor['bytes']
                and count == descriptor['excluded_record_count'] and alternates == descriptor['excluded_alternate_count']
                and digest.hexdigest() == descriptor['sha256'])
    require(_identity(path) == before and _identity(manifest) == manifest_identity,
            'FILTER exclusion source changed during reading.')


def check_exclusions(root, descriptor, *, selected_record_ids=()):
    selected=set(selected_record_ids)
    for row in iter_exclusions(root, descriptor):
        require(row['source_record_id'] not in selected, 'Excluded records entered scientific selection.')


def copy_exclusions(source, target, descriptor, *, selected_record_ids=()):
    """Copy and authenticate without retaining the complete exclusion list."""
    selected=set(selected_record_ids)
    with ExclusionWriter(target, binding=descriptor['binding'],
                         source_records=descriptor['source_record_count']) as output:
        for row in iter_exclusions(source, descriptor):
            require(row['source_record_id'] not in selected, 'Excluded records entered scientific selection.')
            output.append(row)
        actual = output.finish(source_records=descriptor['source_record_count'],
                               source_alternates=descriptor['source_alternate_count'])
    require(actual == descriptor, 'Copied FILTER exclusion identity differs.')
