"""A disease lookup derived from the completed public allele index.

Preparation streams existing exact records once. It never rebuilds source
projections or changes the published allele index. Reuse checks completion
metadata and immutable file identities, as the published index owner does.
"""
from contextlib import contextmanager
import hashlib
import json
import os
from pathlib import Path
import re
import shutil
import sqlite3
import stat
import sys
import tempfile

from genosapiens.bounded_json import load_bounded_json
from genosapiens.common_allele_public_indexes import (
    CACHE_KIB, MAX_ROW_BYTES, CommonAlleleIndexError, CommonAllelePublicIndexSession,
    retention_disease_ids, verify_common_allele_public_indexes,
)
from genosapiens.bounded_wgs_config import PublicQueryLimits
from genosapiens.track1_setup import atomic_json
from genosapiens.track1_installation import canonical, digest

SCHEMA = 'genosapiens.public-disease-lookup/v1'
PREFIX = 'derived/common-allele-public-indexes'
MANIFEST_PREFIX = 'manifests/common-allele-public-indexes'


def require(value):
    if not value:
        raise CommonAlleleIndexError('The public disease lookup is incomplete or has different bindings.')


def request_for(base_sha):
    # Version 1 uses the same explicit-identifier policy as the embedded table.
    return {'schema': SCHEMA, 'base_manifest_sha256': base_sha, 'disease_lookup_version': 1}


def identity(path):
    info = Path(path).stat()
    return (info.st_dev, info.st_ino, info.st_size, info.st_mtime_ns, info.st_ctime_ns)


def sync_directory(path):
    descriptor = os.open(path, os.O_RDONLY | os.O_DIRECTORY)
    try:
        os.fsync(descriptor)
    finally:
        os.close(descriptor)


def load_manifest(reference, pin, base_sha):
    require(type(pin) is dict and set(pin) == {'manifest_path', 'manifest_sha256'})
    request = request_for(base_sha)
    key = digest(canonical(request))
    require(pin['manifest_path'] == f'{MANIFEST_PREFIX}/disease-lookup-{key}.json')
    value, raw = load_bounded_json(Path(reference)/pin['manifest_path'],
        reference_root=reference, maximum_bytes=65536)
    require(set(value) == {'schema','status','classification','request','counts','database'}
        and value['classification'] == 'public'
        and digest(raw) == pin['manifest_sha256'] and value.get('schema') == SCHEMA
        and value.get('request') == request and value.get('status') == 'complete')
    require(set(value['counts']) == {'source_records','disease_links'}
        and all(type(number) is int and number >= 0 for number in value['counts'].values()))
    database = value['database']
    require(set(database) == {'path','sha256','size_bytes'}
        and type(database['size_bytes']) is int and database['size_bytes'] > 0
        and database['path'] == f'{PREFIX}/disease-lookup-{key}/lookup.sqlite'
        and re.fullmatch('[0-9a-f]{64}', database['sha256']) is not None)
    path = Path(reference)/database['path']
    require(path.resolve() == path.absolute() and stat.S_ISREG(path.stat().st_mode)
        and path.stat().st_mode & 0o222 == 0 and path.stat().st_size == database['size_bytes'])
    return value, path


class LookupBinding:
    def __init__(self, path, descriptor, manifest_sha256):
        self.path, self.descriptor, self.manifest_sha256 = path, descriptor, manifest_sha256
        self.identity = identity(path)

    def assert_stable(self):
        info = os.fstat(self.descriptor)
        require(not self.path.is_symlink() and identity(self.path) == self.identity
            and (info.st_dev, info.st_ino, info.st_size, info.st_mtime_ns, info.st_ctime_ns) == self.identity)


@contextmanager
def attach_disease_lookup(index, reference, pin):
    """Join disease keys to the original source rows in the bounded reader."""
    manifest, path = load_manifest(reference, pin, index.verified.manifest_sha256)
    descriptor = os.open(path, os.O_RDONLY | os.O_NOFOLLOW)
    binding = LookupBinding(path, descriptor, pin['manifest_sha256'])
    attached = False
    try:
        require(index.disease_lookup is None)
        index.connection.execute('ATTACH DATABASE ? AS disease_lookup',
            (f'file:/proc/self/fd/{descriptor}?mode=ro&immutable=1',))
        attached = True
        index.connection.execute(f'PRAGMA disease_lookup.cache_size=-{CACHE_KIB}')
        row = index.connection.execute('SELECT completion FROM disease_lookup.binding').fetchall()
        require(row == [(canonical({'request': manifest['request'], 'counts': manifest['counts']}),)])
        index.connection.execute('SELECT 1 FROM disease_lookup.retention_diseases LIMIT 1').fetchone()
        binding.assert_stable()
        index.disease_lookup = binding
        yield binding
        binding.assert_stable()
    finally:
        index.disease_lookup = None
        if attached:
            try:
                index.connection.execute('DETACH DATABASE disease_lookup')
            finally:
                os.close(descriptor)
        else:
            os.close(descriptor)


def _build(index, destination, request, progress):
    """Bound memory by one source row and fixed SQLite caches, with no replay."""
    source_records = disease_links = 0
    connection = sqlite3.connect(destination)
    try:
        connection.execute(f'PRAGMA cache_size=-{CACHE_KIB}')
        connection.execute('PRAGMA temp_store=FILE')
        connection.execute('PRAGMA journal_mode=DELETE')
        connection.executescript('''
            CREATE TABLE binding (completion BLOB NOT NULL);
            CREATE TABLE retention_diseases (
                disease_id TEXT NOT NULL, sequence INTEGER NOT NULL,
                PRIMARY KEY(disease_id,sequence)
            ) WITHOUT ROWID;
        ''')
        # Stream source sequence order instead of random coordinate lookups.
        # Oversized rows fail before Python receives their complete payload.
        rows = index.connection.execute('''SELECT e.sequence,
            CASE WHEN length(r.canonical_row)<=? THEN r.canonical_row ELSE NULL END
            FROM exact_records e NOT INDEXED CROSS JOIN records r ON r.sequence=e.sequence
            WHERE e.source IN ('clinvar','open_targets','curated_functional')
            ORDER BY e.sequence''', (MAX_ROW_BYTES,))
        for sequence, raw in rows:
            require(raw is not None)
            diseases = retention_disease_ids(json.loads(raw))
            connection.executemany('INSERT INTO retention_diseases VALUES (?,?)',
                ((disease, sequence) for disease in diseases))
            source_records += 1
            disease_links += len(diseases)
            if source_records % 10000 == 0:
                connection.commit()
            if source_records % 250000 == 0:
                index.assert_stable()
                progress(source_records, disease_links)
        index.assert_stable()
        counts = {'source_records': source_records, 'disease_links': disease_links}
        connection.execute('INSERT INTO binding VALUES (?)', (canonical({'request': request, 'counts': counts}),))
        connection.commit()
        return counts
    finally:
        connection.close()


def prepare_disease_lookup(reference, base_pin, store):
    """Caller holds the shared public preparation lock. No case inputs enter."""
    reference, store = Path(reference), Path(store)
    verified = verify_common_allele_public_indexes(reference_root=reference,
        manifest_path=base_pin['manifest_path'], expected_manifest_sha256=base_pin['manifest_sha256'],
        verify_sources=False, reuse_qualified_reconciliation=True)
    with CommonAllelePublicIndexSession(verified, limits=PublicQueryLimits(64, 10, 16*1024**2, 60)) as index:
        if index.retention_coverage().get('disease_lookup_version') == 1:
            return None
        request = request_for(verified.manifest_sha256)
        key = digest(canonical(request))
        manifest_path = f'{MANIFEST_PREFIX}/disease-lookup-{key}.json'
        completion = reference/manifest_path
        status = store/'disease-lookup-progress.json'
        if completion.exists():
            _, raw = load_bounded_json(completion, reference_root=reference, maximum_bytes=65536)
            pin = {'manifest_path': manifest_path, 'manifest_sha256': digest(raw)}
            with attach_disease_lookup(index, reference, pin):
                pass
            atomic_json(status, {'status': 'reused', 'index': pin, 'request': request})
            print('Reuse the completed public disease lookup.', file=sys.stderr, flush=True)
            return pin
        parent = reference/PREFIX
        prefix = '.disease-lookup-'+key+'-'
        # A killed producer cannot publish completion. Its own partial staging
        # directories can be replaced while this shared lock remains held.
        for abandoned in parent.glob(prefix+'*'):
            require(abandoned.is_dir() and not abandoned.is_symlink()
                and not any(p.is_symlink() for p in abandoned.rglob('*')))
            shutil.rmtree(abandoned)
        staging = Path(tempfile.mkdtemp(prefix=prefix, dir=parent))
        atomic_json(status, {'status': 'running', 'request': request})
        print('Prepare the missing public disease lookup from the published index.', file=sys.stderr, flush=True)
        def progress(records, links):
            atomic_json(status, {'status': 'running', 'request': request,
                'source_records': records, 'disease_links': links})
            print(f'Disease lookup: {records} public records processed.', file=sys.stderr, flush=True)
        try:
            database = staging/'lookup.sqlite'
            counts = _build(index, database, request, progress)
            # Seal only the new derivative. Do not rehash the source database.
            digestor = hashlib.sha256()
            with database.open('rb') as stream:
                for block in iter(lambda: stream.read(1024*1024), b''):
                    digestor.update(block)
                os.fsync(stream.fileno())
            database.chmod(0o444)
            target = parent/('disease-lookup-'+key)
            if target.exists():
                require(not target.is_symlink() and not any(p.is_symlink() for p in target.rglob('*')))
                shutil.rmtree(target)  # No completion manifest exists here.
            staging.rename(target)
            sync_directory(target)
            sync_directory(parent)
            manifest = {'schema': SCHEMA, 'classification': 'public', 'status': 'complete', 'request': request, 'counts': counts,
                'database': {'path': f'{PREFIX}/{target.name}/lookup.sqlite',
                    'sha256': digestor.hexdigest(), 'size_bytes': (target/'lookup.sqlite').stat().st_size}}
            atomic_json(completion, manifest)
            completion.chmod(0o444)
            sync_directory(completion.parent)
            pin = {'manifest_path': manifest_path, 'manifest_sha256': digest(canonical(manifest))}
            atomic_json(status, {'status': 'complete', 'request': request, 'index': pin, **counts})
            print('The public disease lookup is complete.', file=sys.stderr, flush=True)
            return pin
        except BaseException:
            atomic_json(status, {'status': 'incomplete', 'request': request})
            raise
        finally:
            if staging.exists():
                shutil.rmtree(staging)
