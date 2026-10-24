"""Authenticated public regulatory membership and bounded source context.

The index establishes explicit source links, never regulatory function or score.
Original tissue rows live in bounded compressed chunks; indexed membership points
to those chunks. Querying never opens or decompresses the original rE2G stream.
"""
from __future__ import annotations

from collections import Counter, OrderedDict
from contextlib import ExitStack
from dataclasses import asdict, dataclass
import hashlib
import json
import os
from pathlib import Path
from types import MappingProxyType
import re
import shutil
import sqlite3
import struct
import sys
import time
import uuid
import zlib

from genosapiens import admission_links as links
from genosapiens import variant_class_resources as vc
from genosapiens.evidence_model import SourcedRecord
from genosapiens.wgs_ledger import canonical
from genosapiens.wgs_shards import read_bounded_metadata
from genosapiens.wgs_read_cache import ReadFileCache, cache_bounded_readonly_connection

SCHEMA = 'genosapiens.track1-regulatory-index'
VERSION = 2
BIN_BP = 32768
# New producers use smaller authenticated retrieval units. Keep the format
# ceiling separate so completed 1-MiB indexes retain their original identities.
DEFAULT_CHUNK_BYTES = 16 * 1024
MAX_CHUNK_BYTES = 1024 * 1024
SQLITE_CACHE_KIB = 512 * 1024
COMMIT_SOURCE_BYTES = 64 * 1024 * 1024
MAX_QUERY_BP = 1024 * 1024
ELEMENT_CACHE_BYTES = 16 * 1024 * 1024
ELEMENT_CACHE_ENTRIES = 128
DISK_RESERVE_BYTES = 100_000_000_000
TARGET_KINDS = ('ensembl', 'missing', 'unmapped')
TARGET_POLICY = 'source-native-target; supported-Ensembl-or-missing-or-unmapped; no-symbol-or-alias-mapping'
_SOURCE_OFFSET = struct.Struct('<III')  # byte offset, byte length, ordinal within chunk


class RegulatoryIndexError(ValueError):
    """Static index/source/limit failures, without caller genomic values."""


def _require(ok):
    if not ok:
        raise RegulatoryIndexError('regulatory index binding, source or limit is invalid')


def _sha(content):
    return hashlib.sha256(content).hexdigest()


def _file_sha(path):
    digest = hashlib.sha256()
    with Path(path).open('rb') as source:
        cache = ReadFileCache(source.fileno())
        try:
            for block in iter(lambda: source.read(1024 * 1024), b''):
                digest.update(block)
                cache.release()
        finally:
            cache.release(force=True)
    return digest.hexdigest()


def _signature(path):
    path = Path(path)
    _require(not any(p.is_symlink() for p in (path, *path.parents)) and path.is_file())
    value = path.stat()
    return value.st_dev, value.st_ino, value.st_size, value.st_mtime_ns, value.st_ctime_ns


def regulatory_adapter_sha256():
    return _sha(canonical({p.name: _file_sha(p) for p in (
        Path(__file__), Path(links.__file__), Path(vc.__file__))}))


def _installed_sources(files, resources):
    _require(type(resources) is tuple and len(resources) == 2)
    installed, bindings = {}, []
    for pin in resources:
        data, manifest = files.source(pin)
        contract = vc.VariantClassBundleContract.from_mapping(data)
        join = contract.regulatory_join
        _require((join.registry_resource_id, join.links_resource_id,
                  join.registry_legacy_id_column, join.links_legacy_id_column) ==
                 (links.ENCODE_CCRE, links.ENCODE_RE2G, 4, 0))
        resource = contract.resource(manifest.resource_id)
        _require(resource.resource_id in (links.ENCODE_CCRE, links.ENCODE_RE2G)
                 and resource.resource_id not in installed and resource.release == 'Registry-V4')
        vc._verify_manifest_contract(manifest, resource, pin.contract_sha256)
        artifact = resource.artifact(resource.primary_artifact_role)
        _require(len(resource.artifacts) == 1)
        path = pin.root / artifact.path
        files.open(path, artifact.sha256, artifact.size_bytes)
        installed[resource.resource_id] = (pin, resource, artifact, path)
        bindings.append(dict(resource_id=resource.resource_id, contract_sha256=pin.contract_sha256,
            manifest_sha256=pin.manifest_sha256, artifact=asdict(artifact)))
    _require(resources[0].contract_sha256 == resources[1].contract_sha256)
    return installed, sorted(bindings, key=lambda row: row['resource_id'])


def _capacity(parent, required_bytes=0):
    _require(shutil.disk_usage(parent).free - required_bytes >= DISK_RESERVE_BYTES)


def _durable(path):
    descriptor = os.open(path, os.O_RDONLY)
    try:
        os.fsync(descriptor)
    finally:
        os.close(descriptor)


def build_regulatory_index(destination_parent, *, resources, reference,
        maximum_chunk_bytes=DEFAULT_CHUNK_BYTES, deadline_seconds=7200, progress=None,
        semantic_acceptance_sha256=None, expected_source_counts=None,
        preserve_ingest_checkpoint=False, resume_checkpoint=None, expected_checkpoint_sha256=None):
    """Stream pinned public sources once into an atomically published index.

    No caller count or source order assumption substitutes for actual parsing.
    New builds default to 16-KiB chunks; an explicit size remains identity-bound.
    Failures clean incomplete ingestion. Explicit recovery retains a durable
    complete-ingest checkpoint and still requires the full prepublication audit.
    """
    _require(type(maximum_chunk_bytes) is int and 256 <= maximum_chunk_bytes <= MAX_CHUNK_BYTES
             and type(deadline_seconds) is int and 1 <= deadline_seconds <= 108000)
    _require(progress is None or callable(progress))
    _require(type(preserve_ingest_checkpoint) is bool
             and ((resume_checkpoint is None and expected_checkpoint_sha256 is None)
                  or (resume_checkpoint is not None and type(expected_checkpoint_sha256) is str
                      and links._SHA.fullmatch(expected_checkpoint_sha256))))
    if semantic_acceptance_sha256 is not None:
        _require(type(semantic_acceptance_sha256) is str
                 and links._SHA.fullmatch(semantic_acceptance_sha256))
    if expected_source_counts is not None:
        _require(type(expected_source_counts) is dict and bool(expected_source_counts))
        expected_source_counts = json.loads(canonical(expected_source_counts))
    parent = Path(destination_parent).resolve(strict=True)
    _require(parent.is_dir())
    _capacity(parent)
    deadline = time.monotonic() + deadline_seconds
    scratch = lock = connection = published = owner = None
    checkpoint_ready = False
    checkpoint_content = None
    original_scratch = None
    stage = 'source_authentication'
    try:
        with ExitStack() as stack:
            files = links._Files(stack)
            installed, source_bindings = _installed_sources(files, resources)
            stage = 'reference_authentication'
            _, fasta, reference_binding = links._reference(files, reference)
            # Conservative capacity forecast, followed by live reserve checks.
            # Compressed source size is not a measured final index size.
            forecast_bytes = sum(item[2].size_bytes * 32 for item in installed.values())
            _capacity(parent, forecast_bytes)
            identity = dict(schema=SCHEMA, version=VERSION, adapter_sha256=regulatory_adapter_sha256(),
                sources=source_bindings, reference=reference_binding, bin_bp=BIN_BP,
                maximum_chunk_bytes=maximum_chunk_bytes, target_policy=TARGET_POLICY,
                semantic_acceptance_sha256=semantic_acceptance_sha256,
                expected_source_counts=expected_source_counts)
            content_id = _sha(canonical(identity))
            final = parent / content_id
            _require(not final.exists())
            lock = parent / ('.build-lock-' + content_id)
            lock.mkdir(mode=0o700)
            if resume_checkpoint is not None:
                stage = 'resume_checkpoint_authentication'
                candidate = Path(resume_checkpoint)
                _require(candidate.parent == parent and not candidate.is_symlink()
                    and candidate.name.startswith('.building-'+content_id+'-')
                    and len(candidate.name) == len('.building-')+64+1+32
                    and all(c in '0123456789abcdef' for c in candidate.name[-32:]))
                _signature(candidate/'index.sqlite')
                checkpoint_content = read_bounded_metadata(candidate/'ingest-checkpoint.json')
                _require(_sha(checkpoint_content) == expected_checkpoint_sha256)
                checkpoint = json.loads(checkpoint_content)
                _require(canonical(checkpoint) == checkpoint_content
                    and set(checkpoint) == {'schema','version','classification','status','identity',
                        'content_id','capacity_forecast_bytes','physical_counts','stream_sha256'}
                    and checkpoint['schema'] == 'genosapiens.regulatory-ingest-checkpoint'
                    and type(checkpoint['version']) is int and checkpoint['version'] == 1 and checkpoint['classification'] == 'public'
                    and checkpoint['status'] == 'ingest_complete'
                    and checkpoint['identity'] == identity and checkpoint['content_id'] == content_id)
                scratch = candidate
                checkpoint_ready = True
            else:
                stage = 'source_ingestion'
                scratch = parent / ('.building-' + content_id + '-' + uuid.uuid4().hex)
                scratch.mkdir(mode=0o700)
                db = scratch / 'index.sqlite'
                connection = sqlite3.connect(db)
                connection.executescript('''
                    PRAGMA journal_mode=DELETE;
                    PRAGMA synchronous=FULL;
                    PRAGMA temp_store=FILE;
                    PRAGMA user_version=2;
                    CREATE TABLE elements (ordinal INTEGER PRIMARY KEY, chromosome TEXT NOT NULL,
                        start INTEGER NOT NULL, end INTEGER NOT NULL, current_id TEXT NOT NULL,
                        legacy_id TEXT NOT NULL, source_text BLOB NOT NULL);
                    CREATE TABLE bins (chromosome TEXT NOT NULL, bin INTEGER NOT NULL, ordinal INTEGER NOT NULL,
                        PRIMARY KEY(chromosome,bin,ordinal)) WITHOUT ROWID;
                    CREATE TABLE memberships (legacy_id TEXT NOT NULL, gene_id TEXT NOT NULL,
                        target_kind TEXT NOT NULL CHECK(target_kind IN ('ensembl','missing','unmapped')),
                        source_rows INTEGER NOT NULL CHECK(source_rows>0),
                        PRIMARY KEY(legacy_id,gene_id)) WITHOUT ROWID;
                    CREATE TABLE chunks (id INTEGER PRIMARY KEY, first_ordinal INTEGER NOT NULL,
                        rows INTEGER NOT NULL, bytes INTEGER NOT NULL, sha256 TEXT NOT NULL, compressed BLOB NOT NULL,
                        ensembl_rows INTEGER NOT NULL, missing_rows INTEGER NOT NULL, unmapped_rows INTEGER NOT NULL);
                    CREATE TABLE chunk_memberships (legacy_id TEXT NOT NULL, gene_id TEXT NOT NULL, chunk_id INTEGER NOT NULL,
                        source_rows INTEGER NOT NULL CHECK(source_rows>0),
                        source_offsets BLOB NOT NULL,
                        PRIMARY KEY(legacy_id,gene_id,chunk_id)) WITHOUT ROWID;
                    CREATE TABLE vocabulary (resource_id TEXT NOT NULL, field TEXT NOT NULL,
                        value TEXT NOT NULL, first_ordinal INTEGER NOT NULL, source_rows INTEGER NOT NULL,
                        PRIMARY KEY(resource_id,field,value)) WITHOUT ROWID;
                    CREATE TABLE completion (singleton INTEGER PRIMARY KEY CHECK(singleton=1), content_id TEXT NOT NULL);
                ''')
                connection.execute(f'PRAGMA cache_size=-{SQLITE_CACHE_KIB}')
                descriptor = os.open(db, os.O_RDONLY | os.O_NOFOLLOW)
                stack.callback(os.close, descriptor)
                write_cache = ReadFileCache(descriptor)
                connection.set_progress_handler(lambda: int(time.monotonic() > deadline), 10000)
                ccre_rows = peak_chunk_bytes = peak_member_keys = 0
                anchors = Counter()
                physical_kinds = Counter()
                ccre_digest, re2g_digest = hashlib.sha256(), hashlib.sha256()
                pin, resource, artifact, path = installed[links.ENCODE_CCRE]
                sentinel_fields = tuple(resource.validation['sentinel']['fields'])
                sentinel_found = False
                for ordinal, text in files.lines(path):
                    _require(time.monotonic() <= deadline)
                    row = links.parse_encode_ccre_row(text)
                    ccre_digest.update(text.encode('utf-8'))
                    anchors[row['current_id'][:5]] += 1
                    sentinel_found = sentinel_found or row['columns'] == sentinel_fields
                    if row['chromosome'] in fasta.index:
                        _require(row['end'] <= fasta.index[row['chromosome']].length)
                    connection.execute('INSERT INTO elements VALUES (?,?,?,?,?,?,?)',
                        (ordinal, row['chromosome'], row['start'], row['end'], row['current_id'],
                         row['legacy_id'], text.encode('utf-8')))
                    connection.executemany('INSERT INTO bins VALUES (?,?,?)',
                        ((row['chromosome'], number, ordinal) for number in range(
                            row['start']//BIN_BP, (row['end']-1)//BIN_BP+1)))
                    ccre_rows += 1
                    connection.execute('''INSERT INTO vocabulary VALUES (?,?,?,?,1)
                        ON CONFLICT(resource_id,field,value) DO UPDATE SET source_rows=source_rows+1''',
                        (links.ENCODE_CCRE, 'class', row['class'], ordinal))
                    if ccre_rows % 10000 == 0:
                        connection.commit()
                        _capacity(parent)
                        if progress is not None:
                            progress(dict(resource_id=links.ENCODE_CCRE, ccre_rows=ccre_rows,
                                re2g_rows=0, provenance_chunks=0, peak_chunk_bytes=0, peak_member_keys=0,
                                scratch_bytes=sum(p.stat().st_size for p in scratch.iterdir() if p.is_file())))
                _require(sentinel_found)
                connection.execute('CREATE INDEX element_legacy ON elements(legacy_id)')
                connection.commit()
                pin, resource, artifact, path = installed[links.ENCODE_RE2G]
                source_rows = chunk_id = chunk_first = chunk_rows = 0
                chunk = bytearray()
                members = Counter()
                offsets = {}
                kinds = Counter()
                vocabulary = {}
                uncommitted_bytes = 0
                sentinel_fields = tuple(resource.validation['sentinel']['fields'])
                sentinel_found = False
    
                def flush():
                    nonlocal chunk_id, chunk_first, chunk_rows, peak_chunk_bytes, peak_member_keys, uncommitted_bytes
                    if not chunk:
                        return
                    _require(time.monotonic() <= deadline)
                    chunk_id += 1
                    connection.execute('INSERT INTO chunks VALUES (?,?,?,?,?,?,?,?,?)',
                        (chunk_id, chunk_first, chunk_rows, len(chunk), _sha(chunk), zlib.compress(chunk, 1),
                         *(kinds[kind] for kind in TARGET_KINDS)))
                    ordered = sorted(members)
                    connection.executemany('''INSERT INTO memberships VALUES (?,?,?,?)
                        ON CONFLICT(legacy_id,gene_id) DO UPDATE SET source_rows=source_rows+excluded.source_rows''',
                        ((element, raw_target, links.classify_encode_target(raw_target), members[element,raw_target])
                         for element, raw_target in ordered))
                    connection.executemany('INSERT INTO chunk_memberships VALUES (?,?,?,?,?)',
                        ((element, raw_target, chunk_id, members[element,raw_target], bytes(offsets[element,raw_target]))
                         for element, raw_target in ordered))
                    connection.executemany('''INSERT INTO vocabulary VALUES (?,?,?,?,?)
                        ON CONFLICT(resource_id,field,value) DO UPDATE SET source_rows=source_rows+excluded.source_rows''',
                        ((links.ENCODE_RE2G, field, value, first, count)
                         for (field, value), (first, count) in sorted(vocabulary.items())))
                    peak_chunk_bytes = max(peak_chunk_bytes, len(chunk))
                    peak_member_keys = max(peak_member_keys, len(members))
                    uncommitted_bytes += len(chunk)
                    if uncommitted_bytes >= COMMIT_SOURCE_BYTES:
                        connection.commit()
                        _capacity(parent)
                        uncommitted_bytes = 0
                    if progress is not None:
                        progress(dict(resource_id=links.ENCODE_RE2G, ccre_rows=ccre_rows,
                            re2g_rows=source_rows, provenance_chunks=chunk_id,
                            peak_chunk_bytes=peak_chunk_bytes, peak_member_keys=peak_member_keys,
                            scratch_bytes=sum(p.stat().st_size for p in scratch.iterdir() if p.is_file())))
                    chunk.clear()
                    members.clear()
                    offsets.clear()
                    kinds.clear()
                    vocabulary.clear()
                    chunk_first = chunk_rows = 0
    
                for ordinal, text in files.lines(path):
                    row = links.parse_encode_re2g_row(text)
                    sentinel_found = sentinel_found or row['columns'] == sentinel_fields
                    content = text.encode('utf-8')
                    re2g_digest.update(content)
                    physical_kinds[row['target_kind']] += 1
                    _require(len(content) <= maximum_chunk_bytes)
                    if len(chunk) + len(content) > maximum_chunk_bytes:
                        flush()
                    chunk_first = ordinal if chunk_first == 0 else chunk_first
                    member = row['legacy_id'], row['source_target_id']
                    offsets.setdefault(member, bytearray()).extend(
                        _SOURCE_OFFSET.pack(len(chunk), len(content), ordinal-chunk_first))
                    chunk.extend(content)
                    chunk_rows += 1
                    members[row['legacy_id'], row['source_target_id']] += 1
                    kinds[row['target_kind']] += 1
                    for field in ('gene_biotype', 'method', 'biosample', 'accession'):
                        key = field, row[field]
                        first, count = vocabulary.get(key, (ordinal, 0))
                        vocabulary[key] = first, count + 1
                    source_rows += 1
                flush()
                _require(sentinel_found)
                stage = 'ingest_commit'
                connection.commit()
                connection.close()
                write_cache.release(force=True)
                connection = None
                checkpoint = dict(schema='genosapiens.regulatory-ingest-checkpoint', version=1,
                    classification='public', status='ingest_complete', identity=identity,
                    content_id=content_id, capacity_forecast_bytes=forecast_bytes,
                    physical_counts=dict(ccre_rows=ccre_rows, re2g_rows=source_rows,
                        anchor_family_counts=dict(anchors),
                        rows_by_target_kind={kind: physical_kinds[kind] for kind in TARGET_KINDS}),
                    stream_sha256=dict(ccre=ccre_digest.hexdigest(), re2g=re2g_digest.hexdigest()))
        # Original immutable source owners must close before marking ingest complete.
        db = scratch/'index.sqlite'
        original_scratch = scratch
        stage = 'ingest_checkpoint_seal'
        _durable(db)
        if preserve_ingest_checkpoint or resume_checkpoint is not None:
            if checkpoint_content is None:
                checkpoint_content = canonical(checkpoint)
                with (scratch/'ingest-checkpoint.json').open('xb') as target:
                    target.write(checkpoint_content)
                    target.flush()
                    os.fsync(target.fileno())
            _durable(scratch)
            checkpoint_ready = True
            if progress is not None:
                progress(dict(event='ingest_checkpoint',path=str(scratch),
                    checkpoint_sha256=_sha(checkpoint_content),resumed=resume_checkpoint is not None,
                    physical_counts=checkpoint['physical_counts']))
        physical = checkpoint['physical_counts']
        ccre_rows, source_rows = physical['ccre_rows'], physical['re2g_rows']
        connection = sqlite3.connect(db)
        connection.executescript('PRAGMA journal_mode=DELETE; PRAGMA synchronous=FULL; PRAGMA temp_store=FILE;')
        connection.execute(f'PRAGMA cache_size=-{SQLITE_CACHE_KIB}')
        connection.set_progress_handler(lambda: int(time.monotonic() > deadline), 10000)
        descriptor = os.open(db, os.O_RDONLY | os.O_NOFOLLOW)
        write_cache = ReadFileCache(descriptor)
        try:
            stage = 'index_creation_chunk_id_lookup'
            if not connection.execute("SELECT 1 FROM sqlite_master WHERE type='index' AND name='chunk_id_lookup'").fetchone():
                connection.execute('CREATE INDEX chunk_id_lookup ON chunk_memberships(chunk_id,source_rows)')
            stage = 'index_creation_vocabulary_ordinal'
            if not connection.execute("SELECT 1 FROM sqlite_master WHERE type='index' AND name='vocabulary_ordinal'").fetchone():
                connection.execute('CREATE INDEX vocabulary_ordinal ON vocabulary(resource_id,first_ordinal)')
            stage = 'source_count_reconciliation'
            _capacity(parent)
            _require(connection.execute('SELECT coalesce(sum(rows),0) FROM chunks').fetchone()[0] == source_rows)
            counts = _source_counts(connection)
            _require(counts['ccre_rows'] == ccre_rows and counts['re2g_rows'] == source_rows)
            inventory = dict(collection='single-validation-ingest-pass',
                physical_counts=physical,
                stream_sha256=checkpoint['stream_sha256'],
                vocabulary_path='index.sqlite#vocabulary',
                vocabulary_counts=_vocabulary_counts(connection),
                target_resolution_complete=counts['unresolved_memberships'] == 0)
            connection.execute('INSERT OR REPLACE INTO completion VALUES (1,?)', (content_id,))
            connection.commit()
            _require(connection.execute('PRAGMA quick_check').fetchall() == [('ok',)])
            connection.close()
            connection = None
            _capacity(parent)
        finally:
            try:
                if connection is not None:
                    connection.close()
                    connection = None
            finally:
                try:
                    write_cache.release(force=True)
                finally:
                    os.close(descriptor)
        # Stable original source descriptors must close successfully before seal.
        _require(time.monotonic() <= deadline)
        manifest = dict(**identity, classification='public', status='complete', content_id=content_id,
            counts=counts, inventory=inventory, capacity_forecast_bytes=forecast_bytes,
            index=dict(path='index.sqlite', sha256=_file_sha(db), size_bytes=db.stat().st_size))
        stage = 'prepublication_validation'
        verified = sqlite3.connect(db.as_uri()+'?mode=ro&immutable=1', uri=True)
        verified.set_progress_handler(lambda: int(time.monotonic() > deadline), 10000)
        try:
            _validate_database_closure(verified, manifest)
        finally:
            verified.close()
        stage = 'durable_publication'
        content = canonical(manifest)
        with (scratch / 'manifest.json').open('wb' if resume_checkpoint is not None else 'xb') as target:
            target.write(content)
            target.flush()
            os.fsync(target.fileno())
        _durable(db)
        _durable(scratch)
        _require(time.monotonic() <= deadline)
        _require(not final.exists())
        os.rename(scratch, final)
        scratch = None
        published = final
        _durable(parent)
        stage = 'returned_reader_authentication'
        owner = open_regulatory_index(final, expected_manifest_sha256=_sha(content))
        lock.rmdir()
        lock = None
        if (final/'ingest-checkpoint.json').exists():
            (final/'ingest-checkpoint.json').unlink()
            _durable(final)
        published = None
        return owner
    except (Exception, KeyboardInterrupt) as error:
        if isinstance(error, RegulatoryIndexError):
            failure = error
        else:
            detail = getattr(error, 'sqlite_errorname', None)
            failure = RegulatoryIndexError('regulatory index build failed during '+stage
                + (' (SQLite '+detail+')' if detail else ' ('+type(error).__name__+')'))
        failure.stage = stage
        failure.sqlite_errorcode = getattr(error, 'sqlite_errorcode', None)
        failure.sqlite_errorname = getattr(error, 'sqlite_errorname', None)
        failure.os_errno = getattr(error, 'errno', None)
        raise failure from None
    finally:
        if published is not None:
            if owner is not None:
                owner.close()
            if checkpoint_ready and preserve_ingest_checkpoint:
                os.rename(published, original_scratch)
                scratch = original_scratch
            else:
                shutil.rmtree(published)
        if connection is not None:
            connection.close()
        if scratch is not None and not (checkpoint_ready and preserve_ingest_checkpoint):
            shutil.rmtree(scratch)
        if lock is not None and lock.exists():
            lock.rmdir()


@dataclass(frozen=True, slots=True)
class RegulatorySourceRow:
    ordinal: int
    text: str
    sha256: str


@dataclass(frozen=True, slots=True)
class RegulatoryElementLinks:
    element_ordinal: int
    element_text: str
    targets: tuple[str, ...]
    source_rows: tuple[RegulatorySourceRow, ...]
    unresolved_targets: tuple[str, ...] = ()
    # Raw target, kind, original source count and first original ordinal.
    excluded_memberships: tuple[tuple[str,str,int,int], ...] = ()


class RegulatoryIndex:
    def __init__(self, root, manifest, manifest_sha256, connection):
        self.root, self.manifest, self.manifest_sha256 = root, manifest, manifest_sha256
        self.connection = connection
        self._manifest_content = canonical(manifest)
        self._signatures = {p: _signature(p) for p in (root/'manifest.json', root/'index.sqlite')}
        self._reference_stack = None
        self._reference_pin = None
        self._reference_context = None
        self._closed = False
        self._elements = OrderedDict()
        self._element_bytes = 0

    def _remember_element(self, key, value):
        # Frozen records contain only strings, integers and tuples. Count their
        # Python storage conservatively (including repeated objects), plus map
        # overhead, rather than bounding only compressed/text bytes.
        size = 512 + sys.getsizeof(key) + sum(sys.getsizeof(v) for v in key)
        size += sys.getsizeof(value) + sys.getsizeof(value.element_text)
        for values in (value.targets, value.unresolved_targets):
            size += sys.getsizeof(values) + sum(sys.getsizeof(v) for v in values)
        size += sys.getsizeof(value.excluded_memberships)
        for membership in value.excluded_memberships:
            size += sys.getsizeof(membership)+sum(sys.getsizeof(v) for v in membership)
        size += sys.getsizeof(value.source_rows)
        for row in value.source_rows:
            size += sum(sys.getsizeof(v) for v in (row, row.ordinal, row.text, row.sha256))
        if size > ELEMENT_CACHE_BYTES:
            return
        while self._elements and (self._element_bytes + size > ELEMENT_CACHE_BYTES
                or len(self._elements) >= ELEMENT_CACHE_ENTRIES):
            _, (_, previous_size) = self._elements.popitem(last=False)
            self._element_bytes -= previous_size
        self._elements[key] = value, size
        self._element_bytes += size

    def _reference(self, pin):
        """One authenticated, stable reference lifetime per immutable index owner."""
        self.require_stable()
        if self._reference_context is None:
            stack = ExitStack()
            try:
                files = links._Files(stack)
                maps, fasta, proof = links._reference(files, pin)
                _require(proof == self.manifest['reference'])
                # Keep the requested-base source stable too; no whole-FASTA
                # rehash or sequence-equivalence assertion is introduced.
                files.open(fasta.fasta_path)
                signatures = {path: _signature(path) for path in files.opened}
                for path, stream in files.opened.items():
                    value = os.fstat(stream.fileno())
                    _require(signatures[path] == (value.st_dev, value.st_ino,
                        value.st_size, value.st_mtime_ns, value.st_ctime_ns))
                maps = MappingProxyType({key: (digest, MappingProxyType(aliases))
                    for key, (digest, aliases) in maps.items()})
                fasta.index = MappingProxyType(fasta.index)
                self._signatures.update(signatures)
                self._reference_pin = pin
                self._reference_context = maps, fasta, proof
                self._reference_stack = stack
            except BaseException:
                stack.close()
                raise
        _require(pin == self._reference_pin)
        self.require_stable()
        return self._reference_context

    def require_stable(self):
        _require(not self._closed and canonical(self.manifest) == self._manifest_content)
        _require(all(_signature(p) == signature for p, signature in self._signatures.items()))
        self.connection._read_cache.release()

    def close(self):
        if not self._closed:
            self._closed = True
            self._elements.clear()
            self._element_bytes = 0
            try:
                if self._reference_stack is not None:
                    self._reference_stack.close()
            finally:
                self.connection.close()

    def __enter__(self):
        self.require_stable()
        return self

    def __exit__(self, *_):
        try:
            self.require_stable()
        finally:
            self.close()


def _source_counts(connection):
    """Count exact source identity levels from persisted relations and chunks."""
    scalar = lambda sql: connection.execute(sql).fetchone()[0]
    return dict(ccre_rows=scalar('SELECT count(*) FROM elements'),
        anchor_family_counts=dict(connection.execute(
            'SELECT substr(current_id,1,5),count(*) FROM elements GROUP BY substr(current_id,1,5)')),
        re2g_rows=scalar('SELECT coalesce(sum(rows),0) FROM chunks'),
        provenance_chunks=scalar('SELECT count(*) FROM chunks'),
        unique_legacy_elements=scalar("SELECT count(DISTINCT legacy_id) FROM elements WHERE legacy_id!='.'"),
        unique_source_memberships=scalar('SELECT count(*) FROM memberships'),
        distinct_raw_targets=scalar('SELECT count(DISTINCT gene_id) FROM memberships'),
        stable_memberships=scalar("SELECT count(*) FROM memberships WHERE target_kind='ensembl'"),
        unresolved_memberships=scalar("SELECT count(*) FROM memberships WHERE target_kind!='ensembl'"),
        linked_memberships=scalar("SELECT count(*) FROM memberships m WHERE target_kind='ensembl' AND EXISTS (SELECT 1 FROM elements e WHERE e.legacy_id=m.legacy_id)"),
        unjoined_memberships=scalar('SELECT count(*) FROM memberships m WHERE NOT EXISTS (SELECT 1 FROM elements e WHERE e.legacy_id=m.legacy_id)'),
        unjoined_rows=scalar('SELECT coalesce(sum(source_rows),0) FROM memberships m WHERE NOT EXISTS (SELECT 1 FROM elements e WHERE e.legacy_id=m.legacy_id)'),
        provenance_memberships=scalar('SELECT count(*) FROM chunk_memberships'),
        rows_by_target_kind={kind:scalar('SELECT coalesce(sum('+kind+'_rows),0) FROM chunks') for kind in TARGET_KINDS},
        relations_by_kind={kind:scalar("SELECT count(*) FROM memberships WHERE target_kind='"+kind+"'") for kind in TARGET_KINDS},
        distinct_targets_by_kind={kind:scalar("SELECT count(DISTINCT gene_id) FROM memberships WHERE target_kind='"+kind+"'") for kind in TARGET_KINDS})


def _vocabulary_counts(connection):
    return {resource + ':' + field: dict(distinct_values=count, source_rows=rows)
            for resource, field, count, rows in connection.execute(
                'SELECT resource_id,field,count(*),sum(source_rows) FROM vocabulary GROUP BY resource_id,field')}


def _decode_chunk(compressed, size, digest, maximum):
    _require(type(size) is int and 0 < size <= maximum)
    decoder = zlib.decompressobj()
    content = decoder.decompress(compressed, size+1)
    _require(len(content) == size and decoder.eof and not decoder.unused_data
             and not decoder.unconsumed_tail and _sha(content) == digest)
    return content


def _validate_provenance(connection, manifest):
    """Verify generated artifacts, with at most one bounded chunk in memory.

    Every original physical row must have exactly one correct locator. Original
    streams are never reopened. Parser validation happened at intake; this pass
    verifies preserved bytes and their indexed source-native relations.
    """
    ccre_digest, re2g_digest = hashlib.sha256(), hashlib.sha256()
    for expected, (ordinal, chrom, start, end, anchor, element, content) in enumerate(
            connection.execute('SELECT * FROM elements ORDER BY ordinal'), 1):
        _require(ordinal == expected)
        row = links.parse_encode_ccre_row(content.decode('utf-8'))
        _require((chrom, start, end, anchor, element) ==
                 (row['chromosome'], row['start'], row['end'], row['current_id'], row['legacy_id']))
        ccre_digest.update(content)
    _require(ccre_digest.hexdigest() == next(source['artifact']['sha256']
        for source in manifest['sources'] if source['resource_id'] == links.ENCODE_CCRE))
    for chunk_id, first, rows, size, digest, compressed, *class_counts in connection.execute(
            'SELECT * FROM chunks ORDER BY id'):
        content = _decode_chunk(compressed, size, digest, manifest['maximum_chunk_bytes'])
        physical_rows = content.splitlines(keepends=True)
        _require(len(physical_rows) == rows and b''.join(physical_rows) == content)
        positions = []
        for element, target, count, locators, kind in connection.execute('''
                SELECT c.legacy_id,c.gene_id,c.source_rows,c.source_offsets,m.target_kind
                FROM chunk_memberships c JOIN memberships m
                ON c.legacy_id=m.legacy_id AND c.gene_id=m.gene_id WHERE c.chunk_id=?''', (chunk_id,)):
            _require(len(locators) == count * _SOURCE_OFFSET.size and count <= rows)
            positions.extend((*locator, element, target, kind)
                             for locator in _SOURCE_OFFSET.iter_unpack(locators))
            _require(len(positions) <= rows)
        _require(len(positions) == rows)
        actual_classes = Counter()
        offset = 0
        for relative, (physical, locator) in enumerate(zip(physical_rows, sorted(positions), strict=True)):
            begin, length, relative_ordinal, element, target, kind = locator
            _require((begin, length, relative_ordinal) == (offset, len(physical), relative))
            columns = physical.decode('utf-8').rstrip('\r\n').split('\t')
            _require(len(columns) == 9 and (columns[0], columns[1]) == (element, target))
            actual_classes[kind] += 1
            offset += len(physical)
        _require([actual_classes[kind] for kind in TARGET_KINDS] == class_counts)
        for field, value, ordinal in connection.execute('''SELECT field,value,first_ordinal
                FROM vocabulary WHERE resource_id=? AND first_ordinal>=? AND first_ordinal<?''',
                (links.ENCODE_RE2G, first, first+rows)):
            columns = physical_rows[ordinal-first].decode('utf-8').rstrip('\r\n').split('\t')
            _require(field in links._RE2G_COLUMNS and columns[links._RE2G_COLUMNS.index(field)] == value)
        re2g_digest.update(content)
    _require(manifest['inventory']['stream_sha256'] ==
             dict(ccre=ccre_digest.hexdigest(), re2g=re2g_digest.hexdigest()))


def _validate_database_closure(connection, manifest):
    connection.execute(f'PRAGMA cache_size=-{SQLITE_CACHE_KIB}')
    _require(connection.execute('PRAGMA user_version').fetchone() == (VERSION,)
             and connection.execute('SELECT content_id FROM completion WHERE singleton=1').fetchone()
                == (manifest['content_id'],))
    # Reclassify persisted identities independently of stored counters. Stream
    # distinct raw labels; never materialize every physical tissue row.
    for raw_target, kind in connection.execute(
            'SELECT DISTINCT gene_id,target_kind FROM memberships ORDER BY gene_id,target_kind'):
        _require(links.classify_encode_target(raw_target) == kind)
    counts = _source_counts(connection)
    _require(canonical(counts) == canonical(manifest['counts']))
    inventory = manifest['inventory']
    _require(inventory['collection'] == 'single-validation-ingest-pass'
             and inventory['vocabulary_path'] == 'index.sqlite#vocabulary'
             and inventory['target_resolution_complete'] == (counts['unresolved_memberships'] == 0)
             and inventory['physical_counts'] == {key: counts[key] for key in
                 ('ccre_rows', 're2g_rows', 'anchor_family_counts', 'rows_by_target_kind')}
             and inventory['vocabulary_counts'] == _vocabulary_counts(connection))
    expected = manifest['expected_source_counts']
    if expected is not None:
        _require(set(expected) <= set(counts) and all(
            canonical(counts[key]) == canonical(value) for key, value in expected.items()))
    _require(sum(counts['rows_by_target_kind'].values()) == counts['re2g_rows']
             == connection.execute('SELECT coalesce(sum(source_rows),0) FROM memberships').fetchone()[0]
             == connection.execute('SELECT coalesce(sum(source_rows),0) FROM chunk_memberships').fetchone()[0])
    for kind in TARGET_KINDS:
        _require(counts['rows_by_target_kind'][kind] == connection.execute(
            'SELECT coalesce(sum(source_rows),0) FROM memberships WHERE target_kind=?', (kind,)).fetchone()[0])
    _require(not connection.execute('''SELECT 1 FROM chunks c WHERE
        c.rows != c.ensembl_rows+c.missing_rows+c.unmapped_rows OR c.rows<=0 OR c.bytes<=0
        OR c.rows != (SELECT coalesce(sum(source_rows),0) FROM chunk_memberships m WHERE m.chunk_id=c.id)
        LIMIT 1''').fetchone())
    _require(not connection.execute('''SELECT 1 FROM memberships m WHERE m.source_rows !=
        (SELECT coalesce(sum(source_rows),0) FROM chunk_memberships c
         WHERE c.legacy_id=m.legacy_id AND c.gene_id=m.gene_id) LIMIT 1''').fetchone())
    _require(not connection.execute('SELECT 1 FROM chunk_memberships WHERE length(source_offsets)!=source_rows*12 LIMIT 1').fetchone())
    _require(not connection.execute('''SELECT chunk_id FROM chunk_memberships GROUP BY chunk_id
        EXCEPT SELECT id FROM chunks LIMIT 1''').fetchone())
    _require(not connection.execute('''SELECT legacy_id,gene_id FROM chunk_memberships GROUP BY legacy_id,gene_id
        EXCEPT SELECT legacy_id,gene_id FROM memberships LIMIT 1''').fetchone())
    _require(not connection.execute('''SELECT 1 FROM (
        SELECT id,first_ordinal,row_number() OVER (ORDER BY id) AS expected_id,
        coalesce(sum(rows) OVER (ORDER BY id ROWS BETWEEN UNBOUNDED PRECEDING AND 1 PRECEDING),0)+1 AS expected_first
        FROM chunks) WHERE id!=expected_id OR first_ordinal!=expected_first LIMIT 1''').fetchone())
    _require(connection.execute('SELECT count(*) FROM bins').fetchone()[0] == connection.execute(
        'SELECT coalesce(sum((end-1)/? - start/? + 1),0) FROM elements', (BIN_BP, BIN_BP)).fetchone()[0])
    _require(not connection.execute('''SELECT 1 FROM bins b LEFT JOIN elements e ON e.ordinal=b.ordinal
        WHERE e.ordinal IS NULL OR b.chromosome!=e.chromosome OR b.bin<e.start/? OR b.bin>(e.end-1)/?
        LIMIT 1''', (BIN_BP, BIN_BP)).fetchone())
    for resource, field, total in connection.execute(
            'SELECT resource_id,field,sum(source_rows) FROM vocabulary GROUP BY resource_id,field'):
        _require(resource in (links.ENCODE_CCRE, links.ENCODE_RE2G)
                 and total == counts['ccre_rows' if resource == links.ENCODE_CCRE else 're2g_rows'])
    _require(set(inventory['vocabulary_counts']) == {
        links.ENCODE_CCRE+':class', *(links.ENCODE_RE2G+':'+field for field in
                                    ('gene_biotype', 'method', 'biosample', 'accession'))})
    _validate_provenance(connection, manifest)


def _validate_database_header(connection, manifest):
    """Check the sealed reader interface; exhaustive closure belongs to build."""
    connection.execute(f'PRAGMA cache_size=-{SQLITE_CACHE_KIB}')
    _require(connection.execute('PRAGMA user_version').fetchone() == (VERSION,)
             and connection.execute('SELECT content_id FROM completion WHERE singleton=1').fetchone()
                == (manifest['content_id'],))
    objects = set(connection.execute("SELECT type,name FROM sqlite_master WHERE type IN ('table','index')"))
    _require({('table', name) for name in ('elements','bins','memberships','chunks',
        'chunk_memberships','vocabulary','completion')} <= objects
        and {('index', name) for name in ('element_legacy','chunk_id_lookup','vocabulary_ordinal')} <= objects)


def open_regulatory_index(root, *, expected_manifest_sha256, expected_sources=None):
    """Reuse producer completion; validate bounded metadata and queried evidence."""
    connection = None
    try:
        root = Path(root).resolve(strict=True)
        manifest_signature = _signature(root/'manifest.json')
        content = read_bounded_metadata(root/'manifest.json')
        _require(_sha(content) == expected_manifest_sha256)
        manifest = json.loads(content)
        _require(canonical(manifest) == content and manifest['schema'] == SCHEMA
                 and type(manifest['version']) is int and manifest['version'] == VERSION
                 and manifest['classification'] == 'public' and manifest['status'] == 'complete'
                 # Compatible completed producers retain the same parser,
                 # schema and source/provenance semantics. Reader/default-size
                 # changes preserve their actual immutable producer identities.
                 and manifest['adapter_sha256'] in (regulatory_adapter_sha256(),
                     'cd6389ee3678efc6f72fe242a36d47e465db36f835e6b699aa040bfd49a70ce9',
                     '6513eeb7088d265852ff32e2308beb502cc89d4cada598342f0bba8dea440d1e',
                     'f012c1f7e32bbcfecd9b5c51e59d03d4f701e957fe4a9557bebafd0d2b715433',
                     'a238bce92cfad1bad04d5b7a0054ab1ec6d07605ab808b4ce4ddbc75dc39e3b9')
                 and manifest['bin_bp'] == BIN_BP and manifest['target_policy'] == TARGET_POLICY
                 and 256 <= manifest['maximum_chunk_bytes'] <= MAX_CHUNK_BYTES)
        identity = {name: manifest[name] for name in ('schema','version','adapter_sha256',
            'sources','reference','bin_bp','maximum_chunk_bytes','target_policy',
            'semantic_acceptance_sha256','expected_source_counts')}
        if identity['semantic_acceptance_sha256'] is not None:
            _require(type(identity['semantic_acceptance_sha256']) is str
                     and links._SHA.fullmatch(identity['semantic_acceptance_sha256']))
        if identity['expected_source_counts'] is not None:
            _require(type(identity['expected_source_counts']) is dict and bool(identity['expected_source_counts']))
        _require(manifest['content_id'] == _sha(canonical(identity)) and root.name == manifest['content_id'])
        if expected_sources is not None:
            _require(manifest['sources'] == expected_sources)
        output = manifest['index']
        _require(output['path'] == 'index.sqlite')
        path = root / output['path']
        signature = _signature(path)
        _require(signature[2] == output['size_bytes'] and type(output['sha256']) is str
                 and links._SHA.fullmatch(output['sha256']))
        connection = cache_bounded_readonly_connection(path)
        connection.execute('PRAGMA query_only=ON')
        # This newly opened reader has a release-only callback. Keep file pages
        # during reader header checks; ordinary queries retain release
        # and their own deadline callback. Restore/release even on interruption.
        validation_error = None
        try:
            connection.set_progress_handler(None, 0)
            _validate_database_header(connection, manifest)
        except BaseException as error:
            validation_error = error
            raise
        finally:
            try:
                try:
                    connection.set_progress_handler(
                        lambda: (connection._read_cache.release(), 0)[1], 100000)
                except BaseException as error:
                    if validation_error is None or isinstance(validation_error, Exception):
                        validation_error = error
                    raise
                finally:
                    connection._read_cache.release(force=True)
            except BaseException:
                if validation_error is not None and not isinstance(validation_error, Exception):
                    raise validation_error
                raise
        _require(_signature(path) == signature
                 and _signature(root/'manifest.json') == manifest_signature)
        result = RegulatoryIndex(root, manifest, expected_manifest_sha256, connection)
        _require(result._signatures[path] == signature
                 and result._signatures[root/'manifest.json'] == manifest_signature)
        return result
    except BaseException as error:
        if connection is not None:
            try:
                connection.close()
            except BaseException:
                # Preserve the original failure, especially an interruption;
                # the shared close owner still closes SQLite and its file fd.
                pass
        if not isinstance(error, Exception):
            raise
        raise RegulatoryIndexError('required regulatory index is missing, stale or corrupt') from None


def _source_row_gap(content, left, right):
    """Count complete LF/CRLF/CR boundaries between requested row starts.

    Both ends are verified record starts. Native bytes.count avoids allocating
    positions for thousands of unrelated rows in every compressed chunk.
    """
    return (content.count(b'\n', left, right)+content.count(b'\r', left, right)
            -content.count(b'\r\n', left, right))


def query_regulatory_links(index, chromosome, start, end, *, maximum_elements=65536,
        maximum_source_rows=65536, deadline_seconds=30, selected_targets=None):
    """Half-open bounded window, complete source rows and deduplicated targets.

    An overlapping unjoined/missing-target element is retained with no invented
    target. An empty tuple means complete no element overlap in this window only.
    """
    _require(type(index) is RegulatoryIndex and type(chromosome) is str and chromosome
             and type(start) is int and type(end) is int and 0 <= start < end
             and end-start <= MAX_QUERY_BP and type(maximum_elements) is int and maximum_elements > 0
             and type(maximum_source_rows) is int and maximum_source_rows > 0
             and type(deadline_seconds) is int and 1 <= deadline_seconds <= 300)
    _require(selected_targets is None or type(selected_targets) is tuple
        and len(selected_targets)<=maximum_source_rows and selected_targets==tuple(sorted(set(selected_targets)))
        and all(type(v) is str and links.classify_encode_target(v)=='ensembl' for v in selected_targets))
    selected_set=None if selected_targets is None else frozenset(selected_targets)
    scope_sha=None if selected_targets is None else _sha(canonical(selected_targets))
    index.require_stable()
    deadline = time.monotonic()+deadline_seconds
    db = index.connection
    db.set_progress_handler(lambda: (db._read_cache.release(), int(time.monotonic() > deadline))[1], 1000)
    try:
        elements = db.execute('''SELECT DISTINCT e.ordinal,e.legacy_id,e.source_text FROM bins b
            JOIN elements e ON e.ordinal=b.ordinal WHERE b.chromosome=? AND b.bin BETWEEN ? AND ?
            AND e.start<? AND e.end>? ORDER BY e.ordinal LIMIT ?''',
            (chromosome,start//BIN_BP,(end-1)//BIN_BP,end,start,maximum_elements+1)).fetchall()
        _require(len(elements) <= maximum_elements)
        result, total_rows = [], 0
        for ordinal, element, element_text in elements:
            _require(time.monotonic() <= deadline)
            key = ordinal, element, element_text, scope_sha
            cached = index._elements.get(key)
            if cached is not None:
                hit, _ = cached
                total_rows += len(hit.source_rows)+len(hit.excluded_memberships)
                _require(total_rows <= maximum_source_rows)
                index._elements.move_to_end(key)
                result.append(hit)
                continue
            records, targets, unresolved_targets = [], set(), set()
            locators={}; excluded=[]
            members=db.execute('''SELECT gene_id,target_kind,source_rows FROM memberships
                WHERE legacy_id=? ORDER BY gene_id LIMIT ?''',(element,maximum_source_rows+1)).fetchall()
            _require(len(members)<=maximum_source_rows)
            for target, kind, count in members:
                _require(time.monotonic()<=deadline and kind==links.classify_encode_target(target)
                    and type(count) is int and count>0)
                if selected_set is not None and (kind!='ensembl' or target.split('.')[0] not in selected_set):
                    # The complete producer certifies membership and these row
                    # locators. No compressed biosample payload is needed for
                    # an explicitly omitted association or coarse gate fact.
                    first=db.execute('''SELECT c.first_ordinal,c.rows,m.source_offsets
                        FROM chunk_memberships m JOIN chunks c ON c.id=m.chunk_id
                        WHERE m.legacy_id=? AND m.gene_id=? ORDER BY m.chunk_id LIMIT 1''',
                        (element,target)).fetchone()
                    _require(first is not None and len(first[2])>=_SOURCE_OFFSET.size
                        and len(first[2])%_SOURCE_OFFSET.size==0)
                    offset,length,relative=_SOURCE_OFFSET.unpack_from(first[2])
                    _require(length>0 and relative<first[1])
                    excluded.append((target,kind,count,first[0]+relative))
                    total_rows+=1
                    _require(total_rows<=maximum_source_rows)
                    continue
                seen_count=0
                for chunk_id,number,positions in db.execute('''SELECT chunk_id,source_rows,source_offsets
                        FROM chunk_memberships WHERE legacy_id=? AND gene_id=? ORDER BY chunk_id''',
                        (element,target)):
                    _require(number>0 and len(positions)==number*_SOURCE_OFFSET.size)
                    seen_count+=number
                    _require(total_rows+seen_count<=maximum_source_rows)
                    locators.setdefault(chunk_id,[]).append((target,number,positions))
                _require(seen_count==count)
                total_rows+=count
            for chunk_id in sorted(locators):
                chunk=db.execute('SELECT first_ordinal,rows,bytes,sha256,compressed FROM chunks WHERE id=?',
                    (chunk_id,)).fetchone()
                _require(chunk is not None)
                first,rows,size,digest,compressed=chunk
                _require(time.monotonic() <= deadline and 0 < size <= index.manifest['maximum_chunk_bytes'])
                decoder = zlib.decompressobj()
                content = decoder.decompress(compressed, size+1)
                _require(len(content) == size and decoder.eof and not decoder.unused_data
                         and not decoder.unconsumed_tail and _sha(content) == digest)
                # Only requested offsets need ordinal validation. Count record
                # boundaries in disjoint byte ranges, preserving original mixed
                # newline semantics without walking every unrelated row.
                selected = locators[chunk_id]
                requested=[]
                for target, count, positions in selected:
                    _require(len(positions) == count*_SOURCE_OFFSET.size)
                    requested.extend((offset,length,relative,target)
                        for offset,length,relative in _SOURCE_OFFSET.iter_unpack(positions))
                seen = set(); previous_offset = previous_relative = 0
                for offset, length, relative, target in sorted(requested):
                    _require(time.monotonic() <= deadline and length>0 and offset+length<=size
                             and relative<rows and relative not in seen
                             and (offset==0 or content[offset-1:offset] in (b'\n',b'\r'))
                             and content[offset-1:offset+1] != b'\r\n')
                    observed = previous_relative+_source_row_gap(content,previous_offset,offset)
                    _require(observed==relative)
                    previous_offset,previous_relative=offset,relative
                    seen.add(relative)
                    raw = content[offset:offset+length]
                    _require(offset+length==size or raw.endswith((b'\n',b'\r')))
                    text = raw.decode('utf-8')
                    row = links.parse_encode_re2g_row(text)
                    _require(row['legacy_id']==element and row['source_target_id']==target)
                    records.append(RegulatorySourceRow(first+relative, text, _sha(raw)))
                    if row['ensembl_gene_id'] is not None:
                        targets.add(row['ensembl_gene_id'])
                    else:
                        unresolved_targets.add(row['source_target_id'])
            hit = RegulatoryElementLinks(ordinal, element_text.decode('utf-8'),
                tuple(sorted(targets)), tuple(sorted(records,key=lambda row:row.ordinal)), tuple(sorted(unresolved_targets)),
                tuple(excluded))
            index._remember_element(key, hit)
            result.append(hit)
        _require(time.monotonic() <= deadline)
        index.require_stable()
        return tuple(result)
    except RegulatoryIndexError:
        raise
    except Exception:
        raise RegulatoryIndexError('bounded regulatory query failed') from None
    finally:
        db.set_progress_handler(lambda: (db._read_cache.release(), 0)[1], 100000)


def indexed_regulatory_report(alleles, *, index, reference, selected_targets=None, metadata_only=False):
    """The original association/geometry semantics over the authenticated index.

    Tissue repetitions contribute retained source records to one membership,
    without multiplying members or providing independent regulatory function.
    Exhaust this iterator successfully before accepting any batch.
    """
    try:
        _require(type(index) is RegulatoryIndex)
        alleles = links._inputs(alleles)
        _require(type(metadata_only) is bool and not (metadata_only and selected_targets is not None))
        if metadata_only:selected_targets={a.allele_id:() for a in alleles}
        _require(selected_targets is None or type(selected_targets) is dict
            and set(selected_targets)=={a.allele_id for a in alleles})
        maps, fasta, ref_proof = index._reference(reference)
        _require(ref_proof == index.manifest['reference'])
        geometries, _, unresolved, windows = links.encode_query_geometries(alleles, maps, fasta)
        bindings = {row['resource_id']: row for row in index.manifest['sources']}
        element_binding = bindings[links.ENCODE_CCRE]
        link_binding = bindings[links.ENCODE_RE2G]
        source_bindings = [*index.manifest['sources'], dict(resource_id=SCHEMA,
            manifest_sha256=index.manifest_sha256, content_id=index.manifest['content_id'],
            index_sha256=index.manifest['index']['sha256'])]
        provenance = links._provenance('encode-registry-v4-re2g-index', source_bindings, ref_proof)
        overlapping, emitted = set(), set()
        for chromosome, queries in windows.items():
            for start, end, key in queries:
                allele, geometry = geometries[key]
                geometry_digest = links._digest(geometry)
                geometry_record = SourcedRecord('genosapiens.admission-geometry', 'current',
                    geometry_digest, geometry_digest, None, geometry)
                for left in range(start, end, MAX_QUERY_BP):
                    for hit in query_regulatory_links(index, chromosome, left, min(end, left+MAX_QUERY_BP),
                            selected_targets=None if selected_targets is None else selected_targets[allele.allele_id]):
                        marker = (key, hit.element_ordinal)
                        if marker in emitted:
                            continue
                        emitted.add(marker)
                        overlapping.add(allele.allele_id)
                        element = links.parse_encode_ccre_row(hit.element_text)
                        element_record = links._source_record(links.ENCODE_CCRE, 'Registry-V4',
                            element_binding['artifact']['path'], hit.element_ordinal,
                            element_binding['artifact']['sha256'], element_binding['manifest_sha256'], element)
                        by_target = {}
                        for source in hit.source_rows:
                            row = links.parse_encode_re2g_row(source.text)
                            by_target.setdefault((row['ensembl_gene_id'],row['source_target_id']), []).append(source)
                        def source_records(rows):
                            for source in rows:
                                row=links.parse_encode_re2g_row(source.text)
                                yield links._source_record(links.ENCODE_RE2G, 'Registry-V4',
                                    link_binding['artifact']['path'], source.ordinal,
                                    link_binding['artifact']['sha256'], link_binding['manifest_sha256'], row)
                        if not by_target and not hit.excluded_memberships:
                            yield links._proof(allele, 'regulatory_link', 'unresolved',
                                ('overlapping-element-without-matched-target-link',),
                                (element_record,), provenance)
                        for (target, _), rows in sorted(by_target.items(), key=lambda item: item[0][1]):
                            if selected_targets is not None:
                                from genosapiens.track1_association_scope import select_parsed_regulatory_sources
                                records,_=select_parsed_regulatory_sources(
                                    ((source.ordinal,links.parse_encode_re2g_row(source.text)) for source in rows),
                                    artifact_path=link_binding['artifact']['path'],
                                    source_sha256=link_binding['artifact']['sha256'],
                                    manifest_sha256=link_binding['manifest_sha256'])
                            else:
                                records=source_records(rows)
                            yield links._proof(allele, 'regulatory_link', 'admitted' if target else 'unresolved',
                                ('legacy-ccre-explicit-target-overlap' if target else 'source-target-unresolved',
                                 'context-only'), (element_record, *records, geometry_record), provenance, target)
                        for raw_target, kind, count, first_ordinal in hit.excluded_memberships:
                            from genosapiens.track1_association_scope import _summary
                            target=raw_target if kind=='ensembl' else None
                            membership=_summary('regulatory-membership-metadata', dict(
                                index_manifest_sha256=index.manifest_sha256,
                                element_ordinal=hit.element_ordinal, raw_target=raw_target,
                                target_kind=kind, source_rows=count, first_source_ordinal=first_ordinal))
                            yield links._proof(allele,'regulatory_link','admitted' if target else 'unresolved',
                                ('legacy-ccre-explicit-target-overlap' if target else 'source-target-unresolved',
                                 'context-only'),(element_record,membership,geometry_record),provenance,target)
        for allele in alleles:
            if allele.allele_id in unresolved:
                reasons, sources = unresolved[allele.allele_id]
                yield links._proof(allele, 'regulatory_link', 'unresolved', reasons, sources, provenance)
            elif allele.allele_id not in overlapping:
                yield links._proof(allele, 'regulatory_link', 'evaluated_no_match',
                    ('no-element-overlap-in-verified-registry-scope',), (), provenance)
        index.require_stable()
    except Exception:
        raise RegulatoryIndexError('indexed regulatory association report failed') from None
