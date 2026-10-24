"""Bounded, content-addressed logical evidence; no scientific decisions.

Version 1 stores JSON values as an ordered DAG. A reference is useful only in a
fully authenticated owner: all nodes precede their parents, every root is typed,
and every node must be reachable. SQLite owns uniqueness and closure, not an
in-memory dictionary. Repeated array entries are occurrences, not duplicates to
discard. This codec does not evaluate or certify scientific outcomes.
"""
from __future__ import annotations

from collections.abc import Mapping
from collections import OrderedDict
from contextlib import contextmanager
from dataclasses import dataclass, fields, is_dataclass
from enum import Enum
import hashlib
import json
import math
import os
from pathlib import Path
import re
import sqlite3
import stat
import tempfile
from types import MappingProxyType

from genosapiens.evidence_model import SourcedRecord
from genosapiens.evidence_serialization import _IDENTITIES, _loads
from genosapiens.allele_identity import identity_to_dict
from genosapiens.wgs_failures import WgsFailure
from genosapiens.wgs_ledger import MAX_MANIFEST_BYTES, MAX_SHARD_BYTES, MAX_SHARD_ROWS
from genosapiens.wgs_shards import ScientificShard

SCHEMA = "genosapiens.wgs-normalized-evidence"
VERSION = 1
# Independent downstream representation bounds. Historical execution-config
# bytes and their producer semantics are deliberately not changed.
FANOUT = 128
MAX_NODE_BYTES = 256 * 1024
STRING_CHARS = 16 * 1024
MAX_DEPTH = 64
MAX_LOGICAL_BYTES = 1024**4
MAX_OBJECT_FIELDS = 100_000
MAX_NODES = 2**31 - 1
MAX_ROOTS = 2**31 - 1
WRITER_CACHE_BYTES = 4 * 1024**2
WRITER_CACHE_ENTRIES = 8192
# A content key is built only for bounded source records. It never grants
# authority to an object identity or skips examination of the supplied tree.
SOURCE_KEY_BYTES = 1024 * 1024
SOURCE_KEY_ITEMS = 65536
LIVE_SOURCE_KEY_BYTES = 64 * 1024 * 1024
_SOURCE_FIELDS = frozenset(f.name for f in fields(SourcedRecord))
_LIVE_ARTIFACTS = OrderedDict()
_LIVE_ARTIFACT_LIMIT = 128
# Open handles only; node bytes and decoded values are never cached. Every
# scoped access fetches fresh bytes, checks their digest and decodes them.
READER_OPEN_HANDLES = 16
MAX_RECONSTRUCT_BYTES = 128 * 1024**2
# Phenotype audit ledgers retain exact integer numerators/denominators beyond
# int64. Keep them lossless within an explicit decimal bound below Python's
# default JSON integer digit limit, without changing process-global settings.
MAX_INTEGER_DIGITS = 4096
_MAX_INTEGER_ABS = 10**MAX_INTEGER_DIGITS
_SHA = re.compile(r"[0-9a-f]{64}")
_KIND = re.compile(r"[a-z][a-z0-9_.-]{0,95}")


class EvidenceStoreError(WgsFailure):
    def __init__(self, message="invalid normalized evidence", *, capacity=False,
                 operation="bounded_product_validate"):
        super().__init__("bounded_products", "response_capacity" if capacity else "malformed_evidence",
            operation=operation, exception_kind="project_failure")
        RuntimeError.__init__(self, message)


def _require(ok, message="invalid normalized evidence", *, capacity=False):
    if not ok:
        raise EvidenceStoreError(message, capacity=capacity)


def _bytes(value):
    return (json.dumps(value, sort_keys=True, separators=(",", ":"),
                       ensure_ascii=True, allow_nan=False) + "\n").encode("ascii")


def _digest(value):
    return type(value) is str and _SHA.fullmatch(value) is not None


def _integer(value, maximum=MAX_LOGICAL_BYTES):
    return type(value) is int and 0 <= value <= maximum


def _source_content_key(value, active, depth):
    """Exact normalized contents, or ordinary traversal for unsupported shapes.

    Inspect all values on every call, including read-only proxies' backing
    contents. Only exact built-in JSON shapes participate. The small accounting
    budget bounds the temporary plain tree and key; misses retain the ordinary
    writer's type, depth, cycle, integer and node-capacity checks.
    """
    budget = [SOURCE_KEY_BYTES, SOURCE_KEY_ITEMS]
    visiting = set(active)

    class Fallback(Exception):
        pass

    def plain(current, level):
        budget[1] -= 1
        budget[0] -= 64
        if level > MAX_DEPTH or budget[1] < 0 or budget[0] < 0:
            raise Fallback
        kind = type(current)
        if current is None or kind is bool:
            return current
        if kind is str:
            budget[0] -= 12 * len(current)
            if budget[0] < 0:
                raise Fallback
            return current
        if kind is int:
            budget[0] -= current.bit_length()
            if budget[0] < 0 or abs(current) >= _MAX_INTEGER_ABS:
                raise Fallback
            return current
        if kind is float:
            if not math.isfinite(current):
                raise Fallback
            return current
        if id(current) in visiting:
            raise Fallback
        visiting.add(id(current))
        try:
            if kind is SourcedRecord:
                entries = ((f.name, getattr(current, f.name)) for f in fields(current))
            elif kind in (dict, MappingProxyType):
                if len(current) > MAX_OBJECT_FIELDS:
                    raise Fallback
                entries = current.items()
            elif kind in (list, tuple):
                return [plain(child, level + 1) for child in current]
            else:
                raise Fallback
            result = {}
            for key, child in entries:
                if type(key) is not str or len(key) > STRING_CHARS or key in result:
                    raise Fallback
                budget[0] -= 128 + 12 * len(key)
                result[key] = plain(child, level + 1)
            return result
        finally:
            visiting.remove(id(current))

    try:
        content = _bytes(plain(value, depth))
        return content if len(content) <= SOURCE_KEY_BYTES else None
    except Fallback:
        return None


def _signature(value):
    return (value.st_dev, value.st_ino, value.st_size, value.st_mtime_ns, value.st_ctime_ns)


def _connect(directory):
    temporary = tempfile.TemporaryDirectory(prefix=".normalized-evidence-", dir=directory)
    connection = sqlite3.connect(Path(temporary.name) / "index.sqlite", uri=True)
    connection.execute("PRAGMA cache_size=-65536")
    connection.execute("PRAGMA temp_store=FILE")
    connection.execute("PRAGMA mmap_size=0")
    # This private index is deleted when its owner closes and is never saved
    # evidence. Erasing freed pages during intermediate table cleanup adds a
    # full write pass over large public-resource projections. Keep ordinary
    # transactions and rollback, but avoid that redundant temporary-page I/O.
    connection.execute("PRAGMA secure_delete=OFF")
    connection.execute("CREATE TABLE nodes (digest TEXT PRIMARY KEY, kind TEXT NOT NULL, "
        "size INTEGER NOT NULL, count INTEGER NOT NULL, depth INTEGER NOT NULL, "
        "first_key TEXT, last_key TEXT, shard INTEGER, offset INTEGER, length INTEGER) WITHOUT ROWID")
    connection.execute("CREATE TABLE edges (parent TEXT NOT NULL, child TEXT NOT NULL, "
        "PRIMARY KEY(parent,child)) WITHOUT ROWID")
    connection.execute("CREATE TABLE roots (ordinal INTEGER PRIMARY KEY, name TEXT UNIQUE NOT NULL, "
        "kind TEXT NOT NULL, digest TEXT NOT NULL)")
    connection.execute("CREATE INDEX roots_by_content ON roots(kind,digest)")
    return temporary, connection


def _borrow_prefix_index(source, destination, *, extra_tables=()):
    """Read an immutable live index with independent local additions.

    Only the attachment is shared. All inserts target the private database.
    The source owner must outlive this connection and remain unchanged.
    """
    source.commit()
    database = next(row[2] for row in source.execute('PRAGMA database_list') if row[1] == 'main')
    destination.execute('ATTACH DATABASE ? AS live_prefix',
        (Path(database).absolute().as_uri() + '?mode=ro',))
    for table in ('nodes', 'edges', 'roots'):
        columns = [row[1] for row in destination.execute('PRAGMA table_info(' + table + ')')]
        local = 'private_' + table
        destination.execute('ALTER TABLE ' + table + ' RENAME TO ' + local)
        destination.execute('CREATE TEMP VIEW ' + table + ' AS SELECT * FROM main.' + local
            + ' UNION ALL SELECT * FROM live_prefix.' + table)
        values = ','.join('NEW.' + name for name in columns)
        if table == 'edges':
            statement = ('INSERT OR IGNORE INTO ' + local + ' SELECT ' + values
                + ' WHERE NOT EXISTS(SELECT 1 FROM live_prefix.edges WHERE parent=NEW.parent AND child=NEW.child);')
        else:
            key = 'digest' if table == 'nodes' else 'name'
            statement = ("SELECT RAISE(ABORT,'duplicate shared index key') WHERE EXISTS(SELECT 1 FROM live_prefix."
                + table + ' WHERE ' + key + '=NEW.' + key + ');'
                + 'INSERT INTO ' + local + ' VALUES (' + values + ');')
        destination.execute('CREATE TEMP TRIGGER prefix_' + table + '_insert INSTEAD OF INSERT ON '
            + table + ' BEGIN ' + statement + ' END')
    for table in extra_tables:
        if table.startswith('science_resource_source_'):
            destination.execute('CREATE TEMP VIEW ' + table + ' AS SELECT * FROM live_prefix.' + table)
        else:
            schema = source.execute("SELECT sql FROM sqlite_master WHERE type='table' AND name=?", (table,)).fetchone()[0]
            destination.execute(schema)
            destination.execute('INSERT INTO ' + table + ' SELECT * FROM live_prefix.' + table)
    destination.commit()


def _copy_codec_index(source, destination):
    """Materialize a layered reader for the legacy independent-copy API."""
    layered = source.execute("SELECT 1 FROM sqlite_temp_master WHERE type='view' AND name='nodes'").fetchone()
    if not layered:
        source.backup(destination)
        return
    for table, width in (('nodes', 10), ('edges', 2), ('roots', 4)):
        cursor = source.execute('SELECT * FROM ' + table)
        try:
            destination.executemany('INSERT INTO ' + table + ' VALUES (' + ','.join('?' for _ in range(width)) + ')', cursor)
        finally:
            cursor.close()
    destination.commit()


@dataclass(frozen=True, slots=True, repr=False)
class EvidenceReference:
    """Typed root identity, bound to one complete manifest, never a trust flag."""
    name: str
    kind: str
    sha256: str
    manifest_sha256: str

    def __post_init__(self):
        _require(type(self.name) is str and 0 < len(self.name) <= 512
                 and type(self.kind) is str and _KIND.fullmatch(self.kind) is not None
                 and _digest(self.sha256) and _digest(self.manifest_sha256))


@dataclass(frozen=True, slots=True)
class SharedEvidencePrefix:
    sha256: str
    node_shards: int
    root_shards: int

    def __post_init__(self):
        _require(_digest(self.sha256) and _integer(self.node_shards, MAX_NODES)
            and _integer(self.root_shards, MAX_ROOTS) and self.node_shards > 0 and self.root_shards > 0)


@dataclass(frozen=True, slots=True, repr=False)
class EvidenceManifest:
    node_shards: tuple[ScientificShard, ...]
    root_shards: tuple[ScientificShard, ...]
    node_count: int
    root_count: int
    schema: str = SCHEMA
    version: int = VERSION
    shared_prefix: SharedEvidencePrefix | None = None

    def __post_init__(self):
        _require(self.schema == SCHEMA and type(self.version) is int
            and ((self.version == VERSION and self.shared_prefix is None)
                or (self.version == 2 and type(self.shared_prefix) is SharedEvidencePrefix
                    and self.shared_prefix.node_shards <= len(self.node_shards)
                    and self.shared_prefix.root_shards <= len(self.root_shards)))
            and _integer(self.node_count, MAX_NODES) and _integer(self.root_count, MAX_ROOTS)
            and self.node_count == sum(s.row_count for s in self.node_shards)
            and self.root_count == sum(s.row_count for s in self.root_shards)
            and bool(self.node_count) == bool(self.root_count))
        for kind, shards in (("evidence_nodes", self.node_shards), ("evidence_roots", self.root_shards)):
            _require(type(shards) is tuple and all(type(s) is ScientificShard for s in shards))
            for ordinal, shard in enumerate(shards, 1):
                _require(shard.kind == kind and shard.relative_path == f"{kind}-{ordinal:08d}.jsonl")
        self.as_bytes()

    def as_bytes(self):
        value = {field.name: getattr(self, field.name) for field in fields(self)}
        for key in ("node_shards", "root_shards"):
            value[key] = [{f.name: getattr(s, f.name) for f in fields(s)} for s in value[key]]
        if self.shared_prefix is None:
            del value['shared_prefix']
        else:
            value['shared_prefix'] = {f.name: getattr(self.shared_prefix, f.name) for f in fields(self.shared_prefix)}
        content = _bytes({"classification": "controlled", "llm_policy": "no-llm", "manifest": value})
        _require(len(content) <= MAX_MANIFEST_BYTES, "normalized manifest exceeds capacity", capacity=True)
        return content

    @property
    def sha256(self):
        return hashlib.sha256(self.as_bytes()).hexdigest()


def manifest_from_bytes(content, *, expected_sha256):
    _require(type(content) is bytes and len(content) <= MAX_MANIFEST_BYTES and _digest(expected_sha256)
             and hashlib.sha256(content).hexdigest() == expected_sha256)
    try:
        data = _loads(content)
        _require(set(data) == {"classification", "llm_policy", "manifest"}
                 and data["classification"] == "controlled" and data["llm_policy"] == "no-llm")
        value = data["manifest"]
        expected = {f.name for f in fields(EvidenceManifest)}
        if value.get('version') == VERSION:
            expected.remove('shared_prefix')
        _require(set(value) == expected)
        if 'shared_prefix' in value:
            value['shared_prefix'] = SharedEvidencePrefix(**value['shared_prefix'])
        for key in ("node_shards", "root_shards"):
            value[key] = tuple(ScientificShard(**s) for s in value[key])
        result = EvidenceManifest(**value)
        _require(result.as_bytes() == content)
        return result
    except (TypeError, ValueError, KeyError, AttributeError):
        raise EvidenceStoreError() from None


def _metadata(node, lookup):
    """Validate a single bounded node and derive exact expansion bounds bottom-up."""
    _require(type(node) is dict and set(node) == {"kind", "items"})
    kind, items = node["kind"], node["items"]
    children, first, last = [], None, None
    if kind == "scalar":
        _require(items is None or type(items) in (str, int, float, bool))
        _require(type(items) is not float or math.isfinite(items))
        _require(type(items) is not int or abs(items) < _MAX_INTEGER_ABS,
                 "normalized integer exceeds capacity",capacity=True)
        _require(type(items) is not str or len(items) <= STRING_CHARS)
        return ("string" if type(items) is str else "scalar", len(_bytes(items)) - 1, 1, 1, None, None), ()
    _require(kind in ("array", "object", "arrays", "objects", "strings")
             and type(items) is list and len(items) <= FANOUT)
    if kind in ("arrays", "objects", "strings"):
        _require(1 < len(items) <= FANOUT)
    size, count, depth = (2 if kind != "strings" else 2), 0, 1
    prior = None
    for item in items:
        if kind == "object":
            _require(type(item) is list and len(item) == 2 and type(item[0]) is str
                     and len(item[0]) <= STRING_CHARS)
            key, child = item
            _require(prior is None or prior < key)
            first = key if first is None else first
            last = prior = key
        else:
            child = item
        _require(_digest(child))
        row = lookup(child)
        _require(row is not None, "missing, forward or cyclic evidence dependency")
        child_kind, child_size, child_count, child_depth, child_first, child_last = row
        children.append(child)
        depth = max(depth, child_depth + 1)
        if kind == "object":
            size += len(_bytes(key)) - 1 + 1 + child_size + (1 if count else 0)
            count += 1
        elif kind == "array":
            size += child_size + (1 if count else 0)
            count += 1
        else:
            expected = {"arrays": "array", "objects": "object", "strings": "string"}[kind]
            _require(child_kind == expected and child_count > 0)
            if kind == "objects":
                _require(child_first is not None and (prior is None or prior < child_first))
                first = child_first if first is None else first
                last = prior = child_last
            size += child_size - 2 + (1 if count and kind != "strings" else 0)
            count += child_count
    _require(depth <= MAX_DEPTH and size <= MAX_LOGICAL_BYTES, "normalized expansion exceeds capacity", capacity=True)
    if kind in ("object", "objects"):
        _require(count <= MAX_OBJECT_FIELDS, "normalized object exceeds capacity", capacity=True)
    logical_kind = {"arrays": "array", "objects": "object", "strings": "string"}.get(kind, kind)
    return (logical_kind, size, count, depth, first, last), tuple(children)


class EvidenceWriter:
    """Append complete evaluated values; retain only bounded pages and DAG metadata."""
    def __init__(self, directory, *, scratch_directory=None, progress=None, _checkpoint_prefix=None,
                 _skip_verification=False):
        from genosapiens.wgs_shards import _ShardWriter
        self.directory = Path(directory)
        self.directory.mkdir(parents=True, exist_ok=True)
        self.temporary, self.connection = _connect(scratch_directory)
        self.nodes = _ShardWriter(self.directory, "evidence_nodes", max_rows=MAX_SHARD_ROWS, max_bytes=MAX_SHARD_BYTES)
        self.roots = _ShardWriter(self.directory, "evidence_roots", max_rows=MAX_SHARD_ROWS, max_bytes=MAX_SHARD_BYTES)
        self.node_count = self.root_count = 0
        self.manifest = None
        self.failed = False
        self._skip_verification = bool(_skip_verification)
        self._copied_prefix_source = None
        self._checkpoint_prefix = _checkpoint_prefix
        self._shared_prefix = None
        self._prefetched_metadata = None
        self._produced_indexes = []
        self._pending_reads = 0
        self._building = 0
        # Disposable context workers retain DAG metadata after retiring bytes.
        # Such an owner can never expose another read view or a saved manifest.
        self._fragmented = False
        self.progress = progress
        self._metadata_cache = OrderedDict()
        self._scalar_cache = OrderedDict()
        self._source_cache = OrderedDict()
        self._metadata_cache_bytes = self._scalar_cache_bytes = 0
        self._source_cache_bytes = 0
        self.metadata_cache_hits = self.scalar_cache_hits = 0
        self.source_cache_hits = 0

    def _cache(self, name, key, value, size):
        cache = getattr(self, name)
        byte_name = name + '_bytes'
        used = getattr(self, byte_name)
        if size > WRITER_CACHE_BYTES or WRITER_CACHE_ENTRIES == 0:
            return
        if key in cache:
            cache.move_to_end(key)
            return
        while cache and (used + size > WRITER_CACHE_BYTES or len(cache) >= WRITER_CACHE_ENTRIES):
            _, (_, old_size) = cache.popitem(last=False)
            used -= old_size
        cache[key] = (value, size)
        setattr(self, byte_name, used + size)

    def _remember_metadata(self, digest, row):
        if self._prefetched_metadata is not None and digest in self._prefetched_metadata:
            self._prefetched_metadata[digest] = row
        # Conservative resident accounting includes tuple/dict/string overhead;
        # first/last object keys can be long Unicode strings.
        self._cache('_metadata_cache', digest, row,
            1024 + 4 * sum(len(v) for v in row if type(v) is str))

    def _lookup(self, digest):
        if self._prefetched_metadata is not None and digest in self._prefetched_metadata:
            return self._prefetched_metadata[digest]
        cached = self._metadata_cache.get(digest)
        if cached is not None:
            self._metadata_cache.move_to_end(digest)
            self.metadata_cache_hits += 1
            return cached[0]
        row = self.connection.execute("SELECT kind,size,count,depth,first_key,last_key FROM nodes WHERE digest=?", (digest,)).fetchone()
        if row is not None:
            self._remember_metadata(digest, row)
        return row

    def _scalar_tuple_batch(self, rows):
        """Batch codec-index lookups for bounded three-string source rows.

        This changes no source authority or canonical node encoding. All new
        nodes still pass _node and _metadata. The temporary lookup covers only
        this synchronous batch and clears on every exit, including failure.
        """
        _require(type(rows) in (tuple, list) and len(rows) <= FANOUT
            and self._prefetched_metadata is None and not self.failed)
        if FANOUT < 3 or not all(type(row) is tuple and len(row) == 3
                and all(type(value) is str and len(value) <= STRING_CHARS for value in row) for row in rows):
            return [self._value(row, set(), 2) for row in rows]
        if self._skip_verification and not self._fragmented:
            return self._prepared_scalar_tuple_batch(rows)
        digests = {}
        for row in rows:
            children = [hashlib.sha256(_bytes({'kind': 'scalar', 'items': value})).hexdigest() for value in row]
            digests.update((digest, None) for digest in children)
            digests[hashlib.sha256(_bytes({'kind': 'array', 'items': children})).hexdigest()] = None
        if not digests:
            return []
        query = 'SELECT digest,kind,size,count,depth,first_key,last_key FROM nodes WHERE digest IN ('
        query += ','.join('?' for _ in digests) + ')'
        for digest, *metadata in self.connection.execute(query, tuple(digests)):
            digests[digest] = tuple(metadata)
        self._prefetched_metadata = digests
        try:
            return [self._node('array', [self._node('scalar', value) for value in row]) for row in rows]
        except BaseException:
            self.failed = True
            raise
        finally:
            self._prefetched_metadata = None

    def _prepared_scalar_tuple_batch(self, rows):
        """Encode each bounded node once and batch its private index writes."""
        from genosapiens.wgs_shards import CanonicalRow, _CANONICAL_ROW_SEAL
        nodes, result = {}, []
        def prepare(kind, items):
            node = dict(kind=kind, items=items)
            content = _bytes(node)
            _require(len(content) <= MAX_NODE_BYTES, 'normalized node exceeds capacity', capacity=True)
            digest = hashlib.sha256(content).hexdigest()
            nodes.setdefault(digest, (node, content))
            return digest
        for row in rows:
            result.append(prepare('array', [prepare('scalar', value) for value in row]))
        if not nodes:
            return result
        known = {digest: tuple(metadata) for digest, *metadata in self.connection.execute(
            'SELECT digest,kind,size,count,depth,first_key,last_key FROM nodes WHERE digest IN ('
            + ','.join('?' for _ in nodes) + ')', tuple(nodes))}
        additions, edges = [], []
        try:
            for digest, (node, content) in nodes.items():
                if digest in known:
                    continue
                metadata, children = _metadata(node, known.get)
                _require(self.node_count < MAX_NODES, 'normalized node count exceeds capacity', capacity=True)
                row = dict(sha256=digest, node=node)
                encoded = CanonicalRow(row,
                    b'{"node":' + content[:-1] + b',"sha256":"' + digest.encode('ascii') + b'"}\n',
                    _CANONICAL_ROW_SEAL)
                self.nodes.write(row, key=digest, canonical_bytes=encoded)
                additions.append((digest, *metadata, len(self.nodes.descriptors),
                    self.nodes._size-len(encoded.content), len(encoded.content)))
                edges.extend((digest, child) for child in children)
                known[digest] = metadata
                self._remember_metadata(digest, metadata)
                self.node_count += 1
                self._metadata_budget()
                if self.progress is not None and self.node_count % 4096 == 0:
                    self.progress(self.node_count)
            self.connection.executemany('INSERT INTO nodes VALUES (?,?,?,?,?,?,?,?,?,?)', additions)
            self.connection.executemany('INSERT OR IGNORE INTO edges VALUES (?,?)', edges)
            return result
        except BaseException:
            self.failed = True
            raise

    def _node(self, kind, items):
        _require(not self._fragmented or self._fragmented.active is not None,
                 "disposable owner has no active fragment")
        node = {"kind": kind, "items": items}
        content = _bytes(node)
        _require(len(content) <= MAX_NODE_BYTES, "normalized node exceeds capacity", capacity=True)
        digest = hashlib.sha256(content).hexdigest()
        if self._lookup(digest) is None:
            metadata, children = self._node_metadata(node)
            _require(self.node_count < MAX_NODES, "normalized node count exceeds capacity", capacity=True)
            from genosapiens.wgs_shards import CanonicalRow, _CANONICAL_ROW_SEAL
            row = {"sha256": digest, "node": node}
            # Reuse this owner's already canonical, validated plain node bytes.
            # Sorted wrapper fields are node then sha256; no second recursive
            # plain conversion/JSON encoding of the node is necessary.
            encoded = CanonicalRow(row,
                b'{"node":' + content[:-1] + b',"sha256":"' + digest.encode('ascii') + b'"}\n',
                _CANONICAL_ROW_SEAL)
            if self._fragmented:
                self._fragmented.charge(len(encoded.content))
            self.nodes.write(row, key=digest, canonical_bytes=encoded)
            self.connection.execute("INSERT INTO nodes VALUES (?,?,?,?,?,?,?,?,?,?)",
                (digest, *metadata, len(self.nodes.descriptors),
                    self.nodes._size-len(encoded.content), len(encoded.content)))
            self.connection.executemany("INSERT OR IGNORE INTO edges VALUES (?,?)", ((digest, child) for child in children))
            self._remember_metadata(digest, metadata)
            self.node_count += 1
            self._metadata_budget()
            if self.progress is not None and self.node_count % 4096 == 0:
                self.progress(self.node_count)
        return digest

    def _node_metadata(self, node):
        """Resolve wide nodes with one bounded metadata query.

        Ordinary metadata validation still checks every child and derives all
        expansion bounds. The temporary map cannot outlive this node build.
        Small nodes and nested prepared batches keep their existing lookup.
        """
        kind, items = node['kind'], node['items']
        if (self._prefetched_metadata is not None or type(items) is not list
                or not 16 <= len(items) <= FANOUT
                or kind not in ('array', 'object', 'arrays', 'objects', 'strings')):
            return _metadata(node, self._lookup)
        if kind == 'object':
            if not all(type(item) is list and len(item) == 2 for item in items):
                return _metadata(node, self._lookup)
            children = [item[1] for item in items]
        else:
            children = items
        if not all(_digest(child) for child in children):
            return _metadata(node, self._lookup)
        available, missing = {}, set()
        for child in children:
            cached = self._metadata_cache.get(child)
            if cached is None:
                missing.add(child)
            else:
                available[child] = cached[0]
        if missing:
            for digest, *metadata in self.connection.execute(
                    'SELECT digest,kind,size,count,depth,first_key,last_key FROM nodes WHERE digest IN ('
                    + ','.join('?' for _ in missing) + ')', tuple(missing)):
                available[digest] = tuple(metadata)
        return _metadata(node, available.get)

    def _metadata_budget(self):
        _require(self.nodes.descriptor_bytes + self.roots.descriptor_bytes + 2048 <= MAX_MANIFEST_BYTES,
                 "normalized manifest exceeds capacity", capacity=True)

    def _chunks(self, values, kind):
        # At most FANOUT nodes at each depth; no list proportional to the array.
        levels = [[]]
        def carry(value, level=0):
            if level == len(levels):
                levels.append([])
            levels[level].append(value)
            if len(levels[level]) == FANOUT:
                merged = self._node(kind, levels[level])
                levels[level] = []
                carry(merged, level + 1)
        for value in values:
            carry(value)
        result = None
        for level in range(len(levels)):
            values = levels[level]
            if result is not None:
                values = [*values, result]
            if values:
                result = values[0] if len(values) == 1 else self._node(kind, values)
        _require(result is not None)
        return result

    def _value(self, value, active, depth):
        # Intermediate scientific tapes contain plain decoded source records;
        # their complete contents normalize identically to the typed form.
        source = (type(value) is SourcedRecord or type(value) is dict
            and 'source_sha256' in value and value.keys() == _SOURCE_FIELDS)
        content = None
        if source and WRITER_CACHE_BYTES and WRITER_CACHE_ENTRIES:
            if self._skip_verification:
                # Fresh scientific records already belong to their producer.
                # Avoid another recursive validation walk and retain only a
                # fixed-size content fingerprint, not a large canonical key.
                from genosapiens.evidence_serialization import _owned_canonical_bytes
                content = _owned_canonical_bytes(value)
                if len(content) > LIVE_SOURCE_KEY_BYTES:
                    content = None
            else:
                content = _source_content_key(value, active, depth)
        if content is None:
            return self._value_uncached(value, active, depth)
        key = (len(content), hashlib.sha256(content).digest()) if self._skip_verification else content
        cached = self._source_cache.get(key)
        if cached is not None:
            _require(not self._fragmented or self._fragmented.active is not None,
                     "disposable owner has no active fragment")
            self._source_cache.move_to_end(key)
            self.source_cache_hits += 1
            return cached[0]
        # Build from the same captured contents used by the key. A proxy may
        # wrap a custom mapping whose items change between traversals; reading
        # the original a second time could associate the key with another DAG.
        digest = self._value_uncached(_loads(content), active, depth)
        self._cache('_source_cache', key, digest, 1088 if self._skip_verification else 1024 + len(content))
        return digest

    def _value_uncached(self, value, active, depth):
        _require(depth <= MAX_DEPTH, "normalized input depth exceeds capacity", capacity=True)
        if isinstance(value, Enum):
            value = value.value
        if type(value) in _IDENTITIES.values():
            value = identity_to_dict(value)
        if value is None or type(value) in (str, bool, int, float):
            _require(type(value) is not int or abs(value) < _MAX_INTEGER_ABS,
                     "normalized integer exceeds capacity",capacity=True)
            if type(value) is str and len(value) > STRING_CHARS:
                return self._chunks((self._node("scalar", value[start:start + STRING_CHARS])
                    for start in range(0, len(value), STRING_CHARS)), "strings")
            # Exact immutable scalar keys only. Floats deliberately use the
            # uncached path: Python equality conflates -0.0 and +0.0, whereas
            # their canonical evidence bytes differ. No mutable-tree cache.
            key = (type(value), value) if type(value) is not float else None
            cached = None if key is None else self._scalar_cache.get(key)
            if cached is not None:
                self._scalar_cache.move_to_end(key)
                self.scalar_cache_hits += 1
                return cached[0]
            digest = self._node("scalar", value)
            if key is not None:
                size = 1024 + (4 * len(value) if type(value) is str else
                    (value.bit_length() + 7)//8 if type(value) is int else 0)
                self._cache('_scalar_cache', key, digest, size)
            return digest
        _require(id(value) not in active, "cyclic normalized input")
        active.add(id(value))
        try:
            if is_dataclass(value) and not isinstance(value, type):
                entries = [(f.name, getattr(value, f.name)) for f in fields(value)]
            elif isinstance(value, Mapping):
                _require(len(value) <= MAX_OBJECT_FIELDS, "normalized object exceeds capacity", capacity=True)
                entries = list(value.items())
            else:
                _require(type(value) in (list, tuple), "unsupported normalized evidence type")
                return self._sequence(iter(value), active, depth)
            _require(all(type(key) is str and len(key) <= STRING_CHARS for key, _ in entries))
            entries.sort(key=lambda pair: pair[0])
            def chunks():
                for start in range(0, len(entries), FANOUT):
                    yield self._node("object", [[key, self._value(child, active, depth + 1)]
                        for key, child in entries[start:start + FANOUT]])
            return self._chunks(chunks(), "objects") if entries else self._node("object", [])
        finally:
            active.remove(id(value))

    def _sequence(self, values, active, depth):
        def chunks():
            batch = []
            for value in values:
                batch.append(self._value(value, active, depth + 1))
                if len(batch) == FANOUT:
                    yield self._node("array", batch)
                    batch = []
            if batch:
                yield self._node("array", batch)
        iterator = iter(chunks())
        first = next(iterator, None)
        if first is None:
            return self._node("array", [])
        from itertools import chain
        return self._chunks(chain((first,), iterator), "arrays")

    def add(self, name, kind, value):
        return self._add(name, kind, lambda: self._value(value, set(), 1))

    def add_sequence(self, name, kind, values):
        """Consume a complete ordered stream; exhaustion precedes root publication."""
        return self._add(name, kind, lambda: self._sequence(iter(values), set(), 1))

    def add_alias(self, name, kind, original):
        """Bind another name to a completed root in this producing owner."""
        row = self.connection.execute('SELECT kind,digest FROM roots WHERE name=?', (original,)).fetchone()
        _require(row is not None and row[0] == kind, 'shared evidence root is foreign')
        return self._add(name, kind, lambda: row[1])

    def add_sequence_prefix(self, name, kind, values, *, prefix=None, prefix_values=None):
        """Append values to a shared nonempty array without walking its values.

        Both components remain ordinary v1 array DAGs. A newly built prefix is
        reached by this root, so no extra root or unreachable node is emitted.
        The caller owns logical ordering; this codec preserves every occurrence.
        """
        from itertools import chain
        _require((prefix is None) != (prefix_values is None))
        parts = []
        def build():
            base = prefix if prefix is not None else self._sequence(iter(prefix_values), set(), 1)
            metadata = self._lookup(base)
            _require(metadata is not None and metadata[0] == 'array' and metadata[2] > 0)
            parts.append(base)
            iterator = iter(values)
            sentinel = object()
            first = next(iterator, sentinel)
            if first is not sentinel:
                parts.append(self._sequence(chain((first,), iterator), set(), 1))
            return parts[0] if len(parts) == 1 else self._node('arrays', parts)
        digest = self._add(name, kind, build)
        return digest, tuple(parts)

    def _add(self, name, kind, build):
        _require(self.connection is not None and self.manifest is None and not self.failed)
        _require(not self._fragmented or self._fragmented.active is not None,
                 "disposable owner has no active fragment")
        self._building += 1
        try:
            EvidenceReference(name, kind, "0" * 64, "0" * 64)
            _require(self.root_count < MAX_ROOTS, "normalized root count exceeds capacity", capacity=True)
            digest = build()
            ordinal = self.root_count + 1
            self.connection.execute("INSERT INTO roots VALUES (?,?,?,?)", (ordinal, name, kind, digest))
            row = {"ordinal": ordinal, "name": name, "kind": kind, "sha256": digest}
            if self._fragmented:
                self._fragmented.charge(len(_bytes(row)))
            self.roots.write(row, key=name)
            self.root_count = ordinal
            self._metadata_budget()
            return digest
        except BaseException:
            self.failed = True
            raise
        finally:
            self._building -= 1

    def copy_from(self, reader):
        """Copy an authenticated DAG in bounded node pages into a fresh owner."""
        _require(type(reader) is EvidenceReader and reader.verified
                 and self.node_count == self.root_count == 0)
        with reader.read_into(self):
            for number, shard in enumerate(reader.manifest.node_shards):
                for row, _, _ in reader._rows(shard, number):
                    _require(self._node(row["node"]["kind"], row["node"]["items"]) == row["sha256"])
            for _, name, kind, digest in reader.connection.execute("SELECT * FROM roots ORDER BY ordinal"):
                self._add(name, kind, lambda: digest)
            _require((self.node_count, self.root_count) == (reader.manifest.node_count, reader.manifest.root_count))

    def merge_from(self, reader, *, exclude_prefix=None):
        """Merge authenticated roots, rejecting conflicting shared names.

        A context artifact may also hold local case projections. The coordinator
        regenerates those projections for its own domain. Excluding their root
        prefix copies only nodes reachable from the remaining roots, so unused
        projection nodes cannot enter the new evidence manifest.
        """
        _require(type(reader) is EvidenceReader and reader.verified
            and reader.connection is not self.connection and not self._fragmented
            and (exclude_prefix is None or type(exclude_prefix) is str and 0 < len(exclude_prefix) <= 256))
        created = False
        try:
            with reader.read_into(self):
                reader.connection.execute('CREATE TEMP TABLE merge_nodes(digest TEXT PRIMARY KEY) WITHOUT ROWID')
                created = True
                def included(name):
                    return exclude_prefix is None or not name.startswith(exclude_prefix)
                shared = (self._copied_prefix_source is not None
                    and self._copied_prefix_source is reader._live_prefix_source
                    and reader.connection.execute("SELECT 1 FROM sqlite_temp_master WHERE type='view' AND name='nodes'").fetchone())
                node_table = 'main.private_nodes' if shared else 'nodes'
                edge_table = 'main.private_edges' if shared else 'edges'
                for name, digest in reader.connection.execute('SELECT name,digest FROM roots ORDER BY ordinal'):
                    if included(name) and (not shared or reader.connection.execute(
                            'SELECT 1 FROM main.private_nodes WHERE digest=?', (digest,)).fetchone()):
                        reader.connection.execute('INSERT OR IGNORE INTO merge_nodes VALUES (?)', (digest,))
                for _ in range(MAX_DEPTH):
                    before = reader.connection.total_changes
                    reader.connection.execute('INSERT OR IGNORE INTO merge_nodes SELECT e.child FROM '
                        + edge_table + ' e JOIN merge_nodes m ON e.parent=m.digest JOIN '
                        + node_table + ' n ON n.digest=e.child')
                    if reader.connection.total_changes == before:
                        break
                # The reader authenticated every source byte before this scope.
                # Its boundaries still hash all shards, including excluded roots.
                # Only new reached nodes need decoding and writer insertion.
                # Physical order preserves the codec's child-before-parent rule.
                rows = reader.connection.execute('SELECT n.digest FROM ' + node_table + ' n '
                    'JOIN merge_nodes m ON n.digest=m.digest ORDER BY n.shard,n.offset')
                try:
                    while batch := rows.fetchmany(512):
                        # Existing nodes need an exact membership lookup, not
                        # their full metadata or another LRU insertion. Keep
                        # query parameters and retained rows independently
                        # bounded while preserving physical insertion order.
                        existing = {row[0] for row in self.connection.execute(
                            'SELECT digest FROM nodes WHERE digest IN ('
                            + ','.join('?' for _ in batch) + ')',
                            tuple(row[0] for row in batch))}
                        for (digest,) in batch:
                            if digest not in existing:
                                node = reader._node(digest)
                                _require(self._node(node['kind'], node['items']) == digest)
                finally:
                    rows.close()
                for _, name, kind, digest in reader.connection.execute('SELECT * FROM roots ORDER BY ordinal'):
                    if not included(name):
                        continue
                    previous = self.connection.execute('SELECT kind,digest FROM roots WHERE name=?', (name,)).fetchone()
                    if previous is None:
                        self._add(name, kind, lambda: digest)
                    else:
                        _require(previous == (kind, digest), 'conflicting shared scientific evidence root')
        except BaseException:
            self.failed = True
            raise
        finally:
            if created and reader.connection is not None:
                reader.connection.execute('DROP TABLE IF EXISTS merge_nodes')

    def copy_authenticated_prefix(self, reader):
        """Copy a live authenticated prefix as one complete input operation."""
        with self.authenticated_prefix_scope(reader):
            pass

    def copy_produced_prefix(self, source, *, extra_tables=()):
        """Copy this operation's live immutable resource prefix, not a saved index.

        Source hashes bracket copying. Each derivative owns separate shard
        files and an independent SQLite index, including its policy lookup
        tables. A later derivative failure cannot poison the producing owner.
        """
        from genosapiens.wgs_shards import _bytes as shard_bytes
        import re
        _require(type(source) is EvidenceWriter and source is not self
            and source.connection is not None and not source.failed and source.manifest is None
            and not source._fragmented and source._pending_reads == source._building == 0
            and self.connection is not None and not self.failed and self.manifest is None
            and not self._fragmented and self.node_count == self.root_count == 0
            and self._pending_reads == self._building == 0
            and not self.nodes.descriptors and not self.roots.descriptors
            and type(extra_tables) is tuple and len(extra_tables) <= 4
            and all(type(name) is str and (re.fullmatch(r'population_[0-9a-f]{64}_(contexts|exceptions)', name)
                    or name in ('science_resource_source_owners', 'science_resource_source_records'))
                    for name in extra_tables))
        tables = lambda c: {r[0] for r in c.execute("SELECT name FROM sqlite_master WHERE type='table'")}
        _require(tables(source.connection) == {'nodes', 'edges', 'roots', *extra_tables}
            and tables(self.connection) == {'nodes', 'edges', 'roots'})
        source.connection.commit()
        self._copied_prefix_source = source
        if self._checkpoint_prefix is not None:
            try:
                from genosapiens.wgs_evidence_shared import CheckpointResourcePrefix
                _require(type(self._checkpoint_prefix) is CheckpointResourcePrefix)
                self._checkpoint_prefix.require(source, self.directory)
                self._shared_prefix = self._checkpoint_prefix.descriptor
                # The bundle already contains this live producer's exact
                # prefix. Reuse its stream states and descriptors instead of
                # reading every shared byte again for each new context.
                if self._skip_verification:
                    _borrow_prefix_index(source.connection, self.connection, extra_tables=extra_tables)
                else:
                    source.connection.backup(self.connection)
                for output, original in ((self.nodes, source.nodes), (self.roots, source.roots)):
                    output.descriptors.extend(original.descriptors)
                    output.descriptor_bytes = original.descriptor_bytes
                    output.metrics['shards'] = len(original.descriptors)
                    output.stream = original.stream.copy()
                self.node_count, self.root_count = source.node_count, source.root_count
                self._metadata_budget()
                return
            except BaseException:
                self.failed = True
                raise
        self._pending_reads += 1
        try:
            copy_error = None
            with source.read_view(_producer_owned=source._skip_verification) as reader:
                try:
                    for output, shards in ((self.nodes, reader.manifest.node_shards),
                                           (self.roots, reader.manifest.root_shards)):
                        for shard in shards:
                            digest, size = hashlib.sha256(), 0
                            with reader._open(shard) as handle:
                                before = _signature(os.fstat(handle.fileno()))
                                from contextlib import nullcontext
                                if self._shared_prefix is None:
                                    fd = os.open(self.directory / shard.relative_path,
                                        os.O_WRONLY | os.O_CREAT | os.O_EXCL | os.O_NOFOLLOW, 0o600)
                                    target_scope = os.fdopen(fd, 'wb')
                                else:
                                    _require(not (self.directory / shard.relative_path).exists())
                                    target_scope = nullcontext(None)
                                with target_scope as target:
                                    for block in iter(lambda: handle.read(1024 * 1024), b''):
                                        size += len(block)
                                        _require(size <= shard.size_bytes)
                                        digest.update(block)
                                        output.stream.update(block)
                                        if target is not None:
                                            target.write(block)
                                    _require(size == shard.size_bytes and digest.hexdigest() == shard.sha256
                                        and before == _signature(os.fstat(handle.fileno())))
                                    if target is not None:
                                        target.flush()
                                        os.fsync(target.fileno())
                            output.descriptors.append(shard)
                            output.descriptor_bytes += len(shard_bytes(shard))
                            output.metrics['shards'] += 1
                    source.connection.backup(self.connection)
                    self.node_count, self.root_count = source.node_count, source.root_count
                    self._metadata_budget()
                except BaseException as error:
                    copy_error = error
            if copy_error is not None:
                raise copy_error
        except BaseException:
            self.failed = True
            raise
        finally:
            self._pending_reads -= 1

    @contextmanager
    def authenticated_prefix_scope(self, reader):
        """Copy a live prefix and keep its input scope open for the derivative.

        This accepts no saved index or verification marker. The ordinary
        reader has just authenticated the complete DAG; its committed index
        must contain only codec tables. Full source hashes bracket independent
        shard copies, and any failure poisons this new append-only owner.
        """
        from genosapiens.wgs_shards import _bytes as shard_bytes
        _require(type(reader) is EvidenceReader and reader.verified
                 and reader.connection is not None and reader._source_writer is None
                 and reader._read_scope_depth == 0 and not reader._root_coverage
                 and not reader.connection.in_transaction
                 and self.connection is not None and self.connection is not reader.connection
                 and not self.connection.in_transaction and not self.failed and not self._fragmented
                 and self.manifest is None and self.node_count == self.root_count == 0
                 and self._pending_reads == self._building == 0
                 and not self.nodes.descriptors and not self.roots.descriptors)
        expected_tables = {"nodes", "edges", "roots"}
        for connection in (reader.connection, self.connection):
            _require({row[0] for row in connection.execute(
                "SELECT name FROM sqlite_master WHERE type='table'")} == expected_tables,
                "normalized prefix owner contains foreign projections")
        with reader.read_into(self):
            for output, shards in ((self.nodes, reader.manifest.node_shards),
                                   (self.roots, reader.manifest.root_shards)):
                for shard in shards:
                    digest, size = hashlib.sha256(), 0
                    with reader._open(shard) as source:
                        before = _signature(os.fstat(source.fileno()))
                        fd = os.open(self.directory / shard.relative_path,
                            os.O_WRONLY | os.O_CREAT | os.O_EXCL | os.O_NOFOLLOW, 0o600)
                        with os.fdopen(fd, "wb") as target:
                            for block in iter(lambda: source.read(1024 * 1024), b""):
                                size += len(block)
                                _require(size <= shard.size_bytes)
                                digest.update(block)
                                output.stream.update(block)
                                target.write(block)
                            _require(size == shard.size_bytes and digest.hexdigest() == shard.sha256
                                     and before == _signature(os.fstat(source.fileno())))
                            target.flush()
                            os.fsync(target.fileno())
                    output.descriptors.append(shard)
                    output.descriptor_bytes += len(shard_bytes(shard))
                    output.metrics["shards"] += 1
            _copy_codec_index(reader.connection, self.connection)
            self.node_count, self.root_count = reader.manifest.node_count, reader.manifest.root_count
            self._metadata_budget()
            # Publication remains blocked until the caller's derived work and
            # the final full source-byte authentication both succeed.
            yield self

    @contextmanager
    def read_view(self, *, _producer_owned=False):
        """Private read scope over this operation's already completed roots.

        Close current shards so later additions cannot change the scope's
        files. Reuse the producer's computed DAG index, not a saved cache or a
        wire validation marker. Every accessed node is checked against its
        digest; full shard bytes are authenticated on both scope boundaries.
        Publication remains blocked until the scope completes successfully.
        The ordinary public reader still authenticates the final full artifact.
        """
        _require(self.connection is not None and self.manifest is None and not self.failed
                 and not self._fragmented and self._building == 0 and self._pending_reads == 0)
        self.nodes._close()
        self.roots._close()
        snapshot = self._snapshot_manifest()
        reader = EvidenceReader(snapshot,self.directory)
        reader._producer_owned_read = bool(_producer_owned)
        reader._skip_verification = self._skip_verification
        reader.verification_performed = not self._skip_verification
        reader.connection,reader.temporary,reader._source_writer = self.connection,self.temporary,self
        try:
            for shard in (*snapshot.node_shards,*snapshot.root_shards):
                with reader._open(shard) as handle:
                    reader.signatures.append((shard,_signature(os.fstat(handle.fileno()))))
            reader.verified = True
            with reader.read_into(self):
                yield reader
        except BaseException:
            self.failed = True
            raise
        finally:
            reader.__exit__()

    def produced_prefix(self, *, extra_tables=()):
        """Retain this live producer's index for byte-checked context readers."""
        result = _ProducedEvidencePrefix(self, extra_tables=extra_tables)
        self._produced_indexes.append(result)
        return result

    def finish(self):
        _require(self.connection is not None and self.manifest is None and not self.failed)
        _require(not self._fragmented, "disposable fragment owner cannot publish a manifest")
        _require(self._pending_reads == 0 and self._building == 0, "normalized input verification is incomplete")
        try:
            self.nodes.close()
            self.roots.close()
            self.manifest = self._snapshot_manifest()
            if self._skip_verification:
                key = (os.getpid(), str(self.directory.absolute()), self.manifest.sha256)
                _LIVE_ARTIFACTS[key] = None
                _LIVE_ARTIFACTS.move_to_end(key)
                while len(_LIVE_ARTIFACTS) > _LIVE_ARTIFACT_LIMIT:
                    _LIVE_ARTIFACTS.popitem(last=False)
            return self.manifest
        finally:
            self.close()

    def _snapshot_manifest(self):
        return EvidenceManifest(tuple(self.nodes.descriptors), tuple(self.roots.descriptors),
            self.node_count, self.root_count, version=1 if self._shared_prefix is None else 2,
            shared_prefix=self._shared_prefix)

    def close(self):
        self.nodes._close()
        self.roots._close()
        for prefix in self._produced_indexes:
            prefix.close()
        self._produced_indexes.clear()
        if self.connection is not None:
            self.connection.close()
            self.connection = None
            self.temporary.cleanup()
        self._metadata_cache.clear()
        self._scalar_cache.clear()
        self._source_cache.clear()
        self._metadata_cache_bytes = self._scalar_cache_bytes = 0
        self._source_cache_bytes = 0

    def __enter__(self):
        return self

    def __exit__(self, *_):
        self.close()


def _require_complete_dag_roots(connection, node_count, *, retain_unreferenced=False):
    """Check root coverage after unique, child-before-parent insertion.

    In this finite DAG, every node must be a root or somebody's child. If an
    unreachable component existed, its last inserted node would have neither
    a root nor a later parent. This avoids revisiting all reached edges at
    every depth. It relies on _metadata rejecting missing/forward children
    before insertion and on root targets and node uniqueness being checked.
    """
    connection.execute("CREATE TABLE root_coverage (digest TEXT PRIMARY KEY) WITHOUT ROWID")
    try:
        connection.execute("INSERT OR IGNORE INTO root_coverage SELECT digest FROM roots")
        connection.execute("INSERT OR IGNORE INTO root_coverage SELECT child FROM edges")
        if retain_unreferenced:
            # A live prefix can contain prepared nodes that acquire roots only
            # in the context suffix. Preserve that exact unresolved frontier.
            connection.execute('CREATE TABLE suffix_unreferenced(digest TEXT PRIMARY KEY) WITHOUT ROWID')
            connection.execute('INSERT INTO suffix_unreferenced SELECT n.digest FROM nodes n '
                'LEFT JOIN root_coverage c ON c.digest=n.digest WHERE c.digest IS NULL')
        else:
            _require(connection.execute("SELECT count(*) FROM root_coverage").fetchone()[0] == node_count,
                     "unreferenced normalized evidence")
    finally:
        connection.execute("DROP TABLE root_coverage")


class _ProducedEvidencePrefix:
    """Private live authority, never serialized or reconstructed from a cache."""
    def __init__(self, writer, *, extra_tables):
        import re
        _require(type(writer) is EvidenceWriter and writer.connection is not None
            and not writer.failed and writer.manifest is None and not writer._fragmented
            and writer._pending_reads == writer._building == 0
            and type(extra_tables) is tuple and len(extra_tables) <= 4
            and all(type(name) is str and (re.fullmatch(r'population_[0-9a-f]{64}_(contexts|exceptions)', name)
                or name in ('science_resource_source_owners', 'science_resource_source_records'))
                for name in extra_tables))
        _require({row[0] for row in writer.connection.execute("SELECT name FROM sqlite_master WHERE type='table'")}
            == {'nodes', 'edges', 'roots', *extra_tables})
        writer.nodes._close()
        writer.roots._close()
        writer.connection.commit()
        self.writer, self.connection, self.extra_tables = writer, writer.connection, extra_tables
        self.manifest = EvidenceManifest(tuple(writer.nodes.descriptors), tuple(writer.roots.descriptors),
            writer.node_count, writer.root_count)
        self.version = self._version()
        self._frontier_prepared = False
        self._codec_temporary = self._codec_connection = self._codec_version = None
        self.closed = False

    def _version(self):
        return (self.connection.total_changes,
            self.connection.execute('PRAGMA schema_version').fetchone()[0],
            self.connection.execute('PRAGMA data_version').fetchone()[0])

    def require(self):
        source = self.writer
        _require(not self.closed and source.connection is self.connection and not source.failed and source.manifest is None
            and source._pending_reads == source._building == 0 and not source._fragmented
            and (source.node_count, source.root_count) == (self.manifest.node_count, self.manifest.root_count)
            and tuple(source.nodes.descriptors) == self.manifest.node_shards
            and tuple(source.roots.descriptors) == self.manifest.root_shards
            and self._version() == self.version, 'live evidence prefix changed')
        if self._codec_connection is not None:
            _require(self._index_version(self._codec_connection) == self._codec_version
                and self._codec_connection.execute('PRAGMA query_only').fetchone()[0] == 1,
                'live codec index changed')

    @staticmethod
    def _index_version(connection):
        return (connection.total_changes, connection.execute('PRAGMA schema_version').fetchone()[0],
            connection.execute('PRAGMA data_version').fetchone()[0])

    def _prepare_codec(self):
        temporary, connection = _connect(self.writer.temporary.name)
        try:
            for table, width in (('nodes', 10), ('edges', 2), ('roots', 4)):
                rows = self.connection.execute('SELECT * FROM ' + table)
                try:
                    connection.executemany('INSERT INTO ' + table + ' VALUES ('
                        + ','.join('?' for _ in range(width)) + ')', rows)
                finally:
                    rows.close()
            if self.writer._skip_verification:
                connection.execute('CREATE TABLE suffix_unreferenced(digest TEXT PRIMARY KEY) WITHOUT ROWID')
            else:
                _require_complete_dag_roots(connection, self.manifest.node_count, retain_unreferenced=True)
            connection.commit()
            connection.execute('PRAGMA query_only=ON')
        except BaseException:
            connection.close()
            temporary.cleanup()
            raise
        self._codec_temporary, self._codec_connection = temporary, connection
        self._codec_version = self._index_version(connection)
        self._frontier_prepared = True

    def close(self):
        if not self.closed:
            self.closed = True
            if self._codec_connection is not None:
                self._codec_connection.close()
                self._codec_temporary.cleanup()
                self._codec_temporary = self._codec_connection = self._codec_version = None

    def copy_index(self, reader):
        self.require()
        source, expected, actual = self.writer, self.manifest, reader.manifest
        nodes, roots = len(expected.node_shards), len(expected.root_shards)
        if (actual.node_shards[:nodes] != expected.node_shards
                or actual.root_shards[:roots] != expected.root_shards):
            # A compatible older checkpoint can have a different shard layout.
            # It receives the ordinary complete decoder, not a partial match.
            return 0, 0, 0
        if reader._skip_verification:
            reader._live_prefix_source = source
            _borrow_prefix_index(source.connection, reader.connection)
            return nodes, roots, expected.root_count
        error = None
        with source.read_view(_producer_owned=source._skip_verification) as view:
            _require(view.manifest == expected)
            try:
                # Build one private codec-only index from the unchanged live
                # producer. Each reader receives an independent SQLite copy.
                # Resource projections never enter this index or its readers.
                if self._codec_connection is None:
                    self._prepare_codec()
                self._codec_connection.backup(reader.connection)
            except BaseException as failure:
                error = failure
        self.require()
        if error is not None:
            raise error
        return nodes, roots, expected.root_count


class EvidenceReader:
    """Build an evidence index, with explicit omitted checks for live production."""
    def __init__(self, manifest, directory, *, scratch_directory=None, _produced_prefix=None):
        _require(type(manifest) is EvidenceManifest
            and (_produced_prefix is None or type(_produced_prefix) is _ProducedEvidencePrefix))
        self._produced_prefix = _produced_prefix
        self._live_prefix_source = None
        self.manifest = manifest
        self.manifest_sha256 = manifest.sha256
        self.directory = Path(directory)
        self._skip_verification = ((os.getpid(), str(self.directory.absolute()), manifest.sha256)
            in _LIVE_ARTIFACTS or _produced_prefix is not None and _produced_prefix.writer._skip_verification)
        self.verification_performed = not self._skip_verification
        self._shared_directory = None
        self._shared_metadata_signature = None
        if manifest.shared_prefix is not None:
            from genosapiens.wgs_evidence_shared import shared_directory
            self._shared_directory, self._shared_metadata_signature = shared_directory(manifest, self.directory)
        self.scratch_directory = scratch_directory
        self.connection = None
        self.signatures = []
        self.verified = False
        self._active_derivative = None
        self._read_scope_depth = 0
        self._root_coverage = False
        self._source_writer = None
        self._node_handles = OrderedDict()
        self._node_paths = {}
        self.node_handle_hits = self.node_file_opens = 0

    def _close_node_handles(self):
        while self._node_handles:
            _, handle = self._node_handles.popitem(last=False)
            handle.close()
        self._node_paths.clear()

    def _lookup(self, digest):
        return self.connection.execute("SELECT kind,size,count,depth,first_key,last_key FROM nodes WHERE digest=?", (digest,)).fetchone()

    def _shard_path(self, shard):
        path = self.directory / shard.relative_path
        if self._shared_directory is not None:
            prefix = self.manifest.shared_prefix
            count = prefix.node_shards if shard.kind == 'evidence_nodes' else prefix.root_shards
            number = int(shard.relative_path.removesuffix('.jsonl').rsplit('-', 1)[1])
            if number <= count:
                path = self._shared_directory / shard.relative_path
        return path

    def _open(self, shard):
        path = self._shard_path(shard)
        _require(not any(p.is_symlink() for p in (path, *path.parents)))
        fd = os.open(path, os.O_RDONLY | os.O_NOFOLLOW | os.O_NONBLOCK)
        handle = os.fdopen(fd, "rb")
        info = os.fstat(handle.fileno())
        if not stat.S_ISREG(info.st_mode) or info.st_size != shard.size_bytes:
            handle.close()
            raise EvidenceStoreError()
        return handle

    def _rows(self, shard, number):
        with self._open(shard) as handle:
            signature = _signature(os.fstat(handle.fileno()))
            count, digest, first, last = 0, hashlib.sha256(), None, None
            while True:
                offset = handle.tell()
                line = handle.readline(MAX_NODE_BYTES + 256)
                if not line:
                    break
                _require(len(line) <= MAX_NODE_BYTES + 128 and line.endswith(b"\n"))
                row = _loads(line)
                if not self._skip_verification:
                    _require(_bytes(row) == line)
                count += 1
                _require(count <= shard.row_count)
                if not self._skip_verification:
                    digest.update(line)
                key = row["sha256"] if shard.kind == "evidence_nodes" else row["name"]
                first = key if first is None else first
                last = key
                yield row, offset, len(line)
            _require(count == shard.row_count and (self._skip_verification or digest.hexdigest() == shard.sha256)
                and (first, last) == (shard.first_key, shard.last_key)
                and signature == _signature(os.fstat(handle.fileno())))
        self.signatures.append((shard, signature))

    def _check_prefix_shard(self, shard):
        """Authenticate every copied byte without repeating its proven decode."""
        with self._open(shard) as handle:
            signature = _signature(os.fstat(handle.fileno()))
            if self._skip_verification:
                self.signatures.append((shard, signature))
                return
            digest = hashlib.sha256()
            for block in iter(lambda: handle.read(1024 * 1024), b''):
                digest.update(block)
            _require(digest.hexdigest() == shard.sha256
                and signature == _signature(os.fstat(handle.fileno())))
        self.signatures.append((shard, signature))

    def __enter__(self):
        _require(self.connection is None and not self.verified)
        self.temporary, self.connection = _connect(self.scratch_directory)
        try:
            prefix_nodes = prefix_roots = ordinal = 0
            if self._produced_prefix is not None:
                prefix_nodes, prefix_roots, ordinal = self._produced_prefix.copy_index(self)
            # A copied live index includes its unresolved coverage frontier.
            # Old nodes cannot reference new nodes in this append-only DAG.
            # Suffix edges and roots must cover every remaining prefix node
            # and every new node before the complete reader is accepted.
            for number, shard in enumerate(self.manifest.node_shards):
                if number < prefix_nodes:
                    self._check_prefix_shard(shard)
                    continue
                for row, offset, length in self._rows(shard, number):
                    if not self._skip_verification:
                        _require(set(row) == {"sha256", "node"} and _digest(row["sha256"])
                            and hashlib.sha256(_bytes(row["node"])).hexdigest() == row["sha256"])
                    metadata, children = _metadata(row["node"], self._lookup)
                    self.connection.execute("INSERT INTO nodes VALUES (?,?,?,?,?,?,?,?,?,?)",
                        (row["sha256"], *metadata, number, offset, length))
                    self.connection.executemany("INSERT OR IGNORE INTO edges VALUES (?,?)",
                        ((row["sha256"], child) for child in children))
                    if prefix_nodes and not self._skip_verification:
                        self.connection.execute('INSERT INTO suffix_unreferenced VALUES (?)', (row['sha256'],))
                        self.connection.executemany('DELETE FROM suffix_unreferenced WHERE digest=?',
                            ((child,) for child in children))
            for number, shard in enumerate(self.manifest.root_shards):
                if number < prefix_roots:
                    self._check_prefix_shard(shard)
                    continue
                for row, _, _ in self._rows(shard, number):
                    ordinal += 1
                    _require(set(row) == {"ordinal", "name", "kind", "sha256"}
                        and type(row["ordinal"]) is int and row["ordinal"] == ordinal
                        and self._lookup(row["sha256"]) is not None)
                    EvidenceReference(row["name"], row["kind"], row["sha256"], self.manifest_sha256)
                    self.connection.execute("INSERT INTO roots VALUES (?,?,?,?)",
                        (ordinal, row["name"], row["kind"], row["sha256"]))
                    if prefix_nodes and not self._skip_verification:
                        self.connection.execute('DELETE FROM suffix_unreferenced WHERE digest=?', (row['sha256'],))
            if prefix_nodes:
                if not self._skip_verification:
                    _require(self.connection.execute('SELECT 1 FROM suffix_unreferenced LIMIT 1').fetchone() is None,
                        'unreferenced normalized evidence')
                self.connection.execute('DROP TABLE IF EXISTS suffix_unreferenced')
            elif not self._skip_verification:
                _require_complete_dag_roots(self.connection, self.manifest.node_count)
            self.connection.commit()
            self.verified = True
            self._produced_prefix = None
            return self
        except (sqlite3.Error, ValueError, TypeError, KeyError, AttributeError, OverflowError, RecursionError):
            self.__exit__()
            raise EvidenceStoreError() from None
        except BaseException:
            self.__exit__()
            raise

    def _check_files(self):
        _require(self.verified and self.connection is not None)
        if self._active_derivative is not None:
            _require(not self._active_derivative.failed and self._active_derivative.connection is not None)
            return
        if self._read_scope_depth:
            return
        self._check_file_signatures()

    def _check_file_signatures(self, *, content=False):
        _require(self.verified and self.connection is not None)
        if self._skip_verification:
            return
        # This permission exists only on a live writer view. Saved readers
        # cannot reconstruct it from metadata. Fresh node reads still verify
        # their own digests, while the owner omits repeated full shard scans.
        if self._source_writer is not None and getattr(self, '_producer_owned_read', False):
            content = False
        if self._shared_directory is not None:
            from genosapiens.wgs_evidence_shared import check_metadata
            check_metadata(self._shared_directory, self.manifest.shared_prefix,
                self._shared_metadata_signature, content=content)
        if self._source_writer is not None:
            _require(not self._source_writer.failed and self.connection is self._source_writer.connection,
                     "normalized producing owner changed")
        for shard, signature in self.signatures:
            with self._open(shard) as handle:
                _require(_signature(os.fstat(handle.fileno())) == signature, "normalized evidence changed after authentication")
                if content:
                    # Some filesystems expose identical timestamps for rapid
                    # equal-length writes. A derivative must authenticate all
                    # bytes at its boundaries, including roots it never read.
                    digest = hashlib.sha256()
                    for block in iter(lambda: handle.read(1024 * 1024), b""):
                        digest.update(block)
                    _require(digest.hexdigest() == shard.sha256
                             and _signature(os.fstat(handle.fileno())) == signature,
                             "normalized evidence changed after authentication")

    @contextmanager
    def complete_root_scope(self):
        """Require a schema-specific replay to consume every declared root.

        This is an operation-local disk relation. Merely appearing in a
        manifest or asking for a reference is insufficient: the typed reader
        must authenticate the referenced value during replay or explicitly
        cover it through an exact completed-resource owner. Such reuse retains
        the producer's validation; it does not claim fresh corpus replay.
        """
        _require(self.verified and self.connection is not None and not self._root_coverage
                 and self._source_writer is None and self._read_scope_depth > 0)
        self.connection.execute("CREATE TABLE normalized_root_coverage (name TEXT PRIMARY KEY)")
        # Compatible completed science is covered by its producer, separately
        # from roots decoded by this consumer. Neither relation is persisted.
        self.connection.execute("CREATE TABLE normalized_completed_root_coverage (name TEXT PRIMARY KEY)")
        self._root_coverage = True
        try:
            yield self
            if not self._skip_verification:
                missing = self.connection.execute("SELECT 1 FROM roots r LEFT JOIN normalized_root_coverage c "
                    "ON c.name=r.name LEFT JOIN normalized_completed_root_coverage p ON p.name=r.name "
                    "WHERE r.ordinal<=? AND c.name IS NULL AND p.name IS NULL LIMIT 1",(self.manifest.root_count,)).fetchone()
                _require(missing is None,"scientific evidence contains an unconsumed root")
        finally:
            self._root_coverage = False
            self.connection.execute("DROP TABLE normalized_root_coverage")
            self.connection.execute("DROP TABLE normalized_completed_root_coverage")

    @contextmanager
    def read_scope(self):
        """Authenticate a complete read-only replay at both boundaries.

        Each accessed node still checks its file signature and content digest.
        Nested component replay shares the enclosing full-byte scans. Objects
        assembled in the scope must not escape as accepted results until the
        outer boundary succeeds; failure invalidates this reader.
        """
        _require(self.verified and self.connection is not None
                 and (self._source_writer is None or self._active_derivative is self._source_writer))
        # Component replay may nest inside a derivative's authenticated input
        # scope. That owner checks full bytes at both boundaries and cannot
        # publish while any read remains pending.
        outer = self._read_scope_depth == 0 and self._active_derivative is None
        if outer:
            self._check_file_signatures(content=True)
            self._context_public_shared = None
        self._read_scope_depth += 1
        try:
            yield self
            if outer:
                self._check_file_signatures(content=True)
        except BaseException:
            self.verified = False
            raise
        finally:
            self._read_scope_depth -= 1
            if outer:
                self._close_node_handles()
                self.connection.execute('DROP TABLE IF EXISTS verified_shared_public')
                self._context_public_shared = None

    @contextmanager
    def read_into(self, writer):
        """Authenticate a bulk derivative before its writer may publish.

        Full shard hashes are checked at both boundaries; every individual node
        still checks its source signature and fresh bytes against the digest.
        This removes roots-times-shards scans without admitting a derivative
        before the final check, including corruption of an unused/root shard.
        """
        _require(type(writer) is EvidenceWriter and writer.connection is not None
                 and writer.manifest is None and not writer.failed and self._active_derivative is None)
        self._check_file_signatures(content=True)
        writer._pending_reads += 1
        self._active_derivative = writer
        try:
            yield self
            self._check_file_signatures(content=True)
        except BaseException:
            writer.failed = True
            self.verified = False
            raise
        finally:
            self._active_derivative = None
            writer._pending_reads -= 1
            self._close_node_handles()
            self.connection.execute('DROP TABLE IF EXISTS verified_shared_public')
            self._context_public_shared = None

    def reference(self, name, *, kind):
        self._check_files()
        row = self.connection.execute("SELECT kind,digest FROM roots WHERE name=? AND ordinal<=?", (name,self.manifest.root_count)).fetchone()
        _require(row is not None and row[0] == kind, "missing or foreign typed evidence root")
        return EvidenceReference(name, kind, row[1], self.manifest_sha256)

    def _require_reference(self, reference):
        _require(self.verified and self.connection is not None)
        _require(type(reference) is EvidenceReference and reference.manifest_sha256 == self.manifest_sha256)
        row = self.connection.execute("SELECT kind,digest FROM roots WHERE name=? AND ordinal<=?",
                                      (reference.name,self.manifest.root_count)).fetchone()
        _require(row == (reference.kind, reference.sha256), "foreign typed evidence root")
        if self._root_coverage:
            self.connection.execute("INSERT OR IGNORE INTO normalized_root_coverage VALUES (?)",(reference.name,))

    def _node(self, digest):
        return self._read_node(digest)[0]

    def _read_node(self, digest):
        row = self.connection.execute("SELECT shard,offset,length FROM nodes WHERE digest=?", (digest,)).fetchone()
        _require(row is not None)
        number, offset, length = row
        shard, signature = self.signatures[number]
        scoped = self._read_scope_depth > 0 or self._active_derivative is not None
        if scoped and READER_OPEN_HANDLES > 0:
            handle = self._node_handles.get(shard.relative_path)
            if handle is None:
                while len(self._node_handles) >= READER_OPEN_HANDLES:
                    previous_path, previous = self._node_handles.popitem(last=False)
                    previous.close()
                    self._node_paths.pop(previous_path, None)
                handle = self._open(shard)
                self.node_file_opens += 1
                self._node_handles[shard.relative_path] = handle
                self._node_paths[shard.relative_path] = (
                    (self.directory, self._shared_directory), os.fspath(self._shard_path(shard)))
            else:
                self.node_handle_hits += 1
                self._node_handles.move_to_end(shard.relative_path)
            # Bind the current path and descriptor before and after the read.
            # pread bypasses Python's buffered file bytes, so even a mutation
            # with an unchanged signature still reaches the digest check.
            directory, path = self._node_paths[shard.relative_path]
            _require(directory == (self.directory, self._shared_directory), "normalized reader directory changed")
            if not self._skip_verification:
                _require(_signature(os.stat(path, follow_symlinks=False)) == signature
                    and _signature(os.fstat(handle.fileno())) == signature)
            line = os.pread(handle.fileno(), length, offset)
            if not self._skip_verification:
                _require(_signature(os.fstat(handle.fileno())) == signature
                    and _signature(os.stat(path, follow_symlinks=False)) == signature)
        else:
            with self._open(shard) as handle:
                self.node_file_opens += 1
                _require(_signature(os.fstat(handle.fileno())) == signature)
                handle.seek(offset)
                line = handle.read(length)
                _require(_signature(os.fstat(handle.fileno())) == signature)
        # The authenticated index was built from canonical row envelopes.
        # Check the same digest over the fresh canonical node bytes directly,
        # without parsing and then re-encoding the whole wrapper/node. The
        # expected digest binds canonical bytes; changes still fail here even
        # when the file's signature is unchanged. Decode strictly on EVERY read.
        suffix = b',"sha256":"' + digest.encode('ascii') + b'"}\n'
        _require(len(line) == length and line.startswith(b'{"node":') and line.endswith(suffix))
        content = line[8:-len(suffix)]
        if not self._skip_verification:
            _require(hashlib.sha256(content + b'\n').hexdigest() == digest)
        return _loads(content), content

    def _parts(self, digest, *, interior=False):
        node, wire = self._read_node(digest)
        kind, items = node["kind"], node["items"]
        if kind == "scalar":
            content = wire[len(b'{"items":'):-len(b',"kind":"scalar"}')]
            yield content[1:-1] if interior else content
            return
        del wire
        start, end = (b"{", b"}") if kind in ("object", "objects") else (
            (b'"', b'"') if kind == "strings" else (b"[", b"]"))
        if not interior:
            yield start
        for ordinal, item in enumerate(items):
            if ordinal and kind != "strings":
                yield b","
            if kind == "object":
                key, child = item
                yield _bytes(key)[:-1] + b":"
            else:
                child = item
            yield from self._parts(child, interior=kind in ("arrays", "objects", "strings"))
        if not interior:
            yield end

    def iter_json(self, reference, *, max_bytes=MAX_LOGICAL_BYTES):
        """Canonical logical JSON chunks. Capacity is checked before any output."""
        self._check_files()
        self._require_reference(reference)
        _require(_integer(max_bytes))
        size = self._lookup(reference.sha256)[1]
        _require(size + 1 <= max_bytes, "logical evidence exceeds reconstruction capacity", capacity=True)
        emitted = 0
        for part in self._parts(reference.sha256):
            emitted += len(part)
            _require(emitted <= size)
            yield part
        _require(emitted == size)
        self._check_files()
        yield b"\n"

    def reconstruct(self, reference, *, max_bytes=MAX_RECONSTRUCT_BYTES):
        _require(_integer(max_bytes, MAX_RECONSTRUCT_BYTES), "whole evidence reconstruction exceeds resident capacity", capacity=True)
        return _loads(b"".join(self.iter_json(reference, max_bytes=max_bytes)))

    def array_item(self, reference, index, *, max_item_bytes=MAX_RECONSTRUCT_BYTES):
        """Seek one occurrence through authenticated bounded-fanout metadata."""
        self._check_files()
        self._require_reference(reference)
        metadata = self._lookup(reference.sha256)
        _require(metadata[0] == "array" and type(index) is int and 0 <= index < metadata[2])
        _require(_integer(max_item_bytes, MAX_RECONSTRUCT_BYTES))
        digest = reference.sha256
        while True:
            node = self._node(digest)
            if node["kind"] == "array":
                digest = node["items"][index]
                break
            _require(node["kind"] == "arrays")
            for child in node["items"]:
                count = self._lookup(child)[2]
                if index < count:
                    digest = child
                    break
                index -= count
            else:
                raise EvidenceStoreError("normalized array occurrence is missing")
        size = self._lookup(digest)[1]
        _require(size <= max_item_bytes, "logical evidence item exceeds resident capacity", capacity=True)
        content = b"".join(self._parts(digest))
        _require(len(content) == size)
        self._check_files()
        return _loads(content)

    def _entries(self, digest):
        node = self._node(digest)
        _require(node["kind"] in ("array", "arrays"), "evidence root is not an array")
        for child in node["items"]:
            if node["kind"] == "arrays":
                yield from self._entries(child)
            else:
                yield child

    def iter_array(self, reference, *, max_item_bytes=MAX_RECONSTRUCT_BYTES):
        """Expand one logical array member at a time, including repeated entries."""
        self._check_files()
        self._require_reference(reference)
        yield from self._iter_array_values(reference.sha256, max_item_bytes=max_item_bytes)
        self._check_files()

    def _iter_array_values(self, digest, *, max_item_bytes):
        """Private traversal of a subtree reached from an authenticated root."""
        _require(_integer(max_item_bytes, MAX_RECONSTRUCT_BYTES))
        for digest in self._entries(digest):
            size = self._lookup(digest)[1]
            _require(size <= max_item_bytes, "logical evidence item exceeds resident capacity", capacity=True)
            content = b"".join(self._parts(digest))
            _require(len(content) == size)
            yield _loads(content)

    def __exit__(self, *_):
        self._close_node_handles()
        self._produced_prefix = None
        self.verified = False
        if self.connection is not None:
            if self._source_writer is None:
                self.connection.close()
                self.temporary.cleanup()
            self.connection = None
        self.signatures.clear()
