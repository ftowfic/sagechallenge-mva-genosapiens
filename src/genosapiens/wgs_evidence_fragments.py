"""Private live-worker DAG fragments; never saved evidence or recovery authority.

The coordinator must keep its derivative unpublished until every worker closes
its authenticated input scope. Only fixed input roots may be read by a worker;
completed output is consumed in order and retired before its next job starts.
There is deliberately no public artifact reader or manifest in this protocol.
"""
from dataclasses import asdict
import hashlib
import os
from pathlib import Path
import stat

from genosapiens.evidence_serialization import _loads
from genosapiens.wgs_evidence_store import (
    EvidenceReader, EvidenceReference, EvidenceWriter, MAX_DEPTH, MAX_MANIFEST_BYTES,
    MAX_NODE_BYTES, MAX_NODES, MAX_ROOTS, _bytes, _digest, _integer, _require, _signature)
from genosapiens.wgs_ledger import MAX_SHARD_BYTES, MAX_SHARD_ROWS
from genosapiens.wgs_shards import ScientificShard, _ShardWriter


class _FragmentOutput:
    """One disposable producer, with at most one unconsumed job fragment."""
    def __init__(self, writer, reader, *, max_bytes=None):
        _require(type(writer) is EvidenceWriter and type(reader) is EvidenceReader
            and not writer.failed and not writer._fragmented and writer._building == 0
            and writer._pending_reads == 1 and reader.verified
            and reader._source_writer is writer and reader._active_derivative is writer
            and reader.connection is writer.connection
            and (max_bytes is None or type(max_bytes) is int and max_bytes > 0))
        self.writer, self.reader = writer, reader
        self.max_bytes, self.bytes = max_bytes, 0
        self.ordinal = 0
        self.active = self.pending = None
        self.completed = False
        # Irreversible for this producing owner, including after input closure.
        writer._fragmented = self

    def begin(self, ordinal):
        w = self.writer
        _require(not self.completed and not w.failed and w.connection is not None
            and w._pending_reads == 1 and self.reader._active_derivative is w
            and self.active is None and self.pending is None
            and type(ordinal) is int and ordinal > self.ordinal)
        # The fixed reader owns its original manifest/signatures and base-node
        # locations. New locations are disposable and may never be read back.
        for attr, kind in (("nodes", "evidence_nodes"), ("roots", "evidence_roots")):
            previous = getattr(w, attr)
            _require(previous._handle is None)
            setattr(w, attr, _ShardWriter(w.directory, kind, max_rows=MAX_SHARD_ROWS,
                max_bytes=MAX_SHARD_BYTES, filename_prefix=f"fragment-{ordinal:08d}-{kind}"))
        self.ordinal = ordinal
        self.bytes = 0
        self.active = (w.node_count, w.root_count)

    def charge(self, size):
        """Reserve bytes before either node or root bytes reach a shard."""
        try:
            _require(self.active is not None and type(size) is int and size > 0)
            _require(self.max_bytes is None or self.bytes + size <= self.max_bytes,
                'worker fragment output exceeds reserved bytes', capacity=True)
            self.bytes += size
        except BaseException:
            self.writer.failed = True
            raise

    def seal(self):
        w = self.writer
        _require(self.active is not None and self.pending is None and not w.failed
            and w.connection is not None and w._building == 0 and w._pending_reads == 1)
        node_start, root_start = self.active
        value = dict(ordinal=self.ordinal, node_start=node_start, node_end=w.node_count,
            root_start=root_start, root_end=w.root_count,
            nodes=[asdict(s) for s in w.nodes.close()], roots=[asdict(s) for s in w.roots.close()])
        _require(len(_bytes(value)) <= MAX_MANIFEST_BYTES and value['root_end'] > root_start)
        _require(sum(s['size_bytes'] for key in ('nodes', 'roots') for s in value[key]) == self.bytes)
        self.pending, self.active = value, None
        return value

    def retire(self, ordinal):
        w = self.writer
        _require(not w.failed and self.pending is not None and type(ordinal) is int
            and ordinal == self.pending['ordinal'] and w._pending_reads == 1)
        for section in ('nodes', 'roots'):
            for shard in self.pending[section]:
                path = w.directory / shard['relative_path']
                _require(path.is_file() and not path.is_symlink())
                path.unlink()
        # Counts and on-disk uniqueness/dependency metadata remain cumulative;
        # the next job can reference hashes consumed earlier by this parent.
        for shards in (w.nodes, w.roots):
            shards.descriptors.clear()
            shards.descriptor_bytes = 0
        self.pending = None

    def complete(self):
        w = self.writer
        _require(not self.completed and self.active is None and self.pending is None
            and not w.failed and w.connection is not None and w._pending_reads == 0
            and self.reader._active_derivative is None)
        self.completed = True
        w.close()


def _fragment_rows(directory, shard):
    path = Path(directory) / shard.relative_path
    _require(not any(p.is_symlink() for p in (path, *path.parents)))
    fd = os.open(path, os.O_RDONLY | os.O_NOFOLLOW | os.O_NONBLOCK)
    with os.fdopen(fd, 'rb') as handle:
        info = os.fstat(handle.fileno())
        _require(stat.S_ISREG(info.st_mode) and info.st_size == shard.size_bytes)
        signature = _signature(info)
        digest, count, first, last = hashlib.sha256(), 0, None, None
        while True:
            raw = handle.readline(MAX_NODE_BYTES + 256)
            if not raw:
                break
            _require(len(raw) <= MAX_NODE_BYTES + 128 and raw.endswith(b'\n'))
            row = _loads(raw)
            _require(_bytes(row) == raw and type(row) is dict)
            key = row['sha256'] if shard.kind == 'evidence_nodes' else row['name']
            first = key if first is None else first
            last = key
            digest.update(raw)
            count += 1
            _require(count <= shard.row_count)
            yield row
        _require(count == shard.row_count and digest.hexdigest() == shard.sha256
            and (first, last) == (shard.first_key, shard.last_key)
            and signature == _signature(os.fstat(handle.fileno())))


class _FragmentImport:
    """Incremental consumer for one live worker; carries only counts, not trust."""
    def __init__(self, writer, *, node_count, root_count):
        _require(type(writer) is EvidenceWriter and not writer.failed
            and writer.connection is not None and not writer._fragmented
            and _integer(node_count, MAX_NODES) and _integer(root_count, MAX_ROOTS))
        self.writer = writer
        self.counts = (node_count, root_count)
        self.ordinal = 0

    def append(self, directory, value, *, ordinal):
        w = self.writer
        created = False
        try:
            # Existing read_into owns the coordinator publication barrier.
            _require(not w.failed and w.connection is not None and w.manifest is None
                and w._pending_reads > 0 and w._building == 0)
            _require(type(value) is dict and set(value) == {
                'ordinal', 'node_start', 'node_end', 'root_start', 'root_end', 'nodes', 'roots'}
                and len(_bytes(value)) <= MAX_MANIFEST_BYTES
                and type(ordinal) is int and ordinal > self.ordinal
                and type(value['ordinal']) is int and value['ordinal'] == ordinal
                and (value['node_start'], value['root_start']) == self.counts)
            for key in ('node_start', 'node_end', 'root_start', 'root_end'):
                _require(_integer(value[key], MAX_NODES if key.startswith('node') else MAX_ROOTS))
            _require(value['node_end'] >= value['node_start'] and value['root_end'] > value['root_start'])
            for key, kind in (('nodes', 'evidence_nodes'), ('roots', 'evidence_roots')):
                _require(type(value[key]) is list)
                shards = [ScientificShard(**s) for s in value[key]]
                for number, shard in enumerate(shards, 1):
                    _require(shard.kind == kind and shard.relative_path ==
                        f'fragment-{ordinal:08d}-{kind}-{number:08d}.jsonl')
                start, end = ('node_start', 'node_end') if key == 'nodes' else ('root_start', 'root_end')
                _require(sum(s.row_count for s in shards) == value[end] - value[start])
            w.connection.execute('CREATE TEMP TABLE fragment_nodes (digest TEXT PRIMARY KEY, reached INTEGER NOT NULL) WITHOUT ROWID')
            created = True
            for descriptor in value['nodes']:
                for row in _fragment_rows(directory, ScientificShard(**descriptor)):
                    _require(set(row) == {'sha256', 'node'} and _digest(row['sha256'])
                        and type(row['node']) is dict and set(row['node']) == {'kind', 'items'}
                        and hashlib.sha256(_bytes(row['node'])).hexdigest() == row['sha256'])
                    _require(w._node(row['node']['kind'], row['node']['items']) == row['sha256'])
                    w.connection.execute('INSERT INTO fragment_nodes VALUES (?,0)', (row['sha256'],))
            root_ordinal = value['root_start']
            for descriptor in value['roots']:
                for row in _fragment_rows(directory, ScientificShard(**descriptor)):
                    root_ordinal += 1
                    _require(set(row) == {'ordinal', 'name', 'kind', 'sha256'}
                        and type(row['ordinal']) is int and row['ordinal'] == root_ordinal
                        and w._lookup(row['sha256']) is not None)
                    EvidenceReference(row['name'], row['kind'], row['sha256'], '0' * 64)
                    w._add(row['name'], row['kind'], lambda: row['sha256'])
                    w.connection.execute('UPDATE fragment_nodes SET reached=1 WHERE digest=?', (row['sha256'],))
            # New nodes must reach this job's roots. Older dependencies were
            # already authenticated; do not traverse the growing base again.
            for _ in range(MAX_DEPTH):
                before = w.connection.total_changes
                w.connection.execute('UPDATE fragment_nodes SET reached=1 WHERE reached=0 AND digest IN '
                    '(SELECT child FROM edges JOIN fragment_nodes ON parent=digest WHERE reached=1)')
                if w.connection.total_changes == before:
                    break
            _require(w.connection.execute('SELECT count(*) FROM fragment_nodes WHERE reached=0').fetchone()[0] == 0,
                'unreferenced worker fragment evidence')
            self.counts, self.ordinal = (value['node_end'], value['root_end']), ordinal
        except BaseException:
            w.failed = True
            raise
        finally:
            if created and w.connection is not None:
                w.connection.execute('DROP TABLE IF EXISTS fragment_nodes')
