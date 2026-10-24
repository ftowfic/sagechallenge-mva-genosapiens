"""Durable resource prefixes with separate context-specific evidence files.

Only a live producer can publish a shared bundle. Saved readers authenticate
its manifest and every referenced shard without accepting a persisted index.
"""
import hashlib
import os
from pathlib import Path
import stat
import tempfile

from genosapiens.wgs_evidence_store import (
    EvidenceReader, SharedEvidencePrefix, MAX_MANIFEST_BYTES,
    _ProducedEvidencePrefix, _require, _signature, manifest_from_bytes,
)


def _safe(path):
    path = Path(path).absolute()
    _require(not any(p.is_symlink() for p in (path, *path.parents)))
    return path


def _sync(path):
    descriptor = os.open(path, os.O_RDONLY | os.O_DIRECTORY | os.O_NOFOLLOW)
    try:
        os.fsync(descriptor)
    finally:
        os.close(descriptor)


def _metadata(directory):
    path = _safe(directory / 'evidence.json')
    descriptor = os.open(path, os.O_RDONLY | os.O_NOFOLLOW | os.O_NONBLOCK)
    with os.fdopen(descriptor, 'rb') as stream:
        info = os.fstat(stream.fileno())
        _require(stat.S_ISREG(info.st_mode) and 0 < info.st_size <= MAX_MANIFEST_BYTES)
        content = stream.read(MAX_MANIFEST_BYTES + 1)
        _require(len(content) == info.st_size and _signature(info) == _signature(os.fstat(stream.fileno())))
        return content, _signature(info)


def check_metadata(directory, prefix, signature, *, content):
    path = _safe(directory / 'evidence.json')
    _require(_signature(path.stat()) == signature, 'shared resource manifest changed')
    if content:
        raw, actual = _metadata(directory)
        _require(actual == signature and hashlib.sha256(raw).hexdigest() == prefix.sha256)


def shared_directory(manifest, context_directory):
    """Resolve one sibling bundle, never an arbitrary path from saved input."""
    prefix = manifest.shared_prefix
    _require(type(prefix) is SharedEvidencePrefix)
    directory = _safe(Path(context_directory).parent / ('.resources-' + prefix.sha256))
    raw, signature = _metadata(directory)
    shared = manifest_from_bytes(raw, expected_sha256=prefix.sha256)
    _require(shared.shared_prefix is None and shared.version == 1
        and len(shared.node_shards) == prefix.node_shards
        and len(shared.root_shards) == prefix.root_shards
        and manifest.node_shards[:prefix.node_shards] == shared.node_shards
        and manifest.root_shards[:prefix.root_shards] == shared.root_shards)
    return directory, signature


class CheckpointResourcePrefix:
    """One checked bundle under the context checkpoint store's case lock."""
    def __init__(self, root, live):
        _require(type(live) is _ProducedEvidencePrefix)
        live.require()
        _require(live.manifest.shared_prefix is None)
        self.root = _safe(root)
        self.live = live
        self.descriptor = SharedEvidencePrefix(live.manifest.sha256,
            len(live.manifest.node_shards), len(live.manifest.root_shards))
        self.directory = self.root / ('.resources-' + self.descriptor.sha256)
        if not self.directory.exists():
            pending = Path(tempfile.mkdtemp(prefix='.pending-resources-', dir=self.root))
            with live.writer.read_view() as reader:
                _require(reader.manifest == live.manifest)
                for shard in (*reader.manifest.node_shards, *reader.manifest.root_shards):
                    digest, size = hashlib.sha256(), 0
                    with reader._open(shard) as source:
                        signature = _signature(os.fstat(source.fileno()))
                        descriptor = os.open(pending / shard.relative_path,
                            os.O_WRONLY | os.O_CREAT | os.O_EXCL | os.O_NOFOLLOW, 0o600)
                        with os.fdopen(descriptor, 'wb') as target:
                            for block in iter(lambda: source.read(1024 * 1024), b''):
                                size += len(block)
                                _require(size <= shard.size_bytes)
                                digest.update(block)
                                target.write(block)
                            _require(size == shard.size_bytes and digest.hexdigest() == shard.sha256
                                and signature == _signature(os.fstat(source.fileno())))
                            target.flush()
                            os.fsync(target.fileno())
                    os.chmod(pending / shard.relative_path, 0o400)
            live.require()
            descriptor = os.open(pending / 'evidence.json',
                os.O_WRONLY | os.O_CREAT | os.O_EXCL | os.O_NOFOLLOW, 0o600)
            with os.fdopen(descriptor, 'wb') as target:
                target.write(live.manifest.as_bytes())
                target.flush()
                os.fsync(target.fileno())
            os.chmod(pending / 'evidence.json', 0o400)
            _sync(pending)
            _require(not self.directory.exists())
            os.rename(pending, self.directory)
            _sync(self.root)
        self._check_bundle()

    def _check_bundle(self):
        raw, _ = _metadata(_safe(self.directory))
        manifest = manifest_from_bytes(raw, expected_sha256=self.descriptor.sha256)
        _require(manifest == self.live.manifest)
        reader = EvidenceReader(manifest, self.directory)
        reader._skip_verification = self.live.writer._skip_verification
        for shard in (*manifest.node_shards, *manifest.root_shards):
            reader._check_prefix_shard(shard)
        self._metadata_signature = _metadata(self.directory)[1]
        self._signatures = tuple(reader.signatures)

    def require(self, source, context_directory):
        self.live.require()
        _require(source is self.live.writer and _safe(context_directory).parent == self.root)
        if source._skip_verification:
            return
        check_metadata(self.directory, self.descriptor, self._metadata_signature, content=False)
        for shard, signature in self._signatures:
            path = _safe(self.directory / shard.relative_path)
            _require(_signature(path.stat()) == signature, 'live shared resource shard changed')
