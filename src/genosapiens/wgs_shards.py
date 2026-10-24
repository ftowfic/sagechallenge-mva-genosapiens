"""Shared descriptor for bounded scientific record shards."""

from __future__ import annotations

from dataclasses import dataclass, fields
import hashlib
import re
from types import MappingProxyType
import os
from pathlib import Path
import stat
from typing import Literal

from genosapiens.evidence_model import _sha
from genosapiens.evidence_serialization import _canonical
from genosapiens.wgs_failures import WgsFailure
from genosapiens.wgs_ledger import MAX_MANIFEST_BYTES, MAX_SHARD_BYTES, MAX_SHARD_ROWS


def read_bounded_metadata(path, *, max_bytes=MAX_MANIFEST_BYTES):
    """Read stable metadata under an explicit limit; the shared default is 8 MiB."""
    if type(max_bytes) is not int or max_bytes < 1:
        raise ValueError("metadata capacity is invalid")
    path = Path(path)
    if any(item.is_symlink() for item in (path, *path.parents)):
        raise WgsFailure("bounded_products", "malformed_evidence", operation="bounded_product_validate", exception_kind="project_failure")
    descriptor = os.open(path, os.O_RDONLY | os.O_NOFOLLOW | os.O_NONBLOCK)
    with os.fdopen(descriptor, "rb") as handle:
        before = os.fstat(handle.fileno())
        if not stat.S_ISREG(before.st_mode) or before.st_size < 1:
            raise WgsFailure("bounded_products", "malformed_evidence", operation="bounded_product_validate", exception_kind="project_failure")
        if before.st_size > max_bytes:
            raise WgsFailure("bounded_products", "response_capacity", operation="bounded_product_validate", exception_kind="project_failure")
        content = handle.read(max_bytes + 1)
        after = os.fstat(handle.fileno())
        if (len(content) != before.st_size or
                (before.st_dev, before.st_ino, before.st_size, before.st_mtime_ns, before.st_ctime_ns)
                != (after.st_dev, after.st_ino, after.st_size, after.st_mtime_ns, after.st_ctime_ns)):
            raise WgsFailure("bounded_products", "malformed_evidence", operation="bounded_product_validate", exception_kind="project_failure")
    return content


ShardKind = Literal[
    "read_evidence", "intake", "validation_audit", "evidence_nodes", "evidence_roots", "associations", "context_subjects", "profile_decisions",
]

class WgsScientificError(WgsFailure):
    """Static bounded scientific-artifact failure."""

    def __init__(self, message, *, category="malformed_evidence"):
        super().__init__(
            "bounded_products", category,
            operation="bounded_product_validate", exception_kind="project_failure",
        )
        RuntimeError.__init__(self, message)



def _bytes(value):
    return (_canonical(value) + "\n").encode("ascii")


def _release_shard_cache(descriptor: int) -> None:
    """Release clean pages of one closed stream without changing saved bytes."""
    try:
        os.posix_fadvise(descriptor, 0, 0, os.POSIX_FADV_DONTNEED)
    except (AttributeError, OSError):
        raise WgsScientificError("scientific shard cache release failed") from None


_CANONICAL_ROW_SEAL = object()


@dataclass(frozen=True, slots=True, repr=False)
class CanonicalRow:
    """Exact bytes owned by one producing operation, never a saved trust flag."""
    row: object
    content: bytes
    seal: object
    field_spans: object = None


def encode_row(row, *, field_spans=False):
    if not field_spans:
        return CanonicalRow(row,_bytes(row),_CANONICAL_ROW_SEAL)
    parts,spans,offset = [b'{'],{},1
    for number,field in enumerate(sorted(fields(row),key=lambda field:field.name)):
        prefix = (b',' if number else b'') + _bytes(field.name)[:-1] + b':'
        content = _bytes(getattr(row,field.name))[:-1]
        parts.extend((prefix,content))
        offset += len(prefix)
        spans[field.name] = offset,offset+len(content)
        offset += len(content)
    parts.append(b'}\n')
    return CanonicalRow(row,b''.join(parts),_CANONICAL_ROW_SEAL,MappingProxyType(spans))


def row_content(row, encoded):
    if (type(encoded) is not CanonicalRow or encoded.seal is not _CANONICAL_ROW_SEAL
            or encoded.row is not row or type(encoded.content) is not bytes):
        raise WgsScientificError("canonical scientific row owner changed")
    return encoded.content

class _ShardWriter:
    def __init__(
        self,
        root: Path,
        kind: ShardKind,
        *,
        max_rows: int,
        max_bytes: int,
        filename_prefix: str | None = None,
        row_encoder=None,
    ) -> None:
        self.root = root
        self.encode_row = encode_row if row_encoder is None else row_encoder
        self.kind = kind
        self.max_rows = max_rows
        self.max_bytes = max_bytes
        self.filename_prefix = kind if filename_prefix is None else filename_prefix
        if re.fullmatch(r"[a-z0-9_-]+", self.filename_prefix) is None:
            raise WgsScientificError("scientific shard prefix is invalid")
        self.descriptors: list[ScientificShard] = []
        self.descriptor_bytes = 0
        self.stream = hashlib.sha256()
        self._handle = None
        self._digest = None
        self._path = None
        self._relative = None
        self._rows = 0
        self._size = 0
        self._first = None
        self._last = None
        self.metrics=dict(shards=0,shard_close_nanoseconds=0,fsync_nanoseconds=0)

    def _open(self) -> None:
        ordinal = len(self.descriptors) + 1
        relative = f"{self.filename_prefix}-{ordinal:08d}.jsonl"
        path = self.root / relative
        if path.exists():
            raise WgsScientificError("scientific shard path already exists")
        descriptor = os.open(
            path, os.O_WRONLY | os.O_CREAT | os.O_EXCL, 0o600,
        )
        self._handle = os.fdopen(descriptor, "wb")
        self._digest = hashlib.sha256()
        self._path = path
        self._relative = relative

    def _close(self) -> None:
        if self._handle is None:
            return
        import time
        started=time.monotonic_ns()
        self._handle.flush()
        sync_started=time.monotonic_ns()
        try:os.fsync(self._handle.fileno())
        finally:self.metrics['fsync_nanoseconds']+=time.monotonic_ns()-sync_started
        try:
            _release_shard_cache(self._handle.fileno())
        finally:
            self._handle.close()
            self._handle = None
        descriptor = ScientificShard(
            self.kind,
            self._relative,
            self._digest.hexdigest(),
            self._size,
            self._rows,
            self._first,
            self._last,
        )
        self.descriptor_bytes += len(_bytes(descriptor))
        if self.descriptor_bytes > 8 * 1024 * 1024:
            raise WgsScientificError("scientific shard descriptor budget exceeded", category="response_capacity")
        self.descriptors.append(descriptor)
        self.metrics['shards']+=1
        self.metrics['shard_close_nanoseconds']+=time.monotonic_ns()-started
        self._handle = None
        self._digest = None
        self._path = None
        self._relative = None
        self._rows = 0
        self._size = 0
        self._first = None
        self._last = None

    def sync(self) -> None:
        """Durably publish current bytes without closing a bounded shard."""
        if self._handle is None:
            return
        import time
        self._handle.flush()
        started=time.monotonic_ns()
        try:os.fsync(self._handle.fileno())
        finally:self.metrics['fsync_nanoseconds']+=time.monotonic_ns()-started

    def write(self, value: object, *, key: str, canonical_bytes=None, capture=None) -> tuple[str, int]:
        encoded = self.encode_row(value) if canonical_bytes is None else canonical_bytes
        content = row_content(value, encoded)
        if len(content) > self.max_bytes:
            raise WgsScientificError("scientific row exceeds shard capacity", category="response_capacity")
        if (
            self._handle is not None
            and (self._rows >= self.max_rows
                 or self._size + len(content) > self.max_bytes)
        ):
            self._close()
        if self._handle is None:
            self._open()
        self._handle.write(content)
        self._digest.update(content)
        self.stream.update(content)
        self._rows += 1
        self._size += len(content)
        self._first = key if self._first is None else self._first
        self._last = key
        if capture is not None:
            capture(value, canonical_bytes=encoded, location=(self._relative, self._rows))
        return self._relative, self._rows

    def close(self) -> tuple[ScientificShard, ...]:
        self._close()
        return tuple(self.descriptors)


@dataclass(frozen=True, slots=True, repr=False)
class ScientificShard:
    kind: ShardKind
    relative_path: str
    sha256: str
    size_bytes: int
    row_count: int
    first_key: str
    last_key: str

    def __post_init__(self):
        path = Path(self.relative_path)
        _sha(self.sha256)
        if (self.kind not in (
                "read_evidence", "intake", "validation_audit", "evidence_nodes", "evidence_roots", "associations", "context_subjects", "profile_decisions",
            ) or self.relative_path != path.as_posix() or path.is_absolute()
                or not path.parts or any(part in ("", ".", "..") for part in path.parts)
                or type(self.size_bytes) is not int or not 1 <= self.size_bytes <= MAX_SHARD_BYTES
                or type(self.row_count) is not int or not 1 <= self.row_count <= MAX_SHARD_ROWS
                or not self.first_key or not self.last_key):
            raise ValueError("scientific shard descriptor is invalid")
