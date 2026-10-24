"""Generic bounded metadata, shard and private-index I/O.

The obsolete expanded scientific/order/selection codecs are removed. Current
science is exclusively ScientificArtifact plus normalized scientific products.
"""
from collections.abc import Iterator
import ctypes
from fractions import Fraction
import gc
import hashlib
import json
import os
from pathlib import Path
import sqlite3
import stat
import tempfile
from genosapiens.bounded_wgs_config import BoundedWgsExecutionConfig
from genosapiens.evidence_model import _sha, _typed
from genosapiens.evidence_serialization import _canonical
from genosapiens.wgs_shards import ScientificShard, ShardKind, _ShardWriter, WgsScientificError, _release_shard_cache
_ERRORS = (ValueError, TypeError, RuntimeError, AttributeError, KeyError, IndexError,
    OSError, sqlite3.Error, UnicodeError, OverflowError, RecursionError)
_LARGE_ROW_BYTES = 1024**2
try:
    _MALLOC_TRIM = ctypes.CDLL(None).malloc_trim
except (AttributeError, OSError):
    _MALLOC_TRIM = None
else:
    _MALLOC_TRIM.argtypes = [ctypes.c_size_t]
    _MALLOC_TRIM.restype = ctypes.c_int


def _release_large_row_memory() -> None:
    """Return freed parsing buffers between bounded, authenticated large rows."""

    gc.collect()
    if _MALLOC_TRIM is not None:
        _MALLOC_TRIM(0)



def _bytes(value: object) -> bytes:
    return (_canonical(value) + "\n").encode("ascii")



def _hash(value: object) -> str:
    return hashlib.sha256(_bytes(value)).hexdigest()



def _load_manifest(
    content: bytes,
    *,
    expected_sha256: str,
    result_type,
    result_key: str,
    max_bytes: int = 8 * 1024 * 1024,
):
    try:
        _sha(expected_sha256)
        if (
            type(max_bytes) is not int or max_bytes < 1
            or type(content) is not bytes
            or len(content) > max_bytes
            or hashlib.sha256(content).hexdigest() != expected_sha256
        ):
            raise ValueError
        data = json.loads(content)
        expected = {"classification", "llm_policy", result_key}
        if (
            type(data) is not dict
            or set(data) != expected
            or data["classification"] != "controlled"
            or data["llm_policy"] != "no-llm"
        ):
            raise ValueError
        value = _typed(data[result_key], result_type, decode=True)
        if value.as_bytes() != content:
            raise ValueError
        return value
    except _ERRORS:
        raise WgsScientificError("scientific manifest binding is invalid") from None



def _write_manifest(root: Path, name: str, content: bytes, *, max_bytes: int = 8 * 1024 * 1024) -> None:
    if type(max_bytes) is not int or max_bytes < 1:
        raise ValueError("manifest capacity is invalid")
    if len(content) > max_bytes:
        raise WgsScientificError("scientific manifest exceeds capacity", category="response_capacity")
    path = root / name
    descriptor = os.open(
        path, os.O_WRONLY | os.O_CREAT | os.O_EXCL, 0o600,
    )
    with os.fdopen(descriptor, "wb") as handle:
        handle.write(content)
        handle.flush()
        os.fsync(handle.fileno())



def _scientific_shard_lines(root, shard, *, max_row_bytes=None):
    """Bounded stable reads, including authenticated publication hard links."""
    if max_row_bytes is not None and (type(max_row_bytes) is not int or max_row_bytes<1):
        raise WgsScientificError("scientific row capacity is invalid")
    path = root / shard.relative_path
    if (not root.is_absolute() or path.is_symlink()
            or any(parent.is_symlink() for parent in path.parents)
            or not path.resolve(strict=True).is_relative_to(root.resolve(strict=True))):
        raise WgsScientificError("scientific shard path is invalid")
    descriptor = os.open(path, os.O_RDONLY | os.O_NOFOLLOW | os.O_NONBLOCK)
    with os.fdopen(descriptor, "rb") as handle:
        try:
            before = os.fstat(handle.fileno())
            if not stat.S_ISREG(before.st_mode) or before.st_size != shard.size_bytes:
                raise WgsScientificError("scientific shard size mismatch")
            size = 0
            line_limit=shard.size_bytes if max_row_bytes is None else min(shard.size_bytes,max_row_bytes)
            while content := handle.readline(line_limit + 1):
                size += len(content)
                if size > shard.size_bytes or len(content)>line_limit or not content.endswith(b"\n"):
                    raise WgsScientificError("scientific shard row exceeds capacity")
                yield content
                del content
            after = os.fstat(handle.fileno())
            if (before.st_dev, before.st_ino, before.st_size, before.st_mtime_ns,
                before.st_ctime_ns) != (after.st_dev, after.st_ino, after.st_size,
                                        after.st_mtime_ns, after.st_ctime_ns):
                raise WgsScientificError("scientific shard changed during read")
        finally:
            _release_shard_cache(handle.fileno())


def _read_rows(
    root: Path,
    shards: tuple[ScientificShard, ...],
    kind: ShardKind,
    row_type,
    *,
    expected_stream_sha256: str | None = None,
    row_key,
    selected_row: int | None = None,
    max_row_bytes: int | None = None,
) -> Iterator[tuple[ScientificShard, int, object]]:
    # A location lookup still authenticates the entire bounded shard. Only its
    # endpoints and requested row need typed decoding; full consumers use the
    # default path and validate every row before completing their stream.
    if selected_row is not None and (
        type(selected_row) is not int or len(shards) != 1
        or not 1 <= selected_row <= shards[0].row_count
    ):
        raise WgsScientificError("scientific row location is invalid")
    if max_row_bytes is not None and (type(max_row_bytes) is not int or max_row_bytes < 1):
        raise WgsScientificError("scientific row capacity is invalid")


    def decode(content: bytes):
        value = _typed(json.loads(content), row_type, decode=True)
        if _bytes(value) != content:
            raise ValueError
        key = row_key(value)
        return value, key

    stream = hashlib.sha256()
    for shard in shards:
        if shard.kind != kind:
            raise WgsScientificError("scientific shard kind mismatch")
        path = root / shard.relative_path
        digest = hashlib.sha256()
        size = 0
        count = 0
        first = last = None
        last_offset = 0
        with tempfile.TemporaryFile(mode="w+b") as spool:
            try:
                for count, content in enumerate(_scientific_shard_lines(root,shard,max_row_bytes=max_row_bytes), 1):
                    size += len(content)
                    digest.update(content)
                    stream.update(content)
                    last_offset = spool.tell()
                    if spool.write(content) != len(content):
                        raise OSError
                    del content
            except _ERRORS:
                raise WgsScientificError("scientific shard verification failed") from None
            if (
                count != shard.row_count
                or size != shard.size_bytes
                or digest.hexdigest() != shard.sha256
            ):
                raise WgsScientificError("scientific shard closure mismatch")
            # Authenticate all saved bytes before exposing a row. Decode the
            # endpoints for descriptor closure and reuse those typed rows in
            # the consumer pass. This reader holds at most two rows from one
            # bounded shard; descriptor checks still precede the first yield.
            try:
                spool.seek(0)
                endpoint = spool.readline(shard.size_bytes + 1)
                first_value, first = decode(endpoint)
                del endpoint
                last_value = None
                if count == 1:
                    last = first
                else:
                    spool.seek(last_offset)
                    endpoint = spool.readline(shard.size_bytes + 1)
                    last_value, last = decode(endpoint)
                    del endpoint
                if first != shard.first_key or last != shard.last_key:
                    raise WgsScientificError("scientific shard closure mismatch")
                if selected_row is not None:
                    if selected_row != 1:
                        first_value = None
                    if selected_row != count:
                        last_value = None
                spool.seek(0)
                for row_number in range(1, count + 1):
                    content = spool.readline(shard.size_bytes + 1)
                    if selected_row is not None and row_number != selected_row:
                        del content
                        continue
                    large_row = len(content) >= _LARGE_ROW_BYTES
                    if row_number == 1:
                        value, first_value = first_value, None
                    elif row_number == count:
                        value, last_value = last_value, None
                    else:
                        value, _key = decode(content)
                    yield shard, row_number, value
                    del value, content
                    if large_row:
                        _release_large_row_memory()
                if spool.read(1) != b"":
                    raise WgsScientificError("scientific shard spool closure mismatch")
            except _ERRORS:
                raise WgsScientificError("scientific shard replay failed") from None
    if (
        expected_stream_sha256 is not None
        and stream.hexdigest() != expected_stream_sha256
    ):
        raise WgsScientificError("scientific shard stream mismatch")



def _open_index(
    path: Path, execution_config: BoundedWgsExecutionConfig, *,
    existing: bool = False,
) -> sqlite3.Connection:
    if (not path.is_absolute() or path.is_symlink()
            or (not path.is_file() if existing else path.exists())):
        raise WgsScientificError("scientific index path is invalid")
    connection = (
        sqlite3.connect(path.as_uri() + "?mode=rw", uri=True)
        if existing else sqlite3.connect(path)
    )
    cache_kib = max(1, execution_config.sqlite.cache_bytes // 1024)
    connection.execute("PRAGMA foreign_keys=ON")
    connection.execute(f"PRAGMA cache_size=-{cache_kib}")
    connection.execute("PRAGMA temp_store=FILE")
    connection.execute("PRAGMA mmap_size=0")
    connection.execute("PRAGMA trusted_schema=OFF")

    def exact(left: str, right: str) -> int:
        first, second = Fraction(left), Fraction(right)
        return (first > second) - (first < second)

    connection.create_collation("EXACT_FRACTION", exact)
    return connection
