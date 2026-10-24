"""File-scoped clean-page release for immutable bounded input readers."""
from __future__ import annotations

import os
from pathlib import Path
import sqlite3
import stat
import threading
import time


class ReadFileCache:
    """Borrow an owner's read-only descriptor; never change or close its file."""

    def __init__(self, descriptor: int):
        self.descriptor = descriptor
        self._signature = self._stat()
        self._last_release = 0.0
        self._lock = threading.Lock()
        self.releases = 0
        self.release_nanoseconds = 0

    def _stat(self):
        value = os.fstat(self.descriptor)
        if not stat.S_ISREG(value.st_mode):
            raise OSError("read cache requires a regular file")
        return (value.st_dev, value.st_ino, value.st_size, value.st_mtime_ns,
                value.st_ctime_ns, value.st_mode, value.st_uid, value.st_nlink)

    def release(self, *, force=False):
        with self._lock:
            now = time.monotonic()
            if not force and now - self._last_release < 1:
                return
            started = time.monotonic_ns()
            # SQLite's separately bounded pager retains hot pages. Kernel file
            # pages need not retain every previously authenticated disk block.
            os.posix_fadvise(self.descriptor, 0, 0, os.POSIX_FADV_DONTNEED)
            self.release_nanoseconds += time.monotonic_ns() - started
            self.releases += 1
            self._last_release = now


class _CacheBoundedConnection(sqlite3.Connection):
    _read_cache = None

    def close(self):
        cache = self._read_cache
        self._read_cache = None
        try:
            if cache is not None:
                cache.release(force=True)
        finally:
            try:
                super().close()
            finally:
                if cache is not None:
                    os.close(cache.descriptor)


def cache_bounded_readonly_connection(path: Path):
    """Keep release authority on the exact file for one read-only connection."""
    descriptor = os.open(path, os.O_RDONLY | os.O_NOFOLLOW)
    connection = None
    try:
        cache = ReadFileCache(descriptor)
        connection = sqlite3.connect(path.as_uri() + "?mode=ro", uri=True,
                                     factory=_CacheBoundedConnection)
        connection._read_cache = cache
        # Recheck the pathname so the SQLite connection and release descriptor
        # cannot silently select different files after a concurrent replacement.
        metadata = path.stat()
        if (metadata.st_dev, metadata.st_ino) != cache._signature[:2]:
            raise OSError("read cache connection file changed")
        cache.release(force=True)
        connection.set_progress_handler(lambda: (cache.release(), 0)[1], 100000)
        return connection
    except BaseException:
        if connection is None:
            os.close(descriptor)
        elif connection._read_cache is None:
            try:
                connection.close()
            finally:
                os.close(descriptor)
        else:
            connection.close()
        raise
