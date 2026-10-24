"""Incremental scientific publication. Callers own the destination output lock.

The processing unit is a scalar chunk (32 KiB characters), a mapping's keys,
or one sequence member. No complete serialized document is constructed.
Scientific semantics and small control-record admission do not live here.
"""
from collections.abc import Mapping
from dataclasses import dataclass
import errno
import hashlib
import json
import os
from pathlib import Path
import stat
import sys
import tempfile

CHUNK_CHARS = 32 * 1024
READ_BYTES = 1024 * 1024


class ScientificStorageError(ValueError):
    """Static storage errors contain no artifact content or private paths."""


class ScientificCapacityError(ScientificStorageError):
    """A bounded control unit cannot hold its required metadata."""


@dataclass
class Records:
    """An ordered, single-pass iterable serialized as a JSON array."""
    values: object


def json_chunks(value):
    """Yield canonical ASCII JSON, including the established trailing newline."""
    active = set()

    def encode(item):
        if isinstance(item, str):
            yield b'"'
            for start in range(0, len(item), CHUNK_CHARS):
                yield json.dumps(item[start:start + CHUNK_CHARS], ensure_ascii=True)[1:-1].encode('ascii')
            yield b'"'
        elif item is None or isinstance(item, (bool, int, float)):
            yield json.dumps(item, allow_nan=False, separators=(',', ':')).encode('ascii')
        elif isinstance(item, (Mapping, list, tuple, Records)):
            identity = id(item)
            if identity in active:
                raise ScientificStorageError('Scientific JSON contains a circular reference.')
            active.add(identity)
            try:
                if isinstance(item, Mapping):
                    # Scientific contracts use string keys. Reject implicit key conversion.
                    if any(not isinstance(key, str) for key in item):
                        raise ScientificStorageError('Scientific JSON requires string keys.')
                    yield b'{'
                    for index, key in enumerate(sorted(item)):
                        if index:
                            yield b','
                        yield from encode(key)
                        yield b':'
                        yield from encode(item[key])
                    yield b'}'
                else:
                    yield b'['
                    for index, member in enumerate(item.values if isinstance(item, Records) else item):
                        if index:
                            yield b','
                        yield from encode(member)
                    yield b']'
            finally:
                active.remove(identity)
        else:
            raise ScientificStorageError('Scientific JSON contains an unsupported value.')

    yield from encode(value)
    yield b'\n'


def json_identity(value):
    """Hash canonical bytes without constructing a complete serialized copy."""
    digest = hashlib.sha256()
    size = 0
    for chunk in json_chunks(value):
        digest.update(chunk)
        size += len(chunk)
    return {'size': size, 'sha256': digest.hexdigest()}


def _path(path):
    path = Path(path).absolute()
    if any(p.is_symlink() for p in (path, *path.parents)):
        raise ScientificStorageError('A scientific artifact path contains a symbolic link.')
    if path.exists() and not path.is_file():
        raise ScientificStorageError('A scientific artifact has an invalid file type.')
    return path


def _private_parent(path):
    missing = []
    parent = path.parent
    while not parent.exists():
        missing.append(parent)
        parent = parent.parent
    for directory in reversed(missing):
        directory.mkdir(mode=0o700)
    info = path.parent.stat()
    if info.st_uid != os.getuid() or stat.S_IMODE(info.st_mode) & 0o077:
        raise ScientificStorageError('Scientific output requires a private owned directory.')


def publish_chunks(path, chunks):
    """Durably replace a private artifact and return its complete byte identity.

    A temporary hard link retains prior bytes until directory synchronization.
    Failure restores the previous name where the filesystem permits recovery.
    Cleanup failures never replace the primary exception. Leftover temporary
    files remain incomplete and belong to the caller's locked output inventory.
    """
    path = _path(path)
    _private_parent(path)
    if path.exists():
        info = path.stat()
        if info.st_uid != os.getuid() or stat.S_IMODE(info.st_mode) & 0o077:
            raise ScientificStorageError('An existing scientific artifact is not private.')
    directory = os.open(path.parent, os.O_RDONLY | os.O_DIRECTORY | os.O_NOFOLLOW)
    temporary = previous = None
    replaced = False
    rollback_failed = False
    try:
        descriptor, name = tempfile.mkstemp(prefix='.scientific-', dir=path.parent)
        temporary = Path(name)
        digest = hashlib.sha256()
        size = 0
        with os.fdopen(descriptor, 'wb') as stream:
            os.fchmod(stream.fileno(), 0o600)
            for chunk in chunks(temporary) if callable(chunks) else chunks:
                if not isinstance(chunk, bytes):
                    raise ScientificStorageError('Scientific publication requires byte chunks.')
                if stream.write(chunk) != len(chunk):
                    raise OSError(errno.EIO, 'Incomplete scientific write.')
                digest.update(chunk)
                size += len(chunk)
            stream.flush()
            os.fsync(stream.fileno())
        if path.exists():
            previous = temporary.with_suffix('.previous')
            os.link(path, previous, follow_symlinks=False)
        os.replace(temporary, path)
        replaced = True
        os.fsync(directory)
    except BaseException:
        if replaced:
            try:
                if previous is not None:
                    os.replace(previous, path)
                else:
                    path.unlink()
                os.fsync(directory)
            except OSError:
                # No successful return or completion follows a failed rollback.
                # Keep the old inode available for recovery under the lock.
                rollback_failed = True
        raise
    finally:
        for pending in (temporary, None if rollback_failed else previous):
            if pending is not None:
                try:
                    pending.unlink(missing_ok=True)
                except OSError:
                    pass
        primary = sys.exc_info()[1]
        try:
            os.close(directory)
        except OSError:
            if primary is None:
                raise
    return {'size': size, 'sha256': digest.hexdigest()}


def publish_json(path, value):
    """Publish a complete scientific JSON artifact without a control-file limit."""
    return publish_chunks(path, json_chunks(value))


def file_identity(path):
    """Read complete artifact bytes with a fixed-size buffer."""
    path = _path(path)
    digest = hashlib.sha256()
    size = 0
    descriptor = os.open(path, os.O_RDONLY | os.O_NOFOLLOW | os.O_NONBLOCK)
    with os.fdopen(descriptor, 'rb') as stream:
        before = os.fstat(stream.fileno())
        if not stat.S_ISREG(before.st_mode):
            raise ScientificStorageError('A scientific artifact has an invalid file type.')
        while chunk := stream.read(READ_BYTES):
            size += len(chunk)
            digest.update(chunk)
        after = os.fstat(stream.fileno())
        identity = lambda s: (s.st_dev, s.st_ino, s.st_size, s.st_mtime_ns, s.st_ctime_ns)
        if identity(before) != identity(after) or identity(after) != identity(path.stat()):
            raise ScientificStorageError('A scientific artifact changed during reading.')
    return {'size': size, 'sha256': digest.hexdigest()}
