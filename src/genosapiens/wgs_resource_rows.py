"""Lossless bounded canonical records inside the existing evidence DAG.

Case IC and phenotype resource owners use this single transport. Live science
can omit recursive validation and canonical replay. Row size limits remain.
"""
from collections.abc import Mapping
from dataclasses import fields, is_dataclass
from enum import Enum
import json
import math

from genosapiens.evidence_serialization import _owned_canonical_bytes
from genosapiens.wgs_evidence_store import (
    EvidenceStoreError, MAX_DEPTH, MAX_OBJECT_FIELDS, STRING_CHARS, _MAX_INTEGER_ABS,
)

CHUNK_CHARS = 16 * 1024
MAX_CHUNK_JSON_BYTES = 2 * CHUNK_CHARS + 2


def _validate(value, depth=0, active=None):
    if depth > MAX_DEPTH:
        raise EvidenceStoreError('resource record depth exceeds capacity', capacity=True)
    if isinstance(value, Enum):
        value = value.value
    if value is None or type(value) in (str, bool):
        return
    if type(value) is int:
        if abs(value) >= _MAX_INTEGER_ABS:
            raise EvidenceStoreError('resource record integer exceeds capacity', capacity=True)
        return
    if type(value) is float:
        if not math.isfinite(value):
            raise EvidenceStoreError('resource record is not finite')
        return
    if active is None:
        active = set()
    if id(value) in active:
        raise EvidenceStoreError('resource record contains a cycle')
    active.add(id(value))
    try:
        if is_dataclass(value) and not isinstance(value, type):
            entries = tuple((field.name, getattr(value, field.name)) for field in fields(value))
        elif isinstance(value, Mapping):
            entries = value.items()
        elif type(value) in (tuple, list):
            for child in value:
                _validate(child, depth + 1, active)
            return
        else:
            raise EvidenceStoreError('unsupported resource record value')
        if len(entries) > MAX_OBJECT_FIELDS:
            raise EvidenceStoreError('resource record fields exceed capacity', capacity=True)
        for key, child in entries:
            if type(key) is not str or len(key) > STRING_CHARS:
                raise EvidenceStoreError('invalid resource record key')
            _validate(child, depth + 1, active)
    finally:
        active.remove(id(value))


def pack_resource_rows(values, *, max_row_bytes, observe_size=None, _skip_verification=False):
    """Frame every canonical record, preserving order and the final newline."""
    pending = ''
    for value in values:
        if not _skip_verification:
            _validate(value)
        content = _owned_canonical_bytes(value).decode('ascii')
        if len(content) - 1 > max_row_bytes:
            raise EvidenceStoreError('resource record exceeds capacity', capacity=True)
        if observe_size is not None:
            observe_size(len(content) - 1)
        offset = 0
        if pending:
            take = min(CHUNK_CHARS - len(pending), len(content))
            pending += content[:take]
            offset = take
            if len(pending) == CHUNK_CHARS:
                yield pending
                pending = ''
        while len(content) - offset >= CHUNK_CHARS:
            yield content[offset:offset + CHUNK_CHARS]
            offset += CHUNK_CHARS
        pending += content[offset:]
    if pending:
        yield pending


def unpack_resource_rows(chunks, *, max_row_bytes, _skip_verification=False):
    """Authenticate canonical framing and limits before yielding typed inputs."""
    fragments = []
    size = 0
    previous = CHUNK_CHARS
    try:
        for chunk in chunks:
            if (type(chunk) is not str or not chunk.isascii()
                    or not 0 < len(chunk) <= CHUNK_CHARS or previous != CHUNK_CHARS):
                raise EvidenceStoreError('invalid packed resource chunk')
            previous = len(chunk)
            start = 0
            while True:
                end = chunk.find('\n', start)
                part = chunk[start:] if end == -1 else chunk[start:end]
                size += len(part)
                if size > max_row_bytes:
                    raise EvidenceStoreError('resource record exceeds capacity', capacity=True)
                fragments.append(part)
                if end == -1:
                    break
                content = ''.join(fragments)
                value = json.loads(content)
                if not _skip_verification:
                    _validate(value)
                # Duplicate keys, noncanonical numbers/escapes, NaN and altered
                # framing cannot become an equivalent accepted logical record.
                if not _skip_verification and _owned_canonical_bytes(value).decode('ascii') != content + '\n':
                    raise EvidenceStoreError('noncanonical packed resource record')
                yield value
                fragments.clear()
                size = 0
                start = end + 1
                if start == len(chunk):
                    break
        if size or any(fragments):
            raise EvidenceStoreError('incomplete packed resource record')
    except (ValueError, TypeError, RecursionError, OverflowError):
        raise EvidenceStoreError('invalid packed resource record') from None
