"""Operation-owned, lossless context payload factoring (private scratch v1).

Shared fields are canonical JSON fragments, not independently trusted objects.
The worker authenticates the complete typed result before encoding; the parent
checks fragment hashes and row identities before accepting the live result.
Scratch reads reconstruct ordinary canonical bytes for the strict consumers.
"""
from collections import OrderedDict
from dataclasses import fields
import hashlib
import json
import zlib

from genosapiens.evidence_serialization import _plain
from genosapiens.wgs_failures import WgsFailure
from genosapiens.wgs_ledger import canonical

MAX_BYTES = 64 * 1024**2
CACHE_BYTES = 8 * 1024**2
SHARED = frozenset(('support_ledger', 'evidence_refs', 'source_refs', 'context',
                    'phenotype_proof_ids', 'records'))


class ContextPayloadError(RuntimeError):
    pass


class ContextPayloadCapacityError(ContextPayloadError, WgsFailure):
    def __init__(self, operation):
        WgsFailure.__init__(self, 'bounded_products', 'response_capacity', operation=operation,
                            exception_kind='project_failure')


def unpack(blob, size, digest):
    if type(blob) is not bytes or type(size) is not int or not 0 < size <= MAX_BYTES:
        raise ContextPayloadError('context fragment capacity is invalid')
    decoder = zlib.decompressobj()
    try:
        content = decoder.decompress(blob, size + 1)
    except zlib.error:
        raise ContextPayloadError('context fragment compression is invalid') from None
    if (len(content) != size or not decoder.eof or decoder.unused_data
            or decoder.unconsumed_tail or hashlib.sha256(content).hexdigest() != digest):
        raise ContextPayloadError('context fragment integrity is invalid')
    return content


class Encoder:
    """One bounded job's fragments; identity reuse retains its actual owner."""
    def __init__(self):
        self.fragments = {}
        self.cache = OrderedDict()
        self.cache_bytes = 0
        self.bytes = 0

    def fragment(self, value):
        item = self.cache.get(id(value))
        if item is not None and item[0] is value:
            self.cache.move_to_end(id(value))
            return item[1], item[2]
        content = canonical(_plain(value))
        if len(content) > MAX_BYTES:
            raise ContextPayloadCapacityError('context_fragment_encode')
        digest = hashlib.sha256(content).hexdigest()
        if digest not in self.fragments:
            blob = zlib.compress(content, 1)
            self.bytes += len(content)
            if self.bytes > MAX_BYTES:
                raise ContextPayloadCapacityError('context_payload_encode')
            self.fragments[digest] = (len(content), blob)
        if len(content) <= CACHE_BYTES:
            while self.cache and self.cache_bytes + len(content) > CACHE_BYTES:
                _, old = self.cache.popitem(last=False)
                self.cache_bytes -= old[2]
            self.cache[id(value)] = (value, digest, len(content))
            self.cache_bytes += len(content)
        return digest, len(content)

    def record(self, value):
        data, shared, size = {}, {}, 0
        for field in fields(value):
            child = getattr(value, field.name)
            if field.name in SHARED:
                digest, length = self.fragment(child)
                shared[field.name] = digest
                size += length
            else:
                data[field.name] = _plain(child)
        content = canonical(dict(version=1, data=data, shared=shared))
        self.bytes += len(content)
        if self.bytes > MAX_BYTES:
            raise ContextPayloadCapacityError('context_payload_encode')
        if size + len(content) > MAX_BYTES:
            raise ContextPayloadCapacityError('context_record_encode')
        return content

    def rows(self):
        return tuple((digest, size, blob) for digest, (size, blob) in sorted(self.fragments.items()))


def header(content):
    if type(content) is not bytes or not 0 < len(content) <= MAX_BYTES:
        raise ContextPayloadError('context record capacity is invalid')
    data = json.loads(content)
    if (type(data) is not dict or set(data) != {'version', 'data', 'shared'}
            or type(data['version']) is not int or data['version'] != 1
            or type(data['data']) is not dict or type(data['shared']) is not dict
            or set(data['data']) & set(data['shared'])
            or not set(data['shared']) <= SHARED or canonical(data) != content):
        raise ContextPayloadError('context record structure is invalid')
    return data


def expand(content, get):
    value = header(content)
    data = value['data']
    total = len(content)
    for name, digest in value['shared'].items():
        fragment = get(digest)
        total += len(fragment)
        if total > MAX_BYTES:
            raise ContextPayloadError('expanded context record exceeds capacity')
        data[name] = json.loads(fragment)
    result = canonical(data)
    if len(result) > MAX_BYTES:
        raise ContextPayloadError('expanded context record exceeds capacity')
    return result


class Payloads:
    """Small keys in indexes; each compressed payload is stored only once."""
    def __init__(self, connection):
        self.connection = connection
        self.cache = OrderedDict()
        self.cache_bytes = 0
        self.write_cache = OrderedDict()
        self.write_cache_bytes = 0

    def create(self):
        self.connection.execute('CREATE TABLE context_payloads '
            '(sha256 TEXT PRIMARY KEY NOT NULL,size INTEGER NOT NULL,content BLOB NOT NULL)')

    def put_fragments(self, rows):
        # EncodedContextPart owns hash validation. Compare compressed bytes on
        # collisions, including independently encoded duplicate fragments.
        for digest, size, blob in rows:
            previous = self.write_cache.get(digest)
            if previous is not None:
                self.write_cache.move_to_end(digest)
            else:
                previous = self.connection.execute('SELECT size,content FROM context_payloads WHERE sha256=?',
                                                   (digest,)).fetchone()
            if previous is None:
                self.connection.execute('INSERT INTO context_payloads VALUES (?,?,?)', (digest, size, blob))
            elif previous != (size, blob):
                raise ContextPayloadError('context payload identity conflicts')
            if digest not in self.write_cache and len(blob) <= CACHE_BYTES:
                while self.write_cache and self.write_cache_bytes + len(blob) > CACHE_BYTES:
                    _, old = self.write_cache.popitem(last=False)
                    self.write_cache_bytes -= len(old[1])
                self.write_cache[digest] = (size, blob)
                self.write_cache_bytes += len(blob)

    def put(self, content):
        digest = hashlib.sha256(content).hexdigest()
        self.put_fragments(((digest, len(content), zlib.compress(content, 1)),))
        return digest

    def get(self, digest):
        # Cache belongs to this private SQLite owner, never external files or
        # a serialized "validated" marker. Its lifetime ends with the owner.
        if digest in self.cache:
            self.cache.move_to_end(digest)
            return self.cache[digest]
        row = self.connection.execute('SELECT size,content FROM context_payloads WHERE sha256=?', (digest,)).fetchone()
        if row is None:
            raise ContextPayloadError('context payload is missing')
        content = unpack(row[1], row[0], digest)
        if len(content) <= CACHE_BYTES:
            while self.cache and self.cache_bytes + len(content) > CACHE_BYTES:
                _, old = self.cache.popitem(last=False)
                self.cache_bytes -= len(old)
            self.cache[digest] = content
            self.cache_bytes += len(content)
        return content

    def expanded(self, digest):
        return expand(self.get(digest), self.get)
