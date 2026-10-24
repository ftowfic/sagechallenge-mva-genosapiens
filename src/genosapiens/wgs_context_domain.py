"""Explicit ordered context coverage for a selected or completed domain.

Legacy collectors remain consecutive by default. An explicit domain preserves
original ordinals and still rejects missing, repeated, or foreign contexts.
"""
from dataclasses import dataclass
import hashlib
import struct

from genosapiens.wgs_evidence_store import EvidenceStoreError

MAX_CONTEXTS = 1_000_000


@dataclass(frozen=True, slots=True)
class ContextDomain:
    packed: bytes

    def __post_init__(self):
        if type(self.packed) is not bytes or len(self.packed) % 4 or len(self.packed) > 4 * MAX_CONTEXTS:
            raise EvidenceStoreError('invalid scientific context domain')
        previous = 0
        for (ordinal,) in struct.iter_unpack('<I', self.packed):
            if ordinal <= previous:
                raise EvidenceStoreError('scientific context domain must increase')
            previous = ordinal

    @classmethod
    def from_ordinals(cls, ordinals):
        packed = bytearray()
        previous = 0
        for ordinal in ordinals:
            if type(ordinal) is not int or not previous < ordinal < 2**32 or len(packed) >= MAX_CONTEXTS * 4:
                raise EvidenceStoreError('invalid scientific context domain')
            packed.extend(struct.pack('<I', ordinal))
            previous = ordinal
        return cls(bytes(packed))

    def __len__(self):
        return len(self.packed) // 4

    def __iter__(self):
        return (row[0] for row in struct.iter_unpack('<I', self.packed))

    def at(self, completed):
        if type(completed) is not int or not 0 <= completed < len(self):
            raise EvidenceStoreError('scientific context is outside the declared domain')
        return struct.unpack_from('<I', self.packed, completed * 4)[0]

    @property
    def sha256(self):
        return hashlib.sha256(b'genosapiens.context-domain/v1\0' + self.packed).hexdigest()


def next_context_ordinal(owner):
    domain = getattr(owner, 'context_domain', None)
    if domain is None:
        return owner.context_count + 1
    if type(domain) is not ContextDomain:
        raise EvidenceStoreError('foreign scientific context domain')
    return domain.at(owner.context_count)


def require_complete_domain(owner):
    domain = getattr(owner, 'context_domain', None)
    if domain is not None and (type(domain) is not ContextDomain or owner.context_count != len(domain)):
        raise EvidenceStoreError('scientific context domain is incomplete')
