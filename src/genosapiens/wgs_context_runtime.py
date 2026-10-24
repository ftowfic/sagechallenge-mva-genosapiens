"""Strict consumer profile for the single disposable context SQLite writer.

The profile never configures readers, workers or durable publication stores.
NORMAL scratch is not resumable; successful closure restores FULL durability.
"""
from dataclasses import asdict, dataclass
import hashlib
import json

SCHEMA = 'genosapiens.wgs-context-writer-profile'
MIB = 1024**2


class ContextWriterProfileError(ValueError):
    pass


@dataclass(frozen=True)
class ContextWriterProfile:
    schema: str = SCHEMA
    version: int = 1
    cache_bytes: int = 64*MIB
    transaction_bytes: int = 16*MIB
    transaction_changes: int = 4096
    autocheckpoint_bytes: int = 64*MIB
    truncate_bytes: int = 128*MIB
    max_wal_bytes: int = 2048*MIB
    sparse_batch_commit: bool = True
    audit_index: str = 'always'
    synchronous: str = 'FULL'

    def __post_init__(self):
        integers = ('version', 'cache_bytes', 'transaction_bytes', 'transaction_changes',
                    'autocheckpoint_bytes', 'truncate_bytes', 'max_wal_bytes')
        if (self.schema != SCHEMA or any(type(getattr(self, k)) is not int for k in integers)
                or self.version != 1 or type(self.sparse_batch_commit) is not bool
                or self.max_wal_bytes != 2048*MIB
                or self.audit_index not in ('always', 'when_audited')
                or self.synchronous not in ('FULL', 'NORMAL')):
            raise ContextWriterProfileError('invalid context writer profile')
        small = (16*MIB, 4096, 64*MIB, 128*MIB, True)
        large = (64*MIB, 65536, 512*MIB, 1024*MIB, False)
        actual = (self.transaction_bytes, self.transaction_changes,
                  self.autocheckpoint_bytes, self.truncate_bytes, self.sparse_batch_commit)
        if (self.cache_bytes not in (64*MIB, 512*MIB) or actual not in (small, large)
                or actual == large and self.cache_bytes != 512*MIB
                or self.synchronous == 'NORMAL' and
                (actual != large or self.audit_index != 'when_audited')):
            raise ContextWriterProfileError('unapproved context writer settings')

    def as_bytes(self):
        return json.dumps(asdict(self), sort_keys=True, separators=(',', ':'),
                          allow_nan=False).encode()

    @property
    def sha256(self):
        return hashlib.sha256(self.as_bytes()).hexdigest()


def load_context_writer_profile(content, expected_sha256):
    def unique(pairs):
        result = {}
        for key, value in pairs:
            if key in result:
                raise ContextWriterProfileError('duplicate profile field')
            result[key] = value
        return result
    try:
        if type(content) is not bytes or len(content) > 4096:
            raise ContextWriterProfileError('invalid profile bytes')
        if type(expected_sha256) is not str or hashlib.sha256(content).hexdigest() != expected_sha256:
            raise ContextWriterProfileError('profile digest mismatch')
        data = json.loads(content, object_pairs_hook=unique)
        if type(data) is not dict or set(data) != set(ContextWriterProfile.__dataclass_fields__):
            raise ContextWriterProfileError('unknown or missing profile fields')
        profile = ContextWriterProfile(**data)
        if profile.as_bytes() != content:
            raise ContextWriterProfileError('noncanonical profile bytes')
        return profile
    except (TypeError, ValueError, UnicodeError) as error:
        raise ContextWriterProfileError('invalid context writer profile') from error
