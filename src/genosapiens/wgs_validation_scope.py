"""Small, operation-owned validation reuse for bounded scientific members.

Reaching the entry limit disables additional reuse. It never skips validation.
The scientific owners retain their existing input and output byte limits.
"""
from contextlib import contextmanager

from genosapiens.evidence_model import _VALIDATION_WALK, _validation_session

MAX_VALIDATION_ENTRIES = 8192


class _BoundedEntries(dict):
    def __init__(self, limit):
        super().__init__()
        self.limit = limit

    def __setitem__(self, key, value):
        if key in self or len(self) < self.limit:
            super().__setitem__(key, value)


@contextmanager
def member_validation_scope():
    """Release all retained input and owned-result identities after one operation."""
    with _validation_session(isolated=True):
        active, _ = _VALIDATION_WALK.get()
        cache = _BoundedEntries(MAX_VALIDATION_ENTRIES)
        token = _VALIDATION_WALK.set((active, cache))
        try:
            yield
        finally:
            cache.clear()
            active.clear()
            _VALIDATION_WALK.reset(token)
