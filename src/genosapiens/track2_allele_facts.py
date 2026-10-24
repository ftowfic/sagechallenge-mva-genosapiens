"""Bounded invocation-local reuse of identical allele assessment facts."""
from collections import OrderedDict
import json

from genosapiens.track2_published_genotype import canonical


class AlleleFactCache:
    """Keep immutable JSON values, including completed empty source matches.

    The assessment owner keys entries by source bytes, context, and policy.
    Genotype checks and therapeutic decisions do not enter this cache.
    """
    def __init__(self, *, maximum_entries=256, maximum_bytes=16 * 1024 * 1024):
        if (type(maximum_entries) is not int or maximum_entries <= 0
                or type(maximum_bytes) is not int or maximum_bytes <= 0):
            raise ValueError('Allele fact cache limits must be positive integers.')
        self.maximum_entries = maximum_entries
        self.maximum_bytes = maximum_bytes
        self.entries = OrderedDict()
        self.bytes = self.hits = self.misses = 0

    def resolve(self, key, producer):
        if key in self.entries:
            self.hits += 1
            self.entries.move_to_end(key)
            return json.loads(self.entries[key])
        self.misses += 1
        value = producer()  # A failed lookup never creates a cache entry.
        content = canonical(value)
        size = len(key.encode('ascii')) + len(content)
        if size <= self.maximum_bytes:
            self.entries[key] = content
            self.bytes += size
            while len(self.entries) > self.maximum_entries or self.bytes > self.maximum_bytes:
                old_key, old_content = self.entries.popitem(last=False)
                self.bytes -= len(old_key.encode('ascii')) + len(old_content)
        return json.loads(content)
