"""Bounded sharing of identical raw JSON objects within one member decode.

This is parser-local allocation reuse, not scientific validation or authority.
Every wire row is parsed anew and its ordinary typed codec runs afterwards.
"""
from collections import OrderedDict
import json
import math
import sys


MAX_INTERN_BYTES = 8 * 1024**2
MAX_INTERN_ENTRIES = 4096
MAX_KEY_ITEMS = 4096


class _Object(dict):
    __slots__ = ('ordinal', 'resident_size')


class _TooLarge(Exception):
    pass


class _ObjectPool:
    def __init__(self):
        self.cache = OrderedDict()
        self.bytes = self.ordinal = self.hits = 0

    def object(self, pairs):
        value = _Object(pairs)
        if len(value) != len(pairs):
            raise ValueError('duplicate member JSON object key')
        self.ordinal += 1
        value.ordinal = self.ordinal
        value.resident_size = MAX_INTERN_BYTES + 1
        if not MAX_INTERN_BYTES or not MAX_INTERN_ENTRIES:
            return value
        size = sys.getsizeof(value)
        remaining = MAX_KEY_ITEMS

        def token(child):
            nonlocal size, remaining
            remaining -= 1
            kind = type(child)
            size += child.resident_size if kind is _Object else sys.getsizeof(child)
            if (remaining < 0 or 1024 + 4 * size
                    + 128 * (MAX_KEY_ITEMS - remaining) > MAX_INTERN_BYTES):
                raise _TooLarge
            if kind is _Object:
                # Monotonic ordinals are never reused, including after eviction.
                # Equal ordinals designate this parse's exact same raw object.
                return 6, child.ordinal
            if child is None:
                return (0,)
            if kind is bool:
                return 1, child
            if kind is int:
                return 2, child
            if kind is float:
                if not math.isfinite(child):
                    raise _TooLarge
                return 3, child.hex()
            if kind is str:
                return 4, child
            if kind is list:
                return 5, tuple(token(item) for item in child)
            raise TypeError('unsupported member JSON value')

        try:
            keys = []
            for key, child in pairs:
                size += 128 + sys.getsizeof(key)
                keys.append((key, token(child)))
            key = tuple(keys)
        except _TooLarge:
            return value
        value.resident_size = size
        cached = self.cache.get(key)
        if cached is not None:
            self.cache.move_to_end(key)
            self.hits += 1
            return cached[0]
        # Conservative accounting includes the raw tree, token tuples and LRU
        # overhead, charging shared descendants repeatedly rather than omitting
        # their residency. Oversized objects simply keep ordinary allocation.
        # Token tuples and float.hex() strings need their own allowance; raw
        # scalar size alone undercounts long arrays of floating-point values.
        cost = 1024 + 4 * size + 128 * (MAX_KEY_ITEMS - remaining)
        while self.cache and (self.bytes + cost > MAX_INTERN_BYTES
                or len(self.cache) >= MAX_INTERN_ENTRIES):
            self.bytes -= self.cache.popitem(last=False)[1][1]
        if cost <= MAX_INTERN_BYTES:
            self.cache[key] = value, cost
            self.bytes += cost
        return value


def load_member_json(content):
    pool = _ObjectPool()
    try:
        return json.loads(content, object_pairs_hook=pool.object)
    finally:
        # Keys and cache ownership do not survive parsing, even on failure.
        pool.cache.clear()
        pool.bytes = 0
