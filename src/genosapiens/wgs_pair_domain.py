"""Complete distinct-event pair traversal with bounded resident state.

This is a domain owner, not a scientific filter or an outcome certificate. The
owning context must close admission/source membership first. No qualification,
frequency, consequence or score can remove a member here. Native LongTR locus
parents remain in member accounting while their projected events can pair.
"""
from __future__ import annotations

from dataclasses import dataclass
from array import array
from bisect import bisect_left, bisect_right
import hashlib
from pathlib import Path
import sqlite3
import tempfile
import sys

from genosapiens.evidence_model import AlleleEvidence, _typed
from genosapiens.evidence_serialization import _canonical
from genosapiens.genotype_properties import _allele_event
from genosapiens import longtr_admission as longtr
from genosapiens.wgs_evidence_store import EvidenceStoreError

MAX_CONTEXT_MEMBERS = 100_000
PAIR_PAGE_ROWS = 1024
MAX_COMPLETE_PAIRS = MAX_CONTEXT_MEMBERS * (MAX_CONTEXT_MEMBERS - 1) // 2
MAX_POSITION_BYTES = 64 * 1024**2
SCHEMA = "genosapiens.wgs-complete-pair-domain"
VERSION = 1
_PAIR_SQL = ("SELECT a.allele_id,b.allele_id FROM members a NOT INDEXED "
    "CROSS JOIN members b NOT INDEXED ON b.ordinal>a.ordinal "
    "WHERE a.parent=0 AND b.parent=0 AND a.event!=b.event ORDER BY a.ordinal,b.ordinal")


@dataclass(frozen=True, slots=True, repr=False)
class PairDomainSummary:
    member_count: int
    pairable_count: int
    distinct_event_count: int
    pair_count: int
    membership_sha256: str
    schema: str = SCHEMA
    version: int = VERSION

    def __post_init__(self):
        counts = (self.member_count, self.pairable_count, self.distinct_event_count, self.pair_count)
        if (any(type(count) is not int or count < 0 for count in counts)
                or not self.distinct_event_count <= self.pairable_count <= self.member_count <= MAX_CONTEXT_MEMBERS
                or bool(self.distinct_event_count) != bool(self.pairable_count)
                or self.pair_count > self.pairable_count * (self.pairable_count - 1) // 2
                or self.schema != SCHEMA or type(self.version) is not int or self.version != VERSION
                or type(self.membership_sha256) is not str or len(self.membership_sha256) != 64
                or any(c not in "0123456789abcdef" for c in self.membership_sha256)):
            raise EvidenceStoreError("invalid complete pair domain summary")


class CompletePairDomain:
    """An exact, operation-local SQLite owner over all context members.

    Build exhausts its input before any traversal becomes available. The input
    order is canonical allele-id order, preserving the historical combinations
    order. Every bounded page may cross any former admission batch boundary.
    """
    def __init__(self, *, scratch_directory=None):
        self.temporary = tempfile.TemporaryDirectory(prefix=".pair-domain-", dir=scratch_directory)
        self.connection = sqlite3.connect(Path(self.temporary.name) / "pairs.sqlite")
        self.connection.execute("PRAGMA cache_size=-16384")
        self.connection.execute("PRAGMA mmap_size=0")
        self.connection.execute("PRAGMA temp_store=FILE")
        self.connection.execute("CREATE TABLE members (ordinal INTEGER PRIMARY KEY, "
            "allele_id TEXT UNIQUE NOT NULL, native_id TEXT NOT NULL, event TEXT NOT NULL, parent INTEGER NOT NULL)")
        self.connection.execute("CREATE INDEX event_members ON members(event,parent)")
        self.summary = None
        self.failed = False

    def build(self, alleles):
        if self.connection is None or self.summary is not None or self.failed:
            raise EvidenceStoreError("pair domain owner is unavailable")
        try:
            digest, previous, count, pairable = hashlib.sha256(), None, 0, 0
            digest.update(b"[")
            for allele in alleles:
                count += 1
                if count > MAX_CONTEXT_MEMBERS:
                    raise EvidenceStoreError("complete context member capacity exceeded", capacity=True)
                allele = _typed(allele, AlleleEvidence)
                if previous is not None and allele.allele_id <= previous:
                    raise EvidenceStoreError("pair domain membership is duplicated or reordered")
                event = _canonical(_allele_event(allele))
                parent = longtr.is_parent(allele)
                membership = (allele.allele_id, allele.native_identity_id, event, parent)
                digest.update((b"," if previous is not None else b"") + _canonical(membership).encode("ascii"))
                self.connection.execute("INSERT INTO members VALUES (?,?,?,?,?)",
                    (count, allele.allele_id, allele.native_identity_id, event, int(parent)))
                pairable += not parent
                previous = allele.allele_id
            digest.update(b"]\n")
            groups, duplicates = 0, 0
            for (size,) in self.connection.execute("SELECT count(*) FROM members WHERE parent=0 GROUP BY event"):
                groups += 1
                duplicates += size * (size - 1) // 2
            total = pairable * (pairable - 1) // 2 - duplicates
            if not 0 <= total <= MAX_COMPLETE_PAIRS:
                raise EvidenceStoreError("complete pair count capacity exceeded", capacity=True)
            self.connection.commit()
            self.connection.execute("PRAGMA query_only=ON")
            self.summary = PairDomainSummary(count, pairable, groups, total, digest.hexdigest())
            return self.summary
        except BaseException:
            self.failed = True
            raise

    def _require(self):
        if self.connection is None or self.summary is None or self.failed:
            raise EvidenceStoreError("complete pair domain is not closed")

    def members(self):
        self._require()
        for ordinal, key, native, event, parent in self.connection.execute(
                "SELECT ordinal,allele_id,native_id,event,parent FROM members ORDER BY ordinal"):
            yield (ordinal, key, native, event, bool(parent))

    def pages(self):
        """All unordered distinct-event pairs; no truncation/limit argument."""
        self._require()
        # The rowid range on b supplies historical lexical combination order.
        # No all-pairs table, materialized Cartesian result, or sort is needed.
        cursor = self.connection.execute(_PAIR_SQL)
        emitted = 0
        while page := cursor.fetchmany(PAIR_PAGE_ROWS):
            emitted += len(page)
            if emitted > self.summary.pair_count:
                raise EvidenceStoreError("complete pair domain count changed")
            yield tuple(page)
        if emitted != self.summary.pair_count:
            raise EvidenceStoreError("complete pair domain count changed")

    def pairs(self):
        for page in self.pages():
            yield from page

    def positions(self):
        self._require()
        return PairDomainPositions(self)

    def close(self):
        if self.connection is not None:
            self.connection.close()
            self.connection = None
            self.temporary.cleanup()

    def __enter__(self):
        return self

    def __exit__(self, *_):
        self.close()


class _PairPositionQueries:
    def _finish_positions(self, pair_count):
        self.prefix = array("Q",[0])
        total = 0
        for ordinal in range(1,self.member_count+1):
            index = self.indices[ordinal]
            if index != 2**32-1:
                group = self.groups[self.events[ordinal]]
                total += len(self.pairable)-index-1 - (len(group)-bisect_right(group,index))
            self.prefix.append(total)
        if total != pair_count:
            raise EvidenceStoreError("pair position count changed")
        self.array_bytes = sum(sys.getsizeof(v) for v in (self.events,self.indices,self.pairable,self.prefix,self.groups,*self.groups))

    def ordinal(self, left, right):
        """One-based complete domain position for two one-based member ordinals."""
        self._require_open()
        if (type(left) is not int or type(right) is not int or not 1 <= left < right < self.member_count+1
                or self.indices[left] == 2**32-1 or self.indices[right] == 2**32-1
                or self.events[left] == self.events[right]):
            raise EvidenceStoreError("locator is outside the complete pair domain")
        lo,hi = self.indices[left],self.indices[right]
        group = self.groups[self.events[left]]
        excluded = bisect_left(group,hi)-bisect_right(group,lo)
        return self.prefix[left-1] + hi-lo-1-excluded + 1

    def locator(self, ordinal):
        self._require_open()
        if type(ordinal) is not int or not 1 <= ordinal <= self.prefix[-1]:
            raise EvidenceStoreError("pair domain position is unavailable")
        left = bisect_right(self.prefix,ordinal-1)
        index = self.indices[left]
        group = self.groups[self.events[left]]
        desired = ordinal-self.prefix[left-1]
        lower,upper = index+1,len(self.pairable)-1
        before = bisect_right(group,index)
        while lower < upper:
            middle = (lower+upper)//2
            count = middle-index - (bisect_right(group,middle)-before)
            if count < desired:
                lower = middle+1
            else:
                upper = middle
        result = left,self.pairable[lower]
        if self.ordinal(*result) != ordinal:
            raise EvidenceStoreError("pair position does not replay")
        return result


class PairDomainPositions(_PairPositionQueries):
    """Bounded exact rank/select indexes over member and event metadata only."""
    def __init__(self, domain):
        domain._require()
        self.domain = domain
        self._require_open = domain._require
        self.member_count = domain.summary.member_count
        self.keys = [None]
        self.events = array("I",[0])
        self.indices = array("I",[2**32-1])
        self.pairable = array("I")
        groups,event_ids = [],{}
        size = 0
        for ordinal,key,_,event,parent in domain.members():
            if ordinal != len(self.keys):
                raise EvidenceStoreError("pair position membership changed")
            size += sys.getsizeof(key)+64
            self.keys.append(key)
            if event not in event_ids:
                event_ids[event] = len(groups)
                groups.append(array("I"))
                size += sys.getsizeof(event)+256
            code = event_ids[event]
            self.events.append(code)
            self.indices.append(2**32-1 if parent else len(self.pairable))
            if not parent:
                groups[code].append(len(self.pairable))
                self.pairable.append(ordinal)
            if size > MAX_POSITION_BYTES:
                raise EvidenceStoreError("pair position metadata exceeds capacity",capacity=True)
        self.groups = tuple(groups)
        self._finish_positions(domain.summary.pair_count)
        size += self.array_bytes
        if size > MAX_POSITION_BYTES:
            raise EvidenceStoreError("pair position metadata exceeds capacity",capacity=True)
        self.metadata_bytes = size


class CompactPairDomainPositions(_PairPositionQueries):
    """The same rank/select kernel over verified numeric event memberships.

    The owning case projection must establish source and outcome closure before
    supplying these rows. This index computes positions, not scientific facts.
    """
    def __init__(self, rows, *, member_count, pair_count, require_open):
        if (type(member_count) is not int or not 0 <= member_count <= MAX_CONTEXT_MEMBERS
                or type(pair_count) is not int or not 0 <= pair_count <= MAX_COMPLETE_PAIRS or not callable(require_open)):
            raise EvidenceStoreError("invalid compact pair member count")
        self.member_count,self._require_open = member_count,require_open
        self.events,self.indices,self.pairable = array("I",[0]),array("I",[2**32-1]),array("I")
        groups,event_ids = [],{}
        for ordinal,event,parent in rows:
            if (type(ordinal) is not int or ordinal != len(self.events) or ordinal > member_count
                    or type(event) is not int or not 0 <= event < 2**32 or type(parent) is not bool):
                raise EvidenceStoreError("compact pair membership changed")
            if event not in event_ids:
                event_ids[event] = len(groups)
                groups.append(array("I"))
            code = event_ids[event]
            self.events.append(code)
            self.indices.append(2**32-1 if parent else len(self.pairable))
            if not parent:
                groups[code].append(len(self.pairable))
                self.pairable.append(ordinal)
        if len(self.events) != member_count+1:
            raise EvidenceStoreError("compact pair membership is incomplete")
        self.groups = tuple(groups)
        self._finish_positions(pair_count)
        self.metadata_bytes = self.array_bytes+len(event_ids)*144
        if self.metadata_bytes > MAX_POSITION_BYTES:
            raise EvidenceStoreError("compact pair positions exceed capacity",capacity=True)
