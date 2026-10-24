"""Operation-owned public population policy; bounded construction before the DAG.

Keys close globally before policy consumption. Canonical row order and the exact
original logical policy digest survive without assembling a whole policy JSON
document or retaining all rich ContextPolicy objects. This private SQLite input
is never a saved authority or a second scientific evaluator.
"""
import hashlib
import json
from pathlib import Path
import sqlite3
import tempfile

from genosapiens.candidate_admission import _bytes
from genosapiens.evidence_model import _typed
from genosapiens.population_policy import (
    ContextKey, ContextPolicy, FrequencyRule, PenaltyTransform, RescueRules,
    PopulationPolicy, policy_as_json,
)
from genosapiens.wgs_evidence_store import EvidenceStoreError
from genosapiens.wgs_population_policy import _LogicalDigest

MAX_POLICY_KEYS = 1_000_000
MAX_INPUT_BYTES = 64 * 1024
POLICY_ID = "track1-public-development-population-v2"


class GeneratedPopulationPolicy:
    """A complete, immutable live-operation source for normalized publication."""
    def __init__(self, *, scratch_directory, rule, transform, rescues):
        root = Path(scratch_directory)
        if not root.is_dir() or any(p.is_symlink() for p in (root, *root.parents)):
            raise EvidenceStoreError("population input scratch is unavailable")
        self.rule = _typed(rule, FrequencyRule)
        self.transform = _typed(transform, PenaltyTransform)
        self.rescues = _typed(rescues, RescueRules)
        if len(_bytes((self.rule, self.transform, self.rescues))) > MAX_INPUT_BYTES:
            raise EvidenceStoreError("population input template exceeds capacity", capacity=True)
        self.policy_id = POLICY_ID
        self.temporary = tempfile.TemporaryDirectory(prefix=".population-input-", dir=root)
        self.connection = sqlite3.connect(Path(self.temporary.name) / "keys.sqlite")
        self.connection.execute("PRAGMA cache_size=-65536")
        self.connection.execute("PRAGMA temp_store=FILE")
        self.connection.execute("PRAGMA mmap_size=0")
        self.connection.execute("CREATE TABLE keys(content BLOB PRIMARY KEY) WITHOUT ROWID")
        self.context_count = 0
        self.closed = self.sealed = self.failed = False

    def _require(self):
        if self.closed or self.failed:
            raise EvidenceStoreError("population input owner is unavailable")

    def add(self, key):
        self._require()
        if self.sealed:
            raise EvidenceStoreError("population input is already closed")
        try:
            content = _bytes(_typed(key, ContextKey))
            if len(content) > MAX_INPUT_BYTES:
                raise EvidenceStoreError("population input key exceeds capacity", capacity=True)
            inserted = self.connection.execute("INSERT OR IGNORE INTO keys VALUES (?)", (content,)).rowcount
            self.context_count += inserted
            if self.context_count > MAX_POLICY_KEYS:
                raise EvidenceStoreError("population input key count exceeds capacity", capacity=True)
        except BaseException:
            self.failed = True
            raise

    def _context(self, key):
        # Same ContextPolicy constructor and literal settings as the original
        # public builder. Scientific frequency evaluation remains unchanged.
        return ContextPolicy(key, (self.rule,), self.transform, self.rescues,
            "maximum", "maximum", "strict_greater_than",
            "AF_exact_AC_over_AN;FAF_decimal_typed_value", "neutral",
            "neutral_not_excluded", "0", "any_reliable_excess_member_without_rescue")

    def _rows(self):
        self._require()
        digest, count = hashlib.sha256(), 0
        cursor = self.connection.execute("SELECT content FROM keys ORDER BY content")
        try:
            for (content,) in cursor:
                if type(content) is not bytes or len(content) > MAX_INPUT_BYTES:
                    raise EvidenceStoreError("population input key changed")
                key = _typed(json.loads(content), ContextKey, decode=True)
                if _bytes(key) != content:
                    raise EvidenceStoreError("population input is noncanonical")
                digest.update(content)
                count += 1
                yield self._context(key)
            if count != self.context_count or self.sealed and digest.hexdigest() != self.key_sha256:
                raise EvidenceStoreError("population input coverage changed")
        except BaseException:
            self.failed = True
            raise
        finally:
            cursor.close()

    def seal(self):
        self._require()
        if self.sealed:
            raise EvidenceStoreError("population input was already sealed")
        try:
            self.connection.commit()
            self.connection.execute("PRAGMA query_only=ON")
            keys = hashlib.sha256()
            for (content,) in self.connection.execute("SELECT content FROM keys ORDER BY content"):
                keys.update(content)
            self.key_sha256 = keys.hexdigest()
            logical = _LogicalDigest(self.policy_id)
            self.logical_bytes = len(policy_as_json(PopulationPolicy(self.policy_id, (), ())))
            for ordinal, context in enumerate(self._rows()):
                logical.add_context(context)
                self.logical_bytes += len(_bytes(context)) - 1 + int(ordinal != 0)
            self.policy_sha256 = logical.finish()
            self.sealed = True
            return self
        except BaseException:
            self.failed = True
            raise

    def iter_contexts(self):
        self._require()
        if not self.sealed:
            raise EvidenceStoreError("population input has not globally closed")
        yield from self._rows()

    def iter_exceptions(self):
        self._require()
        if not self.sealed:
            raise EvidenceStoreError("population input has not globally closed")
        # The current public builder generates no curated exception rows. Native
        # and clinical rescue semantics still live in the exact shared template.
        return iter(())

    def close(self):
        if not self.closed:
            self.connection.close()
            self.temporary.cleanup()
            self.closed = True

    def __enter__(self):
        self._require()
        return self

    def __exit__(self, *_):
        self.close()
