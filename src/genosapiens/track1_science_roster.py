"""Disk-backed whole-context selection after complete upstream closure.

The caller supplies projections from authenticated owners. This module cannot
authenticate arbitrary source claims. Selection never means scientific exclusion.
"""
from dataclasses import dataclass
import hashlib
import sqlite3

from genosapiens.track1_science_scope import ScienceScopeMode, ScienceScopePolicy
from genosapiens.wgs_ledger import canonical

MAX_CONTEXTS = 1_000_000
MAX_MEMBERS = 100_000
STRONG = {
    'gencc_submissions': {'strong', 'definitive'},
    'evidence_clingen': {'strong', 'definitive'},
    'evidence_gene2phenotype': {'strong', 'definitive'},
    'evidence_genomics_england': {'green'},
    'orphadata_product6': {'validated'},
    'evidence_orphanet': {'validated'},
}
WEAK = {'moderate', 'limited', 'supportive', 'amber'}


def strict_association_category(dataset, category, *, direct, positive, causal, exact, conflict):
    """Apply source-specific categories after the source owner resolves scope."""
    if any(type(v) is not bool for v in (direct, positive, causal, exact, conflict)):
        raise ValueError('invalid scientific association state')
    if type(dataset) is not str or category is not None and type(category) is not str:
        raise ValueError('invalid scientific association source')
    label = '' if category is None else category.strip().casefold()
    if conflict or not all((direct, positive, causal, exact)):
        return 'nonqualifying'
    if dataset in STRONG and label in STRONG[dataset]:
        return 'strong'
    if dataset in STRONG and label in WEAK:
        return 'weak'
    return 'unknown'


@dataclass(frozen=True, slots=True)
class ContextProjection:
    context_id: str
    context_sha256: str
    phenotype_units: int
    phenotype_available: bool
    molecular_protected: bool = False
    protection_uncertain: bool = False
    association: str = 'unknown'


def _sha(value):
    return type(value) is str and len(value) == 64 and all(c in '0123456789abcdef' for c in value)


class ScienceRoster:
    """Merge all fragments before ranking, including late molecular exceptions."""
    def __init__(self, connection, policy):
        if type(connection) is not sqlite3.Connection or type(policy) is not ScienceScopePolicy:
            raise ValueError('invalid scientific roster owner')
        self.connection, self.policy, self.closed, self.failed = connection, policy, False, False
        self.fragments = self.context_count = 0
        connection.executescript('''
            CREATE TABLE science_roster(
                context_id TEXT PRIMARY KEY, context_sha256 TEXT NOT NULL,
                score INTEGER NOT NULL, available INTEGER NOT NULL,
                protected INTEGER NOT NULL, uncertain INTEGER NOT NULL,
                strong INTEGER NOT NULL, weak INTEGER NOT NULL, nonqualifying INTEGER NOT NULL,
                ordinal INTEGER UNIQUE, member_count INTEGER, membership_sha256 TEXT,
                selected INTEGER, reason TEXT);
            CREATE TABLE science_roster_members(
                context_id TEXT NOT NULL, identity TEXT NOT NULL, sha256 TEXT NOT NULL,
                PRIMARY KEY(context_id,identity));
        ''')

    def add(self, projection, members):
        if self.failed:
            raise ValueError('scientific roster has failed')
        try:
            self._add(projection, members)
        except BaseException:
            self.failed = True
            raise

    def _add(self, projection, members):
        if self.closed or type(projection) is not ContextProjection:
            raise ValueError('scientific roster is closed or projection is invalid')
        p = projection
        if (type(p.context_id) is not str or not 0 < len(p.context_id) <= 256
                or not _sha(p.context_sha256) or type(p.phenotype_units) is not int
                or not -(2**63) < p.phenotype_units < 2**63
                or any(type(v) is not bool for v in
                       (p.phenotype_available, p.molecular_protected, p.protection_uncertain))
                or p.association not in {'strong', 'weak', 'unknown', 'nonqualifying'}):
            raise ValueError('invalid scientific context projection')
        prior = self.connection.execute('SELECT context_sha256,score,available FROM science_roster WHERE context_id=?',
                                        (p.context_id,)).fetchone()
        identity = (p.context_sha256, p.phenotype_units, int(p.phenotype_available))
        if prior is not None and prior != identity:
            raise ValueError('scientific context fragments disagree')
        if prior is None:
            if self.context_count >= MAX_CONTEXTS:
                raise ValueError('scientific roster exceeds context capacity')
            self.connection.execute('INSERT INTO science_roster VALUES (?,?,?,?,?,?,?,?,?,NULL,NULL,NULL,NULL,NULL)',
                (p.context_id, *identity, p.molecular_protected, p.protection_uncertain,
                 p.association == 'strong', p.association == 'weak', p.association == 'nonqualifying'))
            self.context_count += 1
        else:
            self.connection.execute('UPDATE science_roster SET protected=max(protected,?),uncertain=max(uncertain,?),'
                'strong=max(strong,?),weak=max(weak,?),nonqualifying=max(nonqualifying,?) WHERE context_id=?',
                (p.molecular_protected, p.protection_uncertain, p.association == 'strong',
                 p.association == 'weak', p.association == 'nonqualifying', p.context_id))
        count = 0
        for identity, digest in members:
            if type(identity) is not str or not 0 < len(identity) <= 1024 or not _sha(digest):
                raise ValueError('invalid scientific membership identity')
            prior = self.connection.execute('SELECT sha256 FROM science_roster_members WHERE context_id=? AND identity=?',
                                            (p.context_id, identity)).fetchone()
            if prior is not None and prior != (digest,):
                raise ValueError('scientific membership identity changed')
            self.connection.execute('INSERT OR IGNORE INTO science_roster_members VALUES (?,?,?)',
                                    (p.context_id, identity, digest))
            count += 1
            if count > MAX_MEMBERS:
                raise ValueError('scientific fragment exceeds member capacity')
        self.fragments += 1

    def close(self):
        if self.closed or self.failed:
            raise ValueError('scientific roster already closed')
        try:
            return self._close()
        except BaseException:
            self.failed = True
            self.closed = False
            raise

    def _close(self):
        total, = self.connection.execute('SELECT count(*) FROM science_roster').fetchone()
        if total > MAX_CONTEXTS:
            raise ValueError('scientific roster exceeds context capacity')
        for ordinal, (key,) in enumerate(self.connection.execute('SELECT context_id FROM science_roster ORDER BY context_id'), 1):
            digest, count = hashlib.sha256(), 0
            for row in self.connection.execute('SELECT identity,sha256 FROM science_roster_members WHERE context_id=? ORDER BY identity', (key,)):
                digest.update(canonical(row))
                count += 1
                if count > MAX_MEMBERS:
                    raise ValueError('scientific context exceeds member capacity')
            if count == 0:
                raise ValueError('scientific context lacks members')
            self.connection.execute('UPDATE science_roster SET ordinal=?,member_count=?,membership_sha256=? WHERE context_id=?',
                (ordinal, count, digest.hexdigest(), key))
        if self.policy.mode is ScienceScopeMode.EXHAUSTIVE:
            self.connection.execute("UPDATE science_roster SET selected=1,reason='exhaustive'")
        else:
            self.connection.execute("UPDATE science_roster SET selected=0,reason='phenotype_deferred'")
            if self.policy.mode is ScienceScopeMode.RAPID_STRICT_ASSOCIATION:
                self.connection.execute("UPDATE science_roster SET reason=CASE WHEN weak=1 THEN 'strict_association_weak' "
                    "WHEN nonqualifying=1 THEN 'strict_association_nonqualifying' ELSE 'strict_association_unknown' END WHERE strong=0")
            self.connection.execute("UPDATE science_roster SET selected=1,reason='phenotype_missing' WHERE available=0 AND reason='phenotype_deferred'")
            self.connection.execute("UPDATE science_roster SET selected=1,reason='protection_uncertain' WHERE uncertain=1")
            self.connection.execute("UPDATE science_roster SET selected=1,reason='protected_molecular' WHERE protected=1")
            protected, = self.connection.execute('SELECT count(*) FROM science_roster WHERE selected=1').fetchone()
            remaining = self.policy.target_contexts - protected
            if remaining > 0:
                boundary = self.connection.execute("SELECT score FROM science_roster WHERE selected=0 AND reason='phenotype_deferred' "
                    'ORDER BY score DESC,context_id LIMIT 1 OFFSET ?', (remaining - 1,)).fetchone()
                if boundary is None:
                    self.connection.execute("UPDATE science_roster SET selected=1,reason='phenotype_selected' WHERE reason='phenotype_deferred'")
                else:
                    self.connection.execute("UPDATE science_roster SET selected=1,reason='phenotype_selected' WHERE reason='phenotype_deferred' AND score>=?", boundary)
        self.connection.commit()
        self.closed = True
        self.counts = dict(original=total, selected=self.connection.execute('SELECT count(*) FROM science_roster WHERE selected=1').fetchone()[0])
        self.counts['deferred'] = total - self.counts['selected']
        self.counts['overflow'] = max(0, self.counts['selected'] - self.policy.target_contexts) if self.policy.rapid else 0
        digest = hashlib.sha256(self.policy.as_bytes())
        for row in self.rows():
            digest.update(canonical(row))
        self.sha256 = digest.hexdigest()
        return self

    def rows(self):
        if not self.closed:
            raise ValueError('scientific roster is not globally closed')
        columns = ('context_id', 'context_sha256', 'phenotype_units', 'phenotype_available',
            'molecular_protected', 'protection_uncertain', 'strong_association', 'weak_association',
            'nonqualifying_association', 'original_ordinal', 'member_count', 'membership_sha256', 'selected', 'reason')
        for values in self.connection.execute('SELECT * FROM science_roster ORDER BY ordinal'):
            yield dict(zip(columns, values, strict=True))

    def selected(self, context_id):
        if not self.closed:
            raise ValueError('scientific roster is not globally closed')
        row = self.connection.execute('SELECT selected FROM science_roster WHERE context_id=?', (context_id,)).fetchone()
        if row is None:
            raise ValueError('scientific context is outside the roster')
        return bool(row[0])
