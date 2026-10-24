"""Replay explicit scope accounting separately from evaluated scientific scores."""
import hashlib
import json

from genosapiens.track1_science_scope import ScienceScopeMode, ScienceScopePolicy
from genosapiens.wgs_case_resources import bound_reference
from genosapiens.wgs_completed_science import child, fields, integer, require
from genosapiens.wgs_evidence_store import EvidenceStoreError
from genosapiens.wgs_failures import FAILURE_CODES, OPERATIONS
from genosapiens.wgs_ledger import canonical
from genosapiens.wgs_science_checkpoints import ContextIdentity, _sha, _verify_membership
from genosapiens.wgs_science_coordinator import STAGES
from genosapiens.wgs_science_runtime import SCHEMA as BUDGET_SCHEMA, CLOCK_SCHEMA, ScienceBudget, ScienceClock

ROSTER_KEYS = ('context_id', 'context_sha256', 'phenotype_units', 'phenotype_available', 'molecular_protected',
    'protection_uncertain', 'strong_association', 'weak_association', 'nonqualifying_association',
    'original_ordinal', 'member_count', 'membership_sha256', 'selected', 'reason')
DIRECTORY_KEYS = ('original_ordinal', 'context_id', 'identity_sha256', 'membership_sha256', 'state', 'reason',
    'checkpoint_sha256', 'stage', 'operation', 'code', 'attempts')
COUNT_KEYS = ('original', 'selected', 'evaluated', 'deferred', 'failed', 'pending', 'checkpointed', 'reused', 'attempted', 'attempts')


class VerifiedScientificScope:
    def __init__(self, reader, artifact, case):
        self.connection, self.closed = reader.connection, False
        self.table = 'verified_scientific_scope'
        self.created = False
        try:
            reference = reader.reference(artifact.scope_name, kind='scientific-scope-v1')
            require(reference.sha256 == artifact.scope_sha256)
            value = self.header = fields(reader.reconstruct(reference), ('schema', 'scope_sha256', 'policy_sha256',
                'budget_sha256', 'counts', 'status', 'stop_reason', 'policy', 'budget', 'computation_sha256',
                'resources', 'protected_overflow', 'roster', 'contexts'))
            require(value['schema'] == 'genosapiens.science-context-status/v1' and value['resources'] == case.resources.reference
                and _sha(value['computation_sha256']) and _sha(value['scope_sha256']))
            wire = value['policy']
            self.policy = ScienceScopePolicy(ScienceScopeMode(wire['mode']),
                1500 if wire['mode'] == 'exhaustive' else wire['target_contexts'])
            require(canonical(wire) == self.policy.as_bytes() and value['policy_sha256'] == self.policy.sha256)
            budget = dict(value['budget'])
            runtime_schema = budget.pop('schema')
            require(runtime_schema in (BUDGET_SCHEMA, CLOCK_SCHEMA))
            self.budget = (ScienceClock if runtime_schema == CLOCK_SCHEMA else ScienceBudget)(**budget)
            require(self.budget.sha256 == value['budget_sha256'])
            counts = fields(value['counts'], COUNT_KEYS)
            for number in counts.values():
                integer(number, maximum=2_000_000)
            self.connection.execute(f'CREATE TEMP TABLE {self.table}(ordinal INTEGER PRIMARY KEY,id TEXT UNIQUE,'
                'member_sha TEXT,score INTEGER,available INTEGER,protected INTEGER,uncertain INTEGER,strong INTEGER,'
                'weak INTEGER,nonqualifying INTEGER,actual_selected INTEGER,actual_reason TEXT,selected INTEGER,reason TEXT,state TEXT)')
            self.created = True
            digest, count, previous = hashlib.sha256(self.policy.as_bytes()), 0, None
            for row in reader.iter_array(child(reader, value['roster'], 'scientific-scope-roster-v1'), max_item_bytes=16 * 1024):
                fields(row, ROSTER_KEYS)
                count += 1
                require(count <= 1_000_000 and integer(row['original_ordinal'], minimum=1) == count
                    and type(row['context_id']) is str and 0 < len(row['context_id']) <= 256
                    and (previous is None or previous < row['context_id'])
                    and _sha(row['context_sha256']) and _sha(row['membership_sha256'])
                    and type(row['phenotype_units']) is int and -(2**63) < row['phenotype_units'] < 2**63)
                previous = row['context_id']
                integer(row['member_count'], minimum=1, maximum=100_000)
                for key in ('phenotype_available', 'molecular_protected', 'protection_uncertain', 'strong_association',
                            'weak_association', 'nonqualifying_association', 'selected'):
                    integer(row[key], maximum=1)
                digest.update(canonical(row))
                self.connection.execute(f'INSERT INTO {self.table} VALUES (?,?,?,?,?,?,?,?,?,?,?,?,0,NULL,NULL)',
                    (count, row['context_id'], row['membership_sha256'], row['phenotype_units'], row['phenotype_available'],
                        row['molecular_protected'], row['protection_uncertain'], row['strong_association'],
                        row['weak_association'], row['nonqualifying_association'], row['selected'], row['reason']))
            require(count == counts['original'] and digest.hexdigest() == value['scope_sha256'])
            self._replay_selection()
            require(self.connection.execute(f'SELECT count(*) FROM {self.table} WHERE selected!=actual_selected OR reason!=actual_reason').fetchone() == (0,))
            selected = self.connection.execute(f'SELECT count(*) FROM {self.table} WHERE selected=1').fetchone()[0]
            require(selected == counts['selected'] and integer(value['protected_overflow']) ==
                (max(0, selected - self.policy.target_contexts) if self.policy.rapid else 0))
            observed = dict.fromkeys(COUNT_KEYS, 0)
            for row in reader.iter_array(child(reader, value['contexts'], 'scientific-context-directory-v1'), max_item_bytes=16 * 1024):
                fields(row, DIRECTORY_KEYS)
                observed['original'] += 1
                ordinal = integer(row['original_ordinal'], minimum=1)
                require(ordinal == observed['original'])
                expected = self.connection.execute(f'SELECT id,member_sha,selected,reason FROM {self.table} WHERE ordinal=?', (ordinal,)).fetchone()
                require(expected is not None and tuple(row[k] for k in ('context_id', 'membership_sha256')) == expected[:2])
                identity = ContextIdentity(row['context_id'], ordinal, row['membership_sha256'], value['computation_sha256'])
                require(row['identity_sha256'] == identity.sha256)
                attempts = integer(row['attempts'], maximum=2)
                require(row['reason'] == expected[3] or row['reason'] == 'optional_filter_retry' and attempts == 2)
                state = row['state']
                require(state in {'evaluated', 'failed', 'pending', 'deferred'} and (state != 'deferred') == bool(expected[2]))
                observed[state] += 1
                observed['selected'] += expected[2]
                observed['attempted'] += attempts > 0
                observed['attempts'] += attempts
                if state == 'evaluated':
                    require(_sha(row['checkpoint_sha256']))
                    completed = case.context(ordinal)
                    require(completed['context_id'] == row['context_id'])
                    context = reader.reconstruct(bound_reference(reader, completed['context'], 'scientific-case-context-v1'))
                    _verify_membership(reader, context, row['membership_sha256'])
                    observed['checkpointed'] += 1
                    observed['reused'] += attempts == 0
                else:
                    require(row['checkpoint_sha256'] is None)
                if state == 'failed':
                    require(attempts > 0 and row['stage'] in STAGES and row['operation'] in OPERATIONS
                        and row['code'] in FAILURE_CODES | {'optional_population_predicate_failed', 'scientific_deadline_expired'})
                else:
                    require(all(row[key] is None for key in ('stage', 'operation', 'code')))
                require(state != 'deferred' or attempts == 0)
                self.connection.execute(f'UPDATE {self.table} SET state=? WHERE ordinal=?', (state, ordinal))
            require(observed == counts and counts['evaluated'] == case.context_count)
            status = 'incomplete' if counts['failed'] or counts['pending'] else 'rapid_complete' if self.policy.rapid else 'exhaustive_complete'
            require(value['status'] == status == artifact.completion_status)
            require(value['stop_reason'] in {None, 'shared_failure', 'unknown_failure_scope', 'repeated_context_failure',
                'scientific_deadline', 'pilot_context_limit', 'checkpoint_authentication_failed'})
        except BaseException:
            self.close()
            raise

    def _replay_selection(self):
        t, c = self.table, self.connection
        if not self.policy.rapid:
            c.execute(f"UPDATE {t} SET selected=1,reason='exhaustive'")
            return
        c.execute(f"UPDATE {t} SET reason='phenotype_deferred'")
        if self.policy.mode is ScienceScopeMode.RAPID_STRICT_ASSOCIATION:
            c.execute(f"UPDATE {t} SET reason=CASE WHEN weak=1 THEN 'strict_association_weak' WHEN nonqualifying=1 "
                "THEN 'strict_association_nonqualifying' ELSE 'strict_association_unknown' END WHERE strong=0")
        c.execute(f"UPDATE {t} SET selected=1,reason='phenotype_missing' WHERE available=0 AND reason='phenotype_deferred'")
        c.execute(f"UPDATE {t} SET selected=1,reason='protection_uncertain' WHERE uncertain=1")
        c.execute(f"UPDATE {t} SET selected=1,reason='protected_molecular' WHERE protected=1")
        remaining = self.policy.target_contexts - c.execute(f'SELECT count(*) FROM {t} WHERE selected=1').fetchone()[0]
        if remaining > 0:
            boundary = c.execute(f"SELECT score FROM {t} WHERE reason='phenotype_deferred' ORDER BY score DESC,id LIMIT 1 OFFSET ?",
                (remaining - 1,)).fetchone()
            if boundary is None:
                c.execute(f"UPDATE {t} SET selected=1,reason='phenotype_selected' WHERE reason='phenotype_deferred'")
            else:
                c.execute(f"UPDATE {t} SET selected=1,reason='phenotype_selected' WHERE reason='phenotype_deferred' AND score>=?", boundary)

    def close(self):
        if self.created:
            self.connection.execute(f'DROP TABLE {self.table}')
            self.created = False
        self.closed = True

    def __enter__(self):
        require(not self.closed)
        return self

    def __exit__(self, *_):
        self.close()


def verify_scientific_scope(reader, artifact, case):
    return VerifiedScientificScope(reader, artifact, case)
