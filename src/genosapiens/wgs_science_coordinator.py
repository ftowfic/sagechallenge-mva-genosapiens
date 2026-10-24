"""Bounded context scheduling, explicit failure scope, and coverage accounting."""
from dataclasses import dataclass
import sqlite3

from genosapiens.population_filter_candidate import OptionalPopulationFilterError
from genosapiens.track1_science_roster import ScienceRoster
from genosapiens.wgs_context_domain import ContextDomain
from genosapiens.wgs_failures import FAILURE_CODES, OPERATIONS
from genosapiens.wgs_science_checkpoints import ContextCheckpointStore, ContextIdentity, _require, _sha
from genosapiens.wgs_science_runtime import ScienceBudget, ScienceClock, ScienceDeadlineExpired, bounded_operation
from genosapiens.wgs_ledger import canonical

STAGES = frozenset({'resources', 'member_enrichment', 'member_read', 'population_query', 'population_bind',
    'clinical_query', 'clinical_bind', 'member_encode', 'member_publication', 'population_work_filter',
    'genotype_facts', 'hypothesis_evidence', 'population', 'phenotype', 'scores', 'representation',
    'eligibility', 'ordering', 'context_publication'})


@dataclass(frozen=True, slots=True)
class ContextFailure:
    stage: str
    operation: str
    code: str
    scope: str

    def __post_init__(self):
        _require(self.stage in STAGES and self.operation in OPERATIONS
            and self.code in FAILURE_CODES | {'optional_population_predicate_failed', 'scientific_deadline_expired'}
            and self.scope in {'context', 'shared', 'unknown'})


class ContextEvaluationFailure(RuntimeError):
    """An independent producer supplies a static failure and its known scope."""
    def __init__(self, failure):
        _require(type(failure) is ContextFailure)
        self.failure = failure
        super().__init__('scientific context evaluation failed')


class ScienceCoordinator:
    def __init__(self, connection, roster, *, checkpoints, computation_sha256, budget, scratch_directory,
                 maximum_contexts=None, context_timeout_seconds=600):
        _require(type(connection) is sqlite3.Connection and type(roster) is ScienceRoster
            and roster.closed and not roster.failed and type(checkpoints) is ContextCheckpointStore
            and _sha(computation_sha256) and type(budget) in (ScienceBudget, ScienceClock)
            and (maximum_contexts is None or type(maximum_contexts) is int and 0 < maximum_contexts <= 1_000_000)
            and type(context_timeout_seconds) is int and 0 < context_timeout_seconds <= 18000)
        self.connection, self.roster = connection, roster
        self.checkpoints, self.computation_sha256 = checkpoints, computation_sha256
        self.budget, self.scratch_directory, self.maximum_contexts = budget, scratch_directory, maximum_contexts
        self.context_timeout_seconds = context_timeout_seconds
        self.started = self.attempts = self.reused = self.consecutive_failures = 0
        self.previous_failure = self.stop_reason = None
        self.active = None
        connection.execute('CREATE TABLE science_context_directory('
            'ordinal INTEGER PRIMARY KEY,context_id TEXT UNIQUE,identity_sha256 TEXT UNIQUE,membership_sha256 TEXT,'
            'state TEXT,reason TEXT,checkpoint_sha256 TEXT,stage TEXT,operation TEXT,code TEXT,attempts INTEGER)')
        for row in roster.rows():
            identity = ContextIdentity(row['context_id'], row['original_ordinal'], row['membership_sha256'], computation_sha256)
            connection.execute('INSERT INTO science_context_directory VALUES (?,?,?,?,?,?,NULL,NULL,NULL,NULL,0)',
                (identity.original_ordinal, identity.context_id, identity.sha256, identity.membership_sha256,
                    'pending' if row['selected'] else 'deferred', row['reason']))
        connection.commit()

    def identity(self, ordinal):
        row = self.connection.execute('SELECT context_id,membership_sha256 FROM science_context_directory WHERE ordinal=?',
            (ordinal,)).fetchone()
        _require(row is not None)
        return ContextIdentity(row[0], ordinal, row[1], self.computation_sha256)

    def _state(self, ordinal):
        row = self.connection.execute('SELECT state FROM science_context_directory WHERE ordinal=?', (ordinal,)).fetchone()
        _require(row is not None)
        return row[0]

    def _complete(self, identity, checkpoint, *, reused):
        _require(checkpoint.identity == identity and checkpoint.case_name is not None
            and checkpoint.directory.parent == self.checkpoints.root and _sha(checkpoint.completion_sha256))
        self.connection.execute("UPDATE science_context_directory SET state='evaluated',checkpoint_sha256=? WHERE ordinal=?",
            (checkpoint.completion_sha256, identity.original_ordinal))
        self.reused += int(reused)
        self.previous_failure, self.consecutive_failures = None, 0
        self.connection.commit()
        return checkpoint

    def _failed(self, identity, failure):
        self.connection.execute("UPDATE science_context_directory SET state='failed',stage=?,operation=?,code=? WHERE ordinal=?",
            (failure.stage, failure.operation, failure.code, identity.original_ordinal))
        key = failure.stage, failure.code
        self.consecutive_failures = self.consecutive_failures + 1 if self.previous_failure == key else 1
        self.previous_failure = key
        if failure.scope != 'context':
            self.stop_reason = 'shared_failure' if failure.scope == 'shared' else 'unknown_failure_scope'
        elif self.consecutive_failures >= 3:
            self.stop_reason = 'repeated_context_failure'
        self.connection.commit()

    def run_one(self, ordinal, evaluate):
        """The callback owns a fresh writer on every call, including retry.

        It receives the bound context identity and whether optional filtering is
        enabled. It returns only after durable context publication.
        """
        _require(self.active is None and self._state(ordinal) == 'pending')
        identity = self.identity(ordinal)
        # Compatible completed work remains usable after the scheduling budget
        # expires. Loading cannot grant another evaluation or reset the clock.
        try:
            checkpoint = self.checkpoints.load(identity, scratch_directory=self.scratch_directory)
        except BaseException:
            self.stop_reason = 'checkpoint_authentication_failed'
            raise
        if checkpoint is not None:
            return self._complete(identity, checkpoint, reused=True)
        if self.stop_reason is not None:
            return None
        if self.budget.remaining() <= 0:
            self.stop_reason = 'scientific_deadline'
            return None
        if self.maximum_contexts is not None and self.started >= self.maximum_contexts:
            self.stop_reason = 'pilot_context_limit'
            return None
        self.started += 1
        self.active = identity
        try:
            for attempt in range(2):
                if self.budget.remaining() <= 0:
                    self.stop_reason = 'scientific_deadline'
                    return None
                self.attempts += 1
                self.connection.execute('UPDATE science_context_directory SET attempts=attempts+1 WHERE ordinal=?', (ordinal,))
                try:
                    with bounded_operation(self.budget, maximum_seconds=self.context_timeout_seconds):
                        checkpoint = evaluate(identity, optional_filter=attempt == 0)
                    # A live store can attest that its own writer just finished
                    # this artifact. A callback's returned object alone cannot.
                    durable = self.checkpoints.accept_fresh(identity, checkpoint)
                    if durable is None:
                        durable = self.checkpoints.load(identity, scratch_directory=self.scratch_directory)
                    _require(durable is not None and durable == checkpoint)
                    return self._complete(identity, durable, reused=False)
                except OptionalPopulationFilterError:
                    if attempt == 0:
                        self.connection.execute("UPDATE science_context_directory SET reason='optional_filter_retry' WHERE ordinal=?", (ordinal,))
                        continue
                    self._failed(identity, ContextFailure('population_work_filter', 'population_assess',
                        'optional_population_predicate_failed', 'context'))
                    return None
                except ContextEvaluationFailure as error:
                    self._failed(identity, error.failure)
                    return None
                except ScienceDeadlineExpired:
                    if self.budget.remaining() <= 0:
                        self.stop_reason = 'scientific_deadline'
                    else:
                        # An unobserved operation timeout has no proven local
                        # scope. It cannot imply an exhausted case deadline.
                        self._failed(identity, ContextFailure('context_publication', 'outermost_unknown',
                            'bounded_products_timeout', 'unknown'))
                    return None
                except BaseException:
                    self._failed(identity, ContextFailure('context_publication', 'outermost_unknown',
                        'unexpected_unexpected', 'unknown'))
                    raise
        finally:
            self.active = None
            self.connection.commit()

    def counts(self):
        counts = dict(evaluated=0, deferred=0, failed=0, pending=0)
        counts.update(dict(self.connection.execute('SELECT state,count(*) FROM science_context_directory GROUP BY state')))
        counts.update(original=self.roster.counts['original'], selected=self.roster.counts['selected'],
            checkpointed=counts['evaluated'], reused=self.reused, attempted=self.started, attempts=self.attempts)
        _require(counts['original'] == sum(counts[k] for k in ('evaluated', 'deferred', 'failed', 'pending'))
            and counts['selected'] == sum(counts[k] for k in ('evaluated', 'failed', 'pending')))
        return counts

    def completed_domain(self):
        return ContextDomain.from_ordinals(row[0] for row in self.connection.execute(
            "SELECT ordinal FROM science_context_directory WHERE state='evaluated' ORDER BY ordinal"))

    def rows(self):
        keys = ('original_ordinal', 'context_id', 'identity_sha256', 'membership_sha256', 'state', 'reason',
            'checkpoint_sha256', 'stage', 'operation', 'code', 'attempts')
        for values in self.connection.execute('SELECT * FROM science_context_directory ORDER BY ordinal'):
            row = dict(zip(keys, values, strict=True))
            _require(len(canonical(row)) <= 16 * 1024)
            yield row

    def status(self):
        counts = self.counts()
        return dict(schema='genosapiens.science-context-status/v1', scope_sha256=self.roster.sha256,
            policy_sha256=self.roster.policy.sha256, budget_sha256=self.budget.sha256, counts=counts,
            status=('incomplete' if counts['failed'] or counts['pending'] else
                'rapid_complete' if self.roster.policy.rapid else 'exhaustive_complete'), stop_reason=self.stop_reason)
