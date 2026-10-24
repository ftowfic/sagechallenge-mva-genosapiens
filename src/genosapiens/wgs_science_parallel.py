"""Ordered coordination of persistent, independently publishing contexts."""
from dataclasses import asdict
from pathlib import Path
import os
import time
from itertools import islice

from genosapiens.wgs_science_checkpoints import ContextIdentity
from genosapiens.wgs_science_coordinator import ContextFailure, ContextEvaluationFailure
from genosapiens.wgs_science_runtime import bounded_operation, ScienceDeadlineExpired
from genosapiens.population_filter_candidate import OptionalPopulationFilterError


def run_parallel(contexts, *, coordinator, initialize, evaluate, workers, root, updated):
    from genosapiens.wgs_science_workers import PersistentScienceWorkers
    from genosapiens.wgs_contexts import _bytes
    source = iter(contexts)
    try:
        with bounded_operation(coordinator.budget, maximum_seconds=coordinator.context_timeout_seconds):
            first = next(source, None)
    except ScienceDeadlineExpired:
        coordinator.stop_reason = ('scientific_deadline' if coordinator.budget.remaining() <= 0
            else 'unknown_failure_scope')
        return dict(workers=workers, startup_seconds=[], max_outstanding=0)
    if first is None:
        return dict(workers=workers, startup_seconds=[], max_outstanding=0)
    # Exhausting the producing generator deletes its context index. Keep that
    # owner suspended until every worker releases its independent read handle.
    source_owner = source
    source = islice(source_owner, coordinator.roster.counts['selected'] - 1)
    input_started = time.monotonic()
    try:
        with bounded_operation(coordinator.budget, maximum_seconds=coordinator.context_timeout_seconds):
            initialize.prepare(first)
    except ScienceDeadlineExpired:
        coordinator.stop_reason = ('scientific_deadline' if coordinator.budget.remaining() <= 0
            else 'unknown_failure_scope')
        return dict(workers=workers, startup_seconds=[], max_outstanding=0)
    input_preparation_seconds = time.monotonic()-input_started
    # The parent publishes the shared bundle before concurrent context writers.
    from genosapiens.wgs_evidence_shared import CheckpointResourcePrefix
    store = coordinator.checkpoints
    if store._resource_bundle is None:
        store._resource_bundle = CheckpointResourcePrefix(store.root, store._produced_prefix)
    descriptors = Path(root) / 'worker-inputs'
    descriptors.mkdir(mode=0o700)

    def worker_setup(stack):
        load_context = initialize(first, stack)
        def job(request):
            identity = ContextIdentity(**request['identity'])
            closed = load_context(request)
            timings, attempts = {}, 0
            start = time.monotonic()
            result = dict(member_count=closed.member_count)
            for attempt in range(2):
                attempts += 1
                try:
                    with bounded_operation(coordinator.budget, maximum_seconds=coordinator.context_timeout_seconds):
                        checkpoint = evaluate(closed, identity, optional_filter=attempt == 0, timings=timings)
                    result['completion_sha256'] = checkpoint.completion_sha256
                    break
                except OptionalPopulationFilterError:
                    if attempt == 0:
                        continue
                    result['failure'] = asdict(ContextFailure('population_work_filter', 'population_assess',
                        'optional_population_predicate_failed', 'context'))
                except ContextEvaluationFailure as error:
                    result['failure'] = asdict(error.failure)
                    break
                except ScienceDeadlineExpired:
                    if coordinator.budget.remaining() <= 0:
                        result['deadline'] = True
                    else:
                        result['failure'] = asdict(ContextFailure('context_publication', 'outermost_unknown',
                            'bounded_products_timeout', 'unknown'))
                    break
                except (KeyboardInterrupt, SystemExit):
                    raise
                except BaseException as error:
                    from genosapiens.wgs_science_failures import observe_setup_failure
                    observe_setup_failure(error)
                    result['failure'] = asdict(ContextFailure('context_publication', 'outermost_unknown',
                        'unexpected_unexpected', 'unknown'))
                    break
            seconds = time.monotonic()-start
            timings['coordinator_checks'] = max(0, seconds-sum(timings.values()))
            result.update(attempts=attempts, phase_seconds=timings, seconds=seconds,
                worker_pid=os.getpid())
            return result
        job.initialization_metrics = getattr(load_context, 'initialization_metrics', {})
        return job

    buffered = first
    exhausted = False
    with PersistentScienceWorkers(workers, worker_setup,
            startup_timeout=min(coordinator.context_timeout_seconds, coordinator.budget.remaining())) as pool:
        while not exhausted or pool.pending:
            pool.poll()
            if pool.fatal is not None:
                coordinator.stop_reason = pool.fatal
            while not exhausted and len(pool.pending) < workers and coordinator.stop_reason is None:
                closed = buffered
                buffered = None
                if closed is None:
                    try:
                        with bounded_operation(coordinator.budget, maximum_seconds=coordinator.context_timeout_seconds):
                            closed = next(source)
                    except StopIteration:
                        exhausted = True
                        break
                    except ScienceDeadlineExpired:
                        coordinator.stop_reason = ('scientific_deadline' if coordinator.budget.remaining() <= 0
                            else 'unknown_failure_scope')
                        break
                identity = coordinator.identity(closed.ordinal)
                checkpoint = store.load(identity, scratch_directory=coordinator.scratch_directory)
                # Drain older work before recording reuse, preserving failure order.
                if checkpoint is not None and pool.pending:
                    buffered = closed
                    break
                if checkpoint is not None:
                    coordinator._complete(identity, checkpoint, reused=True)
                    updated(None)
                    continue
                if coordinator.budget.remaining() <= 0:
                    coordinator.stop_reason = 'scientific_deadline'
                    break
                if coordinator.maximum_contexts is not None and coordinator.started >= coordinator.maximum_contexts:
                    coordinator.stop_reason = 'pilot_context_limit'
                    break
                path = descriptors / (str(closed.ordinal) + '.json')
                content = _bytes(closed.seed)
                if len(content) > closed._config.contexts.max_context_bytes:
                    raise ValueError('scientific context descriptor exceeds capacity')
                with path.open('xb') as stream:
                    stream.write(content)
                coordinator.started += 1
                pool.submit(dict(ordinal=closed.ordinal, seed_path=str(path), identity=asdict(identity)),
                    timeout=2*coordinator.context_timeout_seconds+30)
            if pool.pending:
                try:
                    ordinal, result = pool.take()
                except RuntimeError:
                    coordinator.stop_reason = pool.fatal or 'unknown_failure_scope'
                    break
                identity = coordinator.identity(ordinal)
                coordinator.attempts += result['attempts']
                coordinator.connection.execute('UPDATE science_context_directory SET attempts=? WHERE ordinal=?',
                    (result['attempts'], ordinal))
                if result['attempts'] == 2:
                    coordinator.connection.execute("UPDATE science_context_directory SET reason='optional_filter_retry' WHERE ordinal=?", (ordinal,))
                if 'completion_sha256' in result:
                    durable = store._read(store._destination(identity), identity)
                    if durable.completion_sha256 != result['completion_sha256']:
                        raise ValueError('scientific worker completion differs')
                    store.omitted_fresh_replays += 1
                    coordinator._complete(identity, durable, reused=False)
                elif 'failure' in result:
                    coordinator._failed(identity, ContextFailure(**result['failure']))
                elif result.get('deadline'):
                    coordinator.stop_reason = 'scientific_deadline'
                else:
                    raise ValueError('scientific worker result is incomplete')
                result['state'] = coordinator._state(ordinal)
                updated(result)
                (descriptors / (str(ordinal) + '.json')).unlink()
            if coordinator.stop_reason is not None:
                break
        coordinator.connection.commit()
        metrics = pool.metrics
        metrics['shared_input_preparation_seconds'] = input_preparation_seconds
    if coordinator.stop_reason is None and next(source_owner, None) is not None:
        raise ValueError('scientific context stream exceeds selected domain')
    return metrics
