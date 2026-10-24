"""Independent scientific context production and bounded case aggregation.

Factories are internal producing operations. They must use the ordinary source,
resource, and scientific owners. This module cannot authorize a submission.
"""
from contextlib import ExitStack, contextmanager, nullcontext
from dataclasses import dataclass
import hashlib
import json
import os
import time
from pathlib import Path

from genosapiens.wgs_case_evidence import CompleteScientificCase, CompletedScientificCase
from genosapiens.wgs_case_resources import CompletedScientificCaseResources
from genosapiens.wgs_context_domain import ContextDomain
from genosapiens.wgs_evidence_store import EvidenceReader, EvidenceWriter
from genosapiens.wgs_ledger import canonical
from genosapiens.wgs_science_aggregation import append_completed_context
from genosapiens.wgs_science_checkpoints import _path, _require, _sync_directory, _write
from genosapiens.wgs_science_coordinator import ScienceCoordinator
from genosapiens.wgs_science_runtime import ScienceClock, ScienceDeadlineExpired, bounded_operation


@dataclass(frozen=True, slots=True)
class ScienceExecutionResult:
    status: dict
    evidence: object | None
    case_values: dict | None
    scope_name: str | None
    scope_sha256: str | None


@contextmanager
def _context_phases(timings):
    """Measure disjoint attempt phases, including the phase that fails.

    Retry attempts accumulate in the same bounded sample. These values are
    diagnostics only and never participate in scientific evidence identities.
    Failure cleanup remains charged to the operation that failed.
    """
    phase, started = 'checkpoint_setup', time.monotonic()
    def advance(name):
        nonlocal phase, started
        now = time.monotonic()
        if phase is not None:
            timings[phase] = timings.get(phase, 0.0) + now - started
        phase, started = name, now
    try:
        yield advance
    finally:
        advance(None)


def _status(root, coordinator, publication):
    value = dict(coordinator.status(), publication=publication)
    value['scientific_status'] = value['status']
    if publication != 'complete':
        value['status'] = 'incomplete'
    temporary = root / '.science-status.pending'
    _write(temporary, canonical(value))
    os.replace(temporary, root / 'science-status.json')
    _sync_directory(root)
    return value


def produce_context(closed, expected, *, optional_filter, checkpoints, scratch_directory,
                    resource_factory, context_factory, reference_pin, observer_factory, timings):
    with _context_phases(timings) as phase:
        pending = checkpoints.begin(expected)
        observer = None if observer_factory is None else observer_factory(pending, closed.ordinal)
        with nullcontext() if observer is None else observer:
            with checkpoints.evidence_writer(pending, scratch_directory=scratch_directory) as writer, ExitStack() as stack:
                phase('resource_binding')
                resources, opened = resource_factory(writer, stack)
                phase('context_evaluation')
                case_name = 'checkpoint-case-' + str(closed.ordinal)
                case = stack.enter_context(CompleteScientificCase(writer, resources, name=case_name,
                    context_domain=ContextDomain.from_ordinals((closed.ordinal,))))
                with context_factory(writer, resources, opened, closed, optional_filter=optional_filter,
                        output_directory=pending, failure_observer=observer) as context:
                    case.add_context(context)
                    context_name, context_sha256 = context.name, context.sha256
                phase('context_publication')
                if observer is not None:
                    observer.stage('context_publication')
                case.publish_primary()
                case.publish(opened, reference_pin, {})
                case_sha256 = case.sha256
                stack.close()
                evidence = writer.finish()
                phase('writer_cleanup')
            phase('checkpoint_commit')
            return checkpoints.commit(pending, expected, evidence, context_name=context_name,
                context_sha256=context_sha256, case_name=case_name, case_sha256=case_sha256,
                scratch_directory=scratch_directory)


def execute_context_science(contexts, *, coordinator, resource_factory, context_factory, output_directory,
                            reference_pin, eligibility_root=None, profile_root=None,
                            authenticated_resources=None, progress=None, failure_observer_factory=None,
                            resource_preparation_seconds=0, workers=1, worker_initialize=None, publication_readers=1,
                            publication_worker_initialize=None, resource_phase_seconds=None):
    """Keep each scientific writer independent, then rank committed contexts."""
    _require(type(coordinator) is ScienceCoordinator)
    root = _path(output_directory)
    root.mkdir(mode=0o700, parents=True, exist_ok=True)
    _require(root.is_dir() and not any(root.iterdir()))
    # Retain the compact scope before the context stream closes its source index.
    # Its policy and complete membership digests are never scientific scores.
    scope_path = root / 'science-scope.jsonl'
    scope_digest, scope_count = hashlib.sha256(), 0
    fd = os.open(scope_path, os.O_WRONLY | os.O_CREAT | os.O_EXCL | os.O_NOFOLLOW, 0o600)
    with os.fdopen(fd, 'wb') as stream:
        for row in coordinator.roster.rows():
            raw = canonical(row)
            _require(len(raw) <= 16 * 1024)
            scope_digest.update(raw)
            scope_count += 1
            stream.write(raw)
        stream.flush()
        os.fsync(stream.fileno())
    _require(scope_count == coordinator.roster.counts['original'])
    status = _status(root, coordinator, 'not_started')
    from genosapiens.wgs_run_manifest import atomic_current
    started = time.monotonic()
    samples = []
    aggregation_progress = None
    worker_metrics = None
    context_elapsed_seconds = None
    publication_phases = dict(input_samples=[], merge_samples=[], readers=publication_readers)
    def performance(publication_seconds=None):
        atomic_current(root / 'science-performance.json', dict(
            schema='genosapiens.science-performance/v1', counts=coordinator.counts(),
            elapsed_seconds=time.monotonic() - started, context_samples=samples,
            publication_seconds=publication_seconds, aggregation=aggregation_progress, publication_phases=publication_phases,
            selected_contexts=coordinator.roster.counts['selected'],
            protected_overflow=coordinator.roster.counts['overflow'],
            seconds_remaining=None if type(coordinator.budget) is ScienceClock else coordinator.budget.remaining(),
            report_seconds_remaining=None if type(coordinator.budget) is ScienceClock else coordinator.budget.remaining(reporting=True),
            workers=workers, worker_metrics=worker_metrics, scope_sha256=coordinator.roster.sha256,
            context_elapsed_seconds=context_elapsed_seconds,
            resource_preparation_seconds=resource_preparation_seconds,
            resource_phase_seconds=resource_phase_seconds,
            resource_index='shared_read_only_prefix_with_private_additions',
            verification=dict(fresh_contexts='live_producer_reuse',
                omitted_fresh_replays=coordinator.checkpoints.omitted_fresh_replays,
                restored_contexts='independent_authentication',
                fresh_artifact_checks='omitted',
                omitted_checks=['hpo_typed_store_checksum_replay',
                    'admission_resource_reconstruction', 'phenotype_proof_model_validation',
                    'gene_resource_model_validation', 'disease_subset_model_validation',
                    'resource_row_recursive_validation', 'resource_row_canonical_replay',
                    'full_shard_rehash', 'node_digest_recheck',
                    'canonical_row_reencoding', 'repeated_file_signatures', 'dag_coverage_replay']),
            limits='First 32 attempts. Context size and source costs can differ.'))

    if worker_initialize is not None:
        from genosapiens.wgs_science_parallel import run_parallel
        def evaluate_worker(closed, identity, *, optional_filter, timings):
            return produce_context(closed, identity, optional_filter=optional_filter,
                checkpoints=coordinator.checkpoints, scratch_directory=coordinator.scratch_directory,
                resource_factory=resource_factory, context_factory=context_factory, reference_pin=reference_pin,
                observer_factory=failure_observer_factory, timings=timings)
        def updated(result):
            if result is not None and len(samples) < 32:
                samples.append(result)
            performance()
            _status(root, coordinator, 'not_started')
            from genosapiens.wgs_runtime import analysis_progress
            counts = coordinator.counts()
            analysis_progress(progress, 'scientific_evaluation', 'context_publication',
                counts={'contexts_completed': counts['checkpointed'], 'contexts_total': counts['selected']},
                unit='contexts', completed=counts['checkpointed'], total=counts['selected'], scope=0)
        worker_metrics = run_parallel(contexts, coordinator=coordinator, initialize=worker_initialize,
            evaluate=evaluate_worker, workers=workers, root=root, updated=updated)
    else:
        source = iter(contexts)
        while True:
            try:
                with bounded_operation(coordinator.budget, maximum_seconds=coordinator.context_timeout_seconds):
                    closed = next(source)
            except StopIteration:
                break
            except ScienceDeadlineExpired:
                coordinator.stop_reason = ('scientific_deadline' if coordinator.budget.remaining() <= 0
                    else 'unknown_failure_scope')
                break
            ordinal = closed.ordinal
            identity = coordinator.identity(ordinal)
            _require(identity.context_id == closed.context_id)
            phase_seconds = {}
            def evaluate(expected, *, optional_filter):
                return produce_context(closed, expected, optional_filter=optional_filter,
                    checkpoints=coordinator.checkpoints, scratch_directory=coordinator.scratch_directory,
                    resource_factory=resource_factory, context_factory=context_factory, reference_pin=reference_pin,
                    observer_factory=failure_observer_factory, timings=phase_seconds)
            context_started = time.monotonic()
            attempts_before = coordinator.attempts
            try:
                coordinator.run_one(ordinal, evaluate)
            except Exception:
                # The coordinator records unknown failures and stops scheduling.
                # Completed artifacts remain separately authenticated. A fresh
                # aggregation owner must still authenticate its resources below.
                if coordinator.stop_reason is None:
                    raise
            if len(samples) < 32 and coordinator.attempts != attempts_before:
                context_seconds = time.monotonic() - context_started
                # Includes the independent durable reload after the callback,
                # initial reuse checks, and coordinator bookkeeping.
                phase_seconds['coordinator_checks'] = max(0.0, context_seconds - sum(phase_seconds.values()))
                samples.append(dict(seconds=context_seconds,
                    member_count=closed.member_count, attempts=coordinator.attempts - attempts_before,
                    state=coordinator._state(ordinal), phase_seconds=phase_seconds))
            performance()
            status = _status(root, coordinator, 'not_started')
            from genosapiens.wgs_runtime import analysis_progress
            counts = coordinator.counts()
            analysis_progress(progress, 'scientific_evaluation', 'context_publication',
                counts={'contexts_completed': counts['checkpointed'], 'contexts_total': counts['selected']},
                unit='contexts', completed=counts['checkpointed'], total=counts['selected'], scope=0)
            if coordinator.stop_reason is not None:
                break

    context_elapsed_seconds = time.monotonic()-started
    values = scope_name = scope_sha256 = evidence = None
    publication_started = time.monotonic()
    try:
        remaining = coordinator.budget.remaining(reporting=True)
        if remaining <= 0:
            raise ScienceDeadlineExpired('reporting deadline expired')
        with bounded_operation(coordinator.budget, maximum_seconds=remaining, reporting=True):
            with EvidenceWriter(root, scratch_directory=coordinator.scratch_directory,
                    _skip_verification=True) as writer, ExitStack() as stack:
                resources, opened = resource_factory(writer, stack)
                case = stack.enter_context(CompleteScientificCase(writer, resources, name='case',
                    context_domain=coordinator.completed_domain()))
                aggregated_contexts = 0
                aggregation_started = time.monotonic()
                aggregation_samples = []
                def checkpoint_for(ordinal):
                    identity = coordinator.identity(ordinal)
                    checkpoint = coordinator.checkpoints._read(coordinator.checkpoints._destination(identity), identity)
                    recorded = coordinator.connection.execute('SELECT checkpoint_sha256 FROM science_context_directory WHERE ordinal=?',
                        (ordinal,)).fetchone()
                    _require(recorded == (checkpoint.completion_sha256,) and checkpoint.case_name is not None)
                    return checkpoint
                if publication_worker_initialize is not None:
                    from genosapiens.wgs_science_publication import prepared_inputs, prepare_projection, append_projection
                    # Compact identities are copied before fork. Workers never use
                    # the coordinator's mutable completion database.
                    identities = {ordinal: coordinator.identity(ordinal) for ordinal in case.context_domain}
                    def initialize_inputs(local):
                        live_resources = publication_worker_initialize(local)
                        def prepare(job):
                            before = time.monotonic()
                            identity = identities[job['ordinal']]
                            checkpoint = coordinator.checkpoints._read(coordinator.checkpoints._destination(identity), identity)
                            _require(checkpoint.completion_sha256 == job['completion_sha256'])
                            with coordinator.checkpoints.open_evidence(checkpoint,
                                    scratch_directory=coordinator.scratch_directory) as reader:
                                opened_at = time.monotonic()
                                with CompletedScientificCaseResources(reader, reader.reference(resources.name, kind='scientific-case-resources-v1'),
                                        expected_parent_bindings=resources.parent_bindings, eligibility_root=eligibility_root,
                                        profile_root=profile_root, authenticated_resources=authenticated_resources,
                                        _live_resources=live_resources if reader._live_prefix_source is not None else None) as saved_resources:
                                    with CompletedScientificCase(reader, reader.reference(checkpoint.case_name, kind='complete-scientific-case-v1'),
                                            saved_resources) as saved:
                                        decoded_at = time.monotonic()
                                        result = prepare_projection(reader, saved, job['path'])
                            result.update(artifact_read_seconds=opened_at-before, projection_decode_seconds=decoded_at-opened_at,
                                input_seconds=time.monotonic()-before)
                            return result
                        return prepare
                    inputs = stack.enter_context(__import__('contextlib').closing(prepared_inputs(case.context_domain,
                        workers=publication_readers, initialize=initialize_inputs, checkpoint_for=checkpoint_for,
                        root=root, timeout=coordinator.context_timeout_seconds)))
                else:
                    inputs = ((ordinal, None, None, None) for ordinal in case.context_domain)
                for ordinal, prepared_path, prepared_metadata, input_metrics in inputs:
                    aggregate_started = time.monotonic()
                    checkpoint = checkpoint_for(ordinal)
                    with coordinator.checkpoints.open_evidence(checkpoint,
                            scratch_directory=coordinator.scratch_directory) as reader:
                        if prepared_path is not None:
                            merged = append_projection(case, reader, checkpoint, prepared_path, prepared_metadata)
                            publication_phases['worker_metrics'] = input_metrics
                            if len(publication_phases['input_samples']) < 32:
                                publication_phases['input_samples'].append(prepared_metadata)
                                publication_phases['merge_samples'].append(merged)
                        else:
                            with CompletedScientificCaseResources(reader, reader.reference(resources.name, kind='scientific-case-resources-v1'),
                                    expected_parent_bindings=resources.parent_bindings, eligibility_root=eligibility_root,
                                    profile_root=profile_root, authenticated_resources=authenticated_resources,
                                    _live_resources=(resources if reader._live_prefix_source is not None
                                        and writer._copied_prefix_source is reader._live_prefix_source else None)) as saved_resources:
                                with CompletedScientificCase(reader, reader.reference(checkpoint.case_name, kind='complete-scientific-case-v1'),
                                        saved_resources, progress=progress) as saved:
                                    append_completed_context(case, reader, saved)
                    aggregated_contexts += 1
                    aggregation_samples.append(time.monotonic() - aggregate_started)
                    del aggregation_samples[:-8]
                    remaining_contexts = len(case.context_domain) - aggregated_contexts
                    aggregation_progress = dict(contexts_completed=aggregated_contexts,
                        contexts_total=len(case.context_domain), elapsed_seconds=time.monotonic() - aggregation_started,
                        estimate_scope='context_aggregation_only',
                        estimate_quality='complete' if not remaining_contexts else
                            'heuristic' if len(aggregation_samples) >= 2 else 'insufficient_samples',
                        remaining_seconds_range=([0.0, 0.0] if not remaining_contexts else
                            [remaining_contexts * min(aggregation_samples), remaining_contexts * max(aggregation_samples)]
                            if len(aggregation_samples) >= 2 else None),
                        report_generation_included=False)
                    performance()
                    from genosapiens.wgs_runtime import analysis_progress
                    analysis_progress(progress, 'publication', 'saved_scientific_projections',
                        counts={'selected_contexts_replayed': aggregated_contexts},
                        unit='contexts', completed=aggregated_contexts, total=len(case.context_domain), scope=0)
                ranking_started = time.monotonic()
                case.publish_primary()
                publication_phases['ranking_seconds'] = time.monotonic()-ranking_started
                selection_started = time.monotonic()
                case.publish(opened, reference_pin, {})
                publication_phases['selection_seconds'] = time.monotonic()-selection_started
                scope_name = 'case-scope'
                status = coordinator.status()
                directory_sha = writer.add_sequence(scope_name + '.contexts', 'scientific-context-directory-v1', coordinator.rows())
                def scope_rows():
                    digest, count = hashlib.sha256(), 0
                    with scope_path.open('rb') as stream:
                        for raw in iter(lambda: stream.readline(16 * 1024 + 1), b''):
                            _require(len(raw) <= 16 * 1024)
                            digest.update(raw)
                            count += 1
                            yield json.loads(raw)
                    _require(digest.hexdigest() == scope_digest.hexdigest() and count == scope_count)
                roster_sha = writer.add_sequence(scope_name + '.roster', 'scientific-scope-roster-v1', scope_rows())
                scope_sha256 = writer.add(scope_name, 'scientific-scope-v1', dict(status,
                    policy=json.loads(coordinator.roster.policy.as_bytes()), budget=json.loads(coordinator.budget.as_bytes()),
                    computation_sha256=coordinator.computation_sha256, resources=resources.reference,
                    protected_overflow=coordinator.roster.counts['overflow'],
                    roster=dict(name=scope_name + '.roster', sha256=roster_sha),
                    contexts=dict(name=scope_name + '.contexts', sha256=directory_sha)))
                values = dict(case_name=case.name, case_sha256=case.sha256, resources_name=resources.name,
                    resources_sha256=resources.sha256, parent_bindings=tuple(sorted(resources.parent_bindings.items())),
                    context_count=case.context_count, positive_context_count=case.positive_context_count,
                    member_count=case.member_count, pair_count=case.pair_count, hypothesis_count=case.hypothesis_count,
                    policy_id=resources.fusion_policy.policy_id, policy_sha256=resources.fusion_policy_sha256,
                    derivation_sha256=resources.fusion_policy.derivation_sha256)
                stack.close()
                evidence = writer.finish()
        status = _status(root, coordinator, 'complete')
    except ScienceDeadlineExpired:
        # This status is small and independent of unfinished report evidence.
        status = _status(root, coordinator, 'incomplete_deadline')
        evidence = values = scope_name = scope_sha256 = None
    except Exception:
        _status(root, coordinator, 'incomplete_failure')
        raise
    performance(time.monotonic() - publication_started)
    return ScienceExecutionResult(status, evidence, values, scope_name, scope_sha256)
