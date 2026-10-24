"""Executable public-resource product builder for bounded WGS analysis."""

from __future__ import annotations

from dataclasses import dataclass, field, replace
from contextlib import closing, ExitStack, contextmanager
import hashlib
import json
import os
import re
from itertools import chain
from pathlib import Path
import sqlite3
import stat
from typing import Iterable, Iterator

from genosapiens import population_policy as frequency
from genosapiens.bounded_wgs_config import (
    BoundedWgsExecutionConfig, bounded_wgs_execution_config_as_json,
)
from genosapiens.candidate_admission import (
    CandidateAdmissionError,
    RepeatCatalogInput,
    ShardedAdmissionSources,
    _bytes as admission_bytes,
    _generate_projected_candidate_admission,
    prepare_sharded_candidate_admission,
)
from genosapiens.common_allele_public_indexes import VerifiedPublicIndexes
from genosapiens.evidence_model import (
    DiseaseGeneEvidence,
    DiseaseProfile,
    scoring_view,
    _typed,
)
from genosapiens.phenotype_ingestion import (
    PhenotypeIngestionResult, phenotype_ingestion_as_json,
)
from genosapiens.genotype_hypotheses import hypothesis_policy_bytes
from genosapiens.population_annotation import (
    PopulationQueryLimits,
    PopulationResourceSession,
)
from genosapiens.track1_preflight import SelectedTechnicalAssessment
from genosapiens.track1_publication import FinalArtifact, FinalRole
from genosapiens.track1_submission import _open_submission_reference
from genosapiens.variant_inventory import IncrementalVariantClassSummary
from genosapiens.vcf_qc import VcfQcIncrementalSummary
from genosapiens.wgs_admission import (
    LedgerAdmissionInput,
    NativeSourceContract,
    iter_ledger_admission_inputs,
    replay_sequence_annotation_admission,
)
from genosapiens.wgs_admission_projection import (
    AdmissionProjectionError, AdmissionProjectionReader, AdmissionProjectionWriter, AuthenticatedAdmissionOwner,
    EncodedAdmissionRow, PendingNativeAdmission, ValidatedSequenceAdmission, _AdmissionProjectionPool, admission_worker_setup,
    _produce_native_admission_record,
)
from genosapiens.wgs_allele_query_cache import (
    CaseAlleleQueryCache, adapter_identity as query_adapter_identity,
    population_source_bindings, clinical_source_bindings,
)
from genosapiens.wgs_analysis import (
    BoundedAnalysisProducts,
    PreparedAnalysisContext,
    _products_from_checkpoint_content,
)
from genosapiens.wgs_shards import read_bounded_metadata
from genosapiens.wgs_native_annotation import NativeAnnotationSession, NativeAnnotationError
from genosapiens.wgs_annotation import (
    annotation_coverage_summary,
    AnnotationStageResult,
    BoundedAnnotationSource,
    VepTransportPlan,
    iter_annotation_stage_encoded_records,
    open_retained_projection,
    iter_retained_alleles,
    run_annotation_partitions,
    write_vep_partitions,
)
from genosapiens.wgs_contexts import (
    WgsContextError,
    iter_closed_contexts,
)
from genosapiens.wgs_failures import (
    WgsFailure,
    failure_boundary,
    managed_resource_session,
)
from genosapiens.wgs_clinical_session import WgsClinicalResourceSession
from genosapiens.wgs_ledger import canonical, read_artifact
from genosapiens.wgs_preparation import WgsNormalizationSummary
from genosapiens.wgs_reporting import build_bounded_reports, profile_source_accounting
from genosapiens.wgs_runtime import runtime_phase
from genosapiens.wgs_scientific_artifact import write_scientific_artifact, ScientificArtifactReader
from genosapiens.wgs_scientific_products import scientific_products
from genosapiens.wgs_phenotype_evidence import _hash as scientific_hash
from genosapiens.wgs_sensitivity import _checked_ic_hash

from genosapiens.wgs_validation_audit import ValidationAuditWriter
from genosapiens.track1_public_validation import (
    HPO_MANIFEST_SHA256,
    RARE_DISEASE_BUNDLE_SHA256,
    REPEAT_CONTRACT_PATH,
    REPEAT_MANIFEST_SHA256,
    RUNTIME_IMAGE,
    VEP_BUNDLE_SHA256,
    VEP_CONTRACT_PATH,
    PublicPolicyFamily,
    PublicPolicyRuntime,
    PublicResourceContext,
    _population_template,
    _subset_gene_store,
    _subset_profile_store,
)


class WgsPublicAnalysisError(RuntimeError):
    """Static failure while building bounded public scientific products."""


def _phase_call(context, phase, function, /, *args, **kwargs):
    with runtime_phase(context.phase_recorder, phase):
        return function(*args, **kwargs)


_MODELS = (
    "compound_heterozygous",
    "dominant",
    "homozygous_recessive",
    "incomplete",
    "mitochondrial",
    "mosaic",
    "repeat",
    "source-requirement-context-only",
    "structural",
    "x_linked",
    "y_linked",
)


def _private_directory(path: Path, *, resume: bool = False) -> Path:
    if path.exists() or path.is_symlink():
        metadata = path.lstat()
        if (not resume or path.is_symlink()
                or not stat.S_ISDIR(metadata.st_mode)
                or stat.S_IMODE(metadata.st_mode) != 0o700
                or metadata.st_uid != os.getuid()):
            raise WgsPublicAnalysisError("bounded product directory already exists")
        return path
    path.mkdir(mode=0o700, parents=True)
    return path


def _read_phenotype(context: PreparedAnalysisContext) -> PhenotypeIngestionResult:
    try:
        content = read_bounded_metadata(context.phenotype_ingestion_path)
        if hashlib.sha256(content).hexdigest() != context.phenotype_ingestion_sha256:
            raise ValueError
        phenotype = _typed(json.loads(content), PhenotypeIngestionResult, decode=True)
        if phenotype_ingestion_as_json(phenotype).encode("ascii") != content:
            raise ValueError
        return phenotype
    except MemoryError:
        raise
    except Exception:
        raise WgsPublicAnalysisError(
            "prepared phenotype ingestion could not be replayed"
        ) from None


def _read_source_summaries(context: PreparedAnalysisContext):
    try:
        values = []
        for artifact, value_type in zip(
            context.source_manifest.summary.reports,
            (
                VcfQcIncrementalSummary,
                WgsNormalizationSummary,
                IncrementalVariantClassSummary,
            ),
            strict=True,
        ):
            values.append(_typed(
                json.loads(read_artifact(context.source_root, artifact)),
                value_type,
                decode=True,
            ))
        return tuple(values)
    except MemoryError:
        raise
    except Exception:
        raise WgsPublicAnalysisError(
            "prepared source summaries could not be replayed"
        ) from None


def _pinned_annotation_coverage(fasta_path, contract):
    from genosapiens.wgs_annotation_coverage import load_reference_coverage
    return load_reference_coverage(fasta_path,
        expected_sha256=contract['fasta']['fai_sha256'],
        expected_bytes=contract['fasta']['fai_size_bytes'])


def _resource_session_sha256(
    resources: PublicResourceContext,
    *,
    executable_stack_sha256: str,
    runtime_image_digest: str,
    vep_contract_bytes: bytes,
) -> str:
    from genosapiens.wgs_annotation_prefix import annotation_resource_session_sha256
    return annotation_resource_session_sha256(
        reference_manifest_sha256=resources.submission_reference.manifest_sha256,
        vep_bundle_manifest_sha256=VEP_BUNDLE_SHA256,
        vep_contract_sha256=hashlib.sha256(vep_contract_bytes).hexdigest(),
        executable_stack_sha256=executable_stack_sha256,runtime_image_digest=runtime_image_digest)



def _admission_inputs(context, *, scratch_name):
    if context.read_intake is not None:
        yield from context.read_intake.admissions(context.execution_config)
    else:
        yield from iter_ledger_admission_inputs(
            source_root=context.source_root, source_manifest=context.source_manifest,
            working_root=context.working_root, working_manifest=context.working_manifest,
            source_routing_context=context.source_routing_context,
            scratch_path=context.work_root / scratch_name, execution_config=context.execution_config,
            prepared_index=context.prepared_index,
            prepared_state_root=context.prepared.root.parent / "state",
            retained_projection_path=context.retained_projection_path)


def _intake_binding(context):
    return None if context.read_intake is None else context.read_intake.manifest.annotation_binding


def _admission_projection_jobs(context, *, plan_root, plan, annotation_root,
                               annotation, resource_session_sha256, setup):
    """Authenticate saved shards and route bounded immutable entries to workers."""
    records = iter_annotation_stage_encoded_records(annotation,plan_root=plan_root,plan=plan,
        source_manifest=context.source_manifest,working_manifest=context.working_manifest,
        output_root=annotation_root,resource_session_sha256=resource_session_sha256,
        bundle_manifest_sha256=VEP_BUNDLE_SHA256,execution_config=context.execution_config,
        intake_binding=_intake_binding(context),prepared_index=context.prepared_index,
        prepared_state_root=context.prepared.root.parent / "state")
    range_mode = json.loads(setup)['projection'] is not None
    def keys():
        if range_mode:
            with closing(open_retained_projection(context.retained_projection_path,
                    prepared_index=context.prepared_index,routing=context.source_routing_context,
                    execution_config=context.execution_config,
                    expected_identities=context.working_manifest.summary.retained_unique_identities,
                    expected_occurrences=context.working_manifest.summary.retained_occurrences)) as db:
                # Planning sends every retained identity with sequence geometry
                # to VEP, including native mosaic/mitochondrial identities.
                # Preparation route alone does not identify that transport.
                with closing(db.execute('SELECT i.ordinal,i.identity_id,f.contig '
                        'FROM identities i JOIN facts f USING(identity_id) ORDER BY i.ordinal')) as cursor:
                    for ordinal,identity,contig in cursor:
                        yield ordinal,identity,contig is not None,None
        else:
            # LongTR owns complete source-locus assembly; read rescue owns its
            # authenticated merged intake. Both feed this same worker engine.
            for value in _admission_inputs(context,scratch_name='projection-admission.sqlite'):
                yield value.retained_ordinal,value.identity_id,value.transport == 'sequence_vep',value.as_bytes()
    previous = 0
    with closing(records), closing(keys()) as inputs:
        for ordinal,identity,sequence,content in inputs:
            if ordinal != previous+1:
                raise WgsPublicAnalysisError('admission input replay order changed')
            previous = ordinal
            saved = next(records,None) if sequence else None
            if sequence and saved is None:
                raise WgsPublicAnalysisError('annotation response coverage is incomplete')
            yield ordinal,identity,saved,content
        if previous != plan.total_retained_identities or next(records,None) is not None:
            raise WgsPublicAnalysisError('annotation response coverage changed')


def _admission_projection_binding(context: PreparedAnalysisContext,
                                  plan: VepTransportPlan,
                                  annotation: AnnotationStageResult,
                                  resource_session_sha256: str, *, consumer_executable_sha256: str,
                                  completed_annotation_prefix_sha256: str | None, screening) -> bytes:
    return canonical({
        "schema": "genosapiens.wgs-enriched-admission-projection", "version": 5,
        "screening_context_sha256": screening.sha256,
        "eligibility_policy_sha256": screening.policy.sha256,
        "analysis_method_sha256": screening.policy.method_sha256,
        "consumer_executable_sha256": consumer_executable_sha256,
        "completed_analysis_input_sha256": None if context.completed_analysis_input is None
            else context.completed_analysis_input.descriptor.sha256,
        "completed_annotation_prefix_sha256": completed_annotation_prefix_sha256,
        "codec": "canonical-input-v1/proofs-v2/facts-v1/zlib-1/layout-v5",
        "identities": plan.total_retained_identities,
        "occurrences": context.working_manifest.summary.retained_occurrences if context.read_intake is None
            else context.read_intake.manifest.admitted_source_occurrences,
        "sequence_requests": plan.sequence_requests,
        "prepared_index_sha256": context.prepared_index.sha256,
        "source_manifest_sha256": context.source_manifest.sha256,
        "working_manifest_sha256": context.working_manifest.sha256,
        "source_routing_sha256": context.source_routing_context.contract.sha256,
        "plan_sha256": plan.sha256,
        "annotation_stage_sha256": annotation.manifest.sha256,
        "resource_session_sha256": resource_session_sha256,
        "bundle_manifest_sha256": VEP_BUNDLE_SHA256,
        "execution_limits_sha256": hashlib.sha256(
            bounded_wgs_execution_config_as_json(context.execution_config)).hexdigest(),
        "read_intake_sha256": None if context.read_intake is None else
            context.read_intake.manifest.sha256,
        "read_intake_source_bindings_sha256": None if context.read_intake is None else
            hashlib.sha256(admission_bytes(context.read_intake.manifest.source_bindings)).hexdigest(),
    })


@dataclass(frozen=True, slots=True)
class PublishedAdmissionProjection:
    binding: bytes
    completion: object
    eligibility_root: Path
    eligibility_manifest: object


def _admission_stage_binding(context, binding, *, runtime_image_digest):
    from genosapiens.wgs_ledger import StageRunBinding, ParentIdentity
    parents=json.loads(binding); digest=hashlib.sha256(binding).hexdigest()
    stack=parents['consumer_executable_sha256']
    limits=parents['execution_limits_sha256']
    return StageRunBinding('admission_projection',digest,stack,runtime_image_digest,limits,
        context.prepared.manifest.preparation_run_options.sha256,
        context.source_manifest.source.input_sha256,
        context.source_routing_context.contract.reference_manifest_sha256,
        parents['eligibility_policy_sha256'],(parents['resource_session_sha256'],),(
            ParentIdentity('stage_parent',parents['schema'],parents['version'],digest),
            ParentIdentity('partition_plan','genosapiens.wgs-vep-transport-plan',1,parents['plan_sha256']),
            ParentIdentity('execution_limits',context.execution_config.schema,context.execution_config.version,limits),
            ParentIdentity('executable_stack','genosapiens.executable-stack',1,stack),
            ParentIdentity('runtime_image','genosapiens.runtime-image',1,runtime_image_digest.removeprefix('sha256:'))))


def _publish_admission_projection(path: Path, *, screening_factory, **arguments):
    from genosapiens.wgs_admission_screening import AdmissionScreeningThread
    from genosapiens.wgs_query_limits import _cancellable_queries
    from genosapiens.wgs_admission_diagnostics import Owner, current, preserving_close
    owner = current() or Owner(path.parent/'admission-diagnostic.json')
    owned = current() is None
    try:
        with owner.activate():
            screening_thread=AdmissionScreeningThread(screening_factory)
            with preserving_close(screening_thread, 'screening_close'):
                with _cancellable_queries(screening_thread._cancelled):
                    try:
                        result = _publish_admission_projection_owned(path,screening_thread=screening_thread,**arguments)
                        owner.value['state'] = 'complete'
                        owner.snapshot(force=True)
                        return result
                    except BaseException as error:
                        owner.capture(error)
                        screening_thread.fail(error)
                        raise
    except BaseException as error:
        owner.capture(error)
        if owned:
            owner.terminal()
        raise


def _publish_admission_projection_owned(path: Path, *, context: PreparedAnalysisContext,
                                  plan_root: Path, plan: VepTransportPlan,
                                  annotation_root: Path, annotation: AnnotationStageResult,
                                  resource_session_sha256: str, consumer_executable_sha256: str,
                                  completed_annotation_prefix_sha256: str | None,
                                  screening_thread, eligibility_dependencies, runtime_image_digest,
                                  native_session=None, profile_intake=None) -> PublishedAdmissionProjection:
    """Validate and screen fresh annotation inputs, then journal the complete seal."""
    from genosapiens.wgs_analysis_eligibility import AnalysisEligibilityWriter, iter_eligibility_decisions
    from genosapiens.wgs_admission_projection import (publish_admission_projection_completion,
        read_admission_projection_completion)
    from genosapiens.wgs_ledger import StageCommitState, stage_lock, publish_bytes
    import time
    from genosapiens.wgs_admission_diagnostics import boundary, current, preserving_close
    diagnostic = current()
    screening=screening_thread.setup.screening
    binding=_admission_projection_binding(context,plan,annotation,resource_session_sha256,
        consumer_executable_sha256=consumer_executable_sha256,
        completed_annotation_prefix_sha256=completed_annotation_prefix_sha256,screening=screening)
    diagnostic.bind(binding)
    parents=json.loads(binding)
    stage_root=path.parent/'admission-stage'
    eligibility_root=path.parent.parent/'eligibility'
    # An earlier analysis failure can leave admission entirely unstarted. Use
    # this stage's artifacts to distinguish its first attempt from a replay.
    resuming_admission=context.resume and any(p.exists() or p.is_symlink()
        for p in (path,stage_root,eligibility_root))
    if profile_intake is not None:
        from genosapiens.wgs_profile_projection import FreshAdmissionProfileIntake
        if type(profile_intake) is not FreshAdmissionProfileIntake or resuming_admission:
            raise WgsPublicAnalysisError('fresh profile intake owner is required')
        profile_intake.bind(binding,screening)
    stage_binding=_admission_stage_binding(context,binding,runtime_image_digest=runtime_image_digest)
    stage_root.mkdir(mode=0o700,exist_ok=resuming_admission)
    state_path=stage_root/'state.sqlite'
    with stage_lock(stage_root,'admission_projection') as lock, StageCommitState(
            stage_root,state_path,stage_binding,lock=lock,create=not resuming_admission) as state:
        if resuming_admission:
            completion=read_admission_projection_completion(state)
            with closing(AdmissionProjectionReader(path,expected_binding=binding,
                    expected_completion_sha256=completion.sha256,execution_config=context.execution_config)) as reader:
                manifest=reader.eligibility_manifest
            return PublishedAdmissionProjection(binding,completion,eligibility_root,manifest)
        start=state.begin(partition_id='admission',scientific=False,max_attempts=1)
        started=time.monotonic(); last=0
        native_metrics=dict(completed_inputs=0,completed_occurrences=0,wall_nanoseconds=0)
        progress_path=path.parent/'admission-progress.json'
        def progress(phase,writer,decisions,pool,*,force=False):
            nonlocal last
            now=time.monotonic()
            if not force and now-last<1: return
            screen=screening_thread.snapshot()
            diagnostic.snapshot(counters=dict(screened=screen['screened_windows'],
                pending=screen['screened_windows']-screen['committed_windows'],
                committed=screen['committed_windows'],cancelled=screen['cancelled_windows'],
                sealed=screen['committed_windows'] if writer.finished and decisions.finished else 0),force=force)
            record=dict(schema='genosapiens.admission-progress',version=2,phase=phase,
                binding_sha256=hashlib.sha256(binding).hexdigest(),elapsed_seconds=now-started,
                total_inputs=parents['identities'],total_occurrences=parents['occurrences'],
                validated_inputs=pool.metrics['received_rows']-pool.metrics['received_native_inputs']+native_metrics['completed_inputs'],
                validated_occurrences=pool.metrics['received_occurrences']-pool.metrics['received_native_occurrences']+native_metrics['completed_occurrences'],
                acquired_inputs=pool.metrics['acquired_rows'],acquired_occurrences=pool.metrics['acquired_occurrences'],
                screened_inputs=screen['screened_inputs'],screened_occurrences=screen['screened_occurrences'],
                written_inputs=writer.count,written_occurrences=writer.occurrences,
                eligibility_written_inputs=decisions.rows,eligibility_written_occurrences=decisions.occurrences,
                committed_inputs=screen['committed_inputs'],committed_occurrences=screen['committed_occurrences'],
                input_dispositions=decisions.counts,occurrence_dispositions=decisions.occurrence_counts,
                public_query_batches=screen['public_query_batches'],screening=screen,native=native_metrics,
                profile_intake=None if profile_intake is None else dict(profile_intake.metrics),
                admission_writer={**writer.metrics,'windows':writer.windows,'sealed':writer.finished},
                eligibility_writer={**decisions.metrics,'windows':decisions.windows,'sealed':decisions.finished,
                    'decision_shards':dict(decisions.writer.metrics),'witness_shards':dict(decisions.witness_writer.metrics)},
                worker_pids=pool.worker_pids,metrics=pool.metrics)
            temporary=progress_path.with_suffix('.tmp')
            fd=os.open(temporary,os.O_WRONLY|os.O_CREAT|os.O_TRUNC|os.O_NOFOLLOW,0o600)
            with os.fdopen(fd,'wb') as out: out.write(canonical(record))
            os.replace(temporary,progress_path);last=now
        try:
            context_artifact=publish_bytes(stage_root,'screening-context.json',screening.as_bytes())
            dependencies=tuple(sorted((*eligibility_dependencies,('screening_context',
                str(stage_root/context_artifact.name),context_artifact.sha256,context_artifact.size_bytes))))
            setup=admission_worker_setup(context,binding,screening=screening,
                defer_sequence_proofs=profile_intake is not None)
            with preserving_close(AdmissionProjectionWriter(path,binding=binding,execution_config=context.execution_config), 'writer_close') as writer, \
                    preserving_close(AnalysisEligibilityWriter(eligibility_root,execution_config=context.execution_config,
                        admission_binding_sha256=hashlib.sha256(binding).hexdigest(),
                        eligibility_policy_sha256=screening.policy.sha256,
                        source_manifest_sha256=screening_thread.setup.source_sha256,
                        reference_manifest_sha256=context.source_routing_context.contract.reference_manifest_sha256,
                        expected_inputs=parents['identities'],expected_occurrences=parents['occurrences'],
                        dependencies=dependencies,resource_bindings=screening_thread.setup.resource_bindings), 'writer_close') as decisions:
                with _AdmissionProjectionPool(context.execution_config,setup=setup, total_inputs=parents['identities']) as pool:
                    pool.abort_check=screening_thread.check
                    progress('validation_and_screening',writer,decisions,pool,force=True)
                    in_screening={}
                    def consume():
                        capacity=min(64,context.execution_config.public_queries.batch_rows,
                            context.execution_config.sqlite.cursor_batch_rows)
                        while len(in_screening)<2 and pool.next_acquire in pool.outstanding:
                            screening_thread.check()
                            with boundary('admission_result_collect', blocking=True):
                                result=pool.take_window(max_rows=capacity,
                                    yield_when=screening_thread.ready_to_write if in_screening else None)
                            if result is None:break
                            rows=[]
                            for row in result.rows:
                                screening_thread.check()
                                if type(row) is PendingNativeAdmission:
                                    if native_session is None:
                                        raise WgsPublicAnalysisError('native annotation session is missing')
                                    native_started=time.monotonic_ns()
                                    try:
                                        if profile_intake is None:
                                            row=_produce_native_admission_record(row,native_session,
                                                binding=binding,screening=screening)
                                        else:
                                            row=profile_intake.prepare_native(row,native_session)
                                    finally:native_metrics['wall_nanoseconds']+=time.monotonic_ns()-native_started
                                    native_metrics['completed_inputs']+=1
                                    native_metrics['completed_occurrences']+=(row.source_occurrences
                                        if profile_intake is None else row.row.source_occurrences)
                                rows.append(row)
                            intake_window=tuple(rows)
                            from genosapiens.wgs_profile_projection import _ValidatedNativeAdmission
                            window=tuple(row.row if type(row) in (ValidatedSequenceAdmission,_ValidatedNativeAdmission)
                                         else row for row in rows)
                            token=screening_thread.submit(window)
                            in_screening[token]=(result,intake_window,window)
                        with boundary('admission_screening_receive', blocking=True):
                            token,batch,witnesses=screening_thread.receive(check_wait=pool.check_running)
                            result,intake_window,window=in_screening[token]
                        if profile_intake is not None:
                            window,batch=profile_intake.finalize_window(intake_window,batch)
                        with boundary('admission_window_persist', blocking=True):
                            writer.add_screened_window(window,batch,eligibility_writer=decisions,witnesses=witnesses)
                        with boundary('admission_window_retire', blocking=True):
                            pool.consumed_window(result)
                            del in_screening[token]
                            screening_thread.committed(token)
                        progress('validation_and_screening',writer,decisions,pool)
                    pending=[]
                    try:
                        with closing(_admission_projection_jobs(context,plan_root=plan_root,plan=plan,
                                annotation_root=annotation_root,annotation=annotation,
                                resource_session_sha256=resource_session_sha256,setup=setup)) as entries:
                            for entry in entries:
                                screening_thread.check()
                                if pending and not pool.fits((*pending,entry)):
                                    if pool.full:consume()
                                    pool.submit(tuple(pending));pending=[]
                                if not pool.fits((entry,)):
                                    raise AdmissionProjectionError('single admission job exceeds capacity')
                                pending.append(entry)
                            if pending:
                                if pool.full:consume()
                                pool.submit(tuple(pending))
                        while pool.outstanding:consume()
                        with boundary('admission_complete', blocking=True):
                            observed=screening_thread.snapshot()
                            if (observed['committed_inputs']!=parents['identities']
                                    or observed['committed_occurrences']!=parents['occurrences']
                                    or pool.metrics['received_rows']!=parents['identities']
                                    or pool.metrics['received_occurrences']!=parents['occurrences']
                                    or pool.metrics['retired_rows']!=observed['committed_inputs']
                                    or pool.metrics['retired_occurrences']!=observed['committed_occurrences']
                                    or writer.metrics['committed_inputs']!=observed['committed_inputs']
                                    or decisions.metrics['committed_inputs']!=observed['committed_inputs']):
                                raise AdmissionProjectionError('admission progress closure differs')
                            screening_thread.finish(expected_inputs=parents['identities'],expected_occurrences=parents['occurrences'])
                            pool.finish()
                    except BaseException as error:
                        diagnostic.capture(error)
                        screening_thread.fail(error)
                        screening_thread.check()
                    finally:
                        # Stop queries and reap their groups while the worker
                        # pool's subreaper is still active, then stop workers.
                        with diagnostic.cleanup('screening_close'):
                            screening_thread.close()
                        secondary = getattr(screening_thread, 'cleanup_failure', None)
                        if secondary is not None:
                            if len(diagnostic.value['secondary']) < 4:
                                diagnostic.value['secondary'].append(secondary)
                            else:
                                diagnostic.value['secondary_truncated'] = True
                # No admission worker or screening session survives its seal.
                with boundary('admission_complete', blocking=True):
                    if profile_intake is not None:
                        profile_intake.finish_intake(expected_inputs=parents['identities'],
                            expected_occurrences=parents['occurrences'],input_dispositions=decisions.counts)
                    manifest=decisions.finish()
                    completion=writer.seal(eligibility_manifest=manifest)
                    publish_bytes(path.parent,'admission-workers.json',canonical(dict(
                        schema='genosapiens.wgs-admission-workers',version=2,
                        binding_sha256=hashlib.sha256(binding).hexdigest(),worker_pids=pool.worker_pids,
                        metrics=pool.metrics,validated_inputs=writer.count,screened_inputs=decisions.rows,
                        screening=screening_thread.snapshot(),native=native_metrics,
                        profile_intake=None if profile_intake is None else dict(
                            identity=profile_intake.identity,metrics=profile_intake.metrics),
                        admission_writer=writer.metrics,eligibility_writer=decisions.metrics,
                        input_dispositions=decisions.counts,complete=True)))
            with boundary('admission_complete', blocking=True):
                publish_admission_projection_completion(completion,state=state,start=start)
                progress('complete',writer,decisions,pool,force=True)
            return PublishedAdmissionProjection(binding,completion,eligibility_root,manifest)
        except BaseException as error:
            diagnostic.capture(error)
            with diagnostic.cleanup('progress_publish'):
                if 'writer' in locals() and 'decisions' in locals() and 'pool' in locals():
                    progress('failed',writer,decisions,pool,force=True)
            with diagnostic.cleanup('journal_publish'):
                if state.committed('admission') is None:
                    state.fail(start,failure_code='admission_projection_failed',
                        interrupted=isinstance(error,(KeyboardInterrupt,SystemExit)))
            raise


def _disease_subset(
    resources: PublicResourceContext,
    targets,
) -> tuple[DiseaseProfile | None, tuple[DiseaseGeneEvidence, ...]]:
    gene_store, diseases = _subset_gene_store(
        resources.disease_genes,
        frozenset(value.target for value in targets if value.target is not None
            and value.namespace is not None and not value.namespace.endswith(':locus')),
        frozenset(value.condition for value in targets if value.condition is not None),
    )
    profile_store = _subset_profile_store(
        resources.disease_profiles,
        diseases,
        resources.hpo_store,
    )
    genes = () if gene_store is None else (
        DiseaseGeneEvidence(
            "track1-public-gene-subset-" + gene_store.evidence_sha256[:16],
            RARE_DISEASE_BUNDLE_SHA256,
            gene_store,
        ),
    )
    profile = None if profile_store is None else DiseaseProfile(
        "track1-public-profile-subset-" + hashlib.sha256(
            canonical(sorted(diseases))
        ).hexdigest()[:16],
        RARE_DISEASE_BUNDLE_SHA256,
        profile_store,
    )
    return profile, genes


def _science_disease_subset(
    resources: PublicResourceContext,
    targets,
    *, _skip_verification=False,
) -> tuple[DiseaseProfile | None, tuple[DiseaseGeneEvidence, ...]]:
    gene_store, diseases = _subset_gene_store(
        resources.disease_genes,
        frozenset(value.target for value in targets if value.target is not None
            and value.namespace is not None and not value.namespace.endswith(':locus')),
        frozenset(value.condition for value in targets if value.condition is not None),
    )
    profile_store = _subset_profile_store(
        resources.disease_profiles,
        diseases,
        resources.hpo_store,
    )
    def resource_record(cls, *values):
        if not _skip_verification:
            return cls(*values)
        from genosapiens.wgs_science_graph import _producer_record
        if cls is DiseaseGeneEvidence:
            values = (*values, ())
        return _producer_record(cls, values)
    genes = () if gene_store is None else (
        resource_record(DiseaseGeneEvidence,
            "track1-public-gene-subset-" + gene_store.evidence_sha256[:16],
            RARE_DISEASE_BUNDLE_SHA256,
            gene_store,
        ),
    )
    profile = None if profile_store is None else resource_record(DiseaseProfile,
        "track1-public-profile-subset-" + hashlib.sha256(
            canonical(sorted(diseases))
        ).hexdigest()[:16],
        RARE_DISEASE_BUNDLE_SHA256,
        profile_store,
    )
    return profile, genes




def _selected_technical_assessments(
    context: PreparedAnalysisContext,
    selection,
    *, admission_owner: AuthenticatedAdmissionOwner,
) -> tuple[SelectedTechnicalAssessment, ...]:
    selected: dict[str, SubmissionAlleleTrace] = {}
    for group in selection.selected:
        for trace in group.alleles:
            prior = selected.setdefault(trace.allele_id, trace)
            if prior != trace:
                raise WgsPublicAnalysisError("selected allele trace is inconsistent")
    if not selected:
        return ()
    values: list[SelectedTechnicalAssessment] = []
    reader = admission_owner.reader()
    try:
        for allele_id, trace in sorted(selected.items()):
            expected = {source.source_allele_id: source
                        for source in trace.source_traces}
            if len(expected) != len(trace.source_traces):
                raise WgsPublicAnalysisError("selected source trace is duplicated")
            for item in reader.iter_allele_inputs(allele_id):
                calls = {call.source.source_fields["source_allele_id"]: call
                         for call in item.allele.calls or ()}
                for obligation in item.technical_obligations:
                    source = expected.get(obligation.source_allele_id)
                    if source is None:
                        continue
                    required = obligation.read_validation_required
                    call = calls.get(obligation.source_allele_id)
                    if call is None:
                        raise WgsPublicAnalysisError("selected admission call is missing")
                    assessment = (call.source.source_fields or {}).get("authorized_read_assessment")
                    measurement = None
                    reason = "read-validation-required" if required else "no-mandatory-technical-obligation"
                    status = "unmet" if required else "not_applicable"
                    if source.authority == "authorized_read_evidence":
                        reason = "authorized-read-output"
                        status = "unmet" if required else "satisfied"
                        measurement = source.authorization.stage_output_sha256
                    elif assessment is not None:
                        reason = "read-validation-required"
                        status = "satisfied" if assessment["outcome"] == "validated" else "unmet"
                        measurement = assessment["measurement_artifact_sha256"]
                    quality = (call.source.source_fields or {}).get("quality_acceptance")
                    if quality is not None and quality["filter_value"] != "PASS":
                        # Read evidence cannot retroactively change the supplied FILTER token.
                        status, reason = "unmet", "record-filter-not-PASS"
                    values.append(SelectedTechnicalAssessment(
                        item.allele.allele_id,
                        obligation.source_allele_id,
                        obligation.sample_name,
                        source.source_input_sha256,
                        reason,
                        "selected-source-allele",
                        status,
                        measurement,
                    ))
    except sqlite3.Error:
        raise WgsPublicAnalysisError("selected admission lookup failed") from None
    finally:
        reader.close()
    result = tuple(sorted(values, key=lambda value: (
        value.allele_id,
        value.source_allele_id,
        value.sample_name,
        value.reason,
        value.scope,
    )))
    expected_sources = {
        (
            trace.allele_id,
            source.source_allele_id,
            source.sample_name,
            source.source_input_sha256,
        )
        for trace in selected.values()
        for source in trace.source_traces
    }
    assessed_sources = {
        (
            value.allele_id,
            value.source_allele_id,
            value.sample_name,
            value.source_input_sha256,
        )
        for value in result
    }
    if (
        assessed_sources != expected_sources
        or len(result) != len(expected_sources)
    ):
        raise WgsPublicAnalysisError(
            "selected technical obligation coverage is incomplete"
        )
    return result


def load_analysis_screening_arguments(resources, bindings, indexes, *, dependencies):
    from genosapiens.common_allele_filter import load_common_allele_filter_policy
    from genosapiens.gnomad_faf95_adapter import (
        load_gnomad_faf95_precision_contract, verify_gnomad_faf95_precision_qualification)
    from genosapiens.track1_public_validation import GNOMAD_PRECISION_REPORT_PATH
    common=load_common_allele_filter_policy(read_bounded_metadata(
        resources.repository_root/'config/track1/common-allele-prefilter-v1.json'),
        bindings.common_allele_policy_sha256)
    content=read_bounded_metadata(resources.repository_root/'config/public-resources/gnomad-joint-faf95-precision-v1.json')
    contract=load_gnomad_faf95_precision_contract(content,hashlib.sha256(content).hexdigest())
    precision=verify_gnomad_faf95_precision_qualification(read_bounded_metadata(
        resources.reference_root/GNOMAD_PRECISION_REPORT_PATH),bindings.gnomad_projection_precision_sha256,contract)
    return analysis_screening_arguments(resources,bindings,indexes,common_policy=common,
        precision=precision,dependencies=dependencies)


def analysis_screening_arguments(resources, bindings, indexes, *, common_policy, precision, dependencies):
    """Current screening resources shared by public validation and its fixtures."""
    from genosapiens.wgs_analysis_eligibility import load_analysis_eligibility_policy
    from genosapiens.wgs_context_runtime import load_context_writer_profile
    policy_bytes=read_bounded_metadata(resources.repository_root/'config/track1/analysis-eligibility-v1.json')
    profile_bytes=read_bounded_metadata(resources.repository_root/'config/track1/context-writer-runtime-v1.json')
    policy=load_analysis_eligibility_policy(policy_bytes,hashlib.sha256(policy_bytes).hexdigest())
    profile=load_context_writer_profile(profile_bytes,hashlib.sha256(profile_bytes).hexdigest())
    if policy.common_filter_policy_sha256!=bindings.common_allele_policy_sha256:
        raise WgsPublicAnalysisError('foreign current screening policy')
    return dict(analysis_eligibility_policy=policy,context_writer_profile=profile,
        eligibility_common_policy=common_policy,eligibility_precision=precision,
        eligibility_dependencies=tuple(sorted(dependencies)),
        eligibility_resource_bindings=tuple(sorted(dict(
            population_manifest=resources.population_joint.manifest_sha256,
            precision_manifest=precision.report_sha256,
            clinical_manifest=resources.clinvar.manifest_sha256,
            exception_manifest=bindings.common_allele_exception_manifest_sha256,
            all_sites_manifest=resources.population_all_sites.bundle_manifest_sha256,
            public_index_manifest=indexes.manifest_sha256).items())))


class _ScientificPublicationScope:
    """Private lifetime for one complete product/publication invocation."""
    def __init__(self):
        self.stack = ExitStack()
        self.reader = None
        self.closed = False

    def close(self):
        try:
            self.stack.close()
        finally:
            self.closed = True


@dataclass(frozen=True, slots=True, repr=False)
class PublicBoundedAnalysisBuilder:
    resources: PublicResourceContext
    policies: PublicPolicyFamily
    runtime: PublicPolicyRuntime
    native_source_contract: NativeSourceContract
    executable_stack_sha256: str
    runtime_image_digest: str
    sample_name: str
    clinical_indexes: VerifiedPublicIndexes
    validation_audit: ValidationAuditWriter | None = None
    annotation_source: BoundedAnnotationSource | None = None
    annotation_prefix: object | None = None
    analysis_eligibility_policy: object = field(kw_only=True)
    eligibility_common_policy: object = field(kw_only=True)
    eligibility_precision: object = field(kw_only=True)
    eligibility_resource_bindings: tuple = field(kw_only=True)
    eligibility_dependencies: tuple = field(kw_only=True)
    context_writer_profile: object = field(kw_only=True)
    focused_work_contract: object | None = field(default=None, kw_only=True)
    public_pilot_inputs: object | None = field(default=None, kw_only=True)
    current_run: object | None = field(default=None, kw_only=True, compare=False)
    science_options: object | None = field(default=None, kw_only=True, compare=False)
    science_budget: object | None = field(default=None, kw_only=True, compare=False)
    _publication_scope: object | None = field(default=None, kw_only=True, compare=False)

    def stage_executable(self, stage):
        return self.executable_stack_sha256 if self.current_run is None else self.current_run.code_identity(stage)

    def stage_complete(self, stage, *paths):
        if self.current_run is not None:
            self.current_run.complete(stage, *paths)

    def __post_init__(self) -> None:
        if self.science_options is not None:
            from genosapiens.wgs_science_options import ScienceRunOptions
            from genosapiens.wgs_science_runtime import ScienceBudget, ScienceClock
            if (type(self.science_options) is not ScienceRunOptions or type(self.science_budget) not in (ScienceBudget, ScienceClock)
                    or (type(self.science_budget) is ScienceClock) != self.science_options.no_budget
                    or self.current_run is None
                    or self.current_run.contracts['science'].get('scientific_scope') != self.science_options.policy.sha256):
                raise WgsPublicAnalysisError('scientific execution options have no bound operator owner')
        elif self.science_budget is not None:
            raise WgsPublicAnalysisError('scientific budget has no selected execution options')
        from genosapiens.wgs_analysis_eligibility import AnalysisEligibilityPolicy
        from genosapiens.common_allele_filter import CommonAlleleFilterPolicy
        from genosapiens.gnomad_faf95_adapter import GnomadFaf95PrecisionQualification
        from genosapiens.wgs_context_runtime import ContextWriterProfile
        from genosapiens.track1_analysis_profiles import ProfileCatalog
        from genosapiens.regulatory_index import RegulatoryIndex
        from genosapiens.admission_links import PinnedPublicSource
        if (self._publication_scope is not None and type(self._publication_scope) is not _ScientificPublicationScope
                or type(self.resources) is not PublicResourceContext
                or type(self.resources.profile_catalog) is not ProfileCatalog
                or type(self.resources.regulatory_index) is not RegulatoryIndex
                or type(self.resources.regulatory_sources) is not tuple
                or len(self.resources.regulatory_sources)!=2
                or any(type(pin) is not PinnedPublicSource for pin in self.resources.regulatory_sources)
                or self.resources.profile_catalog.coverage_issues):
            raise WgsPublicAnalysisError('required profile resources are unavailable')
        self.resources.regulatory_index.require_stable()
        if self.focused_work_contract is not None:
            from genosapiens.wgs_work_contract import FocusedWorkContract, load_focused_work_contract
            if (type(self.focused_work_contract) is not FocusedWorkContract
                    or self.focused_work_contract != load_focused_work_contract(self.resources.repository_root)):
                raise WgsPublicAnalysisError('invalid focused work contract')
            if self.public_pilot_inputs is not None:
                raise WgsPublicAnalysisError('focused work requires fresh scientific construction')
        if (type(self.analysis_eligibility_policy) is not AnalysisEligibilityPolicy
                or type(self.eligibility_common_policy) is not CommonAlleleFilterPolicy
                or type(self.eligibility_precision) is not GnomadFaf95PrecisionQualification
                or not self.eligibility_resource_bindings or not self.eligibility_dependencies
                or type(self.context_writer_profile) is not ContextWriterProfile):
            raise WgsPublicAnalysisError('invalid public eligibility resources')
        if self.annotation_prefix is not None:
            from genosapiens.wgs_annotation_prefix import ResolvedAnnotationPrefix
            if type(self.annotation_prefix) is not ResolvedAnnotationPrefix or self.annotation_source is not None:
                raise WgsPublicAnalysisError("public annotation prefix is invalid")
        if (
            type(self.resources) is not PublicResourceContext
            or type(self.policies) is not PublicPolicyFamily
            or type(self.runtime) is not PublicPolicyRuntime
            or type(self.native_source_contract) is not NativeSourceContract
            or len(self.executable_stack_sha256) != 64
            or not self.runtime_image_digest.startswith("sha256:")
            or not self.sample_name
            or type(self.clinical_indexes) is not VerifiedPublicIndexes
            or self.validation_audit is not None and type(self.validation_audit) is not ValidationAuditWriter
            or self.annotation_source is not None and type(self.annotation_source) is not BoundedAnnotationSource
        ):
            raise WgsPublicAnalysisError(
                "public bounded analysis builder is invalid"
            )

    @contextmanager
    def publication_operation(self):
        """One live publication scope; no reader is serialized in checkpoints."""
        if self._publication_scope is not None:
            raise WgsPublicAnalysisError('scientific publication scope is already active')
        scope = _ScientificPublicationScope()
        try:
            yield replace(self, _publication_scope=scope), scope
        except BaseException as error:
            from genosapiens.wgs_failures import run_cleanup
            run_cleanup((('context_store', scope.close),), primary=error)
            raise
        else:
            scope.close()

    def assemble_primary_publication(self, context, products, *, release_manifest_bytes,
            public_qualification_authority=None):
        from genosapiens.wgs_primary_publication import assemble_primary_publication
        scope = self._publication_scope
        if scope is None or scope.closed or scope.reader is None:
            raise WgsPublicAnalysisError('primary publication has no live scientific owner')
        return assemble_primary_publication(context, products, builder=self,
            scientific_reader=scope.reader, release_manifest_bytes=release_manifest_bytes,
            public_qualification_authority=public_qualification_authority)

    def checkpoint_identity(self, context: PreparedAnalysisContext) -> str:
        """Bind saved products to the selected public policy and resource pins."""
        from genosapiens.wgs_profile_projection import profile_adapter_identity
        if type(context) is not PreparedAnalysisContext:
            raise WgsPublicAnalysisError("analysis checkpoint context is invalid")
        return hashlib.sha256(canonical({
            "schema": "genosapiens.wgs-public-product-builder", "version": 2,
            **({"scientific_scope": self.science_options.policy.sha256} if self.science_options is not None else {}),
            "policy_family": self.policies.policy_sha256,
            "policy_runtime": self.runtime.candidate.bindings.sha256,
            "hypothesis_policy": self.policies.hypothesis_scientific_sha256,
            "epcr_policy": self.policies.epcr_policy_sha256,
            "native_source": self.native_source_contract.sha256,
            "annotation_prefix": None if self.annotation_prefix is None else self.annotation_prefix.descriptor.sha256,
            "public_indexes": (
                self.clinical_indexes.index_sha256,
                self.clinical_indexes.manifest_sha256,
                self.clinical_indexes.closure_sha256,
            ),
            "submission_reference": self.resources.submission_reference.manifest_sha256,
            "clinical_open_targets": None if self.resources.clinical_open_targets is None
                else self.resources.clinical_open_targets.manifest_sha256,
            "clinvar": None if self.resources.clinvar is None
                else self.resources.clinvar.manifest_sha256,
            "sample_name": self.sample_name,
            "analysis_eligibility_policy":self.analysis_eligibility_policy.sha256,
            "context_writer_profile":self.context_writer_profile.sha256,
            "focused_work":None if self.focused_work_contract is None else self.focused_work_contract.bindings(),
            "analysis_profile":context.profile_selection.profile,
            "analysis_profile_policy":context.profile_selection.policy.sha256,
            "selected_profile_manifest":context.profile_manifest_sha256,
            "analysis_profile_catalog":self.resources.profile_catalog.sha256,
            "regulatory_index":self.resources.regulatory_index.manifest_sha256,
            "regulatory_sources":tuple((pin.contract_sha256,pin.manifest_sha256)
                for pin in self.resources.regulatory_sources),
            "profile_adapter":profile_adapter_identity(),
            "completed_analysis_input": (None if context.completed_analysis_input is None
                else context.completed_analysis_input.descriptor.sha256),
            "completed_admission_input": (None if context.completed_admission_input is None
                else context.completed_admission_input.binding()),
            "fixed_annotation": None if self.annotation_source is None
                else self.annotation_source.sha256,
            "finite_public_pilot": None if self.public_pilot_inputs is None else {
                "selection_sha256": hashlib.sha256(self.public_pilot_inputs.selection_bytes).hexdigest(),
                "dense_selection_sha256": hashlib.sha256(self.public_pilot_inputs.dense_bytes).hexdigest(),
                "admission_completion_sha256": self.public_pilot_inputs.admission_owner.completion.sha256,
                "completed_context_inputs": (None if self.public_pilot_inputs.completed_context_inputs is None
                    else self.public_pilot_inputs.completed_context_inputs.sha256),
                "completed_assessment_inputs": (None if self.public_pilot_inputs.completed_assessment_inputs is None
                    else self.public_pilot_inputs.completed_assessment_inputs.sha256),
                "completed_scientific_inputs": (None if self.public_pilot_inputs.completed_scientific_inputs is None
                    else self.public_pilot_inputs.completed_scientific_inputs.sha256),
            },
        })).hexdigest()

    def __call__(
        self, context: PreparedAnalysisContext,
    ) -> BoundedAnalysisProducts:
        if (
            type(context) is not PreparedAnalysisContext
            or self.native_source_contract.parent_source_manifest_sha256
            != context.source_manifest.sha256
            or self.native_source_contract
            != context.source_routing_context.contract
            or self.validation_audit is not None and (
                self.validation_audit.parents["prepared_case"] != context.prepared.manifest.sha256
                or self.validation_audit.parents["policy_bindings"] != self.runtime.candidate.bindings.sha256
                or self.validation_audit.config != context.execution_config)
        ):
            raise WgsPublicAnalysisError(
                "public bounded analysis context is invalid"
            )
        try:
            return self._build(context)
        except (WgsPublicAnalysisError, WgsFailure, MemoryError):
            raise
        except Exception:
            raise WgsPublicAnalysisError(
                "public bounded scientific products failed"
            ) from None

    def resume_preflight(self, context: PreparedAnalysisContext) -> None:
        """Refuse an uncheckpointed scientific prefix before opening resources."""
        if type(context) is not PreparedAnalysisContext or not context.resume:
            raise WgsPublicAnalysisError("analysis resume context is invalid")
        scientific_root = context.work_root / "scientific"
        if scientific_root.exists() or scientific_root.is_symlink():
            _private_directory(scientific_root, resume=True)
            if any(scientific_root.iterdir()):
                raise WgsPublicAnalysisError(
                    "scientific product checkpoint is unavailable for resume"
                )
        annotation_root = context.work_root / "annotation"
        if (annotation_root / "annotation-head.json").exists() and not (
                context.work_root / "state/annotation.sqlite").is_file():
            raise WgsPublicAnalysisError("annotation checkpoint state is missing")

    def bind_checkpointed_products(self, context, products, binding):
        """Open one completed owner after the existing journal verifies products."""
        from genosapiens.wgs_scientific_checkpoint import _completed_product_reuse
        from genosapiens.wgs_scientific_products import validate_scientific_products, product_descriptors
        operation = self._publication_scope
        if (not context.resume or operation is None or operation.closed or operation.reader is not None
                or binding.executable_stack_sha256 != self.executable_stack_sha256):
            raise WgsPublicAnalysisError('completed product publication scope differs')
        reuse = _completed_product_reuse(products.scientific, binding)
        reader = operation.stack.enter_context(ScientificArtifactReader(products.scientific,
            directory=products.scientific_directory, expected_parent_bindings=dict(products.scientific.parent_bindings),
            expected_reference_binding=products.selection.reference, scratch_directory=context.work_root,
            eligibility_root=products.eligibility_root, profile_root=products.profile_root,
            progress=context.analysis_progress, authenticated_resources=self.resources.sensitivity_resources,
            completed_reuse=reuse))
        validate_scientific_products(reader, product_descriptors(products), expected_prefilter_manifest_sha256=context.filter_manifest.sha256)
        operation.reader = reader

    def load_completed_products(
        self, context: PreparedAnalysisContext,
    ) -> BoundedAnalysisProducts:
        """Reconstruct only completed manifests after an uncommitted checkpoint write."""
        if type(context) is not PreparedAnalysisContext or not context.resume:
            raise WgsPublicAnalysisError("product checkpoint recovery context is invalid")
        paths = (
            ("primary.json", "primary/primary-manifest.json"),
            ("final-final-evidence-document.json", "reports/diagnostics.json"),
            ("final-scientific-evidence.json", "scientific/scientific-evidence-manifest.json"),
            ("final-hypothesis-generation.json", "scientific/scientific-artifact.json"),
            ("final-ranking.json", "order/order-manifest.json"),
            ("final-submission-selection.json", "selection/selection-manifest.json"),
            ("final-audit-json.json", "reports/audit-manifest.json"),
            ("final-audit-markdown.json", "reports/audit-index.md"),
            ("final-analysis-profile-manifest.json", "analysis/profile/manifest.json"),
        )
        try:
            content = {
                name: read_bounded_metadata(context.work_root / path)
                for name, path in paths
            }
            if context.read_intake is not None:
                content["final-authorized-read-intake.json"] = (
                    context.read_intake.manifest.as_bytes()
                )
            products = _products_from_checkpoint_content(context, content)
            if products.reports.audit_markdown_index != content["final-audit-markdown.json"]:
                raise ValueError
            return products
        except Exception:
            raise WgsPublicAnalysisError(
                "completed product metadata is unavailable for checkpoint recovery"
            ) from None

    def _build(
        self, context: PreparedAnalysisContext,
    ) -> BoundedAnalysisProducts:
        if self.public_pilot_inputs is not None:
            from genosapiens.wgs_profile_projection import PublicProfilePilotInputs
            if (type(self.public_pilot_inputs) is not PublicProfilePilotInputs
                    or context.completed_admission_input is None or context.resume
                    or not context.work_root.is_relative_to('/mnt/reference/scratch/public-build')
                    or not any(re.fullmatch(r'pilot(?:-[0-9]{2})?', p.name) for p in context.work_root.parents)):
                raise WgsPublicAnalysisError('finite cohort requires an explicit fresh public pilot')
            if self.public_pilot_inputs.completed_scientific_inputs is not None:
                return self._reuse_completed_scientific_products(context)
        if self.annotation_prefix is not None and (context.completed_analysis_input is None or context.resume):
            raise WgsPublicAnalysisError("annotation prefix requires fresh completed-input analysis")
        if context.resume:
            self.resume_preflight(context)
        if context.read_intake is not None:
            from genosapiens.wgs_read_rescue import MergedIntakeIndex
            with MergedIntakeIndex(context.read_intake, execution_config=context.execution_config,
                                  scratch_path=context.work_root / "read-intake-verification.sqlite"):
                pass
        if any(value is None for value in (
            self.resources.clinical_open_targets,
            self.resources.clinvar,
        )):
            raise WgsPublicAnalysisError(
                "public clinical resources are unavailable"
            )
        with failure_boundary(
            "resource_session", "resource_session_construct",
        ):
            manager = WgsClinicalResourceSession(
                indexes=self.clinical_indexes,
                open_targets=self.resources.clinical_open_targets,
                clinvar=self.resources.clinvar,
                reference=self.resources.submission_reference,
                execution_config=context.execution_config,
                scratch_path=context.work_root / (
                    "clinical-session" if not context.resume
                    else "clinical-session-resume-" + os.urandom(8).hex()
                ),
            )
        with managed_resource_session(manager) as clinical_session:
            with failure_boundary(
                "bounded_products", "analysis_product_build",
                known_errors=(WgsPublicAnalysisError, AdmissionProjectionError,),
            ):
                try:
                    return self._build_with_session(context, clinical_session)
                except BaseException as error:
                    if self.science_options is not None:
                        from genosapiens.wgs_science_failures import observe_setup_failure
                        observe_setup_failure(error)
                    raise

    def _build_with_session(
        self,
        context: PreparedAnalysisContext,
        clinical_session: WgsClinicalResourceSession,
    ) -> BoundedAnalysisProducts:
        config: BoundedWgsExecutionConfig = context.execution_config
        phenotype = _read_phenotype(context)
        view = scoring_view(phenotype.observations)
        roots = {
            name: _private_directory(context.work_root / name, resume=context.resume)
            for name in (
                "annotation-plan",
                "annotation",
                "scientific",
                "order",
                "primary",
                "selection",
                "state",
            ) if (name != "annotation-plan" or context.completed_analysis_input is None)
                and (name != "annotation" or context.completed_admission_input is None)
        }
        population_limits = PopulationQueryLimits(config.public_queries.batch_rows,
            config.public_queries.batch_rows,65536,1_000_000,
            config.public_queries.max_response_bytes,config.public_queries.timeout_seconds)
        cache_binding = canonical({
            'schema':'genosapiens.wgs-case-allele-query-cache','version':2,
            'consumer_executable_sha256':self.executable_stack_sha256,
            'case_input_sha256':context.prepared.manifest.sha256 if context.read_intake is None
                else context.read_intake.manifest.sha256,
            'reference_role':context.source_routing_context.contract.reference_role.value,
            'reference_identity_id':context.source_routing_context.contract.reference_identity_id,
            'reference_manifest_sha256':context.source_routing_context.contract.reference_manifest_sha256,
            'population_sources':population_source_bindings(joint=self.resources.population_joint,
                all_sites=self.resources.population_all_sites,reference=self.resources.population_reference,
                limits=population_limits),
            'clinical_sources':clinical_source_bindings(clinical_session),
            'query_semantics_sha256':hashlib.sha256(bounded_wgs_execution_config_as_json(config)).hexdigest(),
            'adapter_sha256':query_adapter_identity(),
            'profile':context.profile_selection.profile,
            'profile_policy_sha256':context.profile_selection.policy.sha256,
            'profile_catalog_sha256':self.resources.profile_catalog.sha256,
            'regulatory_index_sha256':self.resources.regulatory_index.manifest_sha256,
        })
        try:
            return self._build_with_query_cache(context,clinical_session,cache_binding,
                roots=roots,phenotype=phenotype,view=view,population_limits=population_limits)
        finally:
            clinical_session.query_cache=None

    def _build_with_query_cache(self,context,clinical_session,cache_binding,*,roots,phenotype,view,population_limits):
        from genosapiens.wgs_population_input import GeneratedPopulationPolicy
        rule, transform, rescues = _population_template(self.policies.policy_bytes)
        with ExitStack() as admission_scope, GeneratedPopulationPolicy(scratch_directory=context.work_root,
                rule=rule, transform=transform, rescues=rescues) as population_input:
            build = self._build_with_population_input if self.science_options is None else self._build_scoped_with_population_input
            return build(context,clinical_session,cache_binding,
                roots=roots,phenotype=phenotype,view=view,population_limits=population_limits,
                population_input=population_input,admission_scope=admission_scope)

    @contextmanager
    def _screening_session(self,context,*,population_limits,cache_path,cache_binding,query_cache=None):
        from genosapiens.wgs_analysis_eligibility import AnalysisEligibilitySession
        from genosapiens.wgs_filter import _FilterState
        from genosapiens.wgs_admission import _source_header
        from genosapiens.population_clinical_session import PopulationClinicalResourceSession
        config=context.execution_config
        with ExitStack() as scope:
            if query_cache is None:
                query_cache=scope.enter_context(CaseAlleleQueryCache(cache_path,binding=cache_binding,execution_config=config))
            source_connection=filter_state=public=None
            merged=None if context.read_intake is None else context.read_intake.manifest
            if merged is None:
                bindings=dict(self.eligibility_resource_bindings)
                filter_state=scope.enter_context(closing(_FilterState.open_sealed(
                    Path(context.prepared_index.preparation_state_root)/'filter.sqlite',
                    expected=context.filter_manifest.summary.total_source_occurrences,
                    source_manifest_sha256=context.source_manifest.sha256,
                    policy_sha256=self.analysis_eligibility_policy.common_filter_policy_sha256,
                    resource_sha256=tuple(bindings[key] for key in (
                        'population_manifest','clinical_manifest','exception_manifest','precision_manifest')),
                    scope_sha256=(context.filter_manifest.scope_policy_sha256, context.filter_manifest.scope_manifest_sha256),
                    reference_sha256=context.source_manifest.source.reference_manifest_sha256,
                    filter_manifest_sha256=context.filter_manifest.sha256)))
                source_connection=scope.enter_context(closing(open_retained_projection(
                    context.retained_projection_path,prepared_index=context.prepared_index,
                    routing=context.source_routing_context,execution_config=config,
                    expected_identities=context.working_manifest.summary.retained_unique_identities,
                    expected_occurrences=context.working_manifest.summary.retained_occurrences)))
                public=scope.enter_context(managed_resource_session(PopulationClinicalResourceSession(
                    joint=self.resources.population_joint,all_sites=self.resources.population_all_sites,
                    reference=self.resources.population_reference,runtime_image=RUNTIME_IMAGE,
                    limits=population_limits,public_query_limits=config.public_queries,clinical_indexes=self.clinical_indexes)))
            yield AnalysisEligibilitySession(policy=self.analysis_eligibility_policy,
                common_policy=self.eligibility_common_policy,prepared_index=context.prepared_index,
                source_manifest=context.source_manifest,source_connection=source_connection,filter_state=filter_state,
                routing=context.source_routing_context,header=_source_header(context.source_root,context.source_manifest),
                selected_sample=self.sample_name if merged is None else merged.sample_id,
                execution_config=config,public_session=public,precision=self.eligibility_precision,
                resource_bindings=self.eligibility_resource_bindings,case_cache=query_cache,merged_intake=merged)

    def _acquire_admission_projection(self,context,*,roots,population_limits,cache_binding,
                                    clinical_session,query_cache,admission_scope):
        """Converge fresh production and the authenticated, immutable stage input."""
        if context.completed_admission_input is not None:
            from genosapiens.wgs_admission_inputs import ResolvedCompletedAdmissionInput
            completed=context.completed_admission_input
            if (type(completed) is not ResolvedCompletedAdmissionInput or context.resume
                    or context.read_intake is not None or self.annotation_source is not None
                    or self.annotation_prefix is None or context.completed_analysis_input is None
                    or self.executable_stack_sha256!=completed.consumer_authority.executable_stack_sha256
                    or self.runtime_image_digest!=completed.consumer_authority.runtime_image_digest
                    or self.analysis_eligibility_policy.sha256!=completed.eligibility_manifest.eligibility_policy_sha256):
                raise WgsPublicAnalysisError("completed admission consumer differs")
            completed.require_consumer(completed.consumer_authority,completed.consumer_release,
                completed_analysis_input=context.completed_analysis_input,annotation_prefix=self.annotation_prefix)
            annotation_source=completed.annotation_source
            roots['annotation-plan']=annotation_source.plan_root
            roots['annotation']=annotation_source.output_root
            projection=PublishedAdmissionProjection(completed.descriptor.original_binding.encode(),
                completed.completion,completed.descriptor.work_root/'eligibility',completed.eligibility_manifest)
            return projection,annotation_source,Path(completed.descriptor.database.path),None
        config=context.execution_config
        if context.completed_analysis_input is not None:
            if self.annotation_source is not None or context.read_intake is not None:
                raise WgsPublicAnalysisError("completed input cannot import annotation or read intake")
            roots["annotation-plan"] = Path(context.completed_analysis_input.descriptor.annotation_plan_root)
            transport = context.completed_analysis_input.transport
        else:
            retained = iter_retained_alleles(
                source_root=context.source_root,
                source_manifest=context.source_manifest,
                working_root=context.working_root,
                working_manifest=context.working_manifest,
                source_routing_context=context.source_routing_context,
                scratch_path=roots["state"] / "retained.sqlite",
                execution_config=config,
                prepared_index=context.prepared_index,
                prepared_state_root=context.prepared.root.parent / "state",
            )
            if context.read_intake is not None:
                retained = context.read_intake.retained(config, source_manifest=context.source_manifest,
                                                       working_manifest=context.working_manifest)
            transport = _phase_call(
                context,
                "annotation_planning",
                write_vep_partitions,
                retained,
                source_manifest=context.source_manifest,
                working_manifest=context.working_manifest,
                output_root=roots["annotation-plan"],
                execution_config=config,
                intake_binding=_intake_binding(context),
                prepared_index=context.prepared_index,
                prepared_state_root=context.prepared.root.parent / "state",
                resume=context.resume,
            )
        vep_contract_bytes = read_bounded_metadata(
            self.resources.repository_root / VEP_CONTRACT_PATH
        )
        contract = json.loads(vep_contract_bytes)
        resource_session_sha256 = _resource_session_sha256(
            self.resources,
            executable_stack_sha256=self.stage_executable('annotation'),
            runtime_image_digest=self.runtime_image_digest,
            vep_contract_bytes=vep_contract_bytes,
        )
        annotation_arguments = dict(
            plan_root=roots["annotation-plan"],
            plan=transport.plan,
            source_manifest=context.source_manifest,
            working_manifest=context.working_manifest,
            output_root=roots["annotation"],
            cache_directory=(
                self.resources.reference_root / contract["cache"]["output_root"]
            ),
            fasta_path=(
                self.resources.reference_root / contract["fasta"]["output_path"]
            ),
            resource_session_sha256=resource_session_sha256,
            bundle_manifest_sha256=VEP_BUNDLE_SHA256,
            executable_stack_sha256=self.stage_executable('annotation'),
            runtime_image_digest=self.runtime_image_digest,
            execution_config=config,
            intake_binding=_intake_binding(context),
            prepared_index=context.prepared_index,
            prepared_state_root=context.prepared.root.parent / "state",
            resume=context.resume,
            run_options_sha256=context.prepared.manifest.preparation_run_options.sha256,
            phase_recorder=context.phase_recorder,
            reference_coverage=_pinned_annotation_coverage(
                self.resources.reference_root / contract['fasta']['output_path'], contract),
        )
        if self.annotation_source is None:
            if self.annotation_prefix is not None:
                annotation_arguments["completed_prefix"] = self.annotation_prefix
            annotation = _phase_call(
                context,
                "annotation_execution",
                run_annotation_partitions,
                **annotation_arguments,
            )
            annotation_source = BoundedAnnotationSource(roots["annotation-plan"], transport.plan,
                                                        roots["annotation"], annotation)
        else:
            annotation_source = self.annotation_source
            _phase_call(
                context,
                "annotation_execution",
                annotation_source.verify,
                expected_plan=transport.plan,
                source_manifest=context.source_manifest,
                working_manifest=context.working_manifest,
                resource_session_sha256=resource_session_sha256,
                bundle_manifest_sha256=VEP_BUNDLE_SHA256,
                executable_stack_sha256=self.stage_executable('annotation'),
                runtime_image_digest=self.runtime_image_digest,
                execution_config=config,
                intake_binding=_intake_binding(context),
                prepared_index=context.prepared_index,
                prepared_state_root=context.prepared.root.parent / "state",
            )
            roots["annotation-plan"] = annotation_source.plan_root
            roots["annotation"] = annotation_source.output_root
            annotation = annotation_source.result

        admission_projection_path = roots["state"] / "admission.sqlite"
        from genosapiens.wgs_profile_projection import FreshAdmissionProfileIntake
        if self.current_run is not None:
            from genosapiens.wgs_run_manifest import atomic_current
            coverage_path = self.current_run.path.parent / 'annotation-coverage.json'
            atomic_current(coverage_path, annotation_coverage_summary(annotation_source, config))
            self.stage_complete('annotation', roots['annotation'] / 'annotation-manifest.json',
                roots['annotation-plan'] / 'manifest.json', coverage_path)
        profile_intake = None
        if not (context.resume and admission_projection_path.is_file()):
            profile_intake=admission_scope.enter_context(FreshAdmissionProfileIntake(
                roots['state']/'profile-intake.sqlite',selection=context.profile_selection,
                catalog=self.resources.profile_catalog,clinical_session=clinical_session,
                regulatory_index=self.resources.regulatory_index,reference=self.resources.submission_reference,
                annotation_manifest_sha256=annotation.manifest.sha256,
                native_source_sha256=(self.native_source_contract.sha256 if context.read_intake is None else
                    hashlib.sha256(admission_bytes(context.read_intake.manifest.source_bindings)).hexdigest()),
                execution_config=config))
        with failure_boundary(
            "bounded_products", "annotation_replay",
            known_errors=(WgsPublicAnalysisError, AdmissionProjectionError, NativeAnnotationError),
        ):
            with NativeAnnotationSession(root=roots["annotation"], resources=self.resources,
                    config=config, runtime_image="genosapiens-track1-pipeline@" + self.runtime_image_digest,
                    contract_path=self.resources.repository_root / VEP_CONTRACT_PATH,
                    bundle_sha256=VEP_BUNDLE_SHA256, plan_sha256=transport.plan.sha256,
                    stage_sha256=annotation.manifest.sha256,
                    fixed=self.annotation_source is not None or self.annotation_prefix is not None,
                    index=annotation_source.native_index) as native_session, \
                    runtime_phase(context.phase_recorder, "admission_projection"):
                # Classify producer failures before resource-session cleanup;
                # cleanup preserves this origin instead of relabeling it unknown.
                with failure_boundary('bounded_products','annotation_replay',
                        known_errors=(WgsPublicAnalysisError,AdmissionProjectionError,NativeAnnotationError)):
                    projection = _publish_admission_projection(
                        admission_projection_path, context=context,
                        plan_root=roots["annotation-plan"], plan=transport.plan,
                        annotation_root=roots["annotation"], annotation=annotation,
                        resource_session_sha256=resource_session_sha256, native_session=native_session,
                        screening_factory=lambda:self._screening_session(context,population_limits=population_limits,
                            # Screening runs concurrently with main-thread
                            # profile intake. Its SQLite connection and facts
                            # therefore require their own exclusive owner.
                            cache_path=clinical_session.scratch_path/'admission-query-cache.sqlite',
                            cache_binding=cache_binding),
                        eligibility_dependencies=self.eligibility_dependencies,
                        runtime_image_digest=self.runtime_image_digest,
                        consumer_executable_sha256=self.stage_executable('selection'),
                        completed_annotation_prefix_sha256=None if self.annotation_prefix is None
                            else self.annotation_prefix.descriptor.sha256,profile_intake=profile_intake)
            annotation_source = replace(annotation_source, native_index=native_session.index)
        return projection,annotation_source,admission_projection_path,profile_intake

    def _reuse_completed_profile_and_discovery(self, context, *, admission_owner,
            pilot_cohort, route_sources, resources, sources):
        """One explicit public recovery branch, before either completed producer."""
        from genosapiens.wgs_profile_projection import (
            _binding, ProfileProjectionReader, stage_completed_profile_projection)
        from genosapiens.wgs_association_projection import (
            _completed_association_metadata, stage_completed_association_projection)
        from genosapiens.candidate_admission import sharded_admission_policy_bytes
        from genosapiens.track1_roster import compatible_context_discovery_policy
        from genosapiens.wgs_runtime import analysis_progress
        progress = getattr(context, 'analysis_progress', None)
        analysis_progress(progress, 'completed_input_reuse', 'metadata',
            reused_stages=(('admission', 'profile') if self.public_pilot_inputs.completed_context_inputs.profile_only
                else ('admission', 'profile', 'association_discovery')))
        inputs = self.public_pilot_inputs.completed_context_inputs
        inputs.require_authority(context.completed_admission_input.consumer_authority,
            consumer_executable_sha256=self.executable_stack_sha256)
        original = inputs.identity
        if inputs.profile_only:
            # A new association domain can borrow only the unchanged coarse
            # profile. Never open or stage the historical discovery/assessments.
            binding = _binding(admission_owner, context.profile_selection, self.resources.profile_catalog,
                self.resources.regulatory_index.manifest_sha256, self.executable_stack_sha256,
                route_sources, pilot_cohort=pilot_cohort, completed_context_inputs=inputs)
            analysis_progress(progress, 'completed_input_reuse', 'metadata',
                reused_stages=('admission', 'profile'))
            with runtime_phase(context.phase_recorder, 'context_formation'):
                reader = stage_completed_profile_projection(context.work_root/'analysis/profile', inputs=inputs,
                    selection=context.profile_selection, catalog=self.resources.profile_catalog,
                    expected_binding=binding, admission_owner=admission_owner, progress=progress)
            return reader, None
        expected_discovery = hashlib.sha256(canonical(dict(
            schema='genosapiens.wgs-association-discovery', version=1,
            admission_completion_sha256=admission_owner.completion.sha256,
            source_policy_sha256=compatible_context_discovery_policy(inputs, sharded_admission_policy_bytes(sources=sources,
                profiles=None, genes=(), ic_store=None,
                clinical_index_manifest_sha256=self.clinical_indexes.manifest_sha256)),
            reference_manifest_sha256=self.native_source_contract.reference_manifest_sha256,
            profile_manifest_sha256=original['profile_manifest']['sha256'],
            consumer_executable_sha256=original['producer_executable_stack_sha256']))).hexdigest()
        if expected_discovery != original['discovery_binding_sha256']:
            raise WgsPublicAnalysisError('completed discovery producer binding differs')
        profile_root = context.work_root/'analysis/profile'
        discovery_root = context.work_root/'associations'
        if profile_root.exists() or discovery_root.exists():
            raise WgsPublicAnalysisError('completed context staging requires fresh destinations')
        binding = _binding(admission_owner, context.profile_selection, self.resources.profile_catalog,
            self.resources.regulatory_index.manifest_sha256, self.executable_stack_sha256,
            route_sources, pilot_cohort=pilot_cohort, completed_context_inputs=inputs)
        source_reader = ProfileProjectionReader(Path(original['profile_manifest']['path']).parent,
            expected_manifest_sha256=original['profile_manifest']['sha256'], selection=context.profile_selection,
            catalog=self.resources.profile_catalog, expected_binding=binding, admission_owner=admission_owner,
            completed_context_inputs=inputs)
        association_arguments = dict(inputs=inputs, admission_owner=admission_owner, resources=resources,
            clinical_index_manifest_sha256=self.clinical_indexes.manifest_sha256,
            reference_role=self.native_source_contract.reference_role,
            reference_identity_id=self.native_source_contract.reference_identity_id,
            reference_manifest_sha256=self.native_source_contract.reference_manifest_sha256)
        metadata = _completed_association_metadata(profile_reader=source_reader, **association_arguments)
        total_files = len(source_reader.shards)+len(metadata.shards)+len(metadata.subject_shards)+2
        profile_bytes = sum(row.size_bytes for row in source_reader.shards)+original['profile_manifest']['size_bytes']
        total_bytes = profile_bytes+sum(row.size_bytes for row in (*metadata.shards,*metadata.subject_shards))+original['association_manifest']['size_bytes']
        offset_files = offset_bytes = 0
        def staging_progress(phase, subphase, *, counts, **_):
            files = offset_files+counts['staging_files_completed']
            size = offset_bytes+counts['staging_bytes_completed']
            analysis_progress(progress, phase, subphase, counts={
                'metadata_completed':2, 'staging_files_completed':files, 'staging_files_total':total_files,
                'staging_bytes_completed':size, 'staging_bytes_total':total_bytes,
                'source_accounting_inputs':admission_owner.row_count,
                'authentication_scans':admission_owner.metrics.get('authentication_scans',0),
                'complete_domain_scans':admission_owner.metrics.get('complete_domain_traversals',0)},
                unit='bytes', completed=size, total=total_bytes)
        with runtime_phase(context.phase_recorder, 'context_formation'):
            try:
                reader = stage_completed_profile_projection(profile_root, inputs=inputs,
                    selection=context.profile_selection, catalog=self.resources.profile_catalog,
                    expected_binding=binding, admission_owner=admission_owner, progress=staging_progress)
                offset_files, offset_bytes = len(source_reader.shards)+1, profile_bytes
                discovered = stage_completed_association_projection(discovery_root,
                    profile_reader=reader, progress=staging_progress, **association_arguments)
                return reader, discovered
            except BaseException as error:
                from genosapiens.wgs_failures import run_cleanup
                import shutil
                run_cleanup(tuple(('context_store', lambda path=path: shutil.rmtree(path))
                    for path in (profile_root, discovery_root) if path.exists()), primary=error)
                raise

    def _completed_assessments(self, context):
        inputs = None if self.public_pilot_inputs is None else self.public_pilot_inputs.completed_assessment_inputs
        if inputs is not None:
            from genosapiens.wgs_assessment_checkpoint import CompletedAssessmentInputs
            if (type(inputs) is not CompletedAssessmentInputs or self.public_pilot_inputs.completed_context_inputs is None
                    or self.public_pilot_inputs.completed_context_inputs.profile_only):
                raise WgsPublicAnalysisError('assessment restart requires completed profile and discovery')
            inputs.require_authority(context.completed_admission_input.consumer_authority,
                consumer_executable_sha256=self.executable_stack_sha256)
        return inputs

    def _build_with_population_input(self,context,clinical_session,cache_binding,*,roots,phenotype,view,population_limits,population_input,admission_scope):
        config=context.execution_config
        # Clinical facts can refer to this session's raw-record cache. A new
        # session must own new facts, including on resume. Keep legacy state
        # caches untouched and reuse durable stages through their own journals.
        query_cache=admission_scope.enter_context(CaseAlleleQueryCache(clinical_session.scratch_path/'allele-query-cache.sqlite',
            binding=cache_binding,execution_config=config))
        clinical_session.query_cache=query_cache
        projection,annotation_source,admission_projection_path,profile_intake=self._acquire_admission_projection(
            context,roots=roots,population_limits=population_limits,cache_binding=cache_binding,
            clinical_session=clinical_session,query_cache=query_cache,admission_scope=admission_scope)
        annotation=annotation_source.result
        admission_projection_binding=projection.binding
        def authentication_progress(record):
            content=canonical(dict(schema='genosapiens.admission-authentication-progress',version=1,
                admission_completion_sha256=projection.completion.sha256,**record))
            temporary=roots['state']/'admission-authentication-progress.tmp'
            fd=os.open(temporary,os.O_WRONLY|os.O_CREAT|os.O_TRUNC|os.O_NOFOLLOW,0o600)
            with os.fdopen(fd,'wb') as handle:handle.write(content)
            os.replace(temporary,roots['state']/'admission-authentication-progress.json')
        with runtime_phase(context.phase_recorder,'prepared_authority_validation'):
            if self.public_pilot_inputs is None:
                admission_owner=admission_scope.enter_context(closing(AuthenticatedAdmissionOwner(admission_projection_path,
                    expected_binding=admission_projection_binding,expected_completion_sha256=projection.completion.sha256,
                    execution_config=config,progress=authentication_progress)))
            else:
                admission_owner=self.public_pilot_inputs.admission_owner
                if (type(admission_owner) is not AuthenticatedAdmissionOwner
                        or admission_owner.path != admission_projection_path
                        or admission_owner.binding != admission_projection_binding
                        or admission_owner.completion.sha256 != projection.completion.sha256
                        or admission_owner.config != config):
                    raise WgsPublicAnalysisError('finite pilot borrowed admission is foreign')
                admission_owner._assert_stable()
        def observe_admission_reads(kind,value,traceback):
            # ExitStack runs this before closing the owner, on success and
            # failure, after discovery, science, selected checks and any audit.
            content=canonical(dict(schema='genosapiens.admission-read-progress',version=1,
                phase='complete' if kind is None else 'failed',
                admission_completion_sha256=admission_owner.completion.sha256,
                metrics=dict(admission_owner.metrics),cache_bytes=admission_owner.cache_bytes,
                query_cache=dict(query_cache.metrics)))
            temporary=roots['state']/'admission-reads.tmp'
            try:
                fd=os.open(temporary,os.O_WRONLY|os.O_CREAT|os.O_TRUNC|os.O_NOFOLLOW,0o600)
                with os.fdopen(fd,'wb') as handle:handle.write(content)
                os.replace(temporary,roots['state']/'admission-reads.json')
            except OSError:
                if kind is None:raise
            return False
        admission_scope.push(observe_admission_reads)
        from genosapiens.wgs_profile_projection import (
            build_profile_projection,ProfileProjectionReader,_binding as profile_binding,
            analysis_profile_component,profile_component_parents)
        repeat_contract_sha256 = hashlib.sha256(
            read_bounded_metadata(self.resources.repository_root / REPEAT_CONTRACT_PATH)
        ).hexdigest()
        sources = ShardedAdmissionSources(
            annotation.manifest.sha256,
            self.resources.regulatory_sources,
            self.resources.clinical_open_targets,
            self.resources.clinvar,
            self.resources.submission_reference,
            RepeatCatalogInput(
                self.resources.repeat_catalog,
                repeat_contract_sha256,
                REPEAT_MANIFEST_SHA256,
                self.resources.repeat_reference,
            ),
            self.resources.regulatory_index.manifest_sha256,
        )
        from genosapiens.candidate_admission import prepare_sharded_resources,sharded_admission_policy_bytes
        from genosapiens.wgs_association_projection import (
            discover_associations,AssociationProjectionReader,AssessmentProjectionOwner,assessment_parents)
        resource_preparation=prepare_sharded_resources(sources)
        pilot_cohort=None
        if self.public_pilot_inputs is not None:
            from genosapiens.wgs_profile_projection import (authenticate_public_profile_pilot_cohort,
                ReconciledPublicProfilePilotFacts,profile_route_sources)
            inputs=self.public_pilot_inputs
            pilot_cohort=authenticate_public_profile_pilot_cohort(self.resources.repository_root,admission_owner,
                selection_bytes=inputs.selection_bytes,dense_bytes=inputs.dense_bytes)
            route_sources=profile_route_sources(annotation_manifest_sha256=admission_owner.parents['annotation_stage_sha256'],
                native_source_sha256=admission_owner.parents['source_routing_sha256'],
                clinical_session=clinical_session,regulatory_index=self.resources.regulatory_index)
            if inputs.completed_context_inputs is None:
                profile_intake=admission_scope.enter_context(ReconciledPublicProfilePilotFacts(
                    inputs.reconciled_directory,repository=self.resources.repository_root,admission_owner=admission_owner,
                    pilot_cohort=pilot_cohort,selection=context.profile_selection,catalog=self.resources.profile_catalog,
                    route_sources=route_sources,regulatory_index_sha256=self.resources.regulatory_index.manifest_sha256))
        elif profile_intake is None:
            if self.current_run is not None and self.current_run.reusable('selection'):
                from types import SimpleNamespace
                from genosapiens.wgs_profile_projection import profile_route_sources
                profile_intake = SimpleNamespace(route_sources=profile_route_sources(
                    annotation_manifest_sha256=admission_owner.parents['annotation_stage_sha256'],
                    native_source_sha256=admission_owner.parents['source_routing_sha256'],
                    clinical_session=clinical_session, regulatory_index=self.resources.regulatory_index))
            else:
                from genosapiens.wgs_profile_projection import AdmissionProfileFactsProvider
                profile_intake = admission_scope.enter_context(AdmissionProfileFactsProvider(
                    admission_owner=admission_owner, annotation_source=annotation_source,
                    source_manifest=context.source_manifest, working_manifest=context.working_manifest,
                    clinical_session=clinical_session, regulatory_index=self.resources.regulatory_index,
                    reference=self.resources.submission_reference, prepared_index=context.prepared_index,
                    prepared_state_root=context.prepared.root.parent / 'state'))
        profile_root=context.work_root/'analysis/profile'
        profile_root.parent.mkdir(mode=0o700,exist_ok=True)
        completed_inputs = (None if self.public_pilot_inputs is None
            else self.public_pilot_inputs.completed_context_inputs)
        assessment_inputs = self._completed_assessments(context)
        if completed_inputs is not None:
            profile_reader, discovered = self._reuse_completed_profile_and_discovery(context,
                admission_owner=admission_owner, pilot_cohort=pilot_cohort, route_sources=route_sources,
                resources=resource_preparation, sources=sources)
            profile_reader = admission_scope.enter_context(profile_reader)
            profile_sha256 = profile_reader.manifest_sha256
        elif self.current_run is not None and self.current_run.reusable('selection'):
            profile_sha256 = hashlib.sha256(read_bounded_metadata(profile_root/'manifest.json')).hexdigest()
            profile_reader = admission_scope.enter_context(ProfileProjectionReader(profile_root,
                expected_manifest_sha256=profile_sha256, selection=context.profile_selection,
                catalog=self.resources.profile_catalog, expected_binding=profile_binding(admission_owner,
                    context.profile_selection, self.resources.profile_catalog,
                    self.resources.regulatory_index.manifest_sha256, self.stage_executable('selection'),
                    profile_intake.route_sources), admission_owner=admission_owner))
        else:
            with runtime_phase(context.phase_recorder,'context_formation'):
                profile_sha256=build_profile_projection(profile_root,admission_owner=admission_owner,
                    selection=context.profile_selection,catalog=self.resources.profile_catalog,
                    regulatory_index_sha256=self.resources.regulatory_index.manifest_sha256,
                    consumer_executable_sha256=self.stage_executable('selection'),
                    route_sources=profile_intake.route_sources,facts_provider=profile_intake,
                    pilot_cohort=pilot_cohort,
                    maximum_shard_rows=config.source_shards.max_uncompressed_rows,
                    maximum_shard_bytes=config.source_shards.max_uncompressed_bytes)
                profile_reader=admission_scope.enter_context(ProfileProjectionReader(profile_root,
                    expected_manifest_sha256=profile_sha256,selection=context.profile_selection,
                    catalog=self.resources.profile_catalog,expected_binding=profile_binding(admission_owner,
                        context.profile_selection,self.resources.profile_catalog,
                        self.resources.regulatory_index.manifest_sha256,self.stage_executable('selection'),
                        profile_intake.route_sources,pilot_cohort=pilot_cohort),admission_owner=admission_owner,
                    verification_parent=context.work_root))
        self.stage_complete('selection', profile_root/'manifest.json',
            roots['state']/'admission-stage/admission-completion.json')
        profile_component=analysis_profile_component(profile_reader)
        context=context.with_profile_projection(profile_reader)
        discovery_root=context.work_root/'associations'
        if completed_inputs is not None and not completed_inputs.profile_only:
            discovery_binding = discovered.binding_sha256
        else:
            discovery_binding=hashlib.sha256(canonical(dict(
                schema='genosapiens.wgs-association-discovery',version=1,
                admission_completion_sha256=projection.completion.sha256,
                source_policy_sha256=hashlib.sha256(sharded_admission_policy_bytes(sources=sources,
                    profiles=None,genes=(),ic_store=None,
                    clinical_index_manifest_sha256=self.clinical_indexes.manifest_sha256)).hexdigest(),
                reference_manifest_sha256=self.native_source_contract.reference_manifest_sha256,
                profile_manifest_sha256=profile_sha256,
                consumer_executable_sha256=self.stage_executable('associations')))).hexdigest()
            if self.current_run is not None and self.current_run.reusable('associations'):
                from genosapiens.wgs_association_projection import load_association_projection
                raw = read_bounded_metadata(discovery_root/'manifest.json')
                discovered = load_association_projection(raw, expected_sha256=hashlib.sha256(raw).hexdigest())
                if discovered.binding_sha256 != discovery_binding:
                    raise WgsPublicAnalysisError('completed discovery dependencies differ')
            else:
                with runtime_phase(context.phase_recorder,'context_formation'):
                    discovered=discover_associations(
                        root=discovery_root,binding_sha256=discovery_binding,
                        admission_completion_sha256=projection.completion.sha256,
                        expected_rows=profile_reader.manifest['counts']['advancing_rows'],
                        resources=resource_preparation,clinical_session=clinical_session,admission_owner=admission_owner,
                        profile_reader=profile_reader,regulatory_index=self.resources.regulatory_index,
                        reference_role=self.native_source_contract.reference_role,
                        reference_identity_id=self.native_source_contract.reference_identity_id,
                        reference_manifest_sha256=self.native_source_contract.reference_manifest_sha256,execution_config=config,
                        progress_callback=context.analysis_progress)
        self.stage_complete('associations', discovery_root/'manifest.json')
        from genosapiens.wgs_runtime import analysis_progress
        analysis_progress(context.analysis_progress, 'context_assembly', 'resource_prepare')
        # These source subsets are also inputs to scientific publication; a
        # snapshot replaces graph/assessment production, not resource ownership.
        with runtime_phase(context.phase_recorder, 'context_formation'):
            profile,genes=_disease_subset(self.resources,discovered.targets)
        if assessment_inputs is not None:
            prepared_admission = None
            admission_policy_sha256 = assessment_inputs.metadata()['bindings']['policy_sha256']
        else:
            with runtime_phase(context.phase_recorder, 'context_formation'):
                prepared_admission = prepare_sharded_candidate_admission(
                    sources=sources,
                    view=view,
                    profiles=profile,
                    genes=genes,
                    ic_store=self.resources.hpo_store,
                    reached_targets=discovered.targets,resource_preparation=resource_preparation,
                    execution_config=config,
                    graph_progress_path=roots["state"] / "graph-support-progress.json",
                    progress=context.analysis_progress, authenticated_resources=self.resources.sensitivity_resources,
                    clinical_index_manifest_sha256=(
                        self.clinical_indexes.manifest_sha256
                    ),
                )
            admission_policy_sha256 = prepared_admission.policy_sha256
        assessment_owner=AssessmentProjectionOwner(discovery_root,discovered,
            policy_sha256=admission_policy_sha256,profile_reader=profile_reader,progress=context.analysis_progress,
            publication_counts=None if assessment_inputs is None else assessment_inputs.metadata()['counts'])

        def admissions(executor):
            with failure_boundary('bounded_products','annotation_replay',
                    known_errors=(WgsPublicAnalysisError,AdmissionProjectionError)):
                with closing(AssociationProjectionReader(discovery_root,expected_sha256=discovered.sha256,
                        binding_sha256=discovery_binding,admission_completion_sha256=projection.completion.sha256,
                        resources=resource_preparation,admission_owner=admission_owner,profile_reader=profile_reader)) as associations, \
                        closing(associations.batches(max_rows=config.public_queries.batch_rows)) as batches:
                    for batch in batches:
                        with failure_boundary('bounded_products','candidate_admission_convert',
                                known_errors=(CandidateAdmissionError,)):
                            result=_generate_projected_candidate_admission(discovered=batch,
                                prepared=prepared_admission,context_executor=executor,
                                expected_policy_sha256=admission_policy_sha256)
                        yield result
        def complete_context_stream():
            nonlocal prepared_admission
            from genosapiens.wgs_context_admission import ContextAdmissionExecutor
            with closing(admission_owner.reader()) as reader:
                if assessment_inputs is not None:
                    from genosapiens.wgs_assessment_checkpoint import iter_recovered_contexts
                    with closing(iter_recovered_contexts(assessment_inputs, admission_reader=reader,
                            profile_reader=profile_reader, assessment_owner=assessment_owner,
                            scratch_path=context.work_root/'contexts.sqlite', execution_config=config,
                            execution_config_sha256=context.prepared.manifest.execution_limits_sha256,
                            native_contract_sha256=self.native_source_contract.sha256,
                            population_input=population_input, models=_MODELS,
                            progress=context.analysis_progress)) as stream:
                        yield from stream
                    return
                with closing(ContextAdmissionExecutor(prepared_admission, execution_config=config,
                        total_inputs=reader.row_count,
                        profile_reader=profile_reader, progress=context.analysis_progress,
                        progress_path=roots['state'] / 'context-formation-progress.json')) as executor:
                    with closing(iter_closed_contexts(admissions(executor),admission_reader=reader,
                            scratch_path=context.work_root / 'contexts.sqlite',execution_config=config,
                            context_executor=executor, validation_audit=self.validation_audit,
                            writer_profile=self.context_writer_profile,assessment_owner=assessment_owner,
                            profile_reader=profile_reader)) as stream:
                        first = next(stream, None)
                        if first is not None:
                            for key in sorted(executor.context_keys):
                                seed = prepared_admission.graph.seed_contexts[key][0]
                                if seed.native_disease_id is not None and seed.inheritance is not None:
                                    for model in _MODELS:
                                        population_input.add(frequency.ContextKey(seed.native_disease_id,
                                            seed.gene_id, seed.locus_id, seed.inheritance, model))
                        executor.release_preparation()
                        prepared_admission = None
                        if first is not None:
                            yield first
                        yield from stream

        with runtime_phase(context.phase_recorder, "context_formation"):
            with failure_boundary(
                "bounded_products", "context_assemble",
                # Context/graph worker errors are WgsFailure and retain their own operation/category.
                known_errors=(WgsContextError,),
            ):
                contexts = complete_context_stream()
                first_context = next(contexts, None)
        prepared_population_policy = population_input.seal()
        parent_bindings = {
            "prepared_case": context.prepared.manifest.sha256,
            "source_manifest": context.source_manifest.sha256,
            "prefilter_manifest": context.filter_manifest.sha256,
            "working_set": context.working_manifest.sha256,
            "execution_config": context.prepared.manifest.execution_limits_sha256,
            "policy_bindings": self.runtime.candidate.bindings.sha256,
            **({"read_intake": context.read_intake.manifest.sha256} if context.read_intake else {}),
        }
        from genosapiens.wgs_analysis_eligibility import analysis_eligibility_component,eligibility_component_parents
        eligibility_component=analysis_eligibility_component(projection.eligibility_root,projection.eligibility_manifest)
        parent_bindings.update(eligibility_component_parents(eligibility_component))
        parent_bindings.update(profile_component_parents(profile_component))
        parent_bindings.update(assessment_parents(discovered,assessment_owner.manifest))
        if self.focused_work_contract is not None:
            parent_bindings.update(self.focused_work_contract.bindings())
        expected_public_identity_id = self.resources.reference_contract.reference_role(
            "public_annotation_benchmark").identity_id
        with ExitStack() as science_scope:
            science_scope.callback(contexts.close)
            # Admission has globally closed before the first context is yielded.
            # Empty contexts still retain their closure, but need no live query session.
            population_manager = PopulationResourceSession(
                joint=self.resources.population_joint, all_sites=self.resources.population_all_sites,
                reference=self.resources.population_reference, runtime_image=RUNTIME_IMAGE,
                limits=population_limits)
            work_arguments = {}
            if self.focused_work_contract is not None:
                work_source = science_scope.enter_context(self._screening_session(context,
                    population_limits=population_limits, cache_path=roots['scientific']/'work-query-cache.sqlite',
                    cache_binding={}, query_cache=query_cache))
                work_arguments = dict(work_filter_policy=self.focused_work_contract.policy,
                    work_filter_source=work_source)

            def enriched_contexts():
                entered = False
                for closed in (() if first_context is None else chain((first_context,), contexts)):
                    if closed.member_count and not entered:
                        if any(value is None for value in (self.resources.population_joint,
                                self.resources.population_all_sites, self.resources.population_reference)):
                            raise WgsPublicAnalysisError("public population resources are unavailable")
                        science_scope.enter_context(managed_resource_session(population_manager))
                        entered = True
                    yield closed

            scientific = _phase_call(context, "scientific_analysis", write_scientific_artifact,
                enriched_contexts(),
                output_directory=roots["scientific"], execution_config=config,
                run_manifest_sha256=hashlib.sha256(canonical({
                    "schema": "genosapiens.wgs-normalized-scientific-run", "version": 1,
                    "parents": parent_bindings, "executable_stack": self.executable_stack_sha256,
                    "runtime_image": self.runtime_image_digest,
                })).hexdigest(),
                analysis_eligibility=eligibility_component,
                analysis_profile=profile_component,
                observations=phenotype.observations, view=view, profile=profile, gene_evidence=genes,
                ic_store=self.resources.hpo_store,
                expected_ic_sha256=_checked_ic_hash(self.resources.hpo_store, self.resources.sensitivity_resources),
                prepared_population_policy=prepared_population_policy,
                hypothesis_policy=self.policies.hypothesis_policy,
                hypothesis_policy_sha256=hashlib.sha256(hypothesis_policy_bytes(self.policies.hypothesis_policy,
                    admission_policy_sha256=admission_policy_sha256)).hexdigest(),
                admission_policy_sha256=admission_policy_sha256,
                fusion_policy=self.runtime.fusion_policy,
                fusion_policy_sha256=self.runtime.candidate.bindings.fusion_policy_sha256,
                prior_policy=self.runtime.prior_policy, diagnosis=self.policies.disclosed_diagnosis,
                sample_name=self.sample_name,
                source_manifest_sha256=(context.source_manifest.sha256 if context.read_intake is None
                    else context.read_intake.manifest.sha256),
                parent_bindings=parent_bindings, reference_pin=self.resources.submission_reference,
                expected_public_identity_id=expected_public_identity_id, population_session=population_manager,
                clinical_annotator=lambda alleles: clinical_session.annotate(alleles), query_cache=query_cache,
                authenticated_resources=self.resources.sensitivity_resources, progress=context.analysis_progress,
                context_total=None if context.analysis_progress is None else context.analysis_progress.counts.get('contexts_obtained'),
                **work_arguments)
        self.stage_complete('science', roots['scientific']/'scientific-artifact.json')
        from genosapiens.wgs_scientific_checkpoint import _completed_producer_reuse
        return self._publish_scientific_products(context, roots=roots, scientific=scientific,
            parent_bindings=parent_bindings, eligibility_root=projection.eligibility_root,
            profile_root=profile_root, admission_owner=admission_owner,
            annotation_source=annotation_source, completed_reuse=_completed_producer_reuse(scientific,
                consumer_executable_sha256=self.executable_stack_sha256))


    # Retain the legacy entry above for exact upstream checkpoint compatibility.
    # The scoped entry has the same upstream prefix and a separate science contract.
    def _build_scoped_with_population_input(self,context,clinical_session,cache_binding,*,roots,phenotype,view,population_limits,population_input,admission_scope):
        import time
        science_setup_started = time.monotonic()
        config=context.execution_config
        # Clinical facts can refer to this session's raw-record cache. A new
        # session must own new facts, including on resume. Keep legacy state
        # caches untouched and reuse durable stages through their own journals.
        query_cache=admission_scope.enter_context(CaseAlleleQueryCache(clinical_session.scratch_path/'allele-query-cache.sqlite',
            binding=cache_binding,execution_config=config))
        clinical_session.query_cache=query_cache
        projection,annotation_source,admission_projection_path,profile_intake=self._acquire_admission_projection(
            context,roots=roots,population_limits=population_limits,cache_binding=cache_binding,
            clinical_session=clinical_session,query_cache=query_cache,admission_scope=admission_scope)
        annotation=annotation_source.result
        admission_projection_binding=projection.binding
        def authentication_progress(record):
            content=canonical(dict(schema='genosapiens.admission-authentication-progress',version=1,
                admission_completion_sha256=projection.completion.sha256,**record))
            temporary=roots['state']/'admission-authentication-progress.tmp'
            fd=os.open(temporary,os.O_WRONLY|os.O_CREAT|os.O_TRUNC|os.O_NOFOLLOW,0o600)
            with os.fdopen(fd,'wb') as handle:handle.write(content)
            os.replace(temporary,roots['state']/'admission-authentication-progress.json')
        with runtime_phase(context.phase_recorder,'prepared_authority_validation'):
            if self.public_pilot_inputs is None:
                admission_owner=admission_scope.enter_context(closing(AuthenticatedAdmissionOwner(admission_projection_path,
                    expected_binding=admission_projection_binding,expected_completion_sha256=projection.completion.sha256,
                    execution_config=config,progress=authentication_progress)))
            else:
                admission_owner=self.public_pilot_inputs.admission_owner
                if (type(admission_owner) is not AuthenticatedAdmissionOwner
                        or admission_owner.path != admission_projection_path
                        or admission_owner.binding != admission_projection_binding
                        or admission_owner.completion.sha256 != projection.completion.sha256
                        or admission_owner.config != config):
                    raise WgsPublicAnalysisError('finite pilot borrowed admission is foreign')
                admission_owner._assert_stable()
        def observe_admission_reads(kind,value,traceback):
            # ExitStack runs this before closing the owner, on success and
            # failure, after discovery, science, selected checks and any audit.
            content=canonical(dict(schema='genosapiens.admission-read-progress',version=1,
                phase='complete' if kind is None else 'failed',
                admission_completion_sha256=admission_owner.completion.sha256,
                metrics=dict(admission_owner.metrics),cache_bytes=admission_owner.cache_bytes,
                query_cache=dict(query_cache.metrics)))
            temporary=roots['state']/'admission-reads.tmp'
            try:
                fd=os.open(temporary,os.O_WRONLY|os.O_CREAT|os.O_TRUNC|os.O_NOFOLLOW,0o600)
                with os.fdopen(fd,'wb') as handle:handle.write(content)
                os.replace(temporary,roots['state']/'admission-reads.json')
            except OSError:
                if kind is None:raise
            return False
        admission_scope.push(observe_admission_reads)
        from genosapiens.wgs_profile_projection import (
            build_profile_projection,ProfileProjectionReader,_binding as profile_binding,
            analysis_profile_component,profile_component_parents)
        repeat_contract_sha256 = hashlib.sha256(
            read_bounded_metadata(self.resources.repository_root / REPEAT_CONTRACT_PATH)
        ).hexdigest()
        sources = ShardedAdmissionSources(
            annotation.manifest.sha256,
            self.resources.regulatory_sources,
            self.resources.clinical_open_targets,
            self.resources.clinvar,
            self.resources.submission_reference,
            RepeatCatalogInput(
                self.resources.repeat_catalog,
                repeat_contract_sha256,
                REPEAT_MANIFEST_SHA256,
                self.resources.repeat_reference,
            ),
            self.resources.regulatory_index.manifest_sha256,
        )
        from genosapiens.candidate_admission import prepare_sharded_resources,sharded_admission_policy_bytes
        from genosapiens.wgs_association_projection import (
            discover_associations,AssociationProjectionReader,AssessmentProjectionOwner,assessment_parents)
        resource_preparation=prepare_sharded_resources(sources)
        pilot_cohort=None
        if self.public_pilot_inputs is not None:
            from genosapiens.wgs_profile_projection import (authenticate_public_profile_pilot_cohort,
                ReconciledPublicProfilePilotFacts,profile_route_sources)
            inputs=self.public_pilot_inputs
            pilot_cohort=authenticate_public_profile_pilot_cohort(self.resources.repository_root,admission_owner,
                selection_bytes=inputs.selection_bytes,dense_bytes=inputs.dense_bytes)
            route_sources=profile_route_sources(annotation_manifest_sha256=admission_owner.parents['annotation_stage_sha256'],
                native_source_sha256=admission_owner.parents['source_routing_sha256'],
                clinical_session=clinical_session,regulatory_index=self.resources.regulatory_index)
            if inputs.completed_context_inputs is None:
                profile_intake=admission_scope.enter_context(ReconciledPublicProfilePilotFacts(
                    inputs.reconciled_directory,repository=self.resources.repository_root,admission_owner=admission_owner,
                    pilot_cohort=pilot_cohort,selection=context.profile_selection,catalog=self.resources.profile_catalog,
                    route_sources=route_sources,regulatory_index_sha256=self.resources.regulatory_index.manifest_sha256))
        elif profile_intake is None:
            if self.current_run is not None and self.current_run.reusable('selection'):
                from types import SimpleNamespace
                from genosapiens.wgs_profile_projection import profile_route_sources
                profile_intake = SimpleNamespace(route_sources=profile_route_sources(
                    annotation_manifest_sha256=admission_owner.parents['annotation_stage_sha256'],
                    native_source_sha256=admission_owner.parents['source_routing_sha256'],
                    clinical_session=clinical_session, regulatory_index=self.resources.regulatory_index))
            else:
                from genosapiens.wgs_profile_projection import AdmissionProfileFactsProvider
                profile_intake = admission_scope.enter_context(AdmissionProfileFactsProvider(
                    admission_owner=admission_owner, annotation_source=annotation_source,
                    source_manifest=context.source_manifest, working_manifest=context.working_manifest,
                    clinical_session=clinical_session, regulatory_index=self.resources.regulatory_index,
                    reference=self.resources.submission_reference, prepared_index=context.prepared_index,
                    prepared_state_root=context.prepared.root.parent / 'state'))
        profile_root=context.work_root/'analysis/profile'
        profile_root.parent.mkdir(mode=0o700,exist_ok=True)
        completed_inputs = (None if self.public_pilot_inputs is None
            else self.public_pilot_inputs.completed_context_inputs)
        assessment_inputs = self._completed_assessments(context)
        if completed_inputs is not None:
            profile_reader, discovered = self._reuse_completed_profile_and_discovery(context,
                admission_owner=admission_owner, pilot_cohort=pilot_cohort, route_sources=route_sources,
                resources=resource_preparation, sources=sources)
            profile_reader = admission_scope.enter_context(profile_reader)
            profile_sha256 = profile_reader.manifest_sha256
        elif self.current_run is not None and self.current_run.reusable('selection'):
            profile_sha256 = hashlib.sha256(read_bounded_metadata(profile_root/'manifest.json')).hexdigest()
            profile_reader = admission_scope.enter_context(ProfileProjectionReader(profile_root,
                expected_manifest_sha256=profile_sha256, selection=context.profile_selection,
                catalog=self.resources.profile_catalog, expected_binding=profile_binding(admission_owner,
                    context.profile_selection, self.resources.profile_catalog,
                    self.resources.regulatory_index.manifest_sha256, self.stage_executable('selection'),
                    profile_intake.route_sources), admission_owner=admission_owner))
        else:
            with runtime_phase(context.phase_recorder,'context_formation'):
                profile_sha256=build_profile_projection(profile_root,admission_owner=admission_owner,
                    selection=context.profile_selection,catalog=self.resources.profile_catalog,
                    regulatory_index_sha256=self.resources.regulatory_index.manifest_sha256,
                    consumer_executable_sha256=self.stage_executable('selection'),
                    route_sources=profile_intake.route_sources,facts_provider=profile_intake,
                    pilot_cohort=pilot_cohort,
                    maximum_shard_rows=config.source_shards.max_uncompressed_rows,
                    maximum_shard_bytes=config.source_shards.max_uncompressed_bytes)
                profile_reader=admission_scope.enter_context(ProfileProjectionReader(profile_root,
                    expected_manifest_sha256=profile_sha256,selection=context.profile_selection,
                    catalog=self.resources.profile_catalog,expected_binding=profile_binding(admission_owner,
                        context.profile_selection,self.resources.profile_catalog,
                        self.resources.regulatory_index.manifest_sha256,self.stage_executable('selection'),
                        profile_intake.route_sources,pilot_cohort=pilot_cohort),admission_owner=admission_owner,
                    verification_parent=context.work_root))
        self.stage_complete('selection', profile_root/'manifest.json',
            roots['state']/'admission-stage/admission-completion.json')
        profile_component=analysis_profile_component(profile_reader)
        context=context.with_profile_projection(profile_reader)
        discovery_root=context.work_root/'associations'
        if completed_inputs is not None and not completed_inputs.profile_only:
            discovery_binding = discovered.binding_sha256
        else:
            discovery_binding=hashlib.sha256(canonical(dict(
                schema='genosapiens.wgs-association-discovery',version=1,
                admission_completion_sha256=projection.completion.sha256,
                source_policy_sha256=hashlib.sha256(sharded_admission_policy_bytes(sources=sources,
                    profiles=None,genes=(),ic_store=None,
                    clinical_index_manifest_sha256=self.clinical_indexes.manifest_sha256)).hexdigest(),
                reference_manifest_sha256=self.native_source_contract.reference_manifest_sha256,
                profile_manifest_sha256=profile_sha256,
                consumer_executable_sha256=self.stage_executable('associations')))).hexdigest()
            if self.current_run is not None and self.current_run.reusable('associations'):
                from genosapiens.wgs_association_projection import load_association_projection
                raw = read_bounded_metadata(discovery_root/'manifest.json')
                discovered = load_association_projection(raw, expected_sha256=hashlib.sha256(raw).hexdigest())
                if discovered.binding_sha256 != discovery_binding:
                    raise WgsPublicAnalysisError('completed discovery dependencies differ')
            else:
                with runtime_phase(context.phase_recorder,'context_formation'):
                    discovered=discover_associations(
                        root=discovery_root,binding_sha256=discovery_binding,
                        admission_completion_sha256=projection.completion.sha256,
                        expected_rows=profile_reader.manifest['counts']['advancing_rows'],
                        resources=resource_preparation,clinical_session=clinical_session,admission_owner=admission_owner,
                        profile_reader=profile_reader,regulatory_index=self.resources.regulatory_index,
                        reference_role=self.native_source_contract.reference_role,
                        reference_identity_id=self.native_source_contract.reference_identity_id,
                        reference_manifest_sha256=self.native_source_contract.reference_manifest_sha256,execution_config=config,
                        progress_callback=context.analysis_progress)
        self.stage_complete('associations', discovery_root/'manifest.json')
        from genosapiens.wgs_runtime import analysis_progress
        analysis_progress(context.analysis_progress, 'context_assembly', 'resource_prepare')
        # These source subsets are also inputs to scientific publication; a
        # snapshot replaces graph/assessment production, not resource ownership.
        with runtime_phase(context.phase_recorder, 'context_formation'):
            profile,genes=_science_disease_subset(self.resources,discovered.targets, _skip_verification=True)
        if assessment_inputs is not None:
            prepared_admission = None
            admission_policy_sha256 = assessment_inputs.metadata()['bindings']['policy_sha256']
        else:
            from genosapiens.wgs_science_graph import prepare_sharded_candidate_admission
            with runtime_phase(context.phase_recorder, 'context_formation'):
                prepared_admission = prepare_sharded_candidate_admission(
                    sources=sources,
                    view=view,
                    profiles=profile,
                    genes=genes,
                    ic_store=self.resources.hpo_store,
                    reached_targets=discovered.targets,resource_preparation=resource_preparation,
                    execution_config=config,
                    graph_progress_path=roots["state"] / "graph-support-progress.json",
                    progress=context.analysis_progress, authenticated_resources=self.resources.sensitivity_resources,
                    clinical_index_manifest_sha256=(
                        self.clinical_indexes.manifest_sha256
                    ),
                )
            admission_policy_sha256 = prepared_admission.policy_sha256
        assessment_owner=AssessmentProjectionOwner(discovery_root,discovered,
            policy_sha256=admission_policy_sha256,profile_reader=profile_reader,progress=context.analysis_progress,
            publication_counts=None if assessment_inputs is None else assessment_inputs.metadata()['counts'])

        def admissions(executor):
            with failure_boundary('bounded_products','annotation_replay',
                    known_errors=(WgsPublicAnalysisError,AdmissionProjectionError)):
                with closing(AssociationProjectionReader(discovery_root,expected_sha256=discovered.sha256,
                        binding_sha256=discovery_binding,admission_completion_sha256=projection.completion.sha256,
                        resources=resource_preparation,admission_owner=admission_owner,profile_reader=profile_reader)) as associations, \
                        closing(associations.batches(max_rows=config.public_queries.batch_rows)) as batches:
                    if science_selector is not None and profile_reader is not None:
                        associations.profile_rows.close()
                        associations.profile_rows = science_selector.collect_profile_rows()
                    for batch in batches:
                        with failure_boundary('bounded_products','candidate_admission_convert',
                                known_errors=(CandidateAdmissionError,)):
                            result=_generate_projected_candidate_admission(discovered=batch,
                                prepared=prepared_admission,context_executor=executor,
                                expected_policy_sha256=admission_policy_sha256)
                        yield result
        science_selector = None
        if self.science_options is not None:
            from genosapiens.wgs_science_selection import ScientificContextSelector
            from genosapiens.phenotype_scoring import prepare_phenotype_scores
            from genosapiens.wgs_sensitivity import _profiles
            from genosapiens.wgs_phenotype_evidence import _hash
            roster_connection = sqlite3.connect('')
            admission_scope.callback(roster_connection.close)
            science_selector = ScientificContextSelector(policy=self.science_options.policy, roster_connection=roster_connection,
                profile_reader=profile_reader, gene_evidence=genes,
                phenotype_scores=prepare_phenotype_scores(view.observations, ic_store=self.resources.hpo_store,
                    profile_store=None if profile is None else profile.store,
                    prior_identity="diagnosis-blind-primary:" + _hash(view), _skip_verification=True,
                    _ic_store_sha256=getattr(self.resources.sensitivity_resources, 'store_sha256', None)))

        def complete_context_stream():
            nonlocal prepared_admission
            from genosapiens.wgs_science_transport import scientific_context_executor as ContextAdmissionExecutor
            with closing(admission_owner.reader()) as reader:
                if assessment_inputs is not None:
                    from genosapiens.wgs_assessment_checkpoint import iter_recovered_contexts
                    with closing(iter_recovered_contexts(assessment_inputs, admission_reader=reader,
                            profile_reader=profile_reader, assessment_owner=assessment_owner,
                            scratch_path=context.work_root/'contexts.sqlite', execution_config=config,
                            execution_config_sha256=context.prepared.manifest.execution_limits_sha256,
                            native_contract_sha256=self.native_source_contract.sha256,
                            population_input=population_input, models=_MODELS,
                            progress=context.analysis_progress, science_selector=science_selector)) as stream:
                        yield from stream
                    return
                with closing(ContextAdmissionExecutor(prepared_admission, execution_config=config,
                        total_inputs=reader.row_count,
                        profile_reader=profile_reader, progress=context.analysis_progress,
                        progress_path=roots['state'] / 'context-formation-progress.json')) as executor:
                    with closing(iter_closed_contexts(admissions(executor),admission_reader=reader,
                            scratch_path=context.work_root / 'contexts.sqlite',execution_config=config,
                            context_executor=executor, validation_audit=self.validation_audit,
                            writer_profile=self.context_writer_profile,assessment_owner=assessment_owner,
                            profile_reader=profile_reader, science_selector=science_selector)) as stream:
                        first = next(stream, None)
                        if first is not None:
                            for key in sorted(executor.context_keys):
                                seed = prepared_admission.graph.seed_contexts[key][0]
                                if seed.native_disease_id is not None and seed.inheritance is not None:
                                    for model in _MODELS:
                                        population_input.add(frequency.ContextKey(seed.native_disease_id,
                                            seed.gene_id, seed.locus_id, seed.inheritance, model))
                        executor.release_preparation()
                        prepared_admission = None
                        if first is not None:
                            yield first
                        yield from stream

        with runtime_phase(context.phase_recorder, "context_formation"):
            with failure_boundary(
                "bounded_products", "context_assemble",
                # Context/graph worker errors are WgsFailure and retain their own operation/category.
                known_errors=(WgsContextError,),
            ):
                contexts = complete_context_stream()
                from contextlib import nullcontext
                from genosapiens.wgs_science_runtime import bounded_operation
                with nullcontext() if self.science_budget is None else bounded_operation(self.science_budget,
                        maximum_seconds=max(0.001, self.science_budget.remaining())):
                    first_context = next(contexts, None)
        from genosapiens.wgs_run_manifest import atomic_current
        atomic_current(roots['state'] / 'context-science-selection.json', dict(
            schema='genosapiens.science-selection-performance/v1', **science_selector.metrics,
            setup_seconds=time.monotonic() - science_setup_started))
        prepared_population_policy = population_input.seal()
        parent_bindings = {
            "prepared_case": context.prepared.manifest.sha256,
            "source_manifest": context.source_manifest.sha256,
            "prefilter_manifest": context.filter_manifest.sha256,
            "working_set": context.working_manifest.sha256,
            "execution_config": context.prepared.manifest.execution_limits_sha256,
            "policy_bindings": self.runtime.candidate.bindings.sha256,
            **({"read_intake": context.read_intake.manifest.sha256} if context.read_intake else {}),
        }
        from genosapiens.wgs_analysis_eligibility import analysis_eligibility_component,eligibility_component_parents
        eligibility_component=analysis_eligibility_component(projection.eligibility_root,projection.eligibility_manifest)
        parent_bindings.update(eligibility_component_parents(eligibility_component))
        parent_bindings.update(profile_component_parents(profile_component))
        parent_bindings.update(assessment_parents(discovered,assessment_owner.manifest))
        if self.focused_work_contract is not None:
            parent_bindings.update(self.focused_work_contract.bindings())
        expected_public_identity_id = self.resources.reference_contract.reference_role(
            "public_annotation_benchmark").identity_id
        computation_sha256 = hashlib.sha256(canonical({
            "schema": "genosapiens.wgs-normalized-scientific-run", "version": 1,
            "parents": parent_bindings, "executable_stack": self.executable_stack_sha256,
            "runtime_image": self.runtime_image_digest,
        })).hexdigest()
        with ExitStack() as science_scope:
            science_scope.callback(contexts.close)
            coordinator = None
            if science_selector is not None:
                from genosapiens.wgs_science_coordinator import ScienceCoordinator
                from genosapiens.wgs_science_checkpoints import ContextCheckpointStore
                directory = sqlite3.connect(roots['state'] / 'context-science-directory.sqlite')
                science_scope.callback(directory.close)
                directory.execute('DROP TABLE IF EXISTS science_context_directory')
                metadata = self.current_run.path.parent
                coordinator = ScienceCoordinator(directory, science_selector.roster,
                    checkpoints=ContextCheckpointStore(metadata / 'science-contexts'),
                    computation_sha256=computation_sha256, budget=self.science_budget,
                    scratch_directory=context.work_root,
                    maximum_contexts=32 if self.science_options.pilot else None)
            # Admission has globally closed before the first context is yielded.
            # Empty contexts still retain their closure, but need no live query session.
            population_manager = PopulationResourceSession(
                joint=self.resources.population_joint, all_sites=self.resources.population_all_sites,
                reference=self.resources.population_reference, runtime_image=RUNTIME_IMAGE,
                limits=population_limits)
            work_arguments = {}
            if self.focused_work_contract is not None:
                work_source = science_scope.enter_context(self._screening_session(context,
                    population_limits=population_limits, cache_path=roots['state']/'context-work-query-cache.sqlite',
                    cache_binding={}, query_cache=query_cache))
                work_arguments = dict(work_filter_policy=self.focused_work_contract.policy,
                    work_filter_source=work_source)

            def enriched_contexts():
                entered = False
                for closed in (() if first_context is None else chain((first_context,), contexts)):
                    if closed.member_count and not entered:
                        if any(value is None for value in (self.resources.population_joint,
                                self.resources.population_all_sites, self.resources.population_reference)):
                            raise WgsPublicAnalysisError("public population resources are unavailable")
                        science_scope.enter_context(managed_resource_session(population_manager))
                        entered = True
                    yield closed

            def worker_sessions(scope):
                import tempfile
                from genosapiens.wgs_science_workers import open_worker_query_cache
                worker_root = Path(tempfile.mkdtemp(prefix='science-worker-', dir=context.work_root))
                population = scope.enter_context(managed_resource_session(PopulationResourceSession(
                    joint=self.resources.population_joint, all_sites=self.resources.population_all_sites,
                    reference=self.resources.population_reference, runtime_image=RUNTIME_IMAGE,
                    limits=population_limits)))
                clinical = scope.enter_context(managed_resource_session(WgsClinicalResourceSession(
                    indexes=self.clinical_indexes, open_targets=self.resources.clinical_open_targets,
                    clinvar=self.resources.clinvar, reference=self.resources.submission_reference,
                    execution_config=config, scratch_path=worker_root / 'clinical')))
                cache = scope.enter_context(open_worker_query_cache(worker_root / 'queries.sqlite',
                    binding=query_cache.binding, clinical=clinical, execution_config=config))
                clinical.query_cache = cache
                source = (None if self.focused_work_contract is None else scope.enter_context(
                    self._screening_session(context, population_limits=population_limits,
                        cache_path=worker_root / 'work.sqlite', cache_binding={}, query_cache=cache)))
                return population, clinical.annotate, cache, source

            scientific = _phase_call(context, "scientific_analysis", write_scientific_artifact,
                enriched_contexts(), worker_sessions=worker_sessions,
                science_workers=self.science_options.workers,
                publication_readers=self.science_options.publication_readers,
                output_directory=roots["scientific"], execution_config=config,
                run_manifest_sha256=computation_sha256, science_coordinator=coordinator,
                analysis_eligibility=eligibility_component,
                analysis_profile=profile_component,
                observations=phenotype.observations, view=view, profile=profile, gene_evidence=genes,
                ic_store=self.resources.hpo_store,
                expected_ic_sha256=_checked_ic_hash(self.resources.hpo_store, self.resources.sensitivity_resources),
                prepared_population_policy=prepared_population_policy,
                hypothesis_policy=self.policies.hypothesis_policy,
                hypothesis_policy_sha256=hashlib.sha256(hypothesis_policy_bytes(self.policies.hypothesis_policy,
                    admission_policy_sha256=admission_policy_sha256)).hexdigest(),
                admission_policy_sha256=admission_policy_sha256,
                fusion_policy=self.runtime.fusion_policy,
                fusion_policy_sha256=self.runtime.candidate.bindings.fusion_policy_sha256,
                prior_policy=self.runtime.prior_policy, diagnosis=self.policies.disclosed_diagnosis,
                sample_name=self.sample_name,
                source_manifest_sha256=(context.source_manifest.sha256 if context.read_intake is None
                    else context.read_intake.manifest.sha256),
                parent_bindings=parent_bindings, reference_pin=self.resources.submission_reference,
                expected_public_identity_id=expected_public_identity_id, population_session=population_manager,
                clinical_annotator=lambda alleles: clinical_session.annotate(alleles), query_cache=query_cache,
                authenticated_resources=self.resources.sensitivity_resources, progress=context.analysis_progress,
                context_total=None if context.analysis_progress is None else context.analysis_progress.counts.get('contexts_obtained'),
                **work_arguments)
        from genosapiens.wgs_scientific_artifact import scientific_artifact_complete
        if scientific_artifact_complete(scientific):
            self.stage_complete('science', roots['scientific']/'scientific-artifact.json')
        from genosapiens.wgs_scientific_checkpoint import _completed_producer_reuse
        return self._publish_scientific_products(context, roots=roots, scientific=scientific,
            parent_bindings=parent_bindings, eligibility_root=projection.eligibility_root,
            profile_root=profile_root, admission_owner=admission_owner,
            annotation_source=annotation_source, completed_reuse=_completed_producer_reuse(scientific,
                consumer_executable_sha256=self.executable_stack_sha256) if scientific_artifact_complete(scientific) else None)

    def reuse_current_scientific(self, context):
        """Render current reports from matching completed scientific outputs."""
        from genosapiens.wgs_scientific_artifact import scientific_artifact_from_bytes
        from genosapiens.wgs_scientific_checkpoint import _completed_producer_reuse
        if self.current_run is None or not self.current_run.reusable('science'):
            raise WgsPublicAnalysisError('scientific dependencies are not reusable')
        work = context.work_root
        scientific_bytes = read_bounded_metadata(work/'scientific/scientific-artifact.json')
        scientific = scientific_artifact_from_bytes(scientific_bytes,
            expected_sha256=hashlib.sha256(scientific_bytes).hexdigest())
        parents = dict(scientific.parent_bindings)
        profile_root = work/'analysis/profile'
        profile_bytes = read_bounded_metadata(profile_root/'manifest.json')
        profile = json.loads(profile_bytes)
        if (parents.get('prepared_case') != context.prepared.manifest.sha256
                or parents.get('analysis_profile_manifest') != hashlib.sha256(profile_bytes).hexdigest()):
            raise WgsPublicAnalysisError('completed scientific parents differ')
        path = work/'state/admission.sqlite'
        with closing(sqlite3.connect(path.as_uri()+'?mode=ro&immutable=1', uri=True)) as db:
            size, = db.execute('SELECT length(canonical_bytes) FROM run_binding').fetchone()
            if not 0 < size <= 2*1024**2:
                raise WgsPublicAnalysisError('completed admission binding is invalid')
            binding, = db.execute('SELECT canonical_bytes FROM run_binding').fetchone()
        if hashlib.sha256(binding).hexdigest() != profile['parent']['binding_sha256']:
            raise WgsPublicAnalysisError('completed admission binding differs')
        with closing(AuthenticatedAdmissionOwner(path, expected_binding=binding,
                expected_completion_sha256=profile['parent']['completion_sha256'],
                execution_config=context.execution_config)) as owner:
            roots = {name:_private_directory(work/name, resume=True)
                for name in ('scientific','order','primary','selection')}
            return self._publish_scientific_products(context, roots=roots, scientific=scientific,
                parent_bindings=parents, eligibility_root=work/'eligibility', profile_root=profile_root,
                admission_owner=owner, annotation_source=None,
                annotation_coverage=json.loads(read_bounded_metadata(
                    self.current_run.path.parent/'annotation-coverage.json')),
                completed_reuse=_completed_producer_reuse(scientific,
                    consumer_executable_sha256=self.executable_stack_sha256))

    def _reuse_completed_scientific_products(self, context):
        """Consume explicitly authorized sealed science, preserving its producer."""
        from genosapiens.wgs_scientific_checkpoint import CompletedScientificInputs
        from genosapiens.wgs_profile_projection import scientific_profile_component
        inputs = self.public_pilot_inputs.completed_scientific_inputs
        if (type(inputs) is not CompletedScientificInputs
                or self.public_pilot_inputs.completed_context_inputs is not None
                or self.public_pilot_inputs.completed_assessment_inputs is not None
                or context.read_intake is not None or self.validation_audit is not None):
            raise WgsPublicAnalysisError('scientific restart inputs are incompatible')
        completed = context.completed_admission_input
        if (self.annotation_source is not None or self.annotation_prefix is None
                or context.completed_analysis_input is None
                or self.executable_stack_sha256 != completed.consumer_authority.executable_stack_sha256
                or self.runtime_image_digest != completed.consumer_authority.runtime_image_digest):
            raise WgsPublicAnalysisError('scientific restart completed-input consumer differs')
        completed.require_consumer(completed.consumer_authority, completed.consumer_release,
            completed_analysis_input=context.completed_analysis_input, annotation_prefix=self.annotation_prefix)
        inputs.require_authority(context.completed_admission_input.consumer_authority,
            consumer_executable_sha256=self.executable_stack_sha256)
        inputs.require_completed_publication()
        owner = self.public_pilot_inputs.admission_owner
        if type(owner) is not AuthenticatedAdmissionOwner:
            raise WgsPublicAnalysisError('scientific restart admission owner is invalid')
        owner._assert_stable()
        from genosapiens.wgs_runtime import analysis_progress
        analysis_progress(context.analysis_progress, 'publication', 'scientific_input_staging')
        scope = context.work_root
        scientific = inputs.stage_work_root(scope)
        completed_reuse = inputs.completed_reuse(completed.consumer_authority,
            consumer_executable_sha256=self.executable_stack_sha256, artifact=scientific)
        parents = dict(scientific.parent_bindings)
        expected = dict(prepared_case=context.prepared.manifest.sha256,
            source_manifest=context.source_manifest.sha256, prefilter_manifest=context.filter_manifest.sha256,
            working_set=context.working_manifest.sha256, execution_config=context.prepared.manifest.execution_limits_sha256,
            policy_bindings=self.runtime.candidate.bindings.sha256,
            analysis_eligibility_policy=self.analysis_eligibility_policy.sha256,
            analysis_profile_policy=context.profile_selection.policy.sha256,
            analysis_profile_catalog=self.resources.profile_catalog.sha256,
            analysis_profile_regulatory_index=self.resources.regulatory_index.manifest_sha256)
        if (any(parents.get(key) != value for key, value in expected.items())
                or context.profile_manifest_sha256 is not None
                    and context.profile_manifest_sha256 != parents.get('analysis_profile_manifest')):
            raise WgsPublicAnalysisError('scientific restart parents differ from consumer')
        profile = scientific_profile_component(scientific, scope / 'scientific')
        identity = profile['identity']
        if (identity['profile'] != context.profile_selection.profile
                or identity['parent']['binding_sha256'] != owner.binding_sha256
                or identity['parent']['completion_sha256'] != owner.completion.sha256
                or owner.config != context.execution_config
                or owner.completion.sha256 != context.completed_admission_input.descriptor.expected_completion.sha256
                or owner.binding != context.completed_admission_input.descriptor.original_binding.encode('ascii')):
            raise WgsPublicAnalysisError('scientific restart profile or admission differs')
        from genosapiens.wgs_profile_projection import authenticate_public_profile_pilot_cohort
        cohort = authenticate_public_profile_pilot_cohort(self.resources.repository_root, owner,
            selection_bytes=self.public_pilot_inputs.selection_bytes, dense_bytes=self.public_pilot_inputs.dense_bytes)
        expected_cohort = dict(cohort.identity)
        # The pinned old profile keeps its original acceptance provenance. The
        # new amendment binds that exact profile/science, not a relabelled copy.
        expected_cohort['acceptance_sha256'] = identity.get('finite_cohort', {}).get('acceptance_sha256')
        if identity.get('finite_cohort') != expected_cohort:
            raise WgsPublicAnalysisError('scientific restart cohort differs')
        context = replace(context, profile_manifest_sha256=profile['manifest_sha256'])
        roots = dict(scientific=scope / 'scientific')
        for name in ('order', 'primary', 'selection'):
            roots[name] = _private_directory(context.work_root / name, resume=False)
        return self._publish_scientific_products(context, roots=roots, scientific=scientific,
            parent_bindings=parents, eligibility_root=scope / 'eligibility', profile_root=scope / 'analysis/profile',
            admission_owner=owner, annotation_source=None, completed_reuse=completed_reuse)

    def _publish_scientific_products(self, context, *, roots, scientific,
            parent_bindings, eligibility_root, profile_root, admission_owner, annotation_source, completed_reuse=None, annotation_coverage=None):
        """Validate sealed science and derive its normal products and reports.

        This private boundary does not authorize saved-input reuse. Callers own
        the original parent bindings and all retained dependencies; the normal
        completed reader authenticates consumed data and reuses producer science.
        """
        from genosapiens.wgs_runtime import analysis_progress
        from genosapiens.wgs_scientific_checkpoint import CompletedScientificReuse
        if type(completed_reuse) is not CompletedScientificReuse:
            raise WgsPublicAnalysisError('publication requires compatible completed scientific authority')
        completed_reuse.require_artifact(scientific)
        if completed_reuse.consumer_binding != self.executable_stack_sha256:
            raise WgsPublicAnalysisError('completed scientific consumer executable changed')
        operation = self._publication_scope
        if operation is not None and (operation.closed or operation.reader is not None):
            raise WgsPublicAnalysisError('scientific publication reader is unavailable')
        config = context.execution_config
        expected_public_identity_id = self.resources.reference_contract.reference_role(
            "public_annotation_benchmark").identity_id
        analysis_progress(context.analysis_progress, 'publication', 'completed_science_binding')
        with _open_submission_reference(self.resources.submission_reference,
                expected_public_identity_id=expected_public_identity_id) as opened:
            reference_binding = opened.binding
        with ExitStack() as local:
            owner = local if operation is None else operation.stack
            scientific_reader = owner.enter_context(ScientificArtifactReader(scientific, directory=roots["scientific"],
                expected_parent_bindings=parent_bindings, expected_reference_binding=reference_binding,
                scratch_directory=context.work_root,eligibility_root=eligibility_root,
                profile_root=profile_root,progress=context.analysis_progress,
                authenticated_resources=self.resources.sensitivity_resources, completed_reuse=completed_reuse))
            if operation is not None:
                operation.reader = scientific_reader
            with scientific_reader.reader.read_scope():
                if hasattr(scientific_reader, 'scientific_scope'):
                    from genosapiens.wgs_science_reporting import publish_scope_report
                    publish_scope_report(scientific_reader, directory=roots['scientific'])
                    if scientific.completion_status == 'incomplete':
                        raise WgsPublicAnalysisError('scientific analysis is incomplete; see the private scientific scope report')
                analysis_progress(context.analysis_progress, 'publication', 'products')
                products = scientific_products(scientific_reader,
                    expected_prefilter_manifest_sha256=context.filter_manifest.sha256)
                order, primary, selection = products.order, products.primary, products.selection
                for directory, name, content in (
                        (roots["scientific"], "scientific-evidence-manifest.json", scientific.evidence.as_bytes()),
                        (roots["order"], "order-manifest.json", order.as_bytes()),
                        (roots["primary"], "primary-manifest.json", primary.as_bytes()),
                        (roots["selection"], "selection-manifest.json", selection.as_bytes())):
                    from genosapiens.wgs_ledger import publish_bytes
                    path = directory / name
                    if path.exists() and completed_reuse is not None:
                        if read_bounded_metadata(path) != content:
                            raise WgsPublicAnalysisError('completed scientific product differs')
                    else:
                        publish_bytes(directory, name, content)
                if self.validation_audit is not None:
                    self.validation_audit.annotations(scientific_reader)
                technical = _selected_technical_assessments(
                    context,
                    selection,
                    admission_owner=admission_owner,
                )
                qc, normalization, classes = _read_source_summaries(context)
                prepared_resources = {
                    value.role: value.sha256 for value in context.prepared.manifest.resources
                }
                analysis_progress(context.analysis_progress, 'publication', 'reports')
                reports = _phase_call(
                    context,
                    "report_publication",
                    build_bounded_reports,
                    source_root=context.source_root,
                    source_manifest=context.source_manifest,
                    source_accounting=context.source_accounting,
                    source_routing_context=context.source_routing_context,
                    qc_summary=qc,
                    normalization_summary=normalization,
                    class_summary=classes,
                    filter_root=context.filter_root,
                    filter_manifest=context.filter_manifest,
                    working_root=context.working_root,
                    working_manifest=context.working_manifest,
                    scope_policy_sha256=context.filter_manifest.scope_policy_sha256,
                    scope_manifest_sha256=context.filter_manifest.scope_manifest_sha256,
                    policy_sha256=context.prepared.manifest.prefilter_policy_sha256,
                    population_manifest_sha256=(
                        prepared_resources["population_manifest"]
                    ),
                    clinical_manifest_sha256=(
                        prepared_resources["clinical_manifest"]
                    ),
                    exception_manifest_sha256=(
                        prepared_resources["exception_manifest"]
                    ),
                    precision_manifest_sha256=(
                        prepared_resources["precision_manifest"]
                    ),
                    scientific=scientific,
                    scientific_reader=scientific_reader,
                    scientific_directory=roots["scientific"],
                    order=order,
                    order_directory=roots["order"],
                    primary=primary,
                    selection=selection,
                    selection_directory=roots["selection"],
                    selected_technical_assessments=technical,
                    credibility_threshold="0",
                    execution_config=config,
                    output_directory=context.work_root / "reports",
                    maximum_total_json_bytes=config.child_process.max_total_output_bytes,
                    maximum_total_markdown_bytes=config.child_process.max_total_output_bytes,
                    read_intake=context.read_intake,
                    retained_filter_state_path=Path(context.prepared.prepared_index.preparation_state_root) / "filter.sqlite",
                    profile_source_scope=profile_source_accounting(
                        self.resources.regulatory_index, self.resources.profile_catalog),
                    annotation_coverage=(annotation_coverage_summary(annotation_source, config)
                        if annotation_coverage is None else annotation_coverage),
                )
        artifacts = (
            FinalArtifact(FinalRole.EVIDENCE, reports.diagnostic.as_bytes()),
            FinalArtifact(FinalRole.PROVENANCE,
                read_bounded_metadata(roots["scientific"] / "scientific-evidence-manifest.json")),
            FinalArtifact(FinalRole.HYPOTHESES, scientific.as_bytes()),
            FinalArtifact(FinalRole.RANKING, order.as_bytes()),
            FinalArtifact(FinalRole.SELECTION, selection.as_bytes()),
            FinalArtifact(FinalRole.PROFILE, read_bounded_metadata(profile_root/'manifest.json')),
            FinalArtifact(FinalRole.AUDIT_JSON, reports.audit_json_index),
            FinalArtifact(
                FinalRole.AUDIT_MARKDOWN, reports.audit_markdown_index,
            ),
        )
        if context.read_intake is not None:
            artifacts += (FinalArtifact(FinalRole.READ_INTAKE, context.read_intake.manifest.as_bytes()),)
        return BoundedAnalysisProducts(
            scientific,
            roots["scientific"],
            order,
            roots["order"],
            primary,
            selection,
            roots["selection"],
            reports,
            artifacts,
            config,
            read_intake=context.read_intake,
            annotation_source=annotation_source,
            eligibility_root=eligibility_root,
            profile_root=profile_root,
        )
