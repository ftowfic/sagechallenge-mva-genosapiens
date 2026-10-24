"""Normalized complete scientific artifact producer and operation-scoped reader.

The descriptor is written only after the evidence owner has finished primary
order, selection and selected-group replay. Reopening authenticates the complete evidence before
constructing one private case projection shared by that operation's consumers.
The production adapter switches to this contract together with its consumers.
"""
from contextlib import ExitStack, closing, contextmanager
from dataclasses import dataclass
import hashlib
import json
import os
from pathlib import Path
import re
import sys

from genosapiens.bounded_wgs_config import BoundedWgsExecutionConfig
from genosapiens.evidence_model import _sha, _typed
from genosapiens.evidence_serialization import _canonical, _loads
from genosapiens.genotype_evaluation import prepare_gene_evidence
from genosapiens.track1_submission import _open_submission_reference
from genosapiens.track1_submission import SubmissionRepresentationAssessment
from genosapiens.wgs_case_context import ScientificCaseContext
from genosapiens.wgs_case_evidence import CompleteScientificCase, VerifiedScientificCase
from genosapiens.wgs_case_resources import ScientificCaseResources, VerifiedScientificCaseResources, bound_reference
from genosapiens.wgs_context_evidence import EnrichedContextEvidence, VerifiedContextEvidence
from genosapiens.wgs_contexts import ClosedContext
from genosapiens.wgs_evidence_store import EvidenceManifest, EvidenceReference, EvidenceReader, EvidenceWriter, EvidenceStoreError
from genosapiens.wgs_gene_evidence import CaseGeneEvidence
from genosapiens.wgs_phenotype_evidence import CasePhenotypeEvidence
from genosapiens.wgs_population_evidence import CompletePopulationPolicy
from genosapiens.wgs_reference_evidence import (
    METADATA_KIND, SubmissionReferenceMetadata, VerifiedSubmissionReferenceMetadata,
)

SCHEMA = "genosapiens.wgs-normalized-scientific-artifact"
MAX_ARTIFACT_BYTES = 8 * 1024**2
DESCRIPTOR_NAME = "scientific-artifact.json"


class _ScienceFailureObservation:
    """Record code locations without exception text, values, or completion authority."""
    OPERATIONS = {
        'resources': 'analysis_product_build',
        'member_enrichment': 'context_annotate',
        'member_read': 'context_annotation_inputs',
        'population_query': 'context_population_query',
        'population_bind': 'context_population_bind',
        'clinical_query': 'context_clinical_query',
        'clinical_bind': 'context_clinical_bind',
        'member_encode': 'context_payload_encode',
        'member_publication': 'context_annotation_finalize',
        'population_work_filter': 'population_assess',
        'genotype_facts': 'hypothesis_generate',
        'hypothesis_evidence': 'hypothesis_generate',
        'population': 'population_assess',
        'phenotype': 'evidence_adapt',
        'scores': 'evidence_adapt',
        'representation': 'evidence_adapt',
        'eligibility': 'evidence_adapt',
        'ordering': 'evidence_adapt',
        'context_publication': 'scientific_manifest_write',
        'primary_seal': 'scientific_manifest_write',
        'case_seal': 'scientific_manifest_write',
        'evidence_finish': 'scientific_manifest_write',
        'descriptor_publication': 'scientific_manifest_write',
    }

    def __init__(self, root, run_manifest_sha256, progress, *, terminal=True):
        self.root, self.run_manifest_sha256, self.progress = root, run_manifest_sha256, progress
        self.terminal = terminal
        self.subphase, self.context_ordinal, self.primary = 'resources', None, None

    def stage(self, subphase):
        if subphase in self.OPERATIONS:
            self.subphase = subphase

    @staticmethod
    def _frames(error):
        from collections import deque
        library = Path(__file__).parent
        frames = deque(maxlen=16)
        traceback = error.__traceback__
        seen = 0
        while traceback is not None and seen < 128:
            frame = traceback.tb_frame
            module = frame.f_globals.get('__name__')
            if type(module) is str and re.fullmatch(r'genosapiens\.[a-z0-9_]+', module):
                expected = library / (module.split('.')[1] + '.py')
                if frame.f_code.co_filename == str(expected):
                    frames.append((module, traceback.tb_lineno, expected))
            seen += 1
            traceback = traceback.tb_next
        result = []
        for module, line, path in frames:
            # Read only installed code. Never include source lines or frame locals.
            with path.open('rb') as source:
                content = source.read(1024**2 + 1)
            if len(content) > 1024**2:
                continue
            result.append(dict(module=module, line=line, source_sha256=hashlib.sha256(content).hexdigest()))
        return result, seen > 16 or traceback is not None

    def capture(self, error):
        from genosapiens.wgs_failures import static_failure
        from genosapiens.wgs_runtime import AnalysisProgressRecorder
        if self.primary is not None:
            return self.primary
        self.primary = static_failure(error, 'bounded_products', operation=self.OPERATIONS[self.subphase])
        from genosapiens.population_filter_candidate import (
            OptionalPopulationFilterError, PopulationCandidateBindingError,
        )
        local_code = (OptionalPopulationFilterError.failure_code
            if type(error) is OptionalPopulationFilterError else
            PopulationCandidateBindingError.failure_code
            if type(error) is PopulationCandidateBindingError else None)
        recovery = ('clean_context_evaluation_required'
            if type(error) is OptionalPopulationFilterError else 'check_dependencies_before_science_restart')
        # Current invocation authority only. A stale observation file cannot
        # replace admission detail or establish a checkpoint.
        from genosapiens.wgs_admission_diagnostics import current
        owner = current()
        if self.terminal and owner is not None and not hasattr(owner, '_scientific_failure'):
            owner._scientific_failure = (self.subphase, self.primary, local_code, recovery)
        # Freeze the last scientific observation before resource and iterator cleanup.
        # This file is telemetry. It never seals a scientific product or checkpoint.
        if self.terminal and isinstance(self.progress, AnalysisProgressRecorder):
            try:
                self.progress.finish('interrupted' if self.primary.category == 'interrupted' else 'failed')
            except BaseException:
                pass
        temporary = self.root / 'failure-observation.tmp'
        frames = []
        try:
            frames, truncated = self._frames(error)
            attempt = getattr(self.progress, 'attempt_id', None)
            if type(attempt) is not str or re.fullmatch('[0-9a-f]{64}', attempt) is None:
                attempt = None
            run_digest = self.run_manifest_sha256
            if type(run_digest) is not str or re.fullmatch('[0-9a-f]{64}', run_digest) is None:
                run_digest = None
            ordinal = self.context_ordinal
            if ordinal is not None and (type(ordinal) is not int or not 1 <= ordinal < 2**63):
                ordinal = None
            value = dict(schema='genosapiens.scientific-failure-observation', version=1,
                scientific_completion_verified=False, run_manifest_sha256=run_digest,
                attempt_id=attempt, context_ordinal=ordinal, subphase=self.subphase,
                failure_code=self.primary.failure_code, operation=self.primary.operation,
                local_failure_code=local_code, recovery=recovery,
                exception_kind=self.primary.exception_kind, source_frames=frames,
                source_frames_truncated=truncated)
            content = (json.dumps(value, sort_keys=True, separators=(',', ':'), allow_nan=False)+'\n').encode('ascii')
            if len(content) > 16384:
                return self.primary
            descriptor = os.open(temporary, os.O_WRONLY | os.O_CREAT | os.O_EXCL | os.O_NOFOLLOW, 0o600)
            with os.fdopen(descriptor, 'wb') as stream:
                stream.write(content)
                stream.flush()
                os.fsync(stream.fileno())
            os.replace(temporary, self.root / 'failure-observation.json')
        except BaseException:
            # Failed diagnostics must preserve the original scientific failure.
            pass
        try:
            location = ('unavailable' if not frames else
                frames[-1]['module'] + ':' + str(frames[-1]['line']))
            print('Scientific failure observation: step=' + self.subphase
                + ', operation=' + self.primary.operation + ', category=' + self.primary.exception_kind
                + ', source=' + location + ', recovery=' + recovery + '.', file=sys.stderr, flush=True)
        except BaseException:
            pass
        return self.primary

    def __enter__(self):
        return self

    def __exit__(self, kind, error, traceback):
        if error is not None:
            raise self.capture(error) from None


@dataclass(frozen=True, slots=True, repr=False)
class ScientificContextAnnotation:
    """One positive context, with its complete bounded batch-report root."""
    context_ordinal: int
    context_id: str
    member_count: int
    batch_count: int
    clinical_completed: bool
    context_root: EvidenceReference
    reports_root: EvidenceReference

    def __post_init__(self):
        if (any(type(n) is not int or n <= 0 for n in
                (self.context_ordinal, self.member_count, self.batch_count))
                or self.batch_count > self.member_count
                or type(self.context_id) is not str or not self.context_id
                or type(self.clinical_completed) is not bool
                or type(self.context_root) is not EvidenceReference
                or type(self.reports_root) is not EvidenceReference
                or self.context_root.kind != "context-evidence-v1"
                or self.reports_root.kind != "context-annotation-reports-v1"
                or self.context_root.manifest_sha256 != self.reports_root.manifest_sha256):
            raise EvidenceStoreError("invalid complete context annotation binding")


@dataclass(frozen=True,slots=True,repr=False)
class ScientificArtifact:
    evidence: EvidenceManifest
    case_name: str
    case_sha256: str
    resources_name: str
    resources_sha256: str
    reference_name: str
    reference_sha256: str
    parent_bindings: tuple[tuple[str,str],...]
    context_count: int
    positive_context_count: int
    member_count: int
    pair_count: int
    hypothesis_count: int
    policy_id: str
    policy_sha256: str
    derivation_sha256: str
    schema: str = SCHEMA
    version: int = 2

    def __post_init__(self):
        if (type(self.evidence) is not EvidenceManifest or self.schema != SCHEMA
                or type(self.version) is not int or not (self.version == 2 or
                    type(self) is ScopedScientificArtifact and self.version == 3)
                or self.case_name != "case" or self.resources_name != "case-resources"
                or self.reference_name != "case-reference"
                or type(self.policy_id) is not str or not self.policy_id
                or any(type(n) is not int or not 0 <= n < 2**63 for n in (
                    self.context_count,self.positive_context_count,self.member_count,self.pair_count,self.hypothesis_count))
                or self.positive_context_count > self.context_count
                or type(self.parent_bindings) is not tuple or not self.parent_bindings
                or self.parent_bindings != tuple(sorted(self.parent_bindings))
                or len(dict(self.parent_bindings)) != len(self.parent_bindings)
                or any(type(k) is not str or not k for k,_ in self.parent_bindings)):
            raise EvidenceStoreError("invalid normalized scientific descriptor")
        for digest in (self.case_sha256,self.resources_sha256,self.reference_sha256,self.policy_sha256,self.derivation_sha256,
                       *(v for _,v in self.parent_bindings)):
            _sha(digest)
        self.as_bytes()

    def as_bytes(self):
        content = (_canonical({"classification":"controlled","llm_policy":"no-llm","artifact":self})+"\n").encode("ascii")
        if len(content) > MAX_ARTIFACT_BYTES:
            raise EvidenceStoreError("scientific artifact descriptor exceeds capacity",capacity=True)
        return content

    @property
    def sha256(self):
        return hashlib.sha256(self.as_bytes()).hexdigest()


@dataclass(frozen=True, slots=True, repr=False)
class ScopedScientificArtifact(ScientificArtifact):
    """An explicit evaluated domain. Incomplete coverage cannot authorize export."""
    version: int = 3
    scope_name: str = 'case-scope'
    scope_sha256: str = '0' * 64
    completion_status: str = 'incomplete'

    def __post_init__(self):
        super(ScopedScientificArtifact, self).__post_init__()
        if (self.version != 3 or self.scope_name != 'case-scope'
                or self.completion_status not in {'exhaustive_complete', 'rapid_complete', 'incomplete'}):
            raise EvidenceStoreError('invalid scoped scientific descriptor')
        _sha(self.scope_sha256)


def scientific_artifact_complete(artifact):
    return type(artifact) is ScientificArtifact or (type(artifact) is ScopedScientificArtifact
        and artifact.completion_status in {'rapid_complete', 'exhaustive_complete'})


def scientific_artifact_from_bytes(content, *, expected_sha256):
    try:
        _sha(expected_sha256)
        if (type(content) is not bytes or len(content) > MAX_ARTIFACT_BYTES
                or hashlib.sha256(content).hexdigest() != expected_sha256):
            raise EvidenceStoreError("scientific artifact descriptor identity changed")
        value = _loads(content)
        if (type(value) is not dict or set(value) != {"classification","llm_policy","artifact"}
                or value["classification"] != "controlled" or value["llm_policy"] != "no-llm"):
            raise EvidenceStoreError("invalid scientific artifact envelope")
        kind = ScopedScientificArtifact if value['artifact'].get('version') == 3 else ScientificArtifact
        artifact = _typed(value["artifact"],kind,decode=True)
        if artifact.as_bytes() != content:
            raise EvidenceStoreError("noncanonical scientific artifact descriptor")
        return artifact
    except (ValueError,TypeError,KeyError):
        raise EvidenceStoreError("invalid scientific artifact descriptor") from None


def _write_descriptor(root, content):
    partial,final = root/(DESCRIPTOR_NAME+".partial"),root/DESCRIPTOR_NAME
    fd = os.open(partial,os.O_WRONLY|os.O_CREAT|os.O_EXCL|os.O_NOFOLLOW,0o600)
    with os.fdopen(fd,"wb") as handle:
        handle.write(content)
        handle.flush()
        os.fsync(handle.fileno())
    # No replacement of a concurrent or earlier complete descriptor.
    os.link(partial,final,follow_symlinks=False)
    partial.unlink()
    fd = os.open(root,os.O_RDONLY|os.O_DIRECTORY|os.O_NOFOLLOW)
    try:
        os.fsync(fd)
    finally:
        os.close(fd)


def write_scientific_artifact(contexts, *, output_directory, execution_config,
        run_manifest_sha256, observations, view, profile, gene_evidence, ic_store, expected_ic_sha256,
        prepared_population_policy, hypothesis_policy, hypothesis_policy_sha256, admission_policy_sha256,
        fusion_policy, fusion_policy_sha256, prior_policy, diagnosis, sample_name, source_manifest_sha256,
        parent_bindings, reference_pin, expected_public_identity_id, population_session,
        clinical_annotator=None, query_cache=None, authenticated_resources=None, analysis_eligibility=None,
        analysis_profile=None, progress=None, context_total=None,
        work_filter_policy=None, work_filter_source=None, science_coordinator=None,
        science_workers=1, publication_readers=1, worker_sessions=None):
    """Build from the globally closed, assessed context stream."""
    from genosapiens.wgs_runtime import analysis_progress
    started = completed = hypotheses = pairs = 0
    observation = _ScienceFailureObservation(Path(output_directory), run_manifest_sha256, progress)
    def report(subphase, **values):
        observation.stage(subphase)
        analysis_progress(progress, 'scientific_evaluation', subphase,
            counts={'contexts_started':started, 'contexts_completed':completed, 'contexts_total':context_total,
                'hypotheses':hypotheses, **values}, unit='contexts', completed=completed, total=context_total)
    report('resources')
    if (work_filter_policy is None) != (work_filter_source is None):
        raise EvidenceStoreError('population work policy/source ownership incomplete')
    if work_filter_policy is not None:
        from genosapiens.population_filter_candidate import PopulationWorkPolicy
        if type(work_filter_policy) is not PopulationWorkPolicy:
            raise EvidenceStoreError('canonical population work policy required')
        if parent_bindings.get('population_work_policy') != work_filter_policy.sha256:
            raise EvidenceStoreError('population work policy is not bound to scientific parents')
    root = Path(output_directory)
    if not root.is_absolute() or any(p.is_symlink() for p in (root,*root.parents)):
        raise EvidenceStoreError("invalid scientific artifact directory")
    root.mkdir(mode=0o700,parents=True,exist_ok=True)
    if not root.is_dir() or any(root.iterdir()) or type(execution_config) is not BoundedWgsExecutionConfig:
        raise EvidenceStoreError("scientific artifact output must be empty")
    from genosapiens.wgs_association_projection import verify_assessment_projections
    assessments=verify_assessment_projections(root.parent,parent_bindings,full=False,scratch_directory=root)
    if analysis_profile is not None:
        from genosapiens.wgs_profile_projection import verify_analysis_profile_component, verify_profile_assessment_binding
        verify_analysis_profile_component(root.parent,analysis_profile,full=False,scratch_directory=root)
        verify_profile_assessment_binding(assessments,analysis_profile)
    # The case projection and sort indexes can exceed the runtime's 1 GiB
    # tmpfs. Keep them under this operation's explicit writable disk root.
    reference_metadata = None
    def resources_for(writer, stack):
        nonlocal reference_metadata
        import time
        timings = {}
        stamp = time.monotonic()
        opened = stack.enter_context(_open_submission_reference(reference_pin,
            expected_public_identity_id=expected_public_identity_id))
        reference_metadata = SubmissionReferenceMetadata(writer,reference_pin,opened,name="case-reference")
        timings['reference_open_and_publication'] = time.monotonic()-stamp
        stamp = time.monotonic()
        from genosapiens.wgs_science_resources import prepare_live_gene_evidence
        prepare_genes = (prepare_live_gene_evidence if writer._skip_verification
            else prepare_gene_evidence)
        prepared_genes = prepare_genes(gene_evidence)
        timings['gene_index_preparation'] = time.monotonic()-stamp
        stamp = time.monotonic()
        genes = CaseGeneEvidence(writer,prepared_genes,name="case-genes")
        timings['gene_resource_publication'] = time.monotonic()-stamp
        stamp = time.monotonic()
        phenotypes = CasePhenotypeEvidence(writer,run_manifest_sha256=run_manifest_sha256,observations=observations,
            view=view,profile=profile,ic_store=ic_store,expected_ic_sha256=expected_ic_sha256,
            authenticated_resources=authenticated_resources,name="case-phenotypes")
        timings['phenotype_preparation_and_publication'] = time.monotonic()-stamp
        stamp = time.monotonic()
        population = stack.enter_context(CompletePopulationPolicy(writer,prepared_population_policy,name="case-population"))
        timings['population_preparation'] = time.monotonic()-stamp
        stamp = time.monotonic()
        resources = ScientificCaseResources(writer,genes=genes,phenotypes=phenotypes,population=population,
            hypothesis_policy=hypothesis_policy,hypothesis_policy_sha256=hypothesis_policy_sha256,
            admission_policy_sha256=admission_policy_sha256,fusion_policy=fusion_policy,
            fusion_policy_sha256=fusion_policy_sha256,prior_policy=prior_policy,diagnosis=diagnosis,
            sample_name=sample_name,source_manifest_sha256=source_manifest_sha256,reference=opened.binding,
            parent_bindings=parent_bindings,analysis_eligibility=analysis_eligibility,
            analysis_profile=analysis_profile,name="case-resources")
        timings['resource_binding'] = time.monotonic()-stamp
        resources_for.phase_seconds = timings
        return resources, opened

    if science_coordinator is not None:
        from genosapiens.wgs_science_execution import execute_context_science, _status
        from genosapiens.wgs_science_failures import IndependentScienceFailureObservation
        from genosapiens.wgs_science_runtime import ScienceDeadlineExpired, bounded_operation

        @contextmanager
        def independent_context(writer, resources, opened, closed, *, optional_filter, output_directory, failure_observer):
            if type(closed) is not ClosedContext:
                raise EvidenceStoreError('globally closed scientific context required')
            name = 'context-' + str(closed.ordinal)
            with ExitStack() as local:
                selection = None
                if optional_filter and work_filter_policy is not None:
                    from genosapiens.wgs_population_work_filter import ContextPopulationWorkFilter
                    failure_observer.stage('population_work_filter')
                    selection = local.enter_context(ContextPopulationWorkFilter(writer, closed, policy=work_filter_policy,
                        resources=resources, source_access=work_filter_source, population_session=population_session,
                        clinical_annotator=clinical_annotator, query_cache=query_cache, name=name + '.work_filter'))
                failure_observer.stage('member_enrichment')
                members = local.enter_context(EnrichedContextEvidence(writer, closed, population_session=population_session,
                    execution_config=execution_config, clinical_annotator=clinical_annotator, query_cache=query_cache,
                    name=name + '.members', selection=selection, failure_observer=failure_observer))
                with ScientificCaseContext(writer, resources, members, opened, reference_pin, {}, name=name,
                        progress=progress, failure_observer=failure_observer) as result:
                    yield result

        def observer_for(directory, ordinal):
            observer = IndependentScienceFailureObservation(directory, run_manifest_sha256, progress,
                budget=science_coordinator.budget)
            observer.context_ordinal = ordinal
            return observer

        def initialize_worker(first, stack):
            nonlocal population_session, clinical_annotator, query_cache, work_filter_source, progress
            import sqlite3
            from genosapiens.wgs_science_workers import open_worker_reader
            from genosapiens.wgs_context_payloads import Payloads
            from genosapiens.wgs_contexts import _CLOSED_CONTEXT_SEAL, _load_seed
            from genosapiens.bounded_json import read_bounded_bytes
            import time
            timings = {}
            stamp = time.monotonic()
            # Resource values remain copy-on-write. All used connections and
            # query sessions receive explicit child ownership before evaluation.
            template.attach_worker(stack, lambda: _open_submission_reference(reference_pin,
                expected_public_identity_id=expected_public_identity_id))
            timings['shared_resource_attachment'] = time.monotonic()-stamp
            stamp = time.monotonic()
            reader = open_worker_reader(first._reader, execution_config)
            stack.callback(reader.close)
            first._connection.close()
            connection = sqlite3.connect(context_index_path.as_uri() + '?mode=ro', uri=True)
            connection.execute('PRAGMA query_only=ON')
            stack.callback(connection.close)
            payloads = Payloads(connection)
            payloads.public_closed = getattr(first._payloads, 'public_closed', False)
            timings['private_context_reader'] = time.monotonic()-stamp
            stamp = time.monotonic()
            population_session, clinical_annotator, query_cache, work_filter_source = worker_sessions(stack)
            timings['private_query_sessions'] = time.monotonic()-stamp
            progress = None
            def load(request):
                content = read_bounded_bytes(Path(request['seed_path']),
                    maximum_bytes=execution_config.contexts.max_context_bytes)
                seed = _load_seed(content)
                return ClosedContext(request['ordinal'], seed, connection, reader, execution_config,
                    scratch_directory=science_coordinator.scratch_directory, payloads=payloads,
                    _seal=_CLOSED_CONTEXT_SEAL)
            load.initialization_metrics = timings
            return load

        context_index_path = None
        def prepare_workers(first, stack):
            return initialize_worker(first, stack)

        def prepare_worker_inputs(first):
            nonlocal context_index_path
            from copy import copy
            # Baseline contexts share all public source records. Close their
            # common deduplication table once before opening read-only workers.
            if (not getattr(first._payloads, 'public_closed', False)
                    and first._connection.execute('SELECT 1 FROM baseline_seeds LIMIT 1').fetchone()):
                baseline = copy(first)
                baseline._baseline = True
                for _ in baseline.iter_public_records():
                    pass
            first._connection.commit()
            context_index_path = Path(first._connection.execute('PRAGMA database_list').fetchone()[2])
        prepare_workers.prepare = prepare_worker_inputs

        def publication_worker_initialize(stack):
            from contextlib import nullcontext
            template.attach_worker(stack, lambda: nullcontext(None))
            return template.resources

        from genosapiens.wgs_science_resources import ScientificResourceTemplate
        remaining = science_coordinator.budget.remaining()
        try:
            if remaining <= 0:
                raise ScienceDeadlineExpired('scientific resource preparation deadline expired')
            with bounded_operation(science_coordinator.budget, maximum_seconds=remaining):
                template = ScientificResourceTemplate(resources_for,
                    scratch_directory=science_coordinator.scratch_directory)
        except BaseException:
            _status(root, science_coordinator, 'incomplete_resource_initialization')
            raise
        with template, science_coordinator.checkpoints.using_produced_prefix(template.evidence_prefix):
            result = execute_context_science(contexts, coordinator=science_coordinator, resource_factory=template.bind,
                context_factory=independent_context, output_directory=root, reference_pin=reference_pin,
                eligibility_root=root.parent / 'eligibility', profile_root=root.parent / 'analysis/profile',
                authenticated_resources=authenticated_resources, progress=progress, failure_observer_factory=observer_for,
                resource_preparation_seconds=template.elapsed_seconds, workers=science_workers,
                resource_phase_seconds=template.phase_seconds,
                worker_initialize=None if worker_sessions is None else prepare_workers,
                publication_readers=publication_readers,
                publication_worker_initialize=publication_worker_initialize)
        if result.evidence is None:
            raise ScienceDeadlineExpired('scientific report publication remains incomplete')
        artifact = ScopedScientificArtifact(evidence=result.evidence, **result.case_values,
            reference_name=reference_metadata.name, reference_sha256=reference_metadata.sha256,
            scope_name=result.scope_name, scope_sha256=result.scope_sha256, completion_status=result.status['status'])
        remaining = science_coordinator.budget.remaining(reporting=True)
        try:
            if remaining <= 0:
                raise ScienceDeadlineExpired('scientific descriptor deadline expired')
            with bounded_operation(science_coordinator.budget, maximum_seconds=remaining, reporting=True):
                _write_descriptor(root, artifact.as_bytes())
        except ScienceDeadlineExpired:
            _status(root, science_coordinator, 'incomplete_deadline')
            raise
        return artifact

    with EvidenceWriter(root,scratch_directory=root) as writer,ExitStack() as stack,observation:
        resources, opened = resources_for(writer, stack)
        case = stack.enter_context(CompleteScientificCase(writer,resources,name="case"))
        renderers = {}
        for closed in contexts:
            if type(closed) is not ClosedContext:
                raise EvidenceStoreError("globally closed scientific context required")
            started += 1
            observation.context_ordinal = closed.ordinal
            report('member_enrichment', current_context_ordinal=closed.ordinal, members=closed.member_count)
            name = "context-" + str(closed.ordinal)
            with ExitStack() as context_scope,observation:
                selection = None
                if work_filter_policy is not None:
                    from genosapiens.wgs_population_work_filter import ContextPopulationWorkFilter
                    report('population_work_filter', current_context_ordinal=closed.ordinal)
                    selection = context_scope.enter_context(ContextPopulationWorkFilter(writer, closed,
                        policy=work_filter_policy, resources=resources, source_access=work_filter_source,
                        population_session=population_session, clinical_annotator=clinical_annotator,
                        query_cache=query_cache, name=name+'.work_filter'))
                members = context_scope.enter_context(EnrichedContextEvidence(writer,closed,
                    population_session=population_session,execution_config=execution_config,
                    clinical_annotator=clinical_annotator,query_cache=query_cache,name=name + ".members",
                    selection=selection,failure_observer=observation))
                with ScientificCaseContext(writer,resources,members,opened,reference_pin,renderers,
                        name=name,progress=progress,failure_observer=observation) as context:
                    observation.stage('context_publication')
                    case.add_context(context)
                    completed += 1
                    hypotheses += context.hypothesis_count
                    pairs += context.pair_count
                    report('context_publication')
        report('primary_seal')
        case.publish_primary()
        report('case_seal')
        case.publish(opened,reference_pin,renderers)
        values = dict(case_name=case.name,case_sha256=case.sha256,resources_name=resources.name,
            resources_sha256=resources.sha256,parent_bindings=tuple(sorted(resources.parent_bindings.items())),
            reference_name=reference_metadata.name,reference_sha256=reference_metadata.sha256,
            context_count=case.context_count,positive_context_count=case.positive_context_count,member_count=case.member_count,
            pair_count=case.pair_count,hypothesis_count=case.hypothesis_count,
            policy_id=resources.fusion_policy.policy_id,policy_sha256=resources.fusion_policy_sha256,
            derivation_sha256=resources.fusion_policy.derivation_sha256)
        stack.close()
        report('evidence_finish')
        evidence = writer.finish()
    artifact = ScientificArtifact(evidence=evidence,**values)
    report('descriptor_publication')
    _write_descriptor(root,artifact.as_bytes())
    return artifact


class ScientificArtifactReader:
    def __init__(self, artifact, *, directory, expected_parent_bindings, expected_reference_binding, scratch_directory,
                 eligibility_root=None, profile_root=None, progress=None, authenticated_resources=None, completed_reuse=None):
        self.progress = progress
        self.authenticated_resources = authenticated_resources
        self.completed_reuse, self._selected_cache = completed_reuse, None
        if completed_reuse is not None:
            from genosapiens.wgs_scientific_checkpoint import CompletedScientificReuse
            if type(completed_reuse) is not CompletedScientificReuse:
                raise EvidenceStoreError('compatible completed scientific owner required')
            completed_reuse.require_artifact(artifact)
        if type(artifact) not in (ScientificArtifact, ScopedScientificArtifact) or dict(artifact.parent_bindings) != expected_parent_bindings:
            raise EvidenceStoreError("foreign scientific artifact parents")
        self.artifact,self.directory,self.parents = artifact,Path(directory),dict(expected_parent_bindings)
        self.eligibility_root=self.directory.parent/'eligibility' if eligibility_root is None else Path(eligibility_root)
        self.profile_root=self.directory.parent/'analysis/profile' if profile_root is None else Path(profile_root)
        scratch = Path(scratch_directory)
        if (not scratch.is_absolute() or not scratch.is_dir()
                or any(p.is_symlink() for p in (scratch,*scratch.parents))):
            raise EvidenceStoreError("explicit scientific replay scratch directory required")
        self.expected_reference_binding,self.scratch_directory = expected_reference_binding,scratch
        self.stack,self.closed,self.entered = ExitStack(),False,False

    def __enter__(self):
        if self.closed or self.entered:
            raise EvidenceStoreError("scientific artifact reader is unavailable")
        self.entered = True
        try:
            from genosapiens.wgs_association_projection import verify_assessment_projections
            self.assessment_projections=verify_assessment_projections(self.directory.parent,self.parents,full=False,
                scratch_directory=self.scratch_directory)
            reader = EvidenceReader(self.artifact.evidence,self.directory,
                scratch_directory=self.scratch_directory)
            if self.progress is not None:
                from genosapiens.wgs_runtime import analysis_progress
                original_rows = reader._rows
                nodes = len(self.artifact.evidence.node_shards)
                total = nodes + len(self.artifact.evidence.root_shards)
                analysis_progress(self.progress, 'publication', 'scientific_integrity',
                    unit='files', completed=0, total=total)
                def observed_rows(shard, number):
                    yield from original_rows(shard, number)
                    completed = number + 1 + (0 if shard.kind == 'evidence_nodes' else nodes)
                    if number % 250 == 0 or completed == total:
                        analysis_progress(self.progress, 'publication', 'scientific_integrity',
                            unit='files', completed=completed, total=total)
                reader._rows = observed_rows
            self.reader = self.stack.enter_context(reader)
            with self.reader.read_scope(),self.reader.complete_root_scope():
                if self.completed_reuse is None:
                    self._restore_primary()
                else:
                    self._restore_completed_primary()
                    self.completed_reuse.cover_producer_roots(self.reader, self.artifact)
                if type(self.artifact) is ScopedScientificArtifact:
                    from genosapiens.wgs_science_scope_evidence import verify_scientific_scope
                    if self.case.context_domain is None:
                        raise EvidenceStoreError('scoped scientific artifact lacks an explicit domain')
                    self.scientific_scope = self.stack.enter_context(verify_scientific_scope(self.reader, self.artifact, self.case))
                elif self.case.context_domain is not None:
                    raise EvidenceStoreError('explicit scientific domain requires a scoped descriptor')
            return self
        except BaseException:
            self.close()
            raise

    def _restore_completed_primary(self):
        self.completed_reuse.require_artifact(self.artifact)
        self._restore_primary(completed=True)

    def _restore_primary(self, *, completed=False):
        from genosapiens.wgs_runtime import analysis_progress
        analysis_progress(self.progress, 'publication', 'reference_metadata', force=True)
        reference_ref = bound_reference(self.reader,{"name":self.artifact.reference_name,
            "kind":METADATA_KIND,"sha256":self.artifact.reference_sha256},METADATA_KIND)
        metadata = self.stack.enter_context(VerifiedSubmissionReferenceMetadata(self.reader,reference_ref,
            expected_binding=self.expected_reference_binding))
        self.reference,self.reference_pin,self.renderers = metadata.reference,metadata.pin,metadata.renderers
        resources_ref = bound_reference(self.reader,{"name":self.artifact.resources_name,
            "kind":"scientific-case-resources-v1","sha256":self.artifact.resources_sha256},"scientific-case-resources-v1")
        from genosapiens.wgs_case_resources import CompletedScientificCaseResources
        resource_type = CompletedScientificCaseResources if completed else VerifiedScientificCaseResources
        self.resources = self.stack.enter_context(resource_type(self.reader,resources_ref,
            expected_parent_bindings=self.parents,eligibility_root=self.eligibility_root,profile_root=self.profile_root,
            progress=self.progress,authenticated_resources=self.authenticated_resources))
        if self.resources.analysis_profile is not None:
            from genosapiens.wgs_profile_projection import verify_profile_assessment_binding
            verify_profile_assessment_binding(self.assessment_projections,self.resources.analysis_profile)
        if self.resources.reference_binding != self.expected_reference_binding:
            raise EvidenceStoreError("scientific artifact reference parent changed")
        case_ref = bound_reference(self.reader,{"name":self.artifact.case_name,
            "kind":"complete-scientific-case-v1","sha256":self.artifact.case_sha256},"complete-scientific-case-v1")
        if completed:
            from genosapiens.wgs_case_evidence import CompletedScientificCase
            self.case = self.stack.enter_context(CompletedScientificCase(self.reader, case_ref, self.resources, progress=self.progress))
        else:
            self.case = self.stack.enter_context(VerifiedScientificCase(self.reader,case_ref,self.resources,
                self.reference,self.reference_pin,self.renderers,progress=self.progress))
        if (any(getattr(self.case,k) != getattr(self.artifact,k) for k in
                ("context_count","positive_context_count","member_count","pair_count","hypothesis_count"))
                or self.resources.fusion_policy.policy_id != self.artifact.policy_id
                or self.resources.fusion_policy_sha256 != self.artifact.policy_sha256
                or self.resources.fusion_policy.derivation_sha256 != self.artifact.derivation_sha256):
            raise EvidenceStoreError("scientific artifact metadata does not replay")

    def _require_open(self):
        if self.closed or not self.entered:
            raise EvidenceStoreError("scientific artifact reader is closed")
        self.reader._check_files()

    def ordered(self, mode="primary"):
        self._require_open()
        with self.reader.read_scope():
            yield from self.case.ordered(mode)

    def record_at_rank(self, rank, *, mode="primary"):
        self._require_open()
        if self.completed_reuse is not None:
            if mode != 'primary' or rank not in self.case.selection.selected:
                raise EvidenceStoreError('unselected rich records require explicit diagnostic replay')
            for row in self.iter_selected_records():
                if row['ordering']['rank'] == rank:
                    return row
            raise EvidenceStoreError('selected scientific rank is unavailable')
        with self.reader.read_scope():
            return self.case.record_at_rank(self.reader,rank,self.reference,self.reference_pin,self.renderers,mode=mode)

    def iter_scalar_records(self, *, mode="primary"):
        self._require_open()
        with self.reader.read_scope():
            yield from self.case.scalar_records(mode)

    def credible_complete_count(self, threshold):
        self._require_open()
        with self.reader.read_scope():
            return self.case.credible_complete_count(threshold)

    def iter_selected_records(self):
        """Selected primary detail in context/rank order, with owned saved reuse."""
        self._require_open()
        if self.completed_reuse is not None and self._selected_cache is not None:
            with self.reader.read_scope():
                yield from (dict(row) for row in self._selected_cache)
            return
        with self.reader.read_scope():
            selected_cache, selected_bytes = [], 0
            by_context = {}
            for rank in self.case.selection.selected:
                ordered = self.case.ordered_at(rank)
                by_context.setdefault(ordered["locator"][0],[]).append(ordered)
            count = 0
            restored = direct = 0
            from genosapiens.wgs_runtime import analysis_progress
            if self.completed_reuse is not None:
                analysis_progress(self.progress, 'publication', 'selected_detail', unit='contexts', completed=0, total=len(by_context),
                    counts=dict(reused_scientific_contexts=self.case.context_count, selected_contexts_total=len(by_context),
                        selected_contexts_replayed=0, selected_contexts_read=0,
                        scientific_resources_restored=int(self.resources._loaded is not None) if self.completed_reuse is not None
                            and hasattr(self.resources, '_loaded') else None))
            for ordinal in sorted(by_context):
                direct_format = False
                if self.completed_reuse is not None:
                    context_header = self.reader.reconstruct(bound_reference(self.reader,
                        self.case.context(ordinal)['context'], 'scientific-case-context-v1'))
                    direct_format = context_header['version'] == 2
                if self.completed_reuse is not None and (direct_format or all(row['locator'][2] == 0 for row in by_context[ordinal])):
                    for value in self._selected_records_from_saved(ordinal, by_context[ordinal]):
                        selected_bytes += len(_canonical(value))
                        if selected_bytes > 128 * 1024**2:
                            raise EvidenceStoreError('selected scientific detail exceeds capacity', capacity=True)
                        selected_cache.append(value)
                        count += 1
                    direct += 1
                    analysis_progress(self.progress, 'publication', 'selected_detail', unit='contexts', completed=direct+restored,
                        total=len(by_context), counts=dict(selected_contexts_read=direct, selected_contexts_replayed=restored))
                    continue
                with self.case._restore(self.reader,self.case.context(ordinal),self.reference,
                        self.reference_pin,self.renderers) as context:
                    positions = None if context.universe.pairs is None else context.universe.pairs.domain.positions()
                    for ordered in by_context[ordinal]:
                        record = self.case.logical_record(context,ordered["locator"][1:],context.scores,positions=positions)
                        if (record["constructed"].hypothesis.hypothesis_id != ordered["hypothesis_id"]
                                or record["ledger"].total != ordered["total"] or record["eligibility_reasons"]):
                            raise EvidenceStoreError("selected primary detail does not replay")
                        count += 1
                        value = {**record,"ordering":ordered}
                        if self.completed_reuse is not None:
                            selected_bytes += len(_canonical(value))
                            if selected_bytes > 128 * 1024**2:
                                raise EvidenceStoreError('selected scientific detail exceeds capacity', capacity=True)
                            selected_cache.append(value)
                        else:
                            yield value
                restored += 1
                if self.completed_reuse is not None:
                    analysis_progress(self.progress, 'publication', 'selected_detail', unit='contexts', completed=direct+restored,
                        total=len(by_context), counts=dict(selected_contexts_read=direct, selected_contexts_replayed=restored))
            if count != len(self.case.selection.selected):
                raise EvidenceStoreError("selected primary detail coverage changed")
            if self.completed_reuse is not None:
                self._selected_cache = tuple(selected_cache)
                yield from (dict(row) for row in self._selected_cache)

    def _selected_records_from_saved(self, ordinal, ordered_rows):
        from genosapiens.wgs_completed_science import selected_records, selected_singletons
        context = self.reader.reconstruct(bound_reference(self.reader,
            self.case.context(ordinal)['context'], 'scientific-case-context-v1'))
        if context['version'] == 1:
            yield from selected_singletons(self.reader, self.case, self.resources, ordinal, ordered_rows)
        else:
            yield from selected_records(self.reader, self.case, self.resources, ordinal, ordered_rows,
                reference_names=self.reference.public_names)

    def iter_representation_records(self):
        """Read all already replayed member representations without pair expansion."""
        self._require_open()
        count = 0
        with self.reader.read_scope():
            for row in self.case.contexts():
                header = self.reader.reconstruct(bound_reference(self.reader,row["context"],"scientific-case-context-v1"))
                representations = self.reader.reconstruct(bound_reference(self.reader,header["representations"],
                    "complete-submission-representations-v1"))
                value = representations["representations"]
                root = self.reader.reference(value["name"],kind="member-submission-representations-v1")
                if root.sha256 != value["sha256"]:
                    raise EvidenceStoreError("case member representation binding changed")
                observed = 0
                for value in self.reader.iter_array(root):
                    assessment = _typed(value,SubmissionRepresentationAssessment,decode=True)
                    observed += 1
                    count += 1
                    yield {"context_ordinal":row["ordinal"],"context_id":row["context_id"],"assessment":assessment}
                if observed != row["member_count"]:
                    raise EvidenceStoreError("case representation coverage changed")
            if count != self.artifact.member_count:
                raise EvidenceStoreError("case representation count changed")

    def iter_member_records(self):
        """Complete enriched source members; no genotype/pair/score recomputation."""
        self._require_open()
        count = 0
        with self.reader.read_scope():
            resources = self.resources
            if ('population_work_policy' in resources.parent_bindings
                    and hasattr(resources, 'restore_for_selected')):
                resources = resources.restore_for_selected()
            for row in self.case.contexts():
                header = self.reader.reconstruct(bound_reference(self.reader,row["context"],"scientific-case-context-v1"))
                root = bound_reference(self.reader,header["members"],"context-evidence-v1")
                with VerifiedContextEvidence(self.reader,root,resources=resources) as members,closing(members.iter_members()) as values:
                    if (members.ordinal != row["ordinal"] or members.context_id != row["context_id"]
                            or members.member_count != row["member_count"]):
                        raise EvidenceStoreError("case source member binding changed")
                    for member in values:
                        count += 1
                        yield {"context_ordinal":members.ordinal,"context_id":members.context_id,
                            "seed":members.seed,"member":member}
            if count != self.artifact.member_count:
                raise EvidenceStoreError("case source member count changed")

    def annotation_contexts(self):
        """One root binding per positive context, independent of batch size."""
        self._require_open()
        count = 0
        with self.reader.read_scope():
            for row in self.case.contexts():
                if not row["member_count"]:
                    continue
                header = self.reader.reconstruct(bound_reference(self.reader, row["context"], "scientific-case-context-v1"))
                root = bound_reference(self.reader, header["members"], "context-evidence-v1")
                members = self.reader.reconstruct(root)
                reports = self.reader.reference(members["reports"]["name"], kind="context-annotation-reports-v1")
                if (reports.sha256 != members["reports"]["sha256"]
                        or members["ordinal"] != row["ordinal"] or members["member_count"] != row["member_count"]
                        or members["seed"]["context"]["context_id"] != row["context_id"]):
                    raise EvidenceStoreError("case annotation context binding changed")
                count += 1
                yield ScientificContextAnnotation(row["ordinal"], row["context_id"], row["member_count"],
                    members["batch_count"], members["has_clinical"], root, reports)
            if count != self.artifact.positive_context_count:
                raise EvidenceStoreError("case annotation context coverage changed")

    def annotation_reports(self, annotation):
        """Read every bounded report for an authenticated annotation context."""
        from genosapiens.wgs_context_evidence import MAX_RECORD_BYTES
        self._require_open()
        if type(annotation) is not ScientificContextAnnotation:
            raise EvidenceStoreError("complete context annotation binding required")
        with self.reader.read_scope():
            row = self.case.context(annotation.context_ordinal)
            header = self.reader.reconstruct(bound_reference(self.reader, row["context"], "scientific-case-context-v1"))
            root = bound_reference(self.reader, header["members"], "context-evidence-v1")
            members = self.reader.reconstruct(root)
            reports = self.reader.reference(members["reports"]["name"], kind="context-annotation-reports-v1")
            if annotation != ScientificContextAnnotation(row["ordinal"], row["context_id"], row["member_count"],
                    members["batch_count"], members["has_clinical"], root, reports):
                raise EvidenceStoreError("foreign context annotation reports")
            count = 0
            for value in self.reader.iter_array(reports, max_item_bytes=MAX_RECORD_BYTES + 1024):
                count += 1
                if value["ordinal"] != count:
                    raise EvidenceStoreError("context annotation report order changed")
                yield value
            if count != annotation.batch_count:
                raise EvidenceStoreError("context annotation report coverage changed")

    def iter_scientific_records(self, *, mode="primary", include_ordering=False):
        """Explicit complete logical replay in context/ID order.

        Each context is restored once, then every occurrence is reconstructed
        from its computed tapes. This interface is for logical consumption;
        routine publication retains the normalized evidence and complete
        orders instead of writing rich evidence again for every hypothesis.
        """
        self._require_open()
        if self.completed_reuse is not None:
            raise EvidenceStoreError('full scientific replay requires an explicit independent reader')
        if mode != "primary" or type(include_ordering) is not bool:
            raise EvidenceStoreError("invalid scientific replay mode")
        count = 0
        with self.reader.read_scope():
            for row in self.case.contexts():
                with self.case._restore(self.reader,row,self.reference,self.reference_pin,self.renderers,
                        retain_identity_index=True) as context,ExitStack() as stack:
                    scores = self.case.context_scores(self.reader,context,mode,stack)
                    positions = None if context.universe.pairs is None else context.universe.pairs.domain.positions()
                    with closing(context.order.identity_locators()) as locators:
                        for hypothesis_id,locator in locators:
                            record = self.case.logical_record(context,locator,scores,positions=positions)
                            if record["constructed"].hypothesis.hypothesis_id != hypothesis_id:
                                raise EvidenceStoreError("logical scientific identity changed")
                            if include_ordering:
                                total = record["ledger"].total
                                if mode != "primary":
                                    total = scores.ordering_total(total)
                                ties = tuple(key.value for key in record["ledger"].tie_keys[:-1])
                                rank = self.case.rank_of((context.ordinal, *locator), total=total, ties=ties, mode=mode)
                                record = {**record, "ordering": self.case.ordered_at(rank, mode)}
                            count += 1
                            yield record
            if count != self.artifact.hypothesis_count:
                raise EvidenceStoreError("logical scientific replay coverage changed")
    def selection_audit(self):
        self._require_open()
        with self.reader.read_scope():
            yield from self.case.selection_audit()

    @property
    def selected_groups(self):
        self._require_open()
        return self.case.selected_groups

    def close(self):
        try:
            self.stack.close()
        finally:
            self._selected_cache = None
            self.closed = True

    def __exit__(self,*_):
        self.close()
