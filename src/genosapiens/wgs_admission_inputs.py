"""Read-only current-format admission input for public discovery continuations.

The descriptor binds an independently completed stage, not a reusable database
path or a serialized authentication result. Consumers reuse producer completion
through bounded metadata checks; requested records retain their strict codecs.
"""
from __future__ import annotations

from contextlib import closing
from dataclasses import asdict, dataclass
import hashlib
import json
import os
from pathlib import Path
import sqlite3

from genosapiens.evidence_model import _typed
from genosapiens.wgs_admission_projection import AdmissionProjectionCompletion
from genosapiens.wgs_analysis_inputs import AnalysisInputFile, _path, _safe_file
from genosapiens.wgs_ledger import (
    Artifact, JournalEntry, StageRunBinding, canonical, load_model,
)
from genosapiens.wgs_shards import read_bounded_metadata

SCHEMA = 'genosapiens.wgs-completed-admission-discovery-input'
MAX_DESCRIPTOR_BYTES = 64 * 1024
CONFIG_PATH = Path(__file__).resolve().parents[2] / 'config/track1/bounded-wgs-execution-v5.json'
_RESOLVED = object()
_FILES = frozenset((
    'producer_authority', 'producer_release', 'producer_plan', 'producer_start',
    'producer_terminal', 'runtime', 'runtime_sample', 'runtime_phases', 'failure',
    'completed_analysis_input', 'completed_annotation_prefix',
    'admission_head', 'admission_completion', 'screening_context',
    'eligibility_manifest', 'rescue_manifest', 'annotation_manifest',
    'annotation_closure', 'annotation_head', 'annotation_journal', 'annotation_prefix',
))


class AdmissionInputError(RuntimeError):
    """Static rejection of a missing, changed or incomplete admission input."""


def _require(value):
    if not value:
        raise AdmissionInputError('completed admission input is invalid')


def _sha(content):
    return hashlib.sha256(content).hexdigest()


def _digest(value):
    return type(value) is str and len(value) == 64 and all(c in '0123456789abcdef' for c in value)


def _unique(pairs):
    result = {}
    for key, value in pairs:
        _require(key not in result)
        result[key] = value
    return result


def _json(content):
    value = json.loads(content, object_pairs_hook=_unique)
    _require(type(value) is dict and canonical(value) == content)
    return value


def _signature(path):
    s = Path(path).stat(follow_symlinks=False)
    return s.st_dev, s.st_ino, s.st_size, s.st_mtime_ns, s.st_ctime_ns


@dataclass(frozen=True, slots=True)
class AdmissionDatabase:
    path: str
    size_bytes: int
    user_version: int

    def __post_init__(self):
        _path(self.path)
        _require(type(self.size_bytes) is int and self.size_bytes > 0
                 and type(self.user_version) is int and self.user_version in (1, 5))

    def reader_path(self):
        """Validate the real reader's file and parent contract without opening SQL."""
        path = _safe_file(self.path, self.size_bytes)
        if self.user_version == 5:
            # Preflight must enforce the same path contract as the live owner.
            # A file bind alone can leave a root-owned synthetic parent in Docker.
            from genosapiens.wgs_admission_projection import _private_path
            path = _private_path(path, exists=True)
        return path

    def open(self):
        path = self.reader_path()
        db = sqlite3.connect(path.as_uri() + '?mode=ro&immutable=1', uri=True)
        try:
            db.execute('PRAGMA query_only=ON')
            db.execute('PRAGMA cache_size=-8192')
            _require(db.execute('PRAGMA user_version').fetchone() == (self.user_version,))
            return db
        except BaseException:
            db.close()
            raise


@dataclass(frozen=True, slots=True, repr=False)
class CompletedAdmissionInput:
    files: tuple[tuple[str, AnalysisInputFile], ...]
    database: AdmissionDatabase
    stage_state: AdmissionDatabase
    annotation_state: AdmissionDatabase
    original_binding: str
    binding_sha256: str
    stage_head_sha256: str
    producer_executable_sha256: str
    runtime_image_digest: str
    execution_config_sha256: str
    input_dispositions: tuple[tuple[str, int], ...]
    occurrence_dispositions: tuple[tuple[str, int], ...]
    expected_completion: AdmissionProjectionCompletion
    native_index: Artifact | None = None
    schema: str = SCHEMA
    version: int = 1
    classification: str = 'public'

    def __post_init__(self):
        from genosapiens.wgs_admission_projection import _binding
        _require(self.schema == SCHEMA and type(self.version) is int and self.version == 1
                 and self.classification == 'public' and type(self.files) is tuple
                 and all(type(row) is tuple and len(row) == 2 and type(row[0]) is str
                         and type(row[1]) is AnalysisInputFile for row in self.files)
                 and tuple(key for key, _ in self.files) == tuple(sorted(_FILES)))
        _require(all(type(v) is AdmissionDatabase for v in
                     (self.database, self.stage_state, self.annotation_state))
                 and self.database.user_version == 5 and self.stage_state.user_version == 1
                 and self.annotation_state.user_version == 1)
        _require(type(self.original_binding) is str)
        binding = self.original_binding.encode('ascii')
        parents = _binding(binding)
        _require(_sha(binding) == self.binding_sha256
                 and all(_digest(v) for v in (self.binding_sha256, self.stage_head_sha256,
                         self.producer_executable_sha256, self.execution_config_sha256))
                 and self.producer_executable_sha256 == parents['consumer_executable_sha256']
                 and self.execution_config_sha256 == parents['execution_limits_sha256']
                 and type(self.runtime_image_digest) is str
                 and self.runtime_image_digest.startswith('sha256:')
                 and _digest(self.runtime_image_digest[7:])
                 and parents['read_intake_sha256'] is None
                 and parents['read_intake_source_bindings_sha256'] is None)
        for counts, total in ((self.input_dispositions, parents['identities']),
                              (self.occurrence_dispositions, parents['occurrences'])):
            _require(type(counts) is tuple and all(type(v) is tuple and len(v) == 2
                     and type(v[1]) is int and v[1] >= 0 for v in counts)
                     and tuple(k for k, _ in counts) == ('assess', 'defer_refcall', 'exclude_common')
                     and sum(v for _, v in counts) == total)
        completion = self.expected_completion
        _require(type(completion) is AdmissionProjectionCompletion
                 and completion.binding_sha256 == self.binding_sha256
                 and completion.sha256 == self.file('admission_completion').sha256
                 and (completion.rows, completion.occurrences, completion.requests) ==
                 (parents['identities'], parents['occurrences'], parents['sequence_requests'])
                 and completion.eligibility_manifest_sha256 == self.file('eligibility_manifest').sha256)
        work = Path(self.database.path).parent.parent
        _require(Path(self.database.path) == work / 'state/admission.sqlite'
                 and Path(self.stage_state.path) == work / 'state/admission-stage/state.sqlite'
                 and Path(self.annotation_state.path) == work / 'state/annotation.sqlite')
        for role, relative in (
            ('admission_head', 'state/admission-stage/admission_projection-head.json'),
            ('admission_completion', 'state/admission-stage/admission-completion.json'),
            ('screening_context', 'state/admission-stage/screening-context.json'),
            ('eligibility_manifest', 'eligibility/manifest.json'),
            ('rescue_manifest', 'eligibility/rescue-manifest.json'),
            ('annotation_manifest', 'annotation/annotation-manifest.json'),
            ('annotation_closure', 'annotation/annotation-stage-closure.json'),
            ('annotation_head', 'annotation/annotation-head.json'),
            ('annotation_journal', 'annotation/annotation-000001.journal.json'),
            ('annotation_prefix', 'annotation/annotation-prefix.json'),
        ):
            _require(Path(self.file(role).path) == work / relative)
        _require(len(self.as_bytes()) <= MAX_DESCRIPTOR_BYTES)
        _require(self.native_index is None or type(self.native_index) is Artifact
                 and self.native_index.name == 'native-index.jsonl')

    def file(self, role):
        return dict(self.files)[role]

    @property
    def work_root(self):
        return Path(self.database.path).parent.parent

    def as_bytes(self):
        return canonical(asdict(self))

    @property
    def sha256(self):
        return _sha(self.as_bytes())

    def readonly_paths(self):
        return tuple(dict.fromkeys((self.database.path, self.stage_state.path,
            self.annotation_state.path, str(self.work_root / 'state/admission-stage'),
            str(self.work_root / 'annotation'), str(self.work_root / 'eligibility'),
            *(value.path for _, value in self.files))))


def load_completed_admission_input(content, *, expected_sha256):
    try:
        _require(type(content) is bytes and 0 < len(content) <= MAX_DESCRIPTOR_BYTES
                 and _digest(expected_sha256) and _sha(content) == expected_sha256)
        data = _json(content)
        _require(set(data) == set(CompletedAdmissionInput.__dataclass_fields__))
        value = _typed(data, CompletedAdmissionInput, decode=True)
        _require(value.as_bytes() == content)
        return value
    except (KeyError, TypeError, ValueError, UnicodeError):
        raise AdmissionInputError('completed admission input is invalid') from None


def read_completed_admission_stage(descriptor):
    """Reconcile bounded completion metadata, stage SQLite and admission singleton.

    No StageCommitState construction, lock, checkpoint, DDL or rich-table scan.
    The caller must independently authenticate the producer and stage parents.
    """
    from genosapiens.wgs_admission_projection import load_admission_projection_completion
    _require(type(descriptor) is CompletedAdmissionInput)
    root = descriptor.work_root / 'state/admission-stage'
    head = _json(descriptor.file('admission_head').read())
    terminal = load_model(read_bounded_metadata(root/'admission_projection-000003.journal.json'), JournalEntry)
    _require(terminal.sha256 == descriptor.stage_head_sha256)
    # The stage has exactly three bounded control events. Resolve those saved
    # controls without replaying a journal or authenticating its payloads again.
    with closing(descriptor.stage_state.open()) as db:
        saved_events=db.execute('SELECT event_bytes FROM events ORDER BY ordinal LIMIT 4').fetchall()
    entries=tuple(load_model(value,JournalEntry) for value, in saved_events)
    _require(len(entries) == 3 and tuple(e.event_type for e in entries) ==
             ('attempt_started', 'partition_completed', 'stage_completed'))
    _require(entries[-1] == terminal)
    for ordinal,event in enumerate(entries,1):
        _require(event.ordinal == ordinal and event.stage == 'admission_projection'
            and event.stage_parent_sha256 == descriptor.binding_sha256
            and event.previous_entry_sha256 == (None if ordinal == 1 else entries[ordinal-2].sha256))
    first, committed, terminal = entries
    _require(first.partition_id == committed.partition_id == 'admission'
             and all(e.scientific is False for e in entries)
             and first.attempt_number == committed.attempt_number == 1
             and first.attempt_id == committed.attempt_id
             and first.process_id == committed.process_id == terminal.process_id
             and all(e.parents == first.parents for e in entries)
             and len(committed.artifacts) == 1 and committed.artifacts == terminal.artifacts
             and committed.artifacts[0].name == 'admission-completion.json')
    # The descriptor anchors the completed terminal; its pointer must agree.
    _require(head == dict(schema='genosapiens.wgs-stage-journal-head', version=1,
                         stage='admission_projection', ordinal=3, sha256=terminal.sha256))
    file = descriptor.file('admission_completion')
    _require(committed.artifacts == (Artifact('admission-completion.json', file.sha256, file.size_bytes),))
    completion = load_admission_projection_completion(file.read(), expected_sha256=file.sha256)
    _require(completion == descriptor.expected_completion)
    with closing(descriptor.stage_state.open()) as db:
        saved = db.execute('SELECT canonical_bytes FROM run_binding LIMIT 2').fetchall()
        _require(len(saved) == 1)
        stage = load_model(saved[0][0], StageRunBinding)
        _require(stage.stage == 'admission_projection'
                 and stage.stage_parent_sha256 == descriptor.binding_sha256
                 and stage.executable_stack_sha256 == descriptor.producer_executable_sha256
                 and stage.runtime_image_digest == descriptor.runtime_image_digest
                 and stage.execution_config_sha256 == descriptor.execution_config_sha256
                 and stage.parents == first.parents)
        events = db.execute('SELECT ordinal,event_bytes,event_sha256,event_type,partition_id,attempt_id '
                            'FROM events ORDER BY ordinal LIMIT 4').fetchall()
        _require(events == [(e.ordinal, e.as_bytes(), e.sha256, e.event_type, e.partition_id, e.attempt_id)
                            for e in entries])
        counts = dict(rows=completion.rows, occurrences=completion.occurrences,
                      proofs=completion.proofs, requests=completion.requests)
        commits = db.execute('SELECT ordinal,partition_id,cursor_bytes,outputs_bytes,counts_bytes '
                             'FROM commits LIMIT 2').fetchall()
        _require(commits == [(2, 'admission', canonical((completion.rows,)),
                             canonical(('admission-completion.json',)), canonical(counts))])
        outputs = db.execute('SELECT name,commit_ordinal,metadata_bytes FROM committed_outputs LIMIT 2').fetchall()
        _require(len(outputs) == 1 and outputs[0][:2] == ('admission-completion.json', 2))
        metadata = _json(outputs[0][2])
        _require(set(metadata) == {'name', 'sha256', 'size_bytes', 'publication'}
                 and {k: metadata[k] for k in ('name', 'sha256', 'size_bytes')} == asdict(committed.artifacts[0]))
        directory = root.stat()
        _require(metadata['publication'] == [directory.st_dev, directory.st_ino,
                                            *_signature(file.path)])
    with closing(descriptor.database.open()) as db:
        _require(db.execute('SELECT canonical_bytes FROM run_binding LIMIT 2').fetchall() ==
                 [(descriptor.original_binding.encode('ascii'),)]
                 and db.execute('SELECT canonical_bytes FROM completion LIMIT 2').fetchall() == [(completion.as_bytes(),)]
                 and db.execute('SELECT canonical_bytes FROM eligibility LIMIT 2').fetchall() ==
                 [(descriptor.file('eligibility_manifest').read(),)])
    return completion, stage


def _producer(descriptor, execution_config):
    """Authenticate historical controls without interpreting failure as success."""
    from genosapiens.bounded_wgs_config import bounded_wgs_execution_config_as_json
    from genosapiens.track1_release import load_base_release
    from genosapiens.wgs_public_qualification import (
        load_public_qualification_authority, require_public_qualification_release,
    )
    from genosapiens.wgs_runtime import (
        load_stage_runtime_evidence, load_stage_runtime_sample, load_runtime_phase_checkpoint,
    )
    file = descriptor.file
    authority = load_public_qualification_authority(file('producer_authority').read(),
        expected_sha256=file('producer_authority').sha256)
    release = load_base_release(file('producer_release').read(),
        expected_sha256=file('producer_release').sha256)
    require_public_qualification_release(authority, release)
    _require(authority.invocation_id == 'complete-hg002-deepvariant-wgs'
             and authority.preparation_attempts == 0
             and authority.executable_stack_sha256 == descriptor.producer_executable_sha256
             and authority.runtime_image_digest == descriptor.runtime_image_digest
             and authority.execution_config_sha256 == descriptor.execution_config_sha256
             == _sha(bounded_wgs_execution_config_as_json(execution_config))
             and descriptor.work_root == Path(authority.output_root) / 'analysis/analysis-work'
             and authority.plan_sha256 == file('producer_plan').sha256)
    start, terminal, plan = (_json(file(role).read()) for role in
                             ('producer_start', 'producer_terminal', 'producer_plan'))
    _require(start['schema'] == 'genosapiens.track1-public-qualification-attempt-start'
             and type(start['version']) is int and start['version'] == 3
             and terminal['schema'] == 'genosapiens.track1-public-qualification-stage-result'
             and type(terminal['version']) is int and terminal['version'] == 3
             and plan['schema'] == 'genosapiens.track1-public-postannotation-continuation-plan'
             and type(plan['version']) is int and plan['version'] == 3)
    for control in (start, terminal, plan):
        _require(control['candidate_sha256'] == authority.engineering_candidate_sha256
                 and control['classification'] == 'public'
                 and control['controlled_use_authorized'] is False
                 and control['controlled_input_accessed'] is False)
    _require(start['authority_sha256'] == authority.sha256
             and start['candidate_release_sha256'] == release.sha256
             and start['executable_stack_sha256'] == plan['executable_stack_sha256']
             == authority.executable_stack_sha256
             and start['execution_config_sha256'] == plan['execution_config']['sha256']
             == authority.execution_config_sha256
             and start['plan_sha256'] == terminal['plan_sha256'] == authority.plan_sha256
             and start['stage'] == terminal['stage'] == 'analysis'
             and start['container_id'] == terminal['container_id']
             and start['automatic_retries'] == start['resume_attempts'] == 0
             and start['attempt'] == start['attempt_budget'] == 1
             and start['admission_projection_invocations'] == 1
             and start['reused_admission_rows'] == 0
             and start['preparation_invocations'] == start['annotation_child_invocations'] == 0
             and terminal['status'] == 'failed' and terminal['exit_code'] == 1
             and terminal['oom_killed'] is False and terminal['runtime_finalization_returncode'] == 0
             and plan['runtime_image_digest'] == authority.runtime_image_digest
             and Path(plan['output_namespace']) / 'full-wgs' == Path(authority.output_root))
    budget = plan['invocation_budget']
    _require(budget['preparation_attempts'] == budget['automatic_retries'] == budget['resume_attempts'] == 0
             and budget['analysis_attempts'] == budget['admission_projection_invocations'] == 1
             and budget['reused_admission_rows'] == 0
             and budget['source_ingestion_invocations'] == budget['source_audit_invocations'] == 0
             and budget['prefilter_invocations'] == budget['annotation_child_invocations'] == 0
             and budget['annotation_plan_constructors'] == budget['retained_projection_constructors'] == 0)
    for role in ('completed_analysis_input', 'completed_annotation_prefix'):
        _require(start[role] == plan[role] == asdict(file(role))
                 and authority.binding(role) == file(role).sha256)
    command = plan['exact_stage_command_template']
    argv, mounts = command['create_argv_template'], command['mounts']
    _require(command['stage'] == 'analysis' and command['attempt_budget'] == 1
             and command['automatic_retries'] == 0 and command['attempt_id'] == start['runtime_attempt_id']
             and type(argv) is list and all(type(v) is str for v in argv)
             and any(argv[i:i + 4] == ['python', '-m', 'genosapiens.wgs_public_runner', 'analyze']
                     for i in range(len(argv) - 3))
             and type(mounts) is list and all(type(m) is dict and set(m) == {'source', 'target', 'read_only'}
                 and type(m['read_only']) is bool for m in mounts)
             and [m for m in mounts if not m['read_only']] ==
                 [dict(source=authority.output_root, target=authority.output_root, read_only=False)])
    for role in ('completed_analysis_input', 'completed_annotation_prefix'):
        for flag, value in (('--' + role.replace('_', '-'), file(role).path),
                            ('--' + role.replace('_', '-') + '-sha256', file(role).sha256)):
            _require(argv.count(flag) == 1 and argv.index(flag) + 1 < len(argv)
                     and argv[argv.index(flag) + 1] == value)
        _require(any(m['read_only'] and m['source'] == m['target']
                     and Path(file(role).path).is_relative_to(m['target']) for m in mounts))
    runtime_arguments = dict(expected_stage='analysis', expected_attempt_id=start['runtime_attempt_id'],
                             expected_execution_config_sha256=authority.execution_config_sha256)
    evidence = load_stage_runtime_evidence(file('runtime').read(), **runtime_arguments)
    sample = load_stage_runtime_sample(file('runtime_sample').read(), expected_exit_code=1, **runtime_arguments)
    phases = load_runtime_phase_checkpoint(file('runtime_phases').read(), **runtime_arguments)
    _require(evidence.exit_code == 1 and evidence.oom_killed is False
             and evidence.sample_sha256 == file('runtime_sample').sha256
             and evidence.metrics == sample.metrics
             and evidence.metrics.phase_checkpoint_sha256 == file('runtime_phases').sha256
             and evidence.metrics.phases == phases.phases and phases.terminal and not phases.active_phases)
    failure = _json(file('failure').read())
    _require(failure['schema'] == 'genosapiens.wgs-failure-evidence' and failure['version'] == 2
             and failure['stage'] == 'analysis' and failure['exception_kind'] == 'project_failure'
             and failure['operation'] == 'bounded_product_validate'
             and failure['failure_code'].startswith('bounded_products_')
             and not failure['secondary_failures'] and failure['secondary_failures_truncated'] is False)
    return authority, release


@dataclass(frozen=True, slots=True, repr=False)
class ResolvedCompletedAdmissionInput:
    descriptor: CompletedAdmissionInput
    producer: object
    producer_release: object
    completion: AdmissionProjectionCompletion
    eligibility_manifest: object
    annotation_source: object
    consumer_authority: object
    consumer_release: object
    signatures: tuple[tuple[str, tuple[int, ...]], ...]
    _proof: object

    def __post_init__(self):
        _require(self._proof is _RESOLVED)

    def require_stable(self):
        for path, signature in self.signatures:
            _require(_signature(path) == signature)

    def require_consumer(self, authority, release, *, completed_analysis_input, annotation_prefix):
        from genosapiens.wgs_analysis_inputs import require_upstream_dependencies
        from genosapiens.wgs_public_qualification import PublicQualificationAuthority
        from genosapiens.track1_release import BaseRelease
        self.require_stable()
        _require(type(authority) is PublicQualificationAuthority and type(release) is BaseRelease
                 and authority == self.consumer_authority and release == self.consumer_release
                 and authority.binding('completed_admission_input') == self.descriptor.sha256
                 and completed_analysis_input.descriptor.sha256 == self.descriptor.file('completed_analysis_input').sha256
                 and annotation_prefix.descriptor.sha256 == self.descriptor.file('completed_annotation_prefix').sha256)
        completed_analysis_input.require_consumer(authority, release)
        annotation_prefix.require_consumer(authority, release)
        require_upstream_dependencies(self.producer, self.producer_release, authority, release, annotation=True)
        for role in ('analysis_eligibility_policy', 'completed_annotation_prefix'):
            _require(self.producer.binding(role) == authority.binding(role))
        old, new = Path(self.producer.output_root), Path(authority.output_root)
        _require(not new.is_relative_to(old) and not old.is_relative_to(new))

    def binding(self):
        return dict(descriptor_sha256=self.descriptor.sha256,
                    admission_binding_sha256=self.descriptor.binding_sha256,
                    completion_sha256=self.completion.sha256,
                    stage_head_sha256=self.descriptor.stage_head_sha256,
                    producer_authority_sha256=self.producer.sha256,
                    producer_release_sha256=self.producer_release.sha256,
                    producer_executable_sha256=self.producer.executable_stack_sha256)

    def profile_facts_provider(self, *, admission_owner, completed_analysis_input,
                               annotation_prefix, clinical_session, regulatory_index, reference):
        """Borrow current original-input/annotation owners without changing producers."""
        from genosapiens.wgs_admission_projection import AuthenticatedAdmissionOwner
        from genosapiens.wgs_profile_projection import AdmissionProfileFactsProvider
        self.require_consumer(self.consumer_authority, self.consumer_release,
            completed_analysis_input=completed_analysis_input, annotation_prefix=annotation_prefix)
        _require(type(admission_owner) is AuthenticatedAdmissionOwner
                 and admission_owner.binding_sha256 == self.descriptor.binding_sha256
                 and admission_owner.completion == self.completion)
        prepared = completed_analysis_input.prepared
        return AdmissionProfileFactsProvider(admission_owner=admission_owner,
            annotation_source=self.annotation_source,
            source_manifest=prepared.completed_source.manifest,
            working_manifest=completed_analysis_input.working_manifest,
            clinical_session=clinical_session, regulatory_index=regulatory_index, reference=reference,
            prepared_index=prepared.prepared_index, prepared_state_root=prepared.root.parent/'state')


def _authenticate_metadata(descriptor, *, completed_analysis_input, annotation_prefix,
                           execution_config, source_routing_context, _new_control_root=None):
    """Authenticate producer metadata with the already resolved C/P dependencies."""
    from genosapiens.wgs_admission_projection import _binding, _require_screening_context
    from genosapiens.wgs_analysis_eligibility import (
        load_admission_screening_context, load_analysis_eligibility_manifest,
        analysis_eligibility_component,
    )
    from genosapiens.wgs_annotation import (
        BoundedAnnotationSource, AnnotationStageResult, _load_stage_manifest,
    )
    from genosapiens.wgs_annotation_prefix import annotation_resource_session_sha256
    from genosapiens.track1_public_validation import VEP_CONTRACT_PATH
    from genosapiens.wgs_analysis_inputs import require_upstream_dependencies
    from genosapiens.wgs_analysis_inputs import ResolvedCompletedAnalysisInput
    from genosapiens.wgs_annotation_prefix import ResolvedAnnotationPrefix
    from genosapiens.wgs_source_routing import SourceRoutingContext
    from genosapiens.vep_annotation import _readonly
    _require(type(descriptor) is CompletedAdmissionInput
             and type(completed_analysis_input) is ResolvedCompletedAnalysisInput
             and type(annotation_prefix) is ResolvedAnnotationPrefix
             and type(source_routing_context) is SourceRoutingContext)
    descriptor.database.reader_path()
    for path in descriptor.readonly_paths():
        # Only newly written generator control copies are writable. Scientific
        # ancestors and all runtime inputs always require read-only mounts.
        if _new_control_root is None or Path(path).parent != _new_control_root:
            _readonly(Path(path))
    signatures = tuple((path, _signature(path)) for path in dict.fromkeys((
        *(v.path for _, v in descriptor.files), descriptor.database.path,
        descriptor.stage_state.path, descriptor.annotation_state.path,
        *(str(descriptor.work_root / f'state/admission-stage/admission_projection-{i:06d}.journal.json')
          for i in range(1, 4)),
        *((str(descriptor.work_root / 'annotation' / descriptor.native_index.name),)
          if descriptor.native_index is not None else ()))))
    # Read every declared small file, even when its downstream reader only
    # needs a subset. This also enforces the private file/path contract.
    for _, file in descriptor.files:
        file.read()
    producer, release = _producer(descriptor, execution_config)
    require_upstream_dependencies(producer, release, completed_analysis_input.consumer_authority,
                                  completed_analysis_input.consumer_release, annotation=True)
    c, p = completed_analysis_input.descriptor, annotation_prefix.descriptor
    _require(c.as_bytes() == descriptor.file('completed_analysis_input').read()
             and p.as_bytes() == descriptor.file('completed_annotation_prefix').read()
             and p.as_bytes() == descriptor.file('annotation_prefix').read())
    annotation_prefix.require_consumer(producer, release)
    completion, stage = read_completed_admission_stage(descriptor)
    parents = _binding(descriptor.original_binding.encode('ascii'))
    screening = load_admission_screening_context(descriptor.file('screening_context').read())
    _require_screening_context(screening, parents)
    _require(screening.merged_intake is None
             and screening.policy.sha256 == producer.binding('analysis_eligibility_policy'))
    prepared = completed_analysis_input.prepared
    from genosapiens.wgs_admission import _source_header
    from genosapiens.wgs_ledger import read_artifact
    source = prepared.completed_source.manifest
    source_root = prepared.completed_source.root
    header_document = json.loads(read_artifact(source_root, source.header))
    _require(screening.header == _source_header(source_root, source)
             and screening.selected_sample == header_document['selected_sample_name'])
    plan = completed_analysis_input.transport.plan
    contract = read_bounded_metadata(Path(__file__).resolve().parents[2] / VEP_CONTRACT_PATH)
    resource = annotation_resource_session_sha256(reference_manifest_sha256=release.reference_manifest_sha256,
        vep_bundle_manifest_sha256=producer.binding('vep_bundle'), vep_contract_sha256=_sha(contract),
        executable_stack_sha256=producer.executable_stack_sha256, runtime_image_digest=producer.runtime_image_digest)
    manifest = _load_stage_manifest(descriptor.file('annotation_manifest').read())
    expected = dict(prepared_index_sha256=prepared.prepared_index.sha256,
        source_manifest_sha256=source_routing_context.contract.parent_source_manifest_sha256,
        working_manifest_sha256=completed_analysis_input.working_manifest.sha256,
        source_routing_sha256=source_routing_context.contract.sha256, plan_sha256=plan.sha256,
        annotation_stage_sha256=manifest.sha256, resource_session_sha256=resource,
        bundle_manifest_sha256=producer.binding('vep_bundle'),
        completed_analysis_input_sha256=c.sha256, completed_annotation_prefix_sha256=p.sha256,
        identities=plan.total_retained_identities, occurrences=completed_analysis_input.working_manifest.summary.retained_occurrences,
        sequence_requests=plan.sequence_requests)
    _require(all(parents[key] == value for key, value in expected.items())
             and stage.input_sha256 == producer.input_sha256
             and stage.run_options_sha256 == producer.preparation_run_options_sha256
             and stage.reference_sha256 == release.reference_manifest_sha256
             and stage.policy_sha256 == screening.policy.sha256
             and stage.resource_sha256s == (resource,))
    eligibility = load_analysis_eligibility_manifest(descriptor.file('eligibility_manifest').read(),
                                                    completion.eligibility_manifest_sha256)
    _require(eligibility.admission_binding_sha256 == descriptor.binding_sha256
             and eligibility.input_dispositions == descriptor.input_dispositions
             and eligibility.occurrence_dispositions == descriptor.occurrence_dispositions
             and eligibility.input_rows == completion.rows and eligibility.source_occurrences == completion.occurrences
             and eligibility.ordered_original_inputs_sha256 == completion.ordered_input_hashes_sha256
             and eligibility.ordered_original_proofs_sha256 == completion.ordered_proof_hashes_sha256
             and eligibility.ordered_screening_sha256 == completion.screening_sha256
             and eligibility.source_manifest_sha256 == parents['source_manifest_sha256']
             and eligibility.reference_manifest_sha256 == release.reference_manifest_sha256
             and eligibility.eligibility_policy_sha256 == screening.policy.sha256)
    role_map = dict(population_manifest='gnomad_joint_manifest', precision_manifest='gnomad_projection_precision_report',
        clinical_manifest='clinvar_manifest', exception_manifest='common_allele_exception_manifest',
        all_sites_manifest='gnomad_all_sites_manifest', public_index_manifest='common_allele_public_index_manifest')
    _require(dict(eligibility.resource_bindings) == {key: producer.binding(role) for key, role in role_map.items()})
    component = analysis_eligibility_component(descriptor.work_root / 'eligibility', eligibility)
    _require(component['rescue_manifest_sha256'] == descriptor.file('rescue_manifest').sha256)
    annotation_root = descriptor.work_root / 'annotation'
    from genosapiens.wgs_ledger import JournalEntry
    terminal = load_model(descriptor.file('annotation_journal').read(), JournalEntry)
    _require(terminal.event_type == 'stage_completed' and len(terminal.artifacts) == 1
             and terminal.artifacts[0] == Artifact('annotation-stage-closure.json',
                 descriptor.file('annotation_closure').sha256, descriptor.file('annotation_closure').size_bytes))
    _require(_json(descriptor.file('annotation_head').read()) == dict(
        schema='genosapiens.wgs-stage-journal-head', version=1, stage='annotation', ordinal=terminal.ordinal, sha256=terminal.sha256))
    annotation = BoundedAnnotationSource(Path(c.annotation_plan_root), plan, annotation_root,
        AnnotationStageResult(manifest, Artifact('annotation-manifest.json', manifest.sha256,
            descriptor.file('annotation_manifest').size_bytes), terminal, 0, annotation_prefix), descriptor.native_index)
    if descriptor.native_index is not None:
        from genosapiens.wgs_ledger import iter_artifact_lines
        with closing(iter_artifact_lines(annotation_root, descriptor.native_index)) as lines:
            _require(json.loads(next(lines)) == dict(schema='genosapiens.wgs-native-annotation-index',
                version=1, plan_sha256=plan.sha256, stage_sha256=manifest.sha256))
    annotation.verify(expected_plan=plan, source_manifest=source,
        working_manifest=completed_analysis_input.working_manifest, resource_session_sha256=resource,
        bundle_manifest_sha256=producer.binding('vep_bundle'), executable_stack_sha256=producer.executable_stack_sha256,
        runtime_image_digest=producer.runtime_image_digest, execution_config=execution_config,
        prepared_index=prepared.prepared_index, prepared_state_root=prepared.root.parent / 'state')
    for path, signature in signatures:
        _require(_signature(path) == signature)
    return producer, release, completion, eligibility, annotation, signatures


def resolve_completed_admission_input(descriptor, *, consumer_authority, consumer_release,
        completed_analysis_input, annotation_prefix, execution_config, source_routing_context):
    try:
        values = _authenticate_metadata(descriptor, completed_analysis_input=completed_analysis_input,
            annotation_prefix=annotation_prefix, execution_config=execution_config,
            source_routing_context=source_routing_context)
        producer, release, completion, eligibility, annotation, signatures = values
        result = ResolvedCompletedAdmissionInput(descriptor, producer, release, completion, eligibility,
            annotation, consumer_authority, consumer_release, signatures, _RESOLVED)
        result.require_consumer(consumer_authority, consumer_release,
            completed_analysis_input=completed_analysis_input, annotation_prefix=annotation_prefix)
        return result
    except (OSError, sqlite3.Error, KeyError, TypeError, ValueError, UnicodeError):
        raise AdmissionInputError('completed admission input is invalid') from None


_CONTROL_NAMES = dict(
    producer_authority='candidate-authority.json', producer_release='candidate-candidate-release.json',
    producer_plan='candidate-execution-plan.json', producer_start='candidate-analysis-start.json',
    producer_terminal='candidate-analysis-terminal.json', runtime='runtime-analysis.json',
    runtime_sample='.runtime-sample-analysis.json', runtime_phases='.runtime-phases-analysis.json',
    failure='failure-analysis.json',
)


def _write_new(path, content):
    descriptor = os.open(path, os.O_WRONLY | os.O_CREAT | os.O_EXCL | os.O_NOFOLLOW, 0o600)
    with os.fdopen(descriptor, 'wb') as stream:
        stream.write(content)
        stream.flush()
        os.fsync(stream.fileno())
    directory = os.open(path.parent, os.O_RDONLY | os.O_DIRECTORY)
    try:
        os.fsync(directory)
    finally:
        os.close(directory)


def _resolve_original_dependencies(descriptor, execution_config):
    """Use the ordinary C/P readers under the original admission producer."""
    from genosapiens.wgs_analysis_inputs import load_completed_analysis_input, resolve_completed_analysis_input
    from genosapiens.wgs_annotation_prefix import load_completed_annotation_prefix, resolve_completed_annotation_prefix
    from genosapiens.wgs_source_routing import NativeSourceContract, SourceRoutingContext
    from genosapiens.allele_identity import ReferenceRole
    from genosapiens.reference_concordance import ReferenceRoleBinding
    from genosapiens.submission_rendering import load_submission_contig_renderer
    from genosapiens.wgs_public_runner import _load_analysis_resources
    producer, release = _producer(descriptor, execution_config)
    c_file, p_file = descriptor.file('completed_analysis_input'), descriptor.file('completed_annotation_prefix')
    c = load_completed_analysis_input(c_file.read(), expected_sha256=c_file.sha256)
    p = load_completed_annotation_prefix(p_file.read(), expected_sha256=p_file.sha256)
    resources = _load_analysis_resources()
    role = 'public_annotation_benchmark'
    identity = resources.reference_contract.reference_role(role)
    native = NativeSourceContract(json.loads(descriptor.original_binding)['source_manifest_sha256'],
        ReferenceRole(role), identity.identity_id, release.reference_manifest_sha256, ('M',), None)
    reference = resources.population_reference
    renderer = load_submission_contig_renderer(identity_path=reference.identity_path,
        reference_manifest_path=reference.manifest_path, reference_root=reference.reference_root,
        reference_binding=ReferenceRoleBinding(role, identity.identity_id))
    routing = SourceRoutingContext(native, renderer)
    completed = resolve_completed_analysis_input(c, consumer_authority=producer, consumer_release=release,
        execution_config=execution_config, source_routing_context=routing)
    prefix = resolve_completed_annotation_prefix(p, plan=completed.transport.plan, execution_config=execution_config)
    prefix.require_consumer(producer, release)
    return dict(completed_analysis_input=completed, annotation_prefix=prefix,
                source_routing_context=routing, execution_config=execution_config)


def generate_completed_admission_input(inventory_content, *, inventory_sha256, output):
    """Create one private descriptor from explicit, freshly verified inventory.

    Failure leaves diagnostic control copies, never an authenticated handoff.
    The output directory must be new; generator reruns cannot overwrite it.
    """
    from genosapiens.bounded_wgs_config import load_bounded_wgs_execution_config
    from genosapiens.wgs_analysis_eligibility import load_analysis_eligibility_manifest
    from genosapiens.wgs_admission_projection import load_admission_projection_completion
    _require(type(inventory_content) is bytes and len(inventory_content) <= 1024 * 1024
             and _digest(inventory_sha256) and _sha(inventory_content) == inventory_sha256)
    inventory = json.loads(inventory_content, object_pairs_hook=_unique)
    _require(inventory['schema'] == 'genosapiens.wgs-discovery-repair-asset-inventory'
             and type(inventory['version']) is int and inventory['version'] == 1)
    output = _path(str(output))
    _require(os.getuid() == 65532 and not output.parent.exists() and not output.parent.is_symlink()
             and output.name == 'completed-admission-input.json')
    for parent in output.parents[1:]:
        _require(parent.is_dir() and not parent.is_symlink())
    assets = inventory['assets']
    _require(_FILES <= assets.keys())
    original_work = Path(assets['admission_database']['path']).parent.parent
    _require(not output.parent.is_relative_to(original_work)
             and not original_work.is_relative_to(output.parent))
    # Authenticate every requested small source before creating output. Host
    # controls can be readable 0644 files; their new copies use AnalysisInputFile.
    contents = {}
    for role in sorted(_FILES):
        row = assets[role]
        path = _path(row['path'])
        if role in _CONTROL_NAMES:
            _require(path.name == _CONTROL_NAMES[role])
        content = read_bounded_metadata(path)
        _require(len(content) == row['size_bytes'] and _sha(content) == row['sha256'])
        contents[role] = content
    output.parent.mkdir(mode=0o700)
    control_root = output.parent / 'producer'
    control_root.mkdir(mode=0o700)
    files = []
    for role in sorted(_FILES):
        row = assets[role]
        path = Path(row['path'])
        if role in _CONTROL_NAMES:
            path = control_root / _CONTROL_NAMES[role]
            _write_new(path, contents[role])
        files.append((role, AnalysisInputFile(str(path), row['sha256'], row['size_bytes'])))
    def database(role, version):
        row = assets[role]
        return AdmissionDatabase(row['path'], row['size_bytes'], version)
    binding = canonical(inventory['admission_binding'])
    completion = load_admission_projection_completion(contents['admission_completion'],
        expected_sha256=assets['admission_completion']['sha256'])
    eligibility = load_analysis_eligibility_manifest(contents['eligibility_manifest'],
        assets['eligibility_manifest']['sha256'])
    native = assets.get('native_index')
    if native is not None:
        _require(Path(native['path']) == original_work / 'annotation/native-index.jsonl')
    descriptor = CompletedAdmissionInput(tuple(files), database('admission_database', 5),
        database('admission_stage_state', 1), database('annotation_state', 1), binding.decode('ascii'),
        inventory['admission_binding_sha256'], inventory['stage_head_sha256'],
        inventory['producer_executable_sha256'], inventory['runtime_image_digest'],
        inventory['execution_config_sha256'], eligibility.input_dispositions, eligibility.occurrence_dispositions,
        completion, None if native is None else Artifact('native-index.jsonl', native['sha256'], native['size_bytes']))
    config = load_bounded_wgs_execution_config(read_bounded_metadata(CONFIG_PATH), descriptor.execution_config_sha256)
    arguments = _resolve_original_dependencies(descriptor, config)
    _authenticate_metadata(descriptor, **arguments, _new_control_root=control_root)
    _write_new(output, descriptor.as_bytes())
    _require(AnalysisInputFile(str(output), descriptor.sha256, len(descriptor.as_bytes())).read() == descriptor.as_bytes())
    return descriptor


def main(argv=None):
    import argparse
    parser = argparse.ArgumentParser(description='Authenticate completed public layout-v5 admission for discovery reuse')
    parser.add_argument('--inventory', required=True)
    parser.add_argument('--inventory-sha256', required=True)
    parser.add_argument('--output', required=True)
    arguments = parser.parse_args(argv)
    descriptor = generate_completed_admission_input(read_bounded_metadata(Path(arguments.inventory)),
        inventory_sha256=arguments.inventory_sha256, output=arguments.output)
    print(canonical(dict(path=arguments.output, sha256=descriptor.sha256,
                         size_bytes=len(descriptor.as_bytes()))).decode('ascii'), end='')
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
