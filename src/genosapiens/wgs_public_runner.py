"""Docker entry point for one exact public bounded-WGS qualification stage."""

from __future__ import annotations

import argparse
from contextlib import contextmanager
import hashlib
import json
import os
from pathlib import Path
import threading
import time

from genosapiens.allele_identity import ReferenceRole
from genosapiens.bounded_json import MAX_GOVERNING_JSON_BYTES, read_bounded_bytes
from genosapiens.bounded_wgs_config import load_bounded_wgs_execution_config
from genosapiens.common_allele_filter import load_common_allele_filter_policy
from genosapiens.common_allele_public_indexes import (
    CommonAlleleIndexError,
    verify_common_allele_public_indexes,
)
from genosapiens.gnomad_faf95_adapter import (
    load_gnomad_faf95_precision_contract,
    verify_gnomad_faf95_precision_qualification,
)
from genosapiens.population_annotation import PopulationQueryLimits
from genosapiens.reference_concordance import (
    ReferenceConcordanceError,
    load_concordance_report,
    load_public_validation_concordance_report,
    open_installed_reference,
)
from genosapiens.track1_public_validation import (
    HPO_MANIFEST_SHA256,
    REFERENCE_IDENTITY_PATH,
    REFERENCE_MANIFEST_SHA256,
    RUNTIME_IMAGE,
    load_public_policy_family,
    load_public_preparation_resources,
    load_public_resource_context,
)
from genosapiens.track1_public_validation_coordinator import executable_stack_bytes
from genosapiens.track1_release import load_base_release
from genosapiens.wgs_admission import NativeSourceContract
from genosapiens.wgs_analysis import analyze_bounded_case
from genosapiens.wgs_failures import failure_boundary
from genosapiens.wgs_ledger import (
    FilterManifest,
    JournalEntry,
    LedgerManifest,
    MAX_MANIFEST_BYTES,
    PreparedCaseManifest,
    WorkingSetManifest,
    canonical,
    load_model,
    read_journal,
)
from genosapiens.wgs_preparation import (
    PreparationRunOptions,
    SourceAccountingClosure,
    prepare_bounded_case,
    source_worker_reference,
)
from genosapiens.wgs_prepared_case import CompletedJournalBinding
from genosapiens.wgs_public_analysis import PublicBoundedAnalysisBuilder
from genosapiens.wgs_public_filter import PublicEarlyFilterDecisionProvider, PublicFilterWorkerConfig
from genosapiens.wgs_public_qualification import (
    WgsPublicQualificationError,
    load_public_qualification_authority,
    public_qualification_hpo_root,
    require_private_qualification_hpo_input,
    require_public_qualification_release,
)
from genosapiens.wgs_runtime import (
    managed_runtime_phase_recorder,
    runtime_phase,
)


REFERENCE_ROOT = Path("/mnt/reference")
REPOSITORY_ROOT = Path(__file__).resolve().parents[2]
REFERENCE_MANIFEST = REFERENCE_ROOT / (
    "manifests/grch38/"
    "challenge-grch38-masked-v2-plus-hs38d1.a92c50f61c1b.json"
)
CONFIG_PATH = REPOSITORY_ROOT / "config/track1/bounded-wgs-execution-v5.json"
FILTER_POLICY_PATH = REPOSITORY_ROOT / "config/track1/common-allele-prefilter-v1.json"
PRECISION_CONTRACT_PATH = REPOSITORY_ROOT / (
    "config/public-resources/gnomad-joint-faf95-precision-v1.json"
)
PRECISION_REPORT_PATH = REFERENCE_ROOT / (
    "manifests/gnomad/4.1/faf95-precision/"
    "4fb171a7fb5a23ed8ed24dc81b35880f84fa79948b25393048cbcdc991410303.json"
)
INDEX_MANIFEST_PATH = (
    "manifests/common-allele-public-indexes/"
    "c8e7740ccc3b024c8502b035529ad2fc185a5d9dd18a32f9b7cc421d7d8fac4f-"
    "cache-bounded-runtime.json"
)


class PublicRunnerError(RuntimeError):
    pass


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as source:
        while chunk := source.read(1024 * 1024):
            digest.update(chunk)
    return digest.hexdigest()


def _bound(path: Path, expected_sha256: str, *, maximum_bytes=MAX_GOVERNING_JSON_BYTES) -> bytes:
    try:
        if path.is_symlink() or not path.is_file():
            raise ValueError
        content = read_bounded_bytes(path, maximum_bytes=maximum_bytes)
        if hashlib.sha256(content).hexdigest() != expected_sha256:
            raise ValueError
        return content
    except (OSError, ValueError):
        raise PublicRunnerError("public runner binding changed") from None


def _manifest_dependency(role, path, digest):
    # Source/filter/working manifests share the ledger's 8 MiB contract.
    # Governing policy documents retain their separate, smaller default bound.
    content = _bound(path, digest, maximum_bytes=MAX_MANIFEST_BYTES)
    return role, str(path.absolute()), digest, len(content)


def _selected_sample(arguments, authority):
    from genosapiens.wgs_operator_invocation import OperatorCaseBinding
    if type(authority) is OperatorCaseBinding:
        return authority.selected_sample
    from genosapiens.wgs_work_contract import FRESH_INVOCATION, fresh_sample_binding
    supplied = getattr(arguments, 'selected_sample', None)
    if authority.invocation_id == FRESH_INVOCATION:
        if fresh_sample_binding(supplied, authority.input_sha256) != authority.binding('fresh_selected_sample'):
            raise PublicRunnerError('fresh selected sample differs from authority')
        return supplied
    if supplied is not None:
        raise PublicRunnerError('selected sample argument requires fresh qualification authority')
    return 'HG002'


@failure_boundary(
    "resource_authentication", "runner_authority_authentication",
    known_errors=(PublicRunnerError, WgsPublicQualificationError, ValueError),
    category="authentication",
)
def _load_common(arguments):
    authority_bytes = _bound(
        arguments.authority, arguments.authority_sha256,
    )
    authority = load_public_qualification_authority(
        authority_bytes, expected_sha256=arguments.authority_sha256,
    )
    require_private_qualification_hpo_input(authority)
    release_bytes = _bound(arguments.release, arguments.release_sha256)
    release = load_base_release(
        release_bytes, expected_sha256=arguments.release_sha256,
    )
    require_public_qualification_release(authority, release)
    _selected_sample(arguments, authority)
    config_bytes = _bound(CONFIG_PATH, authority.execution_config_sha256)
    config = load_bounded_wgs_execution_config(
        config_bytes, authority.execution_config_sha256,
    )
    stack_sha256 = hashlib.sha256(
        executable_stack_bytes(REPOSITORY_ROOT)
    ).hexdigest()
    if (
        stack_sha256 != authority.executable_stack_sha256
        or authority.runtime_image_digest != release.runtime_image_digest
    ):
        raise PublicRunnerError("public runner executable binding changed")
    return (
        authority_bytes, authority, release_bytes, release,
        config, stack_sha256,
    )


@failure_boundary(
    "resource_authentication", "runner_resource_authentication",
    known_errors=(PublicRunnerError, ValueError), category="authentication",
)
def _load_preparation_resources():
    from genosapiens.track1_early_scope import (
        load_completed_early_scope, load_early_scope_policy, validate_scope_pin)
    from genosapiens.bounded_json import load_bounded_json
    value, _ = load_bounded_json(REPOSITORY_ROOT / "track1-early-scope-pin.json",
        reference_root=REPOSITORY_ROOT, maximum_bytes=65536)
    pin = validate_scope_pin(value)
    manifest = load_completed_early_scope(REFERENCE_ROOT, pin)
    if manifest.document["request"]["policy_sha256"] != load_early_scope_policy(REPOSITORY_ROOT).sha256:
        raise PublicRunnerError("installed early scope policy differs")
    resources = load_public_preparation_resources(
        REPOSITORY_ROOT,
        REFERENCE_ROOT,
        REFERENCE_ROOT / "scratch/public-build",
    )
    return resources, pin, manifest


@failure_boundary(
    "resource_authentication", "runner_resource_authentication",
    known_errors=(PublicRunnerError, ValueError), category="authentication",
)
def _load_analysis_resources():
    from genosapiens.wgs_science_hpo import ScienceHpoLoader
    return load_public_resource_context(
        REPOSITORY_ROOT,
        REFERENCE_ROOT,
        REFERENCE_ROOT / "scratch/public-build",
        analysis_hpo_loader=ScienceHpoLoader(),
    )


@failure_boundary(
    "resource_authentication", "runner_index_reconciliation",
    known_errors=(PublicRunnerError, CommonAlleleIndexError, ValueError),
    category="authentication",
)
def _indexes(authority, *, verify_sources: bool):
    from genosapiens.installed_public_index import load_installed_index_pin
    pin = load_installed_index_pin(REPOSITORY_ROOT)
    expected = authority.binding('common_allele_public_index_manifest')
    if pin and pin['manifest_sha256'] != expected:
        raise PublicRunnerError('installed public index binding differs')
    return verify_common_allele_public_indexes(
        reference_root=REFERENCE_ROOT,
        manifest_path=pin['manifest_path'] if pin else INDEX_MANIFEST_PATH,
        expected_manifest_sha256=expected,
        verify_sources=verify_sources,
        reuse_qualified_reconciliation=not verify_sources,
    )


def _preparation(arguments, phase_recorder=None, *, operator_case=None,
                 executor=None) -> int:
    run_options = PreparationRunOptions(arguments.verify_source_ledger)
    with runtime_phase(phase_recorder, "authority_validation"):
        (
            authority_bytes, authority, release_bytes, _release,
            config, stack_sha256,
        ) = _stage_common(arguments, operator_case)
        if authority.preparation_attempts != 1:
            raise PublicRunnerError("analysis-only authority cannot start preparation")
        if run_options.sha256 != authority.preparation_run_options_sha256:
            raise PublicRunnerError("public preparation audit mode differs from authority")
    with runtime_phase(phase_recorder, "reference_initialization"):
        concordance_content = _bound(
            arguments.concordance, arguments.concordance_sha256,
        )
    try:
        concordance = load_concordance_report(arguments.concordance)
    except ReferenceConcordanceError:
        concordance = load_public_validation_concordance_report(
            arguments.concordance
        )
    if (
        concordance.input_sha256 != authority.input_sha256
        or concordance.validated_input_sha256 != authority.input_sha256
        or concordance.reference_manifest_sha256 != REFERENCE_MANIFEST_SHA256
        or not concordance.accepted
        or concordance.reference_role != (
            "public_annotation_benchmark"
            if authority.invocation_id in {"complete-hg002-deepvariant-wgs", "focused-fresh-subsample"}
            else "challenge_alignment"
        )
        or concordance_content != _bound(
            arguments.concordance, arguments.concordance_sha256,
        )
    ):
        raise PublicRunnerError("public concordance authority changed")
    filter_policy_bytes = _bound(
        FILTER_POLICY_PATH, authority.binding("common_allele_policy"),
    )
    filter_policy = load_common_allele_filter_policy(
        filter_policy_bytes, authority.binding("common_allele_policy"),
    )
    precision_contract_bytes = _bound(
        PRECISION_CONTRACT_PATH,
        authority.binding("gnomad_projection_precision_contract"),
    )
    precision_contract = load_gnomad_faf95_precision_contract(
        precision_contract_bytes,
        authority.binding("gnomad_projection_precision_contract"),
    )
    precision_report = verify_gnomad_faf95_precision_qualification(
        _bound(
            PRECISION_REPORT_PATH,
            authority.binding("gnomad_projection_precision_report"),
        ),
        authority.binding("gnomad_projection_precision_report"),
        precision_contract,
    )
    if authority.completed_source_input is not None:
        from genosapiens.wgs_preparation import require_completed_source_readonly
        require_completed_source_readonly(authority.completed_source_input)
    with runtime_phase(phase_recorder, "resource_initialization"):
        resources, scope_pin, scope_manifest = _load_preparation_resources()
    with runtime_phase(phase_recorder, "public_index_reconciliation"):
        indexes = _indexes(authority, verify_sources=False)
    limits = PopulationQueryLimits(
        config.public_queries.batch_rows,
        config.public_queries.batch_rows,
        65536,
        1_000_000,
        config.public_queries.max_response_bytes,
        config.public_queries.timeout_seconds,
    )

    worker_config = PublicFilterWorkerConfig(
        resources.population_joint, resources.population_all_sites,
        resources.population_reference, RUNTIME_IMAGE, limits,
        config.public_queries, indexes,
        Path(authority.output_root) / "preparation" / "prefilter-staging",
    )

    reference_role = resources.reference_contract.reference_role(
        concordance.reference_role
    )
    from genosapiens.reference_concordance import ReferenceRoleBinding
    from genosapiens.submission_rendering import load_submission_contig_renderer
    contig_renderer = load_submission_contig_renderer(
        identity_path=resources.population_reference.identity_path,
        reference_manifest_path=resources.population_reference.manifest_path,
        reference_root=resources.population_reference.reference_root,
        reference_binding=ReferenceRoleBinding(concordance.reference_role, reference_role.identity_id),
    )

    from genosapiens.wgs_source_routing import SourceRoutingContext

    def source_routing_context_provider(context):
        native = NativeSourceContract(
            context.source_manifest.sha256,
            ReferenceRole(concordance.reference_role),
            reference_role.identity_id,
            REFERENCE_MANIFEST_SHA256,
            ("M",),
            None,
        )
        return SourceRoutingContext(native, contig_renderer)

    def decision_provider(context):
        routing = source_routing_context_provider(context)
        return PublicEarlyFilterDecisionProvider(
            filter_policy,
            precision_report,
            routing,
            worker_config,
            scope_session,
            phase_recorder,
        )(context)

    case_sha256 = hashlib.sha256(canonical({
        "schema": ("genosapiens.public-wgs-qualification-case" if operator_case is None
                   else "genosapiens.operator-wgs-case"),
        "invocation_id": authority.invocation_id,
        "input_sha256": authority.input_sha256,
        "hpo_input_sha256": authority.hpo_input_sha256,
    })).hexdigest()
    selected_sample_name = _selected_sample(arguments, authority)
    selected_sample_sha256 = hashlib.sha256(canonical({
        "schema": ("genosapiens.public-wgs-qualification-sample" if operator_case is None
                   else "genosapiens.operator-wgs-sample"),
        "sample_name": selected_sample_name,
        "input_sha256": authority.input_sha256,
    })).hexdigest()
    def worker_reference(role_name, aliases=()):
        role = resources.reference_contract.reference_role(role_name)
        reference_root = getattr(resources, "reference_root", REFERENCE_ROOT)
        fasta = reference_root / resources.reference_contract.expected_output(
            role.fasta_output_role
        ).path
        fai = reference_root / resources.reference_contract.expected_output(
            role.fai_output_role
        ).path
        return source_worker_reference(fasta, fai, aliases=aliases)

    def run(fetch_reference, reference_binding):
        from genosapiens.track1_early_scope import EarlyScopeSession
        nonlocal scope_session
        with EarlyScopeSession(REFERENCE_ROOT, scope_pin, renderer=contig_renderer) as scope_session:
            return run_with_scope(fetch_reference, reference_binding)

    scope_session = None
    def run_with_scope(fetch_reference, reference_binding):
        result = (prepare_bounded_case if executor is None else executor)(
            output_root=Path(authority.output_root) / "preparation",
            original_input_path=Path(authority.input_path),
            hpo_input_path=Path(authority.hpo_input_path),
            hpo_input_sha256=authority.hpo_input_sha256,
            hpo_input_size_bytes=authority.hpo_input_size_bytes,
            hpo_ontology=resources.hpo_store.graph.ontology,
            hpo_reference_root=(public_qualification_hpo_root(authority)
                if operator_case is None else operator_case.hpo_root),
            concordance_report=concordance,
            fetch_reference=fetch_reference,
            source_worker_reference=reference_binding,
            selected_sample_name=selected_sample_name,
            selected_sample_sha256=selected_sample_sha256,
            case_manifest_sha256=case_sha256,
            ontology_manifest_sha256=HPO_MANIFEST_SHA256,
            executable_stack_sha256=stack_sha256,
            runtime_image_digest=authority.runtime_image_digest,
            execution_config=config,
            run_options=run_options,
            policy_sha256=authority.binding("common_allele_policy"),
            scope_policy_sha256=scope_manifest.document["request"]["policy_sha256"],
            scope_manifest_sha256=scope_manifest.sha256,
            population_manifest_sha256=authority.binding(
                "gnomad_joint_manifest"
            ),
            clinical_manifest_sha256=authority.binding("clinvar_manifest"),
            exception_manifest_sha256=authority.binding(
                "common_allele_exception_manifest"
            ),
            precision_manifest_sha256=authority.binding(
                "gnomad_projection_precision_report"
            ),
            decision_provider=decision_provider,
            source_routing_context_provider=source_routing_context_provider,
            base_release_bytes=release_bytes,
            expected_base_release_sha256=arguments.release_sha256,
            public_qualification_authority_bytes=authority_bytes,
            expected_public_qualification_authority_sha256=(
                arguments.authority_sha256 if operator_case is None else None
            ),
            phase_recorder=phase_recorder,
            **({} if operator_case is None else dict(
                resume=arguments.resume,
                operator_case_binding=operator_case,
                completed_source_input=operator_case.completed_source_input)),
        )
        return result

    if concordance.reference_role == "public_annotation_benchmark":
        aliases = resources.source_to_reference_contigs

        def public_fetch(contig, position, length):
            return resources.reference.fetch(
                aliases[contig], position, length,
            )

        result = run(public_fetch, worker_reference(
            "public_annotation_benchmark",
            tuple(sorted(aliases.items())),
        ))
    else:
        identity_path = REPOSITORY_ROOT / REFERENCE_IDENTITY_PATH
        with open_installed_reference(
            resources.reference_contract,
            identity_path=identity_path,
            reference_manifest_path=REFERENCE_MANIFEST,
            reference_root=REFERENCE_ROOT,
            reuse_published_validation=True,
        ) as (manifest_sha256, fetch_reference, _reference_contigs):
            if manifest_sha256 != REFERENCE_MANIFEST_SHA256:
                raise PublicRunnerError("installed challenge reference changed")
            result = run(fetch_reference, worker_reference("challenge_alignment"))
    print(result.receipt.status, result.receipt.terminal_stage, flush=True)
    return 0 if result.receipt.status == "passed" else 1


@failure_boundary(
    "input_authentication", "runner_input_authentication",
    known_errors=(PublicRunnerError, ValueError), category="authentication",
)
def _load_manifest(path: Path, kind):
    return load_model(read_bounded_bytes(
        path, maximum_bytes=MAX_MANIFEST_BYTES, reference_root=path.parent,
    ), kind)


@failure_boundary(
    "input_authentication", "runner_input_authentication",
    known_errors=(PublicRunnerError, ValueError), category="authentication",
)
def _journal_binding(root: Path, stage: str, expected_head: str):
    first = _load_manifest(root / f"{stage}-000001.journal.json", JournalEntry)
    entries = read_journal(
        root,
        stage=stage,
        parent_sha256=first.stage_parent_sha256,
        expected_head_sha256=expected_head,
    )
    if not entries or entries[-1].event_type != "stage_completed":
        raise PublicRunnerError("prepared journal is incomplete")
    return CompletedJournalBinding(
        stage, root, first.stage_parent_sha256, expected_head,
    )


@contextmanager
def _completed_input_phase(recorder,phase):
    """Bounded diagnostic detail within the existing runtime phase vocabulary."""
    allowed=('completed_analysis_validation','completed_annotation_validation','completed_admission_validation')
    if phase not in allowed:raise PublicRunnerError('invalid completed-input progress phase')
    if recorder is None:
        yield
        return
    timings=getattr(recorder,'_completed_input_timings',None)
    if timings is None:
        timings={};recorder._completed_input_timings=timings
    began=time.monotonic();stop=threading.Event();errors=[]
    def emit(status):
        content=(json.dumps(dict(schema='genosapiens.completed-input-progress',version=1,
            attempt_id=recorder.attempt_id,execution_config_sha256=recorder.execution_config_sha256,
            phase=phase,status=status,phase_elapsed_seconds=time.monotonic()-began,
            completed_seconds=timings),sort_keys=True,separators=(',',':'))+'\n').encode('ascii')
        if len(content)>64*1024:raise PublicRunnerError('completed-input progress exceeds capacity')
        temporary=recorder.root/'completed-input-progress.tmp'
        fd=os.open(temporary,os.O_WRONLY|os.O_CREAT|os.O_TRUNC|os.O_NOFOLLOW,0o600)
        with os.fdopen(fd,'wb') as out:out.write(content)
        os.replace(temporary,recorder.root/'completed-input-progress.json')
    def heartbeat():
        while not stop.wait(1):
            try:emit('running')
            except BaseException as error:errors.append(error);return
    emit('running')
    thread=threading.Thread(target=heartbeat,name='completed-input-progress',daemon=True);thread.start()
    failed=False
    try:
        with runtime_phase(recorder,'prepared_authority_validation'):yield
    except BaseException:
        failed=True
        raise
    finally:
        stop.set();thread.join(timeout=30)
        if not thread.is_alive():
            timings[phase]=time.monotonic()-began
            try:emit('failed' if failed else 'complete')
            except BaseException:
                if not failed:raise
        if not failed and (thread.is_alive() or errors):
            raise PublicRunnerError('completed-input progress failed')


def _stage_common(arguments, operator_case):
    """Keep public candidate and accepted operator authority separate."""
    if operator_case is None:
        return _load_common(arguments)
    from genosapiens.wgs_operator_invocation import OperatorCaseBinding
    if type(operator_case) is not OperatorCaseBinding:
        raise PublicRunnerError('operator case binding is invalid')
    return operator_case.common(arguments)


def _analysis(arguments, phase_recorder=None, *, public_pilot_inputs=None,
              operator_case=None, executor=None) -> int:
    from genosapiens.wgs_analysis_policy import select_analysis_profile
    mode, scope_policy = select_analysis_profile(exhaustive=getattr(arguments, 'exhaustive', False),
        repository_root=REPOSITORY_ROOT)
    with runtime_phase(phase_recorder, "authority_validation"):
        (
            authority_bytes, authority, release_bytes, release,
            config, stack_sha256,
        ) = _stage_common(arguments, operator_case)
    # Only the explicit public candidate harness supplies finite inputs. There
    # is no operator CLI switch or ordinary whole-genome fallback for them.
    pilot_roles={'public_profile_pilot_acceptance','public_profile_pilot_cohort','public_profile_pilot_context_inputs',
        'public_profile_pilot_scientific_inputs','public_profile_pilot_primary_inputs'}
    pilot_bindings={r.role:r.sha256 for r in authority.resource_bindings if r.role in pilot_roles}
    if public_pilot_inputs is None:
        if pilot_bindings: raise PublicRunnerError('finite public pilot inputs are unavailable')
    else:
        from genosapiens.wgs_profile_projection import PublicProfilePilotInputs
        from genosapiens.track1_roster import PROFILE_ACCEPTANCE_SHA256
        if type(public_pilot_inputs) is not PublicProfilePilotInputs:
            raise PublicRunnerError('finite public pilot input type differs')
        expected_pilot = {'public_profile_pilot_acceptance':PROFILE_ACCEPTANCE_SHA256,
            'public_profile_pilot_cohort':public_pilot_inputs.binding_sha256}
        if public_pilot_inputs.completed_context_inputs is not None:
            from genosapiens.wgs_profile_projection import CompletedProfileDiscoveryInputs
            completed = public_pilot_inputs.completed_context_inputs
            if type(completed) is not CompletedProfileDiscoveryInputs:
                raise PublicRunnerError('completed context input type differs')
            completed.require_authority(authority, consumer_executable_sha256=stack_sha256)
            expected_pilot['public_profile_pilot_context_inputs'] = completed.sha256
        if public_pilot_inputs.completed_scientific_inputs is not None:
            from genosapiens.wgs_scientific_checkpoint import CompletedScientificInputs
            completed = public_pilot_inputs.completed_scientific_inputs
            if type(completed) is not CompletedScientificInputs:
                raise PublicRunnerError('completed scientific input type differs')
            completed.require_authority(authority, consumer_executable_sha256=stack_sha256)
            expected_pilot['public_profile_pilot_scientific_inputs'] = completed.sha256
        if (type(public_pilot_inputs) is not PublicProfilePilotInputs
                or authority.invocation_id!='bounded-public-pilot'
                or pilot_bindings != expected_pilot):
            raise PublicRunnerError('finite public pilot authority differs')
    from genosapiens.track1_analysis_profiles import authenticate_analysis_profile
    release.require_analysis_profile(authenticate_analysis_profile(mode, policy=scope_policy))
    from genosapiens.wgs_analysis_policy import load_analysis_selectors, load_analysis_work_contract
    eligibility_policy,writer_profile=load_analysis_selectors(authority,repository_root=REPOSITORY_ROOT)
    focused_work=load_analysis_work_contract(authority,repository_root=REPOSITORY_ROOT)
    analysis_input = None
    completed_input = None
    supplied_path = getattr(arguments, "completed_analysis_input", None)
    supplied_sha = getattr(arguments, "completed_analysis_input_sha256", None)
    bound_sha = next((row.sha256 for row in authority.resource_bindings
                      if row.role == "completed_analysis_input"), None)
    if (supplied_path is None) != (supplied_sha is None) or (bound_sha is None) != (supplied_path is None):
        raise PublicRunnerError("completed analysis input arguments are incomplete")
    if supplied_path is not None:
        from genosapiens.wgs_analysis_inputs import load_completed_analysis_input
        if supplied_sha != bound_sha:
            raise PublicRunnerError("completed analysis input is not authorized")
        analysis_input = load_completed_analysis_input(read_bounded_bytes(supplied_path),
                                                       expected_sha256=supplied_sha)
    prefix_path = getattr(arguments, "completed_annotation_prefix", None)
    prefix_sha = getattr(arguments, "completed_annotation_prefix_sha256", None)
    bound_prefix = next((r.sha256 for r in authority.resource_bindings if r.role == "completed_annotation_prefix"), None)
    if ((prefix_path is None) != (prefix_sha is None) or (bound_prefix is None) != (prefix_path is None)
            or prefix_path is not None and (analysis_input is None or prefix_sha != bound_prefix)):
        raise PublicRunnerError("completed annotation prefix arguments are incomplete")
    admission_path = getattr(arguments, "completed_admission_input", None)
    admission_sha = getattr(arguments, "completed_admission_input_sha256", None)
    bound_admission = next((r.sha256 for r in authority.resource_bindings if r.role == "completed_admission_input"), None)
    if ((admission_path is None) != (admission_sha is None) or (bound_admission is None) != (admission_path is None)
            or admission_path is not None and (prefix_path is None or analysis_input is None or admission_sha != bound_admission)):
        raise PublicRunnerError("completed admission input arguments are incomplete or unbound")
    admission_descriptor = None
    if admission_path is not None:
        from genosapiens.wgs_admission_inputs import load_completed_admission_input, MAX_DESCRIPTOR_BYTES
        admission_descriptor = load_completed_admission_input(read_bounded_bytes(admission_path,
            maximum_bytes=MAX_DESCRIPTOR_BYTES), expected_sha256=admission_sha)
    preparation = (Path(authority.output_root) / "preparation" if analysis_input is None
                   else Path(analysis_input.preparation_root))
    filter_root = preparation / "prefilter"
    working_root = preparation / "working-set"
    publication_root = preparation / "prepared"
    filtered = _load_manifest(filter_root / "manifest.json", FilterManifest)
    working = _load_manifest(
        working_root / "manifest.json", WorkingSetManifest,
    )
    prepared = _load_manifest(
        publication_root / "prepared-case-manifest.json",
        PreparedCaseManifest,
    )
    from genosapiens.wgs_preparation import load_completed_source_input
    descriptor_path = Path(prepared.completed_source_input_path)
    descriptor = load_completed_source_input(read_bounded_bytes(descriptor_path))
    if descriptor.sha256 != prepared.completed_source_input.sha256:
        raise PublicRunnerError("completed source input changed")
    if prepared.source_origin == "completed_input":
        from genosapiens.wgs_preparation import require_completed_source_readonly, CompletedSourceArtifact, COMPLETED_SOURCE_SCHEMA
        require_completed_source_readonly(CompletedSourceArtifact(prepared.completed_source_input_path,
            COMPLETED_SOURCE_SCHEMA, 1, descriptor.sha256, len(descriptor.as_bytes())))
    source_root = Path(descriptor.source_root)
    source = _load_manifest(source_root / "manifest.json", LedgerManifest)
    if operator_case is not None:
        operator_case.require_prepared(source, prepared)
    source_accounting = SourceAccountingClosure(
        source.sha256,
        source.source.input_sha256,
        source.header.sha256,
        source.summary.reports[0].sha256,
        source.summary.reports[1].sha256,
        source.summary.reports[2].sha256,
        source.summary.total_records,
        source.summary.total_physical_alternate_slots,
        source.summary.total_biological_alternates,
    )
    journal_heads = {
        value.role.removesuffix("_journal"): value.sha256
        for value in prepared.journal_heads
    }
    journals = tuple(
        _journal_binding(
            (Path(descriptor.source_journal_root) if stage == "source" else
             Path(descriptor.source_audit_journal_root) if stage == "source_audit" else
             preparation / "journals" / stage),
            stage,
            journal_heads[stage],
        )
        for stage in sorted(journal_heads)
    )
    phenotype_path = preparation / "phenotype/phenotype-ingestion.json"
    phenotype_sha256 = _sha256(phenotype_path)
    early_snapshot_sha256 = _sha256(
        publication_root / "ledger-bound-early-snapshot.json"
    )
    with runtime_phase(phase_recorder, "analysis_resource_initialization"):
        resources = _load_analysis_resources()
        policies = load_public_policy_family(
            REPOSITORY_ROOT, model_loci=resources.model_loci,
        )
    from genosapiens.wgs_analysis_policy import select_analysis_runtime
    runtime = (select_analysis_runtime(authority, release, policies, repository_root=REPOSITORY_ROOT)
               if operator_case is None else operator_case.select_runtime(policies))
    with runtime_phase(phase_recorder, "public_index_reconciliation"):
        indexes = _indexes(authority, verify_sources=False)
    source_reference_role = (
        "public_annotation_benchmark"
        if (authority.invocation_id in {"complete-hg002-deepvariant-wgs", "focused-fresh-subsample"}
            or authority.is_finite_profile_pilot)
        else "challenge_alignment"
    )
    reference_role = resources.reference_contract.reference_role(
        source_reference_role
    )
    native = NativeSourceContract(
        source.sha256,
        ReferenceRole(source_reference_role),
        reference_role.identity_id,
        REFERENCE_MANIFEST_SHA256,
        ("M",),
        None,
    )
    from genosapiens.reference_concordance import ReferenceRoleBinding
    from genosapiens.submission_rendering import load_submission_contig_renderer
    from genosapiens.wgs_source_routing import SourceRoutingContext
    contig_renderer = load_submission_contig_renderer(
        identity_path=resources.population_reference.identity_path,
        reference_manifest_path=resources.population_reference.manifest_path,
        reference_root=resources.population_reference.reference_root,
        reference_binding=ReferenceRoleBinding(
            source_reference_role,
            reference_role.identity_id,
        ),
    )
    source_routing_context = SourceRoutingContext(native, contig_renderer)
    if analysis_input is not None:
        from genosapiens.wgs_analysis_inputs import resolve_completed_analysis_input
        with _completed_input_phase(phase_recorder,"completed_analysis_validation"):
            completed_input = resolve_completed_analysis_input(analysis_input,
                consumer_authority=authority, consumer_release=release,
                execution_config=config, source_routing_context=source_routing_context)
    prefix = None
    if prefix_path is not None:
        from genosapiens.wgs_annotation_prefix import load_completed_annotation_prefix, resolve_completed_annotation_prefix
        descriptor = load_completed_annotation_prefix(read_bounded_bytes(prefix_path), expected_sha256=prefix_sha)
        with _completed_input_phase(phase_recorder,"completed_annotation_validation"):
            prefix = resolve_completed_annotation_prefix(descriptor, plan=completed_input.transport.plan, execution_config=config)
            prefix.require_consumer(authority, release)
    completed_admission = None
    if admission_descriptor is not None:
        from genosapiens.wgs_admission_inputs import resolve_completed_admission_input
        with _completed_input_phase(phase_recorder,"completed_admission_validation"):
            completed_admission = resolve_completed_admission_input(admission_descriptor,
                consumer_authority=authority, consumer_release=release, completed_analysis_input=completed_input,
                annotation_prefix=prefix, execution_config=config, source_routing_context=source_routing_context)
    filter_policy=load_common_allele_filter_policy(_bound(FILTER_POLICY_PATH,
        authority.binding('common_allele_policy')),authority.binding('common_allele_policy'))
    precision_contract=load_gnomad_faf95_precision_contract(_bound(PRECISION_CONTRACT_PATH,
        authority.binding('gnomad_projection_precision_contract')),
        authority.binding('gnomad_projection_precision_contract'))
    precision=verify_gnomad_faf95_precision_qualification(_bound(PRECISION_REPORT_PATH,
        authority.binding('gnomad_projection_precision_report')),
        authority.binding('gnomad_projection_precision_report'),precision_contract)
    dependencies=[]
    for role,path,digest in (('original_source_manifest',source_root/'manifest.json',source.sha256),
            ('original_filter_manifest',filter_root/'manifest.json',filtered.sha256),
            ('original_working_manifest',working_root/'manifest.json',working.sha256)):
        dependencies.append(_manifest_dependency(role,path,digest))
    for role,path,digest in (('completed_analysis_input',supplied_path,supplied_sha),
            ('completed_annotation_prefix',prefix_path,prefix_sha)):
        if path is not None:
            content=_bound(Path(path),digest)
            dependencies.append((role,str(Path(path).absolute()),digest,len(content)))
    eligibility_arguments=dict(analysis_eligibility_policy=eligibility_policy,
        context_writer_profile=writer_profile,eligibility_common_policy=filter_policy,
        eligibility_precision=precision,eligibility_dependencies=tuple(sorted(dependencies)),
        eligibility_resource_bindings=tuple(sorted(dict(
            population_manifest=resources.population_joint.manifest_sha256,
            precision_manifest=precision.report_sha256,
            clinical_manifest=authority.binding('clinvar_manifest'),
            exception_manifest=authority.binding('common_allele_exception_manifest'),
            all_sites_manifest=resources.population_all_sites.bundle_manifest_sha256,
            public_index_manifest=indexes.manifest_sha256).items())))
    builder = PublicBoundedAnalysisBuilder(
        resources,
        policies,
        runtime,
        native,
        stack_sha256,
        authority.runtime_image_digest,
        _selected_sample(arguments, authority),
        indexes,
        annotation_prefix=prefix,
        public_pilot_inputs=public_pilot_inputs,
        focused_work_contract=focused_work,
        current_run=None if operator_case is None else operator_case.current_run,
        science_options=None if operator_case is None else operator_case.science_options,
        science_budget=None if operator_case is None else operator_case.science_budget,
        **eligibility_arguments,
    )
    from genosapiens.wgs_admission_diagnostics import analysis_observation
    with analysis_observation(Path(authority.output_root)/'analysis'):
        result = (analyze_bounded_case if executor is None else executor)(
            output_root=Path(authority.output_root) / "analysis",
            publication_root=publication_root,
            expected_manifest_sha256=prepared.sha256,
            expected_early_snapshot_sha256=early_snapshot_sha256,
            original_input_path=Path(authority.input_path),
            source_root=source_root,
            source_manifest=source,
            source_accounting=source_accounting,
            source_routing_context=source_routing_context,
            filter_root=filter_root,
            filter_manifest=filtered,
            working_root=working_root,
            working_manifest=working,
            policy_sha256=authority.binding("common_allele_policy"),
            population_manifest_sha256=authority.binding("gnomad_joint_manifest"),
            clinical_manifest_sha256=authority.binding("clinvar_manifest"),
            exception_manifest_sha256=authority.binding(
                "common_allele_exception_manifest"
            ),
            precision_manifest_sha256=authority.binding(
                "gnomad_projection_precision_report"
            ),
            case_manifest_sha256=source.source.case_manifest_sha256,
            selected_sample_sha256=source.source.selected_sample_sha256,
            hpo_input_sha256=prepared.hpo_input_sha256,
            phenotype_ingestion_sha256=phenotype_sha256,
            ontology_manifest_sha256=prepared.ontology_manifest_sha256,
            execution_config=config,
            base_release_bytes=release_bytes,
            expected_base_release_sha256=arguments.release_sha256,
            executable_stack_sha256=stack_sha256,
            runtime_image_digest=authority.runtime_image_digest,
            journals=journals,
            product_builder=builder,
            public_qualification_authority_bytes=authority_bytes,
            expected_public_qualification_authority_sha256=(
                arguments.authority_sha256 if operator_case is None else None
            ),
            phase_recorder=phase_recorder,
            completed_analysis_input=completed_input,
            completed_annotation_prefix=prefix,
            completed_admission_input=completed_admission,
            analysis_profile=mode,
            analysis_profile_policy_sha256=scope_policy.sha256,
            **({} if operator_case is None else dict(resume=arguments.resume,
                preparation_versions_bytes=operator_case.preparation_versions,
                report_revision=operator_case.current_run.publication_identity() if operator_case.current_run else None)),
        )
    print(result.receipt.status, result.receipt.terminal_stage, flush=True)
    return 0 if result.receipt.status == "passed" else 1


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("stage", choices=("prepare", "analyze"))
    parser.add_argument("--authority", type=Path, required=True)
    parser.add_argument("--authority-sha256", required=True)
    parser.add_argument("--release", type=Path, required=True)
    parser.add_argument("--release-sha256", required=True)
    parser.add_argument("--concordance", type=Path)
    parser.add_argument("--concordance-sha256")
    parser.add_argument("--selected-sample", help="Exact bound sample for fresh public/synthetic qualification")
    parser.add_argument("--verify-source-ledger", action="store_true")
    parser.add_argument("--exhaustive", action="store_true",
        help="analysis-only broader profile; unavailable until matching exhaustive qualification (default: focused)")
    parser.add_argument("--completed-analysis-input", type=Path)
    parser.add_argument("--completed-analysis-input-sha256")
    parser.add_argument("--completed-annotation-prefix", type=Path)
    parser.add_argument("--completed-annotation-prefix-sha256")
    parser.add_argument("--completed-admission-input", type=Path,
        help="authenticated current layout-v5 admission descriptor for fresh discovery; requires completed analysis and annotation inputs")
    parser.add_argument("--completed-admission-input-sha256",
        help="exact digest authorized under the completed_admission_input resource role")
    arguments = parser.parse_args(argv)
    if (arguments.completed_analysis_input is None) != (arguments.completed_analysis_input_sha256 is None):
        parser.error("completed analysis input path and digest must be supplied together")
    if ((arguments.completed_annotation_prefix is None) != (arguments.completed_annotation_prefix_sha256 is None)
            or arguments.completed_annotation_prefix is not None and arguments.completed_analysis_input is None):
        parser.error("completed annotation prefix requires its digest and completed analysis inputs")
    if ((arguments.completed_admission_input is None) != (arguments.completed_admission_input_sha256 is None)
            or arguments.completed_admission_input is not None and arguments.completed_annotation_prefix is None):
        parser.error("completed admission requires its digest and completed analysis/annotation inputs")
    if arguments.stage == "prepare":
        if arguments.exhaustive:
            parser.error("prepare does not accept --exhaustive")
        if arguments.completed_analysis_input is not None:
            parser.error("prepare does not accept completed analysis inputs")
        if arguments.concordance is None or arguments.concordance_sha256 is None:
            parser.error("prepare requires the exact concordance report")
        with managed_runtime_phase_recorder(
            expected_stage="preparation",
        ) as recorder:
            return _preparation(arguments, recorder)
    if arguments.concordance is not None or arguments.concordance_sha256 is not None:
        parser.error("analyze does not accept a concordance argument")
    if arguments.verify_source_ledger:
        parser.error("analyze does not accept a source audit option")
    from genosapiens.wgs_analysis_policy import select_analysis_profile
    select_analysis_profile(exhaustive=arguments.exhaustive, repository_root=REPOSITORY_ROOT)
    with managed_runtime_phase_recorder(expected_stage="analysis") as recorder:
        return _analysis(arguments, recorder)


if __name__ == "__main__":
    raise SystemExit(main())
