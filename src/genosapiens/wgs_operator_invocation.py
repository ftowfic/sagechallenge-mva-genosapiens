"""Bind a private operator case to an immutable installed package.

This boundary supplies data and authenticated configuration to the existing stage
assembly. It does not accept Python factories, scientific overrides or public
candidate-validation authority.
"""

from dataclasses import dataclass, replace
import argparse
from decimal import Decimal, InvalidOperation
import hashlib
import json
import os
from pathlib import Path
import re
import stat
import sys
from types import SimpleNamespace
from genosapiens.operator_warnings import friendly_operator_warnings

from genosapiens.bounded_json import load_bounded_json, read_bounded_bytes
from genosapiens.vcf_contracts import VcfReadLimits
from genosapiens.bounded_wgs_config import load_bounded_wgs_execution_config
from genosapiens.wgs_ledger import canonical
from genosapiens.wgs_preparation import PreparationRunOptions


from genosapiens.wgs_run_manifest import (SCHEMA as CASE_SCHEMA, VERSION as CASE_VERSION,
    load_current, resolved_inputs, atomic_current, installation_contracts,
    PREPARATION_CONTRACT, LEGACY_PREPARATION_MESSAGE)

_FIELDS = frozenset({
    'schema', 'version', 'classification', 'package_bindings_sha256', 'preparation_contract',
    'input_path', 'input_sha256', 'input_size_bytes', 'selected_sample',
    'hpo_input_path', 'hpo_input_sha256', 'hpo_input_size_bytes',
    'concordance_path', 'concordance_sha256', 'output_root', 'input_identity',
    'runtime_uid', 'runtime_gid', 'inputs', 'stages', 'contracts',
})
_SEAL = object()


def operator_run_parser():
    parser = argparse.ArgumentParser(prog='genosapiens track1 run',
        description='Run focused analysis offline using the installed package.')
    for name in ('vcf', 'hpo', 'output'):
        parser.add_argument('--' + name, type=Path, required=name == 'output')
    parser.add_argument('--sample')
    parser.add_argument('--concordance', type=Path)
    parser.add_argument('--resume', action='store_true')
    parser.add_argument('--allow-low-host-memory', action='store_true',
        help='Accept the risk of memory exhaustion below 26 GiB available host memory. Container limits still apply.')
    parser.add_argument('--max-hours', type=positive_hours,
        help='Optional positive hours for this entire invocation; omitted means no overall limit.')
    parser.add_argument('--science-mode', choices=('exhaustive', 'rapid_phenotype', 'rapid_strict_association'))
    parser.add_argument('--science-target-contexts', type=int, default=1500)
    parser.add_argument('--science-workers', type=int, choices=(1, 2, 4, 8), default=1)
    parser.add_argument('--science-publication-readers', type=int, choices=(1, 2, 4, 8), default=1)
    parser.add_argument('--science-hours', type=positive_hours, default=18000)
    parser.add_argument('--science-report-minutes', type=int, default=30)
    parser.add_argument('--science-pilot', action='store_true', help='Stop after 32 new contexts or the pilot deadline.')
    parser.add_argument('--science-pilot-minutes', type=int, choices=(20, 40), default=20,
        help='Explicitly select a 20-minute or 40-minute pilot. The five-hour deadline remains unchanged.')
    parser.add_argument('--science-no-budget', action='store_true',
        help='Measure elapsed time without an overall scientific deadline. Per-context limits remain active.')
    parser.add_argument('--science-fresh-trial', action='store_true',
        help='Authorize a fresh 30-minute timer for this revised candidate, including a two-minute reporting reserve.')
    return parser


def positive_hours(value):
    try:
        hours = Decimal(value)
        seconds = float(hours * 3600)
        import math
        if not hours.is_finite() or hours <= 0 or not math.isfinite(seconds) or seconds <= 0:
            raise ValueError
        return seconds
    except (ValueError, InvalidOperation, OverflowError):
        raise argparse.ArgumentTypeError('--max-hours requires positive finite hours, for example 24 or 0.5') from None


def validate_operator_inputs(vcf, hpo, sample, *, max_header_bytes, private_output):
    """Syntax-only intake; private choices never enter an exception or receipt."""
    from genosapiens.vcf_qc import read_vcf_samples, VcfQcError
    from genosapiens.phenotype_ingestion import _hpo_lines, PhenotypeIngestionError
    try:
        samples = read_vcf_samples(vcf, max_header_bytes=max_header_bytes)
    except (OSError, ValueError, VcfQcError):
        raise OperatorInvocationError('invalid VCF header; supply a complete supported VCF with unique nonempty sample names within the header size limit') from None
    if (sample is None and len(samples) != 1) or (sample is not None and sample not in samples):
        print('Private sample choices: ' + json.dumps(samples, ensure_ascii=True), file=private_output)
        raise OperatorInvocationError('sample selection required; rerun with --sample NAME from the private list')
    selected = samples[0] if sample is None else sample
    try:
        content = read_bounded_bytes(hpo)
        _hpo_lines(content)
    except PhenotypeIngestionError as error:
        # The shared parser exposes only static syntax and a line number.
        raise OperatorInvocationError(str(error)) from None
    except (OSError, ValueError):
        raise OperatorInvocationError('cannot read HPO text; use a readable UTF-8 file with one HP: identifier and seven digits per line') from None
    return selected, content


class OperatorInvocationError(ValueError):
    """Static diagnostic that must not expose case paths or contents."""


def _require(condition):
    if not condition:
        raise OperatorInvocationError('operator case or package binding is invalid')


def _sha(content):
    return hashlib.sha256(content).hexdigest()


def _digest(value):
    return type(value) is str and re.fullmatch('[0-9a-f]{64}', value) is not None


def _path(value, *, directory=False, private=True):
    _require(type(value) is str)
    path = Path(value)
    _require(path.is_absolute() and '..' not in path.parts
             and path.resolve(strict=True) == path)
    info = path.lstat()
    _require(info.st_uid == os.getuid())
    if directory:
        _require(stat.S_ISDIR(info.st_mode) and stat.S_IMODE(info.st_mode) == 0o700)
    else:
        _require(stat.S_ISREG(info.st_mode))
        if private:
            _require(stat.S_IMODE(info.st_mode) in (0o400, 0o600))
    return path


def _package(root, expected):
    """Read the bounded code/configuration closure once; never reference payloads."""
    from genosapiens.wgs_installation import authenticate_package
    return authenticate_package(root, expected)


@dataclass(frozen=True, slots=True)
class _Binding:
    role: str
    sha256: str


@dataclass(frozen=True, slots=True)
class OperatorCaseBinding:
    """An accepted operator request, never a public qualification authority."""
    input_path: str
    input_sha256: str
    input_size_bytes: int
    hpo_input_path: str
    hpo_input_sha256: str
    hpo_input_size_bytes: int
    output_root: str
    selected_sample: str
    resource_bindings: tuple
    policy_bindings: tuple
    contract: object
    release: object
    release_bytes: bytes
    config: object
    argument_values: tuple
    seal: object
    input_identity: tuple = ()
    invocation_id: str = 'operator-case'
    preparation_attempts: int = 1
    is_finite_profile_pilot: bool = False
    stage_contracts: object = None
    current_run: object = None
    preparation_versions: bytes | None = None
    science_options: object = None
    science_budget: object = None

    def for_stage(self, stage, arguments):
        from genosapiens.track1_installation import InstalledVersions
        _require(stage in ('preparation', 'analysis') and type(self.release) is InstalledVersions)
        key = 'preparation' if stage == 'preparation' else 'science'
        from genosapiens.wgs_science_options import computation_contract
        contract_values = computation_contract(self.stage_contracts) if stage == 'analysis' else self.stage_contracts[key]
        release = replace(self.release,
            executable_stack_sha256=_sha(canonical(contract_values)))
        contract = replace(self.contract, base_release_sha256=release.sha256,
            executable_stack_sha256=release.executable_stack_sha256)
        arguments = SimpleNamespace(**{**vars(arguments), 'release_sha256':release.sha256})
        return replace(self, release=release, release_bytes=release.as_bytes(), contract=contract,
            argument_values=tuple(sorted(vars(arguments).items()))), arguments

    @property
    def hpo_root(self):
        return Path(self.hpo_input_path).parent

    @property
    def completed_source_input(self):
        return None

    @property
    def runtime_image_digest(self):
        return self.contract.runtime_image_digest

    @property
    def preparation_run_options_sha256(self):
        return PreparationRunOptions(dict(self.argument_values)['verify_source_ledger']).sha256

    @property
    def classification_resume_authorized(self):
        return (self.seal is _SEAL
                and dict(self.argument_values).get('resume') is True)

    def binding(self, role):
        matches = tuple(row.sha256 for row in (*self.resource_bindings, *self.policy_bindings)
                        if row.role == role)
        _require(len(matches) == 1)
        return matches[0]

    def common(self, arguments):
        _require(self.seal is _SEAL and tuple(sorted(vars(arguments).items())) == self.argument_values
                )
        return (None, self, self.release_bytes, self.release,
                self.config, self.contract.executable_stack_sha256)

    def select_runtime(self, policies):
        _require(self.seal is _SEAL)
        matches = tuple(runtime for runtime in policies.runtimes
                        if runtime.candidate.bindings == self.release.policy_bindings)
        _require(len(matches) == 1)
        return matches[0]

    def require_prepared(self, source, prepared):
        """A re-pinned descriptor cannot change an already prepared case."""
        _require(self.seal is _SEAL)
        sample = _sha(canonical(dict(schema='genosapiens.operator-wgs-sample',
            sample_name=self.selected_sample, input_sha256=self.input_sha256)))
        case = _sha(canonical(dict(schema='genosapiens.operator-wgs-case',
            invocation_id=self.invocation_id, input_sha256=self.input_sha256,
            hpo_input_sha256=self.hpo_input_sha256)))
        _require(source.source.input_sha256 == self.input_sha256
                 and source.source.input_size_bytes == self.input_size_bytes
                 and source.source.selected_sample_sha256 == sample
                 and source.source.case_manifest_sha256 == case
                 and prepared.hpo_input_sha256 == self.hpo_input_sha256)

    def authenticated_source_identity(self, path, report):
        """Reuse the digest's observed file identity, checking it at consumption."""
        from genosapiens.vcf_qc import _stat_identity
        _require(self.seal is _SEAL and Path(path) == Path(self.input_path)
                 and report.validated_input_sha256 == self.input_sha256
                 and tuple(_stat_identity(Path(path).stat())) == self.input_identity)
        return self.input_identity


def load_operator_case(package_root, case_path, case_sha256, *, resume=False,
                       verify_source_ledger=False):
    """Authenticate only bounded metadata before the normal stage reads inputs."""
    try:
        _require(_digest(case_sha256) and type(resume) is bool
                 and type(verify_source_ledger) is bool)
        root = _path(str(package_root), directory=True)
        descriptor = _path(str(case_path))
        _path(str(descriptor.parent), directory=True)
        data, content = load_bounded_json(descriptor, reference_root=descriptor.parent,
                                         maximum_bytes=65536)
        _require(_sha(content) == case_sha256 and type(data) is dict and set(data) == _FIELDS
                 and data['schema'] == CASE_SCHEMA and type(data['version']) is int
                 and data['version'] == CASE_VERSION and data['classification'] == 'controlled'
                 and data['preparation_contract'] == PREPARATION_CONTRACT
                 and data['runtime_uid'] == os.getuid() > 0
                 and data['runtime_gid'] == os.getgid() > 0
                 and type(data['input_identity']) is list and len(data['input_identity']) == 5
                 and all(type(value) is int and value >= 0 for value in data['input_identity']))
        contract, release, release_bytes, files, read = _package(root, os.environ.get('GENOSAPIENS_PACKAGE_SHA256', data['package_bindings_sha256']))
        # This module must itself come from the authenticated package, not a
        # workspace import that happens to accept another package's manifest.
        _require(Path(__file__).resolve() == root / 'python/genosapiens/wgs_operator_invocation.py')
        from genosapiens.track1_analysis_profiles import authenticate_analysis_profile
        from genosapiens.wgs_analysis_policy import load_analysis_selectors, load_analysis_work_contract
        release.require_analysis_profile(authenticate_analysis_profile())
        for prefix in ('input', 'hpo_input'):
            _require(_digest(data[prefix + '_sha256'])
                     and type(data[prefix + '_size_bytes']) is int
                     and data[prefix + '_size_bytes'] > 0)
            path = _path(data[prefix + '_path'])
            _require(path.stat().st_size == data[prefix + '_size_bytes'])
        _path(str(Path(data['hpo_input_path']).parent), directory=True)
        concordance = _path(data['concordance_path'])
        _require(_digest(data['concordance_sha256']))
        output = _path(data['output_root'], directory=True)
        _require(not output.is_relative_to(root)
                 and type(data['selected_sample']) is str
                 and 0 < len(data['selected_sample']) <= 1024
                 and not any(char.isspace() or ord(char) < 32 for char in data['selected_sample']))
        config = load_bounded_wgs_execution_config(read('bounded-wgs-execution.json'),
                                                  contract.execution_limits_sha256)
        config.runtime.operator_identity(data['runtime_uid'], data['runtime_gid'])
        resources = dict(
            common_allele_policy=contract.common_allele_policy_sha256,
            gnomad_joint_manifest=contract.population_manifest_sha256,
            clinvar_manifest=contract.clinical_manifest_sha256,
            common_allele_exception_manifest=contract.exception_manifest_sha256,
            gnomad_projection_precision_report=contract.precision_manifest_sha256,
            common_allele_public_index_manifest=contract.common_allele_public_index_manifest_sha256,
            hpo_manifest=contract.ontology_manifest_sha256,
            analysis_eligibility_policy=contract.analysis_eligibility_policy_sha256,
            context_writer_profile=contract.context_writer_profile_sha256,
            analysis_profile_policy=contract.analysis_profile_policy_sha256,
            population_work_policy=contract.population_work_policy_sha256,
            focused_work_contract=contract.focused_work_contract_sha256,
            gnomad_projection_precision_contract=files[
                'config/public-resources/gnomad-joint-faf95-precision-v1.json']['sha256'])
        _require(set(resources.values()) <= set(release.resource_manifest_sha256s))
        arguments = SimpleNamespace(concordance=concordance,
            concordance_sha256=data['concordance_sha256'], release_sha256=release.sha256,
            selected_sample=data['selected_sample'], verify_source_ledger=verify_source_ledger,
            exhaustive=False, resume=resume)
        case = OperatorCaseBinding(**{name: data[name] for name in (
            'input_path', 'input_sha256', 'input_size_bytes', 'hpo_input_path',
            'hpo_input_sha256', 'hpo_input_size_bytes', 'output_root', 'selected_sample')},
            resource_bindings=tuple(_Binding(role, digest) for role, digest in sorted(resources.items())),
            policy_bindings=(), contract=contract, release=release, release_bytes=release_bytes,
            config=config, argument_values=tuple(sorted(vars(arguments).items())), seal=_SEAL)
        object.__setattr__(case, 'input_identity', tuple(data['input_identity']))
        contracts = installation_contracts(contract, release, files, read)
        contracts['preparation'].update(input=data['input_sha256'], hpo=data['hpo_input_sha256'],
            sample=_sha(canonical(data['selected_sample'])))
        object.__setattr__(case, 'stage_contracts', contracts)
        from genosapiens.vcf_qc import _stat_identity
        _require(_stat_identity(Path(case.input_path).stat()) == case.input_identity)
        load_analysis_selectors(case, repository_root=root)
        work = load_analysis_work_contract(case, repository_root=root)
        _require(work is not None)
        release.require_focused_work(work)
        _require(release.policy_bindings.bounded_wgs_execution_sha256 == contract.execution_limits_sha256
                 and release.policy_bindings.common_allele_policy_sha256 == contract.common_allele_policy_sha256)
        return case, arguments
    except OperatorInvocationError:
        raise
    except (OSError, ValueError, KeyError, TypeError, RuntimeError):
        raise OperatorInvocationError('operator case or package binding is invalid') from None


def run_operator_stage(stage, *, package_root, parsed, execute, phase_recorder):
    """Dispatch the same producer assembly through the repository-owned driver."""
    _require(stage in ('preparation', 'analysis'))
    case, arguments = load_operator_case(package_root, parsed.case_binding,
        parsed.case_binding_sha256, resume=parsed.resume,
        verify_source_ledger=getattr(parsed, 'verify_source_ledger', False))
    from genosapiens import wgs_public_runner as assembly
    _require(assembly.REPOSITORY_ROOT == Path(package_root))
    def executor(**values):
        if stage == 'preparation':
            values.update(default_preparation_run_options_sha256=case.contract.default_preparation_run_options_sha256,
                          audited_preparation_run_options_sha256=case.contract.audited_preparation_run_options_sha256)
        return execute(values)
    function = assembly._preparation if stage == 'preparation' else assembly._analysis
    return function(arguments, phase_recorder, operator_case=case, executor=executor)


def _save(root, name, content):
    from genosapiens.wgs_ledger import publish_bytes
    return publish_bytes(root, name, content)


def _receipt(path, kind):
    from genosapiens.evidence_model import _typed
    data, content = load_bounded_json(path)
    value = _typed(data, kind, decode=True)
    _require(value.as_bytes() == content)
    return value


def prepare_operator_concordance(vcf, supplied, *, destination, package_root,
                                 read_limits=VcfReadLimits()):
    """Produce required REF evidence once, or authenticate an explicit report."""
    from genosapiens import wgs_public_runner as assembly
    from genosapiens.reference_concordance import (
        ReferenceConcordanceReport, concordance_report_as_json,
        load_reference_identity, open_installed_reference,
        evaluate_vcf_reference_concordance,
    )
    from genosapiens.vcf_qc import _file_sha256, _stat_identity
    before = _stat_identity(vcf.stat())
    try:
        if supplied is None:
            print('Checking every source REF against the installed reference.', flush=True)
            contract = load_reference_identity(package_root / assembly.REFERENCE_IDENTITY_PATH)
            with open_installed_reference(contract,
                    identity_path=package_root / assembly.REFERENCE_IDENTITY_PATH,
                    reference_manifest_path=assembly.REFERENCE_MANIFEST,
                    reference_root=assembly.REFERENCE_ROOT,
                    reuse_published_validation=True) as (manifest, fetch, contigs):
                _require(manifest == assembly.REFERENCE_MANIFEST_SHA256)
                report = evaluate_vcf_reference_concordance(vcf, fetch, contigs,
                    contract=contract, reference_manifest_sha256=manifest, read_limits=read_limits)
            digest = report.validated_input_sha256
            content = concordance_report_as_json(report).encode('ascii')
        else:
            data, content = load_bounded_json(supplied)
            report = ReferenceConcordanceReport.from_mapping(data)
            _require(content == concordance_report_as_json(report).encode('ascii')
                     and report.read_limits == read_limits)
            # A report's declared digest is not authority for a fresh source.
            digest, observed = _file_sha256(vcf)
            _require(observed == before)
        _require(report.accepted and report.input_sha256 == report.validated_input_sha256 == digest
                 and report.reference_role == 'challenge_alignment'
                 and report.reference_manifest_sha256 == assembly.REFERENCE_MANIFEST_SHA256
                 and report.checked_records == report.matching_records == report.total_records
                 and report.unresolved_records == report.mismatching_records == 0
                 and _stat_identity(vcf.stat()) == before)
        _save(destination.parent, destination.name, content)
        return digest, before, _sha(content)
    except (OSError, ValueError, RuntimeError):
        raise OperatorInvocationError(
            'reference concordance failed; check the source reference declaration and complete REF evidence'
            if supplied is None else
            'supplied concordance is invalid or does not match this source/reference; supply its exact accepted complete report') from None


def _operator_case_setup(args, *, package_root, package_sha, hpo_content, selected):
    metadata = args.output / '.operator'
    descriptor = metadata / 'run.json'
    vcf = _path(str(args.vcf))
    if args.resume:
        if not descriptor.is_file():
            raise OperatorInvocationError('intake has no current input binding; supply the original inputs to rebuild intake')
        data, content = load_bounded_json(descriptor)
        case, stage_args = load_operator_case(package_root, descriptor, _sha(content), resume=True)
        _require(case.input_path == str(vcf) and case.selected_sample == selected
                 and case.output_root == str(args.output)
                 and case.hpo_input_sha256 == _sha(hpo_content))
        if args.concordance is not None:
            _require(_sha(read_bounded_bytes(args.concordance)) == data['concordance_sha256'])
        _require(_sha(read_bounded_bytes(Path(data['concordance_path']))) == data['concordance_sha256'])
        host_reference = os.environ.get('GENOSAPIENS_REFERENCE_HOST_PATH')
        if host_reference is not None:
            _require(data['inputs'].get('host_reference_root', host_reference) == host_reference)
            if 'host_reference_root' not in data['inputs']:
                data['inputs']['host_reference_root'] = host_reference
                atomic_current(descriptor, data)
                return load_operator_case(package_root, descriptor, _sha(canonical(data)), resume=True)
        return case, stage_args
    binding, versions, _, files, read = _package(package_root, package_sha)
    config = load_bounded_wgs_execution_config(read('bounded-wgs-execution.json'),
                                              binding.execution_limits_sha256)
    digest, identity, report_sha = prepare_operator_concordance(vcf, args.concordance,
        destination=metadata / 'concordance.json', package_root=package_root,
        read_limits=VcfReadLimits(config.input.max_header_bytes, config.input.max_record_bytes,
                                  config.input.max_record_bytes))
    _save(metadata, 'phenotypes.txt', hpo_content)
    content = canonical(dict(schema=CASE_SCHEMA, version=CASE_VERSION, classification='controlled',
        preparation_contract=PREPARATION_CONTRACT,
        inputs=dict(vcf=str(vcf), hpo=str(getattr(args, 'hpo', metadata / 'phenotypes.txt')),
            sample=selected, concordance=None if args.concordance is None else str(args.concordance),
            host_vcf=os.environ.get('GENOSAPIENS_INPUT_HOST_PATH', str(vcf)),
            host_hpo=os.environ.get('GENOSAPIENS_HPO_HOST_PATH', str(getattr(args, 'hpo', metadata / 'phenotypes.txt'))),
            host_concordance=os.environ.get('GENOSAPIENS_CONCORDANCE_HOST_PATH'),
            host_reference_root=os.environ.get('GENOSAPIENS_REFERENCE_HOST_PATH', '/mnt/reference')),
        stages={}, contracts=installation_contracts(binding, versions, files, read),
        package_bindings_sha256=package_sha, input_path=str(vcf), input_sha256=digest,
        input_size_bytes=identity[2], input_identity=identity, selected_sample=selected,
        hpo_input_path=str(metadata / 'phenotypes.txt'), hpo_input_sha256=_sha(hpo_content),
        hpo_input_size_bytes=len(hpo_content), concordance_path=str(metadata / 'concordance.json'),
        concordance_sha256=report_sha, output_root=str(args.output),
        runtime_uid=os.getuid(), runtime_gid=os.getgid()))
    _save(metadata, descriptor.name, content)
    return load_operator_case(package_root, descriptor, _sha(content))


def _recorded_operation(operation, arguments, recorder, *, stage_case, stage, workspace):
    """Finalize runtime telemetry without replacing a failed operation."""
    from contextlib import nullcontext
    from genosapiens.wgs_admission_diagnostics import Owner, Monitor, finish_terminal, receipt_binding
    owner = None
    monitor = None
    primary = None
    code = None
    if stage == 'analysis':
        owner = Owner(workspace/'analysis-work/state/admission-diagnostic.json')
        monitor = Monitor(owner.path, owner.value['invocation'], emit=lambda text: print(text, flush=True)).start()
    try:
        with owner.activate() if owner is not None else nullcontext():
            try:
                code = operation(arguments, recorder, operator_case=stage_case)
            except BaseException as error:
                primary = error
                raise
            finally:
                if primary is None and code and stage == 'analysis':
                    bound = receipt_binding(workspace)
                    if bound is not None:
                        value = bound[2]
                        primary = static_failure_from_record(value)
                try:
                    recorder.finish()
                except BaseException as error:
                    if primary is None and not code:
                        raise
                    if primary is not None:
                        from genosapiens.wgs_failures import run_cleanup
                        run_cleanup((('runtime_finalize', lambda: (_ for _ in ()).throw(error)),), primary=primary)
                    if owner is not None:
                        owner.secondary(error, 'runtime_finalize')
    finally:
        if monitor is not None:
            try:
                monitor.close()
            except BaseException as error:
                if primary is None and not code:
                    raise
                owner.secondary(error, 'monitor_close')
        if owner is not None:
            def emit_terminal(text):
                observed = getattr(owner, '_scientific_failure', None)
                if observed is not None:
                    subphase, failure, local_code, recovery = observed
                    text = ('Scientific failure: stage=scientific_evaluation, step=' + subphase
                        + ', operation=' + failure.operation + ', code=' + failure.failure_code
                        + ('' if local_code is None else ', local_code=' + local_code)
                        + ', recovery=' + recovery + '.')
                print(text, flush=True)
            finish_terminal(workspace, owner, emit=emit_terminal,
                failed=bool(code or primary))
    return code


def static_failure_from_record(value):
    """Construct only validated static fields from current failure evidence."""
    from genosapiens.wgs_failures import WgsFailure
    # Failure categories can contain underscores. Match the finite code set.
    from genosapiens.wgs_failures import ORIGINS, CATEGORIES
    origin, category = next((o, c) for o in ORIGINS for c in CATEGORIES
        if o + '_' + c == value['failure_code'])
    return WgsFailure(origin, category, operation=value['operation'], exception_kind=value['exception_kind'])


def execute_operator_run(args, *, package_root, package_sha, hpo_content, selected, attempt):
    """Internal orchestration of the existing accepted stage assembly."""
    from genosapiens import wgs_public_runner as assembly
    from genosapiens.wgs_runtime import RuntimePhaseRecorder
    from genosapiens.wgs_preparation import PreparationRunReceipt, PREPARATION_RECEIPT_NAME
    from genosapiens.wgs_analysis import AnalysisRunReceipt, ANALYSIS_RECEIPT_NAME
    from genosapiens.wgs_ledger import load_model, PreparedCaseManifest, MAX_MANIFEST_BYTES
    from genosapiens.track1_publication import check_published_inventory
    case, stage_args = _operator_case_setup(args, package_root=package_root,
        package_sha=package_sha, hpo_content=hpo_content, selected=selected)
    from genosapiens.wgs_run_manifest import CurrentRun, discard_invalid_work
    options = getattr(args, 'science_options', None)
    budget = None
    if options is not None:
        from genosapiens.wgs_science_options import run_identity
        budget = options.open_budget(args.output / '.operator', run_sha256=run_identity(
            input_sha256=case.input_sha256, hpo_sha256=case.hpo_input_sha256, sample=selected),
            started_ns=getattr(args, 'started_ns', None))
        case = replace(case, stage_contracts=options.bind_contracts(case.stage_contracts),
            science_options=options, science_budget=budget)
    run = CurrentRun(args.output / '.operator/run.json', case.stage_contracts)
    run.value['package_bindings_sha256'] = package_sha
    print('Resume stage: ' + (run.plan.restart_stage or 'complete') + '. ' + run.plan.reason + '.', flush=True)
    if args.resume:
        discard_invalid_work(run)
    else:
        run.save()
    case = replace(case, current_run=run)
    for stage, operation, receipt_type, receipt_name in (
            ('preparation', assembly._preparation, PreparationRunReceipt, PREPARATION_RECEIPT_NAME),
            ('analysis', assembly._analysis, AnalysisRunReceipt, ANALYSIS_RECEIPT_NAME)):
        if stage == 'analysis':
            from genosapiens.wgs_ledger import LedgerManifest
            from genosapiens.wgs_preparation import load_completed_source_input
            prepared = load_model(read_bounded_bytes(
                args.output / 'preparation/prepared/prepared-case-manifest.json',
                maximum_bytes=MAX_MANIFEST_BYTES), PreparedCaseManifest)
            source_input = load_completed_source_input(read_bounded_bytes(
                Path(prepared.completed_source_input_path), maximum_bytes=65536))
            _require(source_input.sha256 == prepared.completed_source_input.sha256)
            # Source manifests use the ledger's 8 MiB contract, not the generic
            # 2 MiB governing-JSON limit. The completed artifact also authenticates
            # its recorded size and digest without reopening source shards.
            source = load_model(source_input.source_manifest.read(), LedgerManifest)
            _require(source.sha256 == prepared.source_manifest.sha256)
            case.require_prepared(source, prepared)
        workspace = args.output / stage
        stage_resume = args.resume and workspace.exists() and any(workspace.iterdir())
        operation_args = SimpleNamespace(**{**vars(stage_args), 'resume': stage_resume})
        stage_case, operation_args = case.for_stage(stage, operation_args)
        if stage == 'analysis':
            stage_case = replace(stage_case,
                preparation_versions=read_bounded_bytes(args.output / '.operator/preparation-versions.json'))
        receipt_path = workspace / receipt_name
        if args.resume and ((stage == 'preparation' and run.reusable('preparation'))
                or stage == 'analysis' and run.reusable('reports')):
            receipt = _receipt(receipt_path, receipt_type)
            _require(receipt.status == 'passed')
            print(f'Reusing completed {stage}.', flush=True)
            continue
        if receipt_path.exists():
            previous = _receipt(receipt_path, receipt_type)
            failure_path = workspace / ('failure-' + stage + '.json')
            diagnostic = None
            if failure_path.exists():
                diagnostic, _ = load_bounded_json(failure_path)
            atomic_current(args.output / '.operator/last-error.json', dict(stage=stage,
                previous_status=previous.status, receipt=json.loads(previous.as_bytes()), diagnostic=diagnostic))
            receipt_path.unlink()
            if failure_path.exists():
                failure_path.unlink()
        workspace.mkdir(mode=0o700, exist_ok=args.resume)
        print('Preparing source calls and phenotypes.' if stage == 'preparation'
              else 'Annotating, evaluating focused hypotheses, and publishing reports.', flush=True)
        recorder = RuntimePhaseRecorder(stage=stage, workspace=attempt,
            attempt_id=_sha(attempt.name.encode()), execution_config_sha256=case.contract.execution_limits_sha256)
        from contextlib import nullcontext
        from genosapiens.wgs_science_runtime import bounded_operation
        if budget is not None and stage == 'analysis':
            atomic_current(attempt / 'science-publication-status.json', dict(status='running', budget_sha256=budget.sha256))
        try:
            with nullcontext() if budget is None else bounded_operation(budget,
                    maximum_seconds=max(0.001, budget.remaining(reporting=True)), reporting=True):
                code = _recorded_operation(operation, operation_args, recorder,
                    stage_case=stage_case, stage=stage, workspace=workspace)
        except BaseException:
            if budget is not None and stage == 'analysis':
                atomic_current(attempt / 'science-publication-status.json', dict(
                    status='incomplete_deadline' if budget.remaining(reporting=True) <= 0 else 'incomplete',
                    budget_sha256=budget.sha256))
            raise
        if code:
            if budget is not None and stage == 'analysis':
                atomic_current(attempt / 'science-publication-status.json', dict(
                    status='incomplete_deadline' if budget.remaining(reporting=True) <= 0 else 'incomplete',
                    budget_sha256=budget.sha256))
            return 1
        if stage == 'preparation':
            atomic_current(args.output / '.operator/preparation-versions.json', json.loads(stage_case.release_bytes))
            run.complete('preparation', receipt_path,
                workspace / 'prepared/prepared-case-manifest.json',
                args.output / '.operator/preparation-versions.json')
        else:
            publication = workspace / 'controlled-work/track1/runs' / (
                'bounded-' + prepared.sha256[:20] + '-' + run.publication_identity()[:20])
            outputs = [receipt_path, publication/'manifests/output-manifest.json',
                publication/'reports/preflight.json', publication/'reports/methods.md']
            if (publication/'submission/track1.csv').is_file():
                outputs.append(publication/'submission/track1.csv')
            run.complete('reports', *outputs)

    from contextlib import nullcontext
    from genosapiens.wgs_science_runtime import bounded_operation
    try:
        with nullcontext() if budget is None else bounded_operation(budget,
                maximum_seconds=max(0.001, budget.remaining(reporting=True)), reporting=True):
            prepared = load_model(read_bounded_bytes(
                args.output / 'preparation/prepared/prepared-case-manifest.json',
                maximum_bytes=MAX_MANIFEST_BYTES), PreparedCaseManifest)
            receipt = _receipt(args.output / 'analysis' / ANALYSIS_RECEIPT_NAME, AnalysisRunReceipt)
            _require(receipt.status == 'passed' and receipt.prepared_case_sha256 == prepared.sha256
                     and receipt.base_release_sha256 == stage_case.release.sha256)
            final = args.output / 'analysis/controlled-work/track1/runs' / (
                'bounded-' + prepared.sha256[:20] + '-' + run.publication_identity()[:20])
            check_path = args.output / '.operator/final-check.json'
            if args.resume and check_path.exists():
                check, _ = load_bounded_json(check_path)
                from genosapiens.track1_publication import load_final_output_manifest, OUTPUT_MANIFEST_PATH
                raw = read_bounded_bytes(final / OUTPUT_MANIFEST_PATH)
                manifest = load_final_output_manifest(raw)
                _require(check['output_manifest_sha256'] == _sha(raw) == receipt.output_manifest_sha256
                         and check['status'] == 'passed'
                         and check['competitive_status'] == manifest.competitive_status)
            else:
                print('Checking the final publication inventory and required outputs.', flush=True)
                check = check_published_inventory(final, expected_output_manifest_sha256=receipt.output_manifest_sha256)
                check['output_manifest_sha256'] = receipt.output_manifest_sha256
                _save(check_path.parent, check_path.name, canonical(check))
            _save(attempt, 'results.json', canonical(dict(publication_root=str(final),
                competitive_status=check['competitive_status'], output_manifest_sha256=receipt.output_manifest_sha256)))
            print('Publication complete: ' + check['competitive_status'] + '.', flush=True)
            for name in ('reports/preflight.json', 'reports/methods.md', 'submission/track1.csv'):
                if (final / name).exists():
                    host = os.environ.get('GENOSAPIENS_OUTPUT_HOST_PATH')
                    shown = Path(host) / (final / name).relative_to(args.output) if host else final / name
                    print(str(shown), flush=True)
            if budget is not None:
                atomic_current(attempt / 'science-publication-status.json', dict(status='complete', budget_sha256=budget.sha256))
            return 0 if check['competitive_status'] == 'ready' else 2
    except BaseException:
        if budget is not None:
            atomic_current(attempt / 'science-publication-status.json', dict(
                status='incomplete_deadline' if budget.remaining(reporting=True) <= 0 else 'incomplete',
                budget_sha256=budget.sha256))
        raise



@friendly_operator_warnings()
def run_operator_command(argv=None):
    import time
    import fcntl
    from genosapiens.wgs_runtime import require_bounded_cgroup_limits, supervise_operator_run, WgsRuntimeError
    started = time.monotonic()
    started_ns = time.time_ns()
    parser = operator_run_parser()
    args = parser.parse_args(argv)
    os.umask(0o077)
    try:
        args.output = args.output.absolute()
        if args.resume:
            from genosapiens.wgs_run_manifest import require_resume_archive_complete
            saved = load_current(args.output / '.operator/run.json')
            if saved.get('preparation_contract') != PREPARATION_CONTRACT:
                raise OperatorInvocationError(LEGACY_PREPARATION_MESSAGE)
            require_resume_archive_complete(args.output)
            resolved = resolved_inputs(saved, vars(args))
            args.vcf, args.hpo = Path(resolved['vcf']), Path(resolved['hpo'])
            args.sample = resolved['sample']
            args.concordance = None if resolved['concordance'] is None else Path(resolved['concordance'])
        elif args.vcf is None or args.hpo is None:
            parser.error('--vcf and --hpo are required for a fresh run')
        # Syntax and sample choices precede all input body/resource processing.
        selected, hpo_content = validate_operator_inputs(args.vcf, args.hpo, args.sample,
            max_header_bytes=8 * 1024**2, private_output=sys.stderr)
        args.vcf = args.vcf.resolve(strict=True)
        args.hpo = args.hpo.resolve(strict=True)
        if args.concordance is not None:
            args.concordance = args.concordance.resolve(strict=True)
        args.output = args.output.absolute()
        _require(args.output.resolve() == args.output and os.getuid() > 0 and os.getgid() > 0)
        root = Path(__file__).resolve().parents[2]
        expected = os.environ.get('GENOSAPIENS_PACKAGE_SHA256')
        if expected is None:
            raise OperatorInvocationError('installation is missing; install the delivered operator package before running a case')
        contract, release, _, _, read = _package(root, expected)
        _require(os.environ.get('GENOSAPIENS_RUNTIME_IMAGE', '').split('@')[-1] == contract.runtime_image_digest)
        require_bounded_cgroup_limits('/sys/fs/cgroup')
        import socket
        if (any(name != 'lo' for _, name in socket.if_nameindex())
                or not os.statvfs('/').f_flag & os.ST_RDONLY):
            raise OperatorInvocationError('use the offline Docker launcher with a read-only root filesystem')
        if not args.output.exists():
            if args.resume:
                raise OperatorInvocationError('resume requires an existing compatible output root')
            args.output.mkdir(mode=0o700)
        _path(str(args.output), directory=True)
        metadata = args.output / '.operator'
        # The host launcher owns only the empty scratch mount before intake.
        if not args.resume and any(p.name != '.scratch' for p in args.output.iterdir()):
            raise OperatorInvocationError('output already contains a run; use a fresh root or explicit --resume')
        metadata.mkdir(mode=0o700, exist_ok=args.resume)
        _path(str(metadata), directory=True)
        from genosapiens.bounded_wgs_config import load_bounded_wgs_execution_config
        config = load_bounded_wgs_execution_config(read('bounded-wgs-execution.json'), contract.execution_limits_sha256)
        import importlib.util
        spec = importlib.util.spec_from_file_location('_operator_capacity', root / 'check_bounded_host_capacity.py')
        capacity = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(capacity)
        from genosapiens.wgs_runtime import require_fixed_host_capacity
        host_available_memory = capacity._available_memory(Path('/proc/meminfo'))
        require_fixed_host_capacity(config,
            host_available_memory_bytes=host_available_memory,
            allow_low_host_memory=args.allow_low_host_memory)
        if (args.allow_low_host_memory
                and host_available_memory < config.runtime.host_available_memory_bytes):
            print('Warning: Available host memory is below 26 GiB. '
                  'You accepted the risk of memory exhaustion. The 24 GiB container limit still applies.',
                  file=sys.stderr, flush=True)
        with (metadata / 'invocation.lock').open('a') as lock:
            try:
                fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
            except BlockingIOError:
                raise OperatorInvocationError('this output already has an active invocation') from None
            # The child inherits the locked descriptor. A surviving child keeps
            # ownership even if its supervisor stops before cleanup.
            attempt = metadata / 'current'
            from genosapiens.wgs_run_manifest import preserve_previous_invocation
            # Keep original diagnostic evidence before creating this invocation.
            # A failed archive blocks another resume through its pending marker.
            preserve_previous_invocation(args.output)
            if attempt.exists():
                raise OperatorInvocationError('previous invocation evidence remains unarchived')
            attempt.mkdir(mode=0o700)
            _save(attempt, 'invocation.json', canonical(dict(schema='genosapiens.operator-invocation',
                version=2, runtime_uid=os.getuid(), runtime_gid=os.getgid(),
                maximum_seconds=args.max_hours, resume=args.resume, package_sha256=expected,
                allow_low_host_memory=args.allow_low_host_memory,
                host_available_memory_bytes=host_available_memory)))
            from genosapiens.wgs_science_options import ScienceRunOptions
            options = None
            if args.science_mode is not None:
                _require(float(args.science_hours).is_integer())
                options = ScienceRunOptions(args.science_mode, args.science_target_contexts,
                    int(args.science_hours), args.science_report_minutes * 60, args.science_pilot, args.science_pilot_minutes,
                    args.science_fresh_trial, args.science_no_budget,
                    args.science_workers, args.science_publication_readers)
            elif args.science_pilot or args.science_fresh_trial or args.science_no_budget:
                raise OperatorInvocationError('Scientific timer options require --science-mode')
            maximum_seconds = args.max_hours
            if options is not None:
                remaining = options.supervisor_seconds(metadata, elapsed_seconds=time.monotonic() - started)
                if remaining is not None:
                    bound = time.monotonic() - started + remaining
                    maximum_seconds = bound if maximum_seconds is None else min(maximum_seconds, bound)
            request = canonical(dict(science_options=None if options is None else options.as_dict(),
                started_ns=started_ns, vcf=str(args.vcf), hpo=str(args.hpo), output=str(args.output),
                concordance=None if args.concordance is None else str(args.concordance),
                sample=selected, resume=args.resume, hpo_sha256=_sha(hpo_content)))
            _save(attempt, 'request.json', request)
            child = (sys.executable, '-c',
                'from genosapiens.wgs_operator_invocation import _execute_request; '
                'import sys; sys.exit(_execute_request(sys.argv[1], sys.argv[2]))',
                str(attempt / 'request.json'), _sha(request))
            return supervise_operator_run(child, workspace=attempt, started=started,
                maximum_seconds=maximum_seconds, kill_grace_seconds=config.child_process.kill_grace_seconds,
                lock_fd=lock.fileno())
    except OperatorInvocationError as error:
        print(str(error), file=sys.stderr)
        return 1
    except WgsRuntimeError as error:
        print(str(error) + '. Check Docker limits and available memory.', file=sys.stderr)
        return 1
    except (OSError, ValueError, RuntimeError):
        print('operator prerequisite or binding failed; check the installed package, private file permissions and output root.', file=sys.stderr)
        return 1


@friendly_operator_warnings(suppress_hpo=True)
def _execute_request(path, digest):
    """Private child entry; no operator-facing alternate command or overrides."""
    import signal
    try:
        os.umask(0o077)
        for sig in (signal.SIGINT, signal.SIGTERM):
            signal.signal(sig, signal.default_int_handler)
        request = _path(path)
        data, content = load_bounded_json(request, maximum_bytes=65536)
        _require(_sha(content) == digest and set(data) == {
            'vcf', 'hpo', 'output', 'concordance', 'sample', 'resume', 'hpo_sha256', 'science_options', 'started_ns'})
        selected, hpo = validate_operator_inputs(Path(data['vcf']), Path(data['hpo']), data['sample'],
            max_header_bytes=8 * 1024**2, private_output=sys.stderr)
        _require(_sha(hpo) == data['hpo_sha256'])
        from genosapiens.wgs_science_options import ScienceRunOptions
        options = None if data['science_options'] is None else ScienceRunOptions(**data['science_options'])
        args = SimpleNamespace(science_options=options, started_ns=data['started_ns'], vcf=Path(data['vcf']), hpo=Path(data['hpo']), output=Path(data['output']),
            concordance=None if data['concordance'] is None else Path(data['concordance']),
            resume=data['resume'])
        return execute_operator_run(args, package_root=Path(__file__).resolve().parents[2],
            package_sha=os.environ['GENOSAPIENS_PACKAGE_SHA256'], hpo_content=hpo,
            selected=selected, attempt=request.parent)
    except OperatorInvocationError as error:
        print(str(error), file=sys.stderr, flush=True)
        return 1
    except BaseException:
        print('Case stopped; inspect the private stage status and resume only with matching completed inputs.',
              file=sys.stderr, flush=True)
        return 1
