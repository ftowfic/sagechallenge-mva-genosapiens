"""Current-roster Exomiser comparison with independently replayed raw outputs.

The baseline receives exact blinded public inputs and its previously fixed
analysis. Recovery misses never become execution failures or a superiority gate.
"""
from __future__ import annotations

from dataclasses import asdict
from fractions import Fraction
import hashlib
import json
import os
from pathlib import Path
import time

from genosapiens.track1_exomiser import (
    load_exomiser_stack, render_truth_free_phenopacket, exomiser_command,
    parse_exomiser_outputs, outputs_sha256, pheval_input_as_tsv,
    compare_exomiser_outputs,
)
from genosapiens.track1_public_validation import materialize_runner_view
from genosapiens.track1_public_worker import _engine
from genosapiens.track1_validation_plan import InputKind
from genosapiens.track1_validation_runner import (
    authenticated_held_out_ranker_view, development_ranker_view,
)

SCHEMA = 'genosapiens.track1-exomiser-roster-comparison'
MANIFEST = 'config/public-resources/exomiser-15.1.0-2512.json'
UNSUPPORTED = frozenset({'mixed_class', 'mosaic', 'regulatory', 'repeat', 'sv_cnv'})
REPORT = 'exomiser/current-roster-comparison.json'


def _sha(content):
    return hashlib.sha256(content).hexdigest()


def _bytes(value):
    return (json.dumps(value, sort_keys=True, separators=(',', ':'), allow_nan=False) + '\n').encode('ascii')


def _file_identity(path):
    if path.is_symlink() or any(parent.is_symlink() for parent in path.parents):
        raise ValueError
    digest = hashlib.sha256()
    with path.open('rb') as source:
        before = os.fstat(source.fileno())
        while chunk := source.read(4 * 1024**2):
            digest.update(chunk)
        after = os.fstat(source.fileno())
    if (before.st_dev, before.st_ino, before.st_size, before.st_mtime_ns, before.st_ctime_ns) != (
        after.st_dev, after.st_ino, after.st_size, after.st_mtime_ns, after.st_ctime_ns
    ):
        raise ValueError
    return digest.hexdigest(), before.st_size


def authenticate_installed_baseline(repository):
    stack = load_exomiser_stack(MANIFEST, repository_root=repository)
    configuration = json.loads((repository / MANIFEST).read_bytes())
    manifest = Path(configuration['data']['installed_file_manifest_path'])
    content = manifest.read_bytes()
    if _sha(content) != stack.installed_file_manifest_sha256 or len(content) > 512 * 1024:
        raise ValueError
    root = Path(stack.installed_root).resolve(strict=True)
    identities = set()
    total = 0
    for line in content.decode('ascii').splitlines():
        digest, name = line.split(maxsplit=1)
        path = Path(name)
        if not path.is_absolute() or path != path.resolve(strict=True) or not path.is_relative_to(root) or path in identities:
            raise ValueError
        actual, size = _file_identity(path)
        if actual != digest:
            raise ValueError
        identities.add(path)
        total += size
    if len(identities) != stack.installed_file_count or total != stack.installed_bytes:
        raise ValueError
    return stack


def _blinded_inputs(case, outcome, roster, resources):
    if case.split == 'held_out':
        view = authenticated_held_out_ranker_view(
            outcome.candidate, outcome.development_experiment, roster, case.case_id,
            runner_resources=resources.common_runner_resources,
        )
    else:
        view = development_ranker_view(roster, case.case_id, runner_resources=resources.common_runner_resources)
    materialized = materialize_runner_view(view, resources)
    content = {descriptor.kind: (path, descriptor.sha256, size)
               for descriptor, path, size in zip(view.inputs, materialized.input_paths, materialized.input_sizes, strict=True)}
    return view, content


def baseline_container_spec(stack, runner, output):
    if os.getuid() == 0:
        raise ValueError
    return {
        'Image': stack.image, 'User': '65534:65534',
        'Cmd': list(exomiser_command(output_name='case')),
        'Env': ['EXOMISER_HG38_DATA_VERSION=2512', 'EXOMISER_PHENOTYPE_DATA_VERSION=2512', 'JAVA_TOOL_OPTIONS=-Xmx12g'],
        'Labels': {'org.genosapiens.purpose': 'public-validation-exomiser-case'},
        'HostConfig': {
            'NetworkMode': 'none', 'ReadonlyRootfs': True, 'CapDrop': ['ALL'],
            'SecurityOpt': ['no-new-privileges'], 'Memory': 16 * 1024**3,
            'NanoCpus': 4_000_000_000, 'PidsLimit': 1024,
            'GroupAdd': [str(os.getgid())],
            'Tmpfs': {'/tmp': 'rw,nosuid,nodev,noexec,size=1g'},
            'Mounts': [
                {'Type': 'bind', 'Source': stack.installed_root, 'Target': '/exomiser-data', 'ReadOnly': True},
                {'Type': 'bind', 'Source': str(runner), 'Target': '/inputs', 'ReadOnly': True},
                {'Type': 'bind', 'Source': str(output), 'Target': '/results', 'ReadOnly': False},
            ],
        },
    }


def _adapt(directory):
    genes = tuple(directory.glob('*.genes.tsv'))
    variants = tuple(directory.glob('*.variants.tsv'))
    if len(genes) != 1 or len(variants) != 1:
        raise ValueError
    gene_sha, variant_sha = _file_identity(genes[0])[0], _file_identity(variants[0])[0]
    outputs = parse_exomiser_outputs(genes[0], variants[0])
    if gene_sha != _file_identity(genes[0])[0] or variant_sha != _file_identity(variants[0])[0]:
        raise ValueError
    identities = {'raw_gene_sha256': gene_sha,
                  'raw_variant_sha256': variant_sha,
                  'adapted_sha256': outputs_sha256(outputs),
                  'pheval_gene_sha256': _sha(pheval_input_as_tsv(outputs, kind='gene')),
                  'pheval_variant_sha256': _sha(pheval_input_as_tsv(outputs, kind='variant'))}
    return outputs, identities


def _attempt(stack, runner, output):
    container = None
    started = time.monotonic()
    cpu = memory = None
    status = 'failed'
    exit_code = None
    spec = baseline_container_spec(stack, runner, output)
    try:
        container = _engine('POST', '/containers/create', spec)['Id']
        _engine('POST', f'/containers/{container}/start')
        while True:
            state = _engine('GET', f'/containers/{container}/json')['State']
            if not state['Running']:
                exit_code = state['ExitCode']
                status = 'succeeded' if exit_code == 0 else 'failed'
                break
            if time.monotonic() - started > 900:
                status = 'timed_out'
                break
            stats = _engine('GET', f'/containers/{container}/stats?stream=false', timeout=10)
            observed_cpu = stats.get('cpu_stats', {}).get('cpu_usage', {}).get('total_usage')
            observed_memory = stats.get('memory_stats', {}).get('usage')
            if type(observed_cpu) is int:
                cpu = max(cpu or 0, observed_cpu)
            if type(observed_memory) is int:
                memory = max(memory or 0, observed_memory)
    except Exception:
        status = 'failed'
    finally:
        if container is not None:
            try:
                _engine('POST', f'/containers/{container}/kill')
            except Exception:
                pass
            _engine('DELETE', f'/containers/{container}?force=true')
    record = {'status': status, 'exit_code': exit_code, 'wall_seconds': str(time.monotonic() - started),
              'container_spec_sha256': _sha(_bytes(spec)),
              'sampled_cpu_nanoseconds': cpu, 'sampled_peak_memory_bytes': memory,
              'sampling_semantics': 'observed cgroup samples; CPU is a lower bound, memory is a sampled peak; unobserved values remain null',
              'scratch_bytes': sum(path.stat().st_size for path in output.rglob('*') if path.is_file())}
    outputs = None
    if status == 'succeeded':
        try:
            outputs, identities = _adapt(output)
            record.update(identities)
            record['status'] = 'succeeded' if outputs.genes or outputs.variants else 'succeeded-empty'
        except Exception:
            record['status'] = 'failed-output-validation'
    return record, outputs


def run_current_roster_baseline(outcome, *, repository_root, roster, resources, checkpoint_root):
    from genosapiens.track1_public_validation_coordinator import publish_validation_artifact, publish_validation_file
    repository, journal = Path(repository_root), Path(checkpoint_root)
    if outcome.candidate is None or outcome.status != 'held-out-complete' or outcome.candidate.roster_sha256 != roster.roster_sha256:
        raise ValueError
    report_path = journal / REPORT
    stack = authenticate_installed_baseline(repository)
    if report_path.exists():
        return verify_current_roster_baseline(report_path.read_bytes(), outcome=outcome, roster=roster,
                                             repository_root=repository, checkpoint_root=journal)
    records = []
    for case in roster.cases:
        view, contents = _blinded_inputs(case, outcome, roster, resources)
        record = {'case_id': case.case_id, 'split': case.split, 'scope': case.scope,
                  'input_sha256s': [value.sha256 for value in view.inputs], 'attempts': []}
        if case.scope in UNSUPPORTED:
            record.update(status='unsupported', comparison=None)
            records.append(record)
            continue
        case_root = journal / 'exomiser/cases' / case.case_id
        runner = case_root / 'runner'
        runner.mkdir(parents=True, mode=0o750)
        runner.chmod(0o750)
        path, digest, size = contents[InputKind.VCF]
        publish_validation_file(runner / "case.vcf", path, expected_sha256=digest, expected_size_bytes=size)
        (runner / "case.vcf").chmod(0o640)
        from genosapiens.bounded_json import read_bounded_bytes
        hpo = read_bounded_bytes(contents[InputKind.HPO_TEXT][0], maximum_bytes=2 * 1024 * 1024)
        if _sha(hpo) != contents[InputKind.HPO_TEXT][1]:
            raise ValueError
        for name, content in (
            ('sample.json', render_truth_free_phenopacket(hpo, case_name='public-input')),
            ('analysis.yml', (repository / stack.analysis_path).read_bytes()),
        ):
            publish_validation_artifact(runner / name, content)
            (runner / name).chmod(0o640)
        outputs = []
        for number in (1, 2):
            directory = case_root / ('attempt-' + str(number))
            directory.mkdir(mode=0o770)
            directory.chmod(0o770)
            attempt, adapted = _attempt(stack, runner, directory)
            attempt.update(attempt_number=number, artifact_directory=directory.relative_to(journal).as_posix())
            publish_validation_artifact(case_root / ('attempt-' + str(number) + '.json'), _bytes(attempt))
            record['attempts'].append(attempt)
            outputs.append(adapted)
        succeeded = all(value is not None for value in outputs)
        equal = succeeded and outputs_sha256(outputs[0]) == outputs_sha256(outputs[1])
        record['status'] = 'succeeded' if equal else 'failed'
        record['comparison'] = asdict(compare_exomiser_outputs(outputs[0], truth_members=case.truth_members,
                                                              expected_genes=case.expected_gene_symbols)) if succeeded else None
        records.append(record)
        publish_validation_artifact(case_root / 'comparison.json', _bytes(record))
        print('completed current-roster Exomiser case comparison', flush=True)
    authenticate_installed_baseline(repository)
    report = {'schema': SCHEMA, 'classification': 'public', 'validation_only': True, 'controlled_sources_used': False,
              'candidate_sha256': outcome.candidate.sha256, 'roster_sha256': roster.roster_sha256,
              'executable_stack_sha256': outcome.executable_stack_sha256, 'stack': asdict(stack),
              'records': records, 'full_roster_denominator': len(roster.cases),
              'supported_denominator': sum(case.scope not in UNSUPPORTED for case in roster.cases),
              'status': 'passed' if all(record['status'] in {'succeeded', 'unsupported'} for record in records) else 'failed',
              'comparison_gate': 'report-only; no superiority or non-inferiority release gate',
              'limitations': ['Shared HG002 synthetic probes are not independent patients or clinical calibration.',
                              'Unsupported scopes and empty-output recovery misses retain their complete denominators.',
                              'Each supported case uses its own exact blinded HPO/VCF inputs; no aggregate split run substitutes for a case.']}
    content = _bytes(report)
    publish_validation_artifact(report_path, content)
    return verify_current_roster_baseline(content, outcome=outcome, roster=roster,
                                         repository_root=repository, checkpoint_root=journal)


def verify_current_roster_baseline(content, *, outcome, roster, repository_root, checkpoint_root):
    report = json.loads(content)
    stack = load_exomiser_stack(MANIFEST, repository_root=repository_root)
    if (not isinstance(report, dict) or set(report) != {
            'schema', 'classification', 'validation_only', 'controlled_sources_used',
            'candidate_sha256', 'roster_sha256', 'executable_stack_sha256', 'stack',
            'records', 'full_roster_denominator', 'supported_denominator', 'status',
            'comparison_gate', 'limitations',
        } or report.get('schema') != SCHEMA or report.get('classification') != 'public'
            or report.get('validation_only') is not True or report.get('controlled_sources_used') is not False
            or report.get('candidate_sha256') != outcome.candidate.sha256
            or report.get('roster_sha256') != roster.roster_sha256
            or report.get('executable_stack_sha256') != outcome.executable_stack_sha256
            or report.get('stack') != asdict(stack) or report.get('full_roster_denominator') != len(roster.cases)
            or type(report.get('full_roster_denominator')) is not int
            or type(report.get('supported_denominator')) is not int
            or report.get('supported_denominator') != sum(case.scope not in UNSUPPORTED for case in roster.cases)
            or report.get('comparison_gate') != 'report-only; no superiority or non-inferiority release gate'):
        raise ValueError
    records = report.get('records')
    cases = {case.case_id: case for case in roster.cases}
    if not isinstance(records, list) or len(records) != len(cases) or {record.get('case_id') for record in records} != set(cases):
        raise ValueError
    for record in records:
        if set(record) != {'case_id', 'split', 'scope', 'input_sha256s', 'attempts', 'status', 'comparison'}:
            raise ValueError
        case = cases[record['case_id']]
        if (record['split'], record['scope']) != (case.split, case.scope):
            raise ValueError
        if record['input_sha256s'] != [item.sha256 for item in case.runner_view.inputs]:
            raise ValueError
        if case.scope in UNSUPPORTED:
            if record['status'] != 'unsupported' or record['attempts'] or record['comparison'] is not None:
                raise ValueError
            continue
        if len(record['attempts']) != 2:
            raise ValueError
        case_root = Path(checkpoint_root) / 'exomiser/cases' / case.case_id
        runner = case_root / 'runner'
        source_roster = json.loads((Path(repository_root) / 'tests/fixtures/track1/roster/roster.json').read_bytes())
        source_case = next(value for value in source_roster['cases'] if value['case_id'] == case.case_id)
        expected_inputs = {value['role']: value for value in source_case['runner_inputs']}
        hpo = Path(repository_root) / expected_inputs['hpo_text']['path']
        if _file_identity(hpo)[0] != expected_inputs['hpo_text']['sha256']:
            raise ValueError
        if (_file_identity(runner / 'case.vcf')[0] != expected_inputs['small_variant_vcf']['sha256']
                or _file_identity(runner / 'sample.json')[0] != _sha(render_truth_free_phenopacket(hpo.read_bytes(), case_name='public-input'))
                or _file_identity(runner / 'analysis.yml')[0] != stack.analysis_sha256):
            raise ValueError
        outputs = []
        for number, attempt in enumerate(record['attempts'], 1):
            common_fields = {'status', 'exit_code', 'wall_seconds', 'sampled_cpu_nanoseconds',
                             'sampled_peak_memory_bytes', 'sampling_semantics', 'scratch_bytes',
                             'container_spec_sha256', 'attempt_number', 'artifact_directory'}
            success_fields = {'raw_gene_sha256', 'raw_variant_sha256', 'adapted_sha256',
                              'pheval_gene_sha256', 'pheval_variant_sha256'}
            succeeded_attempt = attempt.get('status') in {'succeeded', 'succeeded-empty'}
            if (set(attempt) != common_fields | (success_fields if succeeded_attempt else set())
                    or attempt.get('status') not in {'succeeded', 'succeeded-empty', 'failed', 'timed_out', 'failed-output-validation'}
                    or type(attempt['wall_seconds']) is not str or Fraction(attempt['wall_seconds']) <= 0
                    or type(attempt['scratch_bytes']) is not int or attempt['scratch_bytes'] < 0
                    or attempt['exit_code'] is not None and type(attempt['exit_code']) is not int
                    or any(value is not None and (type(value) is not int or value < 0)
                           for value in (attempt['sampled_cpu_nanoseconds'], attempt['sampled_peak_memory_bytes']))
                    or succeeded_attempt and attempt['exit_code'] != 0):
                raise ValueError
            relative = Path('exomiser/cases') / case.case_id / ('attempt-' + str(number))
            if attempt['attempt_number'] != number or attempt['artifact_directory'] != relative.as_posix():
                raise ValueError
            directory = Path(checkpoint_root) / relative
            if attempt.get('container_spec_sha256') != _sha(_bytes(baseline_container_spec(stack, runner, directory))):
                raise ValueError
            if attempt['status'] in {'succeeded', 'succeeded-empty'}:
                adapted, identities = _adapt(directory)
                if any(attempt.get(key) != value for key, value in identities.items()):
                    raise ValueError
                if (attempt['status'] == 'succeeded-empty') is not (not adapted.genes and not adapted.variants):
                    raise ValueError
                outputs.append(adapted)
            else:
                outputs.append(None)
        succeeded = all(value is not None for value in outputs)
        equal = succeeded and outputs_sha256(outputs[0]) == outputs_sha256(outputs[1])
        if record['status'] != ('succeeded' if equal else 'failed'):
            raise ValueError
        comparison = asdict(compare_exomiser_outputs(outputs[0], truth_members=case.truth_members,
                                                     expected_genes=case.expected_gene_symbols)) if succeeded else None
        if json.loads(_bytes(comparison)) != record['comparison']:
            raise ValueError
    if report['status'] != ('passed' if all(record['status'] in {'succeeded', 'unsupported'} for record in records) else 'failed'):
        raise ValueError
    return report
