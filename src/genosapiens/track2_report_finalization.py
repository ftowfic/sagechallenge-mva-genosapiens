"""Recover checked synthesis and publish operator completion under the host lock."""
from pathlib import Path

from genosapiens.track2_operator_inputs import small_json
from genosapiens.track2_report_metadata import checked_report_metadata
from genosapiens.track2_run_state import RunState, atomic_record
from genosapiens.track2_scientific_storage import file_identity

PRODUCTION_SCHEMA = 'genosapiens.track2-report-production/v1'
VALIDATORS = frozenset({
    ('genosapiens.track2_synthesis', 'validate_research_report'),
    ('genosapiens.track2_synthesis', 'validate_study_report'),
    ('genosapiens.track2_discovery_report', 'validate_native_report'),
    ('genosapiens.track2_report_metadata', '*'),
    ('genosapiens.track2_report_checks', '*'),
})


class RecoveryError(ValueError):
    """A static recovery reason without source contents or private paths."""
    def __init__(self, reason):
        self.reason = reason
        super().__init__('Report recovery failed: '+reason+'.')


def recovery_check(condition, reason):
    if not condition:
        raise RecoveryError(reason)


def production_binding(root, package, *, native):
    from genosapiens.track2_construction_dependencies import python_dependencies
    root, package = Path(root), Path(package)
    binding = dict(snapshot_sha256=file_identity(root/'preparation/snapshot.json')['sha256'],
        code=python_dependencies(package, [
            ('genosapiens.track2_discovery_report', 'export_native_evidence') if native else
            ('genosapiens.track2_synthesis', 'export_study_evidence'),
            ('genosapiens.track2_report_fresh', 'run_fresh_report')], ignored=VALIDATORS))
    from genosapiens.track2_report_gene_metadata import metadata_binding
    binding['gene_metadata'] = metadata_binding(package)
    binding['overview_renderer'] = file_identity(package/'scripts/render_track2_challenge_overviews.R')['sha256']
    if native:
        from genosapiens.track2_discovery_outcome import checked_outcome
        from genosapiens.track2_lincs_report import native_render_binding
        from genosapiens.track2_scientific_completion import check_completion
        checked_outcome(root/'discovery', plan_path=root/'plan/plan.json')
        check_completion(root/'discovery-figures', stage='figures',
            binding=native_render_binding(root/'discovery'))
        binding.update(discovery_sha256=file_identity(root/'discovery/discovery-results.json')['sha256'],
            figures_completion_sha256=file_identity(root/'discovery-figures/completion.json')['sha256'],
            discovery_completion_sha256=file_identity(root/'discovery/completion.json')['sha256']
                if (root/'discovery/completion.json').exists() else None,
            renderer=file_identity(package/'scripts/render_track2_lincs.R')['sha256'])
    return binding


def production_matches(saved, current):
    """Old routing can contain validation dependencies in addition to production.

    Every current reachable production definition must have its original hash.
    Extra old validation definitions cannot replace a missing production entry.
    """
    return (isinstance(saved, dict) and isinstance(saved.get('code'), dict)
        and {k:v for k,v in saved.items() if k != 'code'} == {k:v for k,v in current.items() if k != 'code'}
        and all(saved['code'].get(k) == v for k,v in current['code'].items()))


def retain_production_routing(root):
    """Preserve existing current-format proof before review replaces routing."""
    root = Path(root)
    report = root/'synthesis/controlled-work/report/manifest.json'
    result = root/'synthesis/controlled-work/result.json'
    route = root/'current-routing.json'
    destination = root/'retained-report-routing.json'
    if not report.is_file() or not result.is_file() or not route.is_file():
        return
    manifest_hash = file_identity(report)['sha256']
    try:
        retained = small_json(destination) if destination.exists() else {}
    except ValueError:
        retained = {}
    if isinstance(retained, dict) and retained.get('manifest_sha256') == manifest_hash:
        return
    try:
        routing = small_json(route)
    except ValueError:
        return
    if isinstance(routing, dict) and 'report_binding' in routing:
        atomic_record(destination, dict(schema=PRODUCTION_SCHEMA, binding=routing['report_binding'],
            manifest_sha256=manifest_hash, result_sha256=file_identity(result)['sha256']))


def save_production(synthesis, binding):
    """Retain original producer dependencies after the synthesis result is durable."""
    synthesis = Path(synthesis)
    atomic_record(synthesis/'production.json', dict(schema=PRODUCTION_SCHEMA, binding=binding,
        manifest_sha256=file_identity(synthesis/'controlled-work/report/manifest.json')['sha256'],
        result_sha256=file_identity(synthesis/'controlled-work/result.json')['sha256']))


def checked_synthesis(root, binding):
    root = Path(root)
    synthesis = root/'synthesis'
    proof_path = synthesis/'production.json'
    if not proof_path.exists():
        proof_path = root/'retained-report-routing.json'
    recovery_check(proof_path.is_file(), 'production_proof_missing')
    proof = small_json(proof_path)
    recovery_check(set(proof) == {'schema', 'binding', 'manifest_sha256', 'result_sha256'}
        and proof['schema'] == PRODUCTION_SCHEMA and production_matches(proof['binding'], binding),
        'production_dependencies_changed_or_incomplete')
    from genosapiens.track2_scientific_reader import read_scientific
    report=synthesis/'controlled-work/report'
    if read_scientific(report/'manifest.json').get('schema')=='genosapiens.track2-challenge-package/v2':
        from genosapiens.track2_report_fresh import coordinator_metadata
        metadata=coordinator_metadata(report,expected_manifest_sha256=proof['manifest_sha256'])
    else:
        metadata = checked_report_metadata(report,expected_manifest_sha256=proof['manifest_sha256'])
    result_path = synthesis/'controlled-work/result.json'
    recovery_check(result_path.is_file(), 'saved_synthesis_result_missing')
    result = small_json(result_path)
    recovery_check(set(result) == set(metadata.readiness) | {'manifest_sha256'} |
            ({'output'} if metadata.kind == 'study' else set())
        and proof['manifest_sha256'] == metadata.manifest_sha256
        and proof['result_sha256'] == file_identity(result_path)['sha256']
        and result.get('manifest_sha256') == metadata.manifest_sha256
        and all(result.get(k) == v for k,v in metadata.readiness.items()),
        'saved_synthesis_result_changed_or_inconsistent')
    if metadata.kind=='study' and metadata.input_snapshot_sha256 is None:
        from genosapiens.track2_freeze import load_frozen_public_envelope
        from genosapiens.track2_scientific_reader import read_scientific
        frozen=load_frozen_public_envelope(root/'preparation/snapshot.json')
        manifest=read_scientific(synthesis/'controlled-work/report/manifest.json')
        if manifest.get('schema')=='genosapiens.track2-challenge-package/v2':
            from genosapiens.track2_report_evidence import scientific_root
            evidence=synthesis/'controlled-work/report/evidence'
            science=scientific_root(evidence,read_scientific(evidence/'manifest.json'))
            manifest=read_scientific(science/('scientific-bindings.json' if (science/'scientific-bindings.json').exists() else 'manifest.json'))
        recovery_check(metadata.classification=='public' and metadata.bindings==frozen['manifest_checksums']
            and manifest['frozen_envelope_sha256']==binding['snapshot_sha256'],
            'public_envelope_identity_changed')
    else:
        snapshot = small_json(root/'preparation/snapshot.json')
        recovery_check(metadata.bindings == snapshot['bindings']
            and metadata.classification == snapshot['classification']
            and metadata.input_snapshot_sha256 == binding['snapshot_sha256'],
            'construction_or_snapshot_identity_changed')
    recovery_check(metadata.bindings['input'] == small_json(root/'construction/run.json')['input_sha256'],
        'construction_identity_changed')
    return metadata


def finalization_binding(package, production):
    from genosapiens.track2_construction_dependencies import python_dependencies
    return dict(production=production, finalization=python_dependencies(package, [
        ('genosapiens.track2_report_finalization', 'checked_synthesis'),
        ('genosapiens.track2_report_finalization', '_publish_finalization')]))


def _publish_finalization(root, package, binding, metadata, *, require_review):
    from genosapiens.track2_report_assembly import presentation_ready
    recovery_check(presentation_ready(Path(root)/'synthesis/controlled-work/report'), 'overview_presentation_incomplete')
    from genosapiens.track2_terminal import assessment_outcome
    state = RunState(root)
    current = finalization_binding(package, binding)
    complete, _ = state.completed('report', current)
    if complete is None:
        state.publish('report', binding=current, result=dict(readiness=metadata.readiness,
            manifest_sha256=metadata.manifest_sha256))
    outcome = assessment_outcome(metadata.readiness, require_review=require_review)
    atomic_record(Path(root)/'current-routing.json', dict(native=metadata.kind == 'native',
        synthesis_complete=True, report_complete=True, report_binding=binding,
        finalization_exit_code=outcome.exit_code))
    atomic_record(Path(root)/'current-diagnostics.json', dict(stage='finalize', state='complete',
        assessment_status=metadata.readiness['assessment_status'], exit_code=outcome.exit_code))
    print('Report checks completed. Final review and submission remain separate.', flush=True)
    return outcome


def finalize_report(root, package, *, require_review=False):
    from genosapiens.track2_terminal import publish_outcome
    root = Path(root)
    try:
        atomic_record(root/'current-diagnostics.json', dict(stage='finalize', state='running'))
        native = (Path(package)/'discovery-resources.json').is_file()
        binding = production_binding(root, package, native=native)
        metadata = checked_synthesis(root, binding)
        return publish_outcome(_publish_finalization(root, package, binding, metadata, require_review=require_review))
    except BaseException as error:
        record_finalization_failure(root, error)
        raise


def record_finalization_failure(root, error):
    import json
    import sys
    from genosapiens.track2_discovery_outcome import safe_diagnostic
    primary = safe_diagnostic('finalize', 'report-completion', error)
    if isinstance(error, TimeoutError): primary['category'] = 'interrupted'
    try:
        route = small_json(Path(root)/'current-routing.json')
        route['report_complete'] = False
        route.pop('finalization_exit_code', None)
        atomic_record(Path(root)/'current-routing.json', route)
    except BaseException as secondary:
        print(json.dumps(dict(primary=primary,
            diagnostic_write=safe_diagnostic('finalize', 'routing', secondary))), file=sys.stderr)
    state = 'interrupted' if primary['category'] == 'interrupted' else 'failed'
    for name, value in [('current-diagnostics.json', dict(state=state, **primary)),
            ('current-outcome.json', dict(schema='genosapiens.track2-finalization-outcome/v1',
                state=state, diagnostic=primary))]:
        try:
            atomic_record(Path(root)/name, value)
        except BaseException as secondary:
            print(json.dumps(dict(primary=primary,
                diagnostic_write=safe_diagnostic('finalize', 'diagnostic-record', secondary))), file=sys.stderr)
    print(json.dumps(primary), file=sys.stderr)


def prepare_synthesis(root, package, *, native, require_review=False):
    """Check existing production before replacing any report artifacts."""
    root = Path(root)
    report=root/'synthesis/controlled-work/report/manifest.json'
    if report.is_file():
        from genosapiens.track2_scientific_reader import read_scientific
        recovery_check(read_scientific(report).get('schema')=='genosapiens.track2-challenge-package/v2',
            'historical_report_requires_separate_revision')
    try:
        binding = production_binding(root, package, native=native)
    except BaseException as error:
        record_finalization_failure(root, error)
        raise
    reason = 'synthesis_missing'
    metadata = None
    if (root/'synthesis/controlled-work/report/manifest.json').is_file() and (
            (root/'synthesis/production.json').exists() and (root/'synthesis/controlled-work/result.json').exists()):
        try:
            metadata = checked_synthesis(root, binding)
        except (ValueError, OSError, KeyError, TypeError) as error:
            from genosapiens.track2_scientific_storage import ScientificCapacityError
            import errno
            if isinstance(error, ScientificCapacityError) or isinstance(error, OSError) and error.errno != errno.ENOENT:
                record_finalization_failure(root, error)
                raise
            record_finalization_failure(root,error)
            raise
    if metadata is None:
        synthesis=root/'synthesis';request=synthesis/'report-request.json'
        if request.exists():
            recovery_check(production_matches(small_json(request)['binding'],binding),
                'incomplete_report_dependencies_changed')
        else:
            recovery_check(not synthesis.exists() or not any(synthesis.iterdir()),
                'unrecognized_incomplete_report_preserved')
            synthesis.mkdir(mode=0o700,parents=True,exist_ok=True)
            atomic_record(request,dict(schema='genosapiens.track2-report-request/v1',binding=binding))
        print('Prepare or resume the shared report workflow.', flush=True)
        atomic_record(root/'current-routing.json', dict(native=native, synthesis_complete=False,
            report_complete=False, report_binding=binding))
    else:
        from genosapiens.track2_report_assembly import presentation_ready
        if not presentation_ready(root/'synthesis/controlled-work/report'):
            atomic_record(root/'current-routing.json', dict(native=native, synthesis_complete=True,
                report_complete=False, report_binding=binding))
            return 0
        print('Reuse checked synthesis. Complete operator finalization.', flush=True)
        try:
            # Reuse this operation's semantic replay, including partial reports.
            _publish_finalization(root, package, binding, metadata, require_review=require_review)
        except BaseException as error:
            record_finalization_failure(root, error)
            raise
    return 0
