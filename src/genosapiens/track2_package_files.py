"""One inventory for copied Track 2 runtime files and preparation fingerprints."""
from pathlib import Path

SCRIPTS = ('track2', 'scripts/run_track2.py', 'scripts/run_track2_public.py', 'scripts/portable_paths.py',
    'scripts/prepare_track2_discovery.py', 'scripts/read_track2_lincs.R',
    'scripts/render_track2_lincs.R', 'scripts/render_track2_overviews.R',
    'scripts/render_track2_challenge_overviews.R', 'scripts/revise_track2_report.py', 'scripts/report_track2_run.py',
    'scripts/diagnose_track2_report_transition.py')
REQUIRED = (*SCRIPTS,
    'src/genosapiens/track2_report_contracts.py',
    'src/genosapiens/track2_report_evidence.py',
    'src/genosapiens/track2_report_sources.py',
    'src/genosapiens/track2_report_calibration_inputs.py',
    'src/genosapiens/track2_report_calibration_stage.py',
    'src/genosapiens/track2_report_worker.py',
    'src/genosapiens/track2_report_workflow.py',
    'src/genosapiens/track2_report_stages.py',
    'src/genosapiens/track2_report_package.py',
    'src/genosapiens/track2_report_fresh.py',
    'src/genosapiens/track2_report_run_source.py',
    'src/genosapiens/track2_report_software.py',
    'src/genosapiens/track2_report_command.py',
    'schemas/track2/report-evidence-v1.schema.json',
    'schemas/track2/challenge-package-v2.schema.json',
    'src/genosapiens/track2_report_discovery_compatibility.py',
    'src/genosapiens/track2_report_producer.py',
    'src/genosapiens/track2_calibration_bank.py',
    'src/genosapiens/track2_challenge_model.py',
    'src/genosapiens/track2_challenge_overviews.py',
    'src/genosapiens/track2_challenge_paper.py',
    'src/genosapiens/track2_challenge_revision.py',
    'src/genosapiens/track2_challenge_selection.py',
    'src/genosapiens/track2_challenge_summary.py',
    'src/genosapiens/track2_checked_evidence.py',
    'src/genosapiens/track2_checked_pages.py',
    'src/genosapiens/track2_compact_ranking.py',
    'src/genosapiens/track2_derived_ranking.py',
    'src/genosapiens/track2_derived_strength.py',
    'src/genosapiens/track2_disk_memo.py',
    'src/genosapiens/track2_figure_workers.py',
    'src/genosapiens/track2_frozen_calibration.py',
    'src/genosapiens/track2_report_calibration.py',
    'src/genosapiens/track2_report_derived.py',
    'src/genosapiens/track2_report_links.py',
    'src/genosapiens/track2_report_tissue_checks.py',
    'src/genosapiens/track2_revision_finalization.py',
    'config/track2/derived-scoring-v1.json',
    'config/track2/report-presentation-v3.json',
    'schemas/track2/derived-strength.schema.json',
    'schemas/track2/frozen-background.schema.json',
    'schemas/track2/report-presentation-v3.schema.json',
 'src/genosapiens/track2_report_gene_metadata.py',
    'src/genosapiens/track2_reader_session.py', 'src/genosapiens/track2_shared_blocks.py',
    'src/genosapiens/track2_scoring_workers.py', 'src/genosapiens/track2_tissue_admission.py',
    'src/genosapiens/track2_speedup_recovery.py', 'config/track2/speedup-compatibility-v1.json',
    'src/genosapiens/track2_concept_supplement.py', 'config/track2/concept-supplement.json',
    'src/genosapiens/track2_protein_loss.py', 'src/genosapiens/track2_protein_supplement.py',
    'src/genosapiens/track2_functional_observations.py', 'config/track2/protein-representation-supplement.json',
    'src/genosapiens/track2_tissue_queries.py','src/genosapiens/track2_tissue_execution.py',
    'src/genosapiens/track2_tissue_aggregation.py','src/genosapiens/track2_tissue_ranking.py',
    'config/track2/exploratory-policy-v2.json','config/track2/report-presentation-v2.json',
    'schemas/track2/tissue-query.schema.json','schemas/track2/tissue-plan.schema.json',
    'schemas/track2/tissue-membership.schema.json','schemas/track2/tissue-condition.schema.json',
    'schemas/track2/tissue-result.schema.json','schemas/track2/tissue-score.schema.json',
    'src/genosapiens/track2_report_finalization.py', 'config/track2/report-presentation-v1.json',
    'config/track2/diagnosis-runtime.json', 'config/track2/discovery-runtime.json',
    'schemas/track2/report-presentation.schema.json', 'schemas/track2/native-figure-page.schema.json',
    'schemas/track2/scientific-completion.schema.json',
    'schemas/track2/scientific-artifact-reference.schema.json')


def runtime_files(repository):
    root = Path(repository).absolute()
    for name in REQUIRED:
        path = root/name
        if not path.is_file() or path.resolve() != path:
            raise ValueError('Required Track 2 runtime input is missing or linked: '+name)
    names = set(SCRIPTS)
    for directory in ('src/genosapiens', 'config', 'schemas'):
        for path in (root/directory).rglob('*'):
            relative = path.relative_to(root)
            if 'experiments' in relative.parts or path.suffix not in {'.py', '.json'}:
                continue
            if path.is_file():
                if path.resolve() != path:
                    raise ValueError('Package source links are not supported: '+str(relative))
                names.add(relative.as_posix())
    return tuple(root/name for name in sorted(names))
