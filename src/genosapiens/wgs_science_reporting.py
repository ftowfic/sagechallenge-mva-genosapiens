"""Private research summaries with explicit scientific coverage."""
from itertools import islice
from pathlib import Path

from genosapiens.wgs_ledger import publish_bytes


def scope_text(header):
    counts = header['counts']
    mode = header['policy']['mode']
    return (f"Scientific status: {header['status']}. Scope mode: {mode}. "
        f"Evaluated contexts: {counts['evaluated']}. Deferred contexts: {counts['deferred']}. "
        f"Failed contexts: {counts['failed']}. Pending contexts: {counts['pending']}. "
        f"Checkpointed contexts: {counts['checkpointed']}. "
        "Ranks apply only to evaluated contexts. The run does not assess deferred contexts. Deferral is not a negative finding. "
        "Failed and pending contexts have no final scientific decision. "
        "These results support research. They do not provide treatment advice.")


def publish_scope_report(reader, *, directory):
    from genosapiens.wgs_run_manifest import atomic_current
    from genosapiens.wgs_scientific_artifact import ScientificArtifactReader, ScopedScientificArtifact
    if type(reader) is not ScientificArtifactReader or type(reader.artifact) is not ScopedScientificArtifact:
        raise ValueError('scoped scientific report requires an authenticated reader')
    reader._require_open()
    header = reader.scientific_scope.header
    findings = list(islice(reader.iter_scalar_records(), 20))
    value = dict(schema='genosapiens.scientific-scope-report/v1', classification='controlled', llm_policy='no-llm',
        scientific_sha256=reader.artifact.sha256, scope=header, findings=findings,
        ranking_domain='evaluated_contexts_only', deferred_are_negative=False,
        complete_science_products_allowed=header['status'] != 'incomplete')
    value['evidence_verification'] = ('performed' if reader.reader.verification_performed
        else 'omitted_under_user_authorization')
    directory = Path(directory)
    atomic_current(directory / 'science-report.json', value)
    lines = ['# Track 1 scientific scope', '', scope_text(header), '',
        ('Evidence verification: ' + value['evidence_verification'] + '.'), '',
        f"Original contexts: {header['counts']['original']}.",
        f"Selected contexts: {header['counts']['selected']}.",
        f"Contexts above the rapid target: {header['protected_overflow']}.", '',
        'The JSON report contains up to 20 ranked findings from evaluated contexts.',
        'The sealed evidence contains all evaluated findings and each context disposition.', '',
        '## Glossary', '',
        '- Context: A disease, gene, locus, and inheritance group with its complete member set.',
        '- Checkpoint: A sealed result that matching later runs can reuse.',
        '- Deferred: Outside the selected scientific scope.',
        '- Pending: Selected work without a final result.', '',
        'Structured identifiers and status values remain unchanged for traceability.', '']
    publish_bytes(directory, 'science-report.md', '\n'.join(lines).encode('utf-8'))
    return value
