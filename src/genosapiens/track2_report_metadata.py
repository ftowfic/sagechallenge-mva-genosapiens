"""Checked identities for the two current scientific report contracts."""
from dataclasses import dataclass
from pathlib import Path
import re

from genosapiens.track2_operator_inputs import require, small_json
from genosapiens.track2_published_genotype import _source_file
from genosapiens.track2_scientific_storage import file_identity
from genosapiens.track2_scientific_reader import read_scientific

STUDY_SCHEMA = 'genosapiens.track2-research-manifest/v1'
NATIVE_SCHEMA = 'genosapiens.track2-scientific-completion/v1'


@dataclass(frozen=True)
class ReportMetadata:
    kind: str
    bindings: dict
    input_snapshot_sha256: str | None
    classification: str
    readiness: dict
    manifest_sha256: str


def checked_report_metadata(output, *, expected_manifest_sha256=None, context=None, for_regeneration=False):
    """Check scientific meaning, with archive validation outside regeneration.

    Regeneration retains manifest and scientific dependency checks. It does not
    require the old presentation or complete source inventory to be unchanged.
    """
    if context is not None:
        from genosapiens.track2_checked_evidence import CheckedEvidenceContext
        require(isinstance(context,CheckedEvidenceContext), 'Report admission requires an invocation-owned context.')
        context._check_active()
    root = Path(output)
    manifest_path = _source_file(root, 'manifest.json')
    identity = file_identity(manifest_path)['sha256']
    if expected_manifest_sha256 is not None:
        require(identity == expected_manifest_sha256, 'The report differs from its independently expected identity.')
    # A study manifest can contain a large output inventory. Its existing
    # scientific contract does not acquire the compact control-record bound.
    manifest = read_scientific(manifest_path)
    require(isinstance(manifest, dict), 'The report manifest contract is invalid.')
    schema = manifest.get('schema')
    if schema == 'genosapiens.track2-challenge-package/v2':
        require(expected_manifest_sha256 is not None,
                'Challenge admission requires an independently expected report manifest identity.')
        from genosapiens.track2_report_package import check_package
        checked=check_package(root,expected_manifest_sha256=expected_manifest_sha256)
        require(file_identity(manifest_path)['sha256']==identity, 'The challenge manifest changed during validation.')
        return ReportMetadata(checked['source_kind'],dict(checked['source_bindings']),
            checked['input_snapshot_sha256'],checked['classification'],checked['readiness'],identity)
    if schema == 'genosapiens.track2-challenge-package/v1':
        require(expected_manifest_sha256 is not None,
                'Challenge admission requires an independently expected report manifest identity.')
        from genosapiens.track2_challenge_revision import check_revision
        checked = check_revision(root,expected_manifest_sha256=expected_manifest_sha256,context=context)
        require(file_identity(manifest_path)['sha256'] == identity,
                'The challenge manifest changed during validation.')
        return ReportMetadata(checked['source_kind'], dict(checked['source_bindings']),
            checked['input_snapshot_sha256'], checked['classification'], checked['readiness'], identity)
    require(schema in {STUDY_SCHEMA, NATIVE_SCHEMA},
        'The report contract is obsolete or unsupported. Preserve the source report and use a supported reader.')
    try:
        if schema == NATIVE_SCHEMA:
            require(manifest.get('stage') == 'native-report', 'The report manifest has an incompatible stage.')
            from genosapiens.track2_discovery_report import validate_native_report,validate_native_science
            pages=context.pages if context is not None else None
            if for_regeneration:
                readiness=validate_native_science(root,dict(bindings=manifest['binding'],**manifest['result']),page_context=pages)
            else:
                readiness=validate_native_report(root,page_context=pages)
            bindings = manifest['binding']
            snapshot_hash = manifest['result']['input_snapshot_sha256']
            snapshot = small_json(_source_file(root, 'input-snapshot.json'))
            classification = snapshot['classification']
            require(readiness['classification'] == classification,
                'The native report classification differs from its readiness evidence.')
            require(snapshot['schema'] == 'genosapiens.track2-private-input-snapshot/v2'
                and snapshot['bindings'] == bindings
                and file_identity(root/'input-snapshot.json')['sha256'] == snapshot_hash,
                'The native report snapshot identity differs.')
            kind = 'native'
        else:
            from genosapiens.track2_synthesis import validate_study_report,validate_study_science
            readiness=(validate_study_science(root,manifest,verify_outputs=False)
                if for_regeneration else validate_study_report(root))
            bindings = manifest['bindings']
            classification = manifest['classification']
            snapshot_hash = manifest.get('input_snapshot_sha256')
            require(('analysis_origins' in manifest) == (snapshot_hash is not None),
                'The study report snapshot identity is missing or unexpected.')
            kind = 'study'
        if for_regeneration:
            from genosapiens.track2_report_checks import check_scientific_supplements
            check_scientific_supplements(root,read_scientific(_source_file(root,'presentation.json')))
        require(isinstance(bindings, dict) and 'input' in bindings and 'plan' in bindings
            and all(v is None and k in {'input_review', 'plan_review'}
                or isinstance(v, str) and re.fullmatch('[0-9a-f]{64}', v) for k, v in bindings.items()),
            'The report input identities are invalid.')
        require(snapshot_hash is None or isinstance(snapshot_hash, str)
            and re.fullmatch('[0-9a-f]{64}', snapshot_hash), 'The report snapshot identity is invalid.')
        require(classification in {'public', 'controlled'}, 'The report classification is invalid.')
        require(readiness['assessment_status'] in {'complete', 'incomplete', 'not_performed'}
            and all(type(readiness[k]) is bool for k in
                ('input_report_ready', 'research_report_ready', 'submission_ready'))
            and not readiness['submission_ready'], 'The report readiness contract is invalid.')
        require(file_identity(manifest_path)['sha256'] == identity, 'The report manifest changed during validation.')
        return ReportMetadata(kind, dict(bindings.items()), snapshot_hash, classification, readiness, identity)
    except (KeyError, TypeError, AttributeError):
        raise ValueError('The report metadata contract is malformed.') from None
