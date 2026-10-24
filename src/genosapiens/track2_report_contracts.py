"""Versioned retained evidence and challenge package envelopes.

These checks establish structure and file identities. Producer adapters must
also replay scientific semantics before admitting evidence to a coordinator.
"""
import json
from pathlib import Path

from genosapiens.json_schema import validate_json_schema
from genosapiens.track2_challenge_model import check_reference
from genosapiens.track2_operator_inputs import require
from genosapiens.track2_scientific_storage import file_identity

EVIDENCE_SCHEMA = 'genosapiens.track2-report-evidence/v1'
PACKAGE_SCHEMA = 'genosapiens.track2-challenge-package/v2'
ROOT = Path(__file__).resolve().parents[2]
RECORDS = ('plan', 'context', 'objective_scores', 'historical_selections',
           'coverage', 'decisions', 'methods', 'reviews', 'source_versions',
           'eligibility', 'figure_index', 'scientific_details')


def _schema(value, name):
    schema = json.loads((ROOT/'schemas/track2'/(name+'.schema.json')).read_bytes())
    validate_json_schema(value, schema)


def _readiness(value, classification, kind):
    require(value.get('classification', classification) == classification,
            'Report readiness has another classification.')
    complete = value['input_report_ready'] and value['assessment_status'] == 'complete'
    require(value['research_report_ready'] == complete if kind == 'native' else
            not value['research_report_ready'] or complete,
            'Report readiness differs from the scientific outcome.')
    require(value['submission_ready'] is False,
            'Report completion cannot establish submission readiness.')


def _bindings(value):
    require({'input', 'plan', 'policy'} <= set(value),
            'Scientific source bindings are incomplete.')
    for key, digest in value.items():
        if digest is None:
            require(key in {'input_review', 'plan_review'}, 'A scientific source binding is missing.')
        else:
            check_reference(dict(path='binding', sha256=digest))


def _inventory(value, references, *, root=None):
    require(bool(value), 'A report inventory is empty.')
    for path, digest in value.items():
        check_reference(dict(path=path, sha256=digest))
        require(path != 'manifest.json', 'A report inventory cannot include its own manifest.')
    for reference in references:
        check_reference(reference)
        require(value.get(reference['path']) == reference['sha256'],
                'A report reference differs from its inventory.')
    if root is not None:
        from genosapiens.track2_challenge_revision import inventory
        require(inventory(root, omit={'manifest.json'}) == value,
                'Retained report files differ from their complete inventory.')


def check_evidence_contract(value, *, root=None):
    """Check the envelope without treating mutable hashes as scientific proof."""
    _schema(value, 'report-evidence-v1')
    _bindings(value['bindings'])
    _readiness(value['readiness'], value['classification'], value['source_kind'])
    references = list(value['records'].values())
    original = value['historical_source']
    if original is not None:
        references.extend(original[key] for key in ('manifest','report'))
    _inventory(value['outputs'], references, root=root)
    return value


def check_package_contract(value, *, root=None):
    """Check v2 structure and scan current files when a root is supplied.

    Public readers never accept an unchecked inventory override. A missing
    root requests structure checks only, as in the original contract API.
    """
    _schema(value, 'challenge-package-v2')
    _bindings(value['source_bindings'])
    _readiness(value['readiness'], value['classification'], value['source_kind'])
    references = [value[key] for key in ('evidence', 'scientific_policy',
        'presentation', 'figures', 'calibration_applicability', 'report_review')]
    references += [value[key] for key in ('calibration', 'derived_rankings', 'producer_admission')
                   if value[key] is not None]
    _inventory(value['outputs'], references, root=root)
    return value


def _complete_package_contract(value, *, root):
    """Calculate and immediately check one inventory for manifest construction.

    This private chain fixes the root and omitted path before scanning.
    It accepts no supplied inventory, callback, or persisted cache result.
    Copy, publication, and public reader boundaries still scan independently.
    """
    from genosapiens.track2_challenge_revision import inventory
    require('outputs' not in value, 'Manifest construction cannot reuse a supplied inventory.')
    value=dict(value,outputs=inventory(root,omit={'manifest.json'}))
    return check_package_contract(value)


def checked_contract_file(root, *, expected_manifest_sha256, package=False):
    """Require an independent identity even before semantic admission."""
    root = Path(root)
    check_reference(dict(path='manifest.json', sha256=expected_manifest_sha256))
    require(file_identity(root/'manifest.json')['sha256'] == expected_manifest_sha256,
            'The report contract differs from its expected identity.')
    from genosapiens.track2_scientific_reader import read_scientific
    value = read_scientific(root/'manifest.json')
    checker = check_package_contract if package else check_evidence_contract
    checker(value, root=root)
    require(file_identity(root/'manifest.json')['sha256'] == expected_manifest_sha256,
            'The report contract changed during validation.')
    return value
