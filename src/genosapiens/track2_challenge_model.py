"""Compact challenge presentation, independent of the original scientific policy."""
import json
from pathlib import Path, PurePosixPath

from genosapiens.json_schema import validate_json_schema
from genosapiens.track2_operator_inputs import require
from genosapiens.track2_scientific_storage import json_identity

SCHEMA = 'genosapiens.track2-report-presentation/v3'
ROOT = Path(__file__).resolve().parents[2]
REFERENCE_NAMES = ('context', 'objective_scores', 'historical_rankings', 'coverage',
                   'figure_index', 'disclosures', 'scientific_details')


def challenge_policy():
    return json.loads((ROOT/'config/track2/report-presentation-v3.json').read_bytes())


def check_reference(value):
    """Check canonical file references. File admission belongs to the package reader."""
    from genosapiens.track2_report_model import relative_destination
    require(isinstance(value, dict) and set(value) == {'path', 'sha256'},
            'A challenge artifact reference is incomplete.')
    import re
    require(isinstance(value['sha256'], str) and re.fullmatch('[0-9a-f]{64}', value['sha256']),
            'A challenge artifact identity is invalid.')
    path = relative_destination(value['path'])
    require('#' not in path and str(PurePosixPath(path)) == path and path != '.',
            'A challenge artifact requires a canonical file path.')
    return value


def build_challenge_model(*, input_origin, classification, source_kind, scientific_policy,
                          sources, artifacts, summary, selections, figures, supplements,
                          derived_policy=None, calibration=None, narrative_version=2):
    require(type(narrative_version) is int and narrative_version in (1, 2),
            'The challenge narrative version is unsupported.')
    policy = challenge_policy()
    model = dict(schema=SCHEMA, policy=policy, policy_sha256=json_identity(policy)['sha256'],
                 input_origin=input_origin, classification=classification, source_kind=source_kind,
                 scientific_policy=scientific_policy, derived_policy=derived_policy,
                 calibration=calibration, sources=sources, artifacts=artifacts, summary=summary,
                 selections=selections, figures=figures, supplements=supplements)
    if narrative_version != 1:
        model['narrative_version'] = narrative_version
    return check_challenge_model(model)


def check_challenge_model(model):
    schema = json.loads((ROOT/'schemas/track2/report-presentation-v3.schema.json').read_bytes())
    validate_json_schema(model, schema)
    require(type(model.get('narrative_version', 1)) is int and model.get('narrative_version', 1) in (1, 2),
            'The challenge narrative version is unsupported.')
    policy = challenge_policy()
    require(model['policy'] == policy and model['policy_sha256'] == json_identity(policy)['sha256'],
            'The challenge presentation policy changed.')
    check_reference(model['scientific_policy'])
    require(set(model['sources']) == {'input', 'plan', 'scientific_manifest'},
            'The challenge source bindings are incomplete.')
    for reference in model['sources'].values():
        check_reference(reference)
    require(set(model['artifacts']) == set(REFERENCE_NAMES) | ({'derived_rankings'} if model['derived_policy'] else set()),
            'The challenge supplementary bindings are incomplete.')
    for reference in model['artifacts'].values():
        check_reference(reference)
    if model['derived_policy'] is not None:
        from genosapiens.track2_derived_strength import policy as scoring_policy
        require(model['derived_policy'] == scoring_policy(), 'The derived scoring policy is unsupported.')
    require(model['calibration'] is None or model['derived_policy'] is not None,
            'Calibration requires a separate derived scoring policy.')
    if model['calibration'] is not None:
        check_reference(model['calibration'])
    require(len(model['selections']) <= policy['main_table_limit'], 'The challenge evidence table is too large.')
    from genosapiens.track2_challenge_selection import check_selection
    check_selection(model['selections'], has_derived_policy=model['derived_policy'] is not None)
    require(len(model['figures']) <= policy['main_figure_limit'], 'The challenge main figure limit is exceeded.')
    labels = set()
    for figure in model['figures']:
        require(figure['label'] not in labels, 'A challenge figure label is duplicated.')
        labels.add(figure['label'])
        for field in ('png', 'svg', 'data'):
            check_reference(figure[field])
    for supplement in model['supplements']:
        check_reference(supplement['artifact'])
    return model
