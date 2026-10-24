"""Shared presentation records. Scientific producers retain decision authority."""
import json
from pathlib import Path, PurePosixPath
import re

from genosapiens.track2_operator_inputs import require
from genosapiens.track2_scientific_storage import json_identity, publish_json

SCHEMA = 'genosapiens.track2-report-presentation/v1'
POLICY_PATH = Path(__file__).resolve().parents[2]/'config/track2/report-presentation-v1.json'
METRICS = ('combined', 'lincs', 'mechanism')


def presentation_policy(*,tissue=False):
    path=POLICY_PATH.with_name('report-presentation-v2.json') if tissue else POLICY_PATH
    policy = json.loads(path.read_bytes())
    require(policy['schema'] == 'genosapiens.track2-report-presentation-policy/'+('v2' if tissue else 'v1'),
            'The report presentation policy is unsupported.')
    return policy


def relative_destination(value):
    """Require package-local destinations, including optional fragment names."""
    require(isinstance(value, str) and bool(value), 'A report destination is missing.')
    path = PurePosixPath(value.split('#', 1)[0])
    require(not path.is_absolute() and '..' not in path.parts and ':' not in value
            and '\\' not in value, 'A report destination is outside the package.')
    return value


def table_caption(title, columns, definitions, *, scope, rules, data):
    """Use the same column dictionary in rendered captions and data exports."""
    require(len(columns) == len(set(columns)) and bool(columns), 'Report columns are not unique.')
    require(set(columns) == set(definitions) and all(definitions.values()),
            'A report table header has no definition.')
    return dict(title=title, columns=list(columns), definitions=dict(definitions),
                scope=scope, rules=list(rules), data=relative_destination(data))


def build_model(*, input_origin, classification, sources, context, scores, selections,
                captions, figures, supplements, coverage, disclosure, tissue=False):
    """Retain presentation inputs without changing their scientific decisions.

    Large context and score arrays may use the scientific reader's disk-backed
    sequences. No serialized copy or per-report byte limit is introduced here.
    """
    policy = presentation_policy(tissue=tissue)
    model = dict(schema='genosapiens.track2-report-presentation/v2' if tissue else SCHEMA, policy=policy, policy_sha256=json_identity(policy)['sha256'],
                 input_origin=input_origin, classification=classification, sources=sources,
                 context=context, scores=scores, selections=selections, captions=captions,
                 figures=figures, supplements=supplements, coverage=coverage, disclosure=disclosure)
    check_model(model)
    return model


def check_model(model):
    """Check identities and navigational contracts before report composition."""
    if model.get('schema') == 'genosapiens.track2-report-presentation/v3':
        from genosapiens.track2_challenge_model import check_challenge_model
        return check_challenge_model(model)
    fields = {'schema', 'policy', 'policy_sha256', 'input_origin', 'classification', 'sources',
              'context', 'scores', 'selections', 'captions', 'figures', 'supplements',
              'coverage', 'disclosure'}
    require(set(model) == fields and model['schema'] in {SCHEMA,'genosapiens.track2-report-presentation/v2'}, 'The report presentation contract is unsupported.')
    require(model['input_origin'] in {'published_patient_genotype', 'diagnosis_guided_vcf'},
            'The report presentation origin is unsupported.')
    policy = presentation_policy(tissue=model['schema']=='genosapiens.track2-report-presentation/v2')
    require(model['policy'] == policy and model['policy_sha256'] == json_identity(policy)['sha256'],
            'The report presentation policy changed.')
    require(bool(model['sources']) and all(isinstance(value, str) and re.fullmatch('[0-9a-f]{64}', value)
            for value in model['sources'].values()), 'The report source identities are incomplete.')
    require(set(model['selections']) == set(METRICS), 'A report selection is missing.')
    for metric in METRICS:
        rows = model['selections'][metric]['rows']
        require(len(rows) <= policy['main_table_limit'] and
                len({r['canonical_active_ingredient_id'] for r in rows}) == len(rows),
                'A report table exceeds its distinct medication limit.')
        require(all(r.get('objective_id') and r.get('scenario_id') and r.get('gene_id') for r in rows),
                'A report representative has no hypothesis identity.')
    for caption in model['captions'].values():
        table_caption(**caption)
    require(sum(bool(f['main']) for f in model['figures']) <= policy['main_figure_limit'],
            'The report exceeds its main figure limit.')
    labels = set()
    for figure in model['figures']:
        require(figure['label'] not in labels and bool(figure['caption']), 'A report figure caption is missing or duplicated.')
        labels.add(figure['label'])
        for field in ('png', 'svg', 'data'):
            relative_destination(figure[field])
        require(bool(figure['provenance']), 'A report figure has no generating provenance.')
    for supplement in model['supplements']:
        relative_destination(supplement['path'])
    return model


def publish_model(path, model):
    check_model(model)
    return publish_json(path, model)
