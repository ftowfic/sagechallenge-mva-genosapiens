"""Structured source versions and the retained evaluated eligibility scope."""
import json
from pathlib import Path

from genosapiens.track2_operator_inputs import require
from genosapiens.track2_scientific_reader import read_scientific
from genosapiens.track2_scientific_storage import json_identity
from genosapiens.track2_challenge_model import check_reference


def source_versions(records, *, origin, kind):
    require(origin in {'producer_structured','legacy_markdown','retained_study'},
            'The source-version adaptation is unsupported.')
    records=list(records)
    selected={}
    for row in records:
        require(isinstance(row,dict) and isinstance(row.get('source'),str),
                'A source-version record has no source identity.')
        if row['source'] in {'fda','lincs'}:
            require(row['source'] not in selected, 'A calibration source identity is duplicated.')
            selected[row['source']]=row
    if kind=='native':
        require(set(selected)=={'fda','lincs'}, 'Native evidence requires original FDA and LINCS source versions.')
        check_reference(dict(path='fda',sha256=selected['fda']['sha256']))
        release=selected['lincs']['release']
        require(isinstance(release,str) and bool(release), 'The retained LINCS release is missing.')
        calibration=dict(state='available',fda_snapshot_sha256=selected['fda']['sha256'],source_release=release,reason=None)
    else:
        require(kind=='study', 'The scientific source kind is unsupported.')
        calibration=dict(state='unavailable',fda_snapshot_sha256=None,source_release=None,
            reason='study_source_has_no_frozen_tissue_calibration_scope')
    return dict(schema='genosapiens.track2-report-source-versions/v1',origin=origin,
        records=records,calibration=calibration)


def historical_source_versions(root, *, kind):
    """Keep Markdown parsing exclusively at the historical adapter boundary."""
    root=Path(root)
    if kind=='study':
        manifest=read_scientific(root/'manifest.json')
        return source_versions([dict(source='study',bindings=manifest['bindings'],runtime=manifest.get('runtime'))],
            origin='retained_study',kind=kind)
    text=(root/'scientific-details.md').read_text()
    require(text.count('## Source versions\n')==1, 'The historical source-version section is missing or duplicated.')
    section=text.split('## Source versions\n',1)[1].split('\n## ',1)[0]
    records=[]
    for line in section.splitlines():
        if line.startswith('`{') and line.endswith('}`'):
            records.append(json.loads(line[1:-1]))
    return source_versions(records,origin='legacy_markdown',kind=kind)


def frozen_eligibility(root, versions, *, kind):
    """Freeze every evaluated ingredient independently of report selections."""
    root=Path(root);plan=read_scientific(root/'plan-source.json')
    context=read_scientific(root/'supplements/context.json')
    nested=plan.get('exploratory_discovery',{})
    hypotheses=[dict(scenario_id=h['scenario_id'],objective_ids=[o['objective_id'] for o in h['objectives']])
        for h in context['hypotheses']]
    base=dict(schema='genosapiens.track2-report-eligibility/v1',
        source_versions_sha256=json_identity(versions)['sha256'],scope='all_retained_evaluated_ingredients',
        hypotheses=hypotheses,fda_snapshot_sha256=versions['calibration']['fda_snapshot_sha256'])
    if kind!='native' or nested.get('schema')!='genosapiens.track2-exploratory-plan/v2':
        return dict(base,state='unavailable',reason='no_retained_tissue_evaluation',ingredients=[],objectives=[])
    require(versions['calibration']['state']=='available', 'Tissue eligibility lacks its frozen source versions.')
    objectives={o['objective_id']:o for o in nested['objectives']}
    evaluated={oid:{} for oid in objectives}
    for row in read_scientific(root/'discovery/tissue-scores.json'):
        oid=row['objective_id'];ingredient=row['medication_entity_key']
        require(oid in objectives and all(row[k]==objectives[oid][k] for k in ('gene_id','scenario_id')),
                'An evaluated ingredient has another objective scope.')
        require(ingredient not in evaluated[oid], 'An evaluated ingredient is duplicated within an objective.')
        evaluated[oid][ingredient]=dict(medication_entity_key=ingredient,
            canonical_active_ingredient_id=row['canonical_active_ingredient_id'],
            nominated=row['nominated'],provisional=row['provisional'])
    rows=[]
    for oid,objective in objectives.items():
        rows.append(dict(objective_id=oid,scenario_id=objective['scenario_id'],gene_id=objective['gene_id'],
            direction=objective['desired_modulation'],
            evaluated=[evaluated[oid][key] for key in sorted(evaluated[oid])]))
    return dict(base,state='available',reason=None,
        ingredients=sorted({key for values in evaluated.values() for key in values}),objectives=rows)


def check_source_records(root, manifest, *, scientific_root=None):
    """Compare retained eligibility with the admitted full scientific export."""
    root=Path(root);science=root if scientific_root is None else Path(scientific_root)
    versions=read_scientific(root/manifest['records']['source_versions']['path'])
    require(versions['schema']=='genosapiens.track2-report-source-versions/v1', 'The source-version contract is unsupported.')
    expected=source_versions(versions['records'],origin=versions['origin'],kind=manifest['source_kind'])
    require(json_identity(expected)==json_identity(versions), 'Structured calibration source identities differ.')
    eligibility=read_scientific(root/manifest['records']['eligibility']['path'])
    require(json_identity(eligibility)==json_identity(frozen_eligibility(science,versions,kind=manifest['source_kind'])),
            'Frozen eligibility differs from the full evaluated evidence.')
    return versions,eligibility
