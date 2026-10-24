"""Source action classes remain separate from exploratory direction assumptions."""
import json

from genosapiens.track2_operator_inputs import require


def classify_action(action_type, policy):
    normalized = ' '.join((action_type or '').split()).upper()
    for category, terms in policy['actions'].items():
        if normalized in terms:
            return category
    return 'unclassified'


def lookup_gene_mechanisms(connection, gene_id):
    require(connection.execute('SELECT 1 FROM genes WHERE gene_id=?', (gene_id,)).fetchone() is not None,
            'The mechanism gene is absent from the catalog. Resolve its stable identifier before searching.')
    cursor=connection.execute('''SELECT m.*,g.group_warning FROM mechanisms m
        JOIN mechanism_genes g USING(mechanism_id) WHERE g.gene_id=? ORDER BY m.mechanism_id''',(gene_id,))
    names=[c[0] for c in cursor.description]
    return [dict(zip(names,row)) for row in cursor]


def _record(row,policy):
    source=json.loads(row['source_json'])
    classification=classify_action(row['action'],policy)
    return dict(mechanism_id=row['mechanism_id'],action=row['action'],direction_class=classification,
        description=row['description'],target_name=row['target_name'],target_type=row['target_type'],
        targets=json.loads(row['targets_json']),references=json.loads(row['references_json']),source=source,
        group_warning=row.get('group_warning'),
        assumption_warning='Positive direction is assumed and remains unconfirmed.' if classification=='assumed_positive' else None)


def nominate_for_objective(objective,rows,identities,policy):
    require(objective['desired_modulation']=='increase','This discovery policy supports increase hypotheses only.')
    candidates={}
    for row in rows:
        record=_record(row,policy)
        if record['direction_class'] not in {'reported_positive','assumed_positive'}:continue
        for identity, in identities.execute('''SELECT c.data FROM mechanism_compounds m JOIN compound_identity c
                ON c.provider='open_targets' AND c.compound_id=m.chembl_id
                WHERE m.mechanism_id=? AND c.resolution='included' ORDER BY c.compound_id''',(row['mechanism_id'],)):
            identity=json.loads(identity);key=identity['ingredient_key']
            candidate=candidates.setdefault(key,dict(ingredient_key=key,records=[],identities=[]))
            if record not in candidate['records']:candidate['records'].append(record)
            if identity not in candidate['identities']:candidate['identities'].append(identity)
    return dict(objective_id=objective['objective_id'],scenario_id=objective['scenario_id'],gene_id=objective['gene_id'],
        nominations=[candidates[k] for k in sorted(candidates)],source_records=[_record(r,policy) for r in rows])


def annotate_candidate(ingredient_key,objective,connection,policy):
    cursor=connection.execute('''SELECT DISTINCT m.* FROM mechanisms m JOIN mechanism_compounds mc USING(mechanism_id)
        JOIN compound_identity c ON c.provider='open_targets' AND c.compound_id=mc.chembl_id
        WHERE c.ingredient_key=? ORDER BY m.mechanism_id''',(ingredient_key,))
    names=[c[0] for c in cursor.description];records=[];conflicts=[]
    for values in cursor:
        row=dict(zip(names,values));targets=json.loads(row['targets_json'])
        row['group_warning']='group_target_applicability_unconfirmed' if len(targets)>1 or (
            row['target_type'] or '').upper() in {'PROTEIN FAMILY','PROTEIN COMPLEX','PROTEIN COMPLEX GROUP'} else None
        record=_record(row,policy);records.append(record)
        same=objective['gene_id'] in targets
        if record['direction_class']=='reported_negative':
            conflicts.append('same_gene_direction_conflict' if same else 'indirect_compensation_unproven')
        if same and record['direction_class'] in {'reported_positive','assumed_positive'} and 'target_product_absent' in objective['known_conflicts']:
            conflicts.append('target_product_absent')
    classes={r['direction_class'] for r in records if objective['gene_id'] in r['targets']}
    direction=('conflicting' if 'reported_negative' in classes else 'reported_positive' if 'reported_positive' in classes
        else 'assumed_positive' if 'assumed_positive' in classes else 'absent')
    return dict(records=records,direction_agreement=direction,conflicts=sorted(set(conflicts)),
        coverage='available' if records else 'mechanism_annotation_unavailable',
        target_group_warnings=sorted({r['group_warning'] for r in records if r['group_warning']}))
