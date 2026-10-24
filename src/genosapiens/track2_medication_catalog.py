"""FDA single-ingredient identity, independent of either discovery route."""
from collections import defaultdict
import csv
from functools import lru_cache
import io
import re
import zipfile

from genosapiens.track2_operator_inputs import require


def norm(text):
    return re.sub(r'[-\s]+', ' ', text.strip().upper())


def zip_rows(path, member, encoding='utf-8-sig'):
    with zipfile.ZipFile(path) as archive:
        with io.TextIOWrapper(archive.open(member), encoding=encoding) as stream:
            yield from csv.DictReader(stream, delimiter='\t')


@lru_cache(maxsize=4)
def _spectators(smiles):
    from rdkit import Chem
    from rdkit.Chem.MolStandardize import rdMolStandardize
    uncharger = rdMolStandardize.Uncharger()
    return frozenset(Chem.MolToInchiKey(uncharger.uncharge(Chem.MolFromSmiles(s))) for s in smiles)


def normalize_structure(value, policy):
    """Preserve the pilot's enumerated spectators, isotopes and stereochemistry."""
    from rdkit import Chem
    from rdkit.Chem.MolStandardize import rdMolStandardize
    if not value or value == '-666':
        return dict(full='', parent='', reason='missing_smiles')
    mol = Chem.MolFromSmiles(value)
    if mol is None:
        return dict(full='', parent='', reason='unparseable_smiles')
    key = Chem.MolToInchiKey
    uncharger = rdMolStandardize.Uncharger()
    full = key(mol)
    fragments = Chem.GetMolFrags(mol, asMols=True, sanitizeFrags=True)
    if len(fragments) == 1:
        return dict(full=full, parent=key(uncharger.uncharge(Chem.Mol(mol))), reason='single_fragment')
    spectators = _spectators(tuple(policy['spectator_smiles']))
    remaining = {key(uncharger.uncharge(Chem.Mol(f))) for f in fragments} - spectators
    if len(remaining) != 1:
        return dict(full=full, parent='', reason='unresolved_multicomponent_structure')
    return dict(full=full, parent=next(iter(remaining)), reason='explicit_counterions_removed')


def build_fda_catalog(source_paths, policy):
    """Build the entire admitted FDA pool before resolving provider compounds."""
    archive = source_paths['products']
    applications = {r['ApplNo']: r for r in zip_rows(archive, 'Applications.txt', 'cp1252')}
    approvals, statuses = defaultdict(list), defaultdict(set)
    for row in zip_rows(archive, 'Submissions.txt', 'cp1252'):
        if row['SubmissionType'] == policy['submission_type'] and row['SubmissionStatus'] == policy['submission_status']:
            approvals[row['ApplNo']].append(row['SubmissionStatusDate'])
    for row in zip_rows(archive, 'MarketingStatus.txt', 'cp1252'):
        statuses[(row['ApplNo'], row['ProductNo'])].add(row['MarketingStatusID'])
    products, exclusions = [], []
    for row in zip_rows(archive, 'Products.txt', 'cp1252'):
        application = applications.get(row['ApplNo'], {})
        status = statuses[(row['ApplNo'], row['ProductNo'])]
        if (application.get('ApplType') not in policy['application_types'] or not approvals[row['ApplNo']]
                or not status.intersection(policy['marketing_status_ids'])):
            continue
        name = norm(row['ActiveIngredient'])
        if not name or policy['mixture_separator'] in name:
            exclusions.append(dict(source='fda', identifier=row['ApplNo']+':'+row['ProductNo'], reason='not_single_ingredient'))
            continue
        products.append(dict(application=application['ApplType']+row['ApplNo'],
            product_number=row['ProductNo'], ingredient_name=name,
            marketing_status=sorted(status), approval_dates=sorted(approvals[row['ApplNo']]), source=row))
    names = {p['ingredient_name'] for p in products}
    records = {r['UNII']: r for r in zip_rows(source_paths['records'], 'UNII_Records_4Aug2026.txt')}
    preferred, aliases = defaultdict(set), defaultdict(set)
    for uid, record in records.items():
        if norm(record['DISPLAY_NAME']) in names:
            preferred[norm(record['DISPLAY_NAME'])].add(uid)
    for row in zip_rows(source_paths['names'], 'UNII_Names_4Aug2026.txt'):
        if row['TYPE'] != 'bn' and norm(row['NAME']) in names:
            aliases[norm(row['NAME'])].add(row['UNII'])
    resolved = {}
    for name in sorted(names):
        matches = preferred.get(name) or aliases.get(name, set())
        if len(matches) != 1 or not matches <= records.keys():
            exclusions.append(dict(source='fda', identifier=name,
                reason='ambiguous_ingredient_name' if matches else 'unresolved_ingredient_name', matches=sorted(matches)))
        else:
            resolved[name] = next(iter(matches))
    medications, full_index, parent_index, uid_keys = {}, defaultdict(set), defaultdict(set), {}
    for uid in sorted(set(resolved.values())):
        record = records[uid]
        structure = normalize_structure(record['SMILES'], policy)
        given = record['INCHIKEY'].strip()
        reason = None
        if given and structure['full'] and given != structure['full']:
            reason = 'source_structure_disagreement'
        elif structure['reason'] == 'unresolved_multicomponent_structure':
            reason = structure['reason']
        key = structure['parent'] or given
        if reason or not key:
            exclusions.append(dict(source='unii', identifier=uid, reason=reason or 'missing_structure'))
            continue
        uid_keys[uid] = key
        medication = medications.setdefault(key, dict(ingredient_key=key, unii=[], names=[], full_keys=[], parent_keys=[]))
        medication['unii'].append(uid)
        medication['names'].extend(n for n, u in resolved.items() if u == uid)
        if given:
            full_index[given].add(key)
            medication['full_keys'].append(given)
        if structure['parent']:
            parent_index[structure['parent']].add(key)
            medication['parent_keys'].append(structure['parent'])
    admitted_products = []
    for product in products:
        uid = resolved.get(product['ingredient_name'])
        if uid in uid_keys:
            admitted_products.append(dict(product, ingredient_key=uid_keys[uid], unii=uid))
    for medication in medications.values():
        for field in ['unii', 'names', 'full_keys', 'parent_keys']:
            medication[field] = sorted(set(medication[field]))
    return dict(schema='genosapiens.track2-fda-catalog/v1', snapshot=policy['snapshot'],
        medications=medications, products=admitted_products, exclusions=exclusions,
        full_index={k: sorted(v) for k, v in full_index.items()},
        parent_index={k: sorted(v) for k, v in parent_index.items()})


def resolve_compound(provider, record, medication_catalog, policy):
    require(provider in {'open_targets', 'lincs'}, 'Unknown medication identity provider.')
    if provider == 'open_targets':
        identifier, smiles, given = record['id'], record.get('canonicalSmiles'), record.get('inchiKey')
        name = record.get('name') or identifier
    else:
        identifier, smiles, given = record['pert_id'], record.get('canonical_smiles'), record.get('inchi_key')
        name = record.get('pert_iname') or identifier
    given = '' if not given or given == '-666' else given.strip()
    structure = normalize_structure(smiles, policy)
    conflict = bool(given and structure['full'] and given != structure['full'])
    exact = set(medication_catalog['full_index'].get(given, [])) if given else set()
    parents = set(medication_catalog['parent_index'].get(structure['parent'], [])) if structure['parent'] else set()
    matches = exact or parents
    reason = 'no_verified_current_FDA_match'
    if conflict:
        reason = 'source_structure_disagreement'
    elif structure['reason'] == 'unresolved_multicomponent_structure':
        reason = structure['reason']
    elif len(matches) > 1 or exact and parents and exact != parents:
        reason = 'ambiguous_structural_identity'
    elif len(matches) == 1:
        reason = 'exact_full_inchikey' if exact else 'stereo_preserved_parent'
    included = reason in {'exact_full_inchikey', 'stereo_preserved_parent'}
    ingredient = next(iter(matches)) if included else None
    return dict(provider=provider, compound_id=identifier, name=name, ingredient_key=ingredient,
        resolution='included' if included else 'unresolved', reason=reason, structure=structure,
        source_inchikey=given, source=record,
        medication=medication_catalog['medications'][ingredient] if included else None,
        products=[p for p in medication_catalog['products'] if p['ingredient_key'] == ingredient] if included else [])


def identity_assertions(resolution):
    """Adapt admitted structural identity to the existing canonical union owner."""
    from genosapiens.track2_drug_identity import DrugIdentityAssertion
    require(resolution['resolution'] == 'included', 'An unresolved compound cannot enter the candidate union.')
    ingredient = resolution['ingredient_key']
    exact = resolution['source_inchikey'] or resolution['structure']['full']
    source = 'chembl' if resolution['provider'] == 'open_targets' else 'lincs'
    identifier = resolution['compound_id'] if source == 'chembl' else 'GSE92742:'+resolution['compound_id']
    assertions = [DrugIdentityAssertion(source, identifier, resolution['name'],
        'salt' if exact != ingredient else 'active_ingredient', exact_substance_key=exact,
        active_moiety_keys=(ingredient,) if exact != ingredient else ())]
    for product in resolution['products']:
        key = product['application']+':'+product['product_number']
        assertions.append(DrugIdentityAssertion('drugs-at-fda', key, product['ingredient_name'],
            'product', ingredient_keys=(ingredient,), regulator_product_id=key))
    return tuple(assertions)
