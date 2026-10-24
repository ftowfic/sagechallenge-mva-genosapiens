"""Bounded local Level 5 reader and FDA discovery filter (not G1 approval)."""
import csv
import json
import math
from pathlib import Path
import re
import subprocess
import tempfile

from genosapiens.json_schema import validate_json_schema
from genosapiens.track2_published_genotype import canonical, sha256, _source_file

ROOT = Path(__file__).resolve().parents[2]
ROUTE = 'lincs_cell_specific_expression'


def read_csv(path):
    with Path(path).open(newline='') as stream:
        return list(csv.DictReader(stream))


def number(value):
    try:
        result = float(value)
        return result if math.isfinite(result) and result != -666 else None
    except (ValueError, TypeError):
        return None


def unique(rows, key):
    result = {r[key]: r for r in rows}
    if len(result) != len(rows) or '' in result:
        raise ValueError('duplicate or missing ' + key)
    return result


def load_source_bundle(root, manifest_path):
    manifest = json.loads(Path(manifest_path).read_bytes())
    validate_json_schema(manifest, json.loads((ROOT / 'schemas/track2/lincs-source.schema.json').read_bytes()))
    paths = {}
    for key, binding in manifest['inputs'].items():
        path = _source_file(Path(root), binding['path'])
        # Only small retained extracts/receipts are admitted, never the GEO matrix.
        if path.stat().st_size > manifest['maximum_input_bytes']:
            raise ValueError('LINCS input exceeds declared small-extract bound')
        if sha256(path.read_bytes()) != binding['sha256']:
            raise ValueError('LINCS input checksum mismatch: ' + key)
        paths[key] = path
    required = {'allowlist', 'products', 'identity_audit', 'landmarks'}
    for cell in manifest['cells']:
        required.update(cell[k] for k in ('matrix', 'drug_metadata', 'genetic_values', 'genetic_metadata'))
    required.update(manifest['provenance_inputs'])
    if required - paths.keys():
        raise ValueError('missing local LINCS source binding')
    if len({c['cell_id'] for c in manifest['cells']}) != len(manifest['cells']):
        raise ValueError('duplicate selected cell')
    return manifest, paths


def read_matrix(path):
    """Exercise the actual pinned R/HDF5 reader, preserving row/column order."""
    with tempfile.TemporaryDirectory() as temp:
        output = Path(temp) / 'values.csv'
        subprocess.run(['Rscript', str(ROOT / 'scripts/read_track2_lincs.R'), str(path), str(output)],
                       check=True, capture_output=True, text=True)
        rows = read_csv(output)
    unique(rows, 'gene_id')
    return rows


def resolve_allowlist(allowlist, products, audit):
    allowed, audited = unique(allowlist, 'pert_id'), unique(audit, 'pert_id')
    resolved, excluded = {}, {}
    for pert_id, row in allowed.items():
        if audited.get(pert_id) != row:
            raise ValueError('allowlist disagrees with retained identity audit')
        key = row['medication_entity_key']
        method = row['identity_match']
        exact = row['calculated_inchikey']
        source = row['source_inchikey']
        if (row['status'] != 'included' or not re.fullmatch(r'[A-Z]{14}-[A-Z]{10}-[A-Z]', key)
                or method not in {'exact_full_inchikey', 'stereo_preserved_parent'}
                or (source not in {'', '-666'} and source != exact)
                or key != (row['normalized_parent_inchikey'] or exact)):
            excluded[pert_id] = 'identity_unresolved_or_conflicting'
            continue
        linked = [p for p in products if p['pert_id'] == pert_id]
        valid = [p for p in linked if re.fullmatch(r'(NDA|ANDA|BLA)\d+', p['application'])
                 and p['marketing_status'] in {'1', '2'} and p['original_approval_date']
                 and p['form_route'] and p['product_number'] and p['approval_url']
                 and p['fda_ingredient']
                 and p['unii'] in row['fda_unii'].split('|')]
        if not valid or len(valid) != len(linked):
            excluded[pert_id] = 'approval_or_active_single_ingredient_product_unresolved'
            continue
        resolved[pert_id] = {**row, 'products': valid}
    for pert_id, row in audited.items():
        if pert_id not in resolved and pert_id not in excluded:
            excluded[pert_id] = row['status']
    return resolved, excluded


def load_cell(cell, paths, allowed):
    landmarks = unique(read_csv(paths['landmarks']), 'pr_gene_id')
    genetic = read_csv(paths[cell['genetic_values']])
    metadata = read_csv(paths[cell['drug_metadata']])
    gmeta = read_csv(paths[cell['genetic_metadata']])
    matrix = read_matrix(paths[cell['matrix']])
    genes = [r['gene_id'] for r in genetic]
    unique(genetic, 'gene_id'); unique(metadata, 'sig_id'); unique(metadata, 'db_column_id')
    unique(gmeta, 'sig_id')
    if (not genes or genes != [r['gene_id'] for r in matrix]
            or set(genes) != set(landmarks)
            or any(r['cell_id'] != cell['cell_id'] for r in metadata + gmeta)):
        raise ValueError('gene/signature/cell alignment mismatch')
    columns = list(matrix[0])[1:]
    if columns != [m['db_column_id'] for m in metadata] or len(columns) > cell['maximum_drug_signatures']:
        raise ValueError('signature columns disagree with metadata or bound')
    if any(r.get('pr_is_lm', '1') not in {'1', '1.0'} for r in landmarks.values()):
        raise ValueError('inferred gene in measured-landmark universe')
    for r in genetic:
        if r['gene_symbol'] != landmarks[r['gene_id']]['pr_gene_symbol']:
            raise ValueError('gene symbol disagrees with landmark identity')
    admitted, rejected = [], []
    for meta in metadata:
        reason = None
        identity = allowed.get(meta['pert_id'])
        if identity is None:
            reason = 'identity_or_FDA_filter_unresolved'
        elif meta.get('medication_entity_key') != identity['medication_entity_key']:
            raise ValueError('drug metadata conflicts with resolved ingredient')
        elif meta['pert_type'] != 'trt_cp' or meta['pert_time_unit'] != 'h' or number(meta['pert_time']) not in {6, 24}:
            reason = 'outside_drug_context'
        elif any(number(meta.get(k)) is None or number(meta[k]) < minimum
                 for k, minimum in [('distil_nsample', 2), ('distil_cc_q75', .2), ('tas', .2)]):
            reason = 'failed_or_missing_drug_QC'
        if reason:
            rejected.append({**meta, 'reason': reason})
        else:
            admitted.append({**meta, 'identity': identity,
                             'values': [number(r[meta['db_column_id']]) for r in matrix]})
    return {'cell_id': cell['cell_id'], 'genes': genes, 'genetic': genetic,
            'genetic_metadata': gmeta, 'conditions': admitted, 'exclusions': rejected}


def freeze_queries(bundle_root, plan_path, source_path, output):
    from genosapiens.track2_synthesis import load_plan
    from genosapiens.track2_codec import decode
    from genosapiens.track2_target_hypothesis import TargetHypothesis
    from genosapiens.track2_lincs_scoring import refine_query
    bundle, plan, _, _ = load_plan(bundle_root, plan_path)
    manifest, paths = load_source_bundle(bundle_root, source_path)
    targets = [decode(TargetHypothesis, r['target']) for r in plan['tables']['target_records']]
    if len(targets) != 1:
        raise ValueError('bounded LINCS query requires exactly one explicitly selected target')
    allowed, _ = resolve_allowlist(read_csv(paths['allowlist']), read_csv(paths['products']), read_csv(paths['identity_audit']))
    queries = {}
    for config in manifest['cells']:
        cell = load_cell(config, paths, allowed)
        queries[config['cell_id']] = refine_query(cell, config, scenario_id=bundle.admitted_origin().scenario_id,
            target=targets[0], source_sha256=sha256(Path(source_path).read_bytes()))
    target = Path(output); target.mkdir(parents=True, exist_ok=False)
    for name, query in queries.items():
        (target / (name + '.json')).write_bytes(canonical(query))
    (target / 'manifest.json').write_bytes(canonical({
        'source_sha256':sha256(Path(source_path).read_bytes()),
        'queries':{name:sha256(canonical(query)) for name,query in queries.items()}}))
    return {name:{k:q[k] for k in ['state','gene_count','query_sha256','reason']} for name,q in queries.items()}
