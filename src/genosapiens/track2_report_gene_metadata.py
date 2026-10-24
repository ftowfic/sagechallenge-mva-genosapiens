"""Compact public gene metadata, independent of scientific discovery catalogs."""
from contextlib import ExitStack
from pathlib import Path

from genosapiens.track2_operator_inputs import require
from genosapiens.track2_scientific_storage import file_identity, json_identity, publish_json

SCHEMA = 'genosapiens.track2-report-gene-metadata/v1'
FILENAME = 'report-gene-metadata.json'


def installed_gene_source(repository, reference):
    """Use the same installed source pin as ordinary Track 2 packages."""
    import json
    from genosapiens.admission_links import PinnedPublicSource
    root, reference = Path(repository), Path(reference)
    installation_inputs(root, reference)
    request = json.loads((root/'config/public-resources/track1-retention-index.json').read_bytes())
    contract = root/request['open_targets_contract_path']
    require(contract.is_relative_to(root) and contract.resolve() == contract.absolute()
            and file_identity(contract)['sha256'] == request['open_targets_contract_sha256'],
            'The installed gene source contract differs.')
    return PinnedPublicSource(contract, request['open_targets_contract_sha256'], reference,
        reference/request['open_targets_manifest_path'], request['open_targets_manifest_sha256'])


def copy_prepared_gene_metadata(repository, reference, destination, *, package_store=None):
    """Reuse the authenticated companion, without reading source partitions."""
    import hashlib
    import json
    import re
    import shutil
    root = Path(repository)
    pin = installed_gene_source(root, reference)
    store = Path(package_store) if package_store is not None else root/'.track2-packages'
    current = store/'current-package.json'
    if not current.is_file():
        return False
    record = json.loads(current.read_bytes())
    name = record['package_name']
    require(re.fullmatch(r'package-[0-9a-f]{64}', name) is not None,
            'The prepared gene metadata package location is invalid.')
    package = store/name
    manifest_path = package/'package.json'
    require(package.resolve() == package.absolute()
            and file_identity(manifest_path)['sha256'] == record['package_sha256'],
            'The prepared gene metadata package binding differs.')
    files = json.loads(manifest_path.read_bytes())['files']
    if FILENAME not in files:
        return False
    module = 'src/genosapiens/track2_report_gene_metadata.py'
    require(file_identity(package/'resources.json')['sha256'] == files['resources.json'],
            'The prepared gene source resource binding differs.')
    source = json.loads((package/'resources.json').read_bytes())['gene_source']
    if (source['contract_sha256'], source['manifest_sha256']) != (pin.contract_sha256, pin.manifest_sha256):
        return False
    require(module in files, 'The prepared gene metadata producer is missing.')
    require(file_identity(package/module)['sha256'] == files[module],
            'The prepared gene metadata producer changed.')
    # Host packaging helpers can change without invalidating the source parser.
    import ast
    def projection_code(path):
        tree = ast.parse(Path(path).read_text())
        selected = {}
        for node in tree.body:
            if isinstance(node, ast.FunctionDef) and node.name in {'metadata_record','prepare_gene_metadata'}:
                selected[node.name] = ast.dump(node, include_attributes=False)
            elif isinstance(node, ast.Assign):
                for target in node.targets:
                    if isinstance(target, ast.Name) and target.id in {'SCHEMA','FILENAME'}:
                        selected[target.id] = ast.dump(node.value, include_attributes=False)
        return selected
    if projection_code(package/module) != projection_code(Path(__file__)):
        return False
    target = Path(destination)
    require(not target.exists(), 'Choose a new gene metadata destination.')
    try:
        with (package/FILENAME).open('rb') as incoming, target.open('xb') as outgoing:
            digest = hashlib.sha256()
            while block := incoming.read(1024*1024):
                digest.update(block)
                outgoing.write(block)
        require(digest.hexdigest() == files[FILENAME], 'The prepared gene metadata file changed.')
    except BaseException:
        target.unlink(missing_ok=True)
        raise
    return True


def installation_inputs(repository, reference):
    """Bind the companion's source manifest before considering package reuse."""
    import hashlib
    import json
    root = Path(repository)
    request = json.loads((root/'config/public-resources/track1-retention-index.json').read_bytes())
    relative = request['open_targets_manifest_path']
    path = Path(reference)/relative
    require(not Path(relative).is_absolute() and '..' not in Path(relative).parts
            and path.resolve() == path.absolute(), 'The gene source manifest path is invalid.')
    with path.open('rb') as stream:
        content = stream.read(2*1024**2+1)
    digest = hashlib.sha256(content).hexdigest()
    require(len(content) <= 2*1024**2 and digest == request['open_targets_manifest_sha256'],
            'The gene source manifest binding differs.')
    return dict(source_manifest_sha256=digest, source_contract_sha256=request['open_targets_contract_sha256'],
                producer='genosapiens.track2_report_gene_metadata.prepare_gene_metadata')


def metadata_record(row, *, release, assembly, source):
    """Keep source-native gene coordinates. Never infer an allele mapping."""
    identifiers = row.get('dbXrefs') or []
    hgnc = sorted({str(x['id']) if str(x['id']).startswith('HGNC:') else 'HGNC:'+str(x['id'])
                   for x in identifiers if x.get('source', '').upper() == 'HGNC' and x.get('id')})
    raw = row.get('genomicLocation')
    locus = None
    if raw and type(raw.get('start')) is int and type(raw.get('end')) is int and 0 < raw['start'] <= raw['end']:
        locus = dict(chromosome=raw.get('chromosome'), start=raw['start'], end=raw['end'],
                     strand=raw.get('strand'), assembly=assembly, coordinate_system='one_based_closed')
    return dict(gene_id=row['id'], symbol=row.get('approvedSymbol'), name=row.get('approvedName'),
        hgnc_id=hgnc[0] if len(hgnc) == 1 else None, hgnc_ids=hgnc,
        identifiers=identifiers, locus=locus, source_locus=raw,
        release=release, assembly=assembly, coordinate_system='one_based_closed',
        source_record_id=row['id'], source=source,
        source_fields={key:row.get(key) for key in ('id','approvedSymbol','approvedName','dbXrefs','genomicLocation')})


def prepare_gene_metadata(pin, destination):
    """Project installed public targets through their existing authenticated reader.

    This public-only package step reads five columns in batches. Runtime reports
    access only the companion in their existing read-only software mount.
    """
    import duckdb
    from genosapiens import admission_links as links
    records = {}
    with ExitStack() as stack:
        files = links._Files(stack, reuse_published_outputs=True)
        release, contract, manifest, artifacts, _ = links._open_targets_metadata(files, pin)
        outputs = {o.path:o for o in manifest.outputs}
        db = stack.enter_context(duckdb.connect(':memory:'))
        db.execute('SET threads=1')
        db.execute("SET memory_limit='256MB'")
        for artifact in artifacts:
            if artifact.dataset != 'target':
                continue
            path = release/artifact.relative_path
            output = outputs[f'{contract.output_prefix}/{artifact.relative_path}']
            files.open(path, output.sha256, output.size_bytes)
            fields = {r[0] for r in db.execute('DESCRIBE SELECT * FROM read_parquet(?)', [str(path)]).fetchall()}
            require('id' in fields, 'Gene display metadata lacks stable identifiers.')
            wanted = [key for key in ('id','approvedSymbol','approvedName','dbXrefs','genomicLocation') if key in fields]
            cursor = db.execute('SELECT '+','.join('"'+key+'"' for key in wanted)+' FROM read_parquet(?)', [str(path)])
            while batch := cursor.fetchmany(1024):
                for values in batch:
                    row = dict(zip(wanted, values))
                    record = metadata_record(row, release=contract.release, assembly=contract.assembly,
                        source=dict(manifest_sha256=pin.manifest_sha256, partition=output.path, partition_sha256=output.sha256))
                    records.setdefault(row['id'], []).append(record)
        data = dict(schema=SCHEMA, classification='public', source_manifest_sha256=pin.manifest_sha256,
                    producer=dict(module_sha256=file_identity(Path(__file__))['sha256'],
                        source_contract_sha256=pin.contract_sha256, duckdb_version=duckdb.__version__),
                    records={key:sorted(value, key=lambda r:json_identity(r)['sha256']) for key,value in sorted(records.items())})
        publish_json(destination, data)
    return file_identity(destination)


def metadata_binding(package):
    path = Path(package)/FILENAME
    return file_identity(path) if path.is_file() else None


def load_gene_metadata(package, gene_ids):
    """Resolve exact stable identifiers only. Ambiguous records remain explicit."""
    path = Path(package)/FILENAME
    if not path.is_file():
        return {}
    from genosapiens.track2_scientific_reader import read_scientific
    data = read_scientific(path)
    require(data['schema'] == SCHEMA and data['classification'] == 'public',
            'The gene display metadata contract is unsupported.')
    result = {}
    for gid in set(gene_ids):
        rows = data['records'].get(gid, [])
        if not rows:
            continue
        require(all(row['gene_id'] == gid and row['source_record_id'] == gid for row in rows),
                'A gene display record has a different stable identifier.')
        if len(rows) != 1:
            result[gid] = dict(gene_id=gid, mapping_state='ambiguous', alternatives=rows,
                               source_manifest_sha256=data['source_manifest_sha256'])
        else:
            result[gid] = dict(rows[0], mapping_state='exact_stable_identifier')
    return result
