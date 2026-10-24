"""Immutable metadata companions preserve the original catalog and matrix receipts."""
import json
from pathlib import Path
import tempfile

from genosapiens.track2_discovery_resources import (file_digest, file_state, load_discovery_resources,
    public_path, query_policy_compatibility)
from genosapiens.track2_operator_inputs import require
from genosapiens.track2_run_state import atomic_record
from genosapiens.track2_tissue_membership import catalog_membership
from genosapiens.track2_tissue_records import identity, validate_tissue_record

SCHEMA = 'genosapiens.track2-discovery-resources/v2'


def prepare_tissue_companion(descriptor, reference_root, source_contract, producer_policy, query_policy):
    root = Path(reference_root).absolute()
    compatibility = query_policy_compatibility(descriptor, root, producer_policy, query_policy)
    database = public_path(root, descriptor['catalog']['path'])
    receipt = json.loads((database.parent/'completion.json').read_bytes())
    require(receipt['bindings']['source_contract_sha256']==identity(source_contract),
        'The tissue preparation source contract differs from the catalog producer.')
    source = source_contract['lincs']
    manifest = next(r for r in descriptor['source_manifests'] if r['source']=='lincs')
    require(manifest['release']==source['release'] and
        file_digest(public_path(root, manifest['path']))==manifest['sha256'],
        'The retained cell metadata provenance changed.')
    metadata_sha256 = source['metadata_sha256'][source['metadata_files']['cell_info']]
    connection = load_discovery_resources(descriptor, root)
    try:
        cells = catalog_membership(connection, source_release=source['release'], metadata_sha256=metadata_sha256)
    finally:
        connection.close()
    membership = dict(schema='genosapiens.track2-tissue-companion/v1', classification='public',
        producer_descriptor_sha256=identity(descriptor), metadata_sha256=metadata_sha256, cells=cells)
    binding = dict(compatibility=compatibility, membership_sha256=identity(membership))
    destination = public_path(root, 'derived/track2/discovery/tissue-v1/'+identity(binding))
    if destination.exists():
        saved = json.loads((destination/'manifest.json').read_bytes())
        load_tissue_resources(saved, root).close()
        require(saved['compatibility']==compatibility and saved['tissue_companion']['sha256']==identity(membership),
            'The immutable tissue companion has a different binding.')
        return saved
    destination.parent.mkdir(parents=True, exist_ok=True)
    staging = Path(tempfile.mkdtemp(prefix='.tissue-', dir=destination.parent))
    atomic_record(staging/'membership.json', membership)
    path = staging/'membership.json'
    path.chmod(0o444)
    require(file_digest(path)==identity(membership), 'Tissue metadata serialization differs from its identity.')
    completion = dict(schema='genosapiens.track2-tissue-preparation-completion/v1', status='complete',
        bindings=binding, artifacts={'membership.json':dict(sha256=file_digest(path), file_state=file_state(path))})
    atomic_record(staging/'completion.json', completion)
    companion = dict(path=(destination/'membership.json').relative_to(root).as_posix(),
        sha256=file_digest(path), completion_sha256=file_digest(staging/'completion.json'),
        metadata_sha256=metadata_sha256, source_release=source['release'])
    result = dict(descriptor, schema=SCHEMA, producer=descriptor, producer_policy=producer_policy,
        query_policy=query_policy, compatibility=compatibility, tissue_companion=companion)
    atomic_record(staging/'manifest.json', result)
    for name in ('completion.json', 'manifest.json'):
        (staging/name).chmod(0o444)
    staging.rename(destination)
    return result


def load_tissue_membership(descriptor, reference_root):
    entry = descriptor['tissue_companion']
    path = public_path(reference_root, entry['path'])
    receipt_path = path.parent/'completion.json'
    require(file_digest(receipt_path)==entry['completion_sha256'], 'The tissue completion identity changed.')
    receipt = json.loads(receipt_path.read_bytes())
    require(receipt['schema']=='genosapiens.track2-tissue-preparation-completion/v1' and receipt['status']=='complete'
        and receipt['bindings']==dict(compatibility=descriptor['compatibility'], membership_sha256=entry['sha256']),
        'The tissue companion is incomplete or has different bindings.')
    require(receipt['artifacts']['membership.json']==dict(sha256=entry['sha256'], file_state=file_state(path))
        and file_digest(path)==entry['sha256'], 'The completed tissue metadata changed.')
    membership = json.loads(path.read_bytes())
    require(membership['schema']=='genosapiens.track2-tissue-companion/v1'
        and membership['classification']=='public'
        and membership['producer_descriptor_sha256']==identity(descriptor['producer']),
        'The tissue companion binds a different producer.')
    cells = membership['cells']
    require(len({r['cell_id'] for r in cells})==len(cells), 'The tissue companion has duplicate cells.')
    for row in cells:
        validate_tissue_record(row, 'membership')
        require(row['metadata_sha256']==membership['metadata_sha256'], 'The cell metadata identity differs.')
    return membership


def planning_membership(descriptor, connection):
    """Reproduce the small companion from an already admitted immutable catalog."""
    entry = descriptor['tissue_companion']
    cells = catalog_membership(connection, source_release=entry['source_release'], metadata_sha256=entry['metadata_sha256'])
    record = dict(schema='genosapiens.track2-tissue-companion/v1', classification='public',
        producer_descriptor_sha256=identity(descriptor['producer']), metadata_sha256=entry['metadata_sha256'], cells=cells)
    require(identity(record)==entry['sha256'], 'Retained catalog annotations differ from the tissue companion.')
    return cells


def load_tissue_resources(descriptor, reference_root, *, matrix_required=True):
    producer = descriptor['producer']
    extras = {'producer', 'producer_policy', 'query_policy', 'compatibility', 'tissue_companion'}
    require(descriptor['schema']==SCHEMA and producer['schema']=='genosapiens.track2-discovery-resources/v1'
        and set(descriptor)==set(producer)|extras and all(descriptor[k]==v for k,v in producer.items() if k!='schema'),
        'The tissue discovery descriptor changed its producer fields.')
    expected = query_policy_compatibility(producer, reference_root, descriptor['producer_policy'],
        descriptor['query_policy'], matrix_required=matrix_required)
    require(descriptor['compatibility']==expected, 'The producer-to-query compatibility record changed.')
    load_tissue_membership(descriptor, reference_root)
    return load_discovery_resources(producer, reference_root, matrix_required=matrix_required)
