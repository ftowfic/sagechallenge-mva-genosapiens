"""Small, authenticated reference-role metadata for scoped construction."""

from genosapiens.reference_concordance import load_reference_identity, ReferenceRoleBinding
from genosapiens.submission_rendering import load_submission_contig_renderer
from genosapiens.track2_published_genotype import sha256


MAPPING_SCHEMA = 'genosapiens.track2-scope-reference/v1'


def scope_reference_mapping(reference):
    """Compose declared source -> canonical -> public coordinate names.

    Only alias and index metadata are read. Sequence equivalence remains a
    separate reference-concordance responsibility. Declared extra sequences
    retain their identity even when the public role has no corresponding name.
    """
    if sha256(reference.contract_path.read_bytes()) != reference.contract_sha256:
        raise ValueError('Reference identity changed.')
    contract = load_reference_identity(reference.contract_path)
    renderers = {}
    for name in ('challenge_alignment', 'public_annotation_benchmark'):
        role = contract.reference_role(name)
        renderers[name] = load_submission_contig_renderer(
            identity_path=reference.contract_path, reference_manifest_path=reference.manifest_path,
            reference_root=reference.root, reference_binding=ReferenceRoleBinding(name, role.identity_id))
    source, target = renderers.values()
    if any(r.reference_manifest_sha256 != reference.manifest_sha256 for r in renderers.values()):
        raise ValueError('Reference manifest changed.')
    return {'schema': MAPPING_SCHEMA, 'reference_sha256': reference.manifest_sha256,
        'contract_sha256': reference.contract_sha256, 'aliases_sha256': source.contig_aliases_sha256,
        'source_role': source.reference_binding.role,
        'source_identity_id': source.reference_binding.identity_id,
        'scope_role': target.reference_binding.role,
        'scope_identity_id': target.reference_binding.identity_id,
        'source_lengths': {name: length for name, (_, length) in source.source_names.items()},
        'scope_lengths': {name: length for name, (_, length) in target.source_names.items()},
        'contig_map': {name: target.reference_aliases.get(canonical)
            for name, (canonical, _) in source.source_names.items()}}


def vep_index_binding(contract_path, reference_root):
    """Resolve only the VEP contract's FASTA index, never the source index."""
    from pathlib import Path
    from genosapiens.vep_bundle import load_vep_bundle_contract
    contract_path = Path(contract_path)
    contract = load_vep_bundle_contract(contract_path)
    return {'role': 'vep_annotation', 'contract_sha256': sha256(contract_path.read_bytes()),
        'path': str(Path(reference_root)/contract.fasta.fai_path),
        'sha256': contract.fasta.fai_sha256, 'size_bytes': contract.fasta.fai_size_bytes}


def read_vep_index(config):
    from pathlib import Path
    from genosapiens.bounded_json import read_bounded_bytes
    expected = vep_index_binding(config['contract_path'], config['reference_root'])
    if config['fai'] != expected:
        raise ValueError('The annotation index role or binding changed.')
    raw = read_bounded_bytes(Path(expected['path']))
    if len(raw) != expected['size_bytes'] or sha256(raw) != expected['sha256']:
        raise ValueError('The annotation index content changed.')
    return raw
