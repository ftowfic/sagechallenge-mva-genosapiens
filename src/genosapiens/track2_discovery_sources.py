"""Native public provenance without an imported study-review assertion."""
from genosapiens.track2_operator_inputs import require
from genosapiens.track2_published_genotype import canonical, sha256


def native_source_envelope(descriptor):
    if descriptor['schema']=='genosapiens.track2-discovery-resources/v2':
        producer=descriptor['producer']
        require(all(descriptor[k]==v for k,v in producer.items() if k!='schema'),
            'Tissue provenance changed an original producer field.')
        original=native_source_envelope(producer)
        return dict(original,schema='genosapiens.track2-native-source-envelope/v2',
            discovery_resources_sha256=sha256(canonical(descriptor)),producer=original,
            tissue_companion=descriptor['tissue_companion'],compatibility=descriptor['compatibility'],
            query_policy_sha256=sha256(canonical(descriptor['query_policy'])))
    require(descriptor['schema']=='genosapiens.track2-discovery-resources/v1'
        and descriptor['classification']=='public', 'Native discovery requires a public resource descriptor.')
    return dict(schema='genosapiens.track2-native-source-envelope/v1',classification='public',
        discovery_resources_sha256=sha256(canonical(descriptor)),build_id=descriptor['build_id'],
        source_manifests=descriptor['source_manifests'],matrix_producer_receipt_sha256=descriptor['matrix']['producer_receipt_sha256'],
        study_review=None,authority='public_resource_provenance')


def native_snapshot(plan, bindings, runtime, descriptor):
    require(set(runtime)=={'images','code_revision'} and runtime['images']==plan['exploratory_discovery']['runtime_images']
        and isinstance(runtime['code_revision'],str) and runtime['code_revision'], 'Native snapshot runtime bindings differ.')
    packet=plan['mechanism_assessment']['automatic_inputs']['public_packet']
    return dict(schema='genosapiens.track2-private-input-snapshot/v2',classification=plan['classification'],
        bindings=bindings,discovery_resources_sha256=plan['exploratory_discovery']['resources']['sha256'],
        policy_sha256=plan['exploratory_discovery']['policy']['sha256'],
        optional_study_sha256=sha256(canonical(packet)) if packet is not None else None,
        source_envelope=native_source_envelope(descriptor),runtime=runtime)
