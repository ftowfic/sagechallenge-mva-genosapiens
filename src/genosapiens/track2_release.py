"""Small default-deny projection for the current diagnosis package/report path.

This stages structured findings for a signed, source-bound release review. It
does not publish remotely, prepare a submission, or relicense source evidence.
Historical report, pitch, PDF and repository-release interfaces are not loaded.
"""

import os
from pathlib import Path
import re
import uuid

from genosapiens.governance import (
    ArtifactClassification, ArtifactRecord, CaseInventory,
    validate_relative_classification, write_inventory,
)
from genosapiens.participant_governance import authorize_controlled_access, authorize_package_review
from genosapiens.incident_response import assert_processing_allowed
from genosapiens.track2_published_genotype import canonical, sha256
from genosapiens.track2_storage import output_storage
from genosapiens.json_schema import validate_json_schema


DOMAINS = frozenset({'disease', 'gene', 'variant', 'genotype', 'age', 'phenotype',
                     'handoff', 'checksum', 'review', 'provenance'})
SCHEMA = 'genosapiens.track2-diagnosis-public-projection/v1'
SELECTION_SCHEMA = Path(__file__).resolve().parents[2] / 'schemas/track2/diagnosis-publication-selection.schema.json'


class Track2ReleaseError(ValueError):
    """Static release rejection, without private source values."""


def _require(condition):
    if not condition:
        raise Track2ReleaseError('diagnosis release projection is invalid or unapproved')


def _validate_selection(selection):
    import json
    try:
        raw = SELECTION_SCHEMA.read_bytes()
        validate_json_schema(selection, json.loads(raw))
        return sha256(raw)
    except ValueError:
        raise Track2ReleaseError('publication selection violates the closed field schema') from None


def _selected_field(bundle, field):
    """Closed typed field set. No narrative, raw identity, checksum or source path."""
    _require(type(field) is dict and set(field) == {'domain', 'scenario_index', 'field'})
    domain, index, name = field['domain'], field['scenario_index'], field['field']
    parts = bundle.components
    if domain == 'disease' and index is None and name == 'resolved_diagnosis':
        value = parts['scope']['resolved_diagnosis']
        _require(re.fullmatch(r'(?:MONDO|OMIM|ORPHA):[0-9]+', value) is not None)
        return {'domain': domain, 'subject': 'scope', 'field': name, 'value': value}
    _require(type(index) is int and 0 <= index < len(parts['scenarios']))
    scenario = parts['scenarios'][index]
    if domain == 'gene' and name == 'gene_id':
        value = scenario['gene_id']
        _require(re.fullmatch(r'ENSG[0-9]{11}(?:\.[0-9]+)?|HGNC:[0-9]+', value) is not None)
    elif domain == 'genotype' and name in {'genotype_completeness', 'phase', 'direction'}:
        value = scenario['mechanism_inputs'][name]
        _require(value in {
            'complete', 'incomplete', 'unknown', 'not_applicable', 'confirmed_trans', 'confirmed_cis',
            'loss_of_function', 'gain_of_function', 'dominant_negative', 'mixed', 'context_dependent',
        })
    else:
        raise Track2ReleaseError('publication field is outside the current projection')
    return {'domain': domain, 'subject': f'scenario-{index + 1}', 'field': name, 'value': value}


def stage_diagnosis_projection(package, report, output, *, selection, controlled_store):
    """Stage only explicitly selected typed findings; every unselected field is denied.

    Variant, age, phenotype, handoff, checksum, review and provenance domains are
    omitted in this minimal projection. Extending its field set is an explicit
    schema change, never an unrestricted string or Markdown passthrough.
    """
    import json
    from genosapiens.track2_diagnosis_package import load_diagnosis_package
    from genosapiens.track2_report_metadata import checked_report_metadata

    package, report = controlled_store.path(package), controlled_store.path(report)
    bundle = load_diagnosis_package(package)
    _require(bundle.data['classification'] == 'controlled')
    metadata = checked_report_metadata(report)
    _require(metadata.classification == 'controlled'
             and metadata.bindings['input'] == bundle.input_sha256
             and metadata.bindings['input_review'] == bundle.review_sha256)
    schema_sha = _validate_selection(selection)
    _require(type(selection) is dict and set(selection) == {'domains', 'fields'}
             and set(selection['domains']) == DOMAINS and type(selection['fields']) is list)
    _require(all(mode in {'omit', 'selected'} for mode in selection['domains'].values()))
    fields = [_selected_field(bundle, field) for field in selection['fields']]
    _require(len({canonical(field) for field in fields}) == len(fields))
    selected_domains = {field['domain'] for field in fields}
    _require(selected_domains == {domain for domain, mode in selection['domains'].items() if mode == 'selected'})
    projection = {'schema': SCHEMA, 'public_id': 'projection-' + uuid.uuid4().hex,
                  'findings': fields}
    raw = canonical(projection)
    binding = {'schema': 'genosapiens.track2-private-projection-binding/v1',
        'classification': 'controlled', 'input_sha256': bundle.input_sha256,
        'input_review_sha256': bundle.review_sha256, 'report_sha256': metadata.manifest_sha256,
        'selection': selection, 'projection_sha256': sha256(raw), 'public_id': projection['public_id']}
    binding['publication_schema_sha256'] = schema_sha
    with output_storage(output, classification='controlled', controlled_store=controlled_store,
                        kind='track2-release-projection') as target:
        target.mkdir(parents=True, exist_ok=False, mode=0o700)
        (target / 'projection.json').write_bytes(raw)
        (target / 'binding.json').write_bytes(canonical(binding))
    return sha256(canonical(binding))


def release_diagnosis_projection(package, report, staged, output, *, controlled_store,
                                 governance, participant_id, reviewer_id, reviewed_at,
                                 incident_root, environment_id):
    """Materialize a local public projection only after exact signed review.

    The authority is the existing authenticated participant-governance manifest.
    Public files contain surrogate IDs and public-output hashes only. The private
    inventory retains the exact case/review binding. External publication and
    final package/license checks remain separate operator work.
    """
    import json
    from genosapiens.track2_diagnosis_package import load_diagnosis_package
    from genosapiens.track2_report_metadata import checked_report_metadata

    authorize_controlled_access(governance, participant_id, os.getuid())
    assert_processing_allowed(incident_root, environment_id)
    staged = controlled_store.path(staged)
    _require({p.name for p in staged.iterdir()} == {'projection.json', 'binding.json'})
    controlled_store.path(staged / 'projection.json')
    controlled_store.path(staged / 'binding.json')
    raw = (staged / 'projection.json').read_bytes()
    binding_raw = (staged / 'binding.json').read_bytes()
    binding, projection = json.loads(binding_raw), json.loads(raw)
    _require(set(binding) == {'schema', 'classification', 'input_sha256', 'input_review_sha256',
        'report_sha256', 'selection', 'projection_sha256', 'public_id', 'publication_schema_sha256'}
        and binding['schema'] == 'genosapiens.track2-private-projection-binding/v1'
        and binding['classification'] == 'controlled' and binding_raw == canonical(binding)
        and binding['projection_sha256'] == sha256(raw) and raw == canonical(projection))
    bundle = load_diagnosis_package(controlled_store.path(package),
        expected_input_sha256=binding['input_sha256'], expected_review_sha256=binding['input_review_sha256'])
    _require(bundle.data['classification'] == 'controlled')
    report = controlled_store.path(report)
    metadata = checked_report_metadata(report)
    _require(metadata.manifest_sha256 == binding['report_sha256']
        and metadata.classification == 'controlled'
        and metadata.bindings['input'] == bundle.input_sha256
        and metadata.bindings['input_review'] == bundle.review_sha256)
    expected = [_selected_field(bundle, field) for field in binding['selection']['fields']]
    selection = binding['selection']
    _require(_validate_selection(selection) == binding['publication_schema_sha256'])
    _require(set(selection) == {'domains', 'fields'} and set(selection['domains']) == DOMAINS
        and all(mode in {'omit', 'selected'} for mode in selection['domains'].values())
        and len({canonical(field) for field in expected}) == len(expected)
        and {field['domain'] for field in expected}
            == {domain for domain, mode in selection['domains'].items() if mode == 'selected'})
    _require(set(projection) == {'schema', 'public_id', 'findings'} and projection['schema'] == SCHEMA
        and projection['public_id'] == binding['public_id']
        and re.fullmatch(r'projection-[0-9a-f]{32}', projection['public_id']) is not None
        and projection['findings'] == expected)
    authorize_package_review(governance, reviewer_id, sha256(binding_raw), reviewed_at)
    target = Path(output).absolute()
    try:
        relative = target.relative_to(controlled_store.reference_root)
        validate_relative_classification(relative.as_posix(), ArtifactClassification.PUBLIC)
        if target.resolve() != target or target.exists():
            raise ValueError
    except (ValueError, OSError):
        raise Track2ReleaseError('public projection destination is invalid') from None
    public_manifest = {'schema': 'genosapiens.track2-public-projection-manifest/v1',
        'classification': 'public', 'public_id': projection['public_id'],
        'outputs': {'projection.json': sha256(raw)}, 'release_review': 'approved',
        'submission_ready': False}
    target.mkdir(parents=True, exist_ok=False, mode=0o700)
    (target / 'projection.json').write_bytes(raw)
    (target / 'manifest.json').write_bytes(canonical(public_manifest))
    inventory = controlled_store.inventory()
    records = tuple(ArtifactRecord.from_file(controlled_store.reference_root, target / name,
        ArtifactClassification.PUBLIC, 'track2-reviewed-public-projection')
        for name in ('projection.json', 'manifest.json'))
    write_inventory(CaseInventory(inventory.case_id, (*inventory.artifacts, *records)),
                    controlled_store.inventory_path)
    return public_manifest
