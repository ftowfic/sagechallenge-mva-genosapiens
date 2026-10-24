"""Complete a retained primary's missing case seal in a fresh producer.

This creates a scientific candidate, not qualified output. The ordinary full
ScientificArtifactReader and publication checks remain mandatory. A caller
must independently bind the failed producer, inventory, primary digest, counts,
parents and reference; this module supplies no general resume authorization.
"""
from contextlib import ExitStack
from pathlib import Path

from genosapiens.genotype_hypotheses import hypothesis_id_for_members
from genosapiens.wgs_case_context import VerifiedScientificCaseContext
from genosapiens.wgs_case_evidence import _replay_selected_groups, MAX_CONTEXT_DIRECTORY_ROW_BYTES
from genosapiens.wgs_case_identities import MAX_IDENTITY_ROW_BYTES
from genosapiens.wgs_case_resources import VerifiedScientificCaseResources, bound_reference
from genosapiens.wgs_context_order import MAX_ORDER_KEY_BYTES
from genosapiens.wgs_evidence_store import EvidenceReader, EvidenceWriter, EvidenceStoreError, _require
from genosapiens.wgs_locator_tape import VerifiedLocatorTape
from genosapiens.wgs_reference_evidence import VerifiedSubmissionReferenceMetadata, METADATA_KIND
from genosapiens.wgs_scientific_artifact import ScientificArtifact, DESCRIPTOR_NAME, _write_descriptor

_COUNTS = ('context_count', 'positive_context_count', 'member_count', 'pair_count', 'hypothesis_count')


def _ref(reader, value, kind):
    _require(type(value) is dict and value.get('kind', kind) == kind)
    return bound_reference(reader, {**value, 'kind': kind}, kind)


def _selected_entries(reader, primary, stack):
    order = reader.reconstruct(_ref(reader, primary['order'], 'complete-case-order-v1'))
    selection = reader.reconstruct(_ref(reader, primary['selection'], 'complete-submission-selection-v1'))
    _require(order['schema'] == 'genosapiens.wgs-complete-case-order' and order['version'] == 1
        and selection['schema'] == 'genosapiens.wgs-submission-selection' and selection['version'] == 1
        and order['context_count'] == primary['context_count']
        and order['hypothesis_count'] == selection['hypothesis_count'] == primary['hypothesis_count']
        and selection['order'] == primary['order'] and selection['inputs'] == primary['selection_inputs'])
    ranks = selection['selected_ranks']
    _require(type(ranks) is list and len(ranks) == selection['selected_count'] <= 10
        and all(type(rank) is int and 1 <= rank <= primary['hypothesis_count'] for rank in ranks)
        and ranks == sorted(set(ranks)))
    tape = stack.enter_context(VerifiedLocatorTape(reader,
        _ref(reader, order['locators'], 'locator-tape-v1'), width=3))
    _require(tape.count == primary['hypothesis_count'])
    entries, start, families = [], 0, 0
    for family in reader.iter_array(_ref(reader, order['families'], 'case-order-families-v1'),
                                   max_item_bytes=MAX_ORDER_KEY_BYTES + 1024):
        families += 1
        _require(family['ordinal'] == families and family['start'] == start
            and type(family['count']) is int and family['count'] > 0)
        for rank in ranks:
            if start <= rank - 1 < start + family['count']:
                locator, = tape.rows(start=rank - 1, count=1)
                _require(1 <= locator[0] <= primary['context_count'])
                entries.append((rank, locator, family['total']))
        start += family['count']
    _require(families == order['family_count'] and start == primary['hypothesis_count']
        and [row[0] for row in entries] == ranks)
    return entries, ranks


def _selected_directories(reader, primary, ordinals):
    contexts, identities = {}, {}
    sums = dict.fromkeys(_COUNTS, 0)
    for row in reader.iter_array(_ref(reader, primary['contexts'], 'scientific-case-contexts-v1'),
                                max_item_bytes=MAX_CONTEXT_DIRECTORY_ROW_BYTES + 1):
        sums['context_count'] += 1
        _require(set(row) == {'ordinal', 'context_id', 'context', 'member_count', 'pair_count', 'hypothesis_count'}
            and row['ordinal'] == sums['context_count']
            and all(type(row[key]) is int and row[key] >= 0
                    for key in ('member_count', 'pair_count', 'hypothesis_count')))
        sums['positive_context_count'] += int(row['member_count'] > 0)
        for key in ('member_count', 'pair_count', 'hypothesis_count'):
            sums[key] += row[key]
        if row['ordinal'] in ordinals:
            _ref(reader, row['context'], 'scientific-case-context-v1')
            contexts[row['ordinal']] = row
    _require(all(sums[key] == primary[key] for key in _COUNTS) and set(contexts) == ordinals,
             'retained primary context coverage changed')
    header = reader.reconstruct(_ref(reader, primary['identities'], 'case-scientific-identities-v1'))
    _require(header['schema'] == 'genosapiens.wgs-case-scientific-identities' and header['version'] == 1
        and all(header[key] == primary[key] for key in ('context_count', 'member_count', 'hypothesis_count')))
    count = 0
    for row in reader.iter_array(_ref(reader, header['contexts'], 'scientific-identity-contexts-v1'),
                                max_item_bytes=MAX_IDENTITY_ROW_BYTES + 1):
        count += 1
        _require(row['ordinal'] == count)
        if count in ordinals:
            _require(all(row[key] == contexts[count][key]
                for key in ('context_id', 'context', 'member_count', 'hypothesis_count')))
            identities[count] = row
    _require(count == primary['context_count'] and set(identities) == ordinals)
    return contexts, identities


def _identity(reader, directories, locator):
    ordinal, left, right = locator
    row = directories[ordinal]
    _require(type(left) is int and type(right) is int and left > 0 and right >= 0)
    if right == 0:
        _require(left <= row['singleton_count'])
        value = reader.array_item(_ref(reader, row['singletons'], 'scientific-singleton-identities-v1'),
                                  left - 1, max_item_bytes=MAX_IDENTITY_ROW_BYTES + 1)
        _require(type(value) is list and len(value) == 3 and value[0] == left)
        return value[1]
    _require(row['pair_context_id'] is not None and left < right <= row['member_count'])
    reference = _ref(reader, row['members'], 'scientific-member-identities-v1')
    members = []
    for member in (left, right):
        value = reader.array_item(reference, member - 1, max_item_bytes=MAX_IDENTITY_ROW_BYTES + 1)
        _require(type(value) is list and len(value) == 3 and value[0] == member)
        members.append(value[1])
    return hypothesis_id_for_members(row['pair_context_id'], tuple(members), row['sample_name'])


def repair_primary_case(evidence, *, directory, output_directory, scratch_directory,
        expected_primary_sha256, expected_context_count, expected_hypothesis_count,
        expected_parent_bindings, expected_reference_binding, eligibility_root=None,
        profile_root=None, progress=None, authenticated_resources=None):
    """Produce missing selected/case roots; never overwrite the failed input."""
    source, output, scratch = map(Path, (directory, output_directory, scratch_directory))
    _require(source.is_absolute() and output.is_absolute() and scratch.is_absolute()
        and not (source / DESCRIPTOR_NAME).exists()
        and not output.exists() and not output.is_symlink()
        and not any(path.is_symlink() for path in (source, *source.parents, *output.parents, scratch, *scratch.parents))
        and source != output and not output.is_relative_to(source)
        and type(expected_context_count) is int and expected_context_count >= 0
        and type(expected_hypothesis_count) is int and expected_hypothesis_count >= 0)
    def report(event, **values):
        if progress is not None:
            progress(event, **values)
    report('primary_authentication')
    reader = EvidenceReader(evidence, source, scratch_directory=scratch)
    if progress is not None:
        original_rows = reader._rows
        def observed_rows(shard, number):
            yield from original_rows(shard, number)
            total = len(evidence.node_shards if shard.kind == 'evidence_nodes' else evidence.root_shards)
            if number % 250 == 0 or number + 1 == total:
                report('primary_authentication', kind=shard.kind, completed=number + 1, total=total)
        reader._rows = observed_rows
    with reader:
        root = reader.reference('case.primary', kind='scientific-primary-evidence-v1')
        _require(root.sha256 == expected_primary_sha256
            and reader.connection.execute('SELECT name FROM roots ORDER BY ordinal DESC LIMIT 1').fetchone() == ('case.primary',)
            and reader.connection.execute("SELECT 1 FROM roots WHERE name IN ('case','case.selected')").fetchone() is None,
            'retained primary is not the unsealed case boundary')
        with EvidenceWriter(output, scratch_directory=scratch) as writer:
            report('primary_copy')
            try:
                with writer.authenticated_prefix_scope(reader), ExitStack() as stack:
                    primary = reader.reconstruct(root)
                    _require(primary['schema'] == 'genosapiens.wgs-scientific-primary' and primary['version'] == 1
                        and primary['context_count'] == expected_context_count
                        and primary['hypothesis_count'] == expected_hypothesis_count
                        and all(type(primary[key]) is int and 0 <= primary[key] < 2**63 for key in _COUNTS))
                    resources_ref = _ref(reader, primary['resources'], 'scientific-case-resources-v1')
                    _require(resources_ref.name == 'case-resources')
                    report('retained_resources')
                    resources = stack.enter_context(VerifiedScientificCaseResources(reader, resources_ref,
                        expected_parent_bindings=expected_parent_bindings,
                        eligibility_root=eligibility_root, profile_root=profile_root,
                        authenticated_resources=authenticated_resources))
                    metadata_ref = reader.reference('case-reference', kind=METADATA_KIND)
                    metadata = stack.enter_context(VerifiedSubmissionReferenceMetadata(reader, metadata_ref,
                        expected_binding=expected_reference_binding))
                    _require(resources.reference_binding == expected_reference_binding)
                    entries, ranks = _selected_entries(reader, primary, stack)
                    ordinals = {row[1][0] for row in entries}
                    contexts, identities = _selected_directories(reader, primary, ordinals)
                    def restore(ordinal):
                        report('selected_context', ordinal=ordinal)
                        return VerifiedScientificCaseContext(reader,
                            _ref(reader, contexts[ordinal]['context'], 'scientific-case-context-v1'),
                            resources, metadata.reference, metadata.pin, metadata.renderers)
                    groups = _replay_selected_groups(entries, restore,
                        lambda locator: _identity(reader, identities, locator))
                    _require([group.source_rank for group in groups] == ranks)
                    primary_ref = dict(name=root.name, kind=root.kind, sha256=root.sha256)
                    selected_sha = writer.add('case.selected', 'scientific-selected-groups-v1', {
                        'primary': primary_ref, 'selection': primary['selection'],
                        'reference': resources.reference_binding, 'groups': groups})
                    header = dict(schema='genosapiens.wgs-complete-scientific-case', version=2,
                        resources=resources.reference, primary=primary_ref,
                        selected=dict(name='case.selected', kind='scientific-selected-groups-v1', sha256=selected_sha),
                        order=primary['order'], assessment_domain_sha256=resources.assessment_domain_sha256,
                        optional_analysis_status='not_run', **{key:primary[key] for key in _COUNTS})
                    case_sha = writer.add('case', 'complete-scientific-case-v1', header)
                    values = dict(case_name='case', case_sha256=case_sha,
                        resources_name=resources.name, resources_sha256=resources.sha256,
                        reference_name=metadata_ref.name, reference_sha256=metadata_ref.sha256,
                        parent_bindings=tuple(sorted(expected_parent_bindings.items())),
                        policy_id=resources.fusion_policy.policy_id, policy_sha256=resources.fusion_policy_sha256,
                        derivation_sha256=resources.fusion_policy.derivation_sha256,
                        **{key:primary[key] for key in _COUNTS})
                report('evidence_finish')
                saved = writer.finish()
            except BaseException:
                writer.failed = True
                raise
    artifact = ScientificArtifact(evidence=saved, **values)
    _write_descriptor(output, artifact.as_bytes())
    report('candidate_descriptor_written', sha256=artifact.sha256)
    return artifact
