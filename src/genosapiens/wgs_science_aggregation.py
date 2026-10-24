"""Aggregate compatible completed context artifacts without hypothesis replay.

Only the checkpoint coordinator may supply a completed producer here. Ordinary
saved-result replay still reconstructs every scientific context independently.
"""
import json

from genosapiens.evidence_model import _typed
from genosapiens.wgs_admission import LedgerAdmissionInput
from genosapiens.wgs_case_evidence import CompleteScientificCase, CompletedScientificCase
from genosapiens.wgs_case_resources import bound_reference
from genosapiens.wgs_completed_science import child, fields, require
from genosapiens.wgs_context_domain import next_context_ordinal
from genosapiens.wgs_context_order import _identity_bytes
from genosapiens.wgs_context_evidence import MAX_MEMBER_BYTES
from genosapiens.wgs_evidence_store import EvidenceReader, EvidenceWriter
from genosapiens.wgs_locator_tape import _decode
from genosapiens.wgs_source_closure import _allele_digests


def _header(owner, ordinal, *, key='ordinal'):
    row = owner.connection.execute(f'SELECT content FROM {owner.prefix}_contexts WHERE {key}=?', (ordinal,)).fetchone()
    require(row is not None)
    return json.loads(row[0])


def _without(value, *keys):
    require(all(key in value for key in keys))
    return {key: item for key, item in value.items() if key not in keys}


def _selection_members(source, ordinal):
    count = 0
    for _, size, content in source.connection.execute(
            f'SELECT ordinal,count,content FROM {source.prefix}_members WHERE context=? ORDER BY ordinal', (ordinal,)):
        _, values = _decode(content, size, 3)
        for offset in range(0, len(values), 3):
            rendered, event, parent = values[offset:offset+3]
            value = None
            if rendered:
                row = source.connection.execute(f'SELECT content FROM {source.prefix}_rendered WHERE id=?', (rendered,)).fetchone()
                require(row is not None)
                value = row[0]
            yield value, event, parent
            count += 1
    require(count == _header(source, ordinal, key='id')['member_count'])


def _context_record_rows(reader, reference):
    """Decode context-specific rows, reusing the case's live resource rows."""
    from itertools import islice
    from genosapiens.wgs_evidence_store import MAX_RECONSTRUCT_BYTES
    entries = iter(reader._entries(reference.sha256))
    while batch := tuple(islice(entries, 512)):
        known = {row[0] for row in reader.connection.execute(
            'SELECT value_digest FROM live_prefix.science_resource_source_records WHERE value_digest IN ('
            + ','.join('?' for _ in batch) + ')', batch)}
        for digest in batch:
            if digest not in known:
                require(reader._lookup(digest)[1] <= MAX_RECONSTRUCT_BYTES)
                yield json.loads(b''.join(reader._parts(digest)))


def _source_events(reader, source, context, *, reuse_resource_records=False):
    closure = reader.reconstruct(bound_reference(reader, source.primary_header['source_closure'],
        'complete-scientific-source-closure-v1'))
    require(closure['context_count'] == 1 and closure['genes'] == source.resources.header['genes']
        and closure['phenotypes'] == source.resources.header['phenotypes'])
    for key, kind, event in (
        ('owners', 'scientific-source-owners-v1', 'owner'),
        ('records', 'scientific-source-records-v1', 'record'),
        ('queries', 'scientific-source-reference-queries-v1', 'query'),
        ('original_members', 'scientific-original-members-v1', 'original_member'),
    ):
        if key not in closure:
            require(key == 'original_members')
            continue
        reference = child(reader, closure[key], kind)
        rows = (_context_record_rows(reader, reference) if key == 'records' and reuse_resource_records
            else reader.iter_array(reference))
        for row in rows:
            yield (event, *row)
    # The portable closure stores each complete allele digest. Its cross-context
    # identity also normalizes only the known population batch envelope. Derive
    # that identity from authenticated member bytes with the ordinary function.
    members = reader.reconstruct(bound_reference(reader, context['members'], 'context-evidence-v1'))
    count = 0
    for row in reader.iter_array(child(reader, members['members'], 'context-enriched-members-v1'),
            max_item_bytes=MAX_MEMBER_BYTES + 5 * 1024**2):
        fields(row, ('admitted', 'enriched', 'batch'))
        enriched = _typed(row['enriched'], LedgerAdmissionInput, decode=True)
        require(enriched.allele is not None)
        yield 'allele', enriched.allele.allele_id, *_allele_digests(enriched.allele)
        count += 1
    require(count == context['member_count'])


def append_completed_context(case, reader, source):
    """Copy one producer's compact projections into a new case domain."""
    require(type(case) is CompleteScientificCase and type(case.writer) is EvidenceWriter
        and type(reader) is EvidenceReader and type(source) is CompletedScientificCase
        and source.connection is reader.connection and source.published and not source.closed
        and source.context_count == 1 and source.context_domain is not None
        and len(source.context_domain) == 1 and source.name.startswith('checkpoint-case-')
        and source.resources.reference == case.resources.reference)
    ordinal = source.context_domain.at(0)
    require(ordinal == next_context_ordinal(case))
    row = source.context(ordinal)
    try:
        # Each consumed source remains checked through the complete append.
        case.writer.merge_from(reader, exclude_prefix=source.name)
        with reader.read_into(case.writer):
            context = reader.reconstruct(bound_reference(reader, row['context'], 'scientific-case-context-v1'))
            require(context['resources'] == case.resources.reference)
            closure = reader.reconstruct(bound_reference(reader, source.primary_header['source_closure'],
                'complete-scientific-source-closure-v1'))
            contexts = iter(reader.iter_array(child(reader, closure['contexts'], 'scientific-source-contexts-v1')))
            closure_header = next(contexts, None)
            require(closure_header is not None and next(contexts, None) is None
                and closure_header['ordinal'] == ordinal and closure_header['context_id'] == row['context_id'])
            case.closure._append_context(closure_header, _source_events(reader, source, context,
                reuse_resource_records=(case.writer._copied_prefix_source is not None
                    and case.writer._copied_prefix_source is reader._live_prefix_source
                    and getattr(case.resources, '_source_prefix', None) is not None)))

            identities = source.identities
            case.identities._append_context(_without(_header(identities, ordinal), 'members', 'singletons'),
                identities.connection.execute(f'SELECT ordinal,id,native FROM {identities.prefix}_members WHERE context=? ORDER BY ordinal', (ordinal,)),
                identities.connection.execute(f'SELECT ordinal,id,subject FROM {identities.prefix}_singletons WHERE context=? ORDER BY ordinal', (ordinal,)))

            inputs = source.inputs
            case.inputs._append_context(_without(_header(inputs, ordinal, key='id'), 'member_map', 'singleton_map', 'root'),
                _selection_members(inputs, ordinal),
                inputs.connection.execute(f'SELECT ordinal,member,code FROM {inputs.prefix}_singletons WHERE context=? ORDER BY ordinal', (ordinal,)),
                ((number, json.loads(value)) for number, value in inputs.connection.execute(
                    f'SELECT id,content FROM {inputs.prefix}_classes WHERE context=? ORDER BY id', (ordinal,))),
                inputs.connection.execute(f'SELECT ordinal,count,content FROM {inputs.prefix}_pairs WHERE context=? ORDER BY ordinal', (ordinal,)))

            facts = source.scalar_facts
            def scalar_rows(kind):
                return ((number, json.loads(value)) for number, value in facts.connection.execute(
                    f'SELECT ordinal,content FROM {facts.prefix}_{kind} WHERE context=? ORDER BY ordinal', (ordinal,)))
            case.scalar_facts._append_context(_without(_header(facts, ordinal), 'singletons', 'members', 'classes'),
                scalar_rows('singletons'), scalar_rows('members'), scalar_rows('classes'),
                facts.connection.execute(f'SELECT ordinal,count,content FROM {facts.prefix}_pairs WHERE context=? ORDER BY ordinal', (ordinal,)))

            order = source.orders['primary']
            def families():
                for family in order.families():
                    locators = order.tape.rows(start=family['start'], count=family['count'])
                    yield family, ((_identity_bytes(source.identities.identity(locator)), *locator) for locator in locators)
            case.orders['primary']._append_context(order.context(ordinal), families())
            case._append_context_directory(row)
    except BaseException:
        case.writer.failed = True
        raise
