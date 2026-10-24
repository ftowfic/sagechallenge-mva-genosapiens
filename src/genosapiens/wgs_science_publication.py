"""Private, bounded publication projections and one ordered final writer.

Spools belong only to this live invocation. They are never checkpoints or an
external replay format. Pickle carries producer-owned scalar tuples and bytes.
"""
import json
from pathlib import Path
import pickle
import sqlite3
import time

from genosapiens.wgs_science_aggregation import (
    _header, _without, _selection_members, _source_events,
)
from genosapiens.wgs_case_resources import bound_reference
from genosapiens.wgs_completed_science import child, require
from genosapiens.wgs_context_order import _identity_bytes
from genosapiens.wgs_context_domain import next_context_ordinal

MAX_SPOOL_BYTES = 256 * 1024**2
MAX_ROW_BYTES = 32 * 1024**2


def prepare_projection(reader, source, path):
    """Decode one committed artifact without rerunning scientific evaluation."""
    path = Path(path)
    require(not path.exists())
    db = sqlite3.connect(path)
    db.execute('PRAGMA cache_size=-8192')
    db.execute('CREATE TABLE items(section TEXT,ordinal INTEGER,content BLOB,PRIMARY KEY(section,ordinal)) WITHOUT ROWID')
    sizes, counts, phases = 0, {}, {}
    def rows(section, values):
        nonlocal sizes
        started = time.monotonic()
        count = 0
        for count, value in enumerate(values, 1):
            raw = pickle.dumps(value, protocol=5)
            sizes += len(raw)
            require(len(raw) <= MAX_ROW_BYTES and sizes <= MAX_SPOOL_BYTES)
            db.execute('INSERT INTO items VALUES (?,?,?)', (section, count, raw))
        counts[section] = count
        phases[section] = time.monotonic()-started
    try:
        ordinal = source.context_domain.at(0)
        row = source.context(ordinal)
        context = reader.reconstruct(bound_reference(reader, row['context'], 'scientific-case-context-v1'))
        closure = reader.reconstruct(bound_reference(reader, source.primary_header['source_closure'],
            'complete-scientific-source-closure-v1'))
        headers = tuple(reader.iter_array(child(reader, closure['contexts'], 'scientific-source-contexts-v1')))
        require(len(headers) == 1 and headers[0]['ordinal'] == ordinal)
        rows('directory', (row,))
        rows('closure_header', headers)
        rows('source_events', _source_events(reader, source, context,
            reuse_resource_records=reader._live_prefix_source is not None))
        identities = source.identities
        rows('identity_header', (_without(_header(identities, ordinal), 'members', 'singletons'),))
        rows('identity_members', identities.connection.execute(
            f'SELECT ordinal,id,native FROM {identities.prefix}_members WHERE context=? ORDER BY ordinal', (ordinal,)))
        rows('identity_singletons', identities.connection.execute(
            f'SELECT ordinal,id,subject FROM {identities.prefix}_singletons WHERE context=? ORDER BY ordinal', (ordinal,)))
        inputs = source.inputs
        rows('input_header', (_without(_header(inputs, ordinal, key='id'), 'member_map', 'singleton_map', 'root'),))
        rows('input_members', _selection_members(inputs, ordinal))
        rows('input_singletons', inputs.connection.execute(
            f'SELECT ordinal,member,code FROM {inputs.prefix}_singletons WHERE context=? ORDER BY ordinal', (ordinal,)))
        rows('input_classes', ((n, json.loads(v)) for n, v in inputs.connection.execute(
            f'SELECT id,content FROM {inputs.prefix}_classes WHERE context=? ORDER BY id', (ordinal,))))
        rows('input_pairs', inputs.connection.execute(
            f'SELECT ordinal,count,content FROM {inputs.prefix}_pairs WHERE context=? ORDER BY ordinal', (ordinal,)))
        facts = source.scalar_facts
        rows('fact_header', (_without(_header(facts, ordinal), 'singletons', 'members', 'classes'),))
        for kind in ('singletons', 'members', 'classes'):
            rows('fact_' + kind, ((n, json.loads(v)) for n, v in facts.connection.execute(
                f'SELECT ordinal,content FROM {facts.prefix}_{kind} WHERE context=? ORDER BY ordinal', (ordinal,))))
        rows('fact_pairs', facts.connection.execute(
            f'SELECT ordinal,count,content FROM {facts.prefix}_pairs WHERE context=? ORDER BY ordinal', (ordinal,)))
        order = source.orders['primary']
        rows('order_header', (order.context(ordinal),))
        def families():
            for index, family in enumerate(order.families()):
                rows('locators-' + str(index), ((_identity_bytes(source.identities.identity(locator)), *locator)
                    for locator in order.tape.rows(start=family['start'], count=family['count'])))
                yield index, family
        rows('families', families())
        db.commit()
        require(path.stat().st_size <= MAX_SPOOL_BYTES)
        return dict(bytes=path.stat().st_size, logical_bytes=sizes, counts=counts, phase_seconds=phases,
            ordinal=ordinal, source_case_sha256=source.sha256)
    finally:
        db.close()


def append_projection(case, reader, checkpoint, path, metadata):
    """Consume only a projection returned by this invocation's private workers."""
    require(metadata['ordinal'] == checkpoint.identity.original_ordinal == next_context_ordinal(case)
        and metadata['source_case_sha256'] == checkpoint.case_sha256
        and Path(path).stat().st_size == metadata['bytes'] <= MAX_SPOOL_BYTES)
    db = sqlite3.connect(Path(path).as_uri() + '?mode=ro', uri=True)
    db.execute('PRAGMA query_only=ON')
    def rows(section):
        count = 0
        for count, (raw,) in enumerate(db.execute(
                'SELECT content FROM items WHERE section=? ORDER BY ordinal', (section,)), 1):
            require(len(raw) <= MAX_ROW_BYTES)
            yield pickle.loads(raw)
        require(count == metadata['counts'][section])
    def one(section):
        values = tuple(rows(section))
        require(len(values) == 1)
        return values[0]
    started = time.monotonic()
    try:
        case.writer.merge_from(reader, exclude_prefix=checkpoint.case_name)
        merged = time.monotonic()
        with reader.read_into(case.writer):
            case.closure._append_context(one('closure_header'), rows('source_events'))
            case.identities._append_context(one('identity_header'), rows('identity_members'), rows('identity_singletons'))
            case.inputs._append_context(one('input_header'), rows('input_members'), rows('input_singletons'),
                rows('input_classes'), rows('input_pairs'))
            case.scalar_facts._append_context(one('fact_header'), rows('fact_singletons'), rows('fact_members'),
                rows('fact_classes'), rows('fact_pairs'))
            case.orders['primary']._append_context(one('order_header'),
                ((family, rows('locators-' + str(index))) for index, family in rows('families')))
            case._append_context_directory(one('directory'))
        return dict(evidence_merge_seconds=merged-started, projection_merge_seconds=time.monotonic()-merged)
    except BaseException:
        case.writer.failed = True
        raise
    finally:
        db.close()


def prepared_inputs(ordinals, *, workers, initialize, checkpoint_for, root, timeout):
    """At most one 256 MiB spool per reader, retained through ordered consume."""
    from genosapiens.wgs_science_workers import PersistentScienceWorkers
    root = Path(root) / 'publication-inputs'
    root.mkdir(mode=0o700)
    with PersistentScienceWorkers(workers, initialize) as pool:
        source = iter(ordinals)
        ended = False
        while not ended or pool.pending:
            while not ended and len(pool.pending) < workers:
                ordinal = next(source, None)
                if ordinal is None:
                    ended = True
                    break
                checkpoint = checkpoint_for(ordinal)
                pool.submit(dict(ordinal=ordinal, path=str(root / (str(ordinal) + '.sqlite')),
                    completion_sha256=checkpoint.completion_sha256), timeout=timeout)
            if pool.pending:
                ordinal, result = pool.take()
                if 'failure' in result:
                    raise ValueError('publication input processing failed')
                path = root / (str(ordinal) + '.sqlite')
                yield ordinal, path, result, pool.metrics
                path.unlink()
