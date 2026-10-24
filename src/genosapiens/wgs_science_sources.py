"""Live resource source registries for independent scientific writers.

The ordinary closure authenticates resource records and references once. This
owner retains no persisted authority. External readers still rebuild closure.
"""
from genosapiens.wgs_context_domain import ContextDomain
from genosapiens.wgs_case_resources import ScientificCaseResources
from genosapiens.wgs_evidence_store import EvidenceStoreError, EvidenceWriter, FANOUT
from genosapiens.wgs_source_closure import CompleteSourceClosure

TABLES = ('science_resource_source_owners', 'science_resource_source_records')
_SEAL = object()


def _require(value):
    if not value:
        raise EvidenceStoreError('scientific resource source owner changed')


def _version(connection):
    return (connection.total_changes,
        connection.execute('PRAGMA schema_version').fetchone()[0],
        connection.execute('PRAGMA data_version').fetchone()[0])


class ScientificSourceTemplate:
    """Prepare bounded SQLite rows and their ordinary normalized node hashes."""
    tables = TABLES

    def __init__(self, writer, resources):
        _require(type(writer) is EvidenceWriter and type(resources) is ScientificCaseResources)
        self.writer, self.resources, self.closed = writer, resources, False
        self.created = []
        resources.require(writer)
        self.owners = self._owners()
        try:
            with CompleteSourceClosure(writer, resources.genes, resources.phenotypes,
                    name='resource-source-preparation', context_domain=ContextDomain.from_ordinals(())) as closure:
                owners, records = self.tables
                writer.connection.execute(f'CREATE TABLE {owners}(id TEXT PRIMARY KEY,manifest TEXT,role TEXT)')
                self.created.append(owners)
                writer.connection.execute(f'INSERT INTO {owners} SELECT * FROM {closure.prefix}_owners')
                writer.connection.execute(f'CREATE TABLE {records}(owner TEXT,id TEXT,digest TEXT,value_digest TEXT NOT NULL,'
                    'PRIMARY KEY(owner,id))')
                self.created.append(records)
                cursor = writer.connection.execute(f'SELECT * FROM {closure.prefix}_records ORDER BY owner,id')
                try:
                    while rows := cursor.fetchmany(FANOUT):
                        digests = writer._scalar_tuple_batch(rows)
                        writer.connection.executemany(f'INSERT INTO {records} VALUES (?,?,?,?)',
                            ((*row, digest) for row, digest in zip(rows, digests, strict=True)))
                finally:
                    cursor.close()
            writer.connection.execute(f'CREATE INDEX science_resource_value_digest ON {records}(value_digest)')
            # These tables are producer-owned constants in every copied writer.
            # Ordinary case collectors use separate mutable tables.
            for table in self.tables:
                for operation in ('INSERT', 'UPDATE', 'DELETE'):
                    writer.connection.execute(f'CREATE TRIGGER {table}_{operation.lower()} BEFORE {operation} ON {table} '
                        "BEGIN SELECT RAISE(ABORT,'scientific resource source table is immutable'); END")
            writer.connection.commit()
            self.version = _version(writer.connection)
        except BaseException:
            self.close()
            raise

    def _owners(self):
        r = self.resources
        return (r, r.genes, r.genes.prepared, r.genes.prepared.evidence,
            r.phenotypes, r.phenotypes.profile, r.phenotypes.observations, r.phenotypes.view)

    def require(self):
        _require(not self.closed and self.writer.connection is not None and not self.writer.failed)
        current = self._owners()
        _require(len(current) == len(self.owners) and
            all(a is b for a, b in zip(current, self.owners, strict=True))
            and _version(self.writer.connection) == self.version)
        self.resources.require(self.writer)
        self.resources.phenotypes.logical_hashes(self.writer)

    def bind(self, writer, resources):
        self.require()
        _require(type(writer) is EvidenceWriter and writer is not self.writer
            and type(resources) is ScientificCaseResources)
        resources.require(writer)
        _require(resources.genes.prepared is self.resources.genes.prepared
            and resources.phenotypes.profile is self.resources.phenotypes.profile
            and resources.phenotypes.view is self.resources.phenotypes.view
            and resources.phenotypes.observations is self.resources.phenotypes.observations)
        return _BoundSourcePrefix(self, writer, resources, _seal=_SEAL)

    def close(self):
        if not self.closed:
            self.closed = True
            for table in reversed(self.created):
                self.writer.connection.execute('DROP TABLE ' + table)
            self.created.clear()

    def __enter__(self):
        return self

    def __exit__(self, *_):
        self.close()


class _BoundSourcePrefix:
    def __init__(self, source, writer, resources, *, _seal):
        _require(_seal is _SEAL and type(source) is ScientificSourceTemplate)
        self.source, self.writer, self.resources = source, writer, resources
        self.closure = None

    def _require(self, closure):
        self.source.require()
        _require(type(closure) is CompleteSourceClosure and closure.writer is self.writer
            and closure.genes is self.resources.genes and closure.phenotypes is self.resources.phenotypes
            and self.writer.connection is not None and not self.writer.failed)
        self.resources.require(self.writer)

    def populate(self, closure):
        self._require(closure)
        _require(self.closure is None)
        self.closure = closure
        owners, records = self.source.tables
        connection = self.writer.connection
        connection.execute(f'INSERT INTO {closure.prefix}_owners SELECT * FROM {owners}')
        connection.execute(f'INSERT INTO {closure.prefix}_records SELECT owner,id,digest FROM {records}')

    def publish_records(self, closure, name):
        self._require(closure)
        _require(self.closure is closure)
        writer = self.writer
        def chunks():
            # Match the full row, so a local record cannot borrow another row's
            # encoded authority. Context records keep ordinary serialization.
            query = (f'SELECT c.owner,c.id,c.digest,b.value_digest FROM {closure.prefix}_records c '
                f'LEFT JOIN {self.source.tables[1]} b ON c.owner=b.owner AND c.id=b.id AND c.digest=b.digest '
                'ORDER BY c.owner,c.id')
            cursor = writer.connection.execute(query)
            try:
                while rows := cursor.fetchmany(FANOUT):
                    ordinary = iter(writer._scalar_tuple_batch([(owner, key, digest)
                        for owner, key, digest, prepared in rows if prepared is None]))
                    batch = [prepared if prepared is not None else next(ordinary)
                        for _, _, _, prepared in rows]
                    yield writer._node('array', batch)
            finally:
                cursor.close()
        def build():
            from itertools import chain
            iterator = iter(chunks())
            first = next(iterator, None)
            return writer._node('array', []) if first is None else writer._chunks(chain((first,), iterator), 'arrays')
        return writer._add(name, 'scientific-source-records-v1', build)
