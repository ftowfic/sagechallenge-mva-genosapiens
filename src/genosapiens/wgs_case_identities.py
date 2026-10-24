"""Bounded identity lookup for all complete ordered hypothesis occurrences.

The private projection retains one member identity and one singleton identity,
not one record per pair. Pair IDs use the unchanged original hash expression
and the evaluated complete pair locator. Each portable context sequence binds
the actual computed source and hypothesis owners.
"""
import hashlib
from genosapiens.wgs_context_domain import next_context_ordinal, require_complete_domain
import json

from genosapiens.evidence_serialization import _canonical
from genosapiens.genotype_hypotheses import hypothesis_id_for_members
from genosapiens.wgs_evidence_store import EvidenceStoreError

MAX_IDENTITY_ROW_BYTES = 64 * 1024


def _member_identity_rows(context):
    for ordinal,member in enumerate(context.members.iter_members(allele_order=True),1):
        yield ordinal,member.allele.allele_id,member.allele.native_identity_id


def _singleton_identity_rows(context):
    for ordinal,value in enumerate(context.universe.singletons(),1):
        subjects = value.hypothesis.evidence_subject_ids
        if len(subjects) != 1:
            raise EvidenceStoreError("scientific singleton subject coverage changed")
        row = ordinal,value.hypothesis.hypothesis_id,subjects[0]
        if context.order.identity((ordinal,0)) != row[1]:
            raise EvidenceStoreError("scientific singleton identity changed")
        yield row


def _identity_header(context):
    return {"ordinal":context.ordinal,"context_id":context.context_id,"context":context.reference,
        "sample_name":context.facts.header.sample_name,
        "pair_context_id":None if context.universe.pairs is None else context.universe.pairs.context.context_id,
        "member_count":context.members.member_count,"singleton_count":context.universe.singleton_count,
        "hypothesis_count":context.hypothesis_count}


class CompleteCaseIdentities:
    def __init__(self, writer, *, name, context_domain=None):
        self.context_domain = context_domain
        self.writer,self.connection,self.name = writer,writer.connection,name
        self.prefix = "case_identities_" + hashlib.sha256(name.encode()).hexdigest()
        self.tables,self.closed,self.published = [],False,False
        self.context_count,self.member_count,self.singleton_count,self.hypothesis_count = 0,0,0,0
        try:
            for suffix,definition in (
                ("contexts","ordinal INTEGER PRIMARY KEY,id TEXT UNIQUE,content TEXT"),
                ("members","context INTEGER,ordinal INTEGER,id TEXT,native TEXT,PRIMARY KEY(context,ordinal),UNIQUE(context,id)"),
                ("singletons","context INTEGER,ordinal INTEGER,id TEXT,subject TEXT,PRIMARY KEY(context,ordinal),UNIQUE(context,id)"),
            ):
                table = self.prefix + "_" + suffix
                self.connection.execute(f"CREATE TABLE {table} ({definition})")
                self.tables.append(table)
        except BaseException:
            writer.failed = True
            self.close()
            raise

    @staticmethod
    def _bounded(value):
        if len(_canonical(value)) > MAX_IDENTITY_ROW_BYTES:
            raise EvidenceStoreError("scientific identity row exceeds capacity",capacity=True)
        return value

    def add_context(self, context):
        self._require_open()
        try:
            if (self.published or context.closed or context.connection is not self.connection
                    or context.ordinal != next_context_ordinal(self) or context.ordinal >= 2**32):
                raise EvidenceStoreError("foreign scientific identity context")
            self._append_context(_identity_header(context),_member_identity_rows(context),_singleton_identity_rows(context))
        except BaseException:
            self.writer.failed = True
            raise

    def _append_context(self, header, member_rows, singleton_rows):
        """Consume producing rows after the caller binds their context owner."""
        self._require_open()
        ordinal = header["ordinal"]
        if self.published or type(ordinal) is not int or ordinal != next_context_ordinal(self) or ordinal >= 2**32:
            raise EvidenceStoreError("foreign scientific identity context")
        try:
            name = self.name + ".context-" + str(ordinal)
            member_count = singleton_count = 0
            def members():
                nonlocal member_count
                for row in member_rows:
                    member_count += 1
                    row = self._bounded(row)
                    if row[0] != member_count:
                        raise EvidenceStoreError("scientific member identity order changed")
                    self.connection.execute(f"INSERT INTO {self.prefix}_members VALUES (?,?,?,?)",(ordinal,*row))
                    yield row
            member_sha = self.writer.add_sequence(name + ".members","scientific-member-identities-v1",members())
            def singletons():
                nonlocal singleton_count
                for row in singleton_rows:
                    singleton_count += 1
                    row = self._bounded(row)
                    if row[0] != singleton_count:
                        raise EvidenceStoreError("scientific singleton identity order changed")
                    self.connection.execute(f"INSERT INTO {self.prefix}_singletons VALUES (?,?,?,?)",(ordinal,*row))
                    yield row
            singleton_sha = self.writer.add_sequence(name + ".singletons","scientific-singleton-identities-v1",singletons())
            if member_count != header["member_count"] or singleton_count != header["singleton_count"]:
                raise EvidenceStoreError("scientific identity coverage changed")
            value = self._bounded({**header,
                "members":{"name":name + ".members","sha256":member_sha},
                "singletons":{"name":name + ".singletons","sha256":singleton_sha}})
            self.connection.execute(f"INSERT INTO {self.prefix}_contexts VALUES (?,?,?)",
                (ordinal,header["context_id"],_canonical(value)))
            self.context_count += 1
            self.member_count += member_count
            self.singleton_count += singleton_count
            self.hypothesis_count += header["hypothesis_count"]
        except BaseException:
            self.writer.failed = True
            raise

    def publish(self):
        require_complete_domain(self)
        self._require_open()
        if self.published:
            raise EvidenceStoreError("case scientific identities are already published")
        try:
            digest = self.writer.add_sequence(self.name + ".contexts","scientific-identity-contexts-v1",
                (json.loads(content) for content, in self.connection.execute(f"SELECT content FROM {self.prefix}_contexts ORDER BY ordinal")))
            self.sha256 = self.writer.add(self.name,"case-scientific-identities-v1",{
                "schema":"genosapiens.wgs-case-scientific-identities","version":1,
                "context_count":self.context_count,"member_count":self.member_count,
                "singleton_count":self.singleton_count,"hypothesis_count":self.hypothesis_count,
                "contexts":{"name":self.name + ".contexts","sha256":digest}})
            self.published = True
        except BaseException:
            self.writer.failed = True
            raise

    @property
    def reference(self):
        self._require_open()
        if not self.published:
            raise EvidenceStoreError("case scientific identities are incomplete")
        return {"name":self.name,"kind":"case-scientific-identities-v1","sha256":self.sha256}

    def identity(self, locator):
        """Caller first authenticates membership through its complete order."""
        self._require_open()
        if (not self.published or type(locator) is not tuple or len(locator) != 3
                or any(type(v) is not int for v in locator)):
            raise EvidenceStoreError("invalid scientific identity locator")
        context,left,right = locator
        if right == 0:
            row = self.connection.execute(f"SELECT id FROM {self.prefix}_singletons WHERE context=? AND ordinal=?",(context,left)).fetchone()
            if row is None:
                raise EvidenceStoreError("scientific singleton locator is missing")
            return row[0]
        row = self.connection.execute(f"SELECT content FROM {self.prefix}_contexts WHERE ordinal=?",(context,)).fetchone()
        if row is None:
            raise EvidenceStoreError("scientific identity context is missing")
        metadata = json.loads(row[0])
        if metadata["pair_context_id"] is None or not 1 <= left < right <= metadata["member_count"]:
            raise EvidenceStoreError("scientific pair identity locator is invalid")
        members = tuple(row[0] for row in self.connection.execute(
            f"SELECT id FROM {self.prefix}_members WHERE context=? AND ordinal IN (?,?) ORDER BY ordinal",(context,left,right)))
        if len(members) != 2:
            raise EvidenceStoreError("scientific pair identities are missing")
        return hypothesis_id_for_members(metadata["pair_context_id"],members,metadata["sample_name"])

    def subject_ids(self, locator):
        """Subjects of an occurrence authenticated by its complete order.

        Locus-call singletons retain their evidence subject even though they
        cannot supply a submitted allele representation.
        """
        self._require_open()
        if (not self.published or type(locator) is not tuple or len(locator) != 3
                or any(type(v) is not int or v < 0 for v in locator)):
            raise EvidenceStoreError("invalid scientific subject locator")
        context, left, right = locator
        if right == 0:
            row = self.connection.execute(f"SELECT subject FROM {self.prefix}_singletons WHERE context=? AND ordinal=?",
                (context, left)).fetchone()
            if row is None:
                raise EvidenceStoreError("scientific singleton subject is missing")
            return row
        if not 0 < left < right:
            raise EvidenceStoreError("invalid scientific pair subject locator")
        values = tuple(row[0] for row in self.connection.execute(
            f"SELECT id FROM {self.prefix}_members WHERE context=? AND ordinal IN (?,?) ORDER BY ordinal",
            (context, left, right)))
        if len(values) != 2:
            raise EvidenceStoreError("scientific pair subjects are missing")
        return values

    def _require_open(self):
        if self.closed or self.connection is None or self.writer.connection is not self.connection:
            raise EvidenceStoreError("case scientific identities are closed")

    def close(self):
        for table in reversed(self.tables):
            self.connection.execute("DROP TABLE " + table)
        self.tables.clear()
        self.closed = True

    def __enter__(self):
        return self

    def __exit__(self,*_):
        self.close()


class CompletedCaseIdentities(CompleteCaseIdentities):
    """Read the producer's compact identities without rebuilding hypotheses."""
    def __init__(self, reader, reference, *, context_at, progress=None, context_domain=None):
        from genosapiens.wgs_completed_science import child, fields, header, integer, require, observed_rows
        from genosapiens.wgs_evidence_replay import EvidenceReplay
        value = header(reader, reference, kind="case-scientific-identities-v1",
            schema="genosapiens.wgs-case-scientific-identities", keys=("context_count", "member_count",
                "singleton_count", "hypothesis_count", "contexts"))
        for key in ("context_count", "member_count", "singleton_count", "hypothesis_count"):
            integer(value[key])
        super().__init__(EvidenceReplay(reader), name=reference.name, context_domain=context_domain)
        try:
            for row in observed_rows(reader.iter_array(child(reader, value["contexts"], "scientific-identity-contexts-v1"),
                    max_item_bytes=MAX_IDENTITY_ROW_BYTES), progress, 'saved_identities', value['context_count']):
                fields(row, ("ordinal", "context_id", "context", "sample_name", "pair_context_id",
                    "member_count", "singleton_count", "hypothesis_count", "members", "singletons"))
                ordinal = integer(row["ordinal"], minimum=1, maximum=2**32-1)
                require(ordinal == next_context_ordinal(self))
                context = context_at(ordinal)
                require(all(row[key] == context[key] for key in
                    ("ordinal", "context_id", "context", "member_count", "hypothesis_count")))
                require(type(row["sample_name"]) is str and bool(row["sample_name"]))
                require(row["pair_context_id"] is None or type(row["pair_context_id"]) is str)
                for kind, root_kind, width in (("members", "scientific-member-identities-v1", 3),
                                               ("singletons", "scientific-singleton-identities-v1", 3)):
                    count = 0
                    for item in reader.iter_array(child(reader, row[kind], root_kind), max_item_bytes=MAX_IDENTITY_ROW_BYTES):
                        require(type(item) is list and len(item) == width)
                        count += 1
                        require(integer(item[0], minimum=1) == count and all(type(v) is str and v for v in item[1:]))
                        self.connection.execute(f"INSERT INTO {self.prefix}_{kind} VALUES (?,?,?,?)", (ordinal, *item))
                    require(count == integer(row["member_count" if kind == "members" else "singleton_count"]))
                self.connection.execute(f"INSERT INTO {self.prefix}_contexts VALUES (?,?,?)",
                    (ordinal, row["context_id"], _canonical(row)))
                self.context_count += 1
                for key in ("member_count", "singleton_count", "hypothesis_count"):
                    setattr(self, key, getattr(self, key) + integer(row[key]))
            require(all(getattr(self, key) == value[key] for key in
                ("context_count", "member_count", "singleton_count", "hypothesis_count")))
            self.sha256, self.published = reference.sha256, True
        except BaseException:
            self.close()
            raise
