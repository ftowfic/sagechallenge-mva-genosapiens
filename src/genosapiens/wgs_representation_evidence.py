"""Complete bounded submission representations from authenticated member facts."""
import hashlib
import json

from genosapiens.evidence_model import _typed
from genosapiens.evidence_serialization import _canonical
from genosapiens.genotype_properties import selected_source_calls
from genosapiens.track1_submission import assess_representation_source_facts, SubmissionRepresentationAssessment
from genosapiens.wgs_evidence_replay import EvidenceReplay
from genosapiens.wgs_evidence_store import EvidenceStoreError
from genosapiens.wgs_reference_evidence import ReferenceQueries, replay_query_source

MAX_REPRESENTATION_BYTES = 64 * 1024**2


class CompleteRepresentationEvidence:
    def __init__(self, writer, members, facts, reference, renderer_cache, reference_pin, *, name):
        self.writer,self.connection,self.name = writer,writer.connection,name
        self.members,self.facts = members,facts
        self.closed,self.created,self.count = False,False,0
        self.queries = None
        self.table = "representations_" + hashlib.sha256(name.encode()).hexdigest()
        try:
            if (members.closed or facts.closed or members.connection is not self.connection or facts.connection is not self.connection
                    or facts.context_id != members.context_id or facts.published is None
                    or facts.input_evidence != {"name":members.name,"sha256":members.sha256}
                    or reference.binding.contract_sha256 != reference_pin.contract_sha256
                    or reference.binding.manifest_sha256 != reference_pin.manifest_sha256):
                raise EvidenceStoreError("foreign representation source owner")
            self.connection.execute(f"CREATE TABLE {self.table} (id TEXT PRIMARY KEY,content TEXT,query_ordinal INTEGER)")
            self.created = True
            self.queries = ReferenceQueries(writer,reference,name=name+".queries",maximum_count=members.member_count)
            def rows():
                for member in members.iter_members():
                    allele = member.allele
                    def source_facts():
                        sources = {v.source_allele_id:v for v in facts.own_sources(allele.allele_id)}
                        calls = selected_source_calls(allele,facts.header,_sources=sources)
                        return calls,sources,facts.header.input_sha256
                    query_start = self.queries.count
                    result = assess_representation_source_facts(allele,source_facts,self.queries.opened,renderer_cache,reference_pin)
                    if self.queries.count - query_start not in (0, 1):
                        raise EvidenceStoreError('representation reference query scope changed')
                    content = _canonical(result)
                    if len(content) > MAX_REPRESENTATION_BYTES:
                        raise EvidenceStoreError("submission representation exceeds capacity",capacity=True)
                    self.connection.execute(f"INSERT INTO {self.table} VALUES (?,?,?)",
                        (allele.allele_id,content,query_start if self.queries.count > query_start else None))
                    self.count += 1
                    yield result
                if self.count != members.member_count:
                    raise EvidenceStoreError("submission member representation coverage changed")
            rows_sha = writer.add_sequence(name + ".members","member-submission-representations-v1",rows())
            self.queries.publish()
            self.sha256 = writer.add(name,"complete-submission-representations-v1",{
                "schema":"genosapiens.wgs-submission-representations","version":1,
                "members":{"name":members.name,"sha256":members.sha256},"genotype_facts":facts.published,
                "reference":reference.binding,"count":self.count,
                "queries":{"name":self.queries.name,"sha256":self.queries.sha256,"count":self.queries.count},
                "representations":{"name":name + ".members","sha256":rows_sha}})
        except BaseException:
            writer.failed = True
            self.close()
            raise

    def get(self, key, default=None):
        if self.closed:
            raise EvidenceStoreError("submission representations are closed")
        row = self.connection.execute(f"SELECT content FROM {self.table} WHERE id=?",(key,)).fetchone()
        return default if row is None else _typed(json.loads(row[0]),SubmissionRepresentationAssessment,decode=True)

    def records(self):
        if self.closed:
            raise EvidenceStoreError("submission representations are closed")
        for content, in self.connection.execute(f"SELECT content FROM {self.table} ORDER BY id"):
            yield _typed(json.loads(content),SubmissionRepresentationAssessment,decode=True)

    def close(self):
        if self.queries is not None:
            self.queries.close()
        if self.created:
            self.connection.execute("DROP TABLE " + self.table)
            self.created = False
        self.closed = True

    def __enter__(self):
        return self

    def __exit__(self,*_):
        self.close()


class VerifiedRepresentationEvidence(CompleteRepresentationEvidence):
    def __init__(self, reader, root, members, facts, reference, renderer_cache, reference_pin):
        if root.kind != "complete-submission-representations-v1":
            raise EvidenceStoreError("foreign submission representation root")
        reader._require_reference(root)
        value = reader.reconstruct(root)
        query = value.get("queries")
        if type(query) is not dict or set(query) != {"name","sha256","count"}:
            raise EvidenceStoreError("invalid retained representation queries")
        source = reader.reference(query["name"],kind="submission-reference-queries-v1")
        if source.sha256 != query["sha256"]:
            raise EvidenceStoreError("retained representation query binding changed")
        opened,finish,rows = replay_query_source(reader,source,reference)
        try:
            super().__init__(EvidenceReplay(reader),members,facts,opened,renderer_cache,reference_pin,name=root.name)
            finish()
        except BaseException:
            if hasattr(self,"created"):
                self.close()
            raise
        finally:
            rows.close()
