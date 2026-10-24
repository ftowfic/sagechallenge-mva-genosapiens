"""Complete singleton score records and the computed complete pair score tape.

This component accepts the source adapter's already computed phenotype result.
The case/publication owner must bind observations, public resource closure and
the primary or sensitivity view that produced it; this is not source admission.
"""
import hashlib
import json

from genosapiens.evidence_fusion import (
    FusionPolicy, HypothesisEvidence, _score, _inheritance_charges, disease_context_charges,
    disease_gene_context, population_assessment_charges, policy_as_json,
)
from genosapiens.evidence_model import (
    EvidenceReference, ScoreLedger, _typed, _construct_owned, _validation_session,
)
from genosapiens.evidence_serialization import _canonical
from genosapiens.wgs_evidence_replay import EvidenceReplay
from genosapiens.wgs_evidence_store import EvidenceStoreError
from genosapiens.wgs_pair_scores import CompletePairScores, VerifiedPairScores
from genosapiens.wgs_score_sources import CompletePairScoreSources, VerifiedPairScoreSources, member_score_charges
from genosapiens.genotype_evaluation import prepare_gene_evidence
from genosapiens.wgs_gene_evidence import CaseGeneEvidence

MAX_SCORE_RECORD_BYTES = 128 * 1024**2


class CompleteContextScores:
    def __init__(self, writer, universe, population, *, policy, policy_sha256, gene_evidence,
                 phenotype_score, phenotype_refs=(), name, gene_resource=None, phenotype_resource=None):
        self.writer,self.connection,self.name = writer,writer.connection,name
        self.universe,self.population = universe,population
        self.closed,self.created,self.sources,self.pairs = False,False,None,None
        self.table = "context_scores_" + hashlib.sha256(name.encode()).hexdigest()
        self.singleton_count = 0
        self.policy = _typed(policy,FusionPolicy)
        self.policy_sha256 = policy_sha256
        self.phenotype_score = phenotype_score
        self.phenotype_resource = phenotype_resource
        self.phenotype_refs = _typed(phenotype_refs,tuple[EvidenceReference,...])
        if (hashlib.sha256(policy_as_json(self.policy)).hexdigest() != policy_sha256
                or universe.connection is not writer.connection or population.connection is not writer.connection
                or population.universe is not universe or universe.closed or population.closed):
            raise EvidenceStoreError("foreign context score owners")
        try:
            if gene_resource is None:
                prepared = (universe.prepared_gene_evidence if gene_evidence is universe.prepared_gene_evidence.evidence
                            else prepare_gene_evidence(gene_evidence))
                gene_resource = CaseGeneEvidence(writer,prepared,name=name + ".genes")
                gene_evidence = prepared.evidence
            gene_resource.require(writer,gene_evidence)
            self.gene_resource = gene_resource
            if phenotype_resource is not None:
                phenotype_resource.require(writer,universe.members.seed.context,phenotype_score,self.phenotype_refs)
                if tuple(sorted(o.observation_id for o in phenotype_resource.resources.view.observations)) != tuple(sorted(universe.observation_ids)):
                    raise EvidenceStoreError("phenotype observations differ from hypothesis universe")
            self.policy.phenotype.component(phenotype_score)
            self.connection.execute(f"CREATE TABLE {self.table} (id TEXT PRIMARY KEY,content TEXT)")
            self.created = True
            shared,audits = disease_context_charges(universe.members.seed.context,gene_evidence,self.policy,
                                                   prepared_gene_evidence=gene_resource.prepared)
            self.shared = shared
            shared_context = disease_gene_context(universe.members.seed.context)
            shared_value = {"charges":shared,"audits":audits,"gene_evidence":gene_resource.reference,
                            "phenotype":phenotype_score.ledger,"phenotype_refs":self.phenotype_refs,
                            "phenotype_evidence":None if phenotype_resource is None else phenotype_resource.reference}
            if len(_canonical(shared_value)) > MAX_SCORE_RECORD_BYTES:
                raise EvidenceStoreError("context score shared evidence exceeds capacity",capacity=True)
            shared_sha = writer.add(name + ".shared","context-score-shared-v1",shared_value)
            def records():
                for constructed in universe.singletons():
                    hypothesis = constructed.hypothesis
                    if disease_gene_context(hypothesis.context) != shared_context:
                        raise EvidenceStoreError("singleton disease score dependency changed")
                    member, = constructed.model.members
                    qualification, = constructed.members
                    key, = hypothesis.evidence_subject_ids
                    member_charges = member_score_charges(universe.members,key,member,qualification,self.policy)
                    allele = universe.members.member(key).allele
                    assessed = population.assessment(hypothesis)
                    # One synchronous hypothesis owns the copied evidence.
                    # The scorer still performs its ordinary canonical and
                    # semantic checks, reusing only this operation's validated
                    # children. Release the walk before yielding any row.
                    with _validation_session(isolated=True):
                        evidence = _construct_owned(HypothesisEvidence,hypothesis,phenotype_score,(
                            *shared,*member_charges,*_inheritance_charges(constructed,self.policy),
                            *population_assessment_charges(assessed,{key:allele},population.policy.policy_sha256)),
                            self.phenotype_refs,None)
                        ledger = _score(evidence,self.policy,policy_sha256)
                    content = _canonical({"hypothesis_id":hypothesis.hypothesis_id,"evidence":evidence,"ledger":ledger})
                    if len(content) > MAX_SCORE_RECORD_BYTES:
                        raise EvidenceStoreError("singleton score record exceeds capacity",capacity=True)
                    self.connection.execute(f"INSERT INTO {self.table} VALUES (?,?)",(hypothesis.hypothesis_id,content))
                    self.singleton_count += 1
                    yield json.loads(content)
            singletons_sha = writer.add_sequence(name + ".singletons","singleton-score-evidence-v1",records())
            if self.singleton_count != universe.singleton_count:
                raise EvidenceStoreError("singleton score coverage changed")
            if universe.pairs is not None:
                self.sources,self.pairs = self._prepare_pairs(writer,universe,population,gene_evidence)
            self.hypothesis_count = self.singleton_count + (0 if self.pairs is None else self.pairs.hypothesis_count)
            if self.hypothesis_count != universe.hypothesis_count:
                raise EvidenceStoreError("complete context score coverage changed")
            self.sha256 = writer.add(name,"complete-context-scores-v1",{
                "schema":"genosapiens.wgs-complete-context-scores","version":1,
                "hypotheses":{"name":universe.name,"sha256":universe.sha256},
                "population":{"name":population.name,"sha256":population.sha256},
                "policy":self.policy,"policy_sha256":policy_sha256,
                "shared":{"name":name + ".shared","sha256":shared_sha},
                "singletons":{"name":name + ".singletons","sha256":singletons_sha},
                "pairs":None if self.pairs is None else {"name":self.pairs.name,"sha256":self.pairs.sha256},
                "singleton_count":self.singleton_count,"hypothesis_count":self.hypothesis_count})
        except BaseException:
            writer.failed = True
            self.close()
            raise

    def _prepare_pairs(self, writer, universe, population, gene_evidence):
        sources = CompletePairScoreSources(writer,universe,policy=self.policy,policy_sha256=self.policy_sha256,
                                          gene_evidence=gene_evidence,name=self.name + ".pair_sources",gene_resource=self.gene_resource)
        try:
            pairs = CompletePairScores(writer,universe,population,sources,
                                       phenotype_score=self.phenotype_score,name=self.name + ".pairs")
        except BaseException:
            sources.close()
            raise
        return sources,pairs

    def singletons(self):
        if self.closed:
            raise EvidenceStoreError("context score owner is closed")
        for content, in self.connection.execute(f"SELECT content FROM {self.table} ORDER BY id"):
            row = json.loads(content)
            yield _typed(row["evidence"],HypothesisEvidence,decode=True),_typed(row["ledger"],ScoreLedger,decode=True)

    def singleton_ledgers(self):
        """Project computed ledgers without rebuilding unused rich evidence.

        This private table is populated only by full scoring above, including
        during independent saved replay. Rich consumers keep singletons().
        """
        if self.closed:
            raise EvidenceStoreError("context score owner is closed")
        for key,content in self.connection.execute(f"SELECT id,content FROM {self.table} ORDER BY id"):
            row = json.loads(content)
            ledger = _typed(row["ledger"],ScoreLedger,decode=True)
            if not key == row["hypothesis_id"] == row["evidence"]["hypothesis"]["hypothesis_id"] == ledger.hypothesis_id:
                raise EvidenceStoreError("singleton ledger identity changed")
            yield ledger

    def record_at(self, locator, *, positions=None):
        """Exact rich source evidence and ledger for one already scored locator."""
        if self.closed:
            raise EvidenceStoreError("context score owner is closed")
        if positions is None and self.universe.pairs is not None:
            positions = self.universe.pairs.domain.positions()
        constructed = self.universe.constructed_at(locator,positions=positions)
        hypothesis = constructed.hypothesis
        assessed = self.population.assessment(hypothesis)
        if locator[1] == 0:
            content, = self.connection.execute(f"SELECT content FROM {self.table} WHERE id=?",(hypothesis.hypothesis_id,)).fetchone()
            row = json.loads(content)
            evidence,ledger = _typed(row["evidence"],HypothesisEvidence,decode=True),_typed(row["ledger"],ScoreLedger,decode=True)
        else:
            if self.pairs is None or self.sources is None:
                raise EvidenceStoreError("pair score evidence is unavailable")
            alleles = {key:self.universe.members.member(key).allele for key in hypothesis.evidence_subject_ids}
            with _validation_session(isolated=True):
                evidence = _construct_owned(HypothesisEvidence,hypothesis,self.phenotype_score,(
                    *self.sources.shared,*(c for key in hypothesis.evidence_subject_ids for c in self.sources.charges(key)),
                    *_inheritance_charges(constructed,self.policy),
                    *population_assessment_charges(assessed,alleles,self.population.policy.policy_sha256)),self.phenotype_refs,None)
                if len(_canonical(evidence)) > MAX_SCORE_RECORD_BYTES:
                    raise EvidenceStoreError("selected pair score record exceeds capacity",capacity=True)
                ledger = _score(evidence,self.policy,self.policy_sha256)
            code = self.pairs.buffer.value_at(positions.ordinal(*locator)-1)-self.pairs.offset
            decision = self.pairs.catalog.result(code)
            if (ledger.total != decision.total or tuple((c.name,c.raw,c.normalized) for c in ledger.components)
                    != tuple((c.name,c.raw,c.normalized) for c in decision.components)):
                raise EvidenceStoreError("selected pair score does not replay computed outcome")
        if evidence.hypothesis != hypothesis:
            raise EvidenceStoreError("selected score hypothesis changed")
        return constructed,evidence,ledger,assessed

    def close(self):
        if self.pairs is not None:
            self.pairs.close()
        if self.sources is not None:
            self.sources.close()
        if self.created:
            self.connection.execute("DROP TABLE " + self.table)
            self.created = False
        self.closed = True

    def __enter__(self):
        return self

    def __exit__(self,*_):
        self.close()


class VerifiedContextScores(CompleteContextScores):
    def __init__(self, reader, reference, universe, population, *, policy, policy_sha256,
                 gene_evidence, phenotype_score, phenotype_refs=(), gene_resource=None, phenotype_resource=None):
        if reference.kind != "complete-context-scores-v1":
            raise EvidenceStoreError("foreign context score root")
        reader._require_reference(reference)
        super().__init__(EvidenceReplay(reader),universe,population,policy=policy,policy_sha256=policy_sha256,
            gene_evidence=gene_evidence,phenotype_score=phenotype_score,phenotype_refs=phenotype_refs,name=reference.name,
            gene_resource=gene_resource,phenotype_resource=phenotype_resource)

    def _prepare_pairs(self, writer, universe, population, gene_evidence):
        reader = writer.reader
        sources = VerifiedPairScoreSources(reader,reader.reference(self.name + ".pair_sources",kind="pair-score-source-evidence-v1"),
            universe,policy=self.policy,policy_sha256=self.policy_sha256,gene_evidence=gene_evidence,gene_resource=self.gene_resource)
        try:
            pairs = VerifiedPairScores(reader,reader.reference(self.name + ".pairs",kind="complete-pair-scores-v1"),
                universe,population,sources,phenotype_score=self.phenotype_score)
        except BaseException:
            sources.close()
            raise
        return sources,pairs
