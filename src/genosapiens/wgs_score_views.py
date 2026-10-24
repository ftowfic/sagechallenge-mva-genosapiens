"""Computed sensitivity views over the fixed complete primary score universe.

Primary pair assignments are reused losslessly. Every distinct score input is
actually evaluated with the changed phenotype; the resulting complete mapping
binds each original outcome to its computed view outcome. Diagnosis changes
only the ordering total through its separately source-validated context prior.
The case sensitivity owner binds ablation/diagnosis selection and resources.
"""
from array import array
from dataclasses import replace
from fractions import Fraction

from genosapiens.evidence_fusion import ranking_ties
from genosapiens.evidence_model import _exact
from genosapiens.evidence_serialization import _canonical
from genosapiens.track1_sensitivity import replace_phenotype_component
from genosapiens.wgs_evidence_replay import EvidenceReplay
from genosapiens.wgs_evidence_store import EvidenceStoreError
from genosapiens.wgs_score_decisions import ScoreDecisionCatalog

MAX_VIEW_RECORD_BYTES = 128 * 1024**2


class _MappedOutcomes:
    def __init__(self, primary, mapping):
        self.primary,self.mapping,self.count = primary,mapping,primary.buffer.count

    def _map(self, code):
        if self.primary.offset and code == 0:
            return 0
        return self.mapping[code-self.primary.offset]+self.primary.offset

    def values(self):
        for code in self.primary.buffer.values():
            yield self._map(code)

    def value_at(self, index):
        return self._map(self.primary.buffer.value_at(index))


class PairScoreView:
    def __init__(self, writer, primary, phenotype_score, *, name):
        self.primary,self.connection,self.name = primary,writer.connection,name
        self.catalog,self.closed,self.mapping = None,False,array("I")
        if primary.closed or primary.connection is not self.connection:
            raise EvidenceStoreError("foreign primary pair score view")
        self.universe,self.sources,self.population = primary.universe,primary.sources,primary.population
        self.offset,self.hypothesis_count = primary.offset,primary.hypothesis_count
        try:
            policy = primary.sources.policy
            policy.phenotype.component(phenotype_score)
            self.catalog = ScoreDecisionCatalog(writer,policy,policy_sha256=primary.sources.policy_sha256,name=name + ".decisions")
            for ordinal,inputs,_ in primary.catalog.records():
                if ordinal != len(self.mapping):
                    raise EvidenceStoreError("primary score view catalog order changed")
                self.mapping.append(self.catalog.evaluate(replace(inputs,phenotype_units=phenotype_score.ledger["raw_evidence_units"])))
            if len(self.mapping) != primary.catalog.count:
                raise EvidenceStoreError("pair score view coverage changed")
            catalog_sha = self.catalog.publish()
            mapping_sha = writer.add_sequence(name + ".mapping","pair-score-view-mapping-v1",iter(self.mapping))
            self.buffer = _MappedOutcomes(primary,self.mapping)
            self.sha256 = writer.add(name,"pair-score-view-v1",{
                "schema":"genosapiens.wgs-pair-score-view","version":1,
                "primary":{"name":primary.name,"sha256":primary.sha256},
                "phenotype":phenotype_score.ledger,
                "decisions":{"name":self.catalog.name,"sha256":catalog_sha},
                "mapping":{"name":name + ".mapping","sha256":mapping_sha},
                "mapping_count":len(self.mapping),"pair_count":self.buffer.count,
                "hypothesis_count":self.hypothesis_count,"offset":self.offset})
        except BaseException:
            writer.failed = True
            self.close()
            raise

    def close(self):
        if self.catalog is not None:
            self.catalog.close()
        self.mapping = array("I")
        self.closed = True


class CompleteContextScoreView:
    def __init__(self, writer, primary, *, phenotype_score, mode, view_binding, prior_contribution="0", name):
        self.primary,self.connection,self.name = primary,writer.connection,name
        self.closed,self.pairs,self.owns_pairs = False,None,False
        self.policy,self.policy_sha256,self.universe = primary.policy,primary.policy_sha256,primary.universe
        self.hypothesis_count,self.singleton_count = primary.hypothesis_count,primary.singleton_count
        self.phenotype_score,self.mode,self.prior_contribution = phenotype_score,mode,prior_contribution
        self.ties = ranking_ties(self.policy)
        if (primary.closed or primary.connection is not self.connection
                or mode not in ("dominant_term_sensitivity","diagnosis_informed_sensitivity")
                or type(prior_contribution) is not str or _exact(prior_contribution) < 0
                or mode == "dominant_term_sensitivity" and _exact(prior_contribution) != 0
                or mode == "diagnosis_informed_sensitivity" and phenotype_score.as_bytes() != primary.phenotype_score.as_bytes()):
            raise EvidenceStoreError("invalid fixed-universe context score view")
        try:
            self.policy.phenotype.component(phenotype_score)
            if len(_canonical((phenotype_score.ledger,view_binding))) > MAX_VIEW_RECORD_BYTES:
                raise EvidenceStoreError("score view evidence exceeds capacity",capacity=True)
            def records():
                count = 0
                for ledger in self.singleton_ledgers():
                    count += 1
                    yield {"hypothesis_id":ledger.hypothesis_id,"ledger":ledger,
                           "ordering_total":self.ordering_total(ledger.total)}
                if count != self.singleton_count:
                    raise EvidenceStoreError("singleton score view coverage changed")
            singletons_sha = writer.add_sequence(name + ".singletons","singleton-score-view-v1",records())
            if primary.pairs is not None:
                if mode == "dominant_term_sensitivity":
                    self.pairs = PairScoreView(writer,primary.pairs,phenotype_score,name=name + ".pairs")
                    self.owns_pairs = True
                else:
                    self.pairs = primary.pairs
            self.sha256 = writer.add(name,"complete-context-score-view-v1",{
                "schema":"genosapiens.wgs-complete-context-score-view","version":1,"mode":mode,
                "primary":{"name":primary.name,"sha256":primary.sha256},"binding":view_binding,
                "phenotype":phenotype_score.ledger,"prior_contribution":prior_contribution,
                "singletons":{"name":name + ".singletons","sha256":singletons_sha},
                "pairs":None if self.pairs is None else {"name":self.pairs.name,"sha256":self.pairs.sha256},
                "singleton_count":self.singleton_count,"hypothesis_count":self.hypothesis_count})
        except BaseException:
            writer.failed = True
            self.close()
            raise

    def _replace(self, evidence, ledger):
        if self.mode == "diagnosis_informed_sensitivity":
            return evidence,ledger
        return replace(evidence,phenotype_score=self.phenotype_score),replace_phenotype_component(
            ledger,self.phenotype_score,transform=self.policy.phenotype,ties=self.ties)

    def singletons(self):
        if self.closed:
            raise EvidenceStoreError("context score view is closed")
        for evidence,ledger in self.primary.singletons():
            yield self._replace(evidence,ledger)

    def singleton_ledgers(self):
        if self.closed:
            raise EvidenceStoreError("context score view is closed")
        for ledger in self.primary.singleton_ledgers():
            yield (ledger if self.mode == "diagnosis_informed_sensitivity" else
                   replace_phenotype_component(ledger,self.phenotype_score,
                       transform=self.policy.phenotype,ties=self.ties))

    def ordering_total(self, total):
        return str(Fraction(total)+Fraction(self.prior_contribution))

    def record_at(self, locator, *, positions=None):
        if self.closed:
            raise EvidenceStoreError("context score view is closed")
        if positions is None and self.universe.pairs is not None:
            positions = self.universe.pairs.domain.positions()
        constructed,evidence,ledger,assessed = self.primary.record_at(locator,positions=positions)
        evidence,ledger = self._replace(evidence,ledger)
        if locator[1]:
            code = self.pairs.buffer.value_at(positions.ordinal(*locator)-1)-self.pairs.offset
            decision = self.pairs.catalog.result(code)
            if (ledger.total != decision.total or tuple((c.name,c.raw,c.normalized) for c in ledger.components)
                    != tuple((c.name,c.raw,c.normalized) for c in decision.components)):
                raise EvidenceStoreError("selected sensitivity score does not replay")
        return constructed,evidence,ledger,assessed

    def close(self):
        if self.owns_pairs:
            self.pairs.close()
        self.closed = True

    def __enter__(self):
        return self

    def __exit__(self,*_):
        self.close()


class VerifiedContextScoreView(CompleteContextScoreView):
    def __init__(self, reader, reference, primary, **kwargs):
        if reference.kind != "complete-context-score-view-v1":
            raise EvidenceStoreError("foreign context score view root")
        reader._require_reference(reference)
        super().__init__(EvidenceReplay(reader),primary,name=reference.name,**kwargs)
