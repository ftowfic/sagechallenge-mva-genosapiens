"""Source-bound case sensitivity choices over the complete primary universe.

The original phenotype and native/equivalence kernels remain authoritative.
Case resources are shared; context evidence binds each computed result to those
resources and the already computed primary scores. No saved validation flag is
accepted in place of replaying the source decision.
"""
from dataclasses import replace

from genosapiens.evidence_model import _exact, _typed
from genosapiens.evidence_serialization import _canonical
from genosapiens.phenotype_scoring import prepare_phenotype_scores
from genosapiens.track1_sensitivity import (
    DisclosedDiagnosis, PriorPolicy, prepare_diagnosis_match_context,
    prior_matches_context, select_dominant_term,
)
from genosapiens.wgs_evidence_store import EvidenceStoreError
from genosapiens.wgs_evidence_replay import EvidenceReplay
from genosapiens.wgs_phenotype_evidence import MAX_CONTEXT_PHENOTYPE_BYTES, _hash


def _primary(writer, scores):
    if (scores.closed or scores.connection is not writer.connection
            or getattr(scores,"mode",None) is not None):
        raise EvidenceStoreError("sensitivity requires the completed primary context")


def _bounded(value):
    if len(_canonical(value)) > MAX_CONTEXT_PHENOTYPE_BYTES:
        raise EvidenceStoreError("sensitivity source evidence exceeds capacity",capacity=True)
    return value


class CaseDiagnosisEvidence:
    def __init__(self, writer, genes, *, diagnosis, prior_policy, name):
        self.writer,self.connection,self.name = writer,writer.connection,name
        self.genes = genes
        try:
            genes.require(writer,genes.prepared.evidence)
            self.diagnosis = _typed(diagnosis,DisclosedDiagnosis | None)
            self.prior_policy = _typed(prior_policy,PriorPolicy)
            self.prepared = None if self.diagnosis is None else prepare_diagnosis_match_context(
                disease_gene_evidence=genes.prepared.evidence,diagnosis=self.diagnosis)
            self.sha256 = writer.add(name,"case-diagnosis-evidence-v1",{
                "schema":"genosapiens.wgs-case-diagnosis-evidence","version":1,
                "genes":genes.reference,"diagnosis":diagnosis,"prior_policy":prior_policy})
        except BaseException:
            writer.failed = True
            raise

    @property
    def reference(self):
        return {"name":self.name,"kind":"case-diagnosis-evidence-v1","sha256":self.sha256}


class ContextDiagnosisEvidence:
    def __init__(self, writer, resources, scores, *, name):
        self.connection,self.name,self.resources = writer.connection,name,resources
        try:
            _primary(writer,scores)
            resources.genes.require(writer,scores.gene_resource.prepared.evidence)
            if resources.connection is not self.connection:
                raise EvidenceStoreError("foreign case diagnosis owner")
            context = scores.universe.members.seed.context
            # Construction always copies exactly these context references into
            # every hypothesis. Model/ID/completeness/phase are not inputs to
            # prior_matches_context; no allele-dependent prior is factored out.
            self.proofs = () if resources.diagnosis is None else prior_matches_context(
                context,context.evidence_refs,resources.genes.prepared.evidence,
                resources.diagnosis,prepared=resources.prepared)
            self.contribution = str(_exact(resources.prior_policy.strength)) if self.proofs else "0"
            self.reason = ("source_native_or_equivalent" if self.proofs else "diagnosis_not_supplied"
                           if resources.diagnosis is None else "no_native_or_equivalent_match")
            self.primary = {"name":scores.name,"sha256":scores.sha256}
            self.sha256 = writer.add(name,"context-diagnosis-evidence-v1",_bounded({
                "schema":"genosapiens.wgs-context-diagnosis-evidence","version":1,
                "resources":resources.reference,"primary":self.primary,"context_sha256":_hash(context),
                "eligible":bool(self.proofs),"proofs":self.proofs,"reason":self.reason,
                "contribution":self.contribution,"hypothesis_count":scores.hypothesis_count}))
        except BaseException:
            writer.failed = True
            raise

    @property
    def reference(self):
        return {"name":self.name,"kind":"context-diagnosis-evidence-v1","sha256":self.sha256}

    def view_arguments(self, writer, scores):
        _primary(writer,scores)
        if self.connection is not writer.connection or self.primary != {"name":scores.name,"sha256":scores.sha256}:
            raise EvidenceStoreError("diagnosis context score binding changed")
        return dict(mode="diagnosis_informed_sensitivity",phenotype_score=scores.phenotype_score,
                    view_binding=self.reference,prior_contribution=self.contribution)


class CaseAblationEvidence:
    def __init__(self, writer, phenotypes, primary_order, leader_scores, *, prior_policy, name):
        self.writer,self.connection,self.name = writer,writer.connection,name
        self.phenotypes,self.order,self.prior_policy = phenotypes,primary_order,_typed(prior_policy,PriorPolicy)
        try:
            if (phenotypes.connection is not self.connection or primary_order.connection is not self.connection
                    or primary_order.closed or not primary_order.published):
                raise EvidenceStoreError("foreign primary ablation owner")
            locator = next(primary_order.tape.rows(start=0,count=min(1,primary_order.hypothesis_count)),None)
            leader = None
            if locator is not None:
                _primary(writer,leader_scores)
                row = primary_order.context(locator[0])
                if (row["scores"] != {"name":leader_scores.name,"sha256":leader_scores.sha256}
                        or row["context_id"] != leader_scores.universe.members.context_id
                        or leader_scores.phenotype_resource is None
                        or leader_scores.phenotype_resource.resources is not phenotypes):
                    raise EvidenceStoreError("primary leader phenotype binding changed")
                leader_scores.phenotype_resource.require(writer,leader_scores.universe.members.seed.context,
                                                        leader_scores.phenotype_score,leader_scores.phenotype_refs)
                constructed,_,_,_ = leader_scores.record_at(locator[1:])
                leader = {"locator":locator,"hypothesis_id":constructed.hypothesis.hypothesis_id,"scores":row["scores"]}
            elif leader_scores is not None:
                raise EvidenceStoreError("empty primary has a foreign leader")
            self.selection = None if leader is None else select_dominant_term(leader_scores.phenotype_score)
            self.term = None if self.selection is None else self.selection["canonical_hpo_id"]
            self.reason = "empty_primary_universe" if leader is None else "no_positive_support" if self.term is None else None
            self.status = "not_applicable" if self.reason else "applied"
            self.removed = ()
            self.view = phenotypes.view
            self.prepared = None
            if self.term is not None:
                self.removed = tuple(sorted(o.observation_id for o in self.view.observations
                    if o.status == "present" and phenotypes.ic_store.graph.resolve(o.source_hpo_id).canonical_id == self.term))
                self.view = replace(self.view,
                    observations=tuple(o for o in self.view.observations if o.observation_id not in self.removed),
                    context_provenance={k:v for k,v in self.view.context_provenance.items() if k not in self.removed})
                self.prepared = prepare_phenotype_scores(self.view.observations,ic_store=phenotypes.ic_store,
                    profile_store=None if phenotypes.profile is None else phenotypes.profile.store,
                    prior_identity="dominant-term:" + self.term + ":" + _hash(prior_policy))
            self.sha256 = writer.add(name,"case-ablation-evidence-v1",_bounded({
                "schema":"genosapiens.wgs-case-ablation-evidence","version":1,
                "resources":phenotypes.reference,"primary_order":primary_order.reference,"leader":leader,
                "prior_policy":prior_policy,"status":self.status,"reason":self.reason,
                "dominant_term":self.selection,"removed_observation_ids":self.removed,"view":self.view}))
        except BaseException:
            writer.failed = True
            raise

    @property
    def reference(self):
        return {"name":self.name,"kind":"case-ablation-evidence-v1","sha256":self.sha256}


class ContextAblationEvidence:
    def __init__(self, writer, resources, scores, *, name):
        self.name,self.connection,self.resources = name,writer.connection,resources
        try:
            _primary(writer,scores)
            if (resources.connection is not self.connection or scores.phenotype_resource is None
                    or scores.phenotype_resource.resources is not resources.phenotypes):
                raise EvidenceStoreError("foreign context ablation resources")
            scores.phenotype_resource.require(writer,scores.universe.members.seed.context,scores.phenotype_score,scores.phenotype_refs)
            self.score = (scores.phenotype_score if resources.prepared is None else
                          resources.prepared.score(scores.universe.members.seed.context.disease_id))
            self.primary = {"name":scores.name,"sha256":scores.sha256}
            self.sha256 = writer.add(name,"context-ablation-evidence-v1",_bounded({
                "schema":"genosapiens.wgs-context-ablation-evidence","version":1,
                "resources":resources.reference,"primary":self.primary,"score":self.score.ledger}))
        except BaseException:
            writer.failed = True
            raise

    @property
    def reference(self):
        return {"name":self.name,"kind":"context-ablation-evidence-v1","sha256":self.sha256}

    def view_arguments(self, writer, scores):
        _primary(writer,scores)
        if self.connection is not writer.connection or self.primary != {"name":scores.name,"sha256":scores.sha256}:
            raise EvidenceStoreError("ablation context score binding changed")
        return dict(mode="dominant_term_sensitivity",phenotype_score=self.score,view_binding=self.reference)


def _replay(reader, reference, kind):
    if reference.kind != kind:
        raise EvidenceStoreError("foreign sensitivity source root")
    reader._require_reference(reference)
    return EvidenceReplay(reader)


class VerifiedCaseDiagnosisEvidence(CaseDiagnosisEvidence):
    def __init__(self, reader, reference, genes, *, diagnosis, prior_policy):
        super().__init__(_replay(reader,reference,"case-diagnosis-evidence-v1"),genes,
                         diagnosis=diagnosis,prior_policy=prior_policy,name=reference.name)


class VerifiedContextDiagnosisEvidence(ContextDiagnosisEvidence):
    def __init__(self, reader, reference, resources, scores):
        super().__init__(_replay(reader,reference,"context-diagnosis-evidence-v1"),resources,scores,name=reference.name)


class VerifiedCaseAblationEvidence(CaseAblationEvidence):
    def __init__(self, reader, reference, phenotypes, primary_order, leader_scores, *, prior_policy):
        super().__init__(_replay(reader,reference,"case-ablation-evidence-v1"),phenotypes,primary_order,leader_scores,
                         prior_policy=prior_policy,name=reference.name)


class VerifiedContextAblationEvidence(ContextAblationEvidence):
    def __init__(self, reader, reference, resources, scores):
        super().__init__(_replay(reader,reference,"context-ablation-evidence-v1"),resources,scores,name=reference.name)
