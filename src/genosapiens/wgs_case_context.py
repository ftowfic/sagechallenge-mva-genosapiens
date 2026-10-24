"""One bounded primary scientific context and its complete dependency closure."""
from contextlib import ExitStack
from genosapiens.wgs_runtime import analysis_progress

from genosapiens.wgs_case_resources import bound_reference, component_reference
from genosapiens.wgs_context_evidence import VerifiedContextEvidence
from genosapiens.wgs_context_order import CompleteContextOrder, VerifiedContextOrder
from genosapiens.wgs_context_population import CompleteContextPopulation, VerifiedContextPopulation
from genosapiens.wgs_context_scores import CompleteContextScores, VerifiedContextScores
from genosapiens.wgs_eligibility_evidence import CompleteEligibilityEvidence, VerifiedEligibilityEvidence
from genosapiens.wgs_evidence_replay import EvidenceReplay
from genosapiens.wgs_evidence_store import EvidenceStoreError
from genosapiens.wgs_genotype_facts import ContextGenotypeFacts, verify_context_genotype_facts
from genosapiens.wgs_hypothesis_evidence import CompleteHypothesisEvidence
from genosapiens.wgs_hypothesis_verification import VerifiedHypothesisEvidence
from genosapiens.wgs_phenotype_evidence import ContextPhenotypeEvidence, VerifiedContextPhenotypeEvidence
from genosapiens.wgs_representation_evidence import CompleteRepresentationEvidence, VerifiedRepresentationEvidence

_KINDS = {
    "members":"context-evidence-v1","facts":"genotype-facts-v1",
    "universe":"complete-hypothesis-universe-v1","population":"context-population-evidence-v1",
    "phenotype":"context-phenotype-evidence-v1","scores":"complete-context-scores-v1",
    "representations":"complete-submission-representations-v1","eligibility":"complete-submission-eligibility-v1",
    "order":"complete-context-order-v1",
}


class ScientificCaseContext:
    def __init__(self, writer, resources, members, reference, reference_pin, renderer_cache, *, name, progress=None,
                 failure_observer=None):
        self.stack,self.closed = ExitStack(),False
        self.failure_observer = failure_observer
        self.writer,self.connection,self.name,self.resources = writer,writer.connection,name,resources
        def report(subphase):
            if failure_observer is not None:
                failure_observer.stage(subphase)
            analysis_progress(progress, 'scientific_evaluation', subphase, scope=members.ordinal)
        try:
            resources.require(writer)
            if reference.binding != resources.reference_binding:
                raise EvidenceStoreError("case submission reference changed")
            if members.closed or members.connection is not self.connection:
                raise EvidenceStoreError("foreign case context members")
            # Ownership is transferred for the duration of this context.
            self.members = self.stack.enter_context(members)
            report('genotype_facts')
            self.facts = self.stack.enter_context(ContextGenotypeFacts(writer,members,
                sample_name=resources.sample_name,source_manifest_sha256=resources.source_manifest_sha256,
                reference_manifest_sha256=resources.reference_binding.manifest_sha256))
            self.facts.publish(name + ".facts")
            report('hypothesis_evidence')
            self.universe = self.stack.enter_context(CompleteHypothesisEvidence(writer,members,self.facts,
                name=name + ".universe",progress=progress,**self._hypothesis_args()))
            report('population')
            self.population = self.stack.enter_context(CompleteContextPopulation(writer,self.universe,
                resources.population,name=name + ".population"))
            report('phenotype')
            self.phenotype = ContextPhenotypeEvidence(writer,resources.phenotypes,members.seed.context,name=name + ".phenotype")
            report('scores')
            self.scores = self.stack.enter_context(CompleteContextScores(writer,self.universe,self.population,
                name=name + ".scores",**self._score_args()))
            report('representation')
            self.representations = self.stack.enter_context(CompleteRepresentationEvidence(writer,members,self.facts,
                reference,renderer_cache,reference_pin,name=name + ".representations"))
            report('eligibility')
            self.eligibility = self.stack.enter_context(CompleteEligibilityEvidence(writer,self.universe,self.population,
                self.representations,name=name + ".eligibility"))
            report('ordering')
            self.order = self.stack.enter_context(CompleteContextOrder(writer,self.scores,name=name + ".order"))
            report('context_publication')
            self._publish(writer)
        except BaseException as error:
            writer.failed = True
            primary = error if failure_observer is None else failure_observer.capture(error)
            from genosapiens.wgs_failures import run_cleanup
            run_cleanup((('context_store', self.close),), primary=primary)
            raise primary from None

    def _hypothesis_args(self):
        r = self.resources
        return dict(policy=r.hypothesis_policy,policy_sha256=r.hypothesis_policy_sha256,
            gene_evidence=r.genes.prepared.evidence,prepared_gene_evidence=r.genes.prepared,
            prepared_policies=r.genotype_policies,
            observation_ids=r.observation_ids,outside_context_partner_uncertainty=r.outside_context_partner_uncertainty)

    def _score_args(self):
        r = self.resources
        return dict(policy=r.fusion_policy,policy_sha256=r.fusion_policy_sha256,
            gene_evidence=r.genes.prepared.evidence,gene_resource=r.genes,
            phenotype_score=self.phenotype.score,phenotype_refs=self.phenotype.profile_refs,phenotype_resource=self.phenotype)

    def _publish(self, writer):
        self.assessment_domain_sha256 = self.resources.assessment_domain_sha256
        self.ordinal,self.context_id = self.members.ordinal,self.members.context_id
        self.hypothesis_count = self.universe.hypothesis_count
        self.pair_count = 0 if self.universe.pairs is None else self.universe.pairs.pair_count
        self.header = {"schema":"genosapiens.wgs-scientific-case-context","version":getattr(self, '_schema_version', 2),
            "assessment_domain_sha256":self.resources.assessment_domain_sha256,
            "resources":self.resources.reference,"ordinal":self.ordinal,"context_id":self.context_id,
            "member_count":self.members.member_count,"pair_count":self.pair_count,"hypothesis_count":self.hypothesis_count}
        for key,kind in _KINDS.items():
            self.header[key] = self.facts.published if key == "facts" else component_reference(getattr(self,key),kind)
        if self.header['version'] == 2:
            from genosapiens.wgs_completed_science import publish_selected_index
            self.header['selected_reconstruction'] = publish_selected_index(writer, self)
        self.sha256 = writer.add(self.name,"scientific-case-context-v1",self.header)

    @property
    def reference(self):
        if self.closed:
            raise EvidenceStoreError("scientific context is closed")
        return component_reference(self,"scientific-case-context-v1")

    def close(self):
        self.stack.close()
        self.closed = True

    def __enter__(self):
        return self

    def __exit__(self,kind,error,traceback):
        if error is None:
            self.close()
            return
        observer = getattr(self, 'failure_observer', None)
        primary = error if observer is None else observer.capture(error)
        from genosapiens.wgs_failures import run_cleanup
        run_cleanup((('context_store', self.close),), primary=primary)
        if primary is not error:
            raise primary from None


class VerifiedScientificCaseContext(ScientificCaseContext):
    def __init__(self, reader, root, resources, reference, reference_pin, renderer_cache, *, retain_identity_index=False):
        self.stack,self.closed = ExitStack(),False
        self.writer,self.connection,self.name,self.resources = reader,reader.connection,root.name,resources
        try:
            resources.require(reader)
            if root.kind != "scientific-case-context-v1" or reference.binding != resources.reference_binding:
                raise EvidenceStoreError("foreign scientific context or submission reference")
            value = reader.reconstruct(root)
            self._schema_version = value.get('version')
            additional = {'selected_reconstruction'} if self._schema_version == 2 else set()
            if (set(value) != {"schema","version","resources","ordinal","context_id","member_count",
                    "pair_count","hypothesis_count","assessment_domain_sha256",*_KINDS,*additional}
                    or value["schema"] != "genosapiens.wgs-scientific-case-context"
                    or type(value["version"]) is not int or value["version"] not in (1, 2)
                    or value["resources"] != resources.reference
                    or value['assessment_domain_sha256']!=resources.assessment_domain_sha256):
                raise EvidenceStoreError("invalid scientific context dependency closure")
            refs = {key:bound_reference(reader,value[key],kind) for key,kind in _KINDS.items()}
            self.members = self.stack.enter_context(VerifiedContextEvidence(reader,refs["members"], resources=resources))
            if self.members.work_filter is not None:
                work = self.members.work_filter_header
                if (work['resources'] != resources.reference
                        or work['assessment_domain_sha256'] != resources.assessment_domain_sha256
                        or work['source_manifest_sha256'] != resources.source_manifest_sha256
                        or work['reference_manifest_sha256'] != resources.reference_binding.manifest_sha256
                        or resources.parent_bindings.get('population_work_policy') != work['policy_sha256']):
                    raise EvidenceStoreError('population work domain resource binding changed')
            self.facts = self.stack.enter_context(verify_context_genotype_facts(reader,self.members,refs["facts"]))
            if (self.facts.header.sample_name != resources.sample_name
                    or self.facts.header.input_sha256 != resources.source_manifest_sha256
                    or self.facts.header.reference_manifest_sha256 != resources.reference_binding.manifest_sha256):
                raise EvidenceStoreError("scientific context genotype parents changed")
            self.universe = self.stack.enter_context(VerifiedHypothesisEvidence(reader,refs["universe"],
                self.members,self.facts,**self._hypothesis_args()))
            self.population = self.stack.enter_context(VerifiedContextPopulation(reader,refs["population"],
                self.universe,resources.population))
            self.phenotype = VerifiedContextPhenotypeEvidence(reader,refs["phenotype"],resources.phenotypes,self.members.seed.context)
            self.scores = self.stack.enter_context(VerifiedContextScores(reader,refs["scores"],self.universe,
                self.population,**self._score_args()))
            self.representations = self.stack.enter_context(VerifiedRepresentationEvidence(reader,refs["representations"],
                self.members,self.facts,reference,renderer_cache,reference_pin))
            self.eligibility = self.stack.enter_context(VerifiedEligibilityEvidence(reader,refs["eligibility"],
                self.universe,self.population,self.representations))
            self.order = self.stack.enter_context(VerifiedContextOrder(reader,refs["order"],self.scores,
                retain_identity_index=retain_identity_index))
            self._publish(EvidenceReplay(reader))
        except BaseException:
            self.close()
            raise
