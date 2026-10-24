"""Complete pair decisions and source-bound phase evidence in one bounded owner.

Every distinct-event pair has an outcome ordinal, including nonconstructible
pairs. Pair/member invariants are evaluated by the existing kernels. This stage
does not publish hypotheses or claim completion of downstream scoring.
"""
from array import array
from collections import OrderedDict
import hashlib
import json
import sys
import zlib

from genosapiens.evidence_model import EnumerationStatus, SourcedRecord, _typed
from genosapiens.evidence_serialization import _canonical
from genosapiens.genotype_hypotheses import _context, _record, construct_evaluated_hypothesis
from genosapiens.genotype_model import (
    assess_genotype_model, model_decision_member, GenotypeModelAssessment, ModelMemberFacts,
    ModelRequirementFacts, POLICY_SHA256 as MODEL_POLICY_SHA256,
)
from genosapiens.model_requirements import assess_model_requirements
from genosapiens.partner_qualification import (
    assess_partner_qualification, qualification_decision_inputs, PartnerQualificationAssessment,
)
from genosapiens import longtr_admission as longtr
from genosapiens.wgs_evidence_store import EvidenceStoreError
from genosapiens.wgs_genotype_evaluation import prepare_bounded_genotype_evaluation
from genosapiens.genotype_evaluation import PreparedGenotypePolicies
from genosapiens.wgs_pair_decisions import PairDecisionCatalog
from genosapiens.wgs_pair_domain import CompletePairDomain
from genosapiens.wgs_pair_outcomes import CHUNK_PAIRS, MAX_CATALOG_ENTRIES, write_pair_outcomes
from genosapiens.wgs_packed_outcomes import outcome_chunk_values
from genosapiens.wgs_phase_patterns import PairPhasePattern, expand_phase_pattern
from genosapiens.wgs_phase_evaluation import CachedPairPhases

MAX_PATTERN_BYTES = 16 * 1024**2
HOT_CLASSES = 4096


from genosapiens.wgs_runtime import analysis_progress


def assemble_evaluated_pair(ordinal, pair, decision, pattern, members, qualifications, *,
        context, header, facts_reference, pair_reference, locus_policy_sha256):
    """Assemble one saved decision and its exact member/phase templates.

    Shared by live computed owners and direct completed-evidence readers. This
    performs no qualification, source query or whole-domain pair evaluation.
    """
    from dataclasses import replace
    from genosapiens.candidate_admission import _hash
    from genosapiens.genotype_hypotheses import admission
    phase = expand_phase_pattern(pattern,allele_ids=pair,
        left_sources=members[0].genotype.event_source_allele_ids,
        right_sources=members[1].genotype.event_source_allele_ids)
    requirements = tuple(ModelRequirementFacts(requirement,
        tuple(r for r in context.evidence_refs if r.record_id in requirement.record_ids),
        result.family,result.required_locus_kind,result.state,result.completeness,result.reasons)
        for requirement,result in zip(context.requirements,decision.requirement_decisions,strict=True))
    binding_fields = {"schema":"genosapiens.wgs-genotype-model-binding","version":1,
        "policy_sha256":MODEL_POLICY_SHA256,"context_sha256":_hash(context),
        "pair_evidence":{**pair_reference,"ordinal":ordinal},
        "genotype_evidence":facts_reference,
        "locus_policy_sha256":locus_policy_sha256,
        "input_sha256":header.input_sha256,"inventory_sha256":header.inventory_sha256,
        "reference_manifest_sha256":header.reference_manifest_sha256,"sample_name":header.sample_name}
    digest = _hash(binding_fields)
    binding = SourcedRecord("genosapiens.wgs-genotype-model-binding","1","model-binding-"+digest,digest,None,binding_fields)
    model = GenotypeModelAssessment(context.context_id,context.model,decision.model_state,
        decision.completeness,decision.model_reasons,tuple(members),requirements,phase.phase,phase,
        decision.read_validation_required,any(m.genotype.screened_low_balance_candidate for m in members),binding)
    qualified = tuple(replace(q,state=state,reasons=reasons,phase=phase,
        binding={**q.binding,"phase_sha256":_hash(phase),"whole_model_sha256":admission._hash(model)})
        for q,state,reasons in zip(qualifications,decision.qualification_states,decision.qualification_reasons,strict=True))
    return ordinal,pair,decision,model,qualified


class CompletePairEvidence:
    def __init__(self, writer, members, facts, *, policy, gene_evidence, name, prepared_gene_evidence=None, prepared_policies=None, progress=None):
        if prepared_policies is not None and (type(prepared_policies) is not PreparedGenotypePolicies
                or policy is not prepared_policies.parent):
            raise EvidenceStoreError('foreign pair policy owner')
        self.writer,self.connection = writer,writer.connection
        self.name = name
        self.facts,self.policy = facts,policy
        self.prefix = "pair_evidence_" + hashlib.sha256(name.encode()).hexdigest()
        self.tables,self.closed = [],False
        self.domain = CompletePairDomain(scratch_directory=writer.temporary.name)
        self.pattern_hot,self.outcome_hot = OrderedDict(),OrderedDict()
        self.pattern_count = self.outcome_count = self.pair_count = self.hypothesis_count = 0
        self.qualifying_partners = set()
        try:
            if members.writer is not writer or facts.writer is not writer:
                raise EvidenceStoreError("foreign complete pair owners")
            if policy.maximum_pairs_per_context is not None:
                raise EvidenceStoreError("complete pair evidence requires an untruncated pair policy")
            self.context = _context(members.seed.context,"compound_heterozygous")
            self.summary = self.domain.build(m.allele for m in members.iter_members(allele_order=True))
            for suffix,definition in (
                ("members","id TEXT PRIMARY KEY,class INTEGER,content TEXT"),
                ("patterns","id INTEGER PRIMARY KEY,digest TEXT UNIQUE,content TEXT"),
                ("outcomes","id INTEGER PRIMARY KEY,decision INTEGER,pattern INTEGER,constructible INTEGER,left_nonqualifying INTEGER,right_nonqualifying INTEGER,UNIQUE(decision,pattern)"),
                ("chunks","ordinal INTEGER PRIMARY KEY,count INTEGER,content BLOB"),
            ):
                table = self.prefix+"_"+suffix
                self.connection.execute(f"CREATE TABLE {table} ({definition})")
                self.tables.append(table)
            self.catalog = PairDecisionCatalog(writer,self.context)
            if facts.published is None:
                facts.publish(name+".genotypes")
            prepared = prepare_bounded_genotype_evaluation(members=members,facts=facts,
                context=members.seed.context,gene_evidence=gene_evidence,
                policies={"loci":policy.loci,"requirements":policy.requirements,"clinical":policy.clinical},
                prepared_gene_evidence=prepared_gene_evidence,prepared_policies=prepared_policies)
            classes = {}
            def member_rows():
                for _,key,_,_,parent in self.domain.members():
                    if parent:
                        continue
                    with prepared.window((key,)):
                        allele = prepared.alleles[key]
                        model = assess_genotype_model(self.context,(allele,),facts.header,loci=policy.loci,prepared=prepared)
                        requirements = assess_model_requirements(self.context,allele,policy=policy.requirements,prepared=prepared)
                        qualification = assess_partner_qualification(self.context,allele,clinical_policy=policy.clinical,
                            requirement_assessments=requirements,genotype_properties=facts.header,
                            disease_gene_evidence=gene_evidence,prepared=prepared)
                        member, = model.members
                        inputs = model_decision_member(member,projected={key} if longtr.is_projection(allele) else set())
                        q_inputs = qualification_decision_inputs(self.context,qualification.genotype,
                            qualification.requirement_assessments,None,support=qualification.clinical_support_refs,
                            contradiction=qualification.clinical_contradiction_refs,conflicts=qualification.clinical_conflict_refs)
                        cls = self.catalog.member(inputs,q_inputs)
                        classes[key] = cls
                        value = {"allele_id":key,"decision_class":cls.ordinal,"model_member":member,
                                 "qualification":qualification}
                        content = _canonical(value)
                        if len(content) > 64*1024**2:
                            raise EvidenceStoreError("pair member evidence exceeds capacity",capacity=True)
                        self.connection.execute(f"INSERT INTO {self.prefix}_members VALUES (?,?,?)",(key,cls.ordinal,content))
                        yield value
            self.members_sha256 = writer.add_sequence(name+".members","pair-member-evidence-v1",member_rows())
            phase_owner = CachedPairPhases(facts)
            pattern_ids = OrderedDict()
            packed = array("I")
            chunks = 0
            def flush():
                nonlocal packed,chunks
                if sys.byteorder != "little":
                    packed.byteswap()
                chunks += 1
                self.connection.execute(f"INSERT INTO {self.prefix}_chunks VALUES (?,?,?)",
                    (chunks,len(packed),zlib.compress(packed.tobytes(),6)))
                packed = array("I")
            def report_pairs():
                analysis_progress(progress, 'scientific_evaluation', 'logical_pairs',
                    counts={'logical_pairs_completed':self.pair_count, 'logical_pairs_total':self.summary.pair_count},
                    unit='pairs', completed=self.pair_count, total=self.summary.pair_count, scope=members.ordinal)
            report_pairs()
            for left,right in self.domain.pairs():
                pattern = phase_owner.pattern(left,right)
                if pattern in pattern_ids:
                    pattern_id = pattern_ids[pattern]
                else:
                    pattern_id = self._pattern(pattern)
                    if len(_canonical(pattern)) <= 4096:
                        self._remember(pattern_ids,pattern,pattern_id)
                decision_id = self.catalog.evaluate(classes[left],classes[right],pattern.phase)
                outcome_id,constructible,left_nonqualifying,right_nonqualifying = self._outcome(decision_id,pattern_id)
                if not right_nonqualifying:
                    self.qualifying_partners.add(left)
                if not left_nonqualifying:
                    self.qualifying_partners.add(right)
                packed.append(outcome_id)
                self.pair_count += 1
                self.hypothesis_count += constructible
                if self.pair_count % 256 == 0:
                    report_pairs()
                if len(packed) == CHUNK_PAIRS:
                    flush()
            report_pairs()
            if packed:
                flush()
            if self.pair_count != self.summary.pair_count:
                raise EvidenceStoreError("complete pair evidence coverage changed")
            self.decision_count = self.catalog.outcomes
            self.decision_evaluations = self.catalog.evaluations
            self.phase_evaluations,self.phase_cache_hits = phase_owner.evaluations,phase_owner.hits
            self.decisions_sha256 = self.catalog.publish(name+".decisions",retain_projection=True)
            self.patterns_sha256 = writer.add_sequence(name+".patterns","pair-phase-patterns-v1",(
                {"ordinal":i,"sha256":digest,"pattern":json.loads(content)} for i,digest,content in
                self.connection.execute(f"SELECT * FROM {self.prefix}_patterns ORDER BY id")))
            self.outcomes_sha256 = writer.add_sequence(name+".outcomes","pair-computed-outcomes-v1",(
                {"ordinal":i,"decision":decision,"pattern":pattern} for i,decision,pattern in
                self.connection.execute(f"SELECT id,decision,pattern FROM {self.prefix}_outcomes ORDER BY id")))
            self.tape = write_pair_outcomes(writer,name+".tape",self._outcomes(),domain=self.summary,
                catalog_sha256=self.outcomes_sha256,catalog_count=self.outcome_count)
            self.domain_sha256 = writer.add_sequence(name+".domain","pair-domain-members-v1",self.domain.members())
            self.availability_sha256 = writer.add_sequence(name+".availability","pair-partner-availability-v1",self.partner_availability())
            self.sha256 = writer.add(name,"complete-pair-evidence-v1",{
                "schema":"genosapiens.wgs-complete-pair-evidence","version":1,"context":self.context,
                "domain_summary":self.summary,"pair_count":self.pair_count,"constructible_count":self.hypothesis_count,
                "pattern_count":self.pattern_count,"outcome_count":self.outcome_count,"decision_count":self.decision_count,
                "genotype_facts":facts.published,
                **{key:{"name":name+"."+key,"sha256":getattr(self,key+"_sha256")}
                   for key in ("members","decisions","patterns","outcomes","domain","availability")},
                "tape_name":name+".tape","tape":self.tape})
        except BaseException:
            writer.failed = True
            self.close()
            raise

    @staticmethod
    def _remember(cache,key,value):
        cache[key] = value
        cache.move_to_end(key)
        if len(cache) > HOT_CLASSES:
            cache.popitem(last=False)

    def _pattern(self, pattern):
        content = _canonical(pattern)
        if len(content) > MAX_PATTERN_BYTES:
            raise EvidenceStoreError("pair phase pattern exceeds capacity",capacity=True)
        digest = hashlib.sha256(content.encode("ascii")).hexdigest()
        if digest in self.pattern_hot:
            return self.pattern_hot[digest]
        row = self.connection.execute(f"SELECT id,content FROM {self.prefix}_patterns WHERE digest=?",(digest,)).fetchone()
        if row is None:
            if self.pattern_count == MAX_CATALOG_ENTRIES:
                raise EvidenceStoreError("pair phase catalog exceeds capacity",capacity=True)
            ordinal = self.pattern_count
            self.connection.execute(f"INSERT INTO {self.prefix}_patterns VALUES (?,?,?)",(ordinal,digest,content))
            self.pattern_count += 1
        else:
            ordinal,prior = row
            if content != prior:
                raise EvidenceStoreError("pair phase pattern digest collision")
        self._remember(self.pattern_hot,digest,ordinal)
        return ordinal

    def _outcome(self, decision, pattern):
        key = (decision,pattern)
        if key in self.outcome_hot:
            return self.outcome_hot[key]
        row = self.connection.execute(f"SELECT id,constructible,left_nonqualifying,right_nonqualifying FROM {self.prefix}_outcomes WHERE decision=? AND pattern=?",key).fetchone()
        if row is None:
            if self.outcome_count == MAX_CATALOG_ENTRIES:
                raise EvidenceStoreError("pair evidence outcome capacity exceeded",capacity=True)
            ordinal = self.outcome_count
            result = self.catalog.result(decision)
            constructible = result.constructible
            left_nonqualifying,right_nonqualifying = (s == "nonqualifying" for s in result.qualification_states)
            self.connection.execute(f"INSERT INTO {self.prefix}_outcomes VALUES (?,?,?,?,?,?)",
                (ordinal,*key,int(constructible),int(left_nonqualifying),int(right_nonqualifying)))
            self.outcome_count += 1
        else:
            ordinal,constructible,left_nonqualifying,right_nonqualifying = row
        result = (ordinal,bool(constructible),bool(left_nonqualifying),bool(right_nonqualifying))
        self._remember(self.outcome_hot,key,result)
        return result

    def _outcomes(self):
        for count,content in self.connection.execute(f"SELECT count,content FROM {self.prefix}_chunks ORDER BY ordinal"):
            decoder = zlib.decompressobj()
            value = decoder.decompress(content,count*4+1)
            if len(value) != count*4 or not decoder.eof or decoder.unused_data or decoder.unconsumed_tail:
                raise EvidenceStoreError("private pair outcome chunk changed")
            values = array("I")
            values.frombytes(value)
            if sys.byteorder != "little":
                values.byteswap()
            yield from values

    def partner_availability(self):
        if self.closed:
            raise EvidenceStoreError("complete pair evidence owner is closed")
        sql = ("SELECT m.allele_id,m.parent,coalesce(g.n,0) FROM members m LEFT JOIN "
               "(SELECT event,count(*) AS n FROM members WHERE parent=0 GROUP BY event) g "
               "ON m.event=g.event ORDER BY m.ordinal")
        for key,parent,event_count in self.domain.connection.execute(sql):
            yield {"allele_id":key,"evaluated_partners":0 if parent else self.summary.pairable_count-event_count,
                "distinct_event_partners":self.summary.pairable_count-event_count,
                "has_not_nonqualifying_partner":key in self.qualifying_partners}

    def member_templates(self):
        """Bounded actual source templates, in complete pair-member order."""
        if self.closed:
            raise EvidenceStoreError("complete pair evidence owner is closed")
        for key, content in self.connection.execute(f"SELECT id,content FROM {self.prefix}_members ORDER BY id"):
            value = json.loads(content)
            yield key, _typed(value["model_member"], ModelMemberFacts, decode=True), \
                _typed(value["qualification"], PartnerQualificationAssessment, decode=True)

    def evaluated_pairs(self):
        """Expand already computed outcomes; no model/qualification/phase replay.

        New parent bindings identify normalized D28 owners explicitly. Full
        member facts, requirements, phase bases and clinical/source proofs are
        reconstructed exactly; obsolete expanded parent hashes are not claimed.
        """
        if self.closed:
            raise EvidenceStoreError("complete pair evidence owner is closed")
        for ordinal,(pair,code) in enumerate(zip(self.domain.pairs(),self._outcomes(),strict=True),1):
            yield self._evaluated_pair(ordinal,pair,code)

    def evaluated_at(self, left, right, *, positions):
        if self.closed or positions.domain is not self.domain:
            raise EvidenceStoreError("foreign pair position owner")
        ordinal = positions.ordinal(left,right)
        index = ordinal-1
        row = self.connection.execute(f"SELECT count,content FROM {self.prefix}_chunks WHERE ordinal=?",
                                      (index//CHUNK_PAIRS+1,)).fetchone()
        if row is None:
            raise EvidenceStoreError("pair outcome position is unavailable")
        code = outcome_chunk_values(*row)[index%CHUNK_PAIRS]
        return self._evaluated_pair(ordinal,(positions.keys[left],positions.keys[right]),code)

    def _evaluated_pair(self, ordinal, pair, code):
        from genosapiens.candidate_admission import _hash
        decision_id,pattern_id = self.connection.execute(
            f"SELECT decision,pattern FROM {self.prefix}_outcomes WHERE id=?",(code,)).fetchone()
        decision = self.catalog.result(decision_id)
        content, = self.connection.execute(f"SELECT content FROM {self.prefix}_patterns WHERE id=?",(pattern_id,)).fetchone()
        pattern = _typed(json.loads(content),PairPhasePattern,decode=True)
        members,qualifications = [],[]
        for key in pair:
            content, = self.connection.execute(f"SELECT content FROM {self.prefix}_members WHERE id=?",(key,)).fetchone()
            row = json.loads(content)
            members.append(_typed(row["model_member"],ModelMemberFacts,decode=True))
            qualifications.append(_typed(row["qualification"],PartnerQualificationAssessment,decode=True))
        return assemble_evaluated_pair(ordinal, pair, decision, pattern, members, qualifications,
            context=self.context, header=self.facts.header, facts_reference=self.facts.published,
            pair_reference={"name":self.name,"sha256":self.sha256},
            locus_policy_sha256=_hash(tuple(sorted(self.policy.loci,key=_canonical))))

    def constructed_pairs(self, *, observation_ids, outside_context_partner_uncertainty=True):
        count = 0
        for evaluated in self.evaluated_pairs():
            value = self.construct_evaluated(evaluated,observation_ids=observation_ids,
                outside_context_partner_uncertainty=outside_context_partner_uncertainty)
            if value is None:
                continue
            count += 1
            yield value
        if count != self.hypothesis_count:
            raise EvidenceStoreError("constructed pair coverage changed")

    def construct_evaluated(self, evaluated, *, observation_ids, outside_context_partner_uncertainty):
        _,_,decision,model,qualified = evaluated
        if not decision.constructible:
            return None
        enumeration_evidence = (_record("partner-search",{"context_id":self.facts.context_id,
            "search_evidence":{"name":self.name,"sha256":self.sha256}}),)
        value = construct_evaluated_hypothesis(context=self.context,admission_context_id=self.facts.context_id,
            references=tuple(m.allele for m in model.members),properties=tuple(m.genotype for m in model.members),
            assessments=qualified,phase=model.phase,phase_evidence=model.phase_evidence,model=model,
            enumeration=EnumerationStatus.COMPLETE,enumeration_evidence=enumeration_evidence,
            sample_name=self.facts.header.sample_name,observation_ids=observation_ids,
            incompleteness_penalty=self.policy.incompleteness_penalty,
            outside_context_partner_uncertainty=outside_context_partner_uncertainty)
        if value is None:
            raise EvidenceStoreError("computed pair construction changed")
        return value

    def close(self):
        self.domain.close()
        if hasattr(self,"catalog"):
            self.catalog.close()
        for table in reversed(self.tables):
            self.connection.execute("DROP TABLE "+table)
        self.tables.clear()
        self.pattern_hot.clear()
        self.outcome_hot.clear()
        self.qualifying_partners.clear()
        self.closed = True

    def __enter__(self):
        return self

    def __exit__(self,*_):
        self.close()
