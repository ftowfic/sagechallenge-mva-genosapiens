"""Complete singleton and pair hypothesis universe over bounded context owners.

Singletons are computed records. Pairs are the complete computed outcome tape
and member templates, with explicit logical reconstruction. This module does
not retain every expanded pair or hypothesis and does not score or prune them.
"""
from dataclasses import replace
import hashlib
import json

from genosapiens import longtr_admission as longtr
from genosapiens.allele_identity import RepeatAllele
from genosapiens.candidate_admission import _hash, _native
from genosapiens.evidence_model import EnumerationStatus, HypothesisPhase, _sha, _typed
from genosapiens.evidence_serialization import _canonical
from genosapiens.genotype_hypotheses import (
    ConstructedHypothesis, EVALUABLE_DOMAIN_RULES, _context, _requires_pairs, _record,
    construct_evaluated_hypothesis as _construct_evaluated_hypothesis, qualify_whole_model, consuming_model_contexts,
)
from genosapiens.genotype_model import assess_genotype_model, GenotypeModelAssessment
from genosapiens.model_requirements import assess_model_requirements
from genosapiens.partner_qualification import assess_partner_qualification, PartnerQualificationAssessment
from genosapiens.partner_relevance import SourceContextRelevance
from genosapiens.wgs_evidence_store import EvidenceStoreError
from genosapiens.wgs_genotype_evaluation import prepare_bounded_genotype_evaluation
from genosapiens.genotype_evaluation import PreparedGenotypePolicies
from genosapiens.wgs_pair_evidence import CompletePairEvidence
from genosapiens.genotype_evaluation import prepare_gene_evidence
from genosapiens.wgs_analysis_eligibility import METHOD_SHA256
from genosapiens.wgs_validation_scope import member_validation_scope


@member_validation_scope()
def construct_evaluated_hypothesis(*args, **kwargs):
    return _construct_evaluated_hypothesis(*args, **kwargs)

MAX_RECORD_BYTES = 64 * 1024**2
MAX_SOURCE_METADATA_BYTES = 64 * 1024**2
SINGLETON_CACHE_BYTES = 8 * 1024**2
SINGLETON_CACHE_ENTRIES = 128


class CompleteHypothesisEvidence:
    def __init__(self, writer, members, facts, *, policy, gene_evidence, policy_sha256,
                 observation_ids, outside_context_partner_uncertainty=True, name, prepared_gene_evidence=None, prepared_policies=None, progress=None):
        self.progress = progress
        self.prepared_policies = prepared_policies
        if prepared_policies is not None and (type(prepared_policies) is not PreparedGenotypePolicies
                or policy is not prepared_policies.parent):
            raise EvidenceStoreError('foreign hypothesis policy owner')
        self.writer,self.connection,self.members,self.facts = writer,writer.connection,members,facts
        self.policy,self.observation_ids = policy,observation_ids
        self.outside_context_partner_uncertainty = outside_context_partner_uncertainty
        self.name,self.prefix = name,"hypotheses_"+hashlib.sha256(name.encode()).hexdigest()
        self.tables,self.closed,self.pairs = [],False,None
        self._singleton_cache_owner = (writer, self.connection)
        self._singleton_cache = {}
        self._singleton_cache_bytes = 0
        self.singleton_count = self.hypothesis_count = 0
        _sha(policy_sha256)
        if type(outside_context_partner_uncertainty) is not bool:
            raise EvidenceStoreError("invalid outside-context uncertainty")
        if policy.maximum_pairs_per_context is not None:
            raise EvidenceStoreError("complete hypothesis evidence requires an untruncated pair policy")
        try:
            self.prepared_gene_evidence = prepared_gene_evidence or prepare_gene_evidence(gene_evidence)
            if prepared_gene_evidence is not None:
                prepared_gene_evidence.require(gene_evidence)
            gene_evidence = self.prepared_gene_evidence.evidence
            for suffix,definition in (
                ("singletons","id TEXT PRIMARY KEY,content TEXT,locus INTEGER"),
                ("constructed","id TEXT PRIMARY KEY,content TEXT"),
                ("availability","id TEXT PRIMARY KEY,evaluated INTEGER,distinct_count INTEGER,available INTEGER"),
                ("unresolved","id TEXT PRIMARY KEY"),
            ):
                table = self.prefix+"_"+suffix
                self.connection.execute(f"CREATE TABLE {table} ({definition})")
                self.tables.append(table)
            if facts.published is None:
                facts.publish(name+".genotypes")
            if _requires_pairs(members.seed.context):
                self.pairs = self._prepare_pairs(writer,members,facts,policy=policy,gene_evidence=gene_evidence,name=name+".pairs")
                self.connection.executemany(f"INSERT INTO {self.prefix}_availability VALUES (?,?,?,?)",
                    ((r["allele_id"],r["evaluated_partners"],r["distinct_event_partners"],int(r["has_not_nonqualifying_partner"]))
                     for r in self.pairs.partner_availability()))
            prepared = prepare_bounded_genotype_evaluation(members=members,facts=facts,context=members.seed.context,
                gene_evidence=gene_evidence,policies={"loci":policy.loci,"requirements":policy.requirements,"clinical":policy.clinical},
                prepared_gene_evidence=self.prepared_gene_evidence,prepared_policies=self.prepared_policies)
            unresolved_bytes = 0
            def qualify(ctx,allele):
                requirements = assess_model_requirements(ctx,allele,policy=policy.requirements,prepared=prepared)
                return assess_partner_qualification(ctx,allele,clinical_policy=policy.clinical,
                    requirement_assessments=requirements,genotype_properties=facts.header,disease_gene_evidence=gene_evidence,prepared=prepared)
            def member_rows():
                nonlocal unresolved_bytes
                for member in members.iter_members(allele_order=True):
                    key = member.allele.allele_id
                    with prepared.window((key,)):
                        allele = prepared.alleles[key]
                        prop = facts.allele(key)
                        if prop.alternate_present is None and not (longtr.is_parent(allele)
                                and longtr.metadata(allele).get("source_alt_scope_complete") is True):
                            for source in prop.source_allele_ids:
                                inserted = self.connection.execute(f"INSERT OR IGNORE INTO {self.prefix}_unresolved VALUES (?)",(source,)).rowcount
                                unresolved_bytes += inserted*len(source)
                                if unresolved_bytes > MAX_SOURCE_METADATA_BYTES:
                                    raise EvidenceStoreError("unresolved source metadata exceeds capacity",capacity=True)
                        base = qualify(members.seed.context,allele)
                        ctx = consuming_model_contexts(members.seed.context, allele, prop)[0]
                        model = assess_genotype_model(ctx,(allele,),facts.header,loci=policy.loci,prepared=prepared)
                        own = qualify(ctx,allele)
                        state,reasons = qualify_whole_model(own.state,own.reasons,model.state)
                        own = replace(own,state=state,reasons=reasons,binding={**own.binding,"whole_model_sha256":_hash(model)})
                        native = _native(allele)
                        interpretation = (native.source.source_fields or {}).get("native_interpretation") or {}
                        locus = isinstance(native.identity,RepeatAllele) and interpretation.get("repeat_representation") == "locus_call"
                        content = _canonical({"context":ctx,"model":model,"qualification":own})
                        if len(content) > MAX_RECORD_BYTES:
                            raise EvidenceStoreError("singleton hypothesis evidence exceeds capacity",capacity=True)
                        self.connection.execute(f"INSERT INTO {self.prefix}_singletons VALUES (?,?,?)",(key,content,int(locus)))
                        yield base
            self.member_qualifications_sha256 = writer.add_sequence(name+".member_qualifications","search-member-qualifications-v1",member_rows())
            self.unresolved_sha256 = writer.add_sequence(name+".unresolved","search-unresolved-sources-v1",
                (r[0] for r in self.connection.execute(f"SELECT id FROM {self.prefix}_unresolved ORDER BY id")))
            unresolved_count, = self.connection.execute(f"SELECT count(*) FROM {self.prefix}_unresolved").fetchone()
            relevance_binding = {"schema":"genosapiens.wgs-context-relevance","version":1,
                "context_sha256":_hash(members.seed.context),"seed_sha256":_hash(members.seed),
                "scope":"globally-closed-retained-context-membership"}
            def relevance_rows():
                for source in facts.source_rows():
                    yield SourceContextRelevance(source.source_record_id,source.source_allele_id,
                        source.derived_allele_ids,"derived",(),"relevant",(),(facts.source_proof(source.source_allele_id),),
                        members.seed.context.evidence_refs,(),(),(),
                        ("globally-closed-admitted-context-member","source-target-association-replayed-upstream"),_hash(relevance_binding))
            self.relevance_sha256 = writer.add_sequence(name+".relevance","search-source-relevance-v1",relevance_rows())
            limits = ["outside-supplied-scope-partners-unknown","qualification-not-population-eligibility",
                      "enumeration-is-source-evaluable-domain-not-exhaustive-discovery"]
            if outside_context_partner_uncertainty:
                limits.append("excluded-or-unresolved-source-partners-not-absence-evidence")
            self.search_sha256 = writer.add(name+".search","complete-partner-search-v1",{
                "schema":"genosapiens.wgs-complete-partner-search","version":1,
                "assessment_method_sha256":METHOD_SHA256,"evaluable_domain_rules":EVALUABLE_DOMAIN_RULES,
                "context":members.seed.context,"member_count":members.member_count,"policy_sha256":policy_sha256,
                "context_evidence":{"name":members.name,"sha256":members.sha256},"genotype_facts":facts.published,
                "pairs":None if self.pairs is None else {"name":self.pairs.name,"sha256":self.pairs.sha256},
                "total_distinct_pairs":0 if self.pairs is None else self.pairs.pair_count,
                "evaluated_pairs":0 if self.pairs is None else self.pairs.pair_count,
                "traversal":"complete","enumeration_status":"complete","source_inventory_status":"assessed-domain-complete",
                "limits":tuple(limits),"outside_context_partner_uncertainty":outside_context_partner_uncertainty,
                "unresolved_relevance_source_allele_ids":(),"relevance_binding":relevance_binding,
                **{key:{"name":name+"."+key,"sha256":getattr(self,key+"_sha256")}
                   for key in ("member_qualifications","relevance","unresolved")}})
            enumeration_evidence = (_record("partner-search",{"context_id":members.context_id,
                "search_evidence":{"name":name+".search","sha256":self.search_sha256}}),)
            for key,content,locus in self.connection.execute(f"SELECT * FROM {self.prefix}_singletons ORDER BY id"):
                row = json.loads(content)
                model = _typed(row["model"],GenotypeModelAssessment,decode=True)
                own = _typed(row["qualification"],PartnerQualificationAssessment,decode=True)
                if own.state == "nonqualifying":
                    continue
                available = self.connection.execute(f"SELECT evaluated,distinct_count,available FROM {self.prefix}_availability WHERE id=?",(key,)).fetchone()
                no_partner = (not unresolved_count and not outside_context_partner_uncertainty and own.state == "qualifying"
                    and available is not None and available[0] == available[1] and not available[2])
                from genosapiens.evidence_model import HypothesisContext
                ctx = _typed(row["context"],HypothesisContext,decode=True)
                value = construct_evaluated_hypothesis(context=ctx,admission_context_id=members.context_id,
                    references=(model.members[0].allele,),properties=(model.members[0].genotype,),assessments=(own,),
                    phase=HypothesisPhase.NOT_APPLICABLE,phase_evidence=None,model=model,enumeration=EnumerationStatus.COMPLETE,
                    enumeration_evidence=enumeration_evidence,sample_name=facts.header.sample_name,observation_ids=observation_ids,
                    incompleteness_penalty=policy.incompleteness_penalty,outside_context_partner_uncertainty=outside_context_partner_uncertainty,
                    conditional=ctx.model == "incomplete" and no_partner and not model.members[0].genotype.read_validation_required,locus_call=bool(locus))
                if value is not None:
                    content = _canonical(value)
                    if len(content) > MAX_RECORD_BYTES:
                        raise EvidenceStoreError("constructed singleton exceeds capacity",capacity=True)
                    self.connection.execute(f"INSERT INTO {self.prefix}_constructed VALUES (?,?)",(value.hypothesis.hypothesis_id,content))
                    self.singleton_count += 1
            self.singletons_sha256 = writer.add_sequence(name+".singletons","constructed-singletons-v1",self.singletons())
            self.hypothesis_count = self.singleton_count + (0 if self.pairs is None else self.pairs.hypothesis_count)
            self.sha256 = writer.add(name,"complete-hypothesis-universe-v1",{
                "schema":"genosapiens.wgs-complete-hypothesis-universe","version":1,"context_id":members.context_id,
                "assessment_method_sha256":METHOD_SHA256,"evaluable_domain_rules":EVALUABLE_DOMAIN_RULES,
                "search":{"name":name+".search","sha256":self.search_sha256},
                "singletons":{"name":name+".singletons","sha256":self.singletons_sha256},
                "pair_evidence":None if self.pairs is None else {"name":self.pairs.name,"sha256":self.pairs.sha256},
                "singleton_count":self.singleton_count,"pair_hypothesis_count":self.hypothesis_count-self.singleton_count,
                "hypothesis_count":self.hypothesis_count,"sample_name":facts.header.sample_name,
                "observation_ids":observation_ids,"incompleteness_penalty":policy.incompleteness_penalty})
        except BaseException:
            writer.failed = True
            self.close()
            raise

    def _prepare_pairs(self, writer, members, facts, *, policy, gene_evidence, name):
        return CompletePairEvidence(writer,members,facts,policy=policy,gene_evidence=gene_evidence,name=name,
                                    prepared_gene_evidence=self.prepared_gene_evidence, prepared_policies=self.prepared_policies, progress=self.progress)

    def _decode_singleton(self, content):
        if (self.closed or self._singleton_cache_owner[0] is not self.writer
                or self._singleton_cache_owner[1] is not self.connection):
            raise EvidenceStoreError("singleton decoding owner changed or closed")
        if type(content) is not str or len(content) > MAX_RECORD_BYTES:
            raise EvidenceStoreError("invalid stored singleton record", capacity=True)
        found = self._singleton_cache.get(content)
        if found is not None:
            return found
        value = _typed(json.loads(content), ConstructedHypothesis, decode=True)
        remaining = SINGLETON_CACHE_BYTES - self._singleton_cache_bytes
        if len(self._singleton_cache) < SINGLETON_CACHE_ENTRIES and len(content) + 2048 <= remaining:
            from genosapiens.wgs_context_evidence import _member_cache_size
            size = _member_cache_size(content, value, limit=remaining)
            if size <= remaining:
                self._singleton_cache[content] = value
                self._singleton_cache_bytes += size
        return value

    def singletons(self):
        if self.closed:
            raise EvidenceStoreError("complete hypothesis evidence owner is closed")
        for (content,) in self.connection.execute(f"SELECT content FROM {self.prefix}_constructed ORDER BY id"):
            yield self._decode_singleton(content)

    def constructed(self):
        """Logical export; downstream production can consume computed families."""
        yield from self.singletons()
        if self.pairs is not None:
            # All hypotheses bind the full search, including singleton/source
            # uncertainty accounting, rather than the pair-only subrecord.
            for value in self.pairs.constructed_pairs(observation_ids=self.observation_ids,
                    outside_context_partner_uncertainty=self.outside_context_partner_uncertainty):
                yield self._bind_pair_search(value)

    def _bind_pair_search(self, value):
        evidence = (_record("partner-search",{"context_id":self.members.context_id,
            "search_evidence":{"name":self.name+".search","sha256":self.search_sha256}}),)
        return replace(value,hypothesis=replace(value.hypothesis,enumeration_evidence=evidence))

    def constructed_at(self, locator, *, positions=None):
        """Reconstruct one evaluated occurrence without scanning earlier pairs."""
        if (self.closed or type(locator) is not tuple or len(locator) != 2
                or any(type(v) is not int for v in locator)):
            raise EvidenceStoreError("invalid hypothesis locator")
        left,right = locator
        if right == 0:
            if not 1 <= left <= self.singleton_count:
                raise EvidenceStoreError("singleton position is unavailable")
            content, = self.connection.execute(f"SELECT content FROM {self.prefix}_constructed ORDER BY id LIMIT 1 OFFSET ?",
                                               (left-1,)).fetchone()
            return self._decode_singleton(content)
        if self.pairs is None or positions is None:
            raise EvidenceStoreError("pair position owner is required")
        value = self.pairs.construct_evaluated(self.pairs.evaluated_at(left,right,positions=positions),
            observation_ids=self.observation_ids,outside_context_partner_uncertainty=self.outside_context_partner_uncertainty)
        if value is None:
            raise EvidenceStoreError("pair position has no constructed hypothesis")
        return self._bind_pair_search(value)

    def close(self):
        try:
            if self.pairs is not None:
                self.pairs.close()
            for table in reversed(self.tables):
                self.connection.execute("DROP TABLE "+table)
            self.tables.clear()
            self.closed = True
        finally:
            self._singleton_cache.clear()
            self._singleton_cache_bytes = 0


    def __enter__(self):
        return self

    def __exit__(self,*_):
        self.close()
