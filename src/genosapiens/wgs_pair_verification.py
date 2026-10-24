"""Source/member replay and complete ordinal closure for saved pair evidence."""
from array import array
from collections import OrderedDict
from dataclasses import replace
import hashlib
from itertools import zip_longest
import json
import zlib

from genosapiens.evidence_model import HypothesisPhase, _typed
from genosapiens.evidence_serialization import _canonical
from genosapiens.genotype_hypotheses import _context
from genosapiens.genotype_model import assess_genotype_model, model_decision_member, ModelMemberFacts
from genosapiens.model_requirements import assess_model_requirements
from genosapiens.partner_qualification import assess_partner_qualification, qualification_decision_inputs, PartnerQualificationAssessment
from genosapiens import longtr_admission as longtr
from genosapiens.wgs_context_evidence import VerifiedContextEvidence
from genosapiens.wgs_evidence_store import EvidenceReader, EvidenceStoreError
from genosapiens.wgs_genotype_evaluation import prepare_bounded_genotype_evaluation
from genosapiens.genotype_evaluation import PreparedGenotypePolicies
from genosapiens.wgs_pair_decisions import MemberDecisionInputs, VerifiedPairDecisions
from genosapiens.wgs_pair_domain import CompletePairDomain, PairDomainSummary
from genosapiens.wgs_pair_evidence import CompletePairEvidence, MAX_PATTERN_BYTES
from genosapiens.wgs_pair_outcomes import MAX_CATALOG_ENTRIES, VerifiedPairOutcomes
from genosapiens.wgs_phase_evaluation import CachedPairPhases
from genosapiens.wgs_phase_patterns import PairPhasePattern

_PHASES = (HypothesisPhase.UNKNOWN,HypothesisPhase.IN_CIS,HypothesisPhase.IN_TRANS,HypothesisPhase.CONFLICTING)
MAX_PATTERN_CACHE_BYTES = 16 * 1024**2


class VerifiedCompletePairEvidence(CompletePairEvidence):
    def __init__(self, reader, reference, members, facts, *, policy, gene_evidence, prepared_gene_evidence=None, prepared_policies=None):
        if (type(reader) is not EvidenceReader or type(members) is not VerifiedContextEvidence
                or members.writer is not reader or facts.writer is not reader
                or reference.kind != "complete-pair-evidence-v1" or policy.maximum_pairs_per_context is not None):
            raise EvidenceStoreError("foreign complete pair evidence")
        if prepared_policies is not None and (type(prepared_policies) is not PreparedGenotypePolicies
                or policy is not prepared_policies.parent):
            raise EvidenceStoreError('foreign pair policy owner')
        self.writer,self.connection,self.name,self.sha256 = reader,reader.connection,reference.name,reference.sha256
        self.facts,self.policy = facts,policy
        self.prefix = "verified_pair_evidence_"+reference.sha256
        self.tables,self.closed = [],False
        self.pattern_hot,self.outcome_hot,self.qualifying_partners = OrderedDict(),OrderedDict(),set()
        self.domain = CompletePairDomain(scratch_directory=reader.temporary.name)
        try:
            header = reader.reconstruct(reference)
            fields = {"schema","version","context","domain_summary","pair_count","constructible_count",
                "pattern_count","outcome_count","decision_count","genotype_facts","members","decisions",
                "patterns","outcomes","domain","availability","tape_name","tape"}
            if (set(header) != fields or header["schema"] != "genosapiens.wgs-complete-pair-evidence"
                    or type(header["version"]) is not int or header["version"] != 1
                    or header["genotype_facts"] != facts.published
                    or any(type(header[k]) is not int or not 0 <= header[k] <= MAX_CATALOG_ENTRIES
                           for k in ("pattern_count","outcome_count","decision_count"))
                    or any(type(header[k]) is not int or header[k] < 0 for k in ("pair_count","constructible_count"))):
                raise EvidenceStoreError("invalid complete pair evidence header")
            self.context = _context(members.seed.context,"compound_heterozygous")
            from genosapiens.evidence_model import HypothesisContext
            if _typed(header["context"],HypothesisContext,decode=True) != self.context:
                raise EvidenceStoreError("foreign complete pair context")
            self.summary = self.domain.build(m.allele for m in members.iter_members(allele_order=True))
            if (PairDomainSummary(**header["domain_summary"]) != self.summary
                    or header["pair_count"] != self.summary.pair_count or header["constructible_count"] > self.summary.pair_count):
                raise EvidenceStoreError("complete pair domain does not replay")
            self.pair_count,self.hypothesis_count = header["pair_count"],header["constructible_count"]
            self.pattern_count,self.outcome_count,self.decision_count = (header[k] for k in ("pattern_count","outcome_count","decision_count"))
            kinds = {"members":"pair-member-evidence-v1","decisions":"pair-decision-catalog-v1",
                     "patterns":"pair-phase-patterns-v1","outcomes":"pair-computed-outcomes-v1",
                     "domain":"pair-domain-members-v1","availability":"pair-partner-availability-v1"}
            refs = {}
            for key,kind in kinds.items():
                value = header[key]
                if set(value) != {"name","sha256"} or value["name"] != self.name+"."+key:
                    raise EvidenceStoreError("foreign complete pair reference")
                ref = reader.reference(value["name"],kind=kind)
                if ref.sha256 != value["sha256"]:
                    raise EvidenceStoreError("complete pair reference changed")
                refs[key] = ref
                setattr(self,key+"_sha256",ref.sha256)
            for expected,actual in zip_longest(self.domain.members(),reader.iter_array(refs["domain"],max_item_bytes=16*1024**2)):
                if expected is None or actual is None or list(expected) != actual:
                    raise EvidenceStoreError("pair member domain coverage changed")
            self.catalog = VerifiedPairDecisions(reader,refs["decisions"],expected_context=self.context)
            if self.catalog.outcome_count != self.decision_count:
                raise EvidenceStoreError("pair decision catalog count changed")
            for suffix,definition in (
                ("members","id TEXT PRIMARY KEY,class INTEGER,content TEXT"),
                ("patterns","id INTEGER PRIMARY KEY,digest TEXT UNIQUE,content TEXT"),
                ("outcomes","id INTEGER PRIMARY KEY,decision INTEGER,pattern INTEGER,constructible INTEGER,left_nonqualifying INTEGER,right_nonqualifying INTEGER,UNIQUE(decision,pattern)"),
                ("chunks","ordinal INTEGER PRIMARY KEY,count INTEGER,content BLOB"),
            ):
                table = self.prefix+"_"+suffix
                self.connection.execute(f"CREATE TABLE {table} ({definition})")
                self.tables.append(table)
            prepared = prepare_bounded_genotype_evaluation(members=members,facts=facts,context=members.seed.context,
                gene_evidence=gene_evidence,policies={"loci":policy.loci,"requirements":policy.requirements,"clinical":policy.clinical},
                prepared_gene_evidence=prepared_gene_evidence,prepared_policies=prepared_policies)
            classes,used_member_classes = {},set()
            expected_keys = (key for _,key,_,_,parent in self.domain.members() if not parent)
            for key,row in zip_longest(expected_keys,reader.iter_array(refs["members"],max_item_bytes=64*1024**2)):
                if (key is None or row is None or set(row) != {"allele_id","decision_class","model_member","qualification"}
                        or row["allele_id"] != key or type(row["decision_class"]) is not int):
                    raise EvidenceStoreError("pair member evidence coverage changed")
                with prepared.window((key,)):
                    allele = prepared.alleles[key]
                    model = assess_genotype_model(self.context,(allele,),facts.header,loci=policy.loci,prepared=prepared)
                    requirements = assess_model_requirements(self.context,allele,policy=policy.requirements,prepared=prepared)
                    qualified = assess_partner_qualification(self.context,allele,clinical_policy=policy.clinical,
                        requirement_assessments=requirements,genotype_properties=facts.header,disease_gene_evidence=gene_evidence,prepared=prepared)
                    member, = model.members
                    if (_typed(row["model_member"],ModelMemberFacts,decode=True) != member
                            or _typed(row["qualification"],PartnerQualificationAssessment,decode=True) != qualified):
                        raise EvidenceStoreError("pair member evidence does not replay")
                    inputs = MemberDecisionInputs(model_decision_member(member,projected={key} if longtr.is_projection(allele) else set()),
                        qualification_decision_inputs(self.context,qualified.genotype,requirements,None,
                            support=qualified.clinical_support_refs,contradiction=qualified.clinical_contradiction_refs,
                            conflicts=qualified.clinical_conflict_refs))
                    if self.catalog.member_inputs(row["decision_class"]) != inputs:
                        raise EvidenceStoreError("pair member decision assignment is foreign")
                classes[key] = row["decision_class"]
                used_member_classes.add(row["decision_class"])
                self.connection.execute(f"INSERT INTO {self.prefix}_members VALUES (?,?,?)",(key,row["decision_class"],_canonical(row)))
            if len(used_member_classes) != self.catalog.member_count:
                raise EvidenceStoreError("unreferenced pair member decision classes")
            count = 0
            for row in reader.iter_array(refs["patterns"],max_item_bytes=MAX_PATTERN_BYTES+1024):
                if set(row) != {"ordinal","sha256","pattern"} or type(row["ordinal"]) is not int or row["ordinal"] != count:
                    raise EvidenceStoreError("pair phase pattern order changed")
                pattern = _typed(row["pattern"],PairPhasePattern,decode=True)
                content = _canonical(pattern)
                if hashlib.sha256(content.encode("ascii")).hexdigest() != row["sha256"]:
                    raise EvidenceStoreError("pair phase pattern digest changed")
                self.connection.execute(f"INSERT INTO {self.prefix}_patterns VALUES (?,?,?)",(count,row["sha256"],content))
                count += 1
            if count != self.pattern_count:
                raise EvidenceStoreError("pair phase pattern count changed")
            count = 0
            # Four uint32 values per outcome avoid retaining a million Python
            # records: left/right class, categorical phase and pattern ordinal.
            assignments,flags = array("I"),array("B")
            used_decisions,used_patterns = bytearray(self.decision_count),bytearray(self.pattern_count)
            for row in reader.iter_array(refs["outcomes"],max_item_bytes=4096):
                if (set(row) != {"ordinal","decision","pattern"} or type(row["ordinal"]) is not int or row["ordinal"] != count
                        or type(row["decision"]) is not int or not 0 <= row["decision"] < self.decision_count
                        or type(row["pattern"]) is not int or not 0 <= row["pattern"] < self.pattern_count):
                    raise EvidenceStoreError("pair computed outcome coverage changed")
                result = self.catalog.result(row["decision"])
                left,right,phase = self.catalog.outcome_inputs(row["decision"])
                nonqualifying = tuple(s == "nonqualifying" for s in result.qualification_states)
                self.connection.execute(f"INSERT INTO {self.prefix}_outcomes VALUES (?,?,?,?,?,?)",
                    (count,row["decision"],row["pattern"],int(result.constructible),*map(int,nonqualifying)))
                assignments.extend((left,right,_PHASES.index(phase),row["pattern"]))
                flags.extend((int(result.constructible),*map(int,nonqualifying)))
                used_decisions[row["decision"]] = used_patterns[row["pattern"]] = 1
                count += 1
            if count != self.outcome_count or not all(used_decisions) or not all(used_patterns):
                raise EvidenceStoreError("pair computed catalog closure changed")
            if header["tape_name"] != self.name+".tape":
                raise EvidenceStoreError("foreign complete pair tape")
            tape = VerifiedPairOutcomes(reader,reader.reference(header["tape_name"],kind="pair-outcome-tape-v1"),
                domain=self.summary,catalog_sha256=refs["outcomes"].sha256,catalog_count=self.outcome_count)
            if _canonical(tape.tape) != _canonical(header["tape"]):
                raise EvidenceStoreError("pair tape descriptor changed")
            self.tape = tape.tape
            for ordinal,(content,values) in enumerate(tape._chunks(),1):
                self.connection.execute(f"INSERT INTO {self.prefix}_chunks VALUES (?,?,?)",(ordinal,len(values),zlib.compress(content,6)))
            phase_owner = CachedPairPhases(facts)
            pattern_cache = OrderedDict()
            cache_bytes = count = constructed = 0
            for pair,code in zip(self.domain.pairs(),self._outcomes(),strict=True):
                left,right = pair
                offset = code*4
                a,b,phase,pattern_id = assignments[offset:offset+4]
                if (a,b) != (classes[left],classes[right]):
                    raise EvidenceStoreError("pair tape assigned foreign member decision inputs")
                if pattern_id in pattern_cache:
                    pattern,size = pattern_cache[pattern_id]
                    pattern_cache.move_to_end(pattern_id)
                else:
                    content, = self.connection.execute(f"SELECT content FROM {self.prefix}_patterns WHERE id=?",(pattern_id,)).fetchone()
                    pattern,size = _typed(json.loads(content),PairPhasePattern,decode=True),len(content)
                    if size <= MAX_PATTERN_CACHE_BYTES:
                        while pattern_cache and (cache_bytes+size > MAX_PATTERN_CACHE_BYTES or len(pattern_cache) >= 4096):
                            _,(_,old_size) = pattern_cache.popitem(last=False)
                            cache_bytes -= old_size
                        pattern_cache[pattern_id] = pattern,size
                        cache_bytes += size
                if pattern.phase is not _PHASES[phase] or pattern != phase_owner.pattern(left,right):
                    raise EvidenceStoreError("pair phase assignment does not replay")
                offset = code*3
                constructed += flags[offset]
                if not flags[offset+2]:
                    self.qualifying_partners.add(left)
                if not flags[offset+1]:
                    self.qualifying_partners.add(right)
                count += 1
            if (count,constructed) != (self.pair_count,self.hypothesis_count):
                raise EvidenceStoreError("complete pair hypothesis count changed")
            for actual,expected in zip_longest(reader.iter_array(refs["availability"],max_item_bytes=4096),self.partner_availability()):
                if actual != expected:
                    raise EvidenceStoreError("pair partner availability does not replay")
            self.phase_evaluations,self.phase_cache_hits = phase_owner.evaluations,phase_owner.hits
        except BaseException:
            self.close()
            raise
