"""Complete submission eligibility assignments over the closed hypothesis domain.

Every singleton and every pair occurrence is assigned. The finite catalog
retains explicit inputs to the original ordered-reason kernel; pair model,
read-validation, population and representation facts remain independent inputs.
"""
from array import array
from collections import OrderedDict
import hashlib
import json
import sys

from genosapiens.evidence_model import Completeness, EnumerationStatus, _typed
from genosapiens.evidence_serialization import _canonical
from genosapiens.genotype_hypotheses import EVALUABLE_DOMAIN_RULES
from genosapiens.wgs_analysis_eligibility import METHOD_SHA256
from genosapiens.population_policy import ContextKey, PopulationOutcome
from genosapiens.track1_submission import (
    SubmissionEligibilityInputs, evaluate_submission_eligibility, submission_representation_reasons,
    requires_partial_member_checks,
)
from genosapiens.wgs_evidence_replay import EvidenceReplay
from genosapiens.wgs_evidence_store import EvidenceStoreError
from genosapiens.wgs_packed_outcomes import PackedOutcomes
from genosapiens.wgs_pair_outcomes import write_pair_outcomes, MAX_CATALOG_ENTRIES

MAX_ELIGIBILITY_BYTES = 64 * 1024
MAX_MEMBER_METADATA_BYTES = 64 * 1024**2
_STATES = ("supported","unresolved","contradicted")


class CompleteEligibilityEvidence:
    def __init__(self, writer, universe, population, representations, *, name):
        self.writer,self.connection,self.name = writer,writer.connection,name
        self.universe,self.population,self.representations = universe,population,representations
        self.prefix = "eligibility_" + hashlib.sha256(name.encode()).hexdigest()
        self.tables,self.closed,self.buffer,self.hot = [],False,None,OrderedDict()
        self.catalog_count,self.hypothesis_count,self.singleton_count,self.metadata_bytes = 0,0,0,0
        try:
            if (any(o.closed or o.connection is not self.connection for o in (universe,population,representations))
                    or population.universe is not universe or representations.members is not universe.members
                    or representations.facts is not universe.facts):
                raise EvidenceStoreError("foreign submission eligibility source owner")
            for suffix,definition in (("catalog","id INTEGER PRIMARY KEY,content TEXT UNIQUE,reasons TEXT"),
                                      ("singletons","ordinal INTEGER PRIMARY KEY,id TEXT UNIQUE,decision INTEGER")):
                table = self.prefix + "_" + suffix
                self.connection.execute(f"CREATE TABLE {table} ({definition})")
                self.tables.append(table)
            def singletons():
                for item in universe.singletons():
                    h = item.hypothesis
                    if item.admission_context_id != universe.members.context_id:
                        raise EvidenceStoreError("singleton eligibility context changed")
                    for ref in h.evidence_subject_refs:
                        if universe.facts.allele(ref.allele_id).native_identity_id != ref.identity_id:
                            raise EvidenceStoreError("singleton eligibility member changed")
                    assessed = population.assessment(h)
                    inputs = SubmissionEligibilityInputs(bool(h.locus_call_refs),len(h.allele_refs),
                        h.enumeration_status is EnumerationStatus.COMPLETE,True,True,item.model.state,
                        any(m.state == "nonqualifying" for m in item.members),item.read_validation_required,True,
                        assessed.outcome is PopulationOutcome.NOT_POPULATION_EXCLUDED,
                        requires_partial_member_checks(h),
                        len(h.allele_refs) == 1 and h.context.model == "incomplete" and item.incompleteness_penalty is not None,
                        tuple(submission_representation_reasons(ref,representations.get(ref.allele_id)) for ref in h.allele_refs))
                    code = self._evaluate(inputs,complete=h.completeness is Completeness.COMPLETE)
                    self.singleton_count += 1
                    self.connection.execute(f"INSERT INTO {self.prefix}_singletons VALUES (?,?,?)",
                                            (self.singleton_count,h.hypothesis_id,code))
                    yield {"ordinal":self.singleton_count,"hypothesis_id":h.hypothesis_id,"decision":code}
            singleton_sha = writer.add_sequence(name + ".singletons","singleton-eligibility-v1",singletons())
            self.hypothesis_count = self.singleton_count
            self.offset = 0
            pair_count,tape = 0,None
            if universe.pairs is not None:
                pairs = universe.pairs
                signatures = bytearray()
                for code,decision_id in self.connection.execute(f"SELECT id,decision FROM {pairs.prefix}_outcomes ORDER BY id"):
                    if code != len(signatures)//6:
                        raise EvidenceStoreError("eligibility pair outcome order changed")
                    decision = pairs.catalog.result(decision_id)
                    signatures.extend((int(decision.constructible),_STATES.index(decision.model_state),
                        int(decision.completeness is Completeness.INCOMPLETE),
                        int(any(s == "nonqualifying" for s in decision.qualification_states)),int(decision.read_validation_required),
                        int(decision.completeness is Completeness.COMPLETE)))
                if len(signatures) != pairs.outcome_count*6:
                    raise EvidenceStoreError("eligibility pair decision coverage changed")
                self.metadata_bytes = sys.getsizeof(signatures)
                members = {}
                context_key = ContextKey.from_context(pairs.context)
                for key,model,_ in pairs.member_templates():
                    reasons = submission_representation_reasons(model.allele,representations.get(key))
                    _,assessment = population.member(context_key,key)
                    if assessment.native_identity_id != model.allele.identity_id:
                        raise EvidenceStoreError("eligibility population member identity changed")
                    excluded = assessment.outcome is PopulationOutcome.POPULATION_EXCLUDED
                    self.metadata_bytes += sys.getsizeof(key)+sys.getsizeof(reasons)+sum(sys.getsizeof(r) for r in reasons)+160
                    if self.metadata_bytes > MAX_MEMBER_METADATA_BYTES:
                        raise EvidenceStoreError("eligibility member metadata exceeds capacity",capacity=True)
                    members[key] = reasons,excluded
                self.offset = int(pairs.pair_count != pairs.hypothesis_count)
                self.buffer = PackedOutcomes(self.connection,name)
                used = bytearray(self.catalog_count+self.offset)
                pair_hot = OrderedDict()
                for (left,right),code in zip(pairs.domain.pairs(),pairs._outcomes(),strict=True):
                    constructible,state,incomplete,nonqualifying,read_required,complete = signatures[6*code:6*code+6]
                    if not constructible:
                        self.buffer.append(0)
                        used[0] = 1
                        continue
                    lr,le = members[left]
                    rr,re = members[right]
                    key = state,incomplete,nonqualifying,read_required,complete,le or re,lr,rr
                    decision = pair_hot.get(key)
                    if decision is None:
                        inputs = SubmissionEligibilityInputs(False,2,True,True,True,_STATES[state],bool(nonqualifying),
                            bool(read_required),True,not (le or re),bool(incomplete),False,(lr,rr))
                        decision = self._evaluate(inputs,complete=bool(complete))
                        pair_hot[key] = decision
                        if len(pair_hot) > 4096:
                            pair_hot.popitem(last=False)
                    code = decision+self.offset
                    if code >= len(used):
                        used.extend(b"\0"*(code+1-len(used)))
                    used[code] = 1
                    self.buffer.append(code)
                    self.hypothesis_count += 1
                pair_count = self.buffer.count
                if pair_count != pairs.pair_count:
                    raise EvidenceStoreError("eligibility pair coverage changed")
                del members,signatures,pair_hot
            if self.hypothesis_count != universe.hypothesis_count or self.singleton_count != universe.singleton_count:
                raise EvidenceStoreError("complete eligibility hypothesis coverage changed")
            catalog_sha = writer.add_sequence(name + ".catalog","submission-eligibility-catalog-v1",self.catalog())
            if self.buffer is not None:
                # Pair catalog indices may refer to a singleton-created class.
                # The generic tape's usage closure requires the actual used
                # classes; persist its separate compact code mapping below.
                mapping = array("I",(0 for _ in used))
                count = 0
                for code,present in enumerate(used):
                    if present:
                        mapping[code] = count
                        count += 1
                if count > MAX_CATALOG_ENTRIES:
                    raise EvidenceStoreError("eligibility pair classes exceed capacity",capacity=True)
                outcome_sha = writer.add_sequence(name + ".pair_codes","pair-eligibility-codes-v1",
                    ({"ordinal":mapping[code],"decision":None if self.offset and code == 0 else code-self.offset}
                     for code,present in enumerate(used) if present))
                tape = write_pair_outcomes(writer,name + ".tape",(mapping[c] for c in self.buffer.values()),
                    domain=universe.pairs.summary,catalog_sha256=outcome_sha,catalog_count=count)
            self.sha256 = writer.add(name,"complete-submission-eligibility-v1",{
                "schema":"genosapiens.wgs-submission-eligibility","version":1,
                "assessment_method_sha256":METHOD_SHA256,"evaluable_domain_rules":EVALUABLE_DOMAIN_RULES,
                "hypotheses":{"name":universe.name,"sha256":universe.sha256},
                "population":{"name":population.name,"sha256":population.sha256},
                "representations":{"name":representations.name,"sha256":representations.sha256},
                "singletons":{"name":name + ".singletons","sha256":singleton_sha},
                "catalog":{"name":name + ".catalog","sha256":catalog_sha},
                "tape_name":None if tape is None else name + ".tape","tape":tape,
                "pair_codes":None if tape is None else {"name":name + ".pair_codes","sha256":outcome_sha},
                "singleton_count":self.singleton_count,"pair_count":pair_count,"hypothesis_count":self.hypothesis_count})
        except BaseException:
            writer.failed = True
            self.close()
            raise

    def _evaluate(self, inputs, *, complete):
        content = _canonical({"inputs":inputs,"complete":complete})
        if len(content) > MAX_ELIGIBILITY_BYTES:
            raise EvidenceStoreError("eligibility decision exceeds capacity",capacity=True)
        code = self.hot.get(content)
        if code is not None:
            return code
        row = self.connection.execute(f"SELECT id FROM {self.prefix}_catalog WHERE content=?",(content,)).fetchone()
        if row is None:
            if self.catalog_count >= MAX_CATALOG_ENTRIES:
                raise EvidenceStoreError("eligibility catalog exceeds capacity",capacity=True)
            code = self.catalog_count
            reasons = evaluate_submission_eligibility(inputs)
            self.connection.execute(f"INSERT INTO {self.prefix}_catalog VALUES (?,?,?)",(code,content,_canonical(reasons)))
            self.catalog_count += 1
        else:
            code = row[0]
        self.hot[content] = code
        while len(self.hot) > 256:
            self.hot.popitem(last=False)
        return code

    def catalog(self):
        for code,content,reasons in self.connection.execute(f"SELECT * FROM {self.prefix}_catalog ORDER BY id"):
            yield {"ordinal":code,**json.loads(content),"reasons":json.loads(reasons)}

    def reasons_at(self, locator, *, positions=None):
        if self.closed or type(locator) is not tuple or len(locator) != 2 or any(type(i) is not int for i in locator):
            raise EvidenceStoreError("invalid eligibility locator")
        if locator[1] == 0:
            row = self.connection.execute(f"SELECT decision FROM {self.prefix}_singletons WHERE ordinal=?",(locator[0],)).fetchone()
            if row is None:
                raise EvidenceStoreError("eligibility singleton is unavailable")
            code, = row
        else:
            if self.buffer is None or positions is None:
                raise EvidenceStoreError("eligibility pair position owner is required")
            code = self.buffer.value_at(positions.ordinal(*locator)-1)-self.offset
            if code < 0:
                raise EvidenceStoreError("eligibility locator has no hypothesis")
        row = self.connection.execute(f"SELECT reasons FROM {self.prefix}_catalog WHERE id=?",(code,)).fetchone()
        if row is None:
            raise EvidenceStoreError("eligibility decision is missing")
        return tuple(json.loads(row[0]))

    def close(self):
        if self.buffer is not None:
            self.buffer.close()
        for table in reversed(self.tables):
            self.connection.execute("DROP TABLE " + table)
        self.tables.clear()
        self.hot.clear()
        self.closed = True

    def __enter__(self):
        return self

    def __exit__(self,*_):
        self.close()


class VerifiedEligibilityEvidence(CompleteEligibilityEvidence):
    def __init__(self, reader, reference, universe, population, representations):
        if reference.kind != "complete-submission-eligibility-v1":
            raise EvidenceStoreError("foreign submission eligibility root")
        reader._require_reference(reference)
        super().__init__(EvidenceReplay(reader),universe,population,representations,name=reference.name)
