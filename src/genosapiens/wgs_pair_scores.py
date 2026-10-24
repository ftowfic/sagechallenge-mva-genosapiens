"""Complete pair score outcomes over authenticated computed source families."""
from collections import OrderedDict
from dataclasses import dataclass
import hashlib
import json
import sys

from genosapiens.evidence_model import Completeness, HypothesisPhase, _typed
from genosapiens.evidence_serialization import _canonical
from genosapiens.population_policy import ContextKey, PopulationOutcome, summarize_population_values
from genosapiens.wgs_evidence_replay import EvidenceReplay
from genosapiens.wgs_evidence_store import EvidenceStoreError
from genosapiens.wgs_packed_outcomes import PackedOutcomes
from genosapiens.wgs_pair_outcomes import write_pair_outcomes
from genosapiens.wgs_score_decisions import ScoreChargeValue, ScoreDecisionCatalog, MAX_SCORE_INPUT_BYTES
from genosapiens.wgs_score_sources import pair_score_inputs, pair_score_numeric_inputs

MAX_MEMBER_METADATA_BYTES = 64 * 1024**2
HOT_FAMILIES = 4096
MAX_FAMILY_CACHE_ROWS = 1_000_000
MAX_HOT_FAMILY_BYTES = 16 * 1024**2
_STATES = ("supported", "unresolved", "contradicted")
_COMPLETENESS = tuple(Completeness)
_PHASES = tuple(HypothesisPhase)


@dataclass(frozen=True, slots=True, repr=False)
class _MemberValues:
    nontechnical: tuple[ScoreChargeValue, ...]
    technical: tuple[str | None, ...]
    population_penalty: str | None
    population_outcome: PopulationOutcome


class CompletePairScores:
    def __init__(self, writer, universe, population, sources, *, phenotype_score, name):
        self.writer, self.connection, self.name = writer, writer.connection, name
        self.universe, self.population, self.sources = universe, population, sources
        self.tables, self.closed = [], False
        self.catalog = self.buffer = None
        self.prefix = "pair_scores_" + hashlib.sha256(name.encode()).hexdigest()
        self.member_values, self.family_hot = {}, OrderedDict()
        self.class_count = self.hypothesis_count = self.family_evaluations = 0
        self.family_rows = self.family_hot_bytes = 0
        self.metadata_bytes = 0
        if (universe.pairs is None or universe.closed or population.closed or sources.closed
                or population.universe is not universe or sources.universe is not universe
                or any(owner.connection is not self.connection for owner in (universe,population,sources))):
            raise EvidenceStoreError("foreign complete pair score owners")
        try:
            # The surrounding source adapter owns phenotype resource replay.
            # This exact shared transform checks scorer identity/raw units.
            sources.policy.phenotype.component(phenotype_score)
            if len(phenotype_score.as_bytes()) > MAX_MEMBER_METADATA_BYTES:
                raise EvidenceStoreError("pair phenotype evidence exceeds capacity",capacity=True)
            self.phenotype_units = phenotype_score.ledger["raw_evidence_units"]
            phenotype_sha = writer.add(name + ".phenotype", "pair-score-phenotype-v1", phenotype_score.ledger)
            for suffix, definition in (
                ("classes", "ordinal INTEGER PRIMARY KEY,digest TEXT UNIQUE,content TEXT"),
                ("technical", "lineage TEXT,key TEXT,ordinal INTEGER,PRIMARY KEY(lineage,key)"),
                ("families", "key TEXT PRIMARY KEY,score INTEGER"),
            ):
                table = self.prefix + "_" + suffix
                self.connection.execute(f"CREATE TABLE {table} ({definition})")
                self.tables.append(table)
            self._prepare_members()
            signatures = self._pair_signatures()
            self.catalog = ScoreDecisionCatalog(writer,sources.policy,policy_sha256=sources.policy_sha256,
                                                name=name + ".decisions")
            self.buffer = PackedOutcomes(self.connection,name)
            self.offset = int(universe.pairs.pair_count != universe.pairs.hypothesis_count)
            for (left,right), code in zip(universe.pairs.domain.pairs(),universe.pairs._outcomes(),strict=True):
                constructible,state,complete,phase = signatures[4*code:4*code+4]
                if not constructible:
                    self.buffer.append(0)
                    continue
                left_class,left_order = self.member_values[left]
                right_class,right_order = self.member_values[right]
                ordering = tuple((scope,position) for _,scope,position in sorted(
                    [(order,0,i) for i,order in enumerate(left_order)] +
                    [(order,1,i) for i,order in enumerate(right_order)] +
                    [(order,2,i) for i,order in enumerate(self.shared_order)]))
                key = (left_class,right_class,state,complete,phase,ordering)
                # The all-disjoint result was computed from every actual source
                # charge. Otherwise check this actual pair even on a cache hit.
                checked_inputs = None
                if not sources.all_owners_disjoint:
                    left_profile,right_profile = sources.profile(left),sources.profile(right)
                    lv,rv = self._class(left_class),self._class(right_class)
                    _,contribution,_ = summarize_population_values(((lv.population_penalty,lv.population_outcome),
                                                                    (rv.population_penalty,rv.population_outcome)))
                    checked_inputs = pair_score_inputs(left_profile,right_profile,sources.shared_profile,
                        completeness=_COMPLETENESS[complete],phase=_PHASES[phase],phenotype_units=self.phenotype_units,
                        model_state=_STATES[state],population_contribution=contribution,policy=sources.policy)
                cached = self.family_hot.get(key)
                score = None if cached is None else cached[0]
                if score is None:
                    encoded_key = _canonical(key)
                    if len(encoded_key) > MAX_SCORE_INPUT_BYTES:
                        raise EvidenceStoreError("pair score family key exceeds capacity",capacity=True)
                    row = self.connection.execute(f"SELECT score FROM {self.prefix}_families WHERE key=?",(encoded_key,)).fetchone()
                    if row is None:
                        inputs = checked_inputs or self._inputs(left_class,right_class,state,complete,phase,ordering)
                        score = self.catalog.evaluate(inputs)
                        if self.family_rows == MAX_FAMILY_CACHE_ROWS:
                            # This is only a private memo. Complete computed
                            # numeric outcomes and pair assignments are retained.
                            self.connection.execute(f"DELETE FROM {self.prefix}_families")
                            self.family_rows = 0
                        self.connection.execute(f"INSERT INTO {self.prefix}_families VALUES (?,?)",(encoded_key,score))
                        self.family_rows += 1
                        self.family_evaluations += 1
                    else:
                        score = row[0]
                    size = len(encoded_key) + sys.getsizeof(key) + sys.getsizeof(ordering) + sum(sys.getsizeof(v) for v in ordering) + 256
                    self.family_hot[key] = score,size
                    self.family_hot_bytes += size
                    while len(self.family_hot) > HOT_FAMILIES or self.family_hot_bytes > MAX_HOT_FAMILY_BYTES:
                        _,(_,removed_size) = self.family_hot.popitem(last=False)
                        self.family_hot_bytes -= removed_size
                self.buffer.append(score + self.offset)
                self.hypothesis_count += 1
            if (self.buffer.count != universe.pairs.pair_count or self.hypothesis_count != universe.pairs.hypothesis_count):
                raise EvidenceStoreError("complete pair score coverage changed")
            catalog_sha = self.catalog.publish()
            def outcomes():
                if self.offset:
                    yield {"ordinal":0,"score":None}
                for i in range(self.catalog.count):
                    yield {"ordinal":i + self.offset,"score":i}
            outcomes_sha = writer.add_sequence(name + ".outcomes", "pair-score-outcomes-v1", outcomes())
            self.tape = write_pair_outcomes(writer,name + ".tape",self.buffer.values(),domain=universe.pairs.summary,
                catalog_sha256=outcomes_sha,catalog_count=self.catalog.count + self.offset)
            self.sha256 = writer.add(name,"complete-pair-scores-v1",{
                "schema":"genosapiens.wgs-complete-pair-scores","version":1,
                "hypotheses":{"name":universe.name,"sha256":universe.sha256},
                "population":{"name":population.name,"sha256":population.sha256},
                "sources":{"name":sources.name,"sha256":sources.sha256},
                "phenotype":{"name":name + ".phenotype","sha256":phenotype_sha},
                "decisions":{"name":self.catalog.name,"sha256":catalog_sha},
                "outcomes":{"name":name + ".outcomes","sha256":outcomes_sha},
                "tape_name":name + ".tape","tape":self.tape,"no_hypothesis_code":0 if self.offset else None,
                "pair_count":self.buffer.count,"hypothesis_count":self.hypothesis_count})
        except BaseException:
            writer.failed = True
            self.close()
            raise

    def _prepare_members(self):
        context_key = ContextKey.from_context(self.universe.pairs.context)
        for _,profile in self.sources.profiles():
            self.connection.executemany(f"INSERT OR IGNORE INTO {self.prefix}_technical VALUES (?,?,NULL)",
                ((lineage,key) for lineage,key,_ in profile.technical))
        self.connection.executemany(f"INSERT OR IGNORE INTO {self.prefix}_technical VALUES (?,?,NULL)",
            ((lineage,key) for lineage,key,_ in self.sources.shared_profile.technical))
        for ordinal,(lineage,key) in enumerate(self.connection.execute(
                f"SELECT lineage,key FROM {self.prefix}_technical ORDER BY lineage,key")):
            self.connection.execute(f"UPDATE {self.prefix}_technical SET ordinal=? WHERE lineage=? AND key=?",(ordinal,lineage,key))
        def order(profile):
            return tuple(self.connection.execute(f"SELECT ordinal FROM {self.prefix}_technical WHERE lineage=? AND key=?",
                (lineage,key)).fetchone()[0] for lineage,key,_ in profile.technical)
        self.shared_order = order(self.sources.shared_profile)
        self.metadata_bytes += sys.getsizeof(self.shared_order) + sum(sys.getsizeof(i) for i in self.shared_order)
        for key,profile in self.sources.profiles():
            _,population = self.population.member(context_key,key)
            value = _MemberValues(profile.nontechnical,tuple(v for _,_,v in profile.technical),
                                  population.penalty,population.outcome)
            content = _canonical(value)
            if len(content) > MAX_SCORE_INPUT_BYTES:
                raise EvidenceStoreError("member score values exceed capacity",capacity=True)
            digest = hashlib.sha256(content.encode("ascii")).hexdigest()
            row = self.connection.execute(f"SELECT ordinal,content FROM {self.prefix}_classes WHERE digest=?",(digest,)).fetchone()
            if row is None:
                ordinal = self.class_count
                self.connection.execute(f"INSERT INTO {self.prefix}_classes VALUES (?,?,?)",(ordinal,digest,content))
                self.class_count += 1
            else:
                ordinal,previous = row
                if previous != content:
                    raise EvidenceStoreError("member score value digest collision")
            indices = order(profile)
            self.metadata_bytes += sys.getsizeof(key) + sys.getsizeof(indices) + sum(sys.getsizeof(i) for i in indices) + 128
            if self.metadata_bytes > MAX_MEMBER_METADATA_BYTES:
                raise EvidenceStoreError("pair score member metadata exceeds capacity",capacity=True)
            self.member_values[key] = ordinal,indices

    def _pair_signatures(self):
        pairs = self.universe.pairs
        result = bytearray()
        for code,decision_id,pattern_id in self.connection.execute(
                f"SELECT id,decision,pattern FROM {pairs.prefix}_outcomes ORDER BY id"):
            if code != len(result)//4:
                raise EvidenceStoreError("pair score source outcomes are incomplete")
            decision = pairs.catalog.result(decision_id)
            pattern, = self.connection.execute(f"SELECT content FROM {pairs.prefix}_patterns WHERE id=?",(pattern_id,)).fetchone()
            phase = HypothesisPhase(json.loads(pattern)["phase"])
            result.extend((int(decision.constructible),_STATES.index(decision.model_state),
                           _COMPLETENESS.index(decision.completeness),_PHASES.index(phase)))
        if len(result) != pairs.outcome_count*4:
            raise EvidenceStoreError("pair score source outcome count changed")
        return result

    def _class(self, ordinal):
        content, = self.connection.execute(f"SELECT content FROM {self.prefix}_classes WHERE ordinal=?",(ordinal,)).fetchone()
        return _typed(json.loads(content),_MemberValues,decode=True)

    def _inputs(self, left_class, right_class, state, complete, phase, ordering):
        left,right = self._class(left_class),self._class(right_class)
        shared = self.sources.shared_profile
        charges = [ScoreChargeValue(c.kind,c.value,slot) for slot,profile in (("0",left),("1",right),(None,shared))
                   for c in profile.nontechnical]
        technical = left.technical,right.technical,tuple(v for _,_,v in shared.technical)
        charges.extend(ScoreChargeValue("technical",technical[scope][position],None if scope==2 else str(scope))
                       for scope,position in ordering)
        _,contribution,_ = summarize_population_values(((left.population_penalty,left.population_outcome),
                                                        (right.population_penalty,right.population_outcome)))
        return pair_score_numeric_inputs(tuple(charges),completeness=_COMPLETENESS[complete],phase=_PHASES[phase],
            phenotype_units=self.phenotype_units,model_state=_STATES[state],population_contribution=contribution,
            policy=self.sources.policy)

    def decisions(self):
        if self.closed:
            raise EvidenceStoreError("complete pair score owner is closed")
        for code in self.buffer.values():
            yield None if self.offset and code == 0 else self.catalog.result(code-self.offset)

    def close(self):
        if self.catalog is not None:
            self.catalog.close()
        if self.buffer is not None:
            self.buffer.close()
        for table in reversed(self.tables):
            self.connection.execute("DROP TABLE " + table)
        self.tables.clear()
        self.member_values.clear()
        self.family_hot.clear()
        self.family_hot_bytes = 0
        self.closed = True

    def __enter__(self):
        return self

    def __exit__(self,*_):
        self.close()


class VerifiedPairScores(CompletePairScores):
    def __init__(self, reader, reference, universe, population, sources, *, phenotype_score):
        if reference.kind != "complete-pair-scores-v1":
            raise EvidenceStoreError("foreign complete pair score root")
        reader._require_reference(reference)
        super().__init__(EvidenceReplay(reader),universe,population,sources,phenotype_score=phenotype_score,name=reference.name)
