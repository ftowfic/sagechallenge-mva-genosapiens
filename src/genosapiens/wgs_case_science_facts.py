"""Bounded scalar access for diagnostics over complete computed science.

One row per singleton, member and evaluated pair class survives a context.
Private pair-code chunks are copied from the authenticated complete tape; no
per-pair JSON or global hypothesis-ID index is introduced. Saved projection
roots are recomputed from scientific owners when a case is opened.
"""
from collections import OrderedDict
import hashlib
from genosapiens.wgs_context_domain import next_context_ordinal, require_complete_domain
import json
import sys

from genosapiens.evidence_serialization import _canonical
from genosapiens.evidence_model import _sha
from genosapiens.population_policy import ContextKey, PopulationOutcome, summarize_population_values
from genosapiens.wgs_evidence_store import EvidenceStoreError
from genosapiens.wgs_packed_outcomes import outcome_chunk_values
from genosapiens.wgs_pair_outcomes import CHUNK_PAIRS

MAX_SCALAR_BYTES = 64 * 1024
MAX_PAIR_CACHE_BYTES = 64 * 1024**2
MAX_CLASS_CACHE_BYTES = 16 * 1024**2


def _population_member(member):
    return (member.penalty,member.outcome.value,
        member.status == "assessed" and any(f.reliable for f in member.frequencies)
        and any(r.status == "unknown" for r in member.rescues))


def _population_values(members):
    penalty,contribution,outcome = summarize_population_values(tuple((p,PopulationOutcome(o)) for p,o,_ in members))
    return {"population_penalty":penalty,"population_contribution":contribution,
        "population_outcome":outcome.value,"unknown_rescue_members":sum(int(u) for _,_,u in members),
        "population_excluded_members":tuple(o == "population_excluded" for _,o,_ in members)}


def constructed_scalar_facts(constructed, population):
    """Direct projection used for singletons and independent rich-row checks."""
    return {"model_state":constructed.model.state,"completeness":constructed.hypothesis.completeness.value,
        "member_states":tuple(m.state for m in constructed.members),
        "clinical_support":tuple(bool(m.clinical_support_refs) for m in constructed.members),
        "alternate_present":tuple(m.genotype.alternate_present for m in constructed.members),
        "conditional_singleton_search_satisfied":constructed.conditional_singleton_search_satisfied,
        "read_validation_required":constructed.read_validation_required,
        **_population_values(tuple(_population_member(m) for m in population.members))}


def _singleton_scalar_rows(context):
    for ordinal,item in enumerate(context.universe.singletons(),1):
        yield ordinal,constructed_scalar_facts(item,context.population.assessment(item.hypothesis))


def _member_scalar_rows(context):
    pairs = context.universe.pairs
    if pairs is not None:
        key = ContextKey.from_context(pairs.context)
        for ordinal,allele_id,_,_,parent in pairs.domain.members():
            if not parent:
                _,assessment = context.population.member(key,allele_id)
                yield ordinal,_population_member(assessment)


def _class_scalar_rows(context):
    pairs = context.universe.pairs
    if pairs is not None:
        for code,decision_id in context.connection.execute(f"SELECT id,decision FROM {pairs.prefix}_outcomes ORDER BY id"):
            decision = pairs.catalog.result(decision_id)
            left,right,_ = pairs.catalog.outcome_inputs(decision_id)
            dependencies = tuple(pairs.catalog.member_inputs(k) for k in (left,right))
            yield code,{"constructible":decision.constructible,"model_state":decision.model_state,
                "completeness":decision.completeness.value,"member_states":decision.qualification_states,
                "clinical_support":tuple(m.qualification.support for m in dependencies),
                "alternate_present":tuple(m.model.genotype.alternate_present for m in dependencies),
                "conditional_singleton_search_satisfied":False,
                "read_validation_required":decision.read_validation_required}


def _scalar_context_header(context):
    from genosapiens.wgs_population_work_filter import work_domain
    return {"ordinal":context.ordinal,"context":context.reference,"hypothesis_count":context.hypothesis_count,
        **work_domain(context.members),
        "assessment_domain_sha256":context.assessment_domain_sha256}


def _scalar_pair_rows(context):
    pairs = context.universe.pairs
    if pairs is not None:
        yield from context.connection.execute(f"SELECT ordinal,count,content FROM {pairs.prefix}_chunks ORDER BY ordinal")


class CompleteCaseScienceFacts:
    def __init__(self, writer, inputs, *, name, assessment_domain_sha256, context_domain=None):
        self.context_domain = context_domain
        self.writer,self.connection,self.inputs,self.name = writer,writer.connection,inputs,name
        self.prefix = "case_science_facts_"+hashlib.sha256(name.encode()).hexdigest()
        self.tables,self.closed,self.published = [],False,False
        self.context_count,self.hypothesis_count = 0,0
        self.cache,self.cache_size = OrderedDict(),{"pairs":0,"classes":0}
        try:
            _sha(assessment_domain_sha256)
            self.assessment_domain_sha256 = assessment_domain_sha256
            if inputs.closed or inputs.connection is not self.connection:
                raise EvidenceStoreError("foreign scalar scientific input owner")
            for suffix,definition in (
                ("contexts","ordinal INTEGER PRIMARY KEY,content TEXT"),
                ("singletons","context INTEGER,ordinal INTEGER,content TEXT,PRIMARY KEY(context,ordinal)"),
                ("members","context INTEGER,ordinal INTEGER,content TEXT,PRIMARY KEY(context,ordinal)"),
                ("classes","context INTEGER,ordinal INTEGER,content TEXT,PRIMARY KEY(context,ordinal)"),
                ("pairs","context INTEGER,ordinal INTEGER,count INTEGER,content BLOB,PRIMARY KEY(context,ordinal)"),
            ):
                table = self.prefix+"_"+suffix
                self.connection.execute(f"CREATE TABLE {table} ({definition})")
                self.tables.append(table)
        except BaseException:
            writer.failed = True
            self.close()
            raise

    def _store(self, kind, context, ordinal, value):
        content = _canonical(value)
        if len(content) > MAX_SCALAR_BYTES:
            raise EvidenceStoreError("scientific scalar row exceeds capacity",capacity=True)
        self.connection.execute(f"INSERT INTO {self.prefix}_{kind} VALUES (?,?,?)",(context,ordinal,content))
        return ordinal,value

    def add_context(self, context):
        self._require_open()
        try:
            if (self.published or context.closed or context.connection is not self.connection
                    or context.assessment_domain_sha256 != self.assessment_domain_sha256
                    or context.ordinal != next_context_ordinal(self) or self.inputs.context_count != self.context_count + 1):
                raise EvidenceStoreError("foreign scalar scientific context")
            self._append_context(_scalar_context_header(context),_singleton_scalar_rows(context),
                _member_scalar_rows(context),_class_scalar_rows(context),_scalar_pair_rows(context))
        except BaseException:
            self.writer.failed = True
            raise

    def _append_context(self, header, singleton_rows, member_rows, class_rows, pair_rows):
        self._require_open()
        ordinal = header["ordinal"]
        if (self.published or header["assessment_domain_sha256"] != self.assessment_domain_sha256
                or ordinal != next_context_ordinal(self) or self.inputs.context_count != self.context_count + 1):
            raise EvidenceStoreError("foreign scalar scientific context")
        try:
            name = self.name+".context-"+str(ordinal)
            def rows(kind, values):
                for number,value in values:
                    yield self._store(kind,ordinal,number,value)
            singleton_sha = self.writer.add_sequence(name+".singletons","scientific-singleton-facts-v1",rows("singletons",singleton_rows))
            members_sha = self.writer.add_sequence(name+".members","scientific-pair-member-facts-v1",rows("members",member_rows))
            classes_sha = self.writer.add_sequence(name+".classes","scientific-pair-class-facts-v1",rows("classes",class_rows))
            self.connection.executemany(f"INSERT INTO {self.prefix}_pairs VALUES (?,?,?,?)",
                ((ordinal,*row) for row in pair_rows))
            value = {**header,
                "singletons":{"name":name+".singletons","sha256":singleton_sha},
                "members":{"name":name+".members","sha256":members_sha},
                "classes":{"name":name+".classes","sha256":classes_sha}}
            self.connection.execute(f"INSERT INTO {self.prefix}_contexts VALUES (?,?)",(ordinal,_canonical(value)))
            self.context_count += 1
            self.hypothesis_count += header["hypothesis_count"]
        except BaseException:
            self.writer.failed = True
            raise

    def publish(self):
        require_complete_domain(self)
        self._require_open()
        if self.published or not self.inputs.published:
            raise EvidenceStoreError("scientific scalar publication is unavailable")
        digest = self.writer.add_sequence(self.name+".contexts","scientific-scalar-contexts-v1",
            (json.loads(content) for content, in self.connection.execute(f"SELECT content FROM {self.prefix}_contexts ORDER BY ordinal")))
        self.sha256 = self.writer.add(self.name,"case-scientific-scalar-facts-v1",{
            "schema":"genosapiens.wgs-case-scientific-scalar-facts","version":1,
            "assessment_domain_sha256":self.assessment_domain_sha256,
            "context_count":self.context_count,"hypothesis_count":self.hypothesis_count,
            "inputs":self.inputs.reference,"contexts":{"name":self.name+".contexts","sha256":digest}})
        self.published = True

    def _row(self, kind, context, ordinal):
        key = kind,context,ordinal
        cached = self.cache.get(key)
        if cached is not None:
            self.cache.move_to_end(key)
            return cached[0]
        fields = "count,content" if kind == "pairs" else "content"
        row = self.connection.execute(f"SELECT {fields} FROM {self.prefix}_{kind} WHERE context=? AND ordinal=?",(context,ordinal)).fetchone()
        if row is None:
            raise EvidenceStoreError("scientific scalar locator is missing")
        value = outcome_chunk_values(*row) if kind == "pairs" else json.loads(row[0])
        if kind not in self.cache_size:
            return value
        size = sys.getsizeof(value)+256 if kind == "pairs" else len(row[0])*8+2048
        limit = MAX_PAIR_CACHE_BYTES if kind == "pairs" else MAX_CLASS_CACHE_BYTES
        if size > limit:
            raise EvidenceStoreError("scientific scalar cache window exceeds capacity",capacity=True)
        while self.cache_size[kind]+size > limit:
            for old in self.cache:
                if old[0] == kind:
                    break
            _,old_size = self.cache.pop(old)
            self.cache_size[kind] -= old_size
        self.cache[key] = value,size
        self.cache_size[kind] += size
        return value

    def get(self, locator):
        self._require_open()
        if (not self.published or type(locator) is not tuple or len(locator) != 3
                or any(type(v) is not int for v in locator)):
            raise EvidenceStoreError("invalid scientific scalar locator")
        context,left,right = locator
        if right == 0:
            value = dict(self._row("singletons",context,left))
        else:
            positions = self.inputs._context(context).positions
            if positions is None:
                raise EvidenceStoreError("scientific scalar pair domain is missing")
            position = positions.ordinal(left,right)-1
            codes = self._row("pairs",context,position//CHUNK_PAIRS+1)
            if position%CHUNK_PAIRS >= len(codes):
                raise EvidenceStoreError("scientific scalar pair position is missing")
            value = dict(self._row("classes",context,codes[position%CHUNK_PAIRS]))
            if not value.pop("constructible"):
                raise EvidenceStoreError("scientific scalar pair has no hypothesis")
            value.update(_population_values(tuple(self._row("members",context,k) for k in (left,right))))
        for key in ("member_states","clinical_support","alternate_present","population_excluded_members"):
            value[key] = tuple(value[key])
        return value

    def pair_records(self, context):
        """All evaluated pairs, including outcomes that emitted no hypothesis."""
        self._require_open()
        if not self.published:
            raise EvidenceStoreError("scientific scalar facts are incomplete")
        positions = self.inputs._context(context).positions
        if positions is None:
            return
        observed = 0
        for chunk, count, content in self.connection.execute(
                f"SELECT ordinal,count,content FROM {self.prefix}_pairs WHERE context=? ORDER BY ordinal", (context,)):
            if chunk != observed // CHUNK_PAIRS + 1:
                raise EvidenceStoreError("scientific pair scalar coverage changed")
            for code in outcome_chunk_values(count, content):
                observed += 1
                value = dict(self._row("classes", context, code))
                for key in ("member_states", "clinical_support", "alternate_present"):
                    value[key] = tuple(value[key])
                yield positions.locator(observed), value
        if observed != positions.prefix[-1]:
            raise EvidenceStoreError("scientific pair scalar count changed")

    def _require_open(self):
        if self.closed or self.writer.connection is not self.connection or self.inputs.closed:
            raise EvidenceStoreError("scientific scalar owner is closed")

    def close(self):
        self.cache.clear()
        for table in reversed(self.tables):
            self.connection.execute("DROP TABLE "+table)
        self.tables.clear()
        self.closed = True

    def __enter__(self):
        return self

    def __exit__(self,*_):
        self.close()


class CompletedCaseScienceFacts(CompleteCaseScienceFacts):
    """Read stored scalar facts and evaluated pair codes, never score science."""
    def __init__(self, reader, reference, inputs, *, assessment_domain_sha256, context_at, progress=None, context_domain=None):
        from genosapiens.wgs_completed_science import child, fields, header, integer, require, outcome_chunks, observed_rows, scalar_value
        from genosapiens.wgs_evidence_replay import EvidenceReplay
        from genosapiens.wgs_case_resources import bound_reference
        value = header(reader, reference, kind='case-scientific-scalar-facts-v1',
            schema='genosapiens.wgs-case-scientific-scalar-facts', keys=('assessment_domain_sha256',
                'context_count', 'hypothesis_count', 'inputs', 'contexts'))
        require(value['assessment_domain_sha256'] == assessment_domain_sha256 and value['inputs'] == inputs.reference)
        super().__init__(EvidenceReplay(reader), inputs, name=reference.name, assessment_domain_sha256=assessment_domain_sha256, context_domain=context_domain)
        try:
            for row in observed_rows(reader.iter_array(child(reader, value['contexts'], 'scientific-scalar-contexts-v1'), max_item_bytes=MAX_SCALAR_BYTES),
                    progress, 'saved_scalar_facts', value['context_count']):
                fields(row, ('ordinal', 'context', 'hypothesis_count', 'assessment_domain_sha256', 'singletons', 'members', 'classes',
                    *(('population_work',) if 'population_work' in row else ())))
                ordinal = integer(row['ordinal'], minimum=1)
                require(ordinal == next_context_ordinal(self) and row['assessment_domain_sha256'] == assessment_domain_sha256)
                context = context_at(ordinal)
                require(row['context'] == context['context'] and row['hypothesis_count'] == context['hypothesis_count'])
                counts = {}
                for key, kind in (('singletons', 'scientific-singleton-facts-v1'), ('members', 'scientific-pair-member-facts-v1'),
                                  ('classes', 'scientific-pair-class-facts-v1')):
                    last = -1 if key == 'classes' else 0
                    counts[key] = 0
                    for item in reader.iter_array(child(reader, row[key], kind), max_item_bytes=MAX_SCALAR_BYTES):
                        require(type(item) is list and len(item) == 2)
                        number = integer(item[0])
                        require(number > last and (key == 'members' or number == last + 1))
                        last = number
                        scalar_value(key, item[1])
                        self._store(key, ordinal, number, item[1])
                        counts[key] += 1
                context_header = reader.reconstruct(bound_reference(reader, row['context'], 'scientific-case-context-v1'))
                from genosapiens.wgs_population_work_filter import verify_saved_work_domain
                verify_saved_work_domain(reader, context_header, row)
                universe = reader.reconstruct(bound_reference(reader, context_header['universe'], 'complete-hypothesis-universe-v1'))
                require(counts['singletons'] == integer(universe['singleton_count']))
                if universe['pair_evidence'] is not None:
                    pair = reader.reconstruct(child(reader, universe['pair_evidence'], 'complete-pair-evidence-v1'))
                    require(pair['pair_count'] == context['pair_count'])
                    require(counts['classes'] == integer(pair['outcome_count']))
                    child(reader, pair['outcomes'], 'pair-computed-outcomes-v1')
                    self.connection.executemany(f'INSERT INTO {self.prefix}_pairs VALUES (?,?,?,?)',
                        ((ordinal, *chunk) for chunk in outcome_chunks(reader, pair['tape_name'], pair['tape'],
                            catalog_sha256=pair['outcomes']['sha256'], catalog_count=pair['outcome_count'],
                            pair_count=pair['pair_count'], domain_sha256=pair['domain_summary']['membership_sha256'])))
                else:
                    require(context['pair_count'] == counts['members'] == counts['classes'] == 0)
                self.connection.execute(f'INSERT INTO {self.prefix}_contexts VALUES (?,?)', (ordinal, _canonical(row)))
                self.context_count += 1
                self.hypothesis_count += integer(row['hypothesis_count'])
            require(all(getattr(self, key) == integer(value[key]) for key in ('context_count', 'hypothesis_count')))
            self.sha256, self.published = reference.sha256, True
        except BaseException:
            self.close()
            raise
